parallel_spec_dirs <- function(spec.dirs, dir, shape, extra = NULL, hv.project.method, fun.init, fun.execute, batch = FALSE, test = 0, node.log.append = TRUE, verbose = FALSE) {
  tryCatch(
    {
      if (verbose) {
        print_function_args()
      }

      catn("Setting up parallel function")

      lock_dir <- paste0(dir, "/locks")
      create_dir_if(lock_dir)

      ram_log <- paste0(dir, "/ram-usage.txt")
      peak_log <- paste0(dir, "/ram-peak.txt")
      stop_file <- paste0(dir, "/stop-file.txt")
      highest_it_file <- paste0(dir, "/highest-iteration.txt")

      create_file_if(ram_log, peak_log)
      create_file_if(highest_it_file, keep = TRUE)

      node_dir <- paste0(dir, "/nodes")
      create_dir_if(node_dir)
      invisible(gc()) # Garbage collect to avoid garbage collection issues when getting peak ram

      # Get max core usage
      used_mem <- get_mem_usage(type = "used", format = "gb")
      remain_mem <- ((config$memory$mem_limit / 1024^3) - used_mem)

      # Set up names for the initiation
      sp_dir <- spec.dirs[1]
      sp_dirs <- spec.dirs[-1]
      if (file.info(sp_dir[[1]])$isdir) {
        sp_filename <- paste0(sp_dir, "/", hv.project.method, ".tif")
      } else {
        sp_filename <- sp_dir
      }

      catn("Processing init function.")
      init_timer <- start_timer("Inititation process")

      init_res <- parallel_init(
        spec.filename = sp_filename,
        shape = shape,
        extra = extra,
        fun.init = fun.init,
        file.log = peak_log,
        file.stop = stop_file,
        verbose = verbose
      )

      if (is.null(init_res)) {
        vebcat("Init process did not return correctly.", color = "fatalError")
        stop("Check init parallel_init function in the parallel_handler")
      }

      if (length(sp_dirs) == 0) {
        catn(highcat(length(sp_dirs)), "species left after the initiation, returning.")
        return(list(init.res = init_res, exec.res = NULL))
      }

      timer_res <- end_timer(init_timer)

      peak_ram <- as.numeric(readLines(peak_log))

      max_cores <- calc_num_cores(
        ram.high = peak_ram,
        verbose = verbose
      )
      max_cores <- max_cores$high

      if (test > 0) { # use test numbers if above 0
        max_cores <- min(test, max_cores, config$memory$total_cores)
      } else {
        max_cores <- min(length(sp_dirs), max_cores, config$memory$total_cores)
      }

      catn("Creating cluster of", max_cores, "core(s).")

      config <- config

      cl <- makeCluster(max_cores)

      # Get the exports that are not null
      export_vars <- c("sp_dirs", "shape", "hv.project.method", "fun.execute", "batch", "node.log.append", "test", "verbose", "node_dir", "ram_log", "extra", "config")

      # export_vars <- export_vars[sapply(export_vars, function(x) !is.null(get(x)))]

      vebprint(export_vars, verbose, "export_vars:")

      clusterExport(cl, export_vars, envir = environment())

      tryCatch(
        {
          clusterEvalQ(cl, {
            global_config <- config
            source("./R/utils/utils.R")
            load_utils(parallel = TRUE)

            config <- global_config

            source("./R/setup/setup_sequence.R")
            source("./R/hypervolume/components/spec_handlers.R")

            source_all("./R/visualize/components")
          })
        },
        error = function(e) {
          vebcat("An error occurred in the parallel process ~ stopping cluster and closing connections.", color = "fatalError")
          stopCluster(cl)
          closeAllConnections()
          stop(e$message)
        }
      )

      current_disk_space <- get_disk_space("/export", units = "GB")

      # Check for batched approach
      if (batch) {
        # Calculate the number of species per batch
        sp_per_core <- ceiling(length(sp_dirs) / max_cores)

        # Split tasks into batches
        sp_dirs <- split(sp_dirs, ceiling(seq_along(sp_dirs) / sp_per_core))

        catn("Directories split into", highcat(sp_per_core), "species per core.")
      }

      # If iterations is not provided, start from the highest saved iteration
      i <- 1
      if (test > 0) end <- test else end <- length(sp_dirs)
      iterations <- i:end

      ram_msg <- FALSE
      # RAM check
      mem_used_gb <- get_mem_usage(type = "used", format = "gb")
      mem_limit_gb <- config$memory$mem_limit / 1024^3

      catn("Running parallel with", highcat(end), "iteration(s).")
      execute_timer <- start_timer("Execute process")
      calculate_etc(timer_res + (1.45 * max_cores), max_cores, length(sp_dirs))

      tryCatch(
        {
          res <- clusterApplyLB(cl, iterations, function(i) {
            repeat {
              if (!mem_check(
                identifier = paste("node", i),
                ram.use = ram_log,
                verbose = verbose
              )) {
                break # Exit loop when memory is OK
              }
            }

            # Setup node log
            node_log <- paste0(node_dir, "/node-", i, ".txt")
            create_file_if(node_log)

            if (node.log.append) {
              node_con <- file(node_log, open = "at")
            } else {
              node_con <- file(node_log, open = "wt")
            }
            if (!inherits(node_con, "connection")) {
              stop("Failed to open log file: ", node_log)
            }

            sink(node_con, type = "output")
            sink(node_con, type = "message")

            sp_dir <- sp_dirs[[i]]
            sp_name <- basename(sp_dir)
            if (file.info(sp_dir[[1]])$isdir) {
              sp_filename <- paste0(sp_dir, "/", hv.project.method, ".tif")
            } else {
              sp_filename <- sp_dir
            }

            catn("Running node for:")

            if (!batch) {
              catn("Iteration:", i)
              catn("Species:", sp_name)
              catn("Species Directory:\n", sp_dir)
              catn("Species Filename:\n", sp_filename)
              catn()
            } else {
              if (test > 0) sp_filename <- sp_filename[1:3]
              catn("Batch", i)
              catn("Number of species:", length(sp_dir))
              catn("Number of filenames:", length(sp_filename))
              catn()
            }

            vebprint(shape, verbose, "Shape:")
            if (!is.null(shape)) {
              region <- load_region(shape)
              region <- handle_region(region)
            } else {
              region <- NULL
            }

            vebprint(extra, verbose, "Extra:")

            catn("\nParallel processing execute function.")
            dt <- track_memory(function() {
              fun.execute(
                spec.filename = sp_filename,
                region = region,
                extra = extra,
                verbose = verbose
              )
            }, identifier = paste0("execute ", sp_name))()

            catn("Parallel execute function completed successfully.")

            sink(type = "message")
            sink(type = "output")
            close(node_con)

            return(data = dt)
          })
        },
        error = function(e) {
          vebcat("An error occurred in the parallel process ~ stopping cluster and closing connections.", color = "fatalError")
          stop(e$message)
        }
      )

      if (batch) {
        return(list(
          init.res = init_res,
          exec.res = res
        ))
      } else {
        return(list(
          init.res = init_res,
          exec.res = res
        ))
      }
    },
    error = function(e) {
      vebcat("Error in visualizing parallel process.", color = "fatalError")
      stop(e$message)
    },
    finally = function() {
      if (exists("cl")) stopCluster(cl)
      if (exists("execute_timer")) end_timer(execute_timer)
      closeAllConnections()
    }
  )
}

parallel_spec_handler <- function(spec.dirs, dir, shape = NULL, extra = NULL, hv.project.method = "0.5-inclusion", col.n = NULL, out.order = NULL, fun, batch = FALSE, test = 0, node.log.append = TRUE, verbose = FALSE) {
  tryCatch(
    {
      if (verbose) {
        print_function_args()
      }
      input_functions <- fun()

      if (!exists("execute", where = input_functions)) {
        vebcat("Missing execute function to run parallel process.", color = "fatalError")
        stop()
      } else {
        vebcat("Found execute function.", veb = verbose)
        fun_execute <- input_functions$execute
      }

      if (!exists("init", where = input_functions)) {
        vebcat("Setting init function be the same as execute.")
        fun_init <- input_functions$execute
      } else {
        vebcat("Found init function.", veb = verbose)
        fun_init <- input_functions$init
      }

      if (!exists("process", where = input_functions)) {
        vebcat("Setting process function to use single function.")
        fun_process <- parallel_process_single
      } else {
        vebcat("Found process function.", veb = verbose)
        fun_process <- input_functions$process
      }

      vebcat("Acquiring", basename(dir), "values for all", hv.project.method, "rasters", color = "funInit")

      create_dir_if(dir)

      logs_sub_dir <- paste0(dir, "/", hv.project.method)
      vebprint(logs_sub_dir, verbose, "logs sub directory:")

      out_dir <- paste0(progressive_dirname(dir, end = 5), "/results")

      create_dir_if(logs_sub_dir, out_dir)

      out_file <- paste0(out_dir, "/", basename(progressive_dirname(dir, end = 7)), "/", basename(dir), ".csv")

      vebprint(out_file, verbose, "out file:")

      if (file.exists(out_file)) {
        post_timer <- start_timer("Post processing")
        catn("Found value table:", out_file)
        process_res <- fread(out_file)
      } else {
        catn("No previous table found, acquiring values for", hv.project.method, "method")

        combined_values <- data.table()

        parallel_res <- parallel_spec_dirs(
          spec.dirs = spec.dirs,
          dir = logs_sub_dir,
          shape = shape,
          extra = extra,
          hv.project.method = hv.project.method,
          fun.init = fun_init,
          fun.execute = fun_execute,
          batch = batch,
          test = test,
          node.log.append = node.log.append,
          verbose = verbose
        )

        catn("Length of init results:", highcat(1))
        catn("Length of execute results:", highcat(length(parallel_res$exec.res)))

        vebprint(parallel_res$init.res, verbose, "Parallel init result:")
        vebprint(head(parallel_res$exec.res, 1), verbose, "Parallel execute head(result, 1):")

        catn("Running post processing.")
        post_timer <- start_timer("Post processing")
        # Run the post process
        if (is.null(parallel_res$exec.res)) {
          process_res <- parallel_res$init.res
        } else {
          process_res <- fun_process(
            parallel.res = parallel_res,
            extra = extra,
            verbose = verbose
          )
        }

        vebprint(head(process_res, 3), verbose, "Processed data:")

        create_dir_if(dirname(out_file))

        fwrite(process_res, out_file, bom = TRUE)
      }

      if (!is.null(out.order)) {
        if (grepl("-", out.order)) {
          process_res <- process_res[order(-process_res[[gsub("-", "", out.order)]]), ]
        } else {
          process_res <- process_res[order(process_res[[out.order]]), ]
        }
      }

      if (is.null(col.n)) {
        catn("col.n is not used, returning the whole table")
      } else {
        # Syntax is: "column-number"
        split_str <- str_split(col.n, "-")[[1]]
        column <- split_str[[1]]
        number <- split_str[[2]]

        process_res <- unique(process_res, by = column)
        process_res <- process_res[1:number, ]
      }

      vebcat("Successfully acquired", basename(dir), "for all", hv.project.method, "rasters.", color = "funSuccess")

      return(process_res)
    },
    error = function(e) {
      vebcat("Error in parallel spec handler.", color = "fatalError")
      stop(e$message)
    },
    finally = function() {
      vebcat("Cleaning up parallel_spec_handler process")
      if (exists("post_timer")) end_timer(post_timer)
      if (exists("parallel_res")) rm(parallel_res)
      rm(list = ls(environment()))
      gc(full = TRUE)
    }
  )
}

parallel_init <- function(spec.filename, shape, extra = NULL, hv.project.method, fun.init, file.log, file.stop, verbose = FALSE) {
  tryCatch(
    {
      # Initiate memory control
      ram_control <- start_mem_tracking(file.log, file.stop)

      if (!is.null(shape)) {
        region <- load_region(shape)
        region <- handle_region(region)
      } else {
        region <- NULL
      }

      vebprint(class(fun.init), veb = verbose)
      vebprint(spec.filename, veb = verbose)
      vebprint(extra, veb = verbose)
      vebprint(region, veb = verbose)
      vebprint(verbose, veb = verbose)

      init_res <- fun.init(
        spec.filename = spec.filename,
        region = region,
        extra = extra,
        verbose = verbose
      )

      # Stop the memory tracker
      peak_mem_usage <- stop_mem_tracking(ram_control)

      return(init_res)
    },
    error = function(e) {
      vebcat("Error occurred:\n", e$message, color = "nonFatalError")
      vebcat("Stopping RAM check process.", color = "nonFatalError")
    },
    finally = function() {
      if (exists("ram_control")) stop_mem_tracking(ram_control)
      if (exists("region")) rm(region)
      rm(list = ls(environment()))
      gc(full = TRUE)
    }
  )
}

parallel_process_single <- function(parallel.res, extra, verbose) {
  tryCatch(
    {
      catn("Using the single processor.")
      merged_dt <- rbindlist(list(parallel.res$init.res, rbindlist(parallel.res$exec.res)), fill = TRUE)

      return(merged_dt)
    },
    error = function(e) {
      vebcat("Error in parallel process single.", color = "fatalError")
      stop(e$message)
    },
    finally = function() {
      if (exists("parallel.res")) rm(parallel.res)
      rm(list = ls(environment()))
      gc(full = TRUE)
    }
  )
}
