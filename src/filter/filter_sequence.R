source_all("./src/filter/components")
if (config$simulation$example) {
  source_all("./example/src/filter")
  src_dir <- "./example/src"
} else {
  source_all("./src/filter/custom_filters")
  src_dir <- "./src"
}

filter_sequence <- function(spec.known = NULL, spec.unknown, approach = "precautionary", column = "scientificName", coord.uncertainty = NULL, cores.max = 1, region = NULL, download.key = NULL, download.doi = NULL, chunk.size = 1e6, chunk.iterations = NULL, force.seq = NULL, verbose = FALSE) {
  vebcat("Initiating filtering sequence", color = "seqInit")

  #------------------------#
  ####       Setup      ####
  #------------------------#

  if (is.null(spec.unknown)) stop("spec.unknown cannot be NULL")

  filter_timer <- start_timer("filter_timer")

  coord_un_file <- paste0(build_climate_path(), "/coordinateUncertainty-m.txt")
  seq_fin_file <- paste0("./outputs/filter/", gsub("_", "-", spec.unknown), "/filter-sequence-completed.txt")
  

  if (!is.null(force.seq) && ("all" %in% force.seq | "filter" %in% force.seq)) {
    catn("Forcing filter sequence.")
    if (file.exists(seq_fin_file)) file.remove(seq_fin_file)
  }

  if (file.exists(seq_fin_file)) {
    catn("filter sequence already run, skipping process.")
    chunk_dir <- readLines(seq_fin_file)
  } else {
    if (grepl("test", spec.unknown)) {
      manual_out <- "./outputs/setup/wrangle/test/manual-check-file.csv"
      manual_edited <- "./resources/data-raw/test/manual-check-file.csv"
    } else {
      manual_out <- "./outputs/setup/wrangle/manual-check-file.csv"
      manual_edited <- list.files("./resources/manual-edit", pattern = ".csv", full.names = TRUE)
    }
    
    man_l <- nrow(fread(manual_out))
    
    if (length(manual_edited) == 0 && (!file.exists(manual_out) || man_l > 0)) {
      vebcat("Need manually handled file to continue", color = "fatalError")
      if (file.exists(manual_out)) {
        catn("Found", highcat(man_l), "species in need of manual handling.")
        catn("Find the file in:", colcat(manual_out, color = "output"))
      } else {
        catn("Could not find manually handled output file, check setup process and input data.")
      }
      stop("Modify manually handled file")
    }
    

    mdwrite(
      config$files$post_seq_md,
      text = "1;Filter Sequence"
    )

    vebcat("Loading dfs.", veb = verbose)

    dts <- select_wfo_column(
      dir.path = "./outputs/setup/wrangle", # change this to outputs path also edit in wfo
      pattern = "wfo-one-clean.csv",
      col.unique = column,
      col.select = NULL,
      verbose = verbose
    )
    
    print(length(manual_edited))
    if (length(manual_edited) > 0) {
      dts <- fix_manual(
        dts = dts,
        manual.edited = manual_edited,
        column = column,
        verbose = verbose
      )
    }
    
    # Change symbols to config need to filename safe
    lapply(dts, function(dt) {
      dt[, `:=` (
        scientificName = {
          tmp <- clean_symbols(scientificName, config$species$standard_symbols, verbose)
          clean_designations(tmp, config$species$standard_infraEpithets, verbose)
        } 
      )]
    })
    
    if (is.null(coord.uncertainty) & file.exists(coord_un_file)) {
      coord.uncertainty <- as.numeric(readLines(coord_un_file))
    }

    vebprint(dts, verbose, "all data tables:")

    #------------------------#
    ####   Filter lists   ####
    #------------------------#
    # get original number of species
    if (grepl("test", spec.unknown)) {
      kpn <- nrow(dts[["test_known"]])
      kan <- 0

      if (grepl("small", spec.unknown)) {
        un <- nrow(dts[["test_small"]])
        download.key <- "0180552-240321170329656"
        download.doi <- "https://doi.org/10.15468/dl.xzxdpx"
      } else if (grepl("big", spec.unknown)) {
        un <- nrow(dts[["test_big"]])
        download.key <- "0180885-240321170329656"
        download.doi <- "https://doi.org/10.15468/dl.sgf54g"
      } else {
        vebcat("Test has to be either 'test_small' or 'test_big'.", color = "fatalError")
        stop("Change the spec.unknown parameter.")
      }
    } else {
      known_names <- basename(list.files(paste0(src_dir, "/setup/custom_setup/wrangle/known_species")))
      known_names <- gsub("^wrangle_|\\.R$", "", known_names)
      pattern <- paste0("^(", paste(known_names, collapse = "|"), ")")
      known_dts <- dts[grep(pattern, names(dts))]

      kn <- list(present = 0, absent = 0)
      for (dt_name in names(known_dts)) {
        dt <- known_dts[[dt_name]]
        n <- nrow(dt)
        if (grepl("present", dt_name)) {
          kn$present <- kn$present + n
        } else {
          kn$absent <- kn$absent + n
        }
      }
      kpn <- sum(unlist(kn$present))
      kan <- sum(unlist(kn$absent))
      un <- nrow(dts[[paste0(spec.unknown, "_absent")]])
    }
    
    spec_known <- get(paste0("filter_", spec.known))
    spec_unknown <- get(paste0("filter_", spec.unknown))

    if (!is.null(spec_known)) {
      known <- spec_known(
        dts = dts,
        column = column,
        verbose = verbose
      )
    } else {
      known <- NULL
      warning("known returns NULL, ignore if intentional, else press ESC and fix spec.known parameter")
      catn("Waiting 5 seconds before continuing...")
      Sys.sleep(5)
    }
    vebprint(known, text = "known output:", veb = verbose)

    if (is.null(spec_unknown)) {
      vebcat("Error, spec_unknown is NULL.", color = "fatalError")
      vebprint(spec_unknown, text = "spec_unknown:")
      return(stop("Check spec.unknown parameter."))
    } else {
      unknown <- spec_unknown(
        known.filtered = known,
        dts = dts,
        column = column,
        verbose = verbose
      )
    }

    vebprint(known, text = "known output:", veb = verbose)

    kfpn <- nrow(known$present)
    if (is.null(known$absent)) {
      kfan <- 0
    } else {
      (
        kfan <- nrow(known$absent)
      )
    }
    ufn <- nrow(unknown$spec)

    md_dt <- data.table(
      present.in = kpn,
      present.out = kfpn,
      absent.in = kan,
      absent.out = kfan,
      unknown.in = un,
      unknown.out = ufn
    )

    mdwrite(
      config$files$post_seq_md,
      text = "2;Filter function summary",
      data = md_dt
    )

    #------------------------#
    ####    Occurrence    ####
    #------------------------#

    occ_name <- paste0(unknown$dir, "/", basename(unknown$dir), "-occ")

    unknown_occ <- get_occurrence(
      spec = unknown$spec,
      file.out = occ_name,
      region = region,
      coord.uncertainty = coord.uncertainty,
      download.key = download.key,
      download.doi = download.doi,
      verbose = verbose
    )

    if (!is.null(unknown_occ$occ)) {
      sp_occ <- unknown_occ$occ
    } else {
      sp_occ <- paste0(occ_name, ".csv")
    }

    #------------------------#
    ####     Chunking     ####
    #------------------------#

    chunk_dir <- paste0(unknown$dir, "/chunk")

    chunk_protocol(
      spec.occ = sp_occ,
      spec.keys = unknown_occ$keys,
      chunk.name = "species",
      chunk.col = "cleanName",
      chunk.dir = chunk_dir,
      chunk.size = chunk.size,
      cores.max = 1,
      iterations = chunk.iterations,
      approach = approach,
      verbose = verbose
    )

    create_file_if(seq_fin_file)
    writeLines(
      chunk_dir,
      con = seq_fin_file
    )
  }

  mdwrite(
    config$files$post_seq_md,
    text = paste0("3;Number of species returned from chunking **", length(list.files(paste0(chunk_dir, "/species"))), "**")
  )

  end_timer(filter_timer)

  vebcat("Filtering sequence completed successfully.", color = "seqSuccess")

  return(paste0(chunk_dir, "/species"))
}
