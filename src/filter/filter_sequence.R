source_all("./src/filter/components")
if (config$run$example) { 
  source_all("./example/src/filter")
  src_dir <- "./example/src"
} else {
  source_all("./src/filter/custom_filters")
  src_dir <- "./src"
}

filter_sequence <- function(spec.known = NULL, spec.unknown = NULL, test = NULL, approach = "precautionary", column = "scientificName", coord.uncertainty = NULL, cores.max = 1, region = NULL, download.key = NULL, download.doi = NULL, chunk.size = 1e6, chunk.iterations = NULL, force.seq = NULL, verbose = FALSE) {
  on.exit(closeAllConnections())
  vebcat("Initiating filtering sequence", color = "seqInit")

  #------------------------#
  ####       Setup      ####
  #------------------------#

  filter_timer <- start_timer("filter_timer")

  coord_un_file <- "./outputs/setup/region/coordinateUncertainty-m.txt"
  seq_fin_file <- "./outputs/filter/filter-sequence-completed.txt"

  if (!is.null(force.seq) && ("all" %in% force.seq | "filter" %in% force.seq)) {
    catn("Forcing filter sequence.")
    if (file.exists(seq_fin_file)) file.remove(seq_fin_file)
  }

  if (file.exists(seq_fin_file)) {
    catn("filter sequence already run, skipping process.")
    chunk_dir <- readLines(seq_fin_file)
  } else {
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
    
    if (!is.null(test)) {
      manual_edited <- "./resources/data-raw/test/manual-edit/wfo-nomatch-editeds.csv"
    } else {
      manual_edited <- "./resources/manual-edit/wfo-nomatch-editeds.csv"
    }
    
    dts <- fix_manual(
      dts = dts,
      manual.edited = ,
      column = column,
      verbose = verbose
    )
    
    
    if (is.null(coord.uncertainty) & file.exists(coord_un_file)) {
      coord.uncertainty <- as.numeric(readLines(coord_un_file))
    }
    
    vebprint(dts, verbose, "all data tables:")
    
    #------------------------#
    ####   Filter lists   ####
    #------------------------#
    # get original number of species
    if (!is.null(test)) {
      kpn <- nrow(dts[["test_known"]])
      kan <- 0
      spec_known <- get("filter_test_known")

      if (test == "small") {
        un <- nrow(dts[["test_small"]])
        spec_unknown <- get("filter_test_small")
        download.key <- "0180552-240321170329656"
        download.doi <- "https://doi.org/10.15468/dl.xzxdpx"
      } else if (test == "big") {
        un <- nrow(dts[["test_big"]])
        spec_unknown <- get("filter_test_big")
        download.key <- "0180885-240321170329656"
        download.doi <- "https://doi.org/10.15468/dl.sgf54g"
      } else {
        vebcat("Test has to be either 'small' or 'big'.", color = "fatalError")
        stop("Change the test parameter.")
      }
    } else {
      print(names(dts))
      known_names <- basename(list.files(paste0(src_dir, "/setup/custom_setup/wrangle/known_species")))
      known_names <- gsub("^wrangle_|\\.R$", "", known_names)
      pattern <- paste0("^(", paste(known_names, collapse="|"), ")")
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
      
      spec_known <- get("filter_known")
      spec_unknown <- get("filter_unknown")
    }
    
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
    kfan <- nrow(known$absent)
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
