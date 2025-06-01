source_all("./R/filter/components")
if (config$simulation$example) {
  source_all("./example/src/filter")
  src_dir <- "./example/src"
} else {
  source_all("./R/filter/custom_filters")
  src_dir <- "./R"
}

filter_sequence <- function(spec.known = NULL, spec.unknown, validation = FALSE, approach = "precautionary", column = "scientificName", coord.uncertainty = NULL, cores.max = 1, region = NULL, download.key = NULL, download.doi = NULL, chunk.size = 1e6, chunk.iterations = NULL, force.seq = NULL, verbose = FALSE) {
  #------------------------#
  ####       Setup      ####
  #------------------------#

  if (is.null(spec.unknown$name)) stop("spec.unknown cannot be NULL")

  filter_timer <- start_timer("filter_timer")

  filter_dir <- "./outputs/filter"
  wrangle_dir <- "./outputs/setup/wrangle"
  unknown_dir <- file.path(filter_dir, gsub("_", "-", spec.unknown$name))
  known_dir <- ifelse(!is.null(spec.known$name), file.path(filter_dir, gsub("_", "-", spec.known$name)), filter_dir)

  absent_final <- file.path(unknown_dir, "absent-final.csv")
  present_final <- file.path(known_dir, "present-final.csv")
  coord_un_file <- file.path(build_climate_path(), "coordinateUncertainty-m.txt")
  if (is.null(coord.uncertainty)) coord.uncertainty <- as.integer(readLines(coord_un_file))

  # Look for file in correct folder based on validation run
  seq_fin_file <- file.path(filter_dir, gsub("_", "-", if (!validation) {
    spec.unknown$name
  } else {
    spec.known$name
  }), "filter-sequence-completed.txt")

  if (!is.null(force.seq) && ("all" %in% force.seq | "filter" %in% force.seq | "validation" %in% force.seq)) {
    catn("Forcing filter sequence.")
    if (file.exists(seq_fin_file)) file.remove(seq_fin_file)
  }

  if (file.exists(seq_fin_file)) {
    catn("filter sequence already run, skipping process.")
    chunk_dir <- readLines(seq_fin_file)
    return(chunk_dir)
  }

  vebcat("Initiating filtering sequence", color = "seqInit")

  if (grepl("test", spec.unknown$name)) {
    manual_out <- file.path(wrangle_dir, "test/manual-check-file.csv")
    manual_edited <- "./resources/data-raw/test/manual-check-file.csv"
  } else {
    manual_out <- "./outputs/setup/wrangle/manual-check-file.csv"
    manual_edited <- "./resources/manual-edit/manual-check-file.csv"
  }

  if (file.exists(manual_out)) {
    man_l <- nrow(fread(manual_out))
    if (man_l > 0 & !file.exists(manual_edited)) {
      catn("Found", highcat(man_l), "species in need of manual handling.")
      catn("Find the file in:", colcat(manual_out, color = "output"))
      stop("Modify manually handled file and save as ", manual_edited)
    } else {
      catn("Manual edits found in:", colcat(manual_edited, color = "output"))
    }
  } else {
    catn("No manual handling file found. Continuing with the process.")
  }

  # Skip if the end file exists
  if (!file.exists(absent_final) || !file.exists(present_final)) {
    mdwrite(
      config$files$post_seq_md,
      text = "1;Filter Sequence"
    )

    vebcat("Loading dfs.", veb = verbose)

    dts <- select_wfo_column(
      dir.path = "./outputs/setup/wrangle",
      pattern = "wfo-one-gbif.csv",
      col.unique = column,
      col.select = NULL,
      verbose = verbose
    )

    if (length(manual_edited) > 0) {
      dts <- fix_manual(
        dts = dts,
        manual.edited = manual_edited,
        column = column,
        verbose = verbose
      )
    }

    vebprint(dts, verbose, "all data tables:")

    # Add sourceDataset
    dts <- setNames(
      lapply(seq_along(dts), function(i) {
        dt <- dts[[i]]
        dt_name <- names(dts)[i]
        dt[, sourceDataset := gsub("_.*", "", dt_name)]
        dt
      }), names(dts)
    )

    #------------------------#
    ####   Filter lists   ####
    #------------------------#
    # get original number of species
    if (grepl("test", spec.unknown$name)) {
      kpn <- nrow(dts[["test_known"]])
      kan <- 0

      if (grepl("small", spec.unknown$name)) {
        un <- nrow(dts[["test_small"]])
      } else if (grepl("big", spec.unknown$name)) {
        un <- nrow(dts[["test_big"]])
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
      un <- nrow(dts[[paste0(spec.unknown$name, "_absent")]])
    }

    spec_known <- get(paste0("filter_", spec.known$name))
    spec_unknown <- get(paste0("filter_", spec.unknown$name))

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
    ufn <- nrow(unknown)

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
  } else {
    known <- list(present = fread(present_final))
    unknown <- fread(absent_final)
  }

  #------------------------#
  #### Filter gbif keys ####
  #------------------------#

  keys_dts <- filter_gbif_keys(
    spec.dts = list(known = known$present, unknown = unknown),
    out.dirs = list(known = known_dir, unknown = unknown_dir),
    verbose = verbose
  )

  #------------------------#
  ####    Occurrence    ####
  #------------------------#

  if (grepl("test", spec.unknown$name)) {
    if (grepl("small", spec.unknown$name)) {
      spec.unknown$download.key <- "0180552-240321170329656"
      spec.unknown$download.doi <- "https://doi.org/10.15468/dl.xzxdpx"
    } else if (grepl("big", spec.unknown$name)) {
      spec.unknown$download.key <- "0038456-240906103802322"
      spec.unknown$download.doi <- "https://doi.org/10.15468/dl.85f3bs"
    } else {
      vebcat("Test has to be either 'test_small' or 'test_big'.", color = "fatalError")
      stop("Change the spec.unknown parameter.")
    }
  }

  if (!validation && (is.null(force.seq) || force.seq != "validation")) {
    occ_name <- file.path(unknown_dir, paste0(gsub("_", "-", spec.unknown$name), "-occ"))
    spec_w_keys <- keys_dts$unknown
    download.key <- spec.unknown$download.key
    download.doi <- spec.unknown$download.doi
    chunk_dir <- file.path(unknown_dir, "chunk")
  } else if (validation || force.seq == "validation") {
    occ_name <- file.path(known_dir, paste0(gsub("_", "-", spec.known$name), "-occ"))
    spec_w_keys <- keys_dts$known
    download.key <- spec.known$download.key
    download.doi <- spec.known$download.doi
    chunk_dir <- file.path(known_dir, "chunk")
  } else {
    vebcat("ERROR: something went wrong when preparing before occurence download.")
    vebprint(validation, text = "validation:")
    vebprint(force.seq, text = "force.seq:")
    stop("Something went")
  }

  occ_data <- get_occ_data(
    species_w_keys = spec_w_keys,
    file.name = occ_name,
    region = region,
    coord.uncertainty = coord.uncertainty,
    download.key = download.key,
    download.doi = download.doi,
    verbose = verbose
  )

  #------------------------#
  ####     Chunking     ####
  #------------------------#

  chunk_protocol(
    spec.occ = occ_data,
    spec.keys = spec_w_keys,
    chunk.name = "species",
    chunk.col = "cleanName",
    chunk.dir = chunk_dir,
    chunk.size = chunk.size,
    iterations = chunk.iterations,
    approach = approach,
    verbose = verbose
  )

  return_dir <- paste0(chunk_dir, "/species")

  create_file_if(seq_fin_file)
  writeLines(
    return_dir,
    con = seq_fin_file
  )

  mdwrite(
    config$files$post_seq_md,
    text = paste0("Number of species in ", basename(dirname(chunk_dir)), " returned from chunking **", length(list.files(return_dir)), "**")
  )

  end_timer(filter_timer)

  vebcat("Filtering sequence completed successfully.", color = "seqSuccess")

  return(return_dir)
}
