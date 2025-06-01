check_syn_wfo <- function(checklist, cols, out.dir, cores.max = 1, verbose = FALSE, counter = 1) {
  if (!"data.table" %in% class(checklist) && !"data.frame" %in% class(checklist)) {
    stop("The input data is not in the 'data.table' or 'data.frame' format.", print(class(checklist)))
  }
  wfo_timer <- start_timer("wfo_match")

  checklist <- as.data.table(checklist)
  out_file <- paste0(out.dir, "/wfo-match.csv")

  if (!file.exists(out_file)) {
    vebcat("Initiating WFO synonym check", color = "funInit")


    input_cols <- unlist(cols)

    # Use the provided columns or default values if not present
    cols$spec.name <- ifelse(!is.null(cols$spec.name), cols$spec.name, "spec.name")
    cols$Genus <- ifelse(!is.null(cols$Genus), cols$Genus, "Genus")
    cols$Species <- ifelse(!is.null(cols$Species), cols$Species, "Species")
    cols$Infraspecific <- ifelse(!is.null(cols$Infraspecific), cols$Infraspecific, "Infraspecific")
    cols$Infraspecific.rank <- ifelse(!is.null(cols$Infraspecific.rank), cols$Infraspecific.rank, "Infraspecific.rank")
    cols$Authorship <- ifelse(!is.null(cols$Authorship), cols$Authorship, "Authorship")

    vebcat("Cols:", cols, veb = verbose)

    catn(
      "Running the WFO synonym check for",
      highcat(nrow(checklist)),
      "species with column(s):\n",
      highcat(paste(input_cols, sep = " = ", collapse = ", "), "\n")
    )

    temp_dt <- checklist[, ..input_cols]
    vebprint(temp_dt[1:3], text = "Table sample:")
    rm(temp_dt)

    if (nrow(checklist) < 10) {
      wfo_result <- WFO.match(
        spec.data = checklist,
        WFO.file = WFO_file,
        spec.name = cols$spec.name,
        Genus = cols$Genus,
        Species = cols$Species,
        Infraspecific.rank = cols$Infraspecific.rank,
        Infraspecific = cols$Infraspecific,
        Authorship = cols$Authorship,
        verbose = verbose,
        counter = counter
      )
    } else {
      custom_evals <- list(
        packages = c(
          "WorldFlora",
          "data.table"
        ),
        source = c(
          "./R/utils/components/handlers/condition_handlers.R",
          "./R/utils/components/handlers/time_handlers.R",
          "./R/utils/components/handlers/file_handlers.R",
          "./R/utils/components/loader.R"
        )
      )

      cores_max <- calc_num_cores(
        ram.high = 2,
        cores.total = config$memory$total_cores,
        verbose = verbose
      )

      wfo_result <- wfo_parallel(
        checklist = checklist,
        cols = cols,
        out.dir = out.dir,
        cores.max = min(nrow(checklist), cores_max$total),
        evals = custom_evals,
        counter = counter,
        verbose = verbose
      )
    }


    if (!is.data.table(wfo_result)) wfo_result <- as.data.table(wfo_result)
    fwrite(wfo_result, paste0(out.dir, "/wfo-match.csv"), bom = T)
  } else {
    catn("Reading existing WFO.match file from:", colcat(out_file, color = "output"))
    wfo_result <- fread(out_file)
  }

  wfo_result <- wfo_mismatch_check(
    wfo.result = wfo_result,
    col.origin = cols[["spec.name"]],
    out.file = paste0(out.dir, "/wfo-match-mismatches.csv"),
    unchecked = FALSE,
    verbose = verbose
  ) # returns list of clean and mismatched info

  end_timer(wfo_timer)

  vebcat("WFO.match synonym check completed", color = "funSuccess")

  return(wfo_result)
}
