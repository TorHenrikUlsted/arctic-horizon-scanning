check_syn_wfo <- function(checklist, column, folder) {
  if (!"data.table" %in% class(checklist) && !"data.frame" %in% class(checklist)) {
    stop("The input data is not in the 'data.table' or 'data.frame' format.", print(class(checklist)))
  }

  cat("Running the WFO synonym check with column:", cc$aquamarine(column), "for table: \n")
  print(head(checklist, 3))
  cat("Number of species to analyse: ", yellow(nrow(checklist)), "\n")
  cat(yellow("Expected waiting time in hours: ", round((((nrow(checklist) * 3.63) / 60) / 60), digits = 2), "hours \n"))
  cat(yellow("Expected waiting time in minutes: ", round(((nrow(checklist) * 3.63) / 60), digits = 2), "minutes \n"))

  wfo_timer <- start_timer("wfo_match")
  wfo_checklist <- WFO.match(spec.data = checklist, spec.name = column, WFO.file = WFO_file, verbose = T, counter = 1)

  wfo_checklist <- set_df_utf8(wfo_checklist)

  fwrite(wfo_checklist, paste0(folder, "wfo-match.csv"), bom = T)

  end_timer(wfo_timer)

  cat(cc$aquamarine("WFO synonym check completed \n"))

  return(wfo_checklist)
}
