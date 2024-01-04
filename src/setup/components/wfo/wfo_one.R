check_syn_wfo_one <- function(wfo_checklist, folder) {
  cat(cc$aquamarine("Starting the WFO.one synonym check \n"))
  
  wfo_one_checklist <- lapply(wfo_checklist, function(df) {
    WFO.one(WFO.result = df, priority = "Accepted", spec.name = "scientificName", verbose = F, counter = 100)
  })
  
  wfo_one_checklist <- rbindlist(wfo_one_checklist)
  
  wfo_one_nomatch <- data.frame(wfo_one_checklist[wfo_one_checklist$One.Reason == "no match found", ])
  
 if (nrow(wfo_one_nomatch) == 0) {
   cat(cc$lightGreen("No missing matches found. \n"))
 } else {
   cat(cc$lightCoral("Missing matches for"), cc$lightSteelBlue(nrow(wfo_one_nomatch)), cc$lightCoral("species, manually check these."), "\n")
 }
  cat("Writing out to:", yellow(paste0(folder, "-wfo-one-nomatch.csv")), "\n")
  
  wfo_one_nomatch <- set_df_utf8(wfo_one_nomatch)

  fwrite(wfo_one_nomatch, paste0(folder, "/wfo-one-nomatch.csv"), bom = T)
  
  wfo_one_match <- data.frame(wfo_one_checklist[wfo_one_checklist$One.Reason != "no match found",])
  
  cat("Removed", cc$lightSteelBlue(nrow(wfo_one_checklist) - nrow(wfo_one_match)), "missing matches from WFO.one result. \n")
  
  wfo_one_dups <- wfo_one_match[duplicated(wfo_one_match$scientificName), ]
  
  cat("Found", cc$lightSteelBlue(nrow(wfo_one_dups)), "duplicated species from the WFO.one result. Writing out to file. \n")
  
  wfo_one_dups <- set_df_utf8(wfo_one_dups)
  
  fwrite(wfo_one_dups, paste0(folder, "/wfo-one-duplicates.csv"), bom = T)
  
  wfo_one_uniq <- wfo_one_match[!duplicated(wfo_one_match$scientificName), ]
  
  cat("Removed", cc$lightSteelBlue(nrow(wfo_one_match) - nrow(wfo_one_uniq)), "species from the WFO.one result. \n")
  
  cat("Writing csv file out to", yellow(paste0(folder, "/wfo-one-uniq.csv")), "\n")

  wfo_one_uniq <- set_df_utf8(wfo_one_uniq)
  
  fwrite(wfo_one_uniq, paste0(folder, "/wfo-one-uniq.csv"), bom = T)
  
  cat(cc$lightGreen("WFO_one results recieved \n"))

  return(list(
    wfo_one_uniq = wfo_one_uniq,
    wfo_one_nomatch = wfo_one_nomatch
  ))
}
