check_syn_wfo_one <- function(wfo_checklist, folder) {
  vebcat("Starting the WFO.one synonym check", color = "funInit")
  
  wfo_one_checklist <- lapply(wfo_checklist, function(df) {
    WFO.one(WFO.result = df, priority = "Accepted", spec.name = "scientificName", verbose = F, counter = 100)
  })
  
  wfo_one_checklist <- rbindlist(wfo_one_checklist)
  
  wfo_one_nomatch <- data.frame(wfo_one_checklist[wfo_one_checklist$One.Reason == "no match found", ])
  
 if (nrow(wfo_one_nomatch) == 0) {
   vebcat("No missing matches found.", color = "proSuccess")
 } else {
   catn(
     colcat("Missing matches for", color = "nonFatalError"), 
     highcat(nrow(wfo_one_nomatch)), colcat("species, manually check these.", color = "nonFatalError")
    )
 }
  catn("Writing out to:", colcat(paste0(folder, "-wfo-one-nomatch.csv"), color = "output"))
  
  wfo_one_nomatch <- set_df_utf8(wfo_one_nomatch)

  fwrite(wfo_one_nomatch, paste0(folder, "/wfo-one-nomatch.csv"), bom = T)
  
  wfo_one_match <- data.frame(wfo_one_checklist[wfo_one_checklist$One.Reason != "no match found",])
  
  catn("Removed", highcat(nrow(wfo_one_checklist) - nrow(wfo_one_match)), "missing matches from WFO.one result.")
  
  wfo_one_dups <- wfo_one_match[duplicated(wfo_one_match$scientificName), ]
  
  catn("Found", highcat(nrow(wfo_one_dups)), "duplicated species from the WFO.one result.")
  
  catn("Writing out to:", colcat(paste0(folder, "/wfo-one-duplicates.csv"), color = "output"))
  
  wfo_one_dups <- set_df_utf8(wfo_one_dups)
  
  fwrite(wfo_one_dups, paste0(folder, "/wfo-one-duplicates.csv"), bom = T)
  
  wfo_one_uniq <- wfo_one_match[!duplicated(wfo_one_match$scientificName), ]
  
  catn("Removed", highcat(nrow(wfo_one_match) - nrow(wfo_one_uniq)), "species from the WFO.one result.")
  
  catn("Writing csv file out to", colcat(paste0(folder, "/wfo-one-uniq.csv"), color = "output"))

  wfo_one_uniq <- set_df_utf8(wfo_one_uniq)
  
  fwrite(wfo_one_uniq, paste0(folder, "/wfo-one-uniq.csv"), bom = T)
  
  vebcat("WFO_one results recieved", color = "funSuccess")

  return(list(
    wfo_one_uniq = wfo_one_uniq,
    wfo_one_nomatch = wfo_one_nomatch
  ))
}
