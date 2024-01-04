get_gbif_names_keys <- function(region, region.name, sp_keys, cores.max = 1, verbose = F) {
  cat(blue("Initiating GBIF name by key acquisition protocol. \n"))

  cl <- makeCluster(cores.max)

  clusterEvalQ(cl, {
    library(rgbif)
  })

  gbif_name_timer <- start_timer("gbif_names_timer")

  gbif_names <- clusterApplyLB(cl, seq_along(sp_keys), function(i) {
    tryCatch(
      {
        progress_msg <- paste(
          "collecting species\n",
          sprintf("%6s | %10s | %8s", "key", "Max keys", "Percent"),
          sprintf("\r%6d | %10d | %8.2f %%", i, length(sp_keys), round(i / length(sp_keys) * 100, 2)),
          collapse = "\n"
        )

        writeLines(progress_msg, con = "./outputs/filter/logs/sp-names-progress.txt")

        scientificName <- occ_count(
          speciesKey = sp_keys[i],
          facet = "scientificName",
          hasCoordinate = T
        )$scientificName

        list(scientificName = scientificName, speciesKey = sp_keys[i])
      },
      error = function(e) {
        e <- conditionMessage(e)
        print(e)
      }
    )
  })


  cat("Finishing up \n")

  stopCluster(cl)

  end_timer(gbif_name_timer)

  cat(cc$lightSteelBlue(length(gbif_names)), "GBIF names found. \n")

  ## Make the species names into strings
  gbif_names <- as.character(do.call(c, gbif_names))
  gbif_names <- data.frame(gbif_names = trimws(gbif_names))

  names_savefile <- paste0("./outputs/filter/gbif-acquisition/gbif-names-", region.name, ".csv")
  cat(yellow("Writing out gbif_names to: \n", names_savefile, "\n"))
  gbif_names <- set_df_utf8(gbif_names)
  fwrite(gbif_names, names_savefile, row.names = F, bom = T)

  cat(cc$lightGreen("GBIF name by key acquisition protocol completed successfully. \n"))

  return(gbif_names)
}
