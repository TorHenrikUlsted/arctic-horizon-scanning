cite_packages <- function(pkgs, formats = c("bibtex")) {
  lapply(formats, function(format) {
    all_citations <- lapply(pkgs, function(pkg) {
      tryCatch(
        {
          cit <- citation(pkg)
          bib <- NULL
          if (format == "bibtex") {
            bib <- toBibtex(cit)
          }
          
          # Add more formats here in the future like this:
          # if (format == "RIS") {
          #   bib <- toRIS(cit)
          # }
          
          # Convert bibentry to character
          bib <- as.character(bib)
          return(bib)
        },
        error = function(e) {
          stop(conditionMessage(e))
        }
      )
    })
    # Filter out NA values
    all_citations <- all_citations[!is.na(all_citations)]
    
    if (!dir.exists("./outputs/utils/references/")) dir.create("./outputs/utils/references/", recursive = T)
    
    # Write to file
    file_name <- paste0("./outputs/utils/references/citations.", format)
    
    writeLines(unlist(all_citations), file_name)
  })
}
