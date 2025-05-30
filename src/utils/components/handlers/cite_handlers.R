cite_packages <- function(pkgs, formats = c("bibtex")) {
  citations <- setNames(lapply(formats, function(format) {
    all_citations <- setNames(lapply(pkgs, function(pkg) {
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
          warning(paste("Error citing package", pkg, ":", conditionMessage(e)))
          return(NULL)
        }
      )
    }), pkgs)
    
    # Remove NULL entries
    all_citations <- all_citations[!sapply(all_citations, is.null)]
    
    return(all_citations)
  }), formats)
  
  return(citations)
}

write_citations <- function(citations, out.dir) {
  # Filter out NA values
  citations <- citations[!is.na(citations)]
  
  if (!dir.exists(out.dir)) dir.create(out.dir, recursive = T)
  
  # Write to file
  for (format in names(citations)) {
    file_name <- paste0(out.dir, "/citations", ".", format)
    
    citation_texts <- sapply(citations, function(x) paste(x$citation, collapse = "\n"))
    
    writeLines(unlist(citations), file_name)
  }
  
  catn("Citations written to directory:", colcat(out.dir, color = "output"))
}
