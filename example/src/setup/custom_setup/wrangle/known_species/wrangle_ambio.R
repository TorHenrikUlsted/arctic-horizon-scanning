wrangle_ambio <- function(name, column, verbose = FALSE) {
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)

  formatted_out <- paste0(dir, "/", name, "-formatted.csv")
  present_out <- paste0(dir, "/", name, "-present.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  
 if (file.exists(formatted_out) && file.exists(present_out) && file.exists(absent_out)) {
   present <- fread(present_out)
   absent <- fread(absent_out)
 } else {
   vebcat("Initiating", name, "wrangling protocol", color = "funInit")
   preformatted <- fread(paste0("./resources/data-raw/", name, ".csv"), header = F)
   
   vebcat("filtering rows.", veb = verbose)
   ## remove the not important rows and columns
   formatted <- preformatted[-c(1, 3:9), ]
   formatted <- formatted[, -25]
   
   ## set the first row to be the header and remove the row from the dataset
   colnames(formatted) <- as.character(unlist(formatted[1, ]))
   formatted <- formatted[-1, ]
   ## give the first column a name
   colnames(formatted)[1] <- "verbatimName"
   formatted <- unique(formatted, by = "verbatimName")
   formatted[, interimName := verbatimName]
   
   formatted[, interimName := { # Clean symbols & designations
     tmp <- clean_string(interimName, verbose)
     tmp <- clean_designations(tmp, config$species$standard_infraEpithets, verbose)
     clean_symbols(tmp, config$species$standard_symbols, verbose)
   }]
   
   formatted[, c("genus", "specificEpithet", "cleanName", "interimAuthorship", "fullName", "structure") := {
     res <- lapply(seq_along(interimName), function(i) {
       cat("\rProcessing species:", i, "of", .N)
       flush.console()
       clean_spec_name(interimName[i], config$species$standard_symbols, config$species$standard_infraEpithets, verbose)
     });catn()
     list(
       vapply(res, function(x) x$genus, character(1)),
       vapply(res, function(x) x$specificEpithet, character(1)),
       vapply(res, function(x) x$cleanName, character(1)),
       vapply(res, function(x) x$other, character(1)),
       vapply(res, function(x) x$fullName, character(1)),
       vapply(res, function(x) x$structure, character(1))
     )
   }]
   
   vebcat("Creating present df.", veb = verbose)
   # Create present and absent lists
   ## create conditions
   ### symbols meaning: Present: ●, IR, IT and Absent: ○, ?, ×
   condition1 <- apply(formatted[, 2:24], 1, function(x) any(x %in% c("●", "IR", "IT")))
   ## Use the condition to create present and absent lists
   present <- subset(formatted, subset = condition1)
   present <- present[, .(verbatimName, cleanName, interimAuthorship)]
   
   setnames(present, old = "cleanName", new = column)
   setnames(present, old = "interimAuthorship", new = paste0(column,"Authorship"))
   
   vebcat("Creating absent df.", veb = verbose)
   ## Only outputs unquie species names
   absent <- subset(formatted, subset = !condition1)
   absent <- absent[, .(verbatimName, cleanName, interimAuthorship)]
   
   setnames(absent, old = "cleanName", new = column)
   setnames(absent, old = "interimAuthorship", new = paste0(column,"Authorship"))
   
   catn("Writing out files to:", colcat(dir, color = "output"))
   
   fwrite(formatted, formatted_out, bom = T)
   
   fwrite(present, present_out, bom = T)
   
   fwrite(absent, absent_out, bom = T)
   
   fl <- nrow(formatted)
   pl <- nrow(present)
   al <- nrow(absent)
   
   md_dt <- data.table(
     formatted = fl,
     present = pl,
     absent = al,
     lost = fl - sum(pl, al)
   )
   
   mdwrite(
     config$files$post_seq_md,
     text = "3;Wrangled AMBIO",
     data = md_dt
   )
   vebcat(name, "wrangling protocol successfully completed.", color = "funSuccess")
 }
  
  present[, sourceDataset := name]
  absent[, sourceDataset := name]
  
  return(list(
    present = present,
    absent = absent
  ))
}
