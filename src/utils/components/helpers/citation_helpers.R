#------------------------#
####    Citations     ####
#------------------------#

format_bibtex <- function(bibtex_string) {
  # Remove leading/trailing whitespace and newline characters
  bibtex_string <- trimws(bibtex_string)
  bibtex_string <- gsub("\n", "", bibtex_string)
  
  # Split the string into individual fields
  fields <- strsplit(bibtex_string, ",(?=[^}]*(?:\\{|$))", perl = TRUE)[[1]]
  
  # Extract the entry type and key
  entry_type <- sub("^@(\\w+)\\{.*", "\\1", fields[1])
  
  # Process each field
  formatted_fields <- c(paste0("@", entry_type, "{,"))
  for (i in 2:length(fields)) {
    field <- fields[i]
    # Extract field name and value
    parts <- strsplit(field, "=")[[1]]
    field_name <- trimws(parts[1])
    field_value <- trimws(paste(parts[-1], collapse = "="))
    
    # Remove surrounding braces or quotes from field value
    field_value <- gsub("^[{\"](.*)[}\"]$", "\\1", field_value)
    
    # Format the field
    if (i == length(fields)) {
      # Last field: remove trailing } if present
      field_value <- sub("\\s*}\\s*$", "", field_value)
      formatted_fields <- c(formatted_fields, paste0("  ", field_name, " = {", field_value, "}"))
    } else {
      formatted_fields <- c(formatted_fields, paste0("  ", field_name, " = {", field_value, "},"))
    }
  }
  
  # Add the closing brace as a separate element
  formatted_fields <- c(formatted_fields, "}")
  
  # Return as a character vector
  return(formatted_fields)
}