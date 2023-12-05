extract_name_after_prefix = function(x, prefixes = c("var", "subsp.")) {
  
  # Create a regular expression pattern to match the specified prefixes
  prefix_pattern = paste0(prefixes, collapse = "|")
  
  # Use str_extract to extract the name after the specified prefixes
  name = str_extract(x, paste0("(?<=", prefix_pattern, ")\\s*\\S+"))
  
  return(name)
} 