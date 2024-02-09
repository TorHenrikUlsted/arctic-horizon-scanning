create_file_if <- function(file_path, keep = F) {
  
  if (!file.exists(file_path)) {
    
    file.create(file_path)
    
  } else {
    
    if (keep == FALSE) {
      
      file.remove(file_path)
      
      file.create(file_path)
      
    }
  }
}