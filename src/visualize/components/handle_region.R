handle_region <- function(region) {
  region_desc <- fread("./resources/region/cavm2003-desc.csv")
  
  index <- match(region$FLOREG, region_desc$FLOREG)
  
  region$country <- region_desc$country[index]
  region$floristicProvince <- region_desc$floristicProvince[index]
  
  # Remove the ice sheet
  region <- region[region$FLOREG != 0, ]
  region <- na.omit(region)
  
  return(region)
}