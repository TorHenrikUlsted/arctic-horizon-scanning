#------------------------#
####      Objects     ####
#------------------------#

get_obj_name <- function(...) {
  sapply(as.list(match.call())[-1], deparse)
}

get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
