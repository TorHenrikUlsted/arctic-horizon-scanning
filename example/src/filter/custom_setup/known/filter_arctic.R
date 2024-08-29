filter_arctic <- function(dts, column, verbose = FALSE) {
  ##############
  # Initialize
  ##############

  mdwrite(
    config$files$post_seq_md,
    text = "3;Arctic"
  )

  aba_absent <- dts$aba_absent
  ambio_absent <- dts$ambio_absent

  aba_present <- dts$aba_present
  ambio_present <- dts$ambio_present

  # Union_dfs merges and removes duplicates while also provide info on how many are removed
  arctic_present <- union_dts(aba_present, ambio_present, verbose = T)

  write_filter_fun(
    file.out = "./outputs/filter/arctic/arctic-present-final.csv",
    spec.in = arctic_present
  )

  # -------------------------------------------------------------------------- #

  arctic_absent <- union_dts(aba_absent, ambio_absent, verbose = T)

  arctic_absent <- write_filter_fun(
    file.out = "./outputs/filter/arctic/arctic-absent-final.csv",
    spec.in = arctic_absent,
    fun = function() {
      # Also remove all arctic_present from absent in case some standard names have changed
      print(column)
      ab <- anti_union(arctic_absent, arctic_present, column)

      return(ab)
    }
  )

  return(list(
    present = arctic_present,
    absent = arctic_absent
  ))
}
