#' Strum Data
#'
#' This function will create a strum data object and fix twin kinships.
#'
#' @param data A list of dataframes that contain a varying number of measured variables in respect to the underlying factor structure.
#' @param type The type of input data: "Pedigree" or "RawData".
#' @param twin The name of the twin zygotic variable as a character string.
#'
#' @return Returns a strum data object
#'
#' @importFrom strum createStrumData
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' strum_data()

strum_data <- function(data, type = "Pedigree", twin = NULL) {
  # Create strum data
  sdata <- createStrumData(inData = data, dType = type)

  # Fix kinship matrix for twin types
  if(!(is.null(twin))) {
    strum_kin <- sdata@phi
    data_family <- split(
      data[, c("family", "id", "father", "mother", twin)],
      data$family
    )
    mz_location <- lapply(data_family, function(z) {z[z[twin] == 1 & !is.na(z[twin]), "id"]}) # Need to remove NAs or else they will be included in mz_location
    strum_kin_names <- names(strum_kin)
    n_strum_kin <- length(strum_kin)
    for (i in 1:n_strum_kin) {
      if (length(mz_location[[i]]) > 1) {
        strum_kin[[strum_kin_names[i]]][mz_location[[strum_kin_names[i]]][1],mz_location[[strum_kin_names[i]]][2]] <- 1
        strum_kin[[strum_kin_names[i]]][mz_location[[strum_kin_names[i]]][2],mz_location[[strum_kin_names[i]]][1]] <- 1
      }
    }
    sdata@phi <- strum_kin
  }

  # Return
  sdata
}
