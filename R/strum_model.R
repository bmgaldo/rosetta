#' Strum Model
#'
#' This function will create a strum model object.
#'
#' @param factor_structure A named list for the hypothesized factor structure. Set the factor names as the list names and measured variables as the list elements.
#' @param strum_object A strum object for the purposes of created a fixed covariance model.
#'
#' @return Returns a strum model object
#'
#' @importFrom strum createStrumModel
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' strum_model()

strum_model <- function(factor_structure, strum_object = NULL) {
  # measurement model
  measurement <- character(length(factor_structure)) # Initialize empty vector
  for(i in 1:length(factor_structure)) {
    measurement[i] <- paste(names(factor_structure)[i], "=~", paste(factor_structure[[i]], collapse = " + "), sep = " ")
  }
  # variance model
  variance <- paste("var(", names(factor_structure), ",e) = 1", sep = "")
  # covariance model
  if (length(names(factor_structure)) > 1) { # If we actually have covariances
    factor_combs <- combn(names(factor_structure), 2) # Choose-2 combinations for factors
    n_combs <- ncol(factor_combs)
    covariance <- character(n_combs) # Initialize empty vector
    for (i in 1:n_combs) {
      covariance[i] <- paste("cov(", paste(factor_combs[,i], collapse = ","), ",e) = NA", sep = "")
    }
  } else { # No covariances
    covariance <- NULL
  }

  # return model text
  return <- paste(
    paste(measurement, collapse = "\n"),
    paste(variance, collapse = "\n"),
    paste(covariance, collapse = "\n"),
    sep = "\n", collapse = ""
  )
  return
}
