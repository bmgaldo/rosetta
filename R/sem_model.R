#' sem_model
#'
#' Create a RAM specified Structural Equation Model for the 'sem' package.
#'
#' @param factor_structure A named list for the hypothesized factor structure. Set the factor names as the list names and measured variables as the list elements.
#' @param sem_object An sem object for the purposes of creating an sem model with fixed covariance values based on sem_object.
#'
#' @return Returns an object of class semmod.
#'
#' @importFrom sem specifyModel
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' sem_model()

sem_model <- function(factor_structure, sem_object = NULL) {
  # Character values and inputs for SEM model
  vars <- stack(factor_structure)[["values"]]
  factors <- as.character(stack(factor_structure)[["ind"]])
  unique_factors <- unique(factors)
  n_factors <- length(factor_structure)
  n_vars <- length(vars)
  coefs <- paste("lam_", vars, sep = "")
  var_params <- paste("e_", vars, sep = "")
  fac_diag_values <- rep(1, n_factors)
  if (n_factors > 1) {
    factor_combs <- combn(unique_factors, 2) # Choose-2 combinations for factors
  } else {
    factor_combs <- matrix(unique_factors)
  }
  n_combs <- ncol(factor_combs)
  fac_cor_params <- character(n_combs) # Initialize empty vector
  for (i in 1:n_combs) {
    fac_cor_params[i] <- paste("f_", paste(factor_combs[,i], collapse = "_"), sep = "")
  }
  if(is.null(sem_object)) { # Set to NA for unconstrained model
    fac_cor_values <- rep(NA, n_combs)
  } else { # Set to estimated values for constrained model
    if (n_factors > 1) {
      fac_cor_values <- fac_cor_estimates(sem_object) # Store dataframe of factor covariance estimates
      fac_cor_values <- fac_cor_values[match(fac_cor_params, fac_cor_values$covariance), ] # Make sure the covariance values are matched/ordered to the correct covariance
      fac_cor_values <- na.omit(fac_cor_values) # Remove any factor covariances that are not need/matched
      fac_cor_values <- fac_cor_values$estimate
      fac_cor_params <- rep(NA, n_combs)
    }
  }

  # Create a tidy dataset.
  ## Columns
  if (n_factors > 1) {
    left <- c(factors, vars, unique_factors, factor_combs[1,])
    right <- c(vars, vars, unique_factors, factor_combs[2,])
    arrow <- c(rep("->", n_vars), rep("<->", n_vars), rep("<->", n_factors), rep("<->", n_combs))
    param <- c(coefs, var_params, rep(NA, n_factors), fac_cor_params)
    value <- c(rep(NA, n_vars), rep(NA, n_vars), fac_diag_values, fac_cor_values)
  } else {
    left <- c(factors, vars, unique_factors)
    right <- c(vars, vars, unique_factors)
    arrow <- c(rep("->", n_vars), rep("<->", n_vars), rep("<->", n_factors))
    param <- c(coefs, var_params, rep(NA, n_factors))
    value <- c(rep(NA, n_vars), rep(NA, n_vars), fac_diag_values)
  }

  ## Dataframe
  data <- data.frame(left = left, arrow = arrow, right = right, param = param,
                     value = value, stringsAsFactors = FALSE)

  # Create Empty character vector that will be filled. This avoids "writing" to disk.
  if (n_factors > 1) {
    text <- character(n_vars*2 + n_factors + n_combs)
  } else {
    text <- character(n_vars*2 + n_factors)
  }
  # Fill in character vector
  for (i in 1:length(text)) { # Specify regression coefficients
    text[i] <- paste(
      paste(as.character(data[i,1:3]), collapse = " "),
      paste(as.character(data[i,4:5]), collapse = ", "),
      sep = ", "
    )
  }

  # Return SEM model
  model_text <- paste(text, collapse = "\n")
  specified_model <- specifyModel(text = model_text, quiet = TRUE)
  specified_model
}
