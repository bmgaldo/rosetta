# Check if correlation or covariance matrices are positive semi definite
is_pos_semidef <- function(matrix) {
  all(eigen(matrix, only.values = TRUE)$values >= 0)
}

# Converts all factor type variables in dataframe to character type.
fac2char <- function(data) {
  factor_index <- sapply(data, is.factor)
  data[factor_index] <- lapply(data[factor_index], as.character)
  data
}

# Returns factor correlation estimates from an sem object.
#' @importFrom sem stdCoef
fac_cor_estimates <- function(sem_object) {
  n_factors <- sem_object$m - sem_object$n
  coef <- fac2char(stdCoef(sem_object))
  n_coef <- nrow(stdCoef(sem_object))
  n_combs <- ncol(combn(n_factors, 2))
  fac_coef <- coef[(n_coef-n_combs+1):n_coef, 1:2]
  colnames(fac_coef) <- c("covariance", "estimate")
  rownames(fac_coef) <- NULL
  fac_coef
}

# Function to estimate missing values of a correlation matrix
#' @importFrom DoE.wrapper lhs.design
#' @importFrom Matrix nearPD
overall_cor <- function(matrix) {
  # optim() function to calculate frobenius norm of difference matrix.
  sm <- function (mat, par) {
    # Store original covariance matrix and matrix that can be modified
    mat_par <- mat
    # Get the missing value locations from the upper triangle
    mat_par[lower.tri(mat_par)] <- 0
    index_na <- which(is.na(mat_par), arr.ind = TRUE)
    # Restore modified matrix and assign values
    mat_par <- mat
    for (i in 1:nrow(index_na)) {
      mat_par[index_na[i,1], index_na[i,2]] <- par[i]
      mat_par[index_na[i,2], index_na[i,1]] <- par[i]
    }
    # Difference between original matrix and nearest positive definite chosen matrix
    matt_diff <- mat - nearPD(mat_par, corr = TRUE)[["mat"]]
    # Calculate Frobenius norm of difference matrix
    frob_norm <- sum(matt_diff^2, na.rm = TRUE)^(1/2)
    frob_norm
  }

  # Initial values
    n_initial <- length(which(is.na(matrix)))/2
    par <- lhs.design(
      nruns = n_initial,
      nfactors = 1,
      default.levels = c(-1, 1)
    )

    # Find values which minimize the frobenius norm
    val <- optim(
      par = par[[1]],
      mat = matrix,
      fn = sm,
      lower = -1,
      upper = 1,
      method = "L-BFGS-B"
    )

    # Put the estimated values back in original matrix
    mat_optim <- matrix
    mat_optim[lower.tri(mat_optim)] <- 0
    index_na <- which(is.na(mat_optim), arr.ind = TRUE)
    mat_optim <- matrix
    for (j in 1:nrow(index_na)) {
      mat_optim[index_na[j,1], index_na[j,2]] <- val[["par"]][j]
      mat_optim[index_na[j,2], index_na[j,1]] <- val[["par"]][j]
    }
    if(!is_pos_semidef(mat_optim)) {
      message("The estimated covariance matrix is not positive semidefinite!\n",
              "It has been replaced by the nearest positive definite matrix.")
      mat_optim <- nearPD(mat_optim, corr = TRUE)$mat
    }
  matrix_est <- mat_optim
  matrix_est
}
