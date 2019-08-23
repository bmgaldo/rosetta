#' Covariance Estimate
#'
#' Estimates missing values in a covariance matrix
#'
#' @param cv_overall See the `cv_overall()` function.
#'
#' @return Returns an estimated CV object.
#'
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' cv_overall()

cov_estimate <- function(cv_overall) {
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

  # Remove mean matrix from cv_overall
  c <- list(cv_overall$mean_mat)
  cv_overall <- cv_overall[!(names(cv_overall) %in% "mean_mat")]

  # Estimate the polygenic and environmental covariance matrices
  cv_est <- list()
  for(i in 1:length(cv_overall)) {
    # initial values
    n_initial <- length(which(is.na(cv_overall[[i]])))/2
    par <- lhs.design(
      nruns = n_initial,
      nfactors = 1,
      default.levels = c(-1, 1)
    )

    # Find values which minimize the frobenius norm
    val <- optim(
      par = par[[1]],
      mat = cv_overall[[i]],
      fn = sm,
      lower = -1,
      upper = 1,
      method = "L-BFGS-B"
    )

    # Put the estimated values back in original matrix
    mat_optim <- cv_overall[[i]]
    mat_optim[lower.tri(mat_optim)] <- 0
    index_na <- which(is.na(mat_optim), arr.ind = TRUE)
    mat_optim <- cv_overall[[i]]
    for (j in 1:nrow(index_na)) {
      mat_optim[index_na[j,1], index_na[j,2]] <- val[["par"]][j]
      mat_optim[index_na[j,2], index_na[j,1]] <- val[["par"]][j]
    }
    if(!is_pos_semi_def(mat_optim)) {
      message("The estimated covariance matrix is not positive semidefinite!\n",
              "It has been replaced by the nearest positive definite matrix.")
      mat_optim <- Matrix::nearPD(mat_optim, corr = TRUE)$mat
    }
    cv_est[[i]] <- mat_optim
  }
  cv_est <- c(c = c, cv_est)
}
