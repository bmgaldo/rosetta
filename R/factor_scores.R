#' Factor Scores
#'
#' Calculate Factor Scores
#'
#' @param model An sem object.
#'
#' @return Returns factor scores for each observation in the dataframes
#'
#' @importFrom sem specifyModel
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' sem_model()

factor_scores <- function(model) {
  m <- model$m
  P <- model$P
  A <- model$A
  var.names <- model$var.names
  observed <- var.names %in% rownames(model$C)
  if (all(observed)) stop("there are no latent variables")
  IAinv <- solve(diag(m) - A)
  Sigma <- IAinv %*% P %*% t(IAinv)
  B <- solve(Sigma[observed, observed]) %*% Sigma[observed, !observed]
  rownames(B) <- var.names[observed]
  colnames(B) <- var.names[!observed]
  B
}

#' @export
factor_score_reliability <- function(lambda, phi, predictor) {
  mdiag <- function(x) {
    diag(diag(x))
  }
  inverse <- function(x) {solve(x)}

  if(is(lambda, "loadings")) {
    lambda <- lambda[, ]
  }

  sigma <- lambda %*% phi %*% t(lambda)
  sigma <- sigma - mdiag(sigma) + diag(nrow(lambda))

  psi <- mdiag(sigma - lambda %*% phi %*% t(lambda))^0.5

  if(round(min(diag(psi))) < 0) {
    stop("The diagonal of psi contains negative values.")
  }

  regress <- inverse(mdiag(phi %*% t(lambda) %*% inverse(sigma) %*% lambda %*% phi))^0.5 %*%
             mdiag(phi %*% t(lambda) %*% inverse(sigma) %*% lambda %*% phi %*% t(lambda) %*% inverse(sigma) %*% lambda %*% phi) %*%
    inverse(mdiag(phi %*% t(lambda) %*% inverse(sigma) %*% lambda %*% phi))^0.5
  diag(regress)
}
