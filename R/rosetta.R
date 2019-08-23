#' Rosetta
#'
#' This function will run a Rosetta analysis based on the sem package.
#'
#' @param data_list A list of dataframes that contain a varying number of measured variables in respect to the underlying factor structure.
#' @param factor_structure A named list for the hypothesized factor structure. Set the factor names as the list names and measured variables as the list elements.
#' @param S A covariance matrix among observed variables.
#'
#' @return Returns factor scores for each observation in the dataframes
#'
#' @importFrom dplyr bind_rows
#' @importFrom sem sem fscores
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' rosetta()

rosetta <- function(data_list, factor_list, S = NULL) {
  # Run unconstrained SEM
  if (is.null(S)) {
    # Unlist data into a single dataframe
    data_bind <- as.data.frame(bind_rows(data_list))

    # Calculate pairwise correlations and estimate correlation matrix if necessary
    overall_cor_matrix <- cor(data_bind, use = "pairwise.complete.obs")
    if (anyNA(overall_cor_matrix)) {
      overall_cor_matrix <- overall_cor(overall_cor_matrix)
    }
  } else {
    overall_cor_matrix <- S
  }

  model <- sem_model(factor_list)
  overall_sem <- sem(
    model = model,
    S = overall_cor_matrix,
    N = ncol(overall_cor_matrix)
  )

  # Run constrained SEM for each dataset
  est_factor_score <- list()

  for (i in 1:length(data_list)) {
    factor_structure <- lapply(factor_list, function(z) {
      z[z %in% colnames(data_list[[i]])]
    })
    model <- sem_model(factor_structure, overall_sem)
    sem <- sem(model, data = data_list[[i]])
    est_factor_score[[i]] <- fscores(sem, data = data_list[[i]])
  }

  if (!is.null(names(data_list))) {
    names(est_factor_score) <- names(data_list)
  }

  # Return factor scores
  est_factor_score
}
