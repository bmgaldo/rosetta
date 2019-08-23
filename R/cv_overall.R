#' CV Overall
#'
#' Creates an overall CV object based on data calculated from the greatest number of observation.
#'
#' @param cv_individual See the `cv_individual()` function.
#' @param factor_structure A named list for the hypothesized factor structure. Set the factor names as the list names and measured variables as the list elements.
#'
#' @return Returns a CV object.
#'
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' cv_overall()

cv_overall <- function(cv_individual, factor_structure) {
  # The rows and columns of the covariance matrix
  rownames <- as.vector(unlist(factor_structure))
  colnames <- as.vector(unlist(factor_structure))

  ## Polygenic V
  for (i in 1:length(cv_individual)) {
    # Creates a 3 column dataframe of the indices and values needed to create a
    # covariance matrix
    dat <- data.frame(
      expand.grid(
        rownames=rownames(cv_individual[[i]]$V[[1]]),
        colnames=colnames(cv_individual[[i]]$V[[1]]),
        stringsAsFactors = FALSE
      ),
      values=as.vector(cv_individual[[i]]$V[[1]])
    )

    # Fill in matrix using values from dataframe
    mat_new <- outer(
      rownames,
      colnames,
      function(x,y) {
        mapply(function(x_sub, y_sub) {
          val <- dat[dat[, 1] == x_sub & dat[, 2] == y_sub, 3]
          if(length(val) == 0L) {
            NA_integer_
          } else {
            val
          }
        },
        x,
        y
        )
      }
    )

    # Create a covariance matrix using values that were calculated from the most
    # amount of data.
    if (i > 1) {
      poly_mat <- ifelse(!is.na(mat_new), mat_new, poly_mat)
    } else {
      poly_mat <- mat_new
    }
  }

  ## Assign row and column names
  rownames(poly_mat) <- rownames
  colnames(poly_mat) <- colnames

  ## Environmental V
  for (i in 1:length(cv_individual)) {
    # Creates a 3 column dataframe of the indices and values needed to create a
    # covariance matrix
    dat <- data.frame(
      expand.grid(
        rownames=rownames(cv_individual[[i]]$V[[2]]),
        colnames=colnames(cv_individual[[i]]$V[[2]]),
        stringsAsFactors = FALSE
      ),
      values=as.vector(cv_individual[[i]]$V[[2]])
    )

    # Fill in matrix using values from dataframe
    mat_new <- outer(
      rownames,
      colnames,
      function(x,y) {
        mapply(function(x_sub, y_sub) {
          val <- dat[dat[, 1] == x_sub & dat[, 2] == y_sub, 3]
          if(length(val) == 0L) {
            NA_integer_
          } else {
            val
          }
        },
        x,
        y
        )
      }
    )

    # Create a covariance matrix using values that were calculated from the most
    # amount of data.
    if (i > 1) {
      env_mat <- ifelse(!is.na(mat_new), mat_new, env_mat)
    } else {
      env_mat <- mat_new
    }
  }

  ## Assign row and column names
  rownames(env_mat) <- rownames
  colnames(env_mat) <- colnames

  # Create a mean matrix using values that were calculated from the most
  # amount of data.
  mean_mat <- matrix(rep(NA_integer_, (length(colnames))), nrow = 1)
  colnames(mean_mat) <- colnames
  for (i in 1:length(cv_individual)) {
    mean_mat[, colnames(mean_mat) %in% colnames(cv_individual[[i]]$C)] <- cv_individual[[i]]$C
  }

  ## Return
  list(mean_mat = mean_mat, polygenic_mat = poly_mat, environmental_mat = env_mat)
}
