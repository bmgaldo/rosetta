#' Rosetta Data
#'
#' Creates a named list of strum dataframes and strum models for each combination of supplied datasets.
#'
#' @param data_list A list of dataframes that contain a varying number of measured variables in respect to the underlying factor structure.
#' @param factor_structure A named list for the hypothesized factor structure. Set the factor names as the list names and measured variables as the list elements.
#' @param twin The name of the twin zygotic variable as a character string.
#'
#' @return Returns named list of strum dataframes and strum models for each combination of supplied datasets.
#'
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' rosetta_data()

rosetta_data <- function(data_list, factor_structure, twin = NULL) {
  n_df <- length(data_list)

  ## List of indices for combinations of dataframes
  comb_index <- 1:n_df
  combs <- lapply(comb_index, function(z) {combn(length(comb_index), z)})

  ## create list store for dataframes and Models
  list_store <- list()
  list_store_index <- 1

  ### Create strum data and model for each individual dataframe
  for (i in 1:n_df) {
    ## Dataframes
    list_item <- list(data =
                        strum_data(
                          data_list[[i]],
                          type = "Pedigree",
                          twin = twin
                        )
    )

    ## Models
    factors <- lapply(factor_structure, function(z) {intersect(z, names(data_list[[i]]))}) # Remove unused variables from factor structure
    factors <- Filter(length, factors) # Remove all empty elements in list
    model_text <- strum_model(factors)
    list_item[["model"]] <- createStrumModel(
      formulas = model_text,
      fixLoadingToOne = FALSE,
      defaultError = '<e,p>'
    )
    list_store[[i]] <- list_item
    names(list_store)[list_store_index] <- paste(i)
    list_store_index <- list_store_index + 1
  }

  if (length(combs) > 1) {
    ### Create strum data and model for each combination of dataframes
    #### Create list of names for merging the dataframes
    col_names <- lapply(data_list, colnames)
    merged_list_names <- list()
    merged_list_names_index <- 1
    for (i in 2:n_df) {
      for (j in 1:ncol(combs[[i]])) {
        merged_list_names[[merged_list_names_index]] <- Reduce(intersect, col_names[combs[[i]][, j]])
        names(merged_list_names)[merged_list_names_index] <- paste(unlist(combs[[i]][, j]), collapse = "_") # method downsides? https://stat.ethz.ch/pipermail/r-help/2001-March/011480.html
        merged_list_names_index <- merged_list_names_index + 1
      }
    }

    #### Create merged dataframes and models
    df_merged <- list()
    df_merged_index <- 1
    for (i in 2:n_df) {
      for (j in 1:ncol(combs[[i]])) {
        ## Dataframes
        df_list_combn <- data_list[combs[[i]][, j]]
        df_list_combn <- lapply(df_list_combn, function(z) {
          z[, names(z) %in% merged_list_names[[df_merged_index]]]
        })
        df_list_combn_merged <- do.call("rbind", df_list_combn)
        list_item <- list(data =
                            strum_data(
                              df_list_combn_merged,
                              type = "Pedigree",
                              twin = twin
                            )
        )

        ## Models
        factors <- lapply(factor_structure, function(z) {intersect(z, names(df_list_combn_merged))})
        factors <- Filter(length, factors) # Remove all empty elements in list
        model_text <- strum_model(factors)

        list_item[["model"]] <- createStrumModel(
          formulas = model_text,
          fixLoadingToOne = FALSE,
          defaultError = '<e,p>'
        )
        list_store[[list_store_index]] <- list_item
        names(list_store)[list_store_index] <- names(merged_list_names)[df_merged_index]
        list_store_index <- list_store_index + 1
        df_merged_index <- df_merged_index + 1
      }
    }
  }
  list_store
}
