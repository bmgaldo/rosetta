#' CV Individual
#'
#' Extracts the environmental and polygenetic variance matrix from a combination of dataframes and models.
#'
#' @param rosetta_data A list of named lists, e.g. list(list(data = x, model = y)). See the `rosetta_data()` function.
#'
#' @return Returns a list of lists in CV format.
#'
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' cv_individual()

cv_individual <- function(rosetta_data) {
  cv <- list()

  for (i in 1:length(rosetta_data)) {
    probands = list()
    aName = ascertainment(rosetta_data[[i]][["model"]])
    if( class(aName) == "character" )
    {
      if( !(aName %in% names(dataVals(rosetta_data[[i]][["data"]]))) )
      {
        stop(paste("Model contains ascertainment but no column of that name exist in data. ",
                   "Either set ascertainment in the model to NULL",
                   "or use a data set with the required proband column.",
                   sep = " "))
      } else
      {
        aValues = dataVals(rosetta_data[[i]][["data"]])[,aName, drop=FALSE]
        if( any(is.na(aValues)) )
        {
          stop(paste("Data contains the missing values(s) in the column '",
                     aName, "'. ",
                     "Complete data is required to run a strum model with ",
                     "ascertainment. Please check your data!",
                     sep = ""))
        }

        aValues = split(aValues, dataVals(rosetta_data[[i]][["data"]])$family)

        probands = lapply(aValues, function(pk) return(pk==1))

        pCount = unlist(lapply(aValues, sum))
        if( any(pCount > 1) )
          warning("Pedigree(s) with more than 1 proband exist!")
        if( any(pCount== 0) )
          warning("Pedigree(s) with no proband exist!")
      }

    } else if( !is.null(aName) )
    {
      warning("Model contains a wrongly specified ascertainment!  Ignoring...")
    }

    allRE = allRandomEffects(rosetta_data[[i]][["model"]]) # Generic method from strum

    y     = .getAnalysisY( rosetta_data[[i]][["model"]], rosetta_data[[i]][["data"]], probands)
    x     = .getAnalysisX( rosetta_data[[i]][["model"]], rosetta_data[[i]][["data"]], probands)
    vcPEC = .getAnalysisVC(rosetta_data[[i]][["model"]], rosetta_data[[i]][["data"]], allRE)

    vcPro = list()
    if( length(probands) > 0 )
      vcPro = mapply(.filterMissingVC, vcPEC, probands, probands, SIMPLIFY=FALSE)

    vc = list(vcAll = vcPEC, vcPro = vcPro)
    filtered = .getValidAnalysisData(y, x, vc)

    step1OptimControl = list(maxit=5500, fnscale=-10)

    cv[[i]] <- .estimateDeltaParameter(filtered$y, filtered$x, filtered$vc, step1OptimControl)

    # Assign names to rows and columns of variance matrices and convert to correlation
    for (j in 1:length(cv[[i]]$V)) {
      colnames(cv[[i]]$V[[j]]) <- rosetta_data[[i]][["model"]]@varList[rosetta_data[[i]][["model"]]@varList$obs == "TRUE", "name"]
      rownames(cv[[i]]$V[[j]]) <- rosetta_data[[i]][["model"]]@varList[rosetta_data[[i]][["model"]]@varList$obs == "TRUE", "name"]
      if(!is_pos_semidef(cv[[i]]$V[[j]])) {
        message("The estimated covariance matrix is not positive semidefinite!\n",
                paste("combination", i, "component", j, "\n", sep = " "),
                "It has been replaced by the nearest positive definite matrix.")
        cv[[i]]$V[[j]] <- Matrix::nearPD(cv[[i]]$V[[j]], corr = FALSE)$mat
      }
      #cv[[i]]$V[[j]] <- cov2cor(cv[[i]]$V[[j]])
    }

    # Assign names to mean vector
    for (j in 1:length(cv[[i]]$V)) {
      colnames(cv[[i]]$C) <- rosetta_data[[i]][["model"]]@varList[rosetta_data[[i]][["model"]]@varList$obs == "TRUE", "name"]
    }
  }
  #return
  cv
}
