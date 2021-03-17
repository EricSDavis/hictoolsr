## Utility functions for bb_plotApa

## Error checking functions --------------------------------------------------------------

## Define function that checks apa for list format
check_apa_list <- function(apa) {

  ## Check that the class of each element is a matrix
  if (!identical(unique(unlist(lapply(apa, class))), "matrix")) {
    stop("apa list must all be class 'matrix'.", call. = FALSE)
  }

  ## Check that all elements are of the same, square dimensions
  if (length(unique(unlist(lapply(apa, dim)))) != 1) {
    stop("apa list must contain matricies of the same dimensions and nrow == ncol.", call. = FALSE)
  }

  ## Check that all elements are numeric
  if (!any(c("numeric", "integer", "double") %in% class(unlist(apa)))) {
    stop("apa list must be a list of numeric matricies.", call. = FALSE)
  }

}

## Define function that checks apa for data.frame format
check_apa_dataFrame <- function(apa) {

  ## Check that format is long (expects 3 column data.frame)
  if (ncol(apa) != 3) {
    stop("apa data.frame expects long format (i.e. bin1, bin2, value)", call. = FALSE)
  }

  ## Check that third column is numeric
  if (!any(c("numeric", "integer", "double") %in% class(apa[[3]]))) {
    stop("apa data.frame values must be numeric.", call. = FALSE)
  }

}

## Define function that checks apa for matrix format
check_apa_matrix <- function(apa) {

  ## Check that format is wide
  if (nrow(apa) != ncol(apa)) {
    stop("apa matrix must have equal numbers of rows and columns (i.e. wide format).", call. = FALSE)
  }

  ## Check that all elements are numeric
  if (!any(c("numeric", "integer", "double") %in% class(as.vector(apa)))) {
    stop("apa matrix must be a numeric.", call. = FALSE)
  }

}

## Define function to check zrange
check_apa_zrange <- function(zrange) {

  ## Ensure properly formatted zrange
  if (!is.null(zrange)){

    ## zrange needs to be a vector
    if (!is.vector(zrange)){

      stop("'zrange' must be a vector of length 2.", call. = FALSE)

    }

    ## zrange vector needs to be length 2
    if (length(zrange) != 2){

      stop("'zrange' must be a vector of length 2.", call. = FALSE)

    }

    ## zrange vector needs to be numbers
    if (!is.numeric(zrange)){

      stop("'zrange' must be a vector of two numbers.", call. = FALSE)

    }

    ## second value should be larger than the first value
    if (zrange[1] >= zrange[2]){

      stop("'zrange' must be a vector of two numbers in which the 2nd value is larger than the 1st.", call. = FALSE)

    }

  }
}


## bb_plotApa functions ------------------------------------------------------------------

## Define a function that sets the zrange if it is null
set_zrange <- function(apa_plot){

  ## Set zrange if it is null
  if (is.null(apa_plot$zrange)) {

    ## Set matrix counts variable for convenience
    cnts <- as.vector(apa_plot$apa)

    ## Only one value
    if (length(unique(cnts)) == 1) {
      apa_plot$zrange <- c(unique(cnts), unique(cnts))
    }

    ## Multiple values
    if (length(unique(cnts)) > 1) {

      ## Divergent data (center on 0 by default)
      if (all(c(-1, 1) %in% unique(sign(cnts))))
      {
        apa_plot$zrange <- c(floor(-max(abs(range(cnts)))),
                             ceiling(max(abs(range(cnts)))))
      }

      ## Positive sequential data
      else if (!-1 %in% unique(sign(cnts)))
      {
        apa_plot$zrange <- c(0, ceiling(max(cnts)))
      }

      ## Negative sequential data
      else
      {
        apa_plot$zrange <- c(floor(min(cnts)), 0)
      }
    }
  }

  return(apa_plot)

}



## Parameter parsing function (until updated) --------------------------------------------

## Define parameter parsing function
parseParams2 <- function(params = params,
                         defaultArgs = formals(eval(match.call()[[1]])),
                         declaredArgs = lapply(match.call()[-1], eval)) {

  ## Place this function inside the parent function
  ##
  ## defaultArgs: List of defaults for each argument of parent function
  ## params: bb_params object to override default arguments of parent function
  ## declaredArgs: List of arguments to override all others

  ## Remove 'params' and '...' from defaultArgs and declaredArgs
  defaultArgs[["params"]] <- NULL
  declaredArgs[["params"]] <- NULL
  if ("..." %in% names(defaultArgs)) defaultArgs[["..."]] <- NULL
  if ("..." %in% names(declaredArgs)) declaredArgs[["..."]] <- NULL

  ## If bb_params are supplied override matching defaultArguments
  if (!is.null(params)) {

    if (class(params) == "bb_params") {

      ## Replace matching defaultArgs with params
      matchedParams <- params[na.omit(sort(match(names(defaultArgs), names(params))))]
      defaultArgs[na.omit(match(names(params), names(defaultArgs)))] <- matchedParams

    } else {

      stop("'params' must be an object of class bb_params.")

    }

  }

  ## Replace default args with declared args
  if (length(declaredArgs) != 0) {
    defaultArgs[names(defaultArgs) %in% names(declaredArgs)] <-  declaredArgs
  }

  ## Set arguments without default to NULL
  unset <- unlist(lapply(defaultArgs, is.name))
  defaultArgs[unset] <- lapply(lapply(defaultArgs[unset], deparse), as.null)

  ## Return final set of overriden arguments
  return(defaultArgs)

}
