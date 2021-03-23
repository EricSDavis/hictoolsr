## Utility functions for bb_plotApa

## Error checking functions --------------------------------------------------------------

## Define function that checks apa for list format
check_apa_list <- function(apa) {

  ## Check that the class of each element is a matrix
  if (!unique(unlist(lapply(apa, function(x) "matrix" %in% class(x))))) {
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

## From BentoBox
## Define a function that converts coordinates/dimensions into default units
defaultUnits <- function(object, default.units){

  if (!(is.null(object$x) & is.null(object$y))){

    if (!"unit" %in% class(object$x)){

      if (!is.numeric(object$x)){

        stop("x-coordinate is neither a unit object nor a numeric value. Cannot place object.", call. = FALSE)

      }

      if (is.null(default.units)){

        stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      object$x <- unit(object$x, default.units)

    }


    if (!"unit" %in% class(object$y)){

      ## Check for "below" y-coord
      if (grepl("b", object$y) == TRUE){
        if (grepl("^[ac-zA-Z]+$", object$y) == TRUE){
          stop("\'below\' y-coordinate detected with additional letters. Cannot parse y-coordinate.", call. = FALSE)
        }

        if(is.na(as.numeric(gsub("b","", object$y)))){
          stop("\'below\' y-coordinate does not have a numeric associated with it. Cannot parse y-coordinate.", call. = FALSE)
        }

        object$y <- plot_belowY(y_coord = object$y)

      } else {

        if (!is.numeric(object$y)){

          stop("y-coordinate is neither a unit object nor a numeric value. Cannot place object.", call. = FALSE)

        }

        if (is.null(default.units)){

          stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

        }

        object$y <- unit(object$y, default.units)

      }


    }

    if (!"unit" %in% class(object$width)){

      if (!is.numeric(object$width)){

        stop("Width is neither a unit object nor a numeric value. Cannot place object.", call. = FALSE)

      }

      if (is.null(default.units)){

        stop("Width detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      object$width <- unit(object$width, default.units)

    }

    if (!"unit" %in% class(object$height)){

      if (!is.numeric(object$height)){

        stop("Height is neither a unit object nor a numeric value. Cannot place object.", call. = FALSE)

      }

      if (is.null(default.units)){

        stop("Height detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      object$height <- unit(object$height, default.units)

    }

  }

  return(object)

}

## From BentoBox
## Define a function that maps a vector to colors
bb_maptocolors <- function(vec, col, num = 100, range = NULL){

  if (is.null(range) == TRUE){
    breaks <- seq(min(vec), max(vec), length.out = num)
  } else {
    vec[which(vec < range[1])] = range[1]
    vec[which(vec > range[2])] = range[2]
    breaks <- seq(range[1], range[2], length.out = num)
  }

  cols <- col(length(breaks) + 1)
  colvec <- as.character(cut(vec, c(-Inf, breaks, Inf), labels = cols))
  return(colvec)


}

## From BentoBox
## Define a function to grab the name of a viewport
viewport_name <- function(viewport){

  return(viewport$name)

}

## From BentoBox
## Define a function to get a list of current viewports
current_viewports <- function(){

  if (!"bb_page" %in% names(lapply(current.vpTree()$children, viewport_name))){

    current <- as.list(names(lapply(current.vpTree()$children, viewport_name)))

  } else {

    ## Check for groups
    page_children <- names(lapply(current.vpTree()$children$bb_page$children, viewport_name))

    if (length(grep(pattern = "bb_group", x = page_children)) > 0){

      group_vps <- as.list(page_children[grep(pattern = "bb_group", x = page_children)])

      group_children <- unlist(lapply(group_vps, vp_children), recursive = F)

      children_vps <- lapply(group_children, viewport_name)

      current <- c(page_children, children_vps)

    } else {

      current <- as.list(page_children)

    }

  }

  return(current)
}

## From BentoBox
## Define a function to make sure a bb_page viewport exists
check_bbpage <- function(error){

  if (!"bb_page" %in% current.vpPath()){

    stop(error, call. = FALSE)

  }

}

## From BentoBox
## Define a function to convert to page units
convert_page <- function(object){



  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = BentoBox:::bbEnv)
  page_units <- get("page_units", envir = BentoBox:::bbEnv)

  ## Convert x and y coordinates and height and width to same page_units
  old_x <- object$x
  old_y <- object$y
  old_height <- object$height
  old_width <- object$width
  new_x <- convertX(old_x, unitTo = page_units)
  new_y <- convertY(unit(page_height, units = page_units) - convertY(old_y, unitTo = page_units), unitTo = page_units)
  new_height <- convertHeight(old_height, unitTo = page_units)
  new_width <- convertWidth(old_width, unitTo = page_units)

  object$x <- new_x
  object$y <- new_y
  object$height <- new_height
  object$width <- new_width


  return(object)

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
