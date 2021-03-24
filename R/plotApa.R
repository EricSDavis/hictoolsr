#' Plots Aggregate Peak Analysis Matrix
#'
#' @param params optional "bb_params" object containing relevant function parameters
#' @param data matrix, list of matricies, or 3 column data.frame of APA results
#' @param x numeric or unit object specifying x-location of plot
#' @param y numeric or unit object specifying y-location of plot
#' @param width numeric or unit object specifying width of plot
#' @param height numeric or unit object specifying height of plot
#' @param just string or numeric vector specifying the justification of the viewport
#'  relative to its (x, y) location
#' @param default.units string indicating the default units to use if x, y, width, or
#'  height are only given as numeric vectors
#' @param draw logical value indicating whether graphics output should be produced
#' @param palette colorRampPalette function to use to map values to colors
#' @param n numeric by which to divide aggregate matrix
#' @param zrange vector of length 2; max and min values to set color scale
#'
#' @return Function will draw a APA matrix and return an object of class "bb_apa"
#'
#'
#' @examples
#'
#' ## Create divergent matrix ####
#' m <- matrix(data = rnorm(n = 21*21, mean = 0, sd = 2), nrow = 21, ncol = 21)
#'
#' ## Define parameters
#' p <- bb_params(width = 3, height = 3, default.units = "inches")
#'
#' ## Create page
#' bb_pageCreate(params = p)
#'
#' ## Plot apa
#' plot <- bb_plotApa(apa = m,
#'                    x = p$width/2, y = p$height/2,
#'                    width = p$width*0.5, height = p$width*0.5, just = c("center", "center"),
#'                    palette = colorRampPalette(c("blue", "white", "red")), zrange = NULL)
#'
#' ## Annotate legend
#' bb_annoHeatmapLegend(plot = plot,
#'                      x = 2.3, y = 0.75, width = 0.1, height = 0.75)
#'
#'
#' ## Create sequential matrix
#' m <- matrix(data = sample(0:100, 21*21, replace = T), nrow = 21, ncol = 21)
#'
#' ## Define parameters
#' p <- bb_params(width = 3, height = 3, default.units = "inches")
#'
#' ## Create page
#' bb_pageCreate(params = p)
#'
#' ## Plot apa
#' plot <- bb_plotApa(apa = m,
#'                    x = p$width/2, y = p$height/2,
#'                    width = p$width*0.5, height = p$width*0.5, just = c("center", "center"),
#'                    palette = colorRampPalette(c("white", "dark red")), zrange = NULL)
#'
#' ## Annotate legend
#' bb_annoHeatmapLegend(plot = plot,
#'                      x = 2.3, y = 0.75, width = 0.1, height = 0.75)
#'
#'
#' @export
#'
plotApa <- function(params = NULL, apa,
                       x = NULL, y = NULL, width = NULL, height = NULL,
                       just = c("left", "top"), default.units = "inches", draw = TRUE,
                       palette = colorRampPalette(c("white", "dark red")),
                       n = NULL, zrange=NULL) {

  ## Parse parameters & Create Object-----------------------------------------------------

  ## Get parsed arguments
  parsedArgs <- parseParams2(
    params = params,
    defaultArgs = formals(eval(match.call()[[1]])),
    declaredArgs = lapply(match.call()[-1], eval)
  )

  ## Evaluate parsed arguments
  parsedArgs <- lapply(parsedArgs, eval)

  ## Initialize object
  apa_plot <- structure(
    .Data = c(parsedArgs,
              list(
                color_palette = parsedArgs$palette,
                grobs = NULL
              )
    ),
    class = "bb_apa"
  )

  ## Set attributes for object
  attr(x = apa_plot, which = "plotted") <- parsedArgs$draw

  ## Parse units
  apa_plot <- defaultUnits(object = apa_plot, default.units = parsedArgs$default.units)


  ## Read in data & check for formatting errors ------------------------------------------

  ## Accept multiple types of input and convert to matrix
  if (class(apa_plot$apa) == "list") {

    ## Check for errors with list format
    check_apa_list(apa_plot$apa)

    ## Convert apa matrix list to a single apa matrix
    apa_plot$apa <- Reduce("+", apa_plot$apa)


  } else if ("data.frame" %in% class(apa_plot$apa)) {

    ## Check for errors with data.frame format
    check_apa_dataFrame(apa_plot$apa)

    ## Rename columns & cast into wide format matrix
    colnames(apa_plot$apa) <- c("var1", "var2", "value")
    apa_plot$apa <- reshape2::acast(apa_plot$apa, var1 ~ var2)

  } else if (is.null(apa_plot$apa)) {

    ## More specific error message if apa is null
    stop("argument 'apa' is missing, with no default.")

  } else if (class(apa_plot$apa) != "matrix") {

    ## Stop for anything that is not a matrix
    stop(sprintf("class %s is not supported.", class(apa_plot$apa)))

  }

  ## Check for matrix errors
  check_apa_matrix(apa_plot$apa)


  ## APA plot specific processing --------------------------------------------------------

  # Divide values by n
  if (!is.null(apa_plot$n)) {
    apa_plot$apa <- apa_plot$apa / n
  }

  ## Check for zrange errors
  check_apa_zrange(apa_plot$zrange)

  ## Set zrange if it is null
  apa_plot <- set_zrange(apa_plot)

  ## Map values to colors
  colv <- bb_maptocolors(vec = as.vector(apa_plot$apa),
                         col = apa_plot$color_palette,
                         num = 1000, range = apa_plot$zrange)

  ## Format color vector back to apa matrix
  m <- matrix(data = colv, nrow = nrow(apa_plot$apa), ncol = ncol(apa_plot$apa))

  ## Viewports ---------------------------------------------------------------------------

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_apa", length(grep(pattern = "bb_apa", x = currentViewports)) + 1)


  ## If placing information is provided but plot == TRUE,
  ## set up it's own viewport separate from bb_makepage

  ## Not translating into page_coordinates
  if (is.null(apa_plot$x) & is.null(apa_plot$y)){

    vp <- viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   height = unit(1, "snpc"), width = unit(1, "snpc"),
                   clip = "on", just = "center", name = vp_name)

    if (apa_plot$draw == TRUE){

      vp$name <- "bb_hic1"
      grid.newpage()

    }

  } else {

    ## Check that BentoBox page exists
    check_bbpage("Use bb_pageCreate() to make a BentoBox page to place a plot.")

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = apa_plot)

    ## Make viewport
    vp <- viewport(x = page_coords$x, y = page_coords$y,
                   height = page_coords$height, width = page_coords$width,
                   clip = "on", just = apa_plot$just,
                   name = vp_name)
  }


  ## Handle graphical objects ------------------------------------------------------------

  ## Initialize gTree for grobs
  assign("apa_grobs", gTree(vp = vp), envir = BentoBox:::bbEnv)

  ## Make grobs
  apaRaster <- rasterGrob(image = m, interpolate = F)

  ## Assign grobs to gTree
  assign(x = "apa_grobs",
         value = addGrob(gTree = get("apa_grobs", envir = BentoBox:::bbEnv),
                         child = apaRaster),
         envir = BentoBox:::bbEnv)

  ## Add grobs to object
  apa_plot$grobs <- get("apa_grobs", envir = BentoBox:::bbEnv)

  ## Plot grobs
  if (apa_plot$draw) {
    grid.draw(apa_plot$grobs)
  }

  ## Return object
  return(apa_plot)

}
