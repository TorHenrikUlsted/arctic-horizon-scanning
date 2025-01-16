#------------------------#
####      ggplot2     ####
#------------------------#

save_ggplot <- function(save.plot, save.name, save.width, save.height, save.dir, save.device = "jpeg", save.unit = "px", vis.title = FALSE, plot.show = FALSE, verbose = FALSE) {
  vebprint(save.plot, veb = plot.show)
  
  fig_out <- paste0(save.dir, "/", save.name, ".", save.device)
  
  create_dir_if(dirname(fig_out))
  
  catn("Saving plot to:", colcat(fig_out, color = "output"))
  ggsave(fig_out, device = save.device, unit = save.unit, width = save.width, height = save.height, plot = save.plot)
}

dynamic_guide_legend <- function(guide_config) {
  args <- list()
  
  for (param_name in names(guide_config)) {
    
    # Add the parameter to the args list
    args[[param_name]] <- guide_config[[param_name]]
  }
  
  # Create and return the guide_legend call
  do.call(guide_legend, args)
}

ggplot.filler <- function(gradient = "viridis-B", scale.type = "fill-c", limits = NULL, breaks = NULL, labels = NULL, begin = NULL, end = NULL, trans = NULL, guide = NULL, na.value = "transparent") {
  tryCatch(
    {
      # Syntax is: "gradient-option"
      split_str <- str_split(gradient, "-")[[1]]
      gradient <- split_str[[1]]
      option <- toupper(split_str[[2]])
      
      # Syntax is: "type-variable"
      split_str <- str_split(scale.type, "-")[[1]]
      scale_type <- split_str[[1]]
      scale_var <- tolower(split_str[[2]])
      
      # Process the guide
      if (!is.null(guide)) {
        if (is.list(guide)) {
          guide_obj <- dynamic_guide_legend(guide)
        } else if (is.character(guide) || is.logical(guide)) {
          guide_obj <- guide
        } else {
          warning("Invalid guide specification. Using default guide.")
          guide_obj <- "legend"
        }
      }
      
      args <- list(
        option = option,
        guide = if (!is.null(guide)) guide_obj else NULL,
        na.value = na.value
      )
      
      if (!is.null(labels)) {
        args$labels <- labels
      }
      
      if (!is.null(limits)) {
        args$limits <- limits
      }
      
      if (!is.null(breaks)) {
        args$breaks <- breaks
      }
      
      if (!is.null(begin)) {
        args$begin <- begin
      }
      
      if (!is.null(end)) {
        args$end <- end
      }
      
      if (!is.null(trans)) {
        args$trans <- trans
      }
      
      if (gradient == "viridis") {
        fun <- paste0("scale_", scale_type, "_viridis_", scale_var)
        return(do.call(fun, args))
      } else if (gradient == "whitebox") {
        fun <- paste0("scale_", scale_type, "_whitebox_", scale_var)
        args$palette <- args$option
        args$option <- NULL
        return(do.call(fun, args))
      }
    },
    error = function(e) {
      vebcat("Error when trying to use custom ggplot.filler function.", color = "fatalError")
      stop(e)
    }
  )
}

ggplot.theme <- function() {
  return(
    ggplot2::theme(
      text = element_text(family = config$ggplot$theme$text$family),
      
      plot.title = element_text(
        color = config$ggplot$theme$plot.title$color,
        vjust = config$ggplot$theme$plot.title$vjust,
        hjust = config$ggplot$theme$plot.title$hjust,
        size = config$ggplot$theme$plot.title$size,
        face = config$ggplot$theme$plot.title$face,
        margin = margin(b = config$ggplot$theme$plot.title$margin$b),
        lineheight = config$ggplot$theme$plot.title$lineheight
      ),
      
      plot.title.position = config$ggplot$theme$plot.title.position,
      plot.margin = do.call(margin, c(config$ggplot$theme$plot.margin[c("t", "r", "b", "l")], unit = config$ggplot$theme$plot.margin$unit)),
      
      axis.text = element_text(size = config$ggplot$theme$axis.text$size),
      axis.title.x = element_text(
        size = config$ggplot$theme$axis.title.x$size,
        hjust = config$ggplot$theme$axis.title.x$hjust
      ),
      
      axis.title.y = element_text(size = config$ggplot$theme$axis.title.y$size),
      legend.text = element_text(size = config$ggplot$theme$legend.text$size),
      legend.title = element_text(
        size = config$ggplot$theme$legend.title$size,
        hjust = config$ggplot$theme$legend.title$hjust
      ),
      
      legend.position = config$ggplot$theme$legend.position
    )
  )
}