#' @title Reading a CSV file and remove spaces from the column names
#' @description Reading a CSV file and remove spaces from the column names
#' @param fileName a CSV file name
#' @importFrom dplyr %>%
#' @importFrom stringr str_replace_all
#' @importFrom magrittr set_colnames
#' @importFrom readr read_csv
#' @export
read_csv_checkColumnNames <- function(fileName){
  readr::read_csv(fileName) %>%
    (function(d){magrittr::set_colnames(d, stringr::str_replace_all(colnames(d), " ", "."))})
}

#' @title Combine multiple ggplots
#' @description \code{multiplot} combines multiple ggplots.
#' @param ... ggplots
#' @param plotlist A list of ggplots
#' @param cols The number of columns
#' @param layout A matrix for layouting the plots
#' @importFrom grid grid.newpage
#' @importFrom grid pushViewport
#' @importFrom grid viewport
#' @export
multiPlot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' @title Plot saving utilities
#' @description
#' \code{savePDF} saves a plot to a PDF file.
#' @param graphicObject The plot object to be saved.
#' @param outputFileName The name of the exported PDF file
#' @param width A plot width
#' @param height A plot height
#' @param pointsize A pointsize
#' @param family A font family
#' @param useDingbats Logical indicating whether points shuold be replaced with the Dingbat font.
#' @importFrom ggplot2 last_plot
#' @importFrom stringr str_detect
#' @importFrom grDevices dev.copy2pdf
#' @export
savePDF <- function(graphicObject=ggplot2::last_plot(), outputFileName, width, height,
                    pointsize=12, family="Helvetica", useDingbats=F){
  windowsFonts("UserDefinedFont"=windowsFont(family))
  out <- ifelse(stringr::str_detect(outputFileName, ".pdf$"), outputFileName, paste0(outputFileName, ".pdf"))
  print(graphicObject)
  grDevices::dev.copy2pdf(file=out, width=width, height=height, pointsize=pointsize, family="UserDefinedFont", useDingbats=useDingbats)
}

#' @title Plot saving utilities
#' @description
#' \code{saveCurrentGraphicPDF} saves the plot currently on the graphics device to a PDF file.
#' @param outputFileName The name of the exported PDF file
#' @param width A plot width
#' @param height A plot height
#' @param pointsize A pointsize
#' @param family A font family
#' @param useDingbats Logical indicating whether points shuold be replaced with the Dingbat font.
#' @importFrom ggplot2 last_plot
#' @importFrom stringr str_detect
#' @importFrom grDevices dev.copy2pdf
#' @export
saveCurrentGraphicPDF <- function(outputFileName, width, height,
                                  pointsize=12, family="Helvetica", useDingbats=F){
  windowsFonts("UserDefinedFont"=windowsFont(family))
  out <- ifelse(stringr::str_detect(outputFileName, ".pdf$"), outputFileName, paste0(outputFileName, ".pdf"))
  grDevices::dev.copy2pdf(file=out, width=width, height=height, pointsize=pointsize, family="UserDefinedFont", useDingbats=useDingbats)
}

#' @title Publication-ready ggplot scales
#' @description Publication-ready ggplot scales.
#' @param ... No extra arguments are supported.
#' @importFrom ggplot2 discrete_scale
#' @importFrom scales manual_pal
#' @export
scale_color_Publication <- function(...) {
  discrete_scale("color", "Publication", manual_pal(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

#' @title Publication-ready ggplot scales
#' @description Publication-ready ggplot scales.
#' @param ... No extra arguments are supported.
#' @importFrom ggplot2 discrete_scale
#' @importFrom scales manual_pal
#' @export
scale_colour_Publication <- function(...) {
  discrete_scale("colour", "Publication", manual_pal(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

#' @title Publication-ready ggplot scales
#' @description Publication-ready ggplot scales.
#' @param ... No extra arguments are supported.
#' @importFrom ggplot2 discrete_scale
#' @importFrom scales manual_pal
#' @export
scale_fill_Publication <- function(...) {
  discrete_scale("fill", "Publication", manual_pal(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

#' @title The publication-ready ggplot theme
#' @description The publication-ready ggplot theme.
#' @param base_size The base font size.
#' @param base_family The base font family.
#' @import ggplot2
#' @importFrom grDevices windowsFonts
#' @importFrom grDevices windowsFont
#' @importFrom grid gpar
#' @importFrom grid polylineGrob
#' @importFrom utils modifyList
#' @importFrom ggthemes theme_foundation
#' @export
theme_Publication <-
  function(base_size=18, base_family="Helvetica") {
    # Translation of a device-independent R graphics font family name to a windows font description
    windowsFonts("UserDefinedFont"=windowsFont(base_family))

    # The preset ggplot parameters
    .pt <- 1 / 0.352777778
    .len0_null <- function(x) {
      if (length(x) == 0)  NULL
      else                 x
    }

    # A ggplot theme for the "L" (left + bottom) border
    theme_L_border <- function(colour = "black", size = 1, linetype = 1) {
      # use with e.g.: ggplot(...) + theme( panel.border=theme_L_border() ) + ...
      structure(
        list(colour = colour, size = size, linetype = linetype),
        class = c("theme_L_border", "element_blank", "element")
      )
    }
    element_grob.theme_L_border <- function(
      element, x = 0, y = 0, width = 1, height = 1,
      colour = NULL, size = NULL, linetype = NULL,
      ...) {
      gp <- gpar(lwd = .len0_null(size * .pt), col = colour, lty = linetype)
      element_gp <- gpar(lwd = .len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
      polylineGrob(
        x = c(x+width, x, x), y = c(y,y,y+height), ..., default.units = "npc",
        gp = modifyList(element_gp, gp)
      )
    }

    # A combined set of ggplot theme options
    (
      theme_foundation(base_size = base_size, base_family = "UserDefinedFont")
      + theme(
        plot.title = element_text(
          face = "bold",
          size = rel(1.2), hjust = 0.5
        ),
        text = element_text(),
        plot.background = element_rect(colour = NA),
        plot.margin = unit(c(10,5,5,5), "mm"),
        panel.background = element_rect(colour = NA),
        panel.border = theme_L_border(),
        panel.grid.major = element_line(size = .5, colour = "#f0f0f0"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(vjust = -0.25),
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(),
        axis.line = element_line(size = .7, colour = "black"),
        axis.ticks = element_line(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.background = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.key.size = unit(0.7, "cm"),
        legend.spacing = unit(0.5, "cm"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
        strip.text = element_text(face = "bold")
      )
    )
}
