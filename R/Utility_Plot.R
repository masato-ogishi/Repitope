#' Plot saving utilities
#' 
#' \code{savePDF} saves a plot to a PDF file.\cr
#' \code{saveCurrentGraphicPDF} saves the plot currently on the graphics device to a PDF file.
#' 
#' @param graphicObject The plot object to be saved.
#' @param outputFileName The name of the exported PDF file
#' @param width A plot width
#' @param height A plot height
#' @param pointsize A pointsize
#' @param useDingbats Logical indicating whether points shuold be replaced with the Dingbat font.
#' @param ghostScriptPath The file path to the GhostScript.exe file.
#' @importFrom ggplot2 last_plot
#' @importFrom stringr str_detect
#' @importFrom grDevices dev.copy2pdf
#' @importFrom grDevices embedFonts
#' @export
#' @rdname Utility_Plot
#' @name Utility_Plot
savePDF <- function(graphicObject=ggplot2::last_plot(), outputFileName, width, height,
                    pointsize=12, useDingbats=F, ghostScriptPath="C:/gs/gs9.16/bin/gswin32c.exe"){
  Sys.setenv(R_GSCMD=ghostScriptPath)
  .gsOption <- "-sFONTPATH=C:/Windows/Fonts -dSubsetFonts=true -dEmbedAllFonts=true"
  windowsFonts(Helvetica=windowsFont("Helvetica"))
  out <- ifelse(stringr::str_detect(outputFileName, ".pdf$"), outputFileName, paste0(outputFileName, ".pdf"))
  print(graphicObject)
  grDevices::dev.copy2pdf(file=out, width=width, height=height, pointsize=pointsize, family="Helvetica", useDingbats=useDingbats)
  grDevices::embedFonts(out, outfile=out, options=.gsOption)
}
#' @export
#' @rdname Utility_Plot
#' @name Utility_Plot
saveCurrentGraphicPDF <- function(outputFileName, width, height,
                                  pointsize=12, useDingbats=F, ghostScriptPath="C:/gs/gs9.16/bin/gswin32c.exe"){
  Sys.setenv(R_GSCMD=ghostScriptPath)
  .gsOption <- "-sFONTPATH=C:/Windows/Fonts -dSubsetFonts=true -dEmbedAllFonts=true"
  windowsFonts(Helvetica=windowsFont("Helvetica"))
  out <- ifelse(stringr::str_detect(outputFileName, ".pdf$"), outputFileName, paste0(outputFileName, ".pdf"))
  grDevices::dev.copy2pdf(file=out, width=width, height=height, pointsize=pointsize, family="Helvetica", useDingbats=useDingbats)
  grDevices::embedFonts(out, outfile=out, options=.gsOption)
}

#' Publication-ready ggplot scales
#' @importFrom ggplot2 discrete_scale
#' @importFrom scales manual_pal
#' @export
#' @rdname Utility_Plot_ggplot-scales
#' @name Utility_Plot_ggplot-scales
scale_color_Publication <- function() {
  discrete_scale("color", "Publication", manual_pal(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}
#' @export
#' @rdname Utility_Plot_ggplot-scales
#' @name Utility_Plot_ggplot-scales
scale_colour_Publication <- function() {
  discrete_scale("colour", "Publication", manual_pal(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}
#' @export
#' @rdname Utility_Plot_ggplot-scales
#' @name Utility_Plot_ggplot-scales
scale_fill_Publication <- function() {
  discrete_scale("fill", "Publication", manual_pal(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

#' The publication-ready ggplot theme
#' @param base_size The base font size.
#' @import ggplot2
#' @importFrom grDevices windowsFonts
#' @importFrom grDevices windowsFont
#' @importFrom grid gpar
#' @importFrom grid polylineGrob
#' @importFrom utils modifyList
#' @importFrom ggthemes theme_foundation
#' @export
#' @rdname Utility_Plot_ggplot-theme
theme_Publication <-
  function(base_size=18) {
    # Translation of a device-independent R graphics font family name to a windows font description
    windowsFonts(Helvetica=windowsFont("Helvetica"))
    
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
      theme_foundation(base_size = base_size, base_family = "Helvetica")
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

