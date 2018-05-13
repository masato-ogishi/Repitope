#' A utility function for creating a Venn diagram.
#'
#' @param x A list of vectors (e.g., integers, chars), with each component corresponding to a separate circle in the Venn diagram. The length should fall into the range between 1 and 5.
#' @param colors A set of colors for the venn diagram.
#' @param show_category_names Logical. Whether the names of \code{x} should be printed?
#' @export
#' @rdname Utility_VennDiagram
#' @name Utility_VennDiagram
vennDiagram <- function(x, colors=ggsci::pal_d3()(length(x)), show_category_names=T){
  try(dev.off(), silent=T) ## remove previous plots if any
  if(show_category_names==T){
    if(identical(names(x), NULL)){
      catNames <- 1:length(x)
    }else{
      catNames <- names(x)
    }
  }else{
    catNames <- rep("", length(x))
  }
  if(length(x)==5){
    cexOption <- c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5)
  }else{
    cexOption <- 1.5
  }
  venn.plot <- VennDiagram::venn.diagram(
    x,
    filename=NULL,
    category.names=catNames,
    col="black",
    fill=colors,
    alpha=0.50,
    cex=cexOption,
    cat.col=colors,
    cat.cex=1.5,
    fontfamily="sans",
    cat.fontfamily="sans",
    margin=0.05
  )
  grid::grid.draw(venn.plot)
  invisible(file.remove(list.files(pattern="^VennDiagram.+log$")))
}
