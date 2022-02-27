
#' plot PCA
#'
#' Plot first two principal components.
#'
#' @param object A Metabolite object.
#' @param color A column in `@sampleData` to show the color of points.
#' @param shape A column in `@sampleData` to show the shape of points.
#' @param size Point size.
#' @export
plot_PCA <- function(object, color = "NEG", shape = "NEG", size = 1.5) {

  if(! all(c("PC1", "PC2") %in%  names(object@sampleData))) {
    object <- run_PCA(object)
  }

  df <- object@sampleData
  df[, (color) := factor(get(color))]
  df[, (shape) := factor(get(shape))]

  p <- ggplot(df, aes_string("PC1", "PC2")) +
    geom_point(aes(color = factor(df[, get(color)]), shape = df[, get(shape)]), size = size) +
    scale_color_discrete(name = color) +
    scale_shape_manual(name = shape, values = point_shape_plate)
  return(p)
}



#' Plot UMAP
#'
#' Plot manifold approximation and projection (UMAP). See more details in \code{\link[M3C]{umap}}.
#'
#' @param object A Metabolite object.
#' @param color A column in `@sampleData` to show the color of points.
#' @param shape A column in `@sampleData` to show the shape of points.
#' @param size Point size.
#' @export
#' @examples
#' \dontrun{
#' p<- plot_UMAP(df_2017_QC_norm_PCA, color = "NEG", shape = "QCSample")
#' p
#' }
plot_UMAP <- function(object, color = "NEG", shape = "NEG", size = 1.5) {

  check_pkg("M3C")
  df <- M3C::umap(t(object@assayData[,-1]))
  df <- cbind(object@assayData[, 1], df$data)

  df_sample <- object@sampleData
  df_sample[, (color) := factor(get(color))]
  df_sample[, (shape) := factor(get(shape))]

  df <- merge(df, df_sample, by = object@sampleID)

  p <- ggplot(df, aes_string("X1", "X2")) +
    geom_point(aes(color = factor(df[, get(color)]), shape = df[, get(shape)]), size = size) +
    scale_color_discrete(name = color) +
    scale_shape_manual(name = shape, values = point_shape_plate) +
    labs(x = "UMAP-1", y = "UMAP-2")
  return(p)
}




#' plot tSNE
#'
#' Plot t-distributed stochastic neighbor embedding. See more details in \code{\link[M3C]{tsne}}.
#'
#' @param object A Metabolite object.
#' @param color A column in `@sampleData` to show the color of points.
#' @param shape A column in `@sampleData` to show the shape of points.
#' @param size Point size.
#' @export
#' @examples
#' \dontrun{
#' p<- plot_tsne(df_2017_QC_norm_PCA, color = "NEG", shape = "QCSample")
#' p
#'
#' }
plot_tsne <- function(object, color = "NEG", shape = "NEG", size = 1.5) {

  check_pkg("M3C")
  df <- M3C::tsne(t(object@assayData[,-1]))
  df <- cbind(object@assayData[, 1], df$data)

  df_sample <- object@sampleData
  df_sample[, (color) := factor(get(color))]
  df_sample[, (shape) := factor(get(shape))]

  df <- merge(df, df_sample, by = object@sampleID)

  p <- ggplot(df, aes_string("X1", "X2")) +
    geom_point(aes(color = factor(df[, get(color)]), shape = df[, get(shape)]), size = size) +
    scale_color_discrete(name = color) +
    scale_shape_manual(name = shape, values = point_shape_plate) +
    labs(x = "tSNE-1", y = "tSNE-2")
  return(p)
}



#' injection order scatterplot
#'
#' Injection order scatterplot. The `@sampleData` should be sorted by injection order, with a new column `ID` from 1 to N.
#'
#' @param object A Metabolite object.
#' @param color A column in `@sampleData` to show the color of points.
#' @param shape A column in `@sampleData` to show the shape of points.
#' @param size Point size.
#' @param ID_order Injection ID order in the `@sampleData`.
#' @param feature_name A vector of selected metabolites to plot. If NULL, will randomly select 16 (default) metabolites to plot.
#' @param random_select An integer, number of randomly selected metabolites to plot.
#' @export
#' @examples
#'\dontrun{
#' p <- plot_injection_order(df_m_PCA, color = "QC_sample")
#' p
#'
#' p <- plot_injection_order(df_m_PCA, color = "QC_sample", feature_name = "X563")
#' p
#'}
#'
plot_injection_order <- function(object, color = "NEG", shape = "NEG", size = 0.6, ID_order = "ID_injection_order", feature_name = NULL, random_select = 16) {

  if(is.null(feature_name)) {
    df_select <- object@assayData[, c(1, sample(2:NCOL(object@assayData), random_select, replace = FALSE)), with = FALSE]
    names(df_select)
  } else {
    df_select <- object@assayData[, c(object@sampleID, feature_name), with = FALSE]
  }

  df_select <- reshape2::melt(df_select, id = object@sampleID)

  df <- merge(object@sampleData, df_select, by = object@sampleID, sort = FALSE, all = TRUE, suffixes = c("", "_"))

  df[, (color) := factor(get(color))]
  df[, (shape) := factor(get(shape))]

  check_pkg("RColorBrewer")
  p <-   ggplot(data = df, aes_string(ID_order, "value")) +
    geom_point(aes( color = df[, get(color)], shape = df[, get(shape)]), size = size) +
    facet_wrap(~variable, scales = "free_y") +
    scale_shape_manual(name = shape, values = point_shape_plate) +
    scale_color_manual(name = color, values = RColorBrewer::brewer.pal(9, "Set1"))
  return(p)
}




#' plot a Metabolite object
#'
#' Plot a Metabolite object including boxplot (more to add.).
#'
#' @param object A Metabolite object.
#' @param plot type of plot
#' @param x The x-axis coordinate.
#' @param y The y-axis coordinate.
#' @param color A column in `@sampleData` to show the color of points.
#' @param shape A column in `@sampleData` to show the shape of points.
#' @param feature_name A vector of selected metabolites to plot. If NULL, will randomly select 16 (default) metabolites to plot.
#' @param random_select An integer, number of randomly selected metabolites to plot.
#' @param size Point size.
#' @export
#' @examples
#'\dontrun{
#' p <- plot_Metabolite(df_m_PCA, plot = "boxplot")
#' p
#'}
#'
plot_Metabolite <- function(object, plot = "boxplot", x = "NEG", y = "value", color = "NEG", shape = "NEG", feature_name = NULL, random_select = 16, size = 0.6) {

  if(is.null(feature_name)) {
    df_select <- object@assayData[, c(1, sample(2:NCOL(object@assayData), random_select, replace = FALSE)), with = FALSE]
    names(df_select)
  } else {
    df_select <- object@assayData[, c(object@sampleID, feature_name), with = FALSE]
  }


  df_select <- reshape2::melt(df_select, id = object@sampleID)

  df <- merge(object@sampleData, df_select, by = object@sampleID, sort = FALSE, all = TRUE, suffixes = c("", "_"))

  df[, (color) := factor(get(color))]
  df[, (shape) := factor(get(shape))]

  check_pkg("RColorBrewer")

  if( plot == "boxplot") {
   p <- ggplot(data = df, aes_string(x = x, y = y, color = df[, get(color)])) +
     # geom_violin() +
      geom_boxplot() +
     facet_wrap(~variable, scales = "free_y") +
     scale_color_manual(name = color, values = RColorBrewer::brewer.pal(9, "Set1")) +
     theme(axis.text.x = element_text(angle = 60, hjust = 1))
   return(p)
  }
}




