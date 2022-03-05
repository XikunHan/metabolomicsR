
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
#' @param plot type of plot, current support `boxplot` and `betweenstats`.
#' @param x The x-axis coordinate.
#' @param feature_name A vector of selected metabolites to plot. If NULL, will randomly select 16 (default) metabolites to plot.
#' @param color A column in `@sampleData` to show the color of points.
#' @param shape A column in `@sampleData` to show the shape of points.
#' @param fill A column in `@sampleData` to show the `fill` for histogram.
#' @param random_select An integer, number of randomly selected metabolites to plot.
#' @param size Point size.
#' @param n_row Number of rows of subfigures for `betweenstats`
#' @param n_col Number of columns of subfigures for `betweenstats`
#' @param ylab Column name to annotate the y-axis in `betweenstats` (eg. "BIOCHEMICAL"), default column: "featureID".
#' @param height Height of the figure.
#' @param width Width of the figure.
#' @param save_to_file Path to save the figure.
#' @export
#' @examples
#'\dontrun{
#' p <- plot_Metabolite(df_m_PCA, plot = "boxplot")
#' p
#'}
#'
plot_Metabolite <- function(object, plot = "boxplot", x = "NEG", feature_name = NULL, color = "NEG", shape = "NEG", fill = "NEG", random_select = 16, size = 0.6, n_row = 1, n_col = 1, ylab = "featureID", height = 10, width = 10, save_to_file = NULL) {

  if(is.null(feature_name)) {
    df_select <- object@assayData[, c(1, sample(2:NCOL(object@assayData), random_select, replace = FALSE)), with = FALSE]

  } else {
    stopifnot(all(feature_name %in% names(object@assayData)))
    df_select <- object@assayData[, c(object@sampleID, feature_name), with = FALSE]
  }

  df_select <- reshape2::melt(df_select, id = object@sampleID)

  df <- merge(object@sampleData, df_select, by = object@sampleID, sort = FALSE, all = TRUE, suffixes = c("", "_"))

  check_pkg("RColorBrewer")

  if( plot == "boxplot") {
    if(! color %in% names(df)) color <- x
    if(! shape %in% names(df)) shape <- x
    
    df[, (color) := factor(get(color))]
    df[, (shape) := factor(get(shape))]
    
    p <- ggplot(data = df, aes_string(x = x, y = "value", color = df[, get(color)])) +
      geom_boxplot() +
      facet_wrap(~variable, scales = "free_y") +
      scale_color_manual(name = color, values = RColorBrewer::brewer.pal(9, "Set1")) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    if(!is.null(save_to_file)) {
      ggsave(save_to_file, p, height = height, width = width)
    }
    return(p)
  }


  if( plot == "histogram") {
    if(! color %in% names(df)) {
      warning("Color not specified")
      df$color_ <- 1
      color <- "color_"
    }
    if(! fill %in% names(df)) fill <- color
    df[, (color) := factor(get(color))]
    df[, (fill) := factor(get(fill))]

    p <- ggplot(data = df, aes_string(x = "value", color = df[, get(color)], fill = df[, get(fill)])) +
      geom_histogram(alpha=0.5) +
      facet_wrap(~variable, scales = "free_y") +
      scale_color_manual(name = color, values = RColorBrewer::brewer.pal(8, "Dark2")) +
      scale_fill_manual(name = fill, values = RColorBrewer::brewer.pal(8, "Dark2")) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    if(!is.null(save_to_file)) {
      ggsave(save_to_file, p, height = height, width = width)
    }
    return(p)
  }

  if(plot == "betweenstats") {
    check_pkg("ggstatsplot")
    v_n <- length(unique(df$variable))
    if(n_row * n_col <  v_n) {
      warning(paste0(n_row, " row and ", n_col, " column for ", v_n, " variables. To increase!"))
      n_col <- ceiling(v_n/n_row)
    }
    p_list <- list()
    if("x_" %in% names(df)) {
      warning("Overwritten `x_` column.")
      df$x <- NULL
    }
    setnames(df, x, "x_") # scope issue

    for(i in seq_along(1L:v_n)) {
      i_variable <- unique(df$variable)[i]
      p <- ggstatsplot::ggbetweenstats(
        data = df[variable == i_variable],
        x = x_,
        y = value,
        pairwise.comparisons = TRUE,
        ylab = object@featureData[featureID == i_variable, get(ylab)]
      )
      p_list[[i]] <- p
    }
    p <- ggstatsplot::combine_plots(
      p_list,
      plotgrid.args = list(nrow = n_row, ncol = n_col)
    )
    if(!is.null(save_to_file)) {
      cowplot::save_plot(save_to_file, p, base_height = height, base_width = width)
    }
    return(p)
  }
  warning("Unkown plot type!")
  return(0)
}



#' volcano plot for regression results
#'
#' @param fit regression summary results.
#' @param x The x-axis column, eg. effect size.
#' @param y The y-axis column, eg. p value.
#' @param p.value_log10 whether to transforme p.value by -log10.
#' @param color A column in fit to show different point colors. Set as NULL to turn off the color argument.
#' @param label A column in fit to label points.
#' @param highlight A column in fit to show the points to highlight. Values as 1 are hightlighted.
#' @param x_lab labels for x-axis.
#' @param y_lab labels for y-axis.
#' @export
#' @examples
#'\dontrun{
#' p <- plot_volcano(fit_lm, color = NULL)
#' p
#'}
#'

plot_volcano <- function(fit, x = "estimate", y = "p.value", p.value_log10 = TRUE,
                    color = "outcome", label = "term", highlight = "significant",
                    x_lab = "Effect size", y_lab = "-log10(P value)"

) {
  stopifnot(is.data.frame(fit))
  fit <- as.data.table(fit)

  if(p.value_log10) {
    fit[, p.value_log10 := -log10(get(y))]
  } else {
    fit[, p.value_log10 :=  get(y)]
  }

  fit$highlight_ <- NA_integer_
  if(! highlight %in% names(fit)) {
    warning(paste0(highlight, " column not in regression data. Bonferroni correction method will be used."))
    fit[, highlight_ := as.integer(get(y) < 0.05/NROW(fit))]
  } else {
    fit[, highlight_ := get(highlight)]
  }

  if(! is.null(color)) {
    fit[, (color) := factor(get(color))]

    p <- ggplot(data= fit, aes_string(x = x, y = "p.value_log10", color = color)) +
      geom_point() +
      ggrepel::geom_label_repel(data = fit[highlight_ == 1], aes(label = fit[highlight_ == 1, get(label)]), max.overlaps = 20, size = 2) +
      theme_minimal()  +
      labs(x = x_lab, y = y_lab)
  } else {
    p <- ggplot(data= fit, aes_string(x = x, y = "p.value_log10")) +
      geom_point() +
      ggrepel::geom_label_repel(data = fit[highlight_ == 1], aes(label = fit[highlight_ == 1, get(label)]), max.overlaps = 20, size = 2) +
      theme_minimal()  +
      labs(x = x_lab, y = y_lab)
  }
  return(p)
}


