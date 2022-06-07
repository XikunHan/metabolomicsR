
#' plot PCA
#'
#' Plot first two principal components.
#'
#' @param object A Metabolite object.
#' @param color A column in `@sampleData` to show the color of points.
#' @param shape A column in `@sampleData` to show the shape of points.
#' @param size Point size.
#' @export
#' @return PCA plot. 
plot_PCA <- function(object, color = "NEG", shape = "NEG", size = 1.5) {

  if(! all(c("PC1", "PC2") %in%  names(object@sampleData))) {
    object <- run_PCA(object)
  }

  df <- object@sampleData
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
#' @return UMAP plot.
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
#' @return tSNE plot.
#'
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
#' @return A scatterplot.
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
#' @return A boxplot of a Metabolite object
plot_Metabolite <- function(object, plot = "boxplot", x = "NEG", feature_name = NULL, color = "NEG", shape = "NEG", fill = "NEG", random_select = 16, size = 0.6, n_row = 1, n_col = 1, ylab = "featureID", height = 10, width = 10, save_to_file = NULL) {

  featureID <- NULL
  
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

    value <- variable <- x_ <- NULL
    
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
#' @return A volcano plot.
plot_volcano <- function(fit, x = "estimate", y = "p.value", p.value_log10 = TRUE,
                    color = "outcome", label = "term", highlight = "significant",
                    x_lab = "Effect size", y_lab = "-log10(P value)"

) {
  highlight_ <- NULL
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




#' ROC
#'
#' Plot Receiver Operating Characteristic (ROC) curve for metabolites with or without covariates
#' @param object A Metabolite object.
#' @param y A column name for the disease (0, 1)
#' @param x One variable name (if x is provided, model_a and model_b should be NULL or vice versa).  
#' @param model_a Column names for model a (one or more covariates, as the first model).
#' @param model_b Column names for model b (one or more covariates, as the second model).
#' @param lab Title (eg. "BIOCHEMICAL"), default value is x.
#' @param text_x Annotation position of text in x-axis. 
#' @param text_y Annotation position of text in y-axis.
#' @export
#' @return ROC. 
plot_ROC <- function(object = NULL, y = NULL, x = NULL, model_a = NULL, model_b = NULL, lab = NULL, text_x = 0.6, text_y = 0.3) {
  
  D <- M <- model_A <- name <- y_magic <- NULL
  
  df <- merge(object@sampleData, object@assayData, by = object@sampleID)
  
  if(length(x) == 0 ) {
    stopifnot(!any(is.null(model_a), is.null(model_b)))
  }
  
  # NULL is removed in c()
  stopifnot(all(c(x, y,model_a, model_b) %in% names(df)))
  
  if(is.null(lab)) lab <- x[1]
  
  df <- df[, c(x, y,model_a, model_b), with = FALSE]
  df <- df[complete.cases(df), ]
  
  df$y_magic <- df[, get(y)]
  
  check_pkg("pROC")
  check_pkg("RColorBrewer")
  
  f_auc_ci <- function(formula = NULL, data = NULL) {
    fit <- pROC::roc(formula = formula, data = data , ci = TRUE)
    auc_ci <- as.numeric(fit$ci)
    auc_ci <- paste0(sprintf("%1.2f",auc_ci[2]),
                     " [",sprintf("%1.2f",auc_ci[1]),
                     ", ",sprintf("%1.2f",auc_ci[3]),"]")
    return(list(fit = fit, auc_ci = auc_ci))
  }
  
  if(length(x) != 0) {
    fit1 <- glm(as.formula(paste(y , " ~ ", paste0(x, collapse = " + "))), 
                data = df,family = binomial())
    df$model_A <- predict(fit1, type=c("response"))
    auc_ci <- f_auc_ci(as.formula(paste(y , " ~ ", "model_A")), data = df)
    df$y_magic <- as.integer(df$y_magic)
    
    p <- ggplot(df, aes(m = model_A, d = y_magic)) + 
      geom_roc(labels=FALSE, n.cuts = 0) +
      geom_abline(intercept = 0, slope = 1,linetype=4) +
      annotate("text",x= text_x, y= text_y, label = auc_ci$auc_ci,
               parse = FALSE,colour = "red", size =5) +
      theme(legend.justification=c(0,0),
            plot.title = element_text(hjust = 0.5),
            legend.position=c(0.18,0.02),
            legend.title = element_blank(),
            legend.background = element_rect(fill=alpha("blue", 0)),
            panel.background = element_blank(),
            axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black"),
            text=element_text(face="bold", size=12)) +
      labs(x="1 - Specificity", y = "Sensitivity",title= lab)
    return(p)
    
  } else {
    # RColorBrewer::display.brewer.all()
    v_color <- RColorBrewer::brewer.pal(9,"Set1")[c(2,5)]
    fit1 <- glm(as.formula(paste(y , " ~ ", paste0(model_a, collapse = " + "))), 
                data = df,family = binomial())
    df$model_A <- predict(fit1, type=c("response"))
    
    fit2 <- glm(as.formula(paste(y , " ~ ", paste0(model_b, collapse = " + "))), 
                data =df,family = binomial())
    df$model_B <- predict(fit2, type=c("response"))
    
    auc_ci1 <- f_auc_ci(as.formula(paste(y , " ~ ", "model_A")), data = df)
    auc_ci2 <- f_auc_ci(as.formula(paste(y , " ~ ", "model_B")), data = df)
    
    fit_p.value <- pROC::roc.test(auc_ci1$fit, auc_ci2$fit)$p.value
    df_roc <- melt_roc(df, y, c("model_A", "model_B"))
    df_roc$D <- as.integer(df_roc$D) # seems a bug in plotROC, need to convert to integer for html
    
    model_A_text <- paste0("Model A: ", auc_ci1$auc_ci)
    model_B_text <- paste0("Model B: ", auc_ci2$auc_ci)
    
    p <- ggplot(df_roc, aes(m = M, d = D, color = name)) + 
      stat_roc(labels=FALSE, n.cuts=0) +
      scale_color_manual(values = v_color) +
      geom_abline(intercept = 0, slope = 1,linetype=4) +
      annotate("text", x= text_x, y= text_y + 0.05,  label = paste0(model_A_text, "\n"), size =4, colour = v_color[1]) +
      annotate("text", x= text_x, y= text_y,  label = paste0(model_B_text, "\n"), size =4, colour = v_color[2]) +
      annotate("text", x= text_x, y= text_y - 0.05,  label = paste0("(P-diff: ", sprintf("%1.3f",fit_p.value), ")\n"), size =4) +
      theme(legend.justification=c(0,0),
            legend.position= "none",
            legend.title = element_blank(),
            legend.background = element_rect(fill=alpha("blue", 0)),
            panel.background = element_blank(),
            axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black"),
            text=element_text(face="bold", size=12)) +
      labs(x="1 - Specificity", y = "Sensitivity",title= lab) 
    return(p)
  }
}






#' quality control visualization
#'
#' This function will plot QC results
#'
#' @param object A Metabolite object.
#' @param nSD Define the N times of the SD as outliers.
#' @export
#' @return QC metrics
#'
plot_QC <- function(object, nSD = 5) {
  df <- calculate_column_constant(object@assayData[, -1], verbose = FALSE)
  df <- data.table(ID = names(object@assayData[, -1]), x = df == 0 | is.na(df))
  column_constant <- df
  p_constant <- ggplot(data= df, aes(x=x)) +
    geom_bar(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
    hrbrthemes::theme_ipsum() +
    labs(x = "", title = "Constant columns") +
    theme(title = element_text(size = 12))
  
  df <- column_missing_rate(object)[-1]
  df <- data.table(ID = names(object@assayData[, -1]), x = df)
  column_missing_rate <- df
  p_column_missing_rate <- ggplot(data = df, aes(x= x)) +
    geom_histogram(bins = 30, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    hrbrthemes::theme_ipsum() +
    labs(x = "", title = "Column missing rate") +
    theme(title = element_text(size = 12))
  
  
  df <- row_missing_rate(object)
  df <- data.table(ID = object@assayData$sampleID, x = df)
  row_missing_rate <- df
  p_row_missing_rate <- ggplot(data = df, aes(x= x)) +
    geom_histogram(bins = 30, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    hrbrthemes::theme_ipsum() +
    labs(x = "", title = "Row missing rate") +
    theme(title = element_text(size = 12))

  
  df <- outlier_rate(object, nSD = nSD)
  df <- data.table(ID = names(object@assayData[, -1]), x = df)
  outlier_rate <- df
  p_outlier_rate <- ggplot(data = df, aes(x= x)) +
    geom_histogram(bins = 30, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    hrbrthemes::theme_ipsum() +
    labs(x = "", title = "Outlier rate") +
    theme(title = element_text(size = 12))
  
  p_list <- list(p_constant, 
                 p_column_missing_rate, 
                 p_row_missing_rate,
                 p_outlier_rate)
  
  
  p <- ggstatsplot::combine_plots(
    p_list,
    plotgrid.args = list(nrow = 2, ncol = 2)
  )
  
  QC_metrics <- list(column_constant = column_constant, 
                     column_missing_rate = column_missing_rate,
                     row_missing_rate = row_missing_rate,
                     outlier_rate = outlier_rate, 
                     p = p
                     )
  return(QC_metrics)
}


