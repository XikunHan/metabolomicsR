


############### QC #######################

#' @rdname column_missing_rate
#' @export
#' @examples
#' \dontrun{
#' # for a data.frame or data.table
#' v <- column_missing_rate(object = df)
#'
#' # to skip the first column (eg. ID)
#' v <- column_missing_rate(object = df[, -1])
#' }
#'
column_missing_rate.default <- function(object) {
  r <- apply(object, 2, function(x) sum(is.na(x)) / length(x))
  names(r) <- names(object)
  return(r)
}


#' @rdname column_missing_rate
#' @export
#' @examples
#' \dontrun{
#'
#' # for a Metabolite object
#' v <- column_missing_rate(object)
#'
#' }
#'
column_missing_rate.Metabolite <- function(object) {
  object <- object@assayData
  r <- apply(object, 2, function(x) sum(is.na(x)) / length(x))
  names(r) <- names(object)
  return(r)
}




#' @rdname filter_column_missing_rate
#' @export
#' @examples
#' \dontrun{
#'
#' d <- filter_column_missing_rate(object)
#'
#' }
#'
filter_column_missing_rate.default <- function(object, threshold = 0.5, verbose = TRUE) {
  object <- as.data.table(object)
  r <- column_missing_rate(object)

  if(verbose) {
    cat(paste0("\n Number of columns with a missing rate >= ", threshold, " : n = ", sum(r >= threshold),  "\n"))
  }
  r <- r[r < threshold]
  return(object[, names(r), with = FALSE])
}




#' @rdname filter_column_missing_rate
#' @export
#'
filter_column_missing_rate.Metabolite <- function(object, threshold = 0.5, verbose = TRUE) {

  ncol_old <- NCOL(object@assayData)
  object@assayData <- filter_column_missing_rate(object@assayData, threshold = threshold, verbose = verbose)
  object@logs <- paste0(object@logs,
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Filter data with a missing rate >= ", threshold,
                        ", removed ", ncol_old -  NCOL(object@assayData), " features. \n")
  if(ncol_old -  NCOL(object@assayData) == 0) {
    return(object)
  }
  return(update_Metabolite(object))
}



#' @rdname row_missing_rate
#' @export
#' @examples
#' \dontrun{
#'
#' # for a data.frame or data.table
#' v <- row_missing_rate(object = df)
#'
#' # to skip the first column (eg. ID)
#' v <- row_missing_rate(object = df[, -1])
#' }
#'
row_missing_rate.default <- function(object) {
  r <- apply(object, 1, function(x) sum(is.na(x)) / length(x))
  names(r) <- unlist(object[, 1])
  return(r)
}


#' @rdname row_missing_rate
#' @export
#' @examples
#' \dontrun{
#'
#' # for a Metabolite object
#' v <- row_missing_rate(object)
#'
#' }
#'
row_missing_rate.Metabolite <- function(object) {
  object <- object@assayData
  # here to skip the first column (sample ID)
  r <- apply(object[, -1], 1, function(x) sum(is.na(x)) / length(x))
  names(r) <- unlist(object[, 1])
  return(r)
}




#' @rdname filter_row_missing_rate
#' @export
#' @examples
#' \dontrun{
#'
#' d <- filter_row_missing_rate(object)
#'
#' }
#'
filter_row_missing_rate.default <- function(object, threshold = 0.5, verbose = TRUE) {
  stopifnot(is.numeric(threshold) & threshold <=1 & threshold >= 0)
  object <- as.data.table(object)
  r <- row_missing_rate(object)
  if(verbose) {
    cat(paste0("\n Number of rows with a missing rate >= ", threshold, " : n = ", sum(r >= threshold),  "\n"))
  }

  if(sum(r >= threshold) == 0) {
    return(object)
  }
  return(object[r < threshold, ])
}




#' @rdname filter_row_missing_rate
#' @export
#'
filter_row_missing_rate.Metabolite <- function(object, threshold = 0.5, verbose = TRUE) {
  stopifnot(is.numeric(threshold) & threshold <=1 & threshold >= 0)
  nrow_old <- NROW(object@assayData)
  object@assayData <- filter_row_missing_rate(object@assayData, threshold = threshold, verbose = verbose)
  object@logs <- paste0(object@logs,
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Filter data with a missing rate >= ", threshold,
                        ", removed ", nrow_old -  NROW(object@assayData), " samples. \n")
  if(nrow_old -  NROW(object@assayData) == 0) {
    return(object)
  }
  return(update_Metabolite(object))
}





#' @rdname filter_column_constant
#' @export
#' @examples
#' \dontrun{
#'
#' # for a data.frame or data.table
#' v <- filter_column_missing_rate(object = df)
#'
#' # if skip the first column (eg. ID)
#' v <- filter_column_missing_rate(object = df[, -1])
#'
#' }
#'
filter_column_constant.default <- function(object, verbose = TRUE) {
  object <- as.data.table(object)
  r <- apply(object, 2, function(x) sd(x, na.rm = TRUE))

  if(verbose) {
    cat(paste0("\nConstant columns n = ", sum(r == 0 | is.na(r)), "\n"))
  }
  r <- r[! (r == 0 | is.na(r))]

  return(object[, names(r), with = FALSE])
}



#' @rdname filter_column_constant
#' @export
#'
filter_column_constant.Metabolite <- function(object, verbose = TRUE) {
  ncol_old <- NCOL(object@assayData)
  # skip the first column (sample ID). and then merge
  object@assayData <- cbind(object@assayData[, 1], filter_column_constant(object = object@assayData[, -1], verbose = verbose))
  object@logs <- paste0(object@logs,
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Filter data with a constant column, removed ",
                        ncol_old -  NCOL(object@assayData), " features. \n")

  if(ncol_old -  NCOL(object@assayData) == 0) {
    return(object)
  }

  return(update_Metabolite(object))
}



#' @rdname subset
#' @export
#'
subset.Metabolite <- function(object, subset, select) {

  if (!missing(subset)) {
    r <-  {
      e <- substitute(subset)
      r <- eval(e, object@sampleData, parent.frame())
      if (!is.logical(r))
        stop("'subset' must be logical")
      r & !is.na(r)
    }

    df_sample <- object@sampleData[r, ]
    object <- update_Metabolite(object, dataset = df_sample[, get(object@sampleID)],  action = "keep_sample")
  }

  if (!missing(select)) {
    vars <- {
      nl <- as.list(seq_along(object@assayData))
      names(nl) <- names(object@assayData)
      eval(substitute(select), nl, parent.frame())
    }
    vars <- names(object@assayData)[vars]
    object <- update_Metabolite(object, dataset = vars,  action = "keep_feature")
  }

  return(object)
}






#' @rdname replace_outlier
#' @export
replace_outlier.default <- function(object, method = "winsorize", nSD = 5) {
  if(is.vector(object)) {
    x <- as.numeric(object)

    stopifnot(method %in% c('as_NA', "winsorize"))

    v_mean <- mean(x, na.rm = TRUE)
    v_sd <- sd(x, na.rm = TRUE)
    r <- x > (v_mean + nSD*v_sd) | x < (v_mean - nSD*v_sd)

    if(method == "as_NA") {
      x <- ifelse(r, NA, x)

    } else if (method == "winsorize") {
      # extreme big and small values replaced with the max and min values.
      r_high <- x > (v_mean + nSD*v_sd)
      r_low <- x < (v_mean - nSD*v_sd)

      v_max <- max(x[!r], na.rm = TRUE)
      v_min <- min(x[!r], na.rm = TRUE)

      x <- ifelse(r_high, v_max, x)
      x <- ifelse(r_low, v_min, x)
    }
    return(x)
  } else stop("Unknown object, convert to vector?", call. = FALSE)

}



#' @rdname replace_outlier
#' @export
#' @examples
#'
#'
#' \dontrun{
#' d <- replace_outlier(object, method = "winsorize", nSD = 5)
#' }
#'
replace_outlier.data.frame <- function(object, method = "winsorize", nSD = 5) {
  object <- apply(object, 2, replace_outlier, method = method, nSD = nSD)
  return(object)
}



#' @rdname replace_outlier
#' @export
#' @examples
#'
#'
#' \dontrun{
#' d <- replace_outlier(object, method = "winsorize", nSD = 5)
#' }
#'
replace_outlier.Metabolite <- function(object, method = "winsorize", nSD = 5) {


  data <- replace_outlier(object@assayData[, -1], method = method, nSD = nSD)
  object@assayData <- cbind(object@assayData[, 1],data)
  object@logs <- paste0(object@logs,
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Replace outliers using `", method, "` method ",nSD," SDs. \n")
  return(object)
}




#' is outlier
#'
#' @param object An object, a vector.
#' @param nSD N times of the SD as outliers.
#' Return TRUE or FALSE for a vector.
#' @export
#' @examples
#'
#'\dontrun{
#' v <- is_outlier(x)
#' }
#'
is_outlier <- function(object, nSD = 5) {
  object <- as.numeric(object)
  v_mean <- mean(object, na.rm = TRUE)
  v_sd <- sd(object, na.rm = TRUE)
  r <- object > (v_mean + nSD*v_sd) | object < (v_mean - nSD*v_sd)
  return(r)
}


#' @rdname outlier_rate
#' @export
#' @examples
#'
#'\dontrun{
#' v <- outlier_rate(x)
#' }
#'
outlier_rate.default <- function(object, nSD = 5) {
  r <- mean(is_outlier(object, nSD = nSD), na.rm = TRUE)
  return(r)
}


#' @rdname outlier_rate
#' @export
outlier_rate.data.frame <- function(object, nSD = 5) {
  r <- apply(object, 2, outlier_rate, nSD = nSD)
  names(r) <- names(object)
  return(r)
}


#' @rdname outlier_rate
#' @export
#' @examples
#' \dontrun{
#'
#' # for a Metabolite object
#' v <- outlier_rate(object)
#'
#' }
#'
outlier_rate.Metabolite <- function(object, nSD = 5) {
  object <- object@assayData[, -1]
  r <- outlier_rate(object, nSD = nSD)
  names(r) <- names(object)
  return(r)
}



#' quality control pipeline
#'
#' This function will run QC steps on a Metabolite object
#'
#' @param object An object, data.frame, data.table or Metabolite.
#' @param filter_column_constant A logical value, whether to filter columns (features) with a constant value.
#' @param filter_column_missing_rate_threshold A numeric threshold to filter columns (features) below a missing rate, default: 0.5. Other values: 0.2, 0.8. If NULL, then skip this step.
#' @param filter_row_missing_rate_threshold A numeric threshold to filter rows (samples) below a missing rate. Default: NULL, to skip this step. Other values: 0.5, 0.2, 0.8.
#' @param replace_outlier_method Method to replace outlier value, see \code{\link{replace_outlier}}.
#' @param nSD Define the N times of the SD as outliers.
#' @param impute_method Imputation method, the default method is half the minimum value (`half-min`) of the metabolite. Currently support 'half-min', "median", "mean", "zero".
#' @param verbose print log information.
#' @export
#'
QC_pipeline <- function(object,
                        filter_column_constant = TRUE,
                        filter_column_missing_rate_threshold = 0.5,
                        filter_row_missing_rate_threshold = NULL,
                        replace_outlier_method = NULL,
                        nSD = 5,
                        impute_method = "half-min",
                        verbose = TRUE
) {

  object@logs <- paste0(object@logs, "\n",
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Run QC pipeline.\n")

  if(filter_column_constant) {
    object <- filter_column_constant(object, verbose = verbose)
  }

  if(!is.null(filter_column_missing_rate_threshold)) {
    object <- filter_column_missing_rate(object, threshold = filter_column_missing_rate_threshold, verbose = verbose)
  }

  if(!is.null(filter_row_missing_rate_threshold)) {
    object <- filter_row_missing_rate(object, threshold = filter_row_missing_rate_threshold, verbose = verbose)
  }

  if(!is.null(replace_outlier_method)) {
    object <- replace_outlier(object, method = replace_outlier_method, nSD = nSD)
  }

  if(!is.null(impute_method)) {
    object <- impute(object, method = impute_method)
  }
  return(object)
}


#' RSD
#'
#' calculate RDS (%)
#'
#' @param x A vector
#' @export
#' @examples
#'\dontrun{
#' v <- RSD(x)
#'}
RSD <- function(x) {
  v_std <- sd(x, na.rm = TRUE)
  v_mean <- mean(x, na.rm = TRUE)
  res <- (v_std/v_mean) * 100
  return(res)
}




#' distinguish genuine untargeted metabolic features without QC samples
#'
#' The makefeature function from genuMet uses a Metabolite object as input. genuMet is an R package used  distinguish genuine untargeted metabolic features without quality control samples.
#'
#' @param object A Metabolite object.
#' @param wsize Window size.
#' @param ssize Slide size.
#' @param defswitch Definition of a switch.
#' @export
#' @references <https://github.com/liucaomics/genuMet>
#' @examples
#'\dontrun{
#' v <- genuMet_makefeature(df)
#'}
#'
genuMet_makefeature <- function(object, wsize=100, ssize= 0.5, defswitch=0.2) {
  if (! requireNamespace("genuMet", quietly = TRUE)) {
    stop(paste0("Please install genuMet package first: ` (\"xyomics/genuMet\") `."), call. = FALSE)
  }
  df <- object@assayData
  df <- t(df[, -1])
  colnames(df) <- object@assayData[, get(object@sampleID)]
  df <- as.data.frame(df)
  metf <- genuMet::makefeature(data=df, wsize=wsize, ssize= ssize, defswitch = defswitch)
  return(metf)
}



############### transformation #######################


#' pareto scale transformation
#'
#' pareto scale transformation
#'
#' @param x A vector
#' @export
#' @examples
#'
#'\dontrun{
#' v <- paretoscale(x)
#'}
#'
pareto_scale <- function(x) {
  v_mean <- mean(x, na.rm = TRUE)
  v_sd <- sd(x, na.rm = TRUE)
  res <- (x - v_mean)/sqrt(v_sd)
  return(res)
}




#' rank-based inverse normal transformation
#'
#' rank-based inverse normal transformation for a metabolite.
#'
#' @param x A vector
#' @export
#' @examples
#'
#'\dontrun{
#' v <- inverse_rank_transform(x)
#'}
#'
inverse_rank_transform <- function(x) {
  stopifnot(is.vector(x))
  transformed <- qnorm((rank(x ,na.last="keep") -0.5) /sum(!is.na(x)))
  return(transformed)
}



#' apply transformation to a Metabolite object
#'
#' Apply transformation to Metabolite  object
#'
#' @param object A Metabolite object.
#' @param method Transform method, eg. "log", "pareto_scale", "scale", "inverse_rank_transform". A User defined method is also supported.
#' @export
#' @examples
#'
#'\dontrun{
#' d <- transformation(x)
#'}
#'
transformation <- function(object, method = "log") {
  stopifnot(inherits(object, "Metabolite"))
  if(!method %in% c("log", "pareto_scale", "scale", "inverse_rank_transform")) {
    warnings(paste0(paste0("A user defined function: ", method, ". \n", c("log", "pareto_scale", "scale", "inverse_rank_transform"), collapse = ", "), " are currently provided in the `metabolomicsR` package."), call. = FALSE)
  }
  object@assayData <- cbind(object@assayData[, 1],
                            apply(object@assayData[, -1], 2, function(x) do.call(method, list(x))))
  object@logs <- paste0(object@logs,
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Transformation using `",method, "` method. \n")
  return(object)
}







#' bridge different data sets based on conversion factors
#'
#' Bridge metabolite data based on a conversion factor file
#'
#' @param object A Metabolite object. In the `featureData`, `conversion_factor_ID` column should be created to match with conversion_factor_data.
#' @param conversion_factor_data A data set with columns `conversion_factor_ID` and `conversion_factor_value`.
#' @param QC_ID_pattern A character pattern to determine QC samples. Default value: "MTRX". Skip QC samples when rescale (median value is already 1).
#' @param verbose print log information.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @return A Metabolite object after multiplying by conversion factor. 
#' 
bridge <- function(object, conversion_factor_data = NULL, QC_ID_pattern = "MTRX", verbose = TRUE) {

  conversion_factor_data <- as.data.table(conversion_factor_data)

  if(!all(c("conversion_factor_ID", "conversion_factor_value") %in% names(conversion_factor_data))) {
    stop(paste0("`conversion_factor_ID`, `conversion_factor_value` do not exist in conversion_factor_data."), call. = FALSE)
  }

  stopifnot(inherits(object, "Metabolite"))

  if(! "conversion_factor_ID" %in% names(object@featureData)) {
    stop(paste0("`conversion_factor_ID` does not exist in @featureData. Please create this new column in @featureData."), call. = FALSE)
  }


  # check how many metabolites have conversion factor.
  v_metabolite <- setdiff(names(object@assayData), object@sampleID)

  v_metabolite_remove <- object@featureData[! object@featureData$conversion_factor_ID %in% conversion_factor_data$conversion_factor_ID , get(object@featureID)]

  if(verbose) {
    cat(paste0("\n Remove ", length(v_metabolite_remove), " metabolites without conversion factors: ", paste0(v_metabolite_remove[seq_len(5)], collapse = ", "),  " ... \n"))
  }

  # check sample IDs
  n_QC_sample <- sum(grepl(QC_ID_pattern, object@assayData[, get(object@sampleID)]))
  sampleIndex <- grep(QC_ID_pattern, object@assayData[, get(object@sampleID)], invert = TRUE)

  if(n_QC_sample == 0) {
    stop(paste0("No QC samples in @assayData (ID: ", object@sampleID,")"), call. = FALSE)
  }
  if(verbose) {
    cat(paste0("\n Number of QC samples n = ", n_QC_sample, "\n"))
  }

  # filter metabolite data.
  object@assayData <- object@assayData[, -(v_metabolite_remove), with = FALSE]
  object <- update_Metabolite(object)

  OutputData <- copy(object@assayData)
  InputData <- copy(object@assayData)


  v_n_col <- dim(object@assayData)[2]
  v_n_col

  pb <- txtProgressBar(min = 0, max = v_n_col, style = 3, file = stderr())


  for(j in 2L:v_n_col) {

    setTxtProgressBar(pb = pb, value = j)
    v_metab <- names(object@assayData)[j]

    InputData_each <- unlist(InputData[, j, with = FALSE])

    # use metabolite ID, get the conversion ID in @featureData
    v_conversion_factor_ID <- unlist(object@featureData[get(object@featureID) == v_metab, "conversion_factor_ID", with = FALSE])

    v_value <- unlist(conversion_factor_data[conversion_factor_data$conversion_factor_ID == v_conversion_factor_ID, "conversion_factor_value", with = FALSE])

    if(length(v_value) == 0) {
      warning(paste0("No conversion factor for  ", v_metab))
      v_value <- NA
    }

    # OutputData[, (v_metab) := object@assayData[, get(v_metab)] * v_value]
    set(OutputData, sampleIndex, j, InputData_each[sampleIndex]* v_value)

  }
  close(con = pb)

  object@assayData <- OutputData
  object@miscData[['No_conversion_factor']] <- v_metabolite_remove

  object@logs <- paste0(object@logs, "\n",
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Rescale by conversion factor.\n")
  return(object)
}



############### normalization #######################

#' QCmatrix normalization
#'
#' Normalization data by the median value of QC samples in each batch. For each metabolite, the values (eg. raw peak area data) were divided by the median value of QC samples in that batch. QC samples and metabolite batches should be specified (see parameters below).
#'
#' @param object A Metabolite object. In the feature annotation slot `feature`, a platform column should be provided for metabolite measurement platform (eg. `PLATFORM`). The values in the `PLATFORM` column (eg. `Neg`, `Polar`, `Pos Early`, and `Pos Late`) are column names in the sample annotation `sample` to determine the batches of samples.
#' @param feature_platform The column name of feature platform for metabolite measurements (eg. `PLATFORM`).
#' @param QC_ID_pattern A character pattern to determine QC samples. Default value: "MTRX".
#' @param test test the function for the first 20 columns.
#' @param verbose print log information.
#' @seealso \code{\link{batch_norm}}
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' d <- QCmatrix_norm(object = df)
#' }
#' @export
#'
#'
QCmatrix_norm <- function(object, feature_platform = "PLATFORM", QC_ID_pattern = "MTRX", test = FALSE, verbose = TRUE) {
  stopifnot(inherits(object, "Metabolite"))
  OutputData <- copy(object@assayData)
  InputData <- copy(object@assayData)

  # check platform column
  if(! feature_platform %in% names(object@featureData)) {
    stop(paste0("`", feature_platform, "` column does no exist in @featureData"), call. = FALSE)
  }

  platforms <- object@featureData[, get(feature_platform)]

  if(verbose) {
    cat("\nPlatform information in @featureData:\n")
    print(table(platforms))
    cat("\n")
  }

  if(all(toupper(unique(platforms)) %in% names(object@sampleData))) {
    if(verbose) {
      cat(paste0("\nSample size in platform ", toupper(sort(unique(platforms))[1])))
      print(table(object@sampleData[, toupper(unique(platforms)[1]), with = FALSE]))
      cat("\n")
    }
  } else {
    stop(paste0("`", paste0(unique(platforms), sep = "", collapse = " "), "` columns do no exist in @sampleData"), call. = FALSE)
  }

  # check sample IDs
  n_QC_sample <- sum(grepl(QC_ID_pattern, object@assayData[, get(object@sampleID)]))

  if(n_QC_sample == 0) {
    stop(paste0("No QC samples in @assayData (ID: ", object@sampleID,")"), call. = FALSE)
  }

  if(verbose) {
    cat(paste0("\n Number of QC samples n = ", n_QC_sample, "\n"))
  }

  v_n_col <- dim(object@assayData)[2]
  if(isTRUE(test)) {
    v_n_col <- 10
  } else if(test >=2 ) {
    v_n_col <- test
  }

  pb <- txtProgressBar(min = 0, max = v_n_col, style = 3, file = stderr())

  sample_IDs <- object@assayData[,get(object@sampleID)]

  list_QC_sample_missing_0.5 <- list()

  for(j in 2L:v_n_col) {

    setTxtProgressBar(pb = pb, value = j)
    v_metab <- names(object@assayData)[j]

    v_platform <- object@featureData[get(object@featureID) == v_metab, get(feature_platform)]
    v_platform <- toupper(v_platform)


    v_batch <- unlist(object@sampleData[, v_platform, with = FALSE])

    InputData_each <- unlist(InputData[, j, with = FALSE])

    for (id_batch in unique(v_batch[! is.na(v_batch)])) {

      Index_each <- which(v_batch == id_batch)

      QCIndex <- grep(QC_ID_pattern, sample_IDs)
      QCIndex <- QCIndex[QCIndex %in% Index_each]

      v_missing_rate <- missing_rate(InputData_each[QCIndex])

      # set as missing if more than half QC are missing.
      # print(cat(v_metab, ": ", v_missing_rate, ", id_batch = ", id_batch, ", v_platform = ", v_platform))
      # print(InputData_each[QCIndex])
      if(v_missing_rate > 0.5) {
        set(OutputData, Index_each, j, NA)
        list_QC_sample_missing_0.5 <- append(list_QC_sample_missing_0.5,
                                             list(list(metabolite = v_metab,
                                                       platform = v_platform,
                                                       batch = id_batch,
                                                       missing_rate = v_missing_rate
                                             )))
        next
      }

      v_median <- median(InputData_each[QCIndex], na.rm = TRUE)
      v_median

      set(OutputData, Index_each, j, InputData_each[Index_each]/v_median)
    }
  }
  close(con = pb)
  object@assayData <- OutputData
  object@miscData[['QC_sample_missing_0.5']] <- rbindlist(lapply(list_QC_sample_missing_0.5, as.data.frame.list))

  object@logs <- paste0(object@logs, "\n",
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": QCmatrix normalization.\n")
  return(object)
}




#' nearest QC sample normalization
#'
#' Normalization data by the median value of the nearest QC samples. For each metabolite, the values (eg. raw peak area data) were divided by the median value of nearest QC samples (eg. the nearest three QC samples). To identify the nearest QC samples, `@assayData` should be ordered by the injection order.
#'
#' @param object A Metabolite object. In the feature annotation slot `feature`, a platform column should be provided for metabolite measurement platform (eg. `PLATFORM`). The values in the `PLATFORM` column (eg. `Neg`, `Polar`, `Pos Early`, and `Pos Late`) are column names in the sample annotation `sample` to determine the batches of samples.
#' @param n_nearest_QCsample Number of nearest QC samples to calculate the median value. The default value is 3 (an outlier QC sample might be used if only n_nearest_QCsample = 1).
#' @param feature_platform The column name of feature platform for metabolite measurements (eg. `PLATFORM`).
#' @param QC_ID_pattern A character pattern to determine QC samples. Default value: "MTRX".
#' @param test test the function for the first 20 columns.
#' @param verbose print log information.
#' @seealso \code{\link{batch_norm}}, \code{\link{QCmatrix_norm}}
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' d <- nearestQC_norm(object = df)
#' }
#' @export
#'
#'
nearestQC_norm <- function(object, n_nearest_QCsample= 3, feature_platform = "PLATFORM", QC_ID_pattern = "MTRX", test = FALSE, verbose = TRUE) {
  stopifnot(inherits(object, "Metabolite"))
  OutputData <- copy(object@assayData)
  InputData <- copy(object@assayData)

  stopifnot(all(object@assayData[,get(object@sampleID)] == object@sampleData[, get(object@sampleID)]))

  # check platform column
  if(! feature_platform %in% names(object@featureData)) {
    stop(paste0("`", feature_platform, "` column does no exist in @featureData"), call. = FALSE)
  }
  platforms <- object@featureData[, get(feature_platform)]
  if(verbose) {
    cat("\nThe top ten sample ID (in injection order):\n")
    cat(paste0(object@assayData[,get(object@sampleID)][seq_len(10)], collapse = ", "),"\n")
    cat("\nPlatform information in @featureData:\n")
    print(table(platforms))
    cat("\n")
  }

  #
  if(all(toupper(unique(platforms)) %in% names(object@sampleData))) {
    if(verbose) {
      cat(paste0("\nSample size in platform ", toupper(sort(unique(platforms))[1])))
      print(table(object@sampleData[, toupper(unique(platforms)[1]), with = FALSE]))
      cat("\n")
    }
  } else {
    stop(paste0("`", paste0(unique(platforms), sep = "", collapse = " "), "` columns do no exist in @sampleData"), call. = FALSE)
  }

  # check sample IDs
  n_QC_sample <- sum(grepl(QC_ID_pattern, object@assayData[, get(object@sampleID)]))

  if(n_QC_sample == 0) {
    stop(paste0("No QC samples in @assayData (ID: ", object@sampleID,")"), call. = FALSE)
  }

  if(verbose) {
    cat(paste0("\n Number of QC samples n = ", n_QC_sample, "\n"))
  }

  v_n_col <- dim(object@assayData)[2]
  v_n_row <- dim(object@assayData)[1]

  if(isTRUE(test)) {
    v_n_col <- 10
  } else if(test >=2 ) {
    v_n_col <- test
  }

  pb <- txtProgressBar(min = 0, max = v_n_col, style = 3, file = stderr())


  sample_IDs <- object@assayData[,get(object@sampleID)]

  v_col_names <- names(object@assayData)


  for(j in 2L:v_n_col) {

    setTxtProgressBar(pb = pb, value = j)

    v_metab <- v_col_names[j]

    v_platform <- object@featureData[get(object@featureID) == v_metab, feature_platform, with = FALSE]
    v_platform <- toupper(v_platform)

    v_batch <- unlist(object@sampleData[, v_platform, with = FALSE])

    InputData_each <- unlist(InputData[, j, with = FALSE]) # first select the values

    for (id_batch in unique(v_batch[! is.na(v_batch)])) {

      Index_each <- which(v_batch == id_batch)

      QCIndex <- grep(QC_ID_pattern, sample_IDs)
      QCIndex <- QCIndex[QCIndex %in% Index_each]
      QCIndex_no_missing <- QCIndex[!is.na(InputData_each[QCIndex])]

      for(i in Index_each) {
        closeQC <- QCIndex_no_missing[which(abs(i-QCIndex_no_missing) %in%  sort(abs(i-QCIndex_no_missing))[seq_len(n_nearest_QCsample)])]
        v_median <- median(InputData_each[closeQC], na.rm = TRUE)
        set(OutputData, i, j, InputData_each[i]/v_median)
      }
    }
  }

  close(con = pb)
  object@assayData <- OutputData

  object@logs <- paste0(object@logs, "\n",
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Nearest QC normalization (n = ",n_nearest_QCsample, ").\n")
  return(object)
}


#' batch normalization
#'
#' Normalization data by the median value of each batch
#'
#' @param object A Metabolite object. In the feature annotation slot `feature`, a platform column should be provided for metabolite measurement platform (eg. `PLATFORM`). The values in the `PLATFORM` column (eg. `Neg`, `Polar`, `Pos Early`, and `Pos Late`) are column names in the sample annotation `sample` to determine the batches of samples.
#' @param feature_platform The column name of feature platform for metabolite measurements (eg. `PLATFORM`).
#' @param QC_ID_pattern A character pattern to determine QC samples. Default value: "MTRX".
#' @param test test the function for the first 20 columns.
#' @param verbose print log information.
#' @seealso \code{\link{QCmatrix_norm}}
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return A Metabolite object after normalization. 
#' @export
#'
#'
#'
batch_norm <- function(object, feature_platform = "PLATFORM", QC_ID_pattern = "MTRX", test = FALSE, verbose = TRUE) {
  stopifnot(inherits(object, "Metabolite"))
  OutputData <- copy(object@assayData)
  InputData <- copy(object@assayData)
  # check platform column
  if(! feature_platform %in% names(object@featureData)) {
    stop(paste0("`", feature_platform, "` column does no exist in @featureData"), call. = FALSE)
  }

  platforms <- object@featureData[, get(feature_platform)]

  if(verbose) {
    cat("\nPlatform information in @featureData:\n")
    print(table(platforms))
    cat("\n")
  }

  #
  if(all(toupper(unique(platforms)) %in% names(object@sampleData))) {
    if(verbose) {
      cat(paste0("\nSample size in platform ", toupper(sort(unique(platforms))[1])))
      print(table(object@sampleData[, toupper(unique(platforms)[1]), with = FALSE]))
      cat("\n")
    }
  } else {
    stop(paste0("`", paste0(unique(platforms), sep = "", collapse = " "), "` columns do no exist in @sampleData"), call. = FALSE)
  }

  # check sample IDs
  n_QC_sample <- sum(grepl("MTRX", object@assayData[, get(object@sampleID)]))

  if(n_QC_sample != 0) {
    warnings(paste0("Possible QC samples (MTRX) in @assayData (ID: ", object@sampleID, " n = ", n_QC_sample,")"), call. = FALSE)
  }

  v_n_col <- dim(object@assayData)[2]
  v_n_col

  if(isTRUE(test)) {
    v_n_col <- 10
  } else if(test >=2 ) {
    v_n_col <- test
  }

  pb <- txtProgressBar(min = 0, max = v_n_col, style = 3, file = stderr())

  sample_IDs <- object@assayData[,get(object@sampleID)]

  for(j in 2L:v_n_col) {
    setTxtProgressBar(pb = pb, value = j)

    v_metab <- names(object@assayData)[j]

    v_platform <- object@featureData[get(object@featureID) == v_metab, get(feature_platform)]
    v_platform <- toupper(v_platform)

    v_batch <- unlist(object@sampleData[, v_platform, with = FALSE])
    InputData_each <- unlist(InputData[, j, with = FALSE])

    for (id_batch in unique(v_batch[! is.na(v_batch)])) {

      Index_each <- (v_batch  == id_batch)
      QCIndex <- grepl(QC_ID_pattern, sample_IDs)

      sampleIndex <- Index_each &  (!QCIndex)
      v_median <- median(InputData_each[sampleIndex], na.rm = TRUE)

      set(OutputData, which(Index_each), j, InputData_each[which(Index_each)]/v_median)
    }
  }

  close(con = pb)
  object@assayData <- OutputData
  object@logs <- paste0(object@logs, "\n",
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Batch normalization.\n")
  return(object)
}





#' LOESS normalization
#'
#' Normalization data by machine learning modelling, eg. locally estimated scatterplot smoothing (LOESS) on QC samples in each batch. For each metabolite, the values (eg. raw peak area data) were divided by the median value of QC samples in that batch. QC samples and metabolite batches should be specified (see parameters below).
#'
#' @param object A Metabolite object. In the feature annotation slot `feature`, a platform column should be provided for metabolite measurement platform (eg. `PLATFORM`). The values in the `PLATFORM` column (eg. `Neg`, `Polar`, `Pos Early`, and `Pos Late`) are column names in the sample annotation `sample` to determine the batches of samples.
#' @param method Modelling method for the normalization, currently support LOESS and KNN. 
#' @param feature_platform The column name of feature platform for metabolite measurements (eg. `PLATFORM`).
#' @param span default value 0.4
#' @param degree default value 2
#' @param k Number of neighbors in KNN modelling (default value 3)
#' @param QC_ID_pattern A character pattern to determine QC samples. Default value: "MTRX".
#' @param test test the function for the first 20 columns.
#' @param verbose print log information.
#' @seealso \code{\link{batch_norm}}
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' d <- QCmatrix_norm(object = df)
#' }
#' @export
#'
#'
modelling_norm <- function(object, method = c("LOESS", "KNN"), feature_platform = "PLATFORM" , QC_ID_pattern = "MTRX", span = 0.75, degree = 2, k = 3, test = FALSE, verbose = TRUE) {

  stopifnot(inherits(object, "Metabolite"))
  
  method <- match.arg(method)
  
  OutputData <- copy(object@assayData)
  InputData <- copy(object@assayData)
  # check platform column
  if(! feature_platform %in% names(object@featureData)) {
    stop(paste0("`", feature_platform, "` column does no exist in @featureData"), call. = FALSE)
  }

  platforms <- object@featureData[, get(feature_platform)]
  if(verbose) {
    cat("\nPlatform information in @featureData:\n")
    print(table(platforms))
    cat("\n")
  }

  if(all(toupper(unique(platforms)) %in% names(object@sampleData))) {
    if(verbose) {
      cat(paste0("\nSample size in platform ", toupper(sort(unique(platforms))[1])))
      print(table(object@sampleData[, toupper(unique(platforms)[1]), with = FALSE]))
      cat("\n")
    }
  } else {
    stop(paste0("`", paste0(unique(platforms), sep = "", collapse = " "), "` columns do no exist in @sampleData"), call. = FALSE)
  }

  # check sample IDs
  n_QC_sample <- sum(grepl(QC_ID_pattern, object@assayData[, get(object@sampleID)]))

  if(n_QC_sample == 0) {
    stop(paste0("No QC samples in @assayData (ID: ", object@sampleID,")"), call. = FALSE)
  }

  if(verbose) {
    cat(paste0("\n Number of QC samples n = ", n_QC_sample, "\n"))
  }


  v_n_col <- dim(object@assayData)[2]
  v_n_col

  if(isTRUE(test)) {
    v_n_col <- 10
  } else if(test >=2 ) {
    v_n_col <- test
  }

  pb <- txtProgressBar(min = 0, max = v_n_col, style = 3, file = stderr())

  sample_IDs <- object@assayData[,get(object@sampleID)]

  list_norm_fail <- c()

  for(j in 2L:v_n_col) {

    setTxtProgressBar(pb = pb, value = j)
    v_metab <- names(object@assayData)[j]

    v_platform <- object@featureData[get(object@featureID) == v_metab, get(feature_platform)]
    v_platform <- toupper(v_platform)

    v_batch <- unlist(object@sampleData[, v_platform, with = FALSE])

    InputData_each <- unlist(InputData[, j, with = FALSE])

    for (id_batch in unique(v_batch[! is.na(v_batch)])) {

      Index_each <- which(v_batch == id_batch)

      QCIndex <- grep(QC_ID_pattern, sample_IDs)
      QCIndex <- QCIndex[QCIndex %in% Index_each]
      
      
      if(method == "LOESS") {
        fit <- tryCatch(
          stats::loess(InputData_each[QCIndex] ~ QCIndex, span = span, degree = degree),
          error = function(e) {
            cat(paste0("Failed to fit model: ", e), "\n")
          })
        
        if(inherits(fit, "loess") && (!is.na(fit$enp))) {
          yfit <- predict(fit, Index_each)
          yfit <- ifelse(yfit <=0, NA, yfit)
          set(OutputData, Index_each, j, InputData_each[Index_each]/yfit)
        } else {
          set(OutputData, Index_each, j, NA)
          list_norm_fail <- c(list_norm_fail, v_metab)
        }
        
      } else if(method == "KNN") {
        # missing values?
        check_pkg("FNN")
        fit <- tryCatch(
          FNN::knn.reg(train = matrix(QCIndex), test = matrix(Index_each), y = InputData_each[QCIndex], k = k),
          error = function(e) {
            cat(paste0("Failed to fit model: ", e), "\n")
          })
        
        if(inherits(fit, "knnReg")) {
          yfit <- fit$pred
          yfit <- ifelse(yfit <=0, NA, yfit)
          set(OutputData, Index_each, j, InputData_each[Index_each]/yfit)
        } else {
          set(OutputData, Index_each, j, NA)
          list_norm_fail <- c(list_norm_fail, v_metab)
        }
      }
    }
  }
  close(con = pb)
  object@assayData <- OutputData
  
  if(method == "LOESS") {
    object@miscData[['list_norm_fail_LOESS']] <- list_norm_fail
    object@logs <- paste0(object@logs, "\n", format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": LOESS normalization.\n")
  } else {
    object@miscData[['list_norm_fail_KNN']] <- list_norm_fail
    object@logs <- paste0(object@logs, "\n", format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": KNN normalization.\n")
  }
  return(object)
}



############### imputation #######################



#' @rdname impute
#' @export
#' @examples
#'
#'
#' \dontrun{
#' d <- impute(object)
#' }
#'
impute.Metabolite <- function(object, method = c('half-min', "median", "mean", "zero", "kNN")) {

  method <- match.arg(method)
  
  if(method == "kNN") {
    return(impute_kNN(object))
  }

  object@assayData <- cbind(object@assayData[, 1],
                            apply(object@assayData[, -1], 2, impute, method = method))

  object@logs <- paste0(object@logs,
                        format(Sys.time(), "%d/%m/%y %H:%M:%OS"),
                        ": Impute data using `",method, "` method. \n")
  return(object)
}

#' @rdname impute
#' @export
#' @note default method is used for a vector
impute.default <- function(object, method = "half-min") {
  if(is.vector(object)) {
    x <- as.numeric(object)

    stopifnot(method %in% c('half-min', "median", "mean", "zero"))
    x_replace <- switch (method,
                         'half-min' = min(x, na.rm = TRUE)/2,
                         "median" = median(x, na.rm = TRUE),
                         "mean" = mean(x, na.rm = TRUE),
                         "zero" = 0
    )
    x <- ifelse(is.na(x), x_replace, x)
    return(x)
  } else stop("Unknown object, convert to vector before impute.", call. = FALSE)

}


#' @rdname impute
#' @export
#' @note `impute_kNN`: Imputation using nearest neighbor averaging (kNN) method, the input is a Metabolite object, assayData was first transposed to row as metabolties and column as samples. 
impute_kNN <- function(object) {
  check_pkg("impute")
  
  x <- transpose(object@assayData, make.names = object@sampleID)
  x <- impute::impute.knn(as.matrix(x))$data
  x <- transpose(as.data.table(x), keep.names = object@sampleID)
  
  setnames(x, names(object@assayData))
  object@assayData <- x
  return(object)
}

############### correlation #######################


#' correlation of features between two Metabolite objects
#'
#' Calculate the correlation of features between two Metabolite objects
#'
#' @param object_X The first Metabolite object.
#' @param object_Y The second Metabolite object.
#' @param method a character string to calculate correlation coefficient. One of "pearson" (default), "kendall", or "spearman".
#' @param verbose print log information.
#' @seealso \code{\link{cor}}
#' @export
#' @return A data.table with correlation coefficients.
#' 
correlation <- function(
  object_X = NULL,
  object_Y = NULL,
  method = "pearson",
  verbose = TRUE
)  {

  v_feature <- setdiff(intersect(names(object_X@assayData), names(object_Y@assayData)),
                       c(object_X@sampleID, object_Y@sampleID))

  v_sample <- intersect(object_X@assayData[, get(object_X@sampleID)], object_Y@assayData[, get(object_Y@sampleID)])

  if(verbose) {
    cat("Identify ", length(object_X@assayData[, get(object_X@sampleID)]), " samples in data A, and ", length(object_Y@assayData[, get(object_Y@sampleID)]), " samples in data B, with ", length(v_sample), " overlap samples.\n")
    cat("Identify ", length(names(object_X@assayData))-1, " features in data A, and ", length(names(object_Y@assayData))-1, " features in data B, with ", length(v_feature), " overlap features.\n")
  }

  df_merge <- merge(object_X@assayData, object_Y@assayData, by.x = object_X@sampleID, by.y = object_Y@sampleID)

  f_get_cor <- function( x = "", y = "", data = NULL, method = NULL) {
    res <- cor.test(unlist(data[, x, with = FALSE]), unlist(data[, y, with = FALSE]), method = method)
    return(c(b = res$estimate, p = res$p.value))
  }

  res <- NULL
  p <- progressr::progressor(along = length(v_feature))
  for(i in seq_along(v_feature)) {
    p(sprintf("i=%g", i))
    v_i <- v_feature[i]

    df_one <- tryCatch(
        f_get_cor(paste0(v_i, ".x"), paste0(v_i, ".y"), data = df_merge, method = method),
        error = function(e) {
          # cat(paste0("Failed : ", v_i, " ", e, "\n"))
          NA
        })
    if(length(df_one) == 2 ) res <- rbind(res, c(term = v_i, df_one))
  }

  res <- data.table(res)
  setnames(res, 2, "r")
  res[, c("r", "p") := lapply(.SD, as.numeric), .SDcols =  c("r", "p")]
  res <- res[complete.cases(res), ]
  res <- res[order(res$r), ]
  return(res)
}

