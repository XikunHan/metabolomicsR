#' @import data.table
#' @import ggplot2
#' @importFrom stats median qnorm sd na.omit cor.test complete.cases as.formula binomial predict glm
#' @importFrom methods new is show validObject
#' @importFrom utils head str
#' @importFrom plotROC geom_roc melt_roc stat_roc
#'
NULL



#' first 10 and last 10 values
#' @param x A vector
#' @noRd
head_tail <- function(x) {
  stopifnot(is.vector(x))
  v_1 <- x[seq_len(10)]
  v_2 <- x[ifelse(length(x) -10 > 0, length(x) -10, 1):length(x)]
  return(unique(na.omit(c(v_1, v_2))))
}


#' paste the first five values
#' @param x A vector
#' @noRd
paste5 <- function(x) {
  if(length(x) >=5) return(paste0(paste0(x[seq_len(5)], collapse = " "), " ...")) else return(paste0(x, collapse = " "))
}


#' calculate missing rate
#' @param x A vector
#' @noRd
missing_rate <- function(x) sum(is.na(x)) / length(x)

#' point shape used in scatter plot
#' @noRd
point_shape_plate <- c(5, 8 , 16, 17,  1, 2, 3, 7, 10,  18,  25)


# utils::globalVariables(c("n"))



#' check R package
#' @param x A character of pacakge name
#' @noRd
check_pkg <- function(x) {
  if (! requireNamespace(x, quietly = TRUE)) {
    stop(paste0("Package \"" , x,"\" must be installed to use this function."), call. = FALSE)
  }
}



#' reorder injection file, internal use
#' @param data A data.table
#' @noRd
f_injection_reorder <- function(data = NULL) {

  # internal function
  data$SUBSET <- data$SUBSET_NEG

  dd <- as.data.table(data[, table(data$SUBSET, data$SUBSET_POLAR)])
  dd <- dd[dd$N>0]
  dd$N <- NULL
  dd[] <- lapply(dd, as.integer)
  setnames(dd, c("SUBSET", "V2"))

  data <- merge(data, dd, by.x = "SUBSET_POLAR", by.y = "V2", all.x = TRUE, suffixes = c("", "_x"))
  data$SUBSET <- ifelse(is.na(data$SUBSET), data$SUBSET_x, data$SUBSET)
  data$SUBSET_x <- NULL


  dd <- as.data.table(data[, table(data$SUBSET, data$`SUBSET_POS EARLY`)])
  dd <- dd[dd$N>0]
  dd$N <- NULL
  dd[] <- lapply(dd, as.integer)
  setnames(dd, c("SUBSET", "V2"))

  data <- merge(data, dd, by.x = "SUBSET_POS EARLY", by.y = "V2", all.x = TRUE, suffixes = c("", "_x"))
  data$SUBSET <- ifelse(is.na(data$SUBSET), data$SUBSET_x, data$SUBSET)
  data$SUBSET_x <- NULL


  dd <- as.data.table(data[, table(data$SUBSET, data$`SUBSET_POS LATE`)])
  dd <- dd[dd$N>0]
  dd$N <- NULL
  dd[] <- lapply(dd, as.integer)
  setnames(dd, c("SUBSET", "V2"))

  data <- merge(data, dd, by.x = "SUBSET_POS LATE", by.y = "V2", all.x = TRUE, suffixes = c("", "_x"))
  data$SUBSET <- ifelse(is.na(data$SUBSET), data$SUBSET_x, data$SUBSET)
  data$SUBSET_x <- NULL

  data$RUN_ORDER <- ifelse(!is.na(data$SAMPLE_RUN_ORDER_NEG), data$SAMPLE_RUN_ORDER_NEG,
                             ifelse(!is.na(data$SAMPLE_RUN_ORDER_POLAR), data$SAMPLE_RUN_ORDER_POLAR,
                                    ifelse(!is.na(data$`SAMPLE_RUN_ORDER_POS EARLY`), data$`SAMPLE_RUN_ORDER_POS EARLY`,data$`SAMPLE_RUN_ORDER_POS LATE` )))

  data <- data[order(data$SUBSET, data$RUN_ORDER)]
  return(data)
}


