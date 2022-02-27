
############### PCA #######################

#' Principal Components Analysis
#'
#' Performs a principal components analysis on the Metabolite object.
#'
#' @param object A Metabolite object.
#' @param nPCs Number of principal components to be calculated. Default value 10.
#' @param impute_method Imputation method, the default method is half the minimum value (`half-min`) of the metabolite. Currently support 'half-min', "median", "mean", "zero". `NULL` without imputation.
#' @param log Performs natural logarithm transformation before PCA analysis.
#' @param scale scale feature in the PCA calculation.
#' @param addPC If TRUE, merge PCs with `@sampleData` and return the `object`, else return `PC`.
#' @importFrom stats prcomp
#' @export
#' @examples
#'
#' # skip the first column (eg. ID) to impute missing values
#' \dontrun{
#' d <- run_PCA(object)
#' }
#'

run_PCA <- function(object, nPCs = 10, impute_method = "half-min", log = TRUE, scale = TRUE, addPC = TRUE) {

  if(! is.null(impute_method)) {
    object <- impute(object, method = impute_method)
  }

  X <- as.matrix(object@assayData[, -1])

  if(log) {
    X[] <- apply(X[], 2, log)
  }

  #   X_center <- colMeans(X, na.rm = TRUE)
  #   X_scale <- apply(X, 2, sd, na.rm = TRUE)
  #
  #   # reference: https://slowkow.com/notes/pca-benchmark/
  # suppressWarnings({
  # res <- irlba::irlba(A = X, nv = nPCs, center = X_center, scale = X_scale)
  # })
  #
  # PCs <- res$u %*% diag(res$d)
  # colnames(PCs) <- paste0("PC", 1: nPCs)
  # PCs <- cbind(object@assayData[, 1], PCs)
  #
  # Variances <- res$d^2/sum(res$d^2) # to check
  # Variances
  # res$d/sum(res$d)

  res <- prcomp(X, center = TRUE, scale = scale, rank. = nPCs)
  PCs <- cbind(object@assayData[, 1], res$x)
  Variances <- (res$sdev^2/sum(res$sdev^2))[1:nPCs]

if (addPC) {
  object@sampleData <- merge(object@sampleData, PCs, by = object@sampleID, all = TRUE, sort = FALSE)
  object@miscData[["Variances_explained_PCs"]] <- Variances
  return(object)

}  else {
  return(list(PCs = PCs, Variances_explained = Variances))

}

}
