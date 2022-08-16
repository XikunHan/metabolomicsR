

####### Metabolite class #####

#' The Metabolite class
#'
#' The Metabolite object is a representation of metabolomic data, metabolomic annotation, and sample annotation.
#' @slot assayData a data.frame or data.table of metabolite measurements (peak area data or normalized data, sample [row] * feature [column]).
#' @slot featureData a data.frame or data.table of metabolite annotation (chemical annotation)
#' @slot sampleData a data.frame or data.table of sample annotation (sample meta data).
#' @slot featureID a character of the metabolite ID column (in feature file and the column names of data), default: CHEM_ID (provided from Metabolon file).
#' @slot sampleID a character of the sample ID column (in sample and the first column of data), default: PARENT_SAMPLE_NAME (provided from Metabolon file).
#' @slot logs Log information of data analysis process.
#' @slot miscData Ancillary data.
#' @name Metabolite-class
#' @export
#' @return A Metabolite class. 
#' @seealso \code{\link{Metabolite}}, \code{\link{load_excel}}, \code{\link{load_data}}
#'
Metabolite <- setClass(
  Class = 'Metabolite',
  slots = list(
    'assayData' = 'data.table',
    'featureData' = 'data.table',
    'sampleData' = 'data.table',
    'featureID' = 'character',
    'sampleID' = 'character',
    'logs' = 'character',
    'miscData' = 'list'
  )
)




#' get assayData
#'
#'  Accessors for Metabolite object. Get the assayData in the Metabolite object.
#' @param object A Metabolite object.
#' 
setGeneric("assayData", function(object) standardGeneric("assayData"))

#' @rdname assayData
#' @docType methods
#' @export
#' @return A data.table of assayData. 
setMethod("assayData", "Metabolite", function(object) object@assayData)




#' set assayData
#'
#'  Accessors for Metabolite object. `assayData<-` will update the assayData in the Metabolite object.
#' @param object A Metabolite object.
#' @param value The new assayData.
#' @rdname assayData_set
#' @export
#' @return A data.table of assayData. 
setGeneric("assayData<-", function(object, value) standardGeneric("assayData<-"))

#' @docType methods
#' @rdname assayData_set
#' @export
setMethod("assayData<-", "Metabolite", function(object, value) {
  stopifnot(inherits(object@assayData, "data.frame"))
  object@assayData <- as.data.table(value)
  update_Metabolite(object)
})


### featureData

#' get featureData
#'
#'  Accessors for Metabolite object. Get the featureData in the Metabolite object.
#' @param object A Metabolite object.
setGeneric("featureData", function(object) standardGeneric("featureData"))

#' @rdname featureData
#' @docType methods
#' @export
#' @return A data.table of featureData. 
setMethod("featureData", "Metabolite", function(object) object@featureData)


#' set featureData
#'
#'  Accessors for Metabolite object. `featureData<-` will update the featureData in the Metabolite object.
#' @param object A Metabolite object.
#' @param value The new featureData.
#' @rdname featureData_set
#' @export
#' @return A data.table of featureData. 
setGeneric("featureData<-", function(object, value) standardGeneric("featureData<-"))

#' @docType methods
#' @rdname featureData_set
#' @export
setMethod("featureData<-", "Metabolite", function(object, value) {
  stopifnot(inherits(object@featureData, "data.frame"))
  object@featureData <- as.data.table(value)
  update_Metabolite(object)
})


### sampleData

#' get sampleData
#'
#'  Accessors for Metabolite object. Get the sampleData in the Metabolite object.
#' @param object A Metabolite object.
setGeneric("sampleData", function(object) standardGeneric("sampleData"))

#' @rdname sampleData
#' @docType methods
#' @export
#' @return A data.table of sampleData. 
setMethod("sampleData", "Metabolite", function(object) object@sampleData)




#' set sampleData
#'
#'  Accessors for Metabolite object. `sampleData<-` will update the sampleData in the Metabolite object.
#' @param object A Metabolite object.
#' @param value The new sampleData.
#' @rdname sampleData_set
#' @export
#' @return A data.table of sampleData. 
setGeneric("sampleData<-", function(object, value) standardGeneric("sampleData<-"))

#' @docType methods
#' @rdname sampleData_set
#' @export
setMethod("sampleData<-", "Metabolite", function(object, value) {
  stopifnot(inherits(object@sampleData, "data.frame"))
  object@sampleData <- as.data.table(value)
  update_Metabolite(object)
})



####### show #####

#' Print a Metabolite class object
#'
#' Print a Metabolite class object
#' @param object A Metabolite object.
#' @docType methods
#' @export
#' @return print a Metabolite object. 
setMethod("show", "Metabolite", function(object) {
  cat("An object of ", is(object), "\n")
  
  validObject(object)
  
  cat("\n***  @assayData (first and last 10 columns [",NROW(object@assayData), " * ", NCOL(object@assayData), "])  ***\n")

  # first 10 and last 10
  print(object@assayData[, head_tail(names(object@assayData)), with = FALSE])

  cat("\n***  @featureData (ID: ",object@featureID, ") ***\n")
  print(object@featureData)

  cat("\n***  @sampleData (ID: ",object@sampleID, ") ***\n")
  print(object@sampleData)

  cat("\n***  @miscData  ***\n")
  print(vapply(object@miscData, length, integer(1)))

  cat("\n***  @logs  ***\n")
  cat(object@logs)
})




#' Create a Metabolite object
#'
#' Create a Metabolite object from three input data sets:
#' 1) metabolite measurements (eg. peak area data or normalized data), and
#' 2) metabolite annotation (eg. chemical annotation)
#' 3) sample annotation (eg. sample meta data)
#'
#' @param assayData a data.frame or data.table of metabolite measurements (peak area data or normalized data, sample [row] * feature [column]).
#' @param featureData a data.frame or data.table of metabolite annotation (chemical annotation)
#' @param sampleData a data.frame or data.table of sample annotation (sample meta data).
#' @param featureID a character of the metabolite ID column (in feature file and the column names of data), default: CHEM_ID (provided from Metabolon file).
#' @param sampleID a character of the sample ID column (in sample and the first column of data), default: PARENT_SAMPLE_NAME (provided from Metabolon file).
#' @param logs Log information.
#' @param miscData Ancillary data.
#' @return A Metabolite object with slots: assayData, featureData, and sampleData.
#' @seealso \code{\link{Metabolite}}, \code{\link{load_excel}}, \code{\link{load_data}}
#' @export
#' @return A Metabolite object. 
#' @examples
#' # df <- create_Metabolite(assayData = df_data, featureData = df_feature, sampleData =  df_sample)
create_Metabolite <- function(
  assayData,
  featureData,
  sampleData,
  featureID = "featureID",
  sampleID = "sampleID",
  logs = "",
  miscData = list()
) {

  setDT(assayData)
  setDT(featureData)
  setDT(sampleData)

  
  featureID_ <- featureID
  featureID <- NULL
  
  assayData[, (setdiff(names(assayData), sampleID)) := lapply(.SD, as.numeric) , .SDcols= setdiff(names(assayData), sampleID)]

  if("featureID" %in% names(featureData)) {
    warnings(paste0("`featureID` column already exists in the feature annotation file. Will be overwritten!"))
  }

  feature_IDs <- setdiff(names(assayData), sampleID)
  sample_IDs <- unlist(assayData[, sampleID, with = FALSE])

  # test feature ID, if a numeric, add X (provided from Metabolon file)
  if(!is.na(as.integer(feature_IDs[1]))) {
    cat(paste0("\n Add X to feature IDs.\n"), file = stderr())
    # featureData[, featureID := paste0("X", get(featureID_))]
    if(! "featureID" %in% names(featureData)) {
      featureData <- cbind(featureID = paste0("X", unlist(featureData[, featureID_, with = FALSE])), featureData)
    } else {
      featureData <- cbind(featureID = paste0("X", unlist(featureData[, featureID_, with = FALSE])), featureData[, -c("featureID"), with = FALSE])
    }
    setnames(assayData, feature_IDs, paste0("X", feature_IDs))
    
  } else {
    # featureData[, featureID := get(featureID_)]
    if(! "featureID" %in% names(featureData)) {
      featureData <- cbind(featureID = unlist(featureData[, featureID_, with = FALSE]), featureData)
    } else {
      featureData <- cbind(featureID = unlist(featureData[, featureID_, with = FALSE]), featureData[, -c("featureID"), with = FALSE])
    }
  }

  setnames(sampleData, sampleID ,"sampleID")
  setnames(assayData, sampleID ,"sampleID")
  
  # sampleID in the first column
  
  stopifnot(names(sampleData)[1] == "sampleID")
  stopifnot(names(assayData)[1] == "sampleID")
  
  # check extra features and samples

  v_feature <- intersect(featureData$featureID, setdiff(names(assayData), "sampleID"))

  v_featureData <- setdiff(featureData$featureID, v_feature)
  if(length(v_featureData) > 0) {
    warning(paste0(length(v_featureData), " extra features in @featureData: ", paste5(v_featureData)))
  }

  featureData <- featureData[featureID %in% v_feature]

  v_assayData <- setdiff(setdiff(names(assayData), "sampleID"), v_feature)

  if(length(v_assayData) > 0) {
    warning(paste0(length(v_assayData), " extra features in @assayData: ", paste5(v_assayData)))
  }
  assayData <- assayData[, c("sampleID", v_feature), with = FALSE]

  v_sample <- intersect(assayData$sampleID, sampleData$sampleID)

  v_sampleData <- setdiff(sampleData$sampleID, v_sample)
  if(length(v_sampleData) > 0) {
    warning(paste0(length(v_sampleData), " extra samples in @sampleData: ", paste5(v_sampleData)))
  }

  sampleData <- sampleData[get("sampleID") %in% v_sample]
  v_assayData <- setdiff(assayData$sampleID,  v_sample)
  if(length(v_assayData) > 0) {
    warning(paste0(length(v_assayData), " extra samples in @assayData: ", paste5(v_assayData)))
  }
  assayData <- assayData[get("sampleID") %in% v_sample]

  new(
    Class = 'Metabolite',
    assayData = assayData,
    featureData = featureData,
    sampleData = sampleData,
    featureID = "featureID",
    sampleID = "sampleID",
    logs = paste0(logs, format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": Initiate data: ", NROW(sampleData), " samples and ", NROW(featureData), " features.\n"),
    miscData = miscData
  )
}



####### setValidity #####

setValidity("Metabolite", function(object) {

  feature_IDs <- setdiff(names(object@assayData), object@sampleID)
  sample_IDs <- object@assayData[, get(object@sampleID)]

  msg <- character()

  if (! object@featureID %in% names(object@featureData)) {
    msg <- append(msg, paste0(object@featureID, " does not exist in @featureData"))
  }
  if (! object@sampleID %in% names(object@assayData)) {
    msg <- append(msg, paste0(object@sampleID, " does not exist in @assayData"))
  }


  if(! all(feature_IDs %in%  object@featureData[, get(object@featureID)])) {
    msg_list <- paste0(feature_IDs[which(! feature_IDs %in%  object@featureData[, get(object@featureID)])[seq_len(5)]], collapse = ", ")
    msg <- append(msg, paste0(msg, "\n", "Some feature IDs (",msg_list, ") are missing in @featureData"))
  }

  if(! all(sample_IDs %in%  object@sampleData[, get(object@sampleID)])) {
    msg_list <- paste0(sample_IDs[which(! sample_IDs %in%  object@sampleData[, get(object@sampleID)])[seq_len(5)]], collapse = ", ")

    msg <- append(msg, paste0(msg, "\n", "Some sample IDs (", msg_list, ") are missing in @sampleData"))
  }
  if (length(msg) > 0) {
    return(paste(msg, sep = "\n"))
  }
  TRUE
})




#' column missing rate
#'
#' Calculate column missing rate -- metabolite missingness.
#'
#' @param object An object, data.frame, data.table or Metabolite.
#' @return Returns a vector of the missing rate for each column
#' @rdname column_missing_rate
#' @export
#' @return A data.table of column missing rate.
#'
column_missing_rate <- function(object) {
  UseMethod(generic = 'column_missing_rate', object = object)
}



#' filter columns using missing rate
#'
#' Remove columns below a specific missing rate threshold.
#'
#' @param object An object, data.frame, data.table or Metabolite.
#' @param threshold missing rate threshold, default is 0.5. Other values: 0.2, 0.8.
#' @param verbose print log information.
#' @export
#' @return An object after filtering column missing rate. 
#' @rdname filter_column_missing_rate
#'
filter_column_missing_rate <- function(object, threshold, verbose) {
  UseMethod(generic = 'filter_column_missing_rate', object = object)
}



#' row missing rate
#'
#' Calculate row missing rate -- sample missingness.
#'
#' @param object An object, data.frame, data.table or Metabolite.
#' @return Returns a vector of the missing rate for each row
#' @rdname row_missing_rate
#' @export
#' @return A data.table of row missing rate.
row_missing_rate <- function(object) {
  UseMethod(generic = 'row_missing_rate', object = object)
}



#' filter rows using missing rate
#'
#' Remove samples below a specific missing rate threshold.
#'
#' @param object An object, data.frame, data.table or Metabolite.
#' @param threshold missing rate threshold, default is 0.5. Other values: 0.2, 0.8.
#' @param verbose print log information.
#' @export
#' @rdname filter_row_missing_rate
#'
filter_row_missing_rate <- function(object, threshold, verbose) {
  UseMethod(generic = 'filter_row_missing_rate', object = object)
}



#' filter columns if values are constant
#'
#' Remove columns if values are constant
#'
#' @param object An object, data.frame, data.table or Metabolite.
#' @param verbose print log information.
#' @export
#' @rdname filter_column_constant
#'
filter_column_constant <- function(object, verbose) {
  UseMethod(generic = 'filter_column_constant', object = object)
}




#' subset a Metabolite object.
#'
#' subset a Metabolite object.
#'
#' @param object An object, data.frame, data.table or Metabolite.
#' @param subset logical expression indicating rows to keep (samples). Expression will be evaluate in the `@sampleData`.
#' @param select expression indicating columns to select (features). See \code{\link[base]{subset}}. Expression will be evaluate in the `@assayData`.
#' @export
#' @return An object after subsetting rows or columns. 
#' @rdname subset
#'
subset <- function(object, subset, select) {
  UseMethod(generic = 'subset', object = object)
}



#' impute missing values
#'
#' impute missing values
#'
#' @param object An object, a vector, data.frame, data.table or Metabolite.
#' @param method Imputation method, the default method is half the minimum value (`half-min`) of the metabolite. Currently support 'half-min', "median", "mean", "zero", "kNN".
#' @export
#' @return An object after imputing missing values. 
#' @rdname impute
#' @note Wei, R., Wang, J., Su, M. et al. Missing Value Imputation Approach for Mass Spectrometry-based Metabolomics Data. Sci Rep 8, 663 (2018). https://doi.org/10.1038/s41598-017-19120-0
#' 
#'
impute <- function(object,  method) {
  UseMethod(generic = 'impute', object = object)
}



#' change outlier values as NA or winsorize
#'
#' @param object An object, a vector, data.frame, data.table or Metabolite.
#' @param method Replace outlier value method, the default method is `winsorize`: replace the outlier values by the maximum and/or minimum values of the remaining values. `as_NA`: set as NA (do not use this method if using half-min imputation).
#' @param nSD Define the N times of the SD as outliers.
#' @export
#' @return An object after replacing outlier values. 
#' @rdname replace_outlier
#'
replace_outlier <- function(object, method, nSD) {
  UseMethod(generic = 'replace_outlier', object = object)
}



#' outlier rate
#'
#' Calculate outlier rate.
#'
#' @param object An object, vector, data.frame, data.table or Metabolite.
#' @param nSD N times of the SD as outliers.
#' @return Returns a vector of the outlier rate.
#' @rdname outlier_rate
#' @export
#' @return A data.table of outlier rate. 
#'
outlier_rate <- function(object, nSD) {
  UseMethod(generic = 'outlier_rate', object = object)
}
