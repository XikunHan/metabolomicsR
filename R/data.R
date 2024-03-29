#' Example data.
#'
#' A dataset containing 356 samples and 758 features.
#' @usage data(df_plasma)
#' 
"df_plasma"





#' Load metabolite data from an excel file
#'
#' @param path Path to the xls/xlsx file.
#' @param data_sheet A integer of xlsx sheet number for metabolite measurements 
#' (peak area data or normalized data, sample [row] * feature [column])
#' @param feature_sheet A integer of xlsx sheet number for metabolite annotation (chemical annotation)
#' @param sample_sheet A integer of xlsx sheet number for sample annotation (sample meta data)
#' @param featureID a character of the metabolite ID column (in feature file and the column names of data file), 
#' default: CHEM_ID (provided from Metabolon file)
#' @param sampleID a character of the sample ID column (in sample file and the first column of data file), 
#' default: PARENT_SAMPLE_NAME (provided from Metabolon file)
#' @return A Metabolite object with slots: assayData, featureData, and sampleData.
#' @export
#' @examples
#'
#' file_path <- system.file("extdata", "QMDiab_metabolomics_OrigScale.xlsx", 
#' package = "metabolomicsR", mustWork = TRUE)
#' 
#' df_plasma <- load_excel(path = file_path, data_sheet = 1, feature_sheet = 4, sample_sheet = 8, 
#' sampleID = "QMDiab-ID", featureID = "BIOCHEMICAL")
#'
load_excel <- function(
  path,
  data_sheet = NULL,
  feature_sheet = NULL,
  sample_sheet = NULL,
  featureID = "CHEM_ID",
  sampleID ="PARENT_SAMPLE_NAME"
) {

  check_pkg("readxl")

  df_data <- setDT(readxl::read_excel(path, sheet = data_sheet))
  # as.numeric
  df_data[, (setdiff(names(df_data), sampleID)) := lapply(.SD, as.numeric) , .SDcols= setdiff(names(df_data), sampleID)]

  df_feature <- setDT(readxl::read_excel(path, sheet = feature_sheet))
  df_sample <- setDT(readxl::read_excel(path, sheet = sample_sheet))

  logs <- paste0(format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": Import data from: ", path, " .\n")

  object <- create_Metabolite(assayData = df_data,
             featureData = df_feature,
             sampleData = df_sample,
             featureID = featureID,
             sampleID = sampleID,
             logs = logs
  )
  return(object)
}



#' Load metabolite data from three separate files
#'
#' Load metabolite data from three separate files (import files using `fread` from data.table).
#'
#' @param data_path Path to the metabolite measurements (peak area data or normalized data, 
#' sample [row] * feature [column])
#' @param feature_path Path to the metabolite annotation (chemical annotation)
#' @param sample_path Path to the sample annotation (sample meta data)
#' @param featureID a character of the metabolite ID column (in feature file and the column names of data file), 
#' default: CHEM_ID (provided from Metabolon file)
#' @param sampleID a character of the sample ID column (in sample file and the first column of data file), 
#' default: PARENT_SAMPLE_NAME (provided from Metabolon file).
#' @return A Metabolite object with slots: assayData, featureData, and sampleData.
#' @export
load_data <- function(
  data_path = NULL,
  feature_path = NULL,
  sample_path = NULL,
  featureID = "CHEM_ID",
  sampleID = "PARENT_SAMPLE_NAME"
)  {

  df_data <- fread(data_path)
  # as.numeric
  df_data[, (setdiff(names(df_data), sampleID)) := lapply(.SD, as.numeric) , .SDcols= setdiff(names(df_data), sampleID)]

  df_feature <- fread(feature_path)
  df_sample <- fread(sample_path)


  logs <- paste0(format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": Import data from: ", 
                 paste5(c(data_path, feature_path, sample_path)), " .\n")

  object <- create_Metabolite(assayData = df_data,
             featureData = df_feature,
             sampleData = df_sample,
             featureID = featureID,
             sampleID = sampleID,
             logs = logs
              )

  return(object)
}


#' Save metabolite data in txt files 
#'
#' Save metabolite data in separate txt files
#'
#' @param object A Metabolite object
#' @param file Output file to save the metabolite measurements 
#' (suffixes: "_assay.txt", "_feature_annotation.txt", "_sample_annotation.txt", "_logs.txt). 
#' @export
#' @return No return value.
#'
#'
save_data_txt <- function(object, file = "")  {
  fwrite(object@assayData, file = paste0(file, "_assay.txt"), sep = "\t", na = "NA", quote = FALSE)
  fwrite(object@featureData, file = paste0(file, "_feature_annotation.txt"), sep = "\t", na = "NA", quote = FALSE)
  fwrite(object@sampleData, file = paste0(file, "_sample_annotation.txt"), sep = "\t", na = "NA", quote = FALSE)
  fwrite(list(object@logs), file = paste0(file, "_logs.txt"), sep = "\t", na = "NA", quote = FALSE)
  }



#' Save metabolite data in an Excel file
#'
#'
#' @param object A Metabolite object
#' @param file Output file to save the metabolite measurements 
#' (suffixes: ".xlsx", "_logs.txt). 
#' @export
#' @return No return value.
#'
#'
save_data_xls <- function(object, file = "")  {
  # xlsx::write.xlsx(object@assayData, file = paste0(file, "xlsx"), sheetName = "Sheet1")
  # xlsx::write.xlsx(object@featureData, file = paste0(file, "xlsx"), sheetName = "Sheet2")
  # xlsx::write.xlsx(object@sampleData, file = paste0(file, "xlsx"), sheetName = "Sheet3")
  # 
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "assayData")
  openxlsx::addWorksheet(wb, "featureData")
  openxlsx::addWorksheet(wb, "sampleData")
  
  openxlsx::writeData(wb, "assayData", object@assayData)
  openxlsx::writeData(wb, "featureData", object@featureData)
  openxlsx::writeData(wb, "sampleData", object@sampleData)
  
  openxlsx::saveWorkbook(wb, file = paste0(file, ".xlsx"), overwrite = TRUE)
  
  fwrite(list(object@logs), file = paste0(file, "_logs.txt"), sep = "\t", na = "NA", quote = FALSE)
}


#' merge two Metabolite objects
#'
#' Merge two Metabolite objects.
#'
#' @param object_X The first Metabolite object.
#' @param object_Y The second Metabolite object.
#' @param all logical; all = TRUE: keep all metabolites; all = FALSE, 
#' keep common metabolites that were present in both datasets.
#' @param verbose print log information.
#' @return A Metabolite object after merging with slots: assayData, featureData, and sampleData.
#' @export
merge_data <- function(
  object_X = NULL,
  object_Y = NULL,
  all = TRUE,
  verbose = TRUE
)  {
  
  stopifnot(object_X@featureID == object_Y@featureID)
  stopifnot(object_X@sampleID == object_Y@sampleID)

  sample_IDs <- c(object_X@assayData[, get(object_X@sampleID)], object_Y@assayData[, get(object_Y@sampleID)])

  if(any(duplicated(sample_IDs))) {
    stop(paste0("\n Find ", sum(duplicated(sample_IDs)), " duplicated sampleID: ", 
                paste0(sample_IDs[duplicated(sample_IDs)][seq_len(5)], collapse = ", "),  " ... \n"))
  }


  merge_assayData <- rbind(object_X@assayData, object_Y@assayData, fill = TRUE)

  merge_featureData <- merge(object_X@featureData, object_Y@featureData, 
                             by.x = object_X@featureID, by.y = object_X@featureID, all = TRUE)

  v_feature_name <- setdiff(intersect(names(object_X@featureData), names(object_Y@featureData)), object_X@featureID)

  # for common names, impute.
  for(i in seq_along(v_feature_name)) {
    v_i <- v_feature_name[i]
    merge_featureData[, (v_i) := ifelse(is.na(get(paste0(v_i, ".x"))),
                                        get(paste0(v_i, ".y")),
                                        get(paste0(v_i, ".x")))]
    merge_featureData[, c(paste0(v_i, ".x"), paste0(v_i, ".y")) := NULL]
  }

  merge_sampleData <- rbind(object_X@sampleData, object_Y@sampleData, fill = TRUE)

    new(
    Class = 'Metabolite',
    assayData = merge_assayData,
    featureData = merge_featureData,
    sampleData = merge_sampleData,
    featureID = object_X@featureID,
    sampleID = object_X@sampleID,
    logs = paste0("******* object_X logs *******", "\n",
                  object_X@logs, "\n",
                  "******* object_Y logs *******", "\n",
                  object_Y@logs, "\n",
                  format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": Merge data, ", 
                  NROW(merge_sampleData), " samples and ", NROW(merge_featureData), " features.\n"),
    miscData = list()
  )

}




#' Update a Metabolite object
#'
#' Update a Metabolite object.
#' @param object A Metabolite object
#' @param dataset A vector or data.table used for a specific action mode.
#' @param action Currently support:
#' \itemize{
#' \item{"injection_order": `@sampleData` will be updated by the order of sampleID that provided 
#' in the injection order data}
#' \item{"keep_feature": feature ID list to keep}
#' \item{"remove_feature": feature ID list to remove}
#' \item{"keep_sample": sample ID list to keep}
#' \item{"remove_sample": sample ID list to remove}
#' \item{"add_sample_annotation": merge data with `@sampleData`}
#' \item{"change_featureID": change the name of featureID (provide the new column name in `@featureData` for dataset)}
#' }
#' @export
#' @return A Metabolite object after updating. 
update_Metabolite <- function(object,
                              dataset = NULL,
                              action = NULL
                              ) {

  df_data <- object@assayData
  df_feature <- object@featureData
  df_sample <- object@sampleData

  if(!is.null(action)) {

    if(action == "keep_feature") {
      if(! is.null(dataset)) {
        stopifnot(is.vector(dataset))
        df_data <- df_data[, unique(c(object@sampleID, dataset)), with = FALSE]
        df_feature <- df_feature[df_feature$featureID %in% dataset]
      } else stop("No dataset provided for action `keep_feature`.")
    }

    if(action == "remove_feature") {
      if(! is.null(dataset)) {
        stopifnot(is.vector(dataset))
        df_data <- df_data[, setdiff(names(df_data), dataset), with = FALSE]
        df_feature <- df_feature[! df_feature$featureID %in% dataset]
      } else stop("No dataset provided for action `remove_feature`.")
    }


    if(action == "keep_sample") {
      if(! is.null(dataset)) {
        stopifnot(is.vector(dataset))
        df_data <- df_data[get(object@sampleID) %in% dataset]
        df_sample <- df_sample[get(object@sampleID) %in% dataset]
      } else stop("No dataset provided for action `keep_sample`.")
    }


    if(action == "remove_sample") {
      if(! is.null(dataset)) {
        stopifnot(is.vector(dataset))
        df_data <- df_data[! get(object@sampleID) %in% dataset]
        df_sample <- df_sample[! get(object@sampleID) %in% dataset]
      } else stop("No dataset provided for action `remove_sample`.")
    }

    if(action == "injection_order") {
      if(! is.null(dataset)) {

        stopifnot(is.data.frame(dataset))
        stopifnot(object@sampleID %in% names(dataset))
        stopifnot(all(df_data[, get(object@sampleID)] %in% dataset[, get(object@sampleID)]))

        df_sample <- merge(dataset[, object@sampleID, with = FALSE], df_sample, by = object@sampleID, sort = FALSE)
        df_data <- merge(dataset[, object@sampleID, with = FALSE], df_data, by = object@sampleID, sort = FALSE)

        if(! "ID_injection_order" %in% names(df_sample)) {

          cat("Creat a new column `ID_injection_order`. \n")
        } else {
          cat("`ID_injection_order` exist, will be overwritten. \n")
        }
        df_sample$ID_injection_order <- seq_len(NROW(df_sample))

        logs = paste0(object@logs, format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": Update injection order. \n")
      } else stop("No dataset provided for action `injection_order`.")
    }

    if(action == "add_sample_annotation") {
      if(! is.null(dataset)) {
        stopifnot(is.data.frame(dataset))
        df_sample_merge <- merge(df_sample, dataset, by =  object@sampleID, sort = FALSE, all.x = TRUE)
        v_sample_name <- setdiff(intersect(names(df_sample), names(dataset)), object@sampleID)

        for(i in seq_along(v_sample_name)) {
          v_i <- v_sample_name[i]
          df_sample_merge[, (v_i) := ifelse(is.na(get(paste0(v_i, ".x"))),
                                              get(paste0(v_i, ".y")),
                                              get(paste0(v_i, ".x")))]
          df_sample_merge[, c(paste0(v_i, ".x"), paste0(v_i, ".y")) := NULL]
        }
        df_sample <- copy(df_sample_merge)
      } else stop("No dataset provided for action `remove_sample`.")
    }

    if(action == "change_featureID") {

      if(! is.null(dataset)) {
        if(!length(dataset) == 1) stop("dataset should be a column name in `@featureData`, 
                                       create the column first if it does not exist.")


        setnames(df_data, df_feature$featureID, df_feature[, get(dataset)])
        df_feature$featureID <- df_feature[, get(dataset)]

      } else stop("No dataset provided for action `change_featureID`.")
    }
  }

  if(is.null(action)) {
    logs = paste0(object@logs, format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": Current data, ", 
                  NROW(df_sample), " samples and ", NROW(df_feature), " features. \n")
  } else {
    logs = paste0(object@logs, format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": Update data, action: ",
                  action, ", ", NROW(df_sample), " samples and ", NROW(df_feature), " features. \n")
  }
  

  new(
    Class = 'Metabolite',
    assayData = df_data,
    featureData = df_feature,
    sampleData = df_sample,
    featureID = object@featureID,
    sampleID = object@sampleID,
    logs = logs,
    miscData = object@miscData
  )
  
}






#' Convert maplet data to a Metabolite object
#' 
#' @param object Data in the maplet package, a SummarizedExperiment data object.
#' @param featureID feature ID column name in the maplet data (-> elementMetadata -> listData). 
#' @param sampleID sample ID column name in the maplet data (-> colData -> listData). 
#' @return A Metabolite object.
#' @export
maplet2Metabolite <- function(object, featureID = "name", sampleID = "sample") {
  
  df_sample <- object@colData@listData
  df_sample <- data.table(as.data.frame(df_sample))
  setnames(df_sample, sampleID, "sampleID")
  
  df_feature <- object@elementMetadata@listData
  df_feature <- data.table(as.data.frame(df_feature))
  setnames(df_feature, featureID, "featureID")
  
  df_data <- object@assays@data@listData[[1]]
  df_data  <- data.table(t(df_data), keep.rownames = TRUE)
  # rename df_data
  setnames(df_data, c("sampleID", df_feature$featureID))
  
  logs <- paste0(format(Sys.time(), "%d/%m/%y %H:%M:%OS"), ": Convert data from maplet data.\n")
  object_convert <- create_Metabolite(assayData = df_data,
                              featureData = df_feature,
                              sampleData = df_sample,
                              featureID = "featureID",
                              sampleID = "sampleID",
                              logs = logs, 
                              miscData = object@metadata
  )
  
  return(object_convert)
}





#' Convert a Metabolite object to maplet data
#' 
#' @param object A Metabolite object
#' @return Data in the maplet package, a SummarizedExperiment data object.
#' @export
Metabolite2maplet <- function(object) {
  save_data_xls(object, file = "tmp_data_convert")
  file_data <- "tmp_data_convert.xlsx"
  
  object_convert <-
    mt_load_checksum(file=file_data) %>%
    mt_load_xls(file=file_data, sheet="assayData", samples_in_row=T, id_col="sampleID") %>%
    mt_anno_xls(file=file_data, sheet="sampleData", anno_type="samples", anno_id_col ="sampleID") %>%
    mt_anno_xls(file=file_data, sheet="featureData",anno_type="features", anno_id_col="featureID", data_id_col = "name") 
  
  file.remove(file_data)
  file.remove("tmp_data_convert_logs.txt")
  
  return(object_convert)
  
}
