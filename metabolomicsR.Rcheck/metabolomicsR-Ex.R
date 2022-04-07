pkgname <- "metabolomicsR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('metabolomicsR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("QCmatrix_norm")
### * QCmatrix_norm

flush(stderr()); flush(stdout())

### Name: QCmatrix_norm
### Title: QCmatrix normalization
### Aliases: QCmatrix_norm

### ** Examples

## Not run: 
##D d <- QCmatrix_norm(object = df)
## End(Not run)



cleanEx()
nameEx("RSD")
### * RSD

flush(stderr()); flush(stdout())

### Name: RSD
### Title: RSD
### Aliases: RSD

### ** Examples

## Not run: 
##D v <- RSD(x)
## End(Not run)



cleanEx()
nameEx("column_missing_rate")
### * column_missing_rate

flush(stderr()); flush(stdout())

### Name: column_missing_rate
### Title: column missing rate
### Aliases: column_missing_rate column_missing_rate.default
###   column_missing_rate.Metabolite

### ** Examples

## Not run: 
##D # for a data.frame or data.table
##D v <- column_missing_rate(object = df)
##D 
##D # to skip the first column (eg. ID)
##D v <- column_missing_rate(object = df[, -1])
## End(Not run)

## Not run: 
##D 
##D # for a Metabolite object
##D v <- column_missing_rate(object)
##D 
## End(Not run)




cleanEx()
nameEx("create_Metabolite")
### * create_Metabolite

flush(stderr()); flush(stdout())

### Name: create_Metabolite
### Title: Create a Metabolite object
### Aliases: create_Metabolite

### ** Examples

## Not run: 
##D 
##D df <- create_Metabolite(assayData = df_data, featureData = df_feature, sampleData =  df_sample)
##D 
## End(Not run)




cleanEx()
nameEx("filter_column_constant")
### * filter_column_constant

flush(stderr()); flush(stdout())

### Name: filter_column_constant
### Title: filter columns if values are constant
### Aliases: filter_column_constant filter_column_constant.default
###   filter_column_constant.Metabolite

### ** Examples

## Not run: 
##D 
##D # for a data.frame or data.table
##D v <- filter_column_missing_rate(object = df)
##D 
##D # if skip the first column (eg. ID)
##D v <- filter_column_missing_rate(object = df[, -1])
##D 
## End(Not run)




cleanEx()
nameEx("filter_column_missing_rate")
### * filter_column_missing_rate

flush(stderr()); flush(stdout())

### Name: filter_column_missing_rate
### Title: filter columns using missing rate
### Aliases: filter_column_missing_rate filter_column_missing_rate.default
###   filter_column_missing_rate.Metabolite

### ** Examples

## Not run: 
##D 
##D d <- filter_column_missing_rate(object)
##D 
## End(Not run)




cleanEx()
nameEx("filter_row_missing_rate")
### * filter_row_missing_rate

flush(stderr()); flush(stdout())

### Name: filter_row_missing_rate
### Title: filter rows using missing rate
### Aliases: filter_row_missing_rate filter_row_missing_rate.default
###   filter_row_missing_rate.Metabolite

### ** Examples

## Not run: 
##D 
##D d <- filter_row_missing_rate(object)
##D 
## End(Not run)




cleanEx()
nameEx("genuMet_makefeature")
### * genuMet_makefeature

flush(stderr()); flush(stdout())

### Name: genuMet_makefeature
### Title: distinguish genuine untargeted metabolic features without QC
###   samples
### Aliases: genuMet_makefeature

### ** Examples

## Not run: 
##D v <- genuMet_makefeature(df)
## End(Not run)




cleanEx()
nameEx("impute")
### * impute

flush(stderr()); flush(stdout())

### Name: impute
### Title: impute missing values
### Aliases: impute impute.Metabolite impute.default impute_kNN

### ** Examples



## Not run: 
##D d <- impute(object)
## End(Not run)




cleanEx()
nameEx("inverse_rank_transform")
### * inverse_rank_transform

flush(stderr()); flush(stdout())

### Name: inverse_rank_transform
### Title: rank-based inverse normal transformation
### Aliases: inverse_rank_transform

### ** Examples


## Not run: 
##D v <- inverse_rank_transform(x)
## End(Not run)




cleanEx()
nameEx("is_outlier")
### * is_outlier

flush(stderr()); flush(stdout())

### Name: is_outlier
### Title: is outlier
### Aliases: is_outlier

### ** Examples


## Not run: 
##D v <- is_outlier(x)
## End(Not run)




cleanEx()
nameEx("load_excel")
### * load_excel

flush(stderr()); flush(stdout())

### Name: load_excel
### Title: Load metabolite data from an excel file
### Aliases: load_excel

### ** Examples


file_path <- system.file("extdata", "QMDiab_metabolomics_OrigScale.xlsx", 
package = "metabolomicsR", mustWork = TRUE)

df_plasma <- load_excel(path = file_path, data_sheet = 1, feature_sheet = 4, sample_sheet = 8, 
sampleID = "QMDiab-ID", featureID = "BIOCHEMICAL")




cleanEx()
nameEx("merge_data")
### * merge_data

flush(stderr()); flush(stdout())

### Name: merge_data
### Title: merge two Metabolite objects
### Aliases: merge_data

### ** Examples

# to merge two Metabolite objects
# df <- merge_data(df_plasma, df_plasma)




cleanEx()
nameEx("modelling_norm")
### * modelling_norm

flush(stderr()); flush(stdout())

### Name: modelling_norm
### Title: LOESS normalization
### Aliases: modelling_norm

### ** Examples

## Not run: 
##D d <- QCmatrix_norm(object = df)
## End(Not run)



cleanEx()
nameEx("nearestQC_norm")
### * nearestQC_norm

flush(stderr()); flush(stdout())

### Name: nearestQC_norm
### Title: nearest QC sample normalization
### Aliases: nearestQC_norm

### ** Examples

## Not run: 
##D d <- nearestQC_norm(object = df)
## End(Not run)



cleanEx()
nameEx("outlier_rate")
### * outlier_rate

flush(stderr()); flush(stdout())

### Name: outlier_rate
### Title: outlier rate
### Aliases: outlier_rate outlier_rate.default outlier_rate.data.frame
###   outlier_rate.Metabolite

### ** Examples


## Not run: 
##D v <- outlier_rate(x)
## End(Not run)

## Not run: 
##D 
##D # for a Metabolite object
##D v <- outlier_rate(object)
##D 
## End(Not run)




cleanEx()
nameEx("pareto_scale")
### * pareto_scale

flush(stderr()); flush(stdout())

### Name: pareto_scale
### Title: pareto scale transformation
### Aliases: pareto_scale

### ** Examples


## Not run: 
##D v <- paretoscale(x)
## End(Not run)




cleanEx()
nameEx("plot_Metabolite")
### * plot_Metabolite

flush(stderr()); flush(stdout())

### Name: plot_Metabolite
### Title: plot a Metabolite object
### Aliases: plot_Metabolite

### ** Examples

## Not run: 
##D p <- plot_Metabolite(df_m_PCA, plot = "boxplot")
##D p
## End(Not run)




cleanEx()
nameEx("plot_UMAP")
### * plot_UMAP

flush(stderr()); flush(stdout())

### Name: plot_UMAP
### Title: Plot UMAP
### Aliases: plot_UMAP

### ** Examples

## Not run: 
##D p<- plot_UMAP(df_2017_QC_norm_PCA, color = "NEG", shape = "QCSample")
##D p
## End(Not run)



cleanEx()
nameEx("plot_injection_order")
### * plot_injection_order

flush(stderr()); flush(stdout())

### Name: plot_injection_order
### Title: injection order scatterplot
### Aliases: plot_injection_order

### ** Examples

## Not run: 
##D p <- plot_injection_order(df_m_PCA, color = "QC_sample")
##D p
##D 
##D p <- plot_injection_order(df_m_PCA, color = "QC_sample", feature_name = "X563")
##D p
## End(Not run)




cleanEx()
nameEx("plot_tsne")
### * plot_tsne

flush(stderr()); flush(stdout())

### Name: plot_tsne
### Title: plot tSNE
### Aliases: plot_tsne

### ** Examples

## Not run: 
##D p<- plot_tsne(df_2017_QC_norm_PCA, color = "NEG", shape = "QCSample")
##D p
##D 
## End(Not run)



cleanEx()
nameEx("plot_volcano")
### * plot_volcano

flush(stderr()); flush(stdout())

### Name: plot_volcano
### Title: volcano plot for regression results
### Aliases: plot_volcano

### ** Examples

## Not run: 
##D p <- plot_volcano(fit_lm, color = NULL)
##D p
## End(Not run)




cleanEx()
nameEx("regression")
### * regression

flush(stderr()); flush(stdout())

### Name: regression
### Title: regression analysis
### Aliases: regression regression_each

### ** Examples

data(df_plasma)
fit_lm <- regression(object = df_plasma, phenoData = NULL, model = "lm", 
outcome = "BMI", covars = c("AGE", "GENDER", "ETHNICITY"), factors = "ETHNICITY")




cleanEx()
nameEx("replace_outlier")
### * replace_outlier

flush(stderr()); flush(stdout())

### Name: replace_outlier
### Title: change outlier values as NA or winsorize
### Aliases: replace_outlier replace_outlier.default
###   replace_outlier.data.frame replace_outlier.Metabolite

### ** Examples



## Not run: 
##D d <- replace_outlier(object, method = "winsorize", nSD = 5)
## End(Not run)



## Not run: 
##D d <- replace_outlier(object, method = "winsorize", nSD = 5)
## End(Not run)




cleanEx()
nameEx("row_missing_rate")
### * row_missing_rate

flush(stderr()); flush(stdout())

### Name: row_missing_rate
### Title: row missing rate
### Aliases: row_missing_rate row_missing_rate.default
###   row_missing_rate.Metabolite

### ** Examples

## Not run: 
##D 
##D # for a data.frame or data.table
##D v <- row_missing_rate(object = df)
##D 
##D # to skip the first column (eg. ID)
##D v <- row_missing_rate(object = df[, -1])
## End(Not run)

## Not run: 
##D 
##D # for a Metabolite object
##D v <- row_missing_rate(object)
##D 
## End(Not run)




cleanEx()
nameEx("run_PCA")
### * run_PCA

flush(stderr()); flush(stdout())

### Name: run_PCA
### Title: Principal Components Analysis
### Aliases: run_PCA

### ** Examples


# skip the first column (eg. ID) to impute missing values
## Not run: 
##D d <- run_PCA(object)
## End(Not run)




cleanEx()
nameEx("transformation")
### * transformation

flush(stderr()); flush(stdout())

### Name: transformation
### Title: apply transformation to a Metabolite object
### Aliases: transformation

### ** Examples


## Not run: 
##D d <- transformation(x)
## End(Not run)




cleanEx()
nameEx("update_Metabolite")
### * update_Metabolite

flush(stderr()); flush(stdout())

### Name: update_Metabolite
### Title: Update a Metabolite object
### Aliases: update_Metabolite

### ** Examples

# df_plasma <- update_Metabolite(df_plasma, dataset = "COMP_IDstr", action = "change_featureID")




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
