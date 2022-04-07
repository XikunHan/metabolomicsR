## ----options, include=FALSE, echo=FALSE---------------------------------------

library(BiocStyle)
knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE)


## -----------------------------------------------------------------------------

library(metabolomicsR)
library(data.table)
library(ggplot2)
library(cowplot)
library(plotROC)

# Load the dataset

file_path <- system.file("extdata", "QMDiab_metabolomics_OrigScale.xlsx", package = "metabolomicsR", mustWork = TRUE)

df_plasma <- load_excel(path = file_path,
                        data_sheet = 1,
                        feature_sheet = 4,
                        sample_sheet = 8,
                        sampleID = "QMDiab-ID",
                        featureID = "BIOCHEMICAL"
                        )

# To save the data
# save_data(df_plasma, file = "~/test")


## -----------------------------------------------------------------------------
df_plasma

## -----------------------------------------------------------------------------
# change the feature ID to the column `COMP_IDstr`
df_plasma <- update_Metabolite(df_plasma, dataset = "COMP_IDstr", action = "change_featureID")


## -----------------------------------------------------------------------------
df_plasma

## -----------------------------------------------------------------------------
# load urine metabolomic data
df_urine <- load_excel(path = file_path,
                        data_sheet = 2,
                        feature_sheet = 5,
                        sample_sheet = 9,
                        sampleID = "QMDiab-ID",
                        featureID = "BIOCHEMICAL"
                        )

df_urine <- update_Metabolite(df_urine, dataset = "COMP_IDstr", action = "change_featureID")

## -----------------------------------------------------------------------------
df_urine

## -----------------------------------------------------------------------------
# load saliva metatabolomic data
df_saliva <- load_excel(path = file_path,
                        data_sheet = 3,
                        feature_sheet = 6,
                        sample_sheet = 10,
                        sampleID = "QMDiab-ID",
                        featureID = "BIOCHEMICAL"
                        )

df_saliva <- update_Metabolite(df_saliva, dataset = "COMP_IDstr", action = "change_featureID")

## -----------------------------------------------------------------------------
df_saliva

## -----------------------------------------------------------------------------
df_plasma_QC <- QC_pipeline(df_plasma, replace_outlier_method = "winsorize", impute_method = NULL)

## -----------------------------------------------------------------------------
df_plasma_QC

## -----------------------------------------------------------------------------
df_urine_QC <- QC_pipeline(df_urine, replace_outlier_method = "winsorize", impute_method = NULL)

## -----------------------------------------------------------------------------
df_urine_QC

## -----------------------------------------------------------------------------
df_saliva_QC <- QC_pipeline(df_saliva, replace_outlier_method = "winsorize", impute_method = NULL)

## -----------------------------------------------------------------------------
df_saliva_QC

## -----------------------------------------------------------------------------
df_plasma_QC <- impute(df_plasma_QC, method = "half-min")
df_plasma_scale <-  transformation(df_plasma_QC, method = "log")
df_plasma_scale <-  transformation(df_plasma_scale, method = "scale")

## -----------------------------------------------------------------------------
df_plasma_scale

## ----  fig.width = 12, fig.height = 7-----------------------------------------
# if no features were selected, randomly show 16 metabolites
plot_Metabolite(df_plasma_QC, plot = "boxplot", x = "T2D", color ="ETHNICITY", shape = "T2D")


## ---- fig.width = 12, fig.height = 6------------------------------------------
# select three metabolites
plot_Metabolite(df_plasma_QC, x = "T2D", plot = "boxplot", feature_name = c("M43027",  "M11953", "M38002"))

## ---- fig.width = 12, fig.height = 7------------------------------------------
# comparisons between groups using `ggbetweenstats`
plot_Metabolite(df_plasma_QC, x = "T2D", plot = "betweenstats",  feature_name = c("M43027",  "M11953"))


## ---- fig.width = 12, fig.height = 7------------------------------------------
plot_Metabolite(df_plasma_scale, plot = "boxplot", x = "T2D", color ="ETHNICITY", shape = "T2D")


## ---- fig.width = 12, fig.height = 7------------------------------------------

plot_Metabolite(df_plasma_scale, plot = "histogram", color = "T2D")

plot_Metabolite(df_plasma_scale, plot = "histogram", color = "ETHNICITY")


## -----------------------------------------------------------------------------
# Take the urine data for example, the batch of the samples are missing. To illustrate the usage of the normalization method, we assume that 120 samples in each batch. 
unique(featureData(df_urine_QC)$PLATFORM)

# add new columns by reference (column name uppercase)
sampleData(df_urine_QC)[, `GC/MS` := rep(1:3, each = 120)[1:359]]
sampleData(df_urine_QC)[, `LC/MS NEG` := rep(1:3, each = 120)[1:359]]
sampleData(df_urine_QC)[, `LC/MS POS` := rep(1:3, each = 120)[1:359]]

# we also added an injection order column (in the current order)
df_urine_QC <- update_Metabolite(df_urine_QC, sampleData(df_urine_QC)[, 1], action = "injection_order")

# and set one QC sample for each 10 samples. 

sampleData(df_urine_QC)[, QCsample := ifelse(1:359%%10 == 0, 1, 0)]
sampleData(df_urine_QC)[, QCsample := factor(QCsample, levels = c(1, 0))]
                       
sampleData(df_urine_QC)[, sampleID := ifelse(1:359 %% 10 == 0, paste0("MTRX-", sampleID), sampleID)]

table(sampleData(df_urine_QC)$QCsample)

assayData(df_urine_QC)$sampleID <- sampleData(df_urine_QC)$sampleID # also change the ID in assayData
df_urine_QC


## ---- fig.width = 12, fig.height = 7------------------------------------------
v_features <- featureData(df_urine_QC)$featureID[1:16]

plot_injection_order(df_urine_QC, color = "QCsample", shape = "GC/MS", feature_name = v_features)

## ---- fig.width = 12, fig.height = 7------------------------------------------
df_urine_QC_norm <- batch_norm(df_urine_QC)
p1 <- plot_injection_order(df_urine_QC_norm, color = "QCsample", shape = "GC/MS", feature_name = v_features)
p1

## ---- fig.width = 12, fig.height = 7------------------------------------------
df_urine_QC_norm <- QCmatrix_norm(df_urine_QC)
p2 <- plot_injection_order(df_urine_QC_norm, color = "QCsample", shape = "GC/MS", feature_name = v_features)
p2

## ---- fig.width = 12, fig.height = 7------------------------------------------
df_urine_QC_norm <- nearestQC_norm(df_urine_QC)
p3 <- plot_injection_order(df_urine_QC_norm, color = "QCsample", shape = "GC/MS", feature_name = v_features)
p3

## ---- fig.width = 12, fig.height = 7------------------------------------------
df_urine_QC_norm_LOESS <- modelling_norm(df_urine_QC, method = "LOESS")
p4 <- plot_injection_order(df_urine_QC_norm, color = "QCsample", shape = "GC/MS", feature_name = v_features)
p4

# to benchmark the performance
df_urine_QC_norm_KNN <- modelling_norm(df_urine_QC, method = "KNN", k = 3) 
p5 <- plot_injection_order(df_urine_QC_norm_KNN, color = "QCsample", shape = "GC/MS", feature_name = v_features)
p5


# to benchmark the performance
df_urine_QC_norm_xgboost <- modelling_norm(df_urine_QC, method = "XGBoost") 
p6 <- plot_injection_order(df_urine_QC_norm_xgboost, color = "QCsample", shape = "GC/MS", feature_name = v_features)
p6


## -----------------------------------------------------------------------------

df_plasma_PCA <- run_PCA(df_plasma_QC)

plot_PCA(df_plasma_PCA, color ="ETHNICITY", shape = "T2D")

plot_UMAP(df_plasma_QC, color ="ETHNICITY", shape = "T2D")

plot_tsne(df_plasma_QC, color ="ETHNICITY", shape = "T2D")


## ---- fig.width = 10, fig.height = 7------------------------------------------

# pairwise correlation of metabolites between two different fluids

# plasma vs urine
dd <- correlation(df_plasma_QC, df_urine_QC, method = "spearman")
dd <- merge(dd, featureData(df_plasma_QC), by.x = "term", by.y = "featureID")

p <- ggplot(dd, aes(x = SUPER_PATHWAY, y = r, fill = SUPER_PATHWAY)) +
  geom_boxplot() +
  geom_jitter( size=0.4, alpha=0.9) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Set1")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  labs(title = "plasma vs urine")

p


## ---- fig.width = 14, fig.height = 7------------------------------------------
# plasma vs saliva
dd <- correlation(df_plasma_QC, df_saliva_QC, method = "spearman")
dd <- merge(dd, featureData(df_plasma_QC), by.x = "term", by.y = "featureID")

p2 <- ggplot(dd, aes(x = SUPER_PATHWAY, y = r, fill = SUPER_PATHWAY)) +
  geom_boxplot() +
  geom_jitter(size=0.4, alpha=0.9) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Set1")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "plasma vs saliva")


# urine vs saliva
dd <- correlation(df_urine_QC, df_saliva_QC, method = "spearman")
dd <- merge(dd, featureData(df_urine_QC), by.x = "term", by.y = "featureID")

p3 <- ggplot(dd, aes(x = SUPER_PATHWAY, y = r, fill = SUPER_PATHWAY)) +
  geom_boxplot() +
  geom_jitter(size=0.4, alpha=0.9) +
   scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Set1")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "urine vs saliva")


p <- cowplot::plot_grid(p2, p3, nrow = 1)
p

## -----------------------------------------------------------------------------
fit_lm <- regression(object = df_plasma_scale, phenoData = NULL, model = "lm", outcome = "BMI",
                          covars = c("AGE", "GENDER", "ETHNICITY"), factors = "ETHNICITY")

head(fit_lm)


## -----------------------------------------------------------------------------

dd <- merge(fit_lm, featureData(df_plasma_scale), by.x = "term", by.y = "featureID")

dd[, sig := ifelse(p.value.adj < 0.1, 1, 0)]
plot_volcano(dd, color = NULL, label = "BIOCHEMICAL")


## -----------------------------------------------------------------------------

fit_glm <- regression(object = df_plasma_scale, phenoData = NULL, model = "logistic", outcome = "T2D",
                         covars = c("AGE","GENDER", "ETHNICITY"), factors = "ETHNICITY")

head(fit_glm)

## -----------------------------------------------------------------------------
dd <- merge(fit_glm, featureData(df_plasma_scale), by.x = "term", by.y = "featureID")

dd[, sig := ifelse(p.value.adj < 0.1, 1, 0)]

plot_volcano(dd, label = "BIOCHEMICAL")


## -----------------------------------------------------------------------------

fit2 <- regression(object = df_plasma_scale, phenoData = NULL, model = "auto", outcome = c("BMI", "T2D"),
                         covars = c("AGE","GENDER", "ETHNICITY"), factors = "ETHNICITY")

head(fit2)

## -----------------------------------------------------------------------------

dd2 <- merge(fit2, featureData(df_plasma_scale), by.x = "term", by.y = "featureID")

dd2[, sig := ifelse(p.value.adj < 0.1, 1, 0)]
plot_volcano(dd2, color = "outcome", label = "BIOCHEMICAL")


## -----------------------------------------------------------------------------

# for example, to fit a Negative Binomial Generalized Linear Model, we first define a `fit_glm.nb` function. 

library(MASS)

fit_glm.nb <- function(data = NULL, formula = NULL, keep = NULL) {
  v_var <- all.vars(formula)
  df <- data[, v_var, with = FALSE]
  fit <- tryCatch(
    do.call("glm.nb", args = list(data = df, formula = formula)),
    error = function(e) {
      cat(paste0("Failed to fit model: ", e), "\n")
    })
  if(inherits(fit, "glm")) {
    res  <- as.data.table(broom::tidy(fit))
    if(! is.null(keep)) res <- res[res$term %in% keep, ]
    res$n <- length(fit$residuals)
  } else {
    res <- list(estimate = NA)
  }
  return(res)
}

# the outcome 'toy_y' was only created to demonstrate the usage to extend regression models

sampleData(df_plasma_scale)[, toy_y := round(BMI, 0)]

# equivalent to:
if(0) {
  dd <- sampleData(df_plasma_scale) 
  dd[, toy_y := round(BMI, 0)]
  sampleData(df_plasma_scale)  <- dd
}

fit <- regression(object = df_plasma_scale, phenoData = NULL, model = "glm.nb", outcome = "toy_y",
                  covars = c("AGE","GENDER", "ETHNICITY"), factors = "ETHNICITY", 
                  feature_name = featureData(df_plasma_scale)$featureID[1:50])


## -----------------------------------------------------------------------------
dd <- merge(fit, featureData(df_plasma_scale), by.x = "term", by.y = "featureID")
dd[, sig := ifelse(p.value.adj < 0.1, 1, 0)]
plot_volcano(dd, color = NULL, label = "BIOCHEMICAL")

## ---- fig.height=8,fig.width=9------------------------------------------------

# ROC for one variable, to add a wrapped for many metabolites separately
plot_ROC(object = df_plasma_scale, y = "T2D", x = "M00527")
plot_ROC(object = df_plasma_scale, y = "T2D", x = c("M00527", "BMI"))


## ---- fig.height=8,fig.width=9------------------------------------------------

plot_ROC(object = df_plasma_scale, y = "T2D", 
         model_a = c("BMI", "AGE", "GENDER"), 
         model_b = c("M00527", "BMI", "AGE", "GENDER"))


## -----------------------------------------------------------------------------
sessionInfo()

