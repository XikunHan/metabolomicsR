#' regression analysis
#'
#' Run regression models with adjusting for covariates. `regression_each` is used for one outcome. In `regression`, several outcomes can be specified to run together.
#' @param object A Metabolite object.
#' @param phenoData A data.table with outcome and covariates. If `phenoData` is NULL, `@sampleData` will be used.
#' @param model Specify a regression model. See \code{\link{fit_lm}} for more details. 'auto' can be used to infer 'lm' or 'logistic' (with only 0, 1, and NA).
#' @param formula A character or formula object to fit model (only used in `regression_each`)
#' @param outcome Column name of the outcome variable.
#' @param covars Column names of covariates.
#' @param factors Variables to be treated as factor.
#' @param feature_name A vector of selected metabolites to run. If both feature_name and random_select are NULL, will run regression for all features.
#' @param time Column name of survival time, used in cox regression, see \code{\link[survival]{coxph}} for more details.
#' @param verbose Print log information.
#' @param ncpus Number of CPUS for parallele job.
#' @param p.adjust.method Adjust for P value method, see \code{\link{p.adjust}}.
#' @param \dots Further arguments passed to regression model.
#' @returns term estimate std.error statistic p.value n outcome p.value.adj.
#' @export
#' @examples
#' data(df_plasma)
#' fit_lm <- regression(object = df_plasma, phenoData = NULL, model = "lm", 
#' outcome = "BMI", covars = c("AGE", "GENDER", "ETHNICITY"), factors = "ETHNICITY")
#'
regression <- function(object, phenoData = NULL, model = NULL, outcome = NULL, covars = NULL, factors = NULL, feature_name = NULL, time = NULL, verbose = TRUE, ncpus = 1, p.adjust.method = "bonferroni", ...) {

  cat(paste0("Regression for ", length(outcome) ," outcome(s). \n\n"))
  res_m <- NULL
  for(i in seq_along(outcome)) {
    v_i <- outcome[i]
    res <- regression_each(object = object, phenoData = phenoData, model = model, outcome = v_i,
                      covars = covars, factors =factors, feature_name = feature_name,
                      time = time, verbose = verbose, ncpus = ncpus, p.adjust.method = p.adjust.method, ...
                      )
    if(NROW(res) == 0 || length(res) == 1) next
    res_m <- rbind(res_m, res, fill = TRUE)
  }
  return(res_m)
}



#' @title available regression methods
#' @description `fit_lm`: linear regression model \code{\link{lm}}.
#'
#' @param data A data.table with all variables to be fitted.
#' @param formula A "formula" object to be fitted.
#' @param keep Variables to keep regression results.
#' @seealso  \code{\link{regression}}
#' @export
#' @returns term estimate std.error statistic p.value n.
fit_lm <- function(data = NULL, formula = NULL, keep = NULL) {
  v_var <- all.vars(formula)
  df <- data[, v_var, with = FALSE]
  fit <- tryCatch(
    do.call("lm", args = list(data = df, formula = formula)),
    error = function(e) {
      cat(paste0("Failed to fit model: ", e), "\n")
    })
  if(inherits(fit, "lm")) {
    res  <- as.data.table(broom::tidy(fit))
    if(! is.null(keep))  res <- res[res$term %in% keep, ]
    res$n <- length(fit$residuals)
  } else {
    res <- list(estimate = NA)
  }
  return(res)
}


#' @description `fit_logistic`: logistic regression model \code{\link{glm}}.
#' @rdname fit_lm
#' @export
fit_logistic <- function(data = NULL, formula = NULL, keep = NULL) {
  v_var <- all.vars(formula)
  df <- data[, v_var, with = FALSE]
  fit <- tryCatch(
    do.call("glm", args = list(data = df, formula = formula,  family = binomial(link = "logit"))),
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



#' @description `fit_poisson`: poisson regression model \code{\link{glm}}.
#' @rdname fit_lm
#' @export
fit_poisson <- function(data = NULL, formula = NULL, keep = NULL) {
  v_var <- all.vars(formula)
  df <- data[, v_var, with = FALSE]
  fit <- tryCatch(
    do.call("glm", args = list(data = df, formula = formula,  family="poisson")),
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


#' @description `fit_cox`: proportional hazards regression model \code{\link[survival]{coxph}}.
#' @rdname fit_lm
#' @export
fit_cox <- function(data = NULL, formula = NULL, keep = NULL) {
  check_pkg("survival")
  
  v_var <- all.vars(formula)
  df <- data[, v_var, with = FALSE]
  fit <- tryCatch(
    do.call(getExportedValue("survival","coxph"), args = list(data = df, formula = formula)),
    error = function(e) {
      cat(paste0("Failed to fit model: ", e), "\n")
    })
  
  if(inherits(fit, "coxph")) {
    res  <- as.data.table(broom::tidy(fit))
    if(! is.null(keep)) res <- res[res$term %in% keep, ]
    res$n <- length(fit$residuals)
  } else {
    res <- list(estimate = NA)
  }
  return(res)
}



#' @description `fit_lme`: linear mixed-effects model \code{\link[nlme]{lme}}.
#' @param \dots Further arguments passed to regression model.
#' @rdname fit_lm
#' @export
fit_lme <- function(data = NULL, formula = NULL, keep = NULL, ...) {
  check_pkg("nlme")
  
  v_var <- all.vars(formula)
  
  v_args <- list(...)
  
  if (! "random" %in% names(v_args)) {
    warnings("`random` formula should be specified for `lme`.")
  } else {
    random <- v_args[["random"]]
  }
  
  if (! "na.action" %in% names(v_args)) {
    na.action <- na.omit
  } else {
    na.action <- v_args[["na.action"]]
  }
  
  if (! "control" %in% names(v_args)) {
    control  <- nlme::lmeControl(opt = "optim")
  } else {
    control <- v_args[["control"]]
  }
  
  v_args <- v_args[! names(v_args) %in% c("random", "na.action", "control", "keep.data")]
  
  if(length(v_args) != 0) {
    warning(paste0("Unknown args in lme: ", paste0(names(v_args), collapse = ", "),"\n"))
  }
  
  
  # select variables.
  df <- data[, c(v_var, all.vars(random)), with = FALSE]
  
  
  
  fit <- tryCatch(
    do.call(getExportedValue("nlme", "lme"), args = list(fixed = formula, data = df,
                                                         random = random,
                                                         na.action = na.action,
                                                         control = control,
                                                         keep.data = FALSE
    )),
    error = function(e) {
      cat(paste0("Failed to fit model: ", e), "\n")
    })
  
  if(inherits(fit, "lme")) {
    res <- as.data.table(summary(fit)$tTable, keep.rownames = TRUE)
    res <- res[, list(term = res$rn, estimate = res$Value, std.error = res$Std.Error,  statistic = res$`t-value`, p.value = res$`p-value`,  n = fit$dims[["N"]])]
    if(! is.null(keep)) res <- res[res$term %in% keep, ]
  } else {
    res <- list(estimate = NA)
  }
  return(res)
}





#' @description `fit_glmer`: logistic linear mixed-effects model \code{\link[lme4]{glmer}}.
#' @param \dots Further arguments passed to regression model.
#' @rdname fit_lm
#' @export
fit_glmer <- function(data = NULL, formula = NULL, keep = NULL, ...) {
  
  
  check_pkg("lme4")
  
  v_var <- all.vars(formula)
  
  v_args <- list(...)
  
  if ("random" %in% names(v_args)) {
    random <- v_args[["random"]]
    formula <- as.formula(paste0(deparse(formula)," ",  gsub("~", " + ", deparse(random))))
    v_var <- c(v_var, all.vars(random))
  }
  
  if (! "nAGQ" %in% names(v_args)) {
    nAGQ  <- 25
  } else {
    nAGQ <- v_args[["nAGQ"]]
  }
  
  if (! "control" %in% names(v_args)) {
    control = lme4::glmerControl()
  } else {
    control <- v_args[["control"]]
  }
  
  
  v_args <- v_args[! names(v_args) %in% c("random", "nAGQ", "control")]
  
  if(length(v_args) != 0) {
    warning(paste0("Unknown args in glmer: ", paste0(names(v_args), collapse = ", "),"\n"))
  }
  
  
  # select variables.
  df <- data[, v_var, with = FALSE]
  
  fit <- tryCatch(
    do.call(getExportedValue("lme4", "glmer"), args = list(formula = formula, data = df,
                                                           family = binomial,
                                                           nAGQ = nAGQ,
                                                           control = control
    )),
    error = function(e) {
      cat(paste0("Failed to fit model: ", e), "\n")
    })
  
  
  if(inherits(fit, "glmerMod")) {
    
    v_messages <- fit@optinfo$conv$lme4$messages
    if(!is.null(v_messages) && grepl("Model failed to converge", v_messages)) {
      
      control <- lme4::glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      
      fit <- tryCatch(
        do.call(getExportedValue("lme4", "glmer"), args = list(formula = formula, data = df,
                                                               family = binomial,
                                                               nAGQ = nAGQ,
                                                               control = control
        )),
        error = function(e) {
          cat(paste0("Failed to fit model: ", e), "\n")
        })
      
    }
    res <- as.data.table(summary(fit)$coefficients, keep.rownames = TRUE)
    
    v_message <- fit@optinfo$conv$lme4$messages
    v_message <- ifelse(is.null(v_message), NA, v_message)
    
    res <- res[, list(term = res$rn, estimate = res$Estimate, std.error = res$`Std. Error`,  statistic = res$`z value`, p.value = res$`Pr(>|z|)`,  n = fit@devcomp$dims[["n"]], message = fit@optinfo$conv$lme4$messages)]
    if(! is.null(keep)) res <- res[res$term %in% keep, ]
  } else {
    res <- list(estimate = NA)
  }
  return(res)
}




#' @description `fit_lmer`: linear mixed-effects model \code{\link[lme4]{lmer}}.
#' @param \dots Further arguments passed to regression model.
#' @rdname fit_lm
#' @export
fit_lmer <- function(data = NULL, formula = NULL, keep = NULL, ...) {
  
  check_pkg("lmerTest")
  
  v_var <- all.vars(formula)
  
  v_args <- list(...)
  
  if ("random" %in% names(v_args)) {
    random <- v_args[["random"]]
    formula <- as.formula(paste0(deparse(formula)," ",  gsub("~", " + ", deparse(random))))
    v_var <- c(v_var, all.vars(random))
  }
  
  v_args <- v_args[! names(v_args) %in% c("random")]
  
  if(length(v_args) != 0) {
    warning(paste0("Unknown args in lmer: ", paste0(names(v_args), collapse = ", "),"\n"))
  }
  
  # select variables.
  df <- data[, v_var, with = FALSE]
  
  check_pkg("lme4")
  
  fit <- tryCatch(
    do.call(getExportedValue("lmerTest", "lmer"), args = list(formula = formula, data = df)),
    error = function(e) {
      cat(paste0("Failed to fit model: ", e), "\n")
    })
  
  if(inherits(fit, "lmerMod")) {
    res <- as.data.table(summary(fit)$coefficients, keep.rownames = TRUE)
    v_message <- fit@optinfo$conv$lme4$messages
    v_message <- ifelse(is.null(v_message), NA, v_message)
    res <- res[, list(term = res$rn, estimate = res$Estimate, std.error = res$`Std. Error`,  statistic = res$`t value`, p.value = res$`Pr(>|t|)`,  n = fit@devcomp$dims[["n"]], message = v_message)]
    if(! is.null(keep)) res <- res[res$term %in% keep, ]
  } else {
    res <- list(estimate = NA)
  }
  return(res)
}




#' @rdname regression
#' @export
regression_each <- function(object, phenoData = NULL, model = NULL, formula = NULL, 
                            outcome = NULL, covars = NULL, factors = NULL, feature_name = NULL, 
                            time = NULL, verbose = TRUE, ncpus = 1, p.adjust.method = "bonferroni", ...) {


  if(is.null(model)) {
    stop("model is not specified.", call. = FALSE)
  }

  v_args <- list(...)

  if(!is.null(feature_name)) {
    x_list <- feature_name
  } else {
    x_list <- names(object@assayData)[2:NCOL(object@assayData)]
  }

  # phenoData, merge with sampleData
  if(is.null(phenoData)) {
    warnings(paste0("`phenoData` is NULL, `@sampleData` will be used. "))
    phenoData <- copy(object@sampleData)
  } else {
    phenoData <- as.data.table(phenoData)
    stopifnot(object@sampleID %in% names(phenoData))
    # phenoData <- merge(phenoData, object@sampleData, by = object@sampleID, suffixes = c("", "_"))
  }

  stopifnot(!any(x_list %in% names(phenoData)))

  # check column names

  if(is.null(formula)) stopifnot(outcome %in% names(phenoData))

  if(is.null(formula)) {
    cat("Build formula using covars and outcomes. \n")
    if(is.null(covars)) stop("Both formula and covars are NULL.")
    v_formula <- paste(outcome, "~", paste(covars, sep="", collapse = " + "))

  } else {
    if(!inherits(formula, "formula")) formula <- as.formula(formula)
    covars <- all.vars(formula)
    outcome <- covars[1]
    v_formula <- deparse(formula)
  }

  stopifnot(all(covars %in% names(phenoData)))

  if (!is.null(factors)) {
    stopifnot(all(factors %in% names(phenoData)))
    stopifnot(all(factors %in% covars))
    phenoData[, (factors) := lapply(.SD, as.factor) , .SDcols = factors]
  }

  df <- merge(phenoData, object@assayData, by = object@sampleID)

  if(model == "auto") {
    cat("model = `auto`, to infer `lm` or `logistic` \n")
    if(all(unique(df[, get(outcome)]) %in% c(0, 1, NA))) {
      model <- "logistic"
    } else {
      model <- "lm"
    }
  }
  
  fit_model <- paste0("fit_", model)
  fit_model <- get(fit_model)

  if(NROW(df) < 3) {
    stop(paste0("Only "), NROW(df), " rows in the merged file.", call. = FALSE)
  }

  if(model == "cox") {
    if(is.null(formula)) {
      if(is.null(time)) stop("`time` is NULL.", call. = FALSE)
      v_formula <- paste('survival::Surv(time = ', time, ', event = ', outcome ,') ~', paste(covars, sep= "", collapse = " + "))
    }
    outcome <- all.vars(as.formula(v_formula))[2]
  } else {
    outcome <- all.vars(as.formula(v_formula))[1]
  }

  if(verbose) {
    cat(paste0("\nRun `", model, "` model for ", length(x_list), " features: \n"))
    cat(v_formula, "+ `feature` \n")
  }
  
  if(ncpus > 1) {
    future::plan(future::multiprocess, workers = ncpus)
  }

  lapply_parallel <- if(future::nbrOfWorkers() > 1) {
    if(verbose) cat(paste0("Parallele runing with ncpus = ", future::nbrOfWorkers()), "\n")
    future.apply::future_lapply
  } else {
      pbapply::pblapply
  }
  
  v_ids <- 1L:length(x_list)
  p <- progressr::progressor(along = v_ids)

  fit_res <- lapply_parallel(
    X = v_ids,
    FUN = function(i) {
      p(sprintf("i=%g", i))
      v_i <- x_list[i]
      v_formula_i <- as.formula(paste0(v_formula, "+", v_i))
      tryCatch(
        do.call(fit_model, args = list(data = df, formula = v_formula_i, keep = v_i, ...)),
        error = function(e) {
          paste0("Failed to fit model: ", e)
          list()
        })
    }
  )
  if(all(is.na(fit_res))) return(NA)
  fit_res <- rbindlist(fit_res, fill=TRUE)
  fit_res$outcome <- outcome
  fit_res$model <- model
  fit_res$p.value.adj <- stats::p.adjust(fit_res$p.value, method = p.adjust.method)
  return(fit_res)
}





#' @rdname regression
#' @note \code{regression_each_as_outcome}: Run linear regression models where each feature is outcome.
#' @export
regression_each_as_outcome <- function(object, phenoData = NULL,
                            exposure = NULL, covars = NULL, factors = NULL, feature_name = NULL, 
                            verbose = TRUE, ncpus = 1, p.adjust.method = "bonferroni", ...) {
  
  
  model <- "lm"
  if(is.null(model)) {
    stop("model is not specified.", call. = FALSE)
  }
  
  v_args <- list(...)
  
  if(!is.null(feature_name)) {
    x_list <- feature_name
  } else {
    x_list <- names(object@assayData)[2:NCOL(object@assayData)]
  }
  
  # phenoData, merge with sampleData
  if(is.null(phenoData)) {
    warnings(paste0("`phenoData` is NULL, `@sampleData` will be used. "))
    phenoData <- copy(object@sampleData)
  } else {
    phenoData <- as.data.table(phenoData)
    stopifnot(object@sampleID %in% names(phenoData))
    # phenoData <- merge(phenoData, object@sampleData, by = object@sampleID, suffixes = c("", "_"))
  }
  
  stopifnot(!any(x_list %in% names(phenoData)))
  
  
    # cat("Build formula using covars and outcomes. \n")
    if(is.null(covars)) warning("Covars are NULL.")
    v_formula <- paste( "~", paste(c(exposure, covars), sep="", collapse = " + "))
  
  # print(covars)
  stopifnot(all(covars %in% names(phenoData)))
  
  if (!is.null(factors)) {
    stopifnot(all(factors %in% names(phenoData)))
    stopifnot(all(factors %in% covars))
    phenoData[, (factors) := lapply(.SD, as.factor) , .SDcols = factors]
  }
  
  df <- merge(phenoData, object@assayData, by = object@sampleID)
  
  
  fit_lm <- function(data = NULL, formula = NULL, keep = NULL) {
    v_var <- all.vars(formula)
    df <- data[, v_var, with = FALSE]
    fit <- tryCatch(
      do.call("lm", args = list(data = df, formula = formula)),
      error = function(e) {
        cat(paste0("Failed to fit model: ", e), "\n")
      })
    
    if(inherits(fit, "lm")) {
      res  <- as.data.table(broom::tidy(fit))
      if(! is.null(keep))  res <- res[res$term %in% keep, ]
      res$n <- length(fit$residuals)
      res$outcome <- v_var[1]
    } else {
      res <- list(estimate = NA)
    }
    return(res)
  }
  
  if(NROW(df) < 3) {
    stop(paste0("Only "), NROW(df), " rows in the merged file.", call. = FALSE)
  }
  
  
  if(verbose) {
    cat(paste0("\nRun `", model, "` model for ", length(x_list), " features: \n"))
    cat("`feature`", v_formula, " \n")
  }
  if(ncpus > 1) {
    future::plan(future::multiprocess, workers = ncpus)
  }
  
  lapply_parallel <- if(future::nbrOfWorkers() > 1) {
    if(verbose) cat(paste0("Parallele runing with ncpus = ", future::nbrOfWorkers()), "\n")
    future.apply::future_lapply
  } else {
    pbapply::pblapply
  }
  v_ids <- 1L:length(x_list)
  # progressr::handlers(global = TRUE)
  p <- progressr::progressor(along = v_ids)
  
  fit_res <- lapply_parallel(
    X = v_ids,
    FUN = function(i) {
      p(sprintf("i=%g", i))
      v_i <- x_list[i]
      v_formula_i <- as.formula(paste0(v_i, v_formula))
      tryCatch(
        do.call(fit_lm, args = list(data = df, formula = v_formula_i, keep = exposure, ...)),
        error = function(e) {
          paste0("Failed to fit model: ", e)
          list()
        })
    }
  )
  
  if(all(is.na(fit_res))) return(NA)
  fit_res <- rbindlist(fit_res, fill=TRUE)
  fit_res$exposure <- exposure
  fit_res$model <- model
  fit_res$p.value.adj <- stats::p.adjust(fit_res$p.value, method = p.adjust.method)
  return(fit_res)
}




