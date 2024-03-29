% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regression.R
\name{regression}
\alias{regression}
\alias{regression_each}
\alias{regression_each_as_outcome}
\title{regression analysis}
\usage{
regression(
  object,
  phenoData = NULL,
  model = NULL,
  outcome = NULL,
  covars = NULL,
  factors = NULL,
  feature_name = NULL,
  time = NULL,
  verbose = TRUE,
  ncpus = 1,
  p.adjust.method = "bonferroni",
  ...
)

regression_each(
  object,
  phenoData = NULL,
  model = NULL,
  formula = NULL,
  outcome = NULL,
  covars = NULL,
  factors = NULL,
  feature_name = NULL,
  time = NULL,
  verbose = TRUE,
  ncpus = 1,
  p.adjust.method = "bonferroni",
  ...
)

regression_each_as_outcome(
  object,
  phenoData = NULL,
  exposure = NULL,
  covars = NULL,
  factors = NULL,
  feature_name = NULL,
  verbose = TRUE,
  ncpus = 1,
  p.adjust.method = "bonferroni",
  ...
)
}
\arguments{
\item{object}{A Metabolite object.}

\item{phenoData}{A data.table with outcome and covariates. If `phenoData` is NULL, `@sampleData` will be used.}

\item{model}{Specify a regression model. See \code{\link{fit_lm}} for more details. 'auto' can be used to infer 'lm' or 'logistic' (with only 0, 1, and NA).}

\item{outcome}{Column name of the outcome variable.}

\item{covars}{Column names of covariates.}

\item{factors}{Variables to be treated as factor.}

\item{feature_name}{A vector of selected metabolites to run. If both feature_name and random_select are NULL, will run regression for all features.}

\item{time}{Column name of survival time, used in cox regression, see \code{\link[survival]{coxph}} for more details.}

\item{verbose}{Print log information.}

\item{ncpus}{Number of CPUS for parallele job.}

\item{p.adjust.method}{Adjust for P value method, see \code{\link{p.adjust}}.}

\item{\dots}{Further arguments passed to regression model.}

\item{formula}{A character or formula object to fit model (only used in `regression_each`)}

\item{exposure}{exposure variables.}
}
\value{
term estimate std.error statistic p.value n outcome p.value.adj.
}
\description{
Run regression models with adjusting for covariates. `regression_each` is used for one outcome. In `regression`, several outcomes can be specified to run together.
}
\note{
\code{regression_each_as_outcome}: Run linear regression models where each feature is outcome.
}
\examples{
data(df_plasma)
fit_lm <- regression(object = df_plasma, phenoData = NULL, model = "lm", 
outcome = "BMI", covars = c("AGE", "GENDER", "ETHNICITY"), factors = "ETHNICITY")

}
