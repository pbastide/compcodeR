% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateRmdCodeDiffExpPhylo.R
\name{DESeq2.length.createRmd}
\alias{DESeq2.length.createRmd}
\title{Generate a \code{.Rmd} file containing code to perform differential expression analysis with DESeq2 with custom model matrix}
\usage{
DESeq2.length.createRmd(
  data.path,
  result.path,
  codefile,
  fit.type,
  test,
  beta.prior = TRUE,
  independent.filtering = TRUE,
  cooks.cutoff = TRUE,
  impute.outliers = TRUE,
  extraDesignFactors = NULL,
  divByLengths = FALSE
)
}
\arguments{
\item{data.path}{The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.}

\item{result.path}{The path to the file where the result object will be saved.}

\item{codefile}{The path to the file where the code will be written.}

\item{fit.type}{The fitting method used to get the dispersion-mean relationship. Possible values are \code{"parametric"}, \code{"local"} and \code{"mean"}.}

\item{test}{The test to use. Possible values are \code{"Wald"} and \code{"LRT"}.}

\item{beta.prior}{Whether or not to put a zero-mean normal prior on the non-intercept coefficients. Default is \code{TRUE}.}

\item{independent.filtering}{Whether or not to perform independent filtering of the data. With independent filtering=TRUE, the adjusted p-values for genes not passing the filter threshold are set to NA.}

\item{cooks.cutoff}{The cutoff value for the Cook's distance to consider a value to be an outlier. Set to Inf or FALSE to disable outlier detection. For genes with detected outliers, the p-value and adjusted p-value will be set to NA.}

\item{impute.outliers}{Whether or not the outliers should be replaced by a trimmed mean and the analysis rerun.}

\item{extraDesignFactors}{A vector containing the extra factors to be passed to the design matrix of \code{DESeq2}. All the factors need to be a \code{sample.annotations} from the \code{\link{compData}} object. It should not contain the "condition" factor column, that will be added automatically.}

\item{divByLengths}{If TRUE, the counts are divided by the sequence lengths. If FALSE, the normalizing method explained in the details section is used. Default to FALSE.}
}
\value{
The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
}
\description{
A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the DESeq2 package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
}
\details{
For more information about the methods and the interpretation of the parameters, see the \code{DESeq2} package and the corresponding publications.


The length matrix are used as a normalization factor and applied to the \code{DESeq2}
model in the way explained in \code{\link[DESeq2]{normalizationFactors}}
(see details and examples of this function).
The provided matrix will be multiplied by the default normalization factor 
obtained through the \code{\link[DESeq2]{estimateSizeFactors}} function.
}
\examples{
try(
if (require(DESeq2)) {
tmpdir <- normalizePath(tempdir(), winslash = "/")
mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
                                    samples.per.cond = 5, n.diffexp = 100, 
                                    id.species = 1:10,
                                    lengths.relmeans = rpois(1000, 1000),
                                    lengths.dispersions = rgamma(1000, 1, 1),
                                    output.file = file.path(tmpdir, "mydata.rds"))
## Add annotations
sample.annotations(mydata.obj)$test_factor <- factor(rep(1:2, each = 5))
sample.annotations(mydata.obj)$test_reg <- 1:10
saveRDS(mydata.obj, file.path(tmpdir, "mydata.rds"))
## Diff Exp
runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "DESeq2", 
           Rmdfunction = "DESeq2.length.createRmd", 
           output.directory = tmpdir, fit.type = "parametric",
           test = "Wald",
           extraDesignFactors = c("test_factor", "test_reg"))
})
}
\references{
Anders S and Huber W (2010): Differential expression analysis for sequence count data. Genome Biology 11:R106
}
\author{
Charlotte Soneson, Paul Bastide, Mélina Gallopin
}
