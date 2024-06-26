% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeTreeSummarizedExperimentFromPhyloseq.R
\name{makeTreeSEFromPhyloseq}
\alias{makeTreeSEFromPhyloseq}
\alias{makeTreeSummarizedExperimentFromPhyloseq}
\alias{makeTreeSummarizedExperimentFromPhyloseq,ANY-method}
\title{Coerce a \code{phyloseq} object to a \code{TreeSummarizedExperiment}}
\usage{
makeTreeSEFromPhyloseq(x)

makeTreeSummarizedExperimentFromPhyloseq(x)

\S4method{makeTreeSummarizedExperimentFromPhyloseq}{ANY}(x)
}
\arguments{
\item{x}{a \code{phyloseq} object}
}
\value{
An object of class \code{TreeSummarizedExperiment}
}
\description{
\code{makeTreeSEFromPhyloseq} converts \code{phyloseq}
objects into \code{TreeSummarizedExperiment} objects.
}
\details{
All data stored in a \code{phyloseq} object is transferred.
}
\examples{
if (requireNamespace("phyloseq")) {
    data(GlobalPatterns, package="phyloseq")
    makeTreeSEFromPhyloseq(GlobalPatterns)
    data(enterotype, package="phyloseq")
    makeTreeSEFromPhyloseq(enterotype)
    data(esophagus, package="phyloseq")
    makeTreeSEFromPhyloseq(esophagus)
}
}
\seealso{
\code{\link[=makeTreeSEFromBiom]{makeTreeSEFromBiom}}
\code{\link[=makeTreeSEFromDADA2]{makeTreeSEFromDADA2}}
\code{\link[=importQIIME2]{importQIIME2}}
\code{\link[=importMothur]{importMothur}}
}
