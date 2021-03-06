% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/MirWalkRetrieval.R
\name{getPredictedTarget_mirwalk}
\alias{getPredictedTarget_mirwalk}
\title{Loads a stored object into workspace}
\usage{
getPredictedTarget_mirwalk(species, geneRegion, identifier, basedir)
}
\arguments{
\item{species}{("hsa","mmu","rno")}

\item{geneRegion}{("utr5","cds","utr3", "promoter")}

\item{identifier}{- (gene_symbol, entrez)}

\item{basedir}{- the base directory where the data is stored}
}
\value{
TRUE if finished correctly
}
\description{
This function loads the prediceted target list into workspace
}

