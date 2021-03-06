% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/MirWalkRetrieval.R
\name{fetchLncInteraction_mirwalk}
\alias{fetchLncInteraction_mirwalk}
\title{Downloads the object from mirwalk}
\usage{
fetchLncInteraction_mirwalk(species = "hsa", identifier = "gene_symbol",
  baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/lncrna/")
}
\arguments{
\item{species}{("hsa","mmu","rno")}

\item{identifier}{- ("gene_symbol","ensembl", rno:"noncode")}

\item{baseurl}{- base url of the resource}
}
\value{
TRUE if finished correctly
}
\description{
This function downloads lncRNA interaction datasets from mirwalk and directly loads the object into workspace
}

