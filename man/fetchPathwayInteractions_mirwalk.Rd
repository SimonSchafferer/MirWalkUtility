% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/MirWalkRetrieval.R
\name{fetchPathwayInteractions_mirwalk}
\alias{fetchPathwayInteractions_mirwalk}
\title{Downloads the object from mirwalk}
\usage{
fetchPathwayInteractions_mirwalk(species = "hsa", identifier = "kegg",
  baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/vtm/")
}
\arguments{
\item{species}{("hsa","mmu","rno")}

\item{identifier}{- ("kegg","wiki","panther","GOBP","GOMF","GOCC","gClass","pClass")}

\item{basedir}{- the base directory where the data is stored}
}
\value{
TRUE if finished correctly
}
\description{
This function downloads validated targets datasets from mirwalk and directly loads the object into workspace
}

