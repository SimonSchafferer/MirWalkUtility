% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/MirWalkRetrieval.R
\name{downloadMirWalkDataBySpecies}
\alias{downloadMirWalkDataBySpecies}
\title{Downloads all Data from Mirwalk}
\usage{
downloadMirWalkDataBySpecies(species = "mmu",
  basedir = "/home/simon/dbsOfflineUse/MirWalk2.0Resources")
}
\arguments{
\item{species}{- a character string (mmu, hsa, rno)}

\item{basedir}{- the root directory where the data should be stored}
}
\value{
TRUE if finished correctly
}
\description{
This function is a wrapper for all sub function to download specific mirwalk resources. I will organise the data in a directory structure in a given base directory
}

