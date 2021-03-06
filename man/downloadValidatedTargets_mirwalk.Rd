% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/MirWalkRetrieval.R
\name{downloadValidatedTargets_mirwalk}
\alias{downloadValidatedTargets_mirwalk}
\title{Download validated miRNA datasets}
\usage{
downloadValidatedTargets_mirwalk(species = "hsa",
  baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/vtm/",
  basedir = basedir)
}
\arguments{
\item{baseurl}{(default from website)}

\item{basedir}{- directory where the data should be stored}
}
\value{
TRUE if finished correctly
}
\description{
This function bulk downloads all validated target interactions for miRNAs of a given species
}

