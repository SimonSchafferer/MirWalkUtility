# MirWalkUtility
Retrieval of data from the mirwalk2 website within R (http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/holistic.html)
The package uses the system command 'unzip'

# Usage

Install and load the package
```{r}
library(devtools)
install_github("SimonSchafferer/MirWalkUtility")
library(MirWalkUtility)
```
### Fetch individual files from MirWalk and load them into the workspace

```{r}
#Setting example values for base directory and species
basedir="/tmp/MirWalk2.0Resources"
species = "mmu" #other possibilities: hsa, rno

test = fetchPredictedTargets_mirwalk( species=species, geneRegion="cds", identifier="gene_symbol")
test = fetchClipDataset_mirwalk( species=species, geneRegion="cds", method="PARCLIP", identifier="refseq") 
test = fetchLncInteraction_mirwalk( species=species, identifier="gene_symbol")
test = fetchValidatedTargets_mirwalk( species=species, identifier="gene_symbol")           
test = fetchPathwayInteractions_mirwalk( species=species, identifier="kegg")
```                                                        

### Bulk download all data for a given species
```{r}
downloadMirWalkDataBySpecies(species=species, basedir=basedir)
```

### After succesfull download, get the files into the workspace by the following functions: 
```{r}
#individual files may be loaded into the workspace e.g.:
test = getPredictedTarget_mirwalk( species=species, "utr3", "gene_symbol", basedir)
test = getClipDataset_mirwalk(  species=species, geneRegion="utr5", method="PARCLIP", identifier="refseq", basedir)
test = getLncInteraction_mirwalk(  species=species, identifier="gene_symbol", basedir)
test = getValidatedTargets_mirwalk( species=species, identifier="gene_symbol", basedir)
test = getPathwayInteractions_mirwalk( species=species, identifier="kegg", basedir)
```
