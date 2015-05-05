#' @title Downloads all Data from Mirwalk
#' @description This function is a wrapper for all sub function to download specific mirwalk resources. I will organise the data in a directory structure in a given base directory
#' @param species - a character string (mmu, hsa, rno)
#' @param basedir - the root directory where the data should be stored
#' @return TRUE if finished correctly
#' @export
downloadMirWalkDataBySpecies = function( species="mmu", basedir="/tmp/MirWalk2.0Resources"){
  
  if(species %in% c("mmu","hsa")){
    downloadPredictedTargets_mirwalk(species=species, basedir=basedir)
    downloadClipDataset_mirwalk(species=species, basedir=basedir)
    downloadLncInteraction_mirwalk(species=species, basedir=basedir)
    downloadValidatedTargets_mirwalk(species=species, basedir=basedir)
    downloadPathwayInteractions_mirwalk(species=species, basedir=basedir)    
  } else if(species %in% c("rno")){
    downloadPredictedTargets_mirwalk(species=species, basedir=basedir)
    downloadLncInteraction_mirwalk(species=species, basedir=basedir)
    downloadValidatedTargets_mirwalk(species=species, basedir=basedir)
    downloadPathwayInteractions_mirwalk(species=species, basedir=basedir)
  }
  return(TRUE)
}

#' @title Download target Predictions
#' @description This function bulk downloads all predicted miRNA targets of a given species
#' @param baseurl (default from website)
#' @param basedir - directory where the data should be stored
#' @return TRUE if finished correctly
#' @export
downloadPredictedTargets_mirwalk = function(species=species, 
                                            baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/r/", 
                                            basedir=basedir){
  dir.create(basedir)
  gene_regions = c("utr5","cds","utr3", "promoter")
  species_supported = c("hsa","mmu","rno")
  ids_supported = c("gene_symbol","entrez")
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  
  tmp=sapply(  gene_regions, function(gene_region){
    sapply( ids_supported, function(id_supported){
      resL = fetchPredictedTargets_mirwalk(species=species, geneRegion=gene_region, identifier=id_supported,baseurl=baseurl )
      setwd( basedir )
      currDir = file.path(basedir, "predictedTargets",species,gene_region, id_supported)
      dir.create( currDir, recursive=TRUE )
      assign( resL$filename, resL$object)
      save( list=resL$filename , file=file.path(currDir, paste0(resL$filename,".rda") ) )
      message( paste0("Saving ",file.path(currDir, paste0(resL$filename,".rda") ), " ... ")  )
      print(currDir)
      
    } )
  } )
  return(TRUE)
}

#' @title Loads a stored object into workspace
#' @description This function loads the prediceted target list into workspace
#' @param species ("hsa","mmu","rno")
#' @param geneRegion ("utr5","cds","utr3", "promoter")
#' @param identifier - (gene_symbol, entrez)
#' @param basedir - the base directory where the data is stored
#' @return TRUE if finished correctly
#' @export
getPredictedTarget_mirwalk = function( species, geneRegion, identifier, basedir){
  gene_regions = c("utr5","cds","utr3", "promoter")
  species_supported = c("hsa","mmu","rno")
  ids_supported = c("gene_symbol","entrez")
  if( !geneRegion %in% gene_regions){ stop( paste0("Allowed gene regions: ", paste0(gene_regions,collapse=", ") ) ) }
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed gene regions: ", paste0(ids_supported,collapse=", ") ) ) }
  fileToLoad = list.files( file.path(basedir, "predictedTargets",species,geneRegion, identifier), pattern="rda") 
  
  res = tryCatch({
    load( file.path(basedir, "predictedTargets",species,geneRegion, identifier, fileToLoad ) )
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
    message( "File is missing,please download files first!")     
  }, finally = {
  })
  
  res = load( file.path(basedir, "predictedTargets",species,geneRegion, identifier, fileToLoad ) )
  return( eval(parse(text=res)) )    
}


#' @title Downloads the object from mirwalk
#' @description This function downloads predicted miRNA targets from mirwalk and directly loads the object into workspace
#' @param species ("hsa","mmu","rno")
#' @param geneRegion ("utr5","cds","utr3", "promoter")
#' @param identifier - (gene_symbol, entrez)
#' @param baseurl - base url of the resource
#' @return TRUE if finished correctly
#' @export
fetchPredictedTargets_mirwalk = function( species="hsa", geneRegion="utr5", identifier="gene_symbol",  baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/r/" ){
  
  gene_regions = c("utr5","cds","utr3", "promoter")
  species_supported = c("hsa","mmu","rno")
  ids_supported = c("gene_symbol","entrez")
  
  if( !geneRegion %in% gene_regions){ stop( paste0("Allowed gene regions: ", paste0(gene_regions,collapse=", ") ) ) }
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed gene regions: ", paste0(ids_supported,collapse=", ") ) ) }
  
  type = switch(identifier, 
                gene_symbol = "g",
                entrez = "e",
  )
  gene_region = switch(geneRegion, 
                       utr5 = "5",
                       cds = "c",
                       utr3 = "3",
                       promoter = "p"
  )
  
  currWD = getwd()
  currTmpFile = tempdir()
  dir.create(currTmpFile)
  
  result = tryCatch({
    setwd(currTmpFile)
    
    filename = paste0( "s",species,gene_region,type,"r.zip" )
    download.file( paste0(baseurl, filename) , file.path(currTmpFile,filename) )
    
    message(paste0("Downloading file: ", paste0(baseurl, filename), " to ", currTmpFile ))
    system( paste0( "unzip -o ", file.path(currTmpFile,filename) ) )
    message(paste0( "Unzipping file ", file.path(currTmpFile,filename) ))
    filenameR = sub("\\.zip","", filename)
    
    allFiles = list.files( currTmpFile )
    filenameR = allFiles[grep( paste0(filenameR,"\\.rd.*$"), allFiles, ignore.case=TRUE )]
    
    tmp = load( file.path(currTmpFile,filenameR ))
    message(paste0("Loaded ", filenameR))
    
    setwd(currWD)
    return( list("filename"=sub("\\.rd.*$","",filenameR, ignore.case=TRUE), "object"=eval(parse(text=tmp)) ))
    
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
  }, finally = {
    setwd(currWD)
  })
  
  return( result )
  
}

#' @title Download CLIP datasets from mirwalk
#' @description This function bulk downloads all CLIP-Seq data for miRNAs of a given species
#' @param baseurl (default from website)
#' @param basedir - directory where the data should be stored
#' @return TRUE if finished correctly
#' @export
downloadClipDataset_mirwalk = function(species="hsa", 
                                       baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/clip/", 
                                       basedir=basedir){
  dir.create(basedir)
  gene_regions = c("utr5","cds","utr3")
  species_supported = c("hsa","mmu")
  if(species == "mmu"){
    method_supported = c("PARCLIP","HITCLIP","CLIPSEQ")  
  } else{
    method_supported = c("CLASH","PARCLIP","HITCLIP","iCLIP","CLIPSEQ")  
  }
  ids_supported = c("refseq")
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  
  tmp=sapply( gene_regions, function(geneRegion){
    sapply( method_supported, function(method) {
      sapply( ids_supported, function(id){
        resL = fetchClipDataset_mirwalk(species=species,geneRegion=geneRegion,method=method, identifier=id, baseurl=baseurl )
        setwd( basedir )
        currDir = file.path(basedir, "clipData",species,method,geneRegion,id)
        dir.create( currDir, recursive=TRUE )
        resL$filename = gsub("-","",resL$filename)
        assign( resL$filename, resL$object)
        save( list=resL$filename , file=file.path(currDir, paste0(resL$filename,".rda") ) )
        message( paste0("Saving ",file.path(currDir, paste0(resL$filename,".rda") ), " ... ")  )
        print(currDir)
      } )
    })    
  } )
  return(TRUE)
}

#' @title Downloads the object from mirwalk
#' @description This function downloads CLIP datasets from mirwalk and directly loads the object into workspace
#' @param species ("hsa","mmu","rno")
#' @param geneRegion ("utr5","cds","utr3")
#' @param identifier - (refseq)
#' @param method - (mmu: PARCLIP","HITCLIP","CLIPSEQ", hsa: "CLASH","PARCLIP","HITCLIP","iCLIP","CLIPSEQ")
#' @param baseurl - base url of the resource
#' @return TRUE if finished correctly
#' @export
fetchClipDataset_mirwalk = function( species="hsa", geneRegion="utr5", method="CLASH", identifier="refseq",  baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/clip/" ){
  
  gene_regions = c("utr5","cds","utr3")
  species_supported = c("hsa","mmu")
  
  if(species == "mmu"){
    method_supported = c("PARCLIP","HITCLIP","CLIPSEQ")  
  } else{
    method_supported = c("CLASH","PARCLIP","HITCLIP","iCLIP","CLIPSEQ")  
  }
  
  ids_supported = c("refseq")
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed Annotation: ", paste0(ids_supported,collapse=", ") ) ) }
  if( !method %in% method_supported){ stop( paste0("Allowed methods: ", paste0(method_supported,collapse=", ") ) ) }
  if( !geneRegion %in% gene_regions){ stop( paste0("Allowed gene region: ", paste0(gene_regions,collapse=", ") ) ) }
  
  gene_region = switch(geneRegion, 
                       utr5 = "5utr",
                       cds = "cds",
                       utr3 = "3utr"
  )
  type = switch(identifier, 
                refseq = "R"
  )
  
  currWD = getwd()
  currTmpFile = file.path( tempdir(), paste0(species,"clip",type ) )
  dir.create(currTmpFile)
  
  result = tryCatch({
    setwd(currTmpFile)
    
    filename = paste0(species,"-clip",type,".zip" )
    download.file( paste0(baseurl, filename) , file.path(currTmpFile,filename) )
    
    message(paste0("Downloading file: ", paste0(baseurl, filename), " to ", currTmpFile ))
    system( paste0( "unzip -o ", file.path(currTmpFile,filename) ) )
    message(paste0( "Unzipping file ", file.path(currTmpFile,filename) ))
    
    system( paste0( "rm ", file.path(currTmpFile,filename) ) )
    filenameOI = paste0(species,"-",gene_region,"-",method,".zip")
    system( paste0( "unzip -o ", file.path(currTmpFile,filenameOI) ) )
    filenameOI = sub("\\.zip","", filenameOI)
    allFiles = list.files( currTmpFile )
    filenameR = allFiles[grep( paste0(filenameOI,"\\.rd.*$"), allFiles, ignore.case=TRUE )]
    
    tmp = load( file.path(currTmpFile,filenameR ))
    message(paste0("Loaded ", filenameR))    
    
    setwd(currWD)
    return( list("filename"=sub("\\.rd.*$","",filenameR, ignore.case=TRUE), "object"=eval(parse(text=tmp)) ))
    
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
  }, finally = {
    setwd(currWD)
  })
  
  return(result)
}

#' @title Loads a stored object into workspace
#' @description This function loads the CLIP data into workspace
#' @param species ("hsa","mmu")
#' @param geneRegion ("utr5","cds","utr3")
#' @param identifier - (refseq)
#' @param method - (mmu: PARCLIP","HITCLIP","CLIPSEQ", hsa: "CLASH","PARCLIP","HITCLIP","iCLIP","CLIPSEQ")
#' @param basedir - the base directory where the data is stored
#' @return TRUE if finished correctly
#' @export
getClipDataset_mirwalk = function(  species="hsa", geneRegion="utr5", method="CLASH", identifier="refseq", basedir){
  
  gene_regions = c("utr5","cds","utr3")
  species_supported = c("hsa","mmu")
  if(species == "mmu"){
    method_supported = c("PARCLIP","HITCLIP","CLIPSEQ")  
  } else{
    method_supported = c("CLASH","PARCLIP","HITCLIP","iCLIP","CLIPSEQ")  
  }
  ids_supported = c("refseq")
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !geneRegion %in% gene_regions){ stop( paste0("Allowed gene regions: ", paste0(gene_regions,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed gene identifier: ", paste0(ids_supported,collapse=", ") ) ) }
  if( !method %in% method_supported){ stop( paste0("Allowed methods: ", paste0(method_supported,collapse=", ") ) ) }
  
  dirToLoad = file.path(basedir, "clipData",species,method,geneRegion,identifier)
  
  fileToLoad = list.files(dirToLoad , pattern="rda") 
  
  res = tryCatch({
    load( file.path(dirToLoad, fileToLoad ) )
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
    message( "File is missing,please download files first!")     
  }, finally = {
  })  
  return( eval(parse(text=res)) )    
}

#' @title Download lncRNA interaction datasets from mirwalk
#' @description This function bulk downloads all lncRNA interaction data for miRNAs of a given species
#' @param baseurl (default from website)
#' @param basedir - directory where the data should be stored
#' @return TRUE if finished correctly
#' @export
downloadLncInteraction_mirwalk = function(species="hsa", 
                                          baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/lncrna/", 
                                          basedir=basedir){
  dir.create(basedir)
  species_supported = c("hsa","mmu","rno")
  if(  species == "rno"){
    ids_supported = c("noncode")  
  }  else{
    ids_supported = c("gene_symbol","ensembl")
  }
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  
  tmp=sapply( ids_supported, function(id_supported){
    resL = fetchLncInteraction_mirwalk(species=species, identifier=id_supported,baseurl=baseurl )
    setwd( basedir )
    currDir = file.path(basedir, "lncRNAInteraction",species,id_supported)
    dir.create( currDir, recursive=TRUE )
    resL$filename = gsub("-","",resL$filename)
    assign( resL$filename, resL$object)
    save( list=resL$filename , file=file.path(currDir, paste0(resL$filename,".rda") ) )
    message( paste0("Saving ",file.path(currDir, paste0(resL$filename,".rda") ), " ... ")  )
    print(currDir)
  } )
  return(TRUE)
}

#' @title Downloads the object from mirwalk
#' @description This function downloads lncRNA interaction datasets from mirwalk and directly loads the object into workspace
#' @param species ("hsa","mmu","rno")
#' @param identifier - ("gene_symbol","ensembl", rno:"noncode")
#' @param baseurl - base url of the resource
#' @return TRUE if finished correctly
#' @export
fetchLncInteraction_mirwalk = function( species="hsa", identifier="gene_symbol",  baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/lncrna/" ){
  
  species_supported = c("hsa","mmu","rno")
  if(  species == "rno"){
    ids_supported = c("noncode")  
  }  else{
    ids_supported = c("gene_symbol","ensembl")
  }
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed identifier: ", paste0(ids_supported,collapse=", ") ) ) }
  
  type = switch(identifier, 
                gene_symbol = "lncRNAgr",
                ensembl = "lncRNAer",
                noncode = "LncRNA-noncodeR"
  )
  
  fileID = switch(identifier, 
                  gene_symbol = "gene",
                  ensembl = "ensembl",
                  noncode = "noncode"
  )
  
  currWD = getwd()
  currTmpFile = file.path(tempdir(), "lncRNA")
  dir.create(currTmpFile)
  
  result = tryCatch({
    setwd(currTmpFile)
    
    if( species == "rno" ){
      filename = paste0( species,type,".zip" )  
    } else{
      filename = paste0( "s",species,type,".zip" )
    }
    
    download.file( paste0(baseurl, filename) , file.path(currTmpFile,filename) )
    
    message(paste0("Downloading file: ", paste0(baseurl, filename), " to ", currTmpFile ))
    system( paste0( "unzip -o ", file.path(currTmpFile,filename) ) )
    message(paste0( "Unzipping file ", file.path(currTmpFile,filename) ))
    filenameR = sub("\\.zip","", filename)
    
    allFiles = list.files( currTmpFile )
    filenameR = allFiles[grep( paste0(".*",species,".*",fileID,".*\\.rd.*$"), allFiles, ignore.case=TRUE )]
    
    tmp = load( file.path(currTmpFile,filenameR ))
    message(paste0("Loaded ", filenameR))
    
    setwd(currWD)
    return( list("filename"=sub("\\.rd.*$","",filenameR, ignore.case=TRUE), "object"=eval(parse(text=tmp)) ))
    
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
  }, finally = {
    setwd(currWD)
  })
  
  return( result )
}

#' @title Loads a stored object into workspace
#' @description This function loads the CLIP data into workspace
#' @param species ("hsa","mmu","rno")
#' @param identifier - ("gene_symbol","ensembl", rno:"noncode")
#' @param basedir - the base directory where the data is stored
#' @return TRUE if finished correctly
#' @export
getLncInteraction_mirwalk = function(  species="hsa", identifier="gene_symbol", basedir){
  
  species_supported = c("hsa","mmu","rno")
  if(  species == "rno"){
    ids_supported = c("noncode")  
  }  else{
    ids_supported = c("gene_symbol","ensembl")
  }
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed gene regions: ", paste0(ids_supported,collapse=", ") ) ) }
  
  dirToLoad = file.path(basedir, "lncRNAInteraction",species,identifier)
  fileToLoad = list.files( dirToLoad, pattern="rda") 
  
  res = tryCatch({
    load( file.path(dirToLoad, fileToLoad ) )
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
    message( "File is missing,please download files first!")     
  }, finally = {
  })
  
  res = load( file.path(dirToLoad, fileToLoad ) )
  return( eval(parse(text=res)) )    
}


#' @title Download validated miRNA datasets
#' @description This function bulk downloads all validated target interactions for miRNAs of a given species
#' @param baseurl (default from website)
#' @param basedir - directory where the data should be stored
#' @return TRUE if finished correctly
#' @export
downloadValidatedTargets_mirwalk = function(species="hsa", 
                                            baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/vtm/", 
                                            basedir=basedir){
  dir.create(basedir)
  species_supported = c("hsa","mmu","rno")
  ids_supported = c("gene_symbol","entrez")
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  
  tmp=sapply( ids_supported, function(id_supported){
    resL = fetchValidatedTargets_mirwalk(species=species, identifier=id_supported,baseurl=baseurl )
    setwd( basedir )
    currDir = file.path(basedir, "validatedTargets",species,id_supported)
    dir.create( currDir, recursive=TRUE )
    resL$filename = gsub("-","",resL$filename)
    assign( resL$filename, resL$object)
    save( list=resL$filename , file=file.path(currDir, paste0(resL$filename,".rda") ) )
    message( paste0("Saving ",file.path(currDir, paste0(resL$filename,".rda") ), " ... ")  )
    print(currDir)
  } )
  return(TRUE)
}

#' @title Downloads the object from mirwalk
#' @description This function downloads validated targets datasets from mirwalk and directly loads the object into workspace
#' @param species ("hsa","mmu","rno")
#' @param identifier - ("gene_symbol","entrez")
#' @param baseurl - base url of the resource
#' @return TRUE if finished correctly
#' @export
fetchValidatedTargets_mirwalk = function( species="hsa", identifier="gene_symbol",  baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/vtm/" ){
  
  species_supported = c("hsa","mmu","rno")
  ids_supported = c("gene_symbol","entrez")
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed identifier: ", paste0(ids_supported,collapse=", ") ) ) }
  
  type = switch(identifier, 
                gene_symbol = "gene",
                entrez = "entrez"
  )
  
  currWD = getwd()
  currTmpFile = file.path(tempdir(), "validatedTargets")
  dir.create(currTmpFile)
  
  result = tryCatch({
    setwd(currTmpFile)
    
    filename = paste0( species,"-vtm-",type,".rdata.zip" )
    
    download.file( paste0(baseurl, filename) , file.path(currTmpFile,filename) )
    
    message(paste0("Downloading file: ", paste0(baseurl, filename), " to ", currTmpFile ))
    system( paste0( "unzip -o ", file.path(currTmpFile,filename) ) )
    message(paste0( "Unzipping file ", file.path(currTmpFile,filename) ))
    filenameR = sub("\\.zip","", filename)
    
    allFiles = list.files( currTmpFile )
    filenameR = allFiles[grep( paste0(filenameR,".*\\.rd.*$"), allFiles, ignore.case=TRUE )]
    
    tmp = load( file.path(currTmpFile,filenameR ))
    message(paste0("Loaded ", filenameR))
    
    setwd(currWD)
    return( list("filename"=sub("\\.rd.*$","",filenameR, ignore.case=TRUE), "object"=eval(parse(text=tmp)) ))
    
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
  }, finally = {
    setwd(currWD)
  })
  
  return( result )
}

#' @title Loads a stored object into workspace
#' @description This function loads the Validated targets into workspace
#' @param species ("hsa","mmu","rno")
#' @param identifier - ("gene_symbol","entrez")
#' @param basedir - the base directory where the data is stored
#' @return TRUE if finished correctly
#' @export
getValidatedTargets_mirwalk = function(  species="hsa", identifier="gene_symbol", basedir){
  
  species_supported = c("hsa","mmu","rno")
  ids_supported = c("gene_symbol","entrez")
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed gene regions: ", paste0(ids_supported,collapse=", ") ) ) }
  
  dirToLoad = file.path(basedir, "validatedTargets",species,identifier)
  fileToLoad = list.files( dirToLoad, pattern="rda") 
  
  res = tryCatch({
    load( file.path(dirToLoad, fileToLoad ) )
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
    message( "File is missing,please download files first!")     
  }, finally = {
  })
  
  res = load( file.path(dirToLoad, fileToLoad ) )
  return( eval(parse(text=res)) )    
}

#' @title Download pathway - miRNA - datasets
#' @description This function bulk downloads all pathway datasets for miRNAs of a given species
#' @param baseurl (default from website)
#' @param basedir - directory where the data should be stored
#' @return TRUE if finished correctly
#' @export
downloadPathwayInteractions_mirwalk = function(species="hsa", 
                                               baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/vtm/", 
                                               basedir=basedir){
  dir.create(basedir)
  species_supported = c("hsa","mmu","rno")
  ids_supported = c("kegg","wiki","panther","GOBP","GOMF","GOCC","gClass","pClass")
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  
  tmp=sapply( ids_supported, function(id_supported){
    resL = fetchPathwayInteractions_mirwalk(species=species, identifier=id_supported,baseurl=baseurl )
    setwd( basedir )
    currDir = file.path(basedir, "pathways",species,id_supported)
    dir.create( currDir, recursive=TRUE )
    resL$filename = gsub("-","",resL$filename)
    assign( resL$filename, as.data.frame( resL$object ) )
    save( list=resL$filename , file=file.path(currDir, paste0(resL$filename,".rda") ) )
    message( paste0("Saving ",file.path(currDir, paste0(resL$filename,".rda") ), " ... ")  )
    print(currDir)
  } )
  return(TRUE)
}

#' @title Downloads the object from mirwalk
#' @description This function downloads validated targets datasets from mirwalk and directly loads the object into workspace
#' @param species ("hsa","mmu","rno")
#' @param identifier - ("kegg","wiki","panther","GOBP","GOMF","GOCC","gClass","pClass")
#' @param basedir - the base directory where the data is stored
#' @return TRUE if finished correctly
#' @export
fetchPathwayInteractions_mirwalk = function( species="hsa", identifier="kegg",  baseurl = "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/vtm/" ){
  
  species_supported = c("hsa","mmu","rno")
  ids_supported = c("kegg","wiki","panther","GOBP","GOMF","GOCC","gClass","pClass")
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed identifier: ", paste0(ids_supported,collapse=", ") ) ) }
  
  type = switch(identifier, 
                kegg = "kegg",
                wiki = "Wiki",
                panther = "PantherPath",
                GOBP = "GOBP",
                GOMF = "GOMF",
                GOCC = "GOCC",
                gClass = "geneclass",
                pClass = "proteinclass"
  )
  
  currWD = getwd()
  currTmpFile = file.path(tempdir(), "pathways")
  dir.create(currTmpFile)
  
  result = tryCatch({
    setwd(currTmpFile)
    
    filename = paste0( species,"-",type,"-vtm",".zip" )
    
    download.file( paste0(baseurl, filename) , file.path(currTmpFile,filename) )
    
    message(paste0("Downloading file: ", paste0(baseurl, filename), " to ", currTmpFile ))
    system( paste0( "unzip -o ", file.path(currTmpFile,filename) ) )
    message(paste0( "Unzipping file ", file.path(currTmpFile,filename) ))
    filenameR = sub("\\.zip","", filename)
    
    allFiles = list.files( currTmpFile )
    filenameR = allFiles[grep( paste0(filenameR,".*\\.txt.*$"), allFiles, ignore.case=TRUE )]
    
    tmp = read.csv( file.path(currTmpFile,filenameR ), sep="\t", header=TRUE)
    
    colnames(tmp)[grep("miRNAname",colnames(tmp),ignore.case=TRUE)] = "miRNAname"
    colnames(tmp)[grep("Gene",colnames(tmp),ignore.case=TRUE)] = "Gene"
    colnames(tmp)[grep("EntrezID",colnames(tmp),ignore.case=TRUE)] = "EntrezID"
    
    message(paste0("Loaded ", filenameR))
    
    setwd(currWD)
    return( list("filename"=sub("\\.rd.*$","",filenameR, ignore.case=TRUE), "object"=tmp ))
    
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
  }, finally = {
    setwd(currWD)
  })
  
  return( result )
}

#' @title Loads a stored object into workspace
#' @description This function loads the pathway interactions into workspace
#' @param species ("hsa","mmu","rno")
#' @param identifier - ("kegg","wiki","panther","GOBP","GOMF","GOCC","gClass","pClass")
#' @param basedir - the base directory where the data is stored
#' @return TRUE if finished correctly
#' @export
getPathwayInteractions_mirwalk = function(  species="hsa", identifier="kegg", basedir){
  
  species_supported = c("hsa","mmu","rno")
  ids_supported = c("kegg","wiki","panther","GOBP","GOMF","GOCC","gClass","pClass")
  
  if( !species %in% species_supported){ stop( paste0("Allowed species: ", paste0(species_supported,collapse=", ") ) ) }
  if( !identifier %in% ids_supported){ stop( paste0("Allowed gene regions: ", paste0(ids_supported,collapse=", ") ) ) }
  
  dirToLoad = file.path(basedir, "pathways",species,identifier)
  fileToLoad = list.files( dirToLoad, pattern="rda") 
  
  res = tryCatch({
    load( file.path(dirToLoad, fileToLoad ) )
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
    message( "File is missing,please download files first!")     
  }, finally = {
  })
  
  res = load( file.path(dirToLoad, fileToLoad ) )
  return( eval(parse(text=res)) )    
}
