######Number of regulated genes#############
#####Based on all genes #############
RegGene1 <- function(maxExpr = 10, pVal= 0.001, FC=1, contrast = 1){
  pvalCol <- colnames(fData(gset)[grep("_PVal", colnames(fData(gset)))])[contrast]
  FCCol   <- colnames(fData(gset)[grep("_FC", colnames(fData(gset)))])[contrast]
  RegGenes <- fData(gset)[ (fData(gset)$maxExpr >maxExpr &
                              fData(gset)[pvalCol] <pVal&
                              abs(fData(gset)[FCCol]) >FC),]
  return (nrow(RegGenes))
}

RegGeneTable <- function(maxExpr = 0, pVal= 0.01, FC=.33) {   
  j=0
  ReguList <- list()
  for (i in colnames(cont.wt1)){
    j = j + 1
    ReguList[[i]] <- RegGene1(maxExpr = maxExpr, pVal= pVal, FC=FC, contrast = j)
  }
  Output <-data.frame(t(data.frame(ReguList)))
  colnames(Output) <- "Num"
  return(Output)
}
