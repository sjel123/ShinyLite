####QC Functions####

NormalizeData <- function(data, design=design){
  require(edgeR)
  dge <- DGEList(counts=data)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=TRUE)
  return(v)
}

######Scatterplot Function
ScatterFunction <- function(data=DFS, xx= "A1", yy= "A2"){
  require(ggplot2)
  x=data[,xx]
  y=data[,yy]
  ggplot(as.data.frame(data), aes(x=x,y=y))+geom_point()
}


######Scatterplot Function by Group
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*0.5, col=c("gray60", "black")[(abs(r)>0.65)+1])
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2],0,1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


ScatterFunctionMatrix <- function(data=DFS, group="EMPTY_0"){
  SFMFunction <- function(group=group1){
    
    index=which(pData(gset)$Group==group)
    data1 <- data[,index]
    if(length(index)>10){index1 <- sample(index, 10, replace=FALSE)
      data1=data[,index1]
    }
    if(length(index)==1){
    data1=data[,c(index, index)]}
    
    if(nrow(data1)>1000){
      data1 <- data1[sample(1:nrow(data1), 1000),]
    }
    TITLE <- ifelse(length(index) > 10, paste0(group, " Subsample of 10"),  paste0(group, " ALL"))
    print(TITLE)
    pairs(data1, upper.panel = panel.cor, diag.panel = panel.hist, main=TITLE)
  }
  for (i in levels(factor(pData(gset)$Group))){
    SFMFunction(i)
  }
  par(mfrow=c(1,1))
}


############# Histogram function
HistFunc <- function (data=DFS){
  DFS1<- stack(as.data.frame(data))
  color = data.frame(ind=colnames(data), Factor)
  DFS1 <- merge(DFS1, color, by = "ind")
  p <- ggplot(DFS1, aes(x=values, color=Factor))+geom_density(aes(group=ind))
  p
}

#Correlation HeatMap
CorFunction1 <- function(data=gset){
  data1 <- data
  if (class(data)=="ExpressionSet") data1 <- exprs(data)
  #data1 <- data1[,order(colnames(data1))]
  CorData <- melt(cor(data1))
  p <- ggplot(CorData, aes(x=X1, y=X2, fill=value))+ geom_tile()+
    scale_fill_gradient(low="white", high="red")+ theme(axis.text.x=element_text(angle=90,hjust=1))+
    scale_x_discrete(limits=(unique(CorData$X1)))+
    scale_y_discrete(limits=(unique(CorData$X2)))
  print(p)
  #p <- ggplotly(p)
  p
}

CorFunction <- function(data=gset, orderby = "Hour", Labels="Donor"){
  #data1 <- data#[,order(colnames(data))]
  require(ggplot2)
  require(reshape2)
  Anno <- pData(data)[c(orderby, Labels)]
  Anno <- Anno[order(Anno[,1]),]
  #which(names(pData(data))==orderby)
  DFS_O <- exprs(data)#[,order(pData(data)[,orderby])]
    
  #colnames(DFS_O)<- pData(data)[Labels][,1]
  dim(DFS_O)
  DFS_OMin <- DFS_O [apply(DFS_O , 1, min)>1,]
    #tDFS_0<- c(data.frame(t(DFS_O)), pData(gset))
    #index <- which (colnames(pData(data))==orderby)
   # colnames(DFS_O)<-make.names(pData(data)[orderby][,1], unique = TRUE)
  CorData <- melt(cor(DFS_O))
  names(CorData)<- c("X1", "X2","value")
    CorData <- merge(CorData, Anno, by.x="X1", by.y=0)
    CorData <- merge(CorData, Anno, by.x="X2", by.y=0)
      CorData$X2 <- factor(CorData$X2, levels=rownames(Anno))
      CorData$X1 <- factor(CorData$X1, levels=rownames(Anno))
        Order_y <- Anno[match(levels(CorData$X1), row.names(Anno)), ]
        Order_x <- Anno[match(levels(CorData$X2), row.names(Anno)), ]
  

         # names(CorData) <- c("X1", "X2", "value")
   # CorData$X11 <-  gsub("\\.[1-9]", "",CorData$X1)
   # CorData$X22 <-  gsub("\\.[1-9]", "",CorData$X2)
   # CorData$X11 <- factor(CorData$X11, level= unique(as.character(pData(data)[Order,orderby])))
   # CorData$X22 <- factor(CorData$X22, level= unique(as.character(pData(data)[Order,orderby])))
   # ylabels <- Order_y[,2]
  p <- ggplot(CorData, aes(x=X1, y=X2, fill=value))+ geom_tile()+
    scale_fill_gradient(low="white", high="red")+ theme(axis.text.x=element_text(angle=90,hjust=1))
  p <- p + scale_y_discrete(labels=Anno[,2])
  p <- p + scale_x_discrete(labels=Anno[,2])
  
  #pData(data)[Order, Labels])
  #p <- p + scale_x_discrete(labels=pData(data)[Order, Labels])
  #p <- p + scale_x_discrete(limits=CorData$X1[order(CorData$X11)]) + 
  #  p <- p +  scale_y_discrete(limits=(x$V2)[order(x$V3)])
  print(p)
  #p <- ggplotly(p)
  #p
}

###Density Plots
#Histograms
DensityPlotFunction <- function(data=gsetHb){
  if (class(data)=="ExpressionSet") data <- exprs(data)
  DFS1<- stack(as.data.frame(data))
  SampleAnnotation <- pData(gset)
  if (class(data)=="ExpressionSet") SampleAnnotation <- pData(data)
  color = SampleAnnotation$Cohort#c(rep('green',5),rep('red',4), rep('blue',4))
  fill = c(rep(1,3),rep(0,2), rep(1,2), rep(0,2), rep(1,2), rep(0,2))
  p <- ggplot(DFS1, aes(x=values))+geom_density(aes(group=ind))
  print(p)
  p
}

###Limma Function
LimmaFunction <- function(){
  design <- model.matrix(~0+factor(pData(gset)$Group))
  colnames(design) <- gsub("factor\\(pData\\(gset\\)\\$Group\\)", '', colnames(design))
  colnames(design)<- make.names(colnames(design))
  #colnames(design) <- gsub("factor\\(pData\\(gset\\)\\$Donor\\)", '', colnames(design))
  
  fit1 <- lmFit(gset, design)
  cont.wt1 <- makeContrasts(contrasts=CONtrasts,
                            levels=design )
  fit1.wt <- (contrasts.fit(fit1, cont.wt1))  
  efit1 <- eBayes(fit1.wt)
  return(list(efit1, cont.wt1))
} 

#KeansClustering
## Prepare Data
AnovaAnalysis <- function(data = gsetHb){
  DFS <- exprs(gsetHb)
  poorQC <- !colnames(DFS)%in%poorqc
  data1 <- DFS[,poorQC]
  AovFactor <- Factor[poorQC]
  WTAnal.1 <- apply(data1, 1, function(x)summary(aov(x~ AovFactor))[[1]][1,5])
  WTAnal.THSD <- apply(data1, 1, function(x){
    a1 <- aov(x~ as.factor(AovFactor))
    posthoc <- TukeyHSD(a1)[[1]][,4]
  })
  
  WTAnal.2 <- apply(data1, 1, function(x){
    SampleAnnotation <- pdata(data)
    Mean <- ddply(data.frame(value=x, Treatment=SampleAnnotation$Cohort[poorQC]),
                  c("Treatment"),summarise,
                  mean =mean(value))
    max(Mean$mean)/min(Mean$mean)
  })
  WTAnal<- data.frame(Symbol=fData(data)$gene_name, P.Value=WTAnal.1, FC=WTAnal.2, t(WTAnal.THSD))
  #Aov.Results1 <- data.frame(as.character(DF[,1]),as.numeric(WTAnal.1))
  #colnames(Aov.Results1) <- c("Symbol", "p.value")
  #head(Aov.Results1[order(Aov.Results1$p.value),])
  return(WTAnal)
}

tukeysHSDLoop <- function(){
  poorQC <- !colnames(DFS)%in%poorqc
  data1 <- DFS[,poorQC]
  AovFactor <- Factor[poorQC]
  WTAnal.THSD <- apply(data1[1:2,], 1, function(x){
    a1 <- aov(x~ as.factor(AovFactor))
    posthoc <- 
      TukeyHSD(a1)[[1]][,4]
  })
}

VarFilter <- function(data, quant=0.95){
  data <- data[(rowSums(data)>10),]
  vardf <- apply(data, 1, var)
  vardfFilt <- vardf>quantile(vardf, quant)
  print(table(vardfFilt))
  data[vardfFilt,]
}

######Number of regulated genes#############
#####Based on all genes #############
RegGene <- function(maxExpr = 10, pVal= 0.001, FC=1, contrast = 1){
  pvalCol <- colnames(fData(gset)[grep("_PVal", colnames(fData(gset)))])[contrast]
  FCCol   <- colnames(fData(gset)[grep("_FC", colnames(fData(gset)))])[contrast]
  #pvalCol <- paste0("Hour_",Time, ".Hour_0_pVal")
  #FCCol <- paste0("Hour_",Time, ".Hour_0.1_FC")
  RegGenes <- fData(gset)[(fData(gset)$gene_type =="protein_coding" & 
                             fData(gset)$maxExpr >maxExpr &
                             fData(gset)[pvalCol] <pVal&
                             abs(fData(gset)[FCCol]) >FC),]
  return (droplevels(RegGenes$gene_name))
}

#####Based on all genes #############
RegGene1 <- function(maxExpr = 10, pVal= 0.001, FC=1, contrast = 1){
  pvalCol <- colnames(fData(gset)[grep("_PVal", colnames(fData(gset)))])[contrast]
  FCCol   <- colnames(fData(gset)[grep("_FC", colnames(fData(gset)))])[contrast]
  #pvalCol <- paste0("Hour_",Time, ".Hour_0_pVal")
  #FCCol <- paste0("Hour_",Time, ".Hour_0.1_FC")
  RegGenes <- fData(gset)[ (fData(gset)$gene_type =="protein_coding"&
                              fData(gset)$maxExpr >maxExpr &
                              fData(gset)[pvalCol] <pVal&
                              abs(fData(gset)[FCCol]) >FC),]
  return (RegGenes$gene_name)
}


###Combine limma output pvalue into single table  
CombineLimmaOutput <- function(efit=efit1, CONTRASTs = cont.wt1) {
  #fit <- contrasts.fit(FIT, CONTRASTs)  
  #efit <- eBayes(fit)
  
  pValueDF <- list()
  for (i in colnames(CONTRASTs)){
    pValueDF[[i]] <- topTable(efit,coef =i,  number="inf", sort.by = "none")[,"P.Value"]
  }
  names(pValueDF) <- paste0(names(pValueDF), "_PVal")
  pValueDF[["Limmaanova"]] <- topTable(efit, number="inf", sort.by = "none")[,"P.Value"] 
  
  ###Combine limma output Fold Change into single table
  FCDF <- list()
  for (i in colnames(CONTRASTs)){
    FCDF[[i]] <- topTable(efit,coef =i,  number="inf", sort.by = "none")[,"logFC"]
  }
  names(FCDF) <- paste0(names(FCDF), "_FC")
  #### Add output to dataframe
  DF <- data.frame(as.data.frame(pValueDF), as.data.frame(FCDF))
  return(DF)
}

