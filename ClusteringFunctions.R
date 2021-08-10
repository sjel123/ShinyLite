##Clustering Functions
#PCA PLOT
PCAFunction <- function(data=gset, labelCol="Cohort", shape1 ="Patient", color1="Cohort", 
                        Continous=FALSE){
  require(ggrepel)
  data1 <- as.data.frame(exprs(data))
  data1 <- data1[apply(data1,1,var)!=0,]
  data.PC = prcomp(t(data1), scale.=TRUE)
  PC1_percent <- summary(data.PC)$importance[2,1]
  PC1_percent <- signif(PC1_percent*100,2)
  PC2_percent <- summary(data.PC)$importance[2,2]
  PC2_percent <- signif(PC2_percent*100,2)
  SampleAnnotation <- pData(data)
  PCAData <- data.frame(data.PC$x[,1:2], SampleAnnotation)#, by.x=0, by.y="sample_id")#cbind(data.PC$x[,1:2], df1Annotation)
  
  if(labelCol == "none") {label1=""}
  if(labelCol != "none") {label1 <- PCAData[,labelCol]}
  PCAData$shape2 <- as.factor(PCAData[,shape1])
  if(Continous==FALSE)PCAData$color2 <- as.factor(PCAData[,color1])
  if(Continous!=FALSE)PCAData$color2 <- (PCAData[,color1])
  #SampleAnnotation$labelCol <- (SampleAnnotation[labelCol])
  p <- ggplot(data=PCAData, aes(x=PC1, y=PC2, 
                                shape = shape2,
                                color = color2,
                                label = label1))+
    geom_point(size = 4)+ geom_text_repel(hjust=1.4, nudge_x =0.05)
  p <- p + scale_shape_discrete(name=shape1)#, labels=PCAData$shape2)+
  if(Continous!=FALSE) p <- p + scale_color_continuous(name=color1)#, labels=PCAData$color2)
  if(Continous==FALSE) p <- p + scale_color_discrete(name=color1)#, labels=PCAData$color2)
  p <- p + theme(legend.text = element_text(size=12))
  p <- p + xlab(paste0("PC1 ", PC1_percent, "%")) + ylab(paste0("PC2 ", PC2_percent,"%"))
  return(p)
}
#PCAFunction(gsetMax,labelCol="Status",  shape1 ="Location_sampled_description", color1 = "Location.1")

PCAEsetFunction <- function(ESET=fgset, labelCol="Cohort"){
  data <- exprs(ESET)
  data <- data[which(apply(data,1, function(v) !all(is.na(v)))),]
  data.PC = prcomp(t(data), scale.=TRUE)
  #SampleAnnotation$labelCol <- (SampleAnnotation[labelCol])
  PCAData <- data.frame(data.PC$x[,1:2], pData(ESET))#, by.x=0, by.y="sample_id")#cbind(data.PC$x[,1:2], df1Annotation)
  p <- ggplot(data=PCAData, aes(x=PC1, y=PC2, 
                                shape = factor(RNA.Batch),
                                size = 2,#df1Annotation$cell,
                                color= Cohort,
                                label = pData(ESET)[labelCol]))+
    geom_point()+ geom_text(hjust=1.4, nudge_x =0.05)
  return(p)
  
}

###KMEANSClustering
#Kmeans clustering
KmeansClusterFun <- function(mydata1=gset, n=8, shade=5, ESET =NULL, orderBy = "Group", PCutoff=0.01){
  mydata3 <- exprs(ESET)[fData(ESET)$Limmaanova < PCutoff,]
  row.names(mydata3) <- make.names(fData(ESET)[fData(ESET)$Limmaanova < PCutoff,"gene_name"], unique = TRUE)
  mydata2 <- as.data.frame(t(scale(t(mydata3), center=TRUE, scale=TRUE)))
  #Find duplicated gene_names
  #dupnames <- as.numeric(duplicated(row.names(mydata2)))
  #dupnames <- gsub(1, "_1", dupnames)
  #dupnames <- gsub (0, "", dupnames)
  rownames(mydata2)  <- make.names(rownames(mydata2))#gsub (0, "", dupnames)
  colnames(mydata2)  <- make.names(pData(ESET)[,"Group"])
  #row.names(mydata2) <-paste0(rownames(mydata2),dupnames)
  
  #Order Data
  #colnames(mydata2) <- colnames(mydata1)
  mydata2 <- (mydata2[,order(pData(ESET)[orderBy])])
  
  wss <- (nrow(mydata2)-1)*sum(apply(mydata2,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata2,
                                       centers=i)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares") 
  
  # K-Means Cluster Analysis
  set.seed(122345)
  fit <- kmeans(mydata2, n) # 10 cluster solution
  
  ClusterLineGraphs (fit) 
  mydata2 <- data.frame(mydata2, fit$cluster)
  
  
  #   if (!is.null(ESET)) {
  #     if (dim(pData(ESET))[1] == length(colnames(mydata2))-1){
  #       colnames(mydata2) <- c(pData(ESET)$Cohortgroup, "fit.cluster")
  #     }
  #   }
  o <- order(mydata2[,"fit.cluster"])
  mydata2 <- mydata2[o,]
  library(gplots)
  library(RColorBrewer)
  #cols <- palette(brewer.pal(14, "Dark2"))
  hmcol <- colorRampPalette(brewer.pal(11, "RdBu"))(shade)
  
  nco <- ncol(mydata2)
  output <- heatmap(as.matrix(mydata2[,1:ncol(mydata2)-1]), col=rev(hmcol), Rowv=NA, Colv=NA)
  return(output)[[5]]
}  

ClusterLineGraphs <- function(fit){
  ClusterData <- data.frame(fit$center)
  #Order by Cluster
  #     o <- order(ClusterData[,"fit.cluster"])
  #     ClusterData <- ClusterData[o,-which(colnames(ClusterData) %in% "fit.cluster")]
  ClusterData$ClusterNumber <- row.names(ClusterData)
  ClusterData.melt <- melt(ClusterData, id = "ClusterNumber")
  ClusterData.melt<- ClusterData.melt[order(ClusterData.melt$ClusterNumber,decreasing = TRUE),]
  p <- ggplot(ClusterData.melt, aes(x=variable, y= value, group=ClusterNumber))+geom_line()+
    facet_wrap(~ClusterNumber, ncol=4)
  print(p)
  #order by sample
}

KmeansClusterFun2 <- function(mydata1=gsetMax, n=8, shade=5, ESET =NULL, orderBy = "Hour"){
  mydata3 <- exprs(mydata1)[fData(mydata1)$pVal <0.0001,]
  row.names(mydata3) <- fData(mydata1)[fData(mydata1)$pVal <0.0001,"gene_name"]
  mydata2 <- as.data.frame(t(scale(t(mydata3), center=TRUE, scale=TRUE)))
  #Find duplicated gene_names
  dupnames <- as.numeric(duplicated(row.names(mydata2)))
  dupnames <- gsub(1, "_1", dupnames)
  dupnames <- gsub (0, "", dupnames)
  row.names(mydata2) <-paste0(rownames(mydata2),dupnames)
  
  #Order Data
  #colnames(mydata2) <- colnames(mydata1)
  mydata2 <- (mydata2[,order(pData(ESET)[orderBy])])
  
  wss <- (nrow(mydata2)-1)*sum(apply(mydata2,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata2,
                                       centers=i)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares") 
  
  # K-Means Cluster Analysis
  set.seed(122345)
  fit <- kmeans(mydata2, n) # 10 cluster solution
  
  ClusterLineGraphs (fit) 
  mydata2 <- data.frame(mydata2, fit$cluster)
  
  
  #   if (!is.null(ESET)) {
  #     if (dim(pData(ESET))[1] == length(colnames(mydata2))-1){
  #       colnames(mydata2) <- c(pData(ESET)$Cohortgroup, "fit.cluster")
  #     }
  #   }
  o <- order(mydata2[,"fit.cluster"])
  mydata2 <- mydata2[o,]
  library(gplots)
  library(RColorBrewer)
  #cols <- palette(brewer.pal(14, "Dark2"))
  hmcol <- colorRampPalette(brewer.pal(11, "RdBu"))(shade)
  
  nco <- ncol(mydata2)
  output <- heatmap(as.matrix(mydata2[,1:ncol(mydata2)-1]), col=rev(hmcol), Rowv=NA, Colv=NA)
  return(output)[[5]]
}  


####HClust
HierClustFunc <- function(data1=gsetMax, pvalcutoff = 0.0001, label="Group"){
  my_palette <- colorRampPalette(c('blue','white','red'))(256)
  scaled <- (exprs(data1))[fData(data1)$Limmaanova< pvalcutoff,]    # scale all but the first column to make information comparable
  scaled.scale <- t(apply(scaled,1,function(x)scale(x)))
  index <- which (colnames(pData(data1))==label)
  colnames(scaled.scale) <- pData(data1)[,index]   #colnames(scaled)
  out.1 <- heatmap.2(scaled.scale,               # specify the (scaled) data to be used in the heatmap
                     cexRow=0.5, 
                     cexCol=0.95,          # decrease font size of row/column labels
                     col = my_palette,     # arguments to read in custom colors
                     #colsep=c(2,4,5),      # Adding on the separators that will clarify plot even more
                     #rowsep = c(6,14,18,25,30,36,42,47), 
                     #sepcolor="black", 
                     #sepwidth=c(0.01,0.01),  
                     scale="none",         # we have already scaled the data 
                     dendrogram="both",#c("both","row","column","None")"None",    # no need to see dendrograms in this one 
                     labRow = fData(data1)$gene_name,
                     trace="none")         # cleaner heatmap
  #return (out.1)
}

SummaryPlot <- function(mydata1=mydata2, hmcol1=hmcol, ESET=NULL){
  outcome <- as.data.frame(apply(mydata1[,-which(names(mydata1) %in% "fit.cluster")], 1, 
                                 function(x)(aggregate(x, by =(pData(ESET)[Label]), FUN=mean)[2])))
  outcome <- as.data.frame(t(outcome))
  
  colnames(outcome) <- (apply(mydata1[1,-which(names(mydata1) %in% "fit.cluster")], 1, 
                              function(x)(aggregate(x, by =(pData(ESET)[Label]), FUN=mean)))[[1]][,1])
  rownames(outcome) <- row.names(mydata1)
  outcome$fit.cluster <- mydata1$fit.cluster
  if (ncol(outcome)==17) outcome <- outcome[,c(7,11:16,1:3, 8:10, 4:6, 17)]
  outputsum <- heatmap(as.matrix(outcome[,1:ncol(outcome)-1]), col=rev(hmcol1), Rowv=NA, Colv=NA)
  
  if (ncol(outcome)==17) outcome <- outcome[,c(1,2,5,8,11,14,3,6,9,12,15,4,7,10,13,16,17)]
  outputsum <- heatmap(as.matrix(outcome[,1:ncol(outcome)-1]), col=rev(hmcol1), Rowv=NA, Colv=NA)
  return(outputsum)
}

individualheatmapfunction <- function(mydata2, hmcol){
  i=0
  for (i in unique(mydata2$fit.cluster)){
    print(sprintf("Cluster %s", i))
    
    # ha_column = c(rep("red", 2), rep("blue", 4),
    #              rep("green", 3), rep("orange", 3))
    
    heatmap(as.matrix(mydata2[grep(i, mydata2$fit.cluster),1:ncol(mydata2)-1]), 
            col=rev(hmcol), Rowv=NA, Colv=NA, main=sprintf("Cluster %s", i), 
            margin=c(4,4))#,
    # ColSideColors=ha_column)
    
    #mtext(side=3, text=c("ILC3" ,"ILC1+/+","ILC+/-", "ILC1" ), 
    #     at=c(0, .2, .4, .6), line=1.25 , col="black", adj=1)
    #legend("topright", legend=c("low","median","high"), 
    #      fill = hmcol[c(shade, shade/2, 1)])
    
  }
}
