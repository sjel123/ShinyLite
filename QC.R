#################
################


###############
##############
library(limma)
library(ggplot2)
library(edgeR)
library(reshape)
#library(plotly)
library(plyr)
library(Biobase)
library(gplots)
library(preprocessCore)


#if expressset already created then skip to set gset to expressionset and skip to
####Study Design
gset <- readRDS("GSE89408/expression.rds")
gset$Group <- gset$condition
fData(gset)$maxExpr <- apply(exprs(gset), 1, max)
Cutoff <- 3.7 #10 used to generate GSETMax and remove low level expression.  For cytoreason cell data set to zero

##########Load Data
DF <- read.table("Data/Gene-count-table.txt", sep="\t",
                 header=TRUE, check.names=FALSE, stringsAsFactors = FALSE)
DF[DF=="N.F."]<-0
DF[is.na(DF)]

SampleAnnotation <- read.table("Data/sample.annotation.txt", 
                               header=TRUE, check.names=TRUE, sep="\t", stringsAsFactors = FALSE)
row.names(SampleAnnotation) <- SampleAnnotation$Sample

SampleAnnotation$Treatment <- SampleAnnotation$biospecimenName
  SampleAnnotation$Treatment <- gsub("Plate-[1:2]|_no stim_|_no stim _|_replicate|_[1-3]0min_stim_|[0-9]+h ","", SampleAnnotation$Treatment)
  SampleAnnotation$Treatment <- gsub("_stim_", "", SampleAnnotation$Treatment)
SampleAnnotation$Plate <- SampleAnnotation$biospecimenName
  SampleAnnotation$Plate <- gsub("_[A-z0-9 ]+", "", SampleAnnotation$Plate)
SampleAnnotation$Donor <- SampleAnnotation$subjectName
SampleAnnotation$Donor <-gsub("Donor |_[0-9]+", "", SampleAnnotation$Donor)
SampleAnnotation$Time <- SampleAnnotation$biospecimenName
  SampleAnnotation$Time <-gsub("Plate-[1:2]|_replicate|_|DMSO|cpd1674| stim|stim","", SampleAnnotation$Time)
  SampleAnnotation$Time <- rep(c(0,10,30,120, 240, 480, 1440),each=8)
SampleAnnotation$Group <- paste(SampleAnnotation$Treatment, SampleAnnotation$Time, sep="_")

GeneCount <- read.table("Data/fc-counting-summary.txt", 
                        header=TRUE, check.names=FALSE)
MapSummary <- read.table("Data/star-mapping-summary.txt", 
                         header=TRUE, check.names=FALSE)
ExpressedGenes <- read.table("Data/expr_count_RPKM.txt", 
                             header=TRUE, check.names=FALSE)

##################
#Create Expression set
##################
gset <- new("ExpressionSet", exprs = ((as.matrix(DF[,9:ncol(DF)]))))
fData(gset) <- DF[,1:8]
#Add max value per gene
fData(gset)$maxExpr <- apply(exprs(gset), 1, max)
pData(gset) <- data.frame(SampleAnnotation)
pData(gset) <- data.frame(pData(gset), GeneCount[,2:5])
pData(gset) <- data.frame(pData(gset), MapSummary[,2:5])
pData(gset) <- data.frame(pData(gset), ExpressedGenes[,2:5])

  row.names(pData(gset)) <- gset$SAMPLE_ID
#Test Row names and column names match  
table(row.names(pData(gset)) == colnames(exprs(gset)))

####If gene_name doesnt exist that add row_name to gene_name variable
Val <- grep("gene_name", names(fData(gset)))
if(identical(Val, integer(0))) {fData(gset)$gene_name <- row.names(fData(gset))}
##################
#Study Design
##################
textplot(table(pData(gset)$Group))
textplot(table(pData(gset)$Group, pData(gset)$TMT))
title("Study Design",xlab = "Time_point", ylab = "Treatment")

##############################
####ToDO ADD ENTREZ IDS#######
##############################


##################
#Normalize Data
##################
Factor <- factor(gset$Group)
design <- model.matrix(~0+Factor)

NormalizeData <- function(data, design=design){
  dge <- DGEList(counts=data)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=TRUE)
  return(v)
}

DF_Norm <- NormalizeData((exprs(gset)), design)
DFS<- DF_Norm$E
DFSMAX <- DFS[apply(DFS, 1, max) >4,]


###########################
########Update GSET with Normalized Values
###########################
assayDataElement(gset, "raw_data")<- exprs(gset)

exprs(gset) <-(DFS)
# For Cytoreason analysis no need for cutoffs
gsetMax <- gset[fData(gset)$maxExpr>Cutoff]



###########################
########QC
###########################
source ("QCFunction.R")
source ("RNASeqQC2/QCFunction.R")
source ("ClusteringFunctions.R")
source ("RNASeqQC2/RegulatedGenesFunction.R")

#########SampleAnnotation
library(knitr)
kable(pData(gset))

###########Scatterplot
Sample1 <- colnames(exprs(gset))[1]
Sample2 <- colnames(exprs(gset))[2]
ScatterFunction (exprs(gset),Sample1,Sample2)
ScatterFunction (exprs(gset)[fData(gset)$maxExpr>Cutoff,],Sample1,Sample2)
ScatterFunctionMatrix (exprs(gset)[fData(gset)$maxExpr>Cutoff,])

##########BoxPlot
boxplot(exprs(gset))
boxplot(exprs(gset)[fData(gset)$maxExpr>Cutoff,])

par(mfrow=c(1,1))
#Histograms
HistFunc(exprs(gset))
HistFunc(exprs(gset)[fData(gset)$maxExpr>Cutoff,])

#PCA PLOT
PCAFunction(gset[fData(gset)$maxExpr>Cutoff,],labelCol="none",  shape1 ="Treatment", color1="Time")
PCAFunction(gset[fData(gset)$"Limmaanova" <0.001],labelCol="none",  shape1 ="Group", color1="Group")

#Correlation heatmap

CorFunction(gset[fData(gset)$maxExpr>Cutoff,], orderby = "Group", Labels="Group")
CorFunction(gset[fData(gset)$"Limmaanova" <0.001&fData(gset)$"gene_type"=="protein_coding",], orderby = "Methyl", Label="Group")
CorFunction(gset[fData(gset)$"Limmaanova" <0.001,], orderby = "Group", Label="Group")

###Density Plots
#Histograms
DensityPlotFunction (gset)
DensityPlotFunction (gsetMax)
DensityPlotFunction (DFSMAX)

## Identification of differential genes
#####Limma Aanlysis to identify differential expressed genes
CONtrasts <- c("DMSO_0 - cpd1674_0",
              "DMSO_10 - cpd1674_10",
              "DMSO_30 - cpd1674_30",
              "DMSO_120 - cpd1674_120",
              "DMSO_240 - cpd1674_240",
              "DMSO_480 - cpd1674_480",
               "DMSO_1440 - cpd1674_1440")





efit1 <- LimmaFunction()[[1]] 
cont.wt1 <- LimmaFunction()[[2]]
save(efit1, cont.wt1, file="Data/efit1.RData")    
############

topTable(efit1, coef=1)
topTable(efit1, coef=2)
topTable(efit1, coef=3)
topTable(efit1, coef=7)
maxExpr<- fData(gset)$maxExpr>0
plot(x=topTable(efit1, coef=1, n="inf")[maxExpr,10],
     y=topTable(efit1, coef=2, n="inf")[maxExpr,10])


#Volanoplot
volcanoplot(efit1, coef=5,  highlight=20, names=efit1$genes$gene_name)
volcanoplot(efit1, coef=1,  highlight=20, names=efit1$genes$parent)
#MvsA plot
limma::plotMA (efit1,coef=5, xlab="Average log-Expression", ylab= "Log Fold Change")
o <- order(efit1$p.value[,1])
x <- efit1$Amean
y <- efit1$coefficients[,1]
G <- efit1$genes$parent

text(x[o[1:20]], y[o[1:20]], labels=G[o[1:20]])


###Combine limma output pvalue into single table  

test1 <- CombineLimmaOutput(efit = efit1, CONTRASTs = cont.wt1)

#### Add output to gset

Index <- which(names(fData(gset))=="maxExpr")
fData(gset) <- data.frame(fData(gset)[1:Index], as.data.frame(test1))
#### Save gset
save(gset, gsetMax, cont.wt1, file="Data/Vav1.RData")


###
#Regulated Genes Table
###
source("RegulatedGenesFunction.R")
RegGeneTable(maxExpr = 3.0, pVal= 0.0001, FC=1.4)

#Kmeans clustering

KmeansOutput <- KmeansClusterFun(gset,n=5, shade=50, ESET=gset, orderBy="Group",PCutoff=0.00005) 
KmeansOutput <- KmeansClusterFun(gsetPC,n=5, shade=50, ESET=gset, orderBy="TMT",PCutoff=0.0005101) 
KmeansOutput <- KmeansClusterFun(gset,n=5, shade=50, ESET=gset, orderBy="Group",PCutoff=0.05) 

###################################################

gsetPC <- gset[fData(gset)$maxExpr>0, ]
HCData <- HierClustFunc(data1=gsetPC, pvalcutoff = 0.000011, label="Group")


HCData <- row.names(topTable(efit1, coef=1, n=240))
HierClustFunc(data1=gset[HCData,], pvalcutoff = 1)


#QC MappingPlots
##Gene Rate
DF <- pData(gset)
#DF$Treatment <- DF$Location_sampled_description
#  DF$Treatment <- gsub("__Intact_Tissue", '', DF$Treatment)
ggplot(DF, aes(x=Group, y=Gene_Rate, color= Treatment))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(position=position_jitterdodge(dodge.width=0.75, jitter.width=.15))+
  ggtitle("Gene_Rate")+
  theme(axis.text.x = element_text(angle = 90, hjust =1))

##Number of Expressed Gene
DF1 <- DF[,c("Group", "X0", "X0.5", "X1", "X10", "Treatment")]
DF1.melt <- melt(DF1)

ggplot(DF1.melt, aes(x=Group, y=value, color= Treatment))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(position=position_jitterdodge(dodge.width=0.75, jitter.width=.15))+
  ggtitle("Expressed Genes")+facet_wrap(~variable, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust =1))

##Gene Rate
DF2 <- DF[,c("Group", "Uniq_Rate", "Multi_Rate", "Unmap_Rate", "Treatment")]
DF2.melt <- melt(DF2)
ggplot(DF2.melt, aes(x=Group, y=value, color= Treatment))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(position=position_jitterdodge(dodge.width=0.75, jitter.width=.15))+
  ggtitle("GeneRate")+facet_wrap(~variable, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust =1))

####Snp
SNP <- read.table("Data/RNASeq-snp-corr.txt", header=T, sep="\t")
SNP$Pair <- make.names(SNP$Pair)
SNP.melt <- melt(SNP)
SNP.melt$variable <- make.names(SNP.melt$variable)
DonorData <- pData(gset)[,c("Sample", "Subject_Number")]
DonorData$Sample <- make.names(DonorData$Sample)  
SNP.melt <- merge(SNP.melt, DonorData, by.x="Pair", by.y="Sample", all.x=TRUE)
SNP.melt <- merge(SNP.melt, DonorData, by.x="variable", by.y="Sample", all.x=TRUE)
SNP.melt <-SNP.melt[order(SNP.melt$Donor.x, SNP.melt$Donor.y),]
SNP.melt$variable <- factor(SNP.melt$variable, levels = unique(SNP.melt$variable[order(SNP.melt$Donor.x)]))
SNP.melt$Pair<- factor(SNP.melt$Pair, levels = unique(SNP.melt$Pair[order(SNP.melt$Donor.y)]))
ggplot(SNP.melt,aes(x=factor(variable), y=factor(Pair), fill=value))+geom_tile()+
  scale_fill_gradient(low="white", high="red")+
  scale_x_discrete(labels=gset$Donor[order(gset$Donor)])+
  scale_y_discrete(labels=gset$Donor[order(gset$Donor)])
