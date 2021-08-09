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



##########Load Data
load("crisa_normalized_eset.rds")
gset <- eset
gset$group <- paste(gset$treatment, gset$Visit, sep="_")
#Test Row names and column names match  
table(row.names(pData(gset)) == colnames(exprs(gset)))

library(hgu133plus2.db)

select(hgu133plus2.db, c("1007_s_at","1053_at"), c("SYMBOL", "GENENAME")) ##  This is just a trying example
PROBES<- as.character(row.names(exprs(gset)))
OUT <- select(hgu133plus2.db,keys= PROBES, columns=c("SYMBOL", "GENENAME"))

table(table(OUT$PROBEID))
OUT <- OUT[!(duplicated(OUT$PROBEID)),]

table(row.names(exprs(gset))==OUT$PROBEID)
fData(gset) <- OUT
colnames(fData(gset))[2] <- "gene_name"

save(gset, file="Crissa.gset")

