library(Biobase)
library(gplots)
#library(plotly)
library(ggplot2)

####Functions

PlotExpressionFunction2 <- function(data2=mydata1, generow=2, ESET=gset, res=res, x=1, col="Time_point",
                                    sha=1, LOG=FALSE){
  data2 <- exprs(ESET)
    plotdata <- as.data.frame((data2[generow,]))
      if(nrow(plotdata)==0){ plotdata = data.frame(t(rep(0, ncol(data2))))
                         names(plotdata) <- names(data2)
      }
    df1Annotation <- pData(ESET)
    
    plotdata2 <- cbind(plotdata, df1Annotation)#, by=0)#), by.x=0, by.y="Sample")
      names(plotdata2)[1]<-"value"
      
      plotdata2 $treatment <- factor( plotdata2 $treatment, levels=c("NONLESIONAL", "Vehicle", "Crisaborole"))
      plotdata2 $Visit <- factor( plotdata2 $Visit, levels=c("DAY1", "DAY8", "DAY15"))

    #############Upsated
      #pd = ggplot2::position_jitterdodge(dodge.width = .75, jitter.width = 0.3, seed = 1)
      plotdata2$Jittervalue <- jitter(as.numeric(plotdata2$Visit ))
    g <- ggplot(plotdata2, aes(x=Jittervalue, y=value,  fill=group, group=Visit))+ 
      geom_boxplot(outlier.shape = NA, alpha=0.2) +geom_point(aes(color=Visit))+ geom_line(aes(group=Subject.Number), alpha=0.3)
    g <- g  +  facet_grid(~treatment)+
        ggtitle(res$gene_name[generow]) 
    g <- g  +  theme (text=element_text(size=16), 
                      axis.text.x = element_text(angle = 90, hjust =1))
    #g <- g + geom_point(position=position_dodge(width=.75))   
    g <- g + guides(colour = guide_legend(override.aes = list(size=30))) + theme(legend.position = "none")
    # Change the x axis name
    g <- g + scale_x_discrete(name ="Visit", 
                         limits=c("Day1","Day8","Day15"))
    g
return(g)
}




