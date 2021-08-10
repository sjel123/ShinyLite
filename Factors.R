library(Biobase)
library(gplots)
#library(plotly)
library(ggplot2)

####Functions

PlotExpressionFunction2 <- function(data2=mydata1, generow=2, ESET=gset, res=res, x=1, col="Time_point",
                                    sha=1, LOG=FALSE, SAmple="BOTH"){
 
  if(SAmple=="R848") ESET <- ESET[,ESET$treatment%in% c("Unstimulated", "R848")]
  if(SAmple=="LPS") ESET <- ESET[,ESET$treatment%in% c("Unstimulated", "LPS")]
  data2 <- exprs(ESET)
    plotdata <- as.data.frame((data2[generow,]))
      if(nrow(plotdata)==0){ plotdata = data.frame(t(rep(0, ncol(data2))))
                         names(plotdata) <- names(data2)
      }
    df1Annotation <- pData(ESET)
    
    plotdata2 <- cbind(plotdata, df1Annotation)#, by=0)#), by.x=0, by.y="Sample")
      names(plotdata2)[1]<-"value"
    #print(head(plotdata2))
  if (LOG==TRUE)plotdata2$value = 2^(plotdata2$value)
      plotdata2$Group <- factor(plotdata2[,x])#, levels=c("Wk1", "Wk5", "Wk9", "Fup"))
      plotdata2$Color2 <- factor(plotdata2[,col])#, levels=c("Wk1", "Wk5", "Wk9", "Fup"))
  #print(sprintf("plotdata2$Color2 is  %s",head(plotdata2$Color2)))
      plotdata2$Shape2 <- factor(plotdata2[,sha])#, levels=c("Wk1", "Wk5", "Wk9", "Fup"))
  
    g <- ggplot(plotdata2, aes(x=Group, y=value,  fill=Color2,shape=Shape2))+ 
           geom_boxplot() + #stat_summary(fun.y=mean, geom="line", position=position_dodge(width=0.9),
                            #      aes(group=Shape2 ))+
  #geom_boxplot(aes(group=Time_point, fill=Group)) +
        geom_point(size=3, position=position_dodge(width=.75))+  ggtitle(res$gene_name[generow])
  
    g <- g  +  theme (text=element_text(size=16), 
                      axis.text.x = element_text(angle = 90, hjust =1))#+stat_summary(fun.y = mean, geom="line", size=2)
    g <- g + geom_point(position=position_dodge(width=.75))   
    g <- g + guides(colour = guide_legend(override.aes = list(size=30)))
  #gg <- ggplotly(g)
    p <- ggplot(plotdata2, aes(x=Group, y=value, group=Group, color=Group), shape=as.factor(Shape2))+
          geom_point(size=3)+facet_wrap(~Group)+
          ggtitle(res$gene_name[generow])
    p <- p  + theme (text=element_text(size=16))  #+stat_summary(fun.y = mean, geom="line", size=2)
    #p <- p + geom_line(aes(group=Donor))
  
    q <-ggplot(plotdata2, aes(x=Time, y=value))+#, group=interaction(Time,Treatment), color=Treatment))+ 
      geom_boxplot() +facet_wrap(~ Condition)+
      geom_point(position=position_dodge(width=.75))+ theme (text=element_text(size=16))#
  
    r <- ggplot(plotdata2, aes(x=Time, y=value, group=Group, color=Conditions))+#interaction(Time,Treatment), color=Treatment))+ 
         geom_boxplot() + stat_summary(fun.y=median, geom="line", aes(group=Conditions),
                                       position=position_dodge(width=0.0))+ #Change width to 0.75 to line up with box
         geom_point(position=position_dodge(width=0.75)) +      
         theme (text=element_text(size=16))
  s <- ggplot(plotdata2, aes(x=Time, y=value, group=Group, color=Conditions))+#interaction(Time,Treatment), color=Treatment))+ 
    geom_boxplot() + #stat_summary(fun.y=median, geom="line", aes(group=Conditions),
                                #  position=position_dodge(width=0.75))+ #Change width to 0.75 to line up with box
    geom_point(position=position_dodge(width=0.75)) + 
    geom_line(aes(x=Time, y=value, group=interaction(Donor,Conditions)), position=position_dodge(width=0.0))+
    theme (text=element_text(size=16))+facet_wrap(~Conditions)
  s
  return(list(g,p,q, r,s))
#return(g)
}

FCGraph <- function(generow=2, ESET=gset){
  data2 <- fData(gset)
  plotdata <- data2[generow,]
  print(plotdata)
  plotdata <- as.data.frame(t(plotdata))
  print(plotdata)
  # if(nrow(plotdata)==0){ plotdata = data.frame(t(rep(0, ncol(data2))))
  # names(plotdata) <- names(data2)
  # }
  plotdata$Group <- row.names(plotdata)
  plotdata <- plotdata[25:38,]#grep("FC",plotdata$Group ),]
  names(plotdata) <- c("value", "Group")
  print(head(plotdata))
  plotdata$value <- as.numeric(as.character(plotdata$value))
  print(head(plotdata))
  plotdata$Group <- gsub("stimulated", "", plotdata$Group)
  print(head(plotdata))
  plotdata$Group <- factor(plotdata$Group, levels=c(plotdata$Group[c(2,1,4,3,6,5,8,7,10,9,12,11,14,13)]))
  print(head(plotdata))
  p <- ggplot(plotdata, aes(x=Group, y=value))+geom_bar(stat='identity')+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(p) 
  return(p)
}

invokeRemote = function(proxy, method, args = list()) {
  if (!inherits(proxy, 'datatableProxy'))
    stop('Invalid proxy argument; table proxy object was expected')
  
  msg = list(id = proxy$id, call = list(method = method, args = args))
  
  sess = proxy$session
  if (proxy$deferUntilFlush) {
    sess$onFlushed(function() {
      sess$sendCustomMessage('datatable-calls', msg)
    }, once = TRUE)
  } else {
    sess$sendCustomMessage('datatable-calls', msg)
  }
  proxy
} 

updateSearch = function(proxy, keywords = list(global = NULL, columns = NULL)) {
  global = keywords$global
  if (is.null(global)) {
    keywords['global'] = list(NULL)
  } else {
    if (!is.character(global) || length(global) != 1)
      stop('keywords$global must be a character string')
  }
  columns = keywords$columns
  if (is.null(columns)) {
    keywords['columns'] = list(NULL)
  } else {
    if (is.character(columns)) {
      if (length(columns) == 0) stop(
        'The length of keywords$columns must be greater than zero if it is a character vector'
      )
    } else if (is.list(columns)) {
      if (any(sapply(columns, length) > 1)) stop(
        'keywords$columns should be a list of NULL or character strings'
      )
    } else stop('keywords$columns must be either a character vector or a list')
  }
  invokeRemote(proxy, 'updateSearch', list(keywords))
}


