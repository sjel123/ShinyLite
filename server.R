

library(shiny)
library(DT)   
#library(lubridate)
library(ggplot2)
#library(plotly)
library(Biobase)
source("./Factors.R")
library(gplots)
library(tmod)
library(reshape)
#source("./HeatmapGeneSets.R")
################################################
#Load Data.  This is run only once when the server.R is called

load("/app/Shiny/Vav1/ShinyLite/Data/Vav1.RData")

# Set variables for down stream analysis
mydata1 <- NULL#exprs(gset)
gsetraw <- gset
  exprs(gsetraw)<-  as.matrix(assayDataElement(gset, "raw_data"))
mydata2 <-  NULL#exprs(gsetraw)

#Create table of results 
res <- fData(gset)#[c(1,2,9:ncol(fData(gset)))]
index <- grep("PVal|anova|FC", colnames(res))

res[,index] <- signif(res[,index],2)

#eliminate some data to make it run load faster
fData(gset)<-fData(gset)[,1:2]
fData(gsetraw)<-fData(gsetraw)[,1:2]

##End Load Data  
##################################################

shinyServer <- function(input, output, session) {  
   
  #Select default variables to shown in the results table
  observe({
    choices2 <-  names(res)
    updateCheckboxGroupInput(session, "show_vars1",
                             label = "Select Columns to Include",
                             choices = choices2,
                             selected = c("gene", "gene_name", "maxExpr", "Limmaanova"                             
                             ) , 
                             inline = TRUE)
  })
    
    #Used to remove columns rows from data table() based on categories selected 
    #from show_vars1 checkbox on client
    Names <- reactive({
      NameLabels <- (input$'show_vars1')
    NameLabels
    }) 
    
  ##Used to select rows based on slider inputs
  ## Currently the only option is All
  ## This is useful to have a small number of default genes in the results table
  RowValue <- reactive({
    Value <- as.numeric(row.names(unique(res[,Names()])))
    if(input$action_selectiontype == "TWO") Value <- c(792, 796, 797, 798)#c(20990, 3385,7582, 10233,8432)
    Value
  })
  
  #Set Default Input value for graphs x axis, Color and Shape Variables
  updateRadioButtons(session, "Order", choices=names(pData(gset)), selected="Time", inline = TRUE)
  Order <- reactive ({
    input$Order
  })
  updateRadioButtons(session, "Color", choices=names(pData(gset)), selected="Time", inline = TRUE)
  Color <- reactive ({
    input$Color
  })
  updateRadioButtons(session, "Shape", choices=names(pData(gset)), selected="Treatment", inline = TRUE)
  Shape <- reactive ({
    input$Shape
  })

output$text1 <- renderText({
  paste("Input$order value is ", input$Order)
  paste("Input$color value is ", input$Color)
  paste("Input$shape value is ", input$Shape)
})

  # turn input selection to reactive
  Proteinnames <- reactive({
    NameLabels <- (input$'protein')
    print(sprintf("Proteinnames %s", NameLabels))
    NameLabels
  })
  
  
  #  if (input$go=="PullDown"){
  #  Make Datatable a reactive input  
  DataForTable <- reactive({
    if (input$go==1){DAT <- res[which(Proteinnames() == res$gene_name),]
                     print(sprintf("DAT is %s", DAT))
                     print(DAT)}
    if (input$go==2) DAT <- res
    DAT[,Names()]
  })

  #Set DataTable options
  output$x1 = DT::renderDataTable(DataForTable(), server = T, escape=T, selection = 'single', options=list(
    lengthMenu = list(c(5, 10, 15, -1), c('5', '10', '15', 'All'))))
  
  SW <- reactive({ 
    s = NULL
    s = as.numeric(rownames(DataForTable()[RowValue(),]))[input$x1_rows_selected]
      print(sprintf("Names() %s", Names()))
      print(sprintf("input$x1_rows_selected %s", input$x1_rows_selected))
    if (is.null(s)) s=2  
    if (length(s)==0) s=2 
      print(sprintf("row selected %s", s))

    s
  })
  
  LOG <- reactive ({
    value = FALSE
    if (input$Log=="Log")  {value=TRUE}
    value
  })
  
  
  updateSelectInput(session, "protein", choices=res$gene_name[1:10], selected="MTOR")
  
  
  output$main_plot <-  renderPlot({
    generow <- as.numeric(SW())
    generow1 <- which(generow==as.numeric(row.names(exprs(gset))))
    print(sprintf("main plot SW() %s", generow))
   # print(sprintf("mydata1 %s", mydata1))
   print(PlotExpressionFunction2(mydata1, generow = generow1, x=Order(), 
                                 ESET=gset, res=res, col=Color(),
                                 sha=Shape(), LOG=LOG())[1])
    print("Complete1")
  })
  
 

output$main_plot6 <-  renderPlot({
  generow <- as.numeric(SW())
  generow1 <- which(generow==row.names(exprs(gset)))
  print(sprintf("main plot for plot 6 SW() %s", generow))
  # print(sprintf("mydata1 %s", mydata1))
  print(PlotExpressionFunction2(mydata1, generow = generow1, x=Order(), 
                                ESET=gset, res=res, col=Color(),
                                sha=Shape(), LOG=TRUE)[1])
  print("Complete6")
  
 
})


dataset <- reactive ({
  x <- which(SW()==as.numeric(row.names(exprs(gset))))
    PlotData <- data.frame(t(exprs(gset[x,])), pData(gset))
  names(PlotData)[1] <- c("Value")
  PlotData
})


output$downloadData <- downloadHandler(
  filename = function() { paste("Download", '.csv', sep='') },
  content = function(file) {
    write.csv(res[input$x1_rows_all,Names()], file)
  })

output$downloadDataraw <- downloadHandler(
  filename = function() { paste("Normalized", '.csv', sep='') },
  content = function(file) {
     data1 <- exprs(gset)
     row.names(data1) <- fData(gset)$gene_name
     colnames(data1)<- gset$Treatment
    write.csv(data1, file)
  })
output$downloadDatameta <- downloadHandler(
  filename = function() { paste("Meta", '.csv', sep='') },
  content = function(file) {
    write.csv(pData(gset), file)
  })


}


