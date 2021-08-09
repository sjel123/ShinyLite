

#Load Data
################################################
#Load Data.  This is run only once when the server.R is called


##End Load Data                                                                                                                
##################################################
ui <- fluidPage(
  mainPanel(  
    fluidRow(
      img(src = "CSI_3.jpg", align="right", height = 72, width = 336))
    ),#'End Main Panle'
 
  tabsetPanel(  
      tabPanel("Graphs",
        titlePanel("Expression of Criaborole"),
          selectInput("protein", "Choose protein", NULL),
            textOutput("protein"),
          radioButtons('go', label= h3("Select Input or Datatable search"),
                     choices = list("PullDown"=1, "DataTable"=2),
                     selected=2),
        
          fluidRow(checkboxGroupInput('show_vars1', 'Columns in Table to show:', choices=NULL, inline=TRUE)),

          fluidRow(column (12, DT::dataTableOutput('x1'),  hr(),fluid=FALSE)
          ),#End Fluid Row

        downloadButton('downloadData', 'Download Table'),
        downloadButton('downloadDataraw', 'Download Normalized Data'),
        downloadButton('downloadDatameta', 'Download Meta Data'),
        bookmarkButton(),
          fluidRow(radioButtons("Log", "Selection type",
                              choices = c("Log", "Normal"), #"647-TYK2", "554-Jak1", "BOTH"),
                              selected = "Normal", inline = TRUE)),
        fluidRow(
          column(12, h3("Expression - edgeR Norm"),
                 plotOutput("main_plot",  height = "600px"))
                  )#End Fluid Row

      ), # end tabPanel

      tabPanel("Study Info",
               tags$iframe(style="height:400px; width:100%; scrolling=yes", src="IL4_13_p40.pdf")           
               #img(src = "TregDesign"))
            )#End tabPanel
    )#End Tabset Panel

  #uiOutput("DataTable"),
	

)#End UI






