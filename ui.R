

#Load Data
################################################
#Load Data.  This is run only once when the server.R is called


##End Load Data                                                                                                                
##################################################
ui <- fluidPage(
  mainPanel(  
   fluidRow(
      img(src = "CSI_3.jpg", align="right", height = 72, width = 336),
      img(src = "2021_I&I_Signature__CMYK.jpg", align="left", height = 72)),
    h5("Link to QUICKRNAseq", a("Link", href=
                                  "http://162.48.173.66:3838/Shiny/IRF5/IRF5/Results/"))
    
    ),#'End Main Panel'
 
  tabsetPanel(  
      tabPanel("Graphs",
        titlePanel("Expression of IRF5KO macrophages"),
          selectInput("protein", "Choose protein", NULL),
            textOutput("protein"),
          radioButtons('go', label= h3("Select Input or Datatable search"),
                     choices = list("PullDown"=1, "DataTable"=2),
                     selected=2),
        
          fluidRow(checkboxGroupInput('show_vars1', 'Columns in Table to show:', choices=NULL, inline=TRUE)),
          fluidRow(radioButtons("action_selectiontype", "Selection type",
                      choices = c("All"),
                      selected = "All", inline = TRUE)),
         h3("Adjusted P Values"),

          fluidRow(column (12, DT::dataTableOutput('x1'),  hr(),fluid=FALSE)
          ),#End Fluid Row

        textOutput ("text1"), 
        downloadButton('downloadData', 'Download Table'),
        downloadButton('downloadDataraw', 'Download Normalized Data'),
        downloadButton('downloadDatameta', 'Download Meta Data'),
        bookmarkButton(),
          fluidRow(radioButtons("Log", "Selection type",
                              choices = c("Log", "Normal"), #"647-TYK2", "554-Jak1", "BOTH"),
                              selected = "Normal", inline = TRUE)),
        fluidRow(
          column(6, h3("Expression - edgeR Norm (Log)"),
                 plotOutput("main_plot",  height = "600px")),
          column(6, h3("Expression - edgrNorm (Standard)"),
                 plotOutput("main_plot6",  height = "600px"))
                 ),#End Fluid Row

          fluidRow(radioButtons("Order", "X Value",
                                choices = c("Time", "Cell", "Treatment"), #"647-TYK2", "554-Jak1", "BOTH"),
                                selected = "Cell", inline = TRUE)),
          fluidRow(radioButtons("Color", "Color type",
                                choices = c("Time", "Cell", "Treatment"), #"647-TYK2", "554-Jak1", "BOTH"),
                                selected = "Cell", inline = TRUE)),
          fluidRow(radioButtons("Shape", "Shape type",
                                choices = c("Time", "Cell", "Treatment"), #"647-TYK2", "554-Jak1", "BOTH"),
                                selected = "Cell", inline = TRUE))
      ), # end tabPanel
    tabPanel("QC"#, 
            # includeMarkdown("./QCmarkdown.md")
             ),#End tabPanel


      tabPanel("Study Info",
             
               tags$iframe(style="height:1400px; width:100%; scrolling=yes", src="IRF5_ESD_Review_0708821.pdf")           
        
            )#End tabPanel
    )#End Tabset Panel

)#End UI






