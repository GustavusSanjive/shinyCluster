library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Interactive Heat Maps"),
  
  fluidRow(),
  
  mainPanel( 
    selectInput(inputId = "datasets",label = "Please select a geneset", 
                choices = names(GeneChoice), selected = NULL), 
    
    radioButtons("checkpatient", "Choose patients", choices=Group_count[ ,1])
    
  ),
  
  
    tabsetPanel(
      tabPanel( "Cluster Data",
                tableOutput("table")
      ),
      
      
      tabPanel("HeatMap",
               d3heatmapOutput("heatmap") 
      ),   
  
      
      tabPanel("HeatMap2",
               plotOutput("heatmap2"),
               
               
               sliderInput("upper", label = h4("upper level of expression"),
                                                 min = 0, max = 5, step = 0.1, value = 2),
               
               br(),
               
               sliderInput("lower", label = h4("lower level of expression"),
                           min = -5, max = 0, step = 0.1, value = -2),
               plotOutput("bwplot_meancent")
               
               
      ),   
      
      tabPanel("Genome Browser",
               
               uiOutput("Gene_select"),
               uiOutput("chooseProbeset"),
               tableOutput("Biomart"),
               tableOutput("chrom_table"),
               uiOutput("slider1"),
               uiOutput("slider2"),
               plotOutput("chrom_plot"),
               uiOutput("chrom_start")
      ),
      
      tabPanel("Statistics", 
              uiOutput("probeset"), 
              plotOutput("checkdat"), 
              plotOutput("Diagnostics"),
              plotOutput("Stdresid"),
              plotOutput("cookdist"),
              plotOutput("detect_out"),
              tableOutput("Summary"), 
              tableOutput("pval"),
              tableOutput("meanvals"),
              tableOutput("samplesize")
              ##plotOutput("Tukey") 
      ),
      
      tabPanel( "Test Power",
                
                sliderInput("mu", label = h4("Mu"),
                            min = 5, max = 7, step = 0.1, value = 6),
                
                br(),
                
                sliderInput("sigma", label = h4("Sigma"),
                            min = 0, max = 4, step = 0.1, value = 2),
                
                br(),
                
                sliderInput("n", label = h4("N sample points"),
                            min = 5, max = 305, step = 5, value = 50),
                
                br(),
                
                plotOutput("mPower")       
      )
      
        )
            
      )
  
  )
  




