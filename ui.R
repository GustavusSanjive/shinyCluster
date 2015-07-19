library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Make Cluster Data and Heat Map"),
  
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
  
      
      tabPanel("Statistics", 
              uiOutput("probeset"), 
              plotOutput("checkdat"), 
              plotOutput("Diagnostics"), 
              tableOutput("Summary"), 
              tableOutput("pval"), 
              plotOutput("Tukey") 
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
  




