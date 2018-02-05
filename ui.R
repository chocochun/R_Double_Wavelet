source("block_functions.R")

testinput <- function(){
  input <- list()
  input$N.sim <- 3
  input$N.subj <- 5
  input$N.dim1 <- 10
  input$N.dim2 <- 10
  input$effect_size <- 1
  input$N.time <- 128
  input$waveP <- "db3"
  input$waveT <- "sym8"
  input$phi <- 0.6
  input$spat_cor <- "Independent"
  input$spat_phi <- 2
  input$phi_sigma <- input$GauSigma <- 1.5
  input$correlation <- 0.2
  input$randomsigma <- 3
  input

}

restinput <- function(){
  input <- list()
  input$rest_N.sim <- 3
  input$rest_N.dim1 <- 10
  input$rest_N.dim2 <- 10
  input$rest_N.time <- 128
  input$rest_waveP <- "db3"
  input$rest_waveT <- "sym8"
  input$rest_phi <- 0.6
  input$rest_spat_cor <- "Independent"
  input$rest_spat_phi <- 2
  input$rest_phi_sigma <- input$rest_GauSigma <- 1.5
  input$rest_correlation <- c(0,0.2)
  input$rest_randomsigma <- 3
  input
  
}

# Define UI for application that draws a histogram
shinyUI(navbarPage("Double-wavelet Transform for Task-induced fMRI Data Analysis",
                 
                     tabPanel("R Shiny Simulaiton Tool",
                            
                   fluidPage(                                    
                     fluidRow(

                     column(3,
                            h4("Simulation Setting") ,
      numericInput("N.sim", 
                   label = h5("Number of Simulation"), 
                   value = 100),
      numericInput("N.subj", 
                   label = h5("Number of subjects per Simulation"), 
                   value = 5),
      hr(),
      numericInput("N.dim1", 
                   label = h5("Roi Length"), 
                   value = 10),
      numericInput("N.dim2", 
                   label = h5("Roi Width"), 
                   value = 10),
      h6("Note: Number of voxels = Length * Width"),
      hr(),
      sliderInput("correlation", label = h5("Correlation between two ROIs"), 
                  min = 0, max = 1, value = 0.2),
      
      numericInput("effect_size", 
                   label = h5("Effect size for Activate ROI"), 
                   value = 0.2),
      
      selectInput("N.time",h5("Length of Time series"),
                  choices = c(128), selected = 128)
                     ),

      column(4,
             actionButton("go", "Start Simulation"), 
             hr(),
             tableOutput("table"),
             verbatimTextOutput("value")
             
      ),
      
      # correlation size 
      column(2,
             
             #downloadButton('download', "Download Report"),
             
             
      h4("Double Wavelet"),
      # double wavelet option
      selectInput("waveP",h5("Spatial wavelet"),
                  choices = allwave, selected = "db3"),
      selectInput("waveT",h5("Temporal wavelet"),
                  choices = allwave, selected = "sym8")
      
      ),
      
      column(3,
             h4("Temporal correlation"),
             
             sliderInput("phi", label = h5("Autoregressive One parameter"), 
                         min = 0, max = 1, value = 0.6),
             
             selectInput("spat_cor",h4("Spatial Correlation"),
                         choices = c( "Independent" ,"Exponential", "Identical"), # "Gaussian", 
                         selected = "Exponential"),

             numericInput("randomsigma", 
                          label = h5("Random Variance"), 
                          value = 3), 
             
             
      numericInput("spat_phi", 
                   label = h5("Exponential Decaying parameter"), 
                   value = 2),
      
      numericInput("phi_sigma", 
                   label = h5("Exponential Decaying Sigma"), 
                   value = 1.5),
      
      numericInput("GauSigma", 
                   label = h5("Gaussian Smoothing Sigma"), 
                   value = 1.5)
      )
      ) 
    )
  ),
  
  tabPanel("Matlab GUI",
           tabsetPanel(
             tabPanel("Download",
                      
                      h4("Gui"),
                      
                      downloadButton('download_gui', label = "Download"),
                      
                      
                      h4("Sample Data"),
                      
                      downloadLink('download_subj1', label = "subject1")
                      

                      ),
             tabPanel("Example"),
             tabPanel("About")
             
           )
           )
  
  
))
