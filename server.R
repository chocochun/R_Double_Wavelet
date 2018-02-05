
source("block_functions.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  analysis <- function(input=input,X=X,con=con){
    shiny::isolate({
      source("block_functions.R")
      
    if (input$spat_cor ==  "Independent" ){
      genfunc <- gen2roi_error_ind
    } else if (input$spat_cor ==  "Exponential" ){
      genfunc <- gen2roi_error_exp
    } else if (input$spat_cor ==  "Gaussian" ){
      genfunc <- gen2roi_error_gau
    } else if (input$spat_cor ==  "Identical" ){
      genfunc <- gen2roi_error_same
    }
      
      data <- genfunc(input,X)
      data <- normalseries(data)
      data_smooth <- spatial_smoothing(data=data, input)
      result_mean <- block_mean(data_smooth, X=X, con=con)
      result_dw <- block_dw(data,input, X=X, con=con, level=1)
      result <- rbind( c(result_mean, "Mean") , c(result_dw, "DW") )
    }
    )
      return(result)
    
  }
  
  output$value <- renderPrint({ print(paste0("Simulation ", input$go))  })
  
  resulttable <- eventReactive(input$go, {

    cl <- makeCluster(8, 'PSOCK')
    clusterExport(
    cl, varlist=c("input","analysis"),
    envir=environment())
    
    mat <- parSapply(cl, 1:(input$N.sim*input$N.subj), function(i) {
      source("block_functions.R")
      analysis(input, X, con)
    },simplify = FALSE)
    stopCluster(cl)
    
    #mat <- replicate( input$N.sim *input$N.subj  , analysis(input, X, con)  
    #                  ,simplify=FALSE )
    
    data <- as.data.frame(do.call("rbind", mat))
    colnames(data) <- c("roi1", "roi2", "method")
    data$roi1 <- as.numeric(as.character(data$roi1)) 
    data$roi2 <- as.numeric(as.character(data$roi2)) 
    
    meandata <- data[ data$method == "Mean" ,]
    dwdata <- data[ data$method == "DW" ,]
    
    errorrate <- function(testdata, Nsimu = input$N.sim, Nsub = input$N.subj){
      resultP <- rep(NA,Nsimu)
      
      for (i in 1: Nsimu){
        test1 <- t.test( testdata[ (1:Nsub) + (i-1)*Nsub  ] )
        resultP[i] <- test1$p.value
      }
      return(resultP)
    }
    
   #  mean_roi1 <-  errorrate(meandata$roi1)
   #  mean_roi2 <-  errorrate(meandata$roi2)
   #  dw_roi1 <-  errorrate(dwdata$roi1)
   #  dw_roi2 <-  errorrate(dwdata$roi2)
    
    meantypeI <- mean(errorrate(meandata$roi2) < 0.05)
    meantypeII <- mean(errorrate(meandata$roi1) > 0.05)
    
    dwtypeI <- mean(errorrate(dwdata$roi2) < 0.05)
    dwtypeII <- mean(errorrate(dwdata$roi1) > 0.05)
    
    result <- rbind( c(meantypeI , meantypeII ), c(dwtypeI, dwtypeII) )
    colnames(result) <- c("Type I Error", "Type II Error")
    rownames(result) <- c("Mean-Voxel", "Double-Wavelet")
    
    result
    
      })
  
  
  output$table <- renderTable({
    resulttable()
    }, rownames = TRUE)
  
  output$download <- downloadHandler(
    filename = function() { 
      paste0("DW_Block_" ,Sys.Date(), '.html') 
    },
    content = function(file) {
      
      out = render('block.Rmd', clean = TRUE)
      file.rename(out, file) # move pdf to file for downloading
      
    },
    contentType = 'application/html'
    
  )
  
  
  rest_analysis <- function(input=input,X=X){
    shiny::isolate({
      source("block_functions.R")
      
      if (input$spat_cor ==  "Independent" ){
        genfunc <- gen2roi_error_ind
      } else if (input$spat_cor ==  "Exponential" ){
        genfunc <- gen2roi_error_exp
      } else if (input$spat_cor ==  "Gaussian" ){
        genfunc <- gen2roi_error_gau
      } else if (input$spat_cor ==  "Identical" ){
        genfunc <- gen2roi_error_same
      }
      
      
      
      data <- genfunc(input,X, block=FALSE)
      data <- normalseries(data)
      data_smooth <- spatial_smoothing(data=data, input)
      result_mean <- cor_mean(data_smooth)
      result_dw <- cor_dw(data,input)
      result <- rbind( c(result_mean, "Mean") , c(result_dw, "DW") )
    }
    )
    return(result)
    
  }
  
  
  output$rest_value <- renderPrint({ print(paste0("Simulation ", input$rest_go))  })
  
  rest_resultplot <- eventReactive(input$rest_go, {
    
    allcor <- seq( input$rest_correlation[1],  input$rest_correlation[2], by =0.1)
    
    result <- matrix(NA, ncol=5, nrow = length(allcor)*2 )
    
    restcount <- 0
    
    for (restcor in allcor){
      
      restinput <- list()
      restinput$N.sim <- input$rest_N.sim
      restinput$N.dim1 <- input$rest_N.dim1
      restinput$N.dim2 <- input$rest_N.dim2
      restinput$N.time <- input$rest_N.time
      restinput$waveP <- input$rest_waveP
      restinput$waveT <- input$rest_waveT
      restinput$phi <- input$rest_phi
      restinput$spat_cor <- input$rest_spat_cor
      restinput$spat_phi <- input$rest_spat_phi
      restinput$phi_sigma <- input$rest_phi_sigma
      restinput$GauSigma <- input$rest_GauSigma
      restinput$correlation <- restcor
      restinput$randomsigma <- input$rest_randomsigma
      
      cl <- makeCluster(2, 'PSOCK')
      clusterExport(
        cl, varlist=c("restinput","rest_analysis"),
        envir=environment())
      
      mat <- parSapply(cl, 1:(restinput$N.sim), function(i) {
        source("block_functions.R")
        rest_analysis(restinput, X)
      },simplify = FALSE)
      stopCluster(cl)
      
      data <- as.data.frame(do.call("rbind", mat))
      colnames(data) <- c("correlation","method")
      data$truth <- restcor
      data$correlation <- as.numeric(as.character(data$correlation))
      
      meandata <- data[  data$method== "Mean"  ,]
      dwdata <- data[  data$method== "DW"  ,]
      
      restcount <- restcount + 1
      
      result[restcount, 1] <- "Mean"
      result[restcount, 2] <- mean(meandata$correlation - meandata$truth)
      result[restcount, 3] <- var(meandata$correlation - meandata$truth)
      result[restcount, 4] <- mean(meandata$correlation - meandata$truth)^2 + 
        var(meandata$correlation - meandata$truth)
      result[restcount,5] <- restcor
      
      restcount <- restcount + 1
      
      result[restcount, 1] <- "Double-Wavelet"
      result[restcount, 2] <- mean(dwdata$correlation - dwdata$truth)
      result[restcount, 3] <- var(dwdata$correlation - dwdata$truth)
      result[restcount, 4] <- mean(dwdata$correlation - dwdata$truth)^2 + 
        var(dwdata$correlation - dwdata$truth)
      result[restcount,5] <- restcor
      
      
      plotresult <- as.data.frame(result)
      
      if (sum(is.na(plotresult[,1])) > 0)    
        plotresult <- plotresult[ - which(is.na(plotresult[,1])) ,]

      
      colnames(plotresult) <- c("Method", "Bias", "Variance", "MSE", "Truth")
      plotresult$MSE <- as.numeric(as.character(plotresult$MSE))
      plotresult$Truth <- as.numeric(as.character(plotresult$Truth))
      
     p1 <-  ggplot(plotresult, aes(Truth, MSE, group=Method, color=Method))+ geom_line() +
        ylab("MSE")  + geom_point() + xlim(input$rest_correlation)

     renderPlot({p1})
     
}
    
})
  
  
  output$rest_plot <- renderUI({
    rest_resultplot()
  })
  
  output$download_gui <- downloadHandler(
    filename <- function() {
      paste("dw_gui", "zip", sep=".")
    },
    
    content <- function(file) {
      file.copy("gui/dw_gui.zip", file)
    },
    contentType = "application/zip"
  )
  
  output$download_subj1 <- downloadHandler(
    filename <- function() {
      paste("subj1", "nii", sep=".")
    },
    
    content <- function(file) {
      file.copy("gui/subj1_run1.nii", file)
    }
  )
  
  
})
