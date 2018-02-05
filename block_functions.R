library(neuRosim)
library(MASS)
library(spatialfil)
library(nlme)
library(waveslim)
library(shiny)
library(parallel)
library(rmarkdown)

allwave <- c("haar" , "sym8", "db3",  "d4", "mb4", "w4",
             "bs3.1", "fk4", "d6", "fk6", "d8", 
             "fk8", "la8", "mb8", "bl14", "fk14", "d16", "la16", "mb16", 
             "la20", "bl20", "fk22", "mb24")
con <- c(1,-1)

## correlation
## function to analyze mean correlation for single subject
cor_mean <- function(data){
   
  roi1 <- apply(data[[1]],3,mean)
  roi2 <- apply(data[[2]],3,mean)
  
  corresult <- cor(roi1, roi2)
  corresult
}

## function to analyze dw correlation for single subject
cor_dw <- function(data, input, level=1){
  
  waveP <- input$waveP
  waveT <- input$waveT
  
  roicoef <- rep(NA,length(data))
  
  sizeroi <- dim(data[[1]])
  
  tempdw2d <- mydwt.2d(matrix(0,sizeroi[1],sizeroi[2] ), waveP, level)
  size2dLL <- dim(tempdw2d$LL1) # low-low
  
  tempdw1d <- mydwt(rep(0, sizeroi[3]), waveP, level)
  size1dLL <- length(tempdw1d$s1) # low
  
  dw_roidata_2d <- array(NA, c(size2dLL[1],size2dLL[2], size1dLL))
  dw_roidata_2d_temp <- array(NA, c(size2dLL[1],size2dLL[2], input$N.time))
  
  
  coefdw <- matrix(NA, size2dLL[1],size2dLL[2])
  
  dwcount <- 0
  waveform <- dwvariance <- matrix(NA, 2,8)
  allcor <- rep(NA, 8) 
   
  for (j in 1:length(tempdw2d)){
    for (k in 1:length(tempdw1d)){
      
      dwtemp <- list()
      dwcount <- dwcount+1
      
      for (i in 1:length(data)){
        
        roidata <- data[[i]]
        
        for ( t1 in 1: sizeroi[3] ){
          tempdw  <- mydwt.2d(roidata[,,t1], waveP, level )
          dw_roidata_2d_temp[,,t1] <- tempdw[[j]]
        }
        
        for (d1 in 1:size2dLL[1]){
          for (d2 in 1: size2dLL[2]){
            tempdw1d <- mydwt(dw_roidata_2d_temp[d1,d2,] , waveP, level)
            dw_roidata_2d[d1,d2,] <- tempdw1d[[k]]
          }
        }
        
        dwtemp[[i]] <- dw_roidata_2d
        waveform[i,dwcount] <- var(dw_roidata_2d^2)
        dwvariance[i,dwcount] <- var(dw_roidata_2d)
      }
      
      allcor[dwcount] <- cor_mean(dwtemp)
      
    }
  }

  waveform_mean <- apply(waveform,2, mean)
  dwvariance_mean <- apply(dwvariance, 2, mean)
  
  # composite correlation based on waveform 
  as.vector(waveform_mean/sum(waveform)) %*% as.vector(allcor)
  
}


# function to generate independent  data
gen2roi_error_ind <- function(input,X=X, block=TRUE){
  
  dim1 <- input$N.dim1
  dim2 <- input$N.dim2
  alltime <- input$N.time
  roi <- 2 
  SigCov  <- matrix( c(1,input$correlation, input$correlation, 1), ncol=2 )*input$randomsigma
  phi <- input$phi
  beta1 <- c(1+input$effect_size, 1)
  beta2 <- c(2,2)
  
  
  roi1 <- roi2 <- array(NA, c(dim1,dim2, alltime))
  
  # generate 2 roi data
  for (d1 in 1:dim1 ){
    for (d2 in 1:dim2 ){
      temp_v <- mvrnorm(1, rep(0,roi), SigCov)
      roi1[d1,d2,1] = temp_v[1]
      roi2[d1,d2,1] = temp_v[2]
      
      for (t1 in 2:alltime){
        temp_v <- mvrnorm(1, rep(0,roi),SigCov)
        roi1[d1,d2,t1] = phi * roi1[d1,d2,t1-1] + temp_v[1]
        roi2[d1,d2,t1] = phi * roi2[d1,d2,t1-1] + temp_v[2]
      }
    }
  }
  roidata <- list(roi1,roi2)
  if (block==TRUE){
    roidata <- addbeta(roidata, X,beta1, beta2)
  }
  return(roidata)
}


# function to generate gaussian  data
gen2roi_error_gau <- function(input, X=X, block=TRUE){

  GauSigma <- input$GauSigma
  
  dim1 <- input$N.dim1
  dim2 <- input$N.dim2
  alltime <- input$N.time
  roi <- 2 
  SigCov  <- matrix( c(1,input$correlation, input$correlation, 1), ncol=2 )*input$randomsigma
  phi <- input$phi
  beta1 <- c(1+input$effect_size, 1)
  beta2 <- c(2,2)
  
  roi1 <- roi2 <- array(NA, c(dim1,dim2, alltime))
  
  # generate 2 roi data
  for (d1 in 1:dim1 ){
    for (d2 in 1:dim2 ){
      temp_v <- mvrnorm(1, rep(0,roi), SigCov)
      roi1[d1,d2,1] = temp_v[1]
      roi2[d1,d2,1] = temp_v[2]
      
      for (t1 in 2:alltime){
        temp_v <- mvrnorm(1, rep(0,roi),SigCov)
        roi1[d1,d2,t1] = phi * roi1[d1,d2,t1-1] + temp_v[1]
        roi2[d1,d2,t1] = phi * roi2[d1,d2,t1-1] + temp_v[2]
      }
    }
  }
  roidata <- list(roi1,roi2)
  
  if (block==TRUE){
    roidata <- addbeta(roidata, X,beta1, beta2)
  }
  
  roidata <- spatial_smoothing(data=roidata, input)
  
  return(roidata)
}


# function to generate independent error data
gen2roi_error_exp <- function(input, X=X, block=TRUE){
  
  phi_sigma <- input$phi_sigma
  spat_phi <- input$spat_phi
  dim1 <- input$N.dim1
  dim2 <- input$N.dim2
  alltime <- input$N.time
  roi <- 2 
  SigCov  <- matrix( c(1,input$correlation, input$correlation, 1), ncol=2 )*input$randomsigma
  phi <- input$phi
  beta1 <- c(1+input$effect_size, 1)
  beta2 <- c(2,2)
  
  pairdistV <- cbind(rep( 1:dim1, each= dim2 ), 1:dim2)
  distMat <- naive_pdist(pairdistV)
  spat_COV_true <-  phi_sigma*exp(-distMat/spat_phi) 
  chol_spat <- chol(spat_COV_true);
  
  Voxel <- dim1*dim2
  beta_temp_spat1 <- mvrnorm(1, rep(0,Voxel),diag(Voxel))
  beta_temp_spat2 <- mvrnorm(1, rep(0,Voxel),diag(Voxel))
  
  beta_spat1 <- t(chol_spat)%*%beta_temp_spat1
  beta_spat2 <- t(chol_spat)%*%beta_temp_spat2
  
  roi1 <- roi2 <- array(NA, c(dim1,dim2, alltime))
  
  tempcount <- 0
  
  # generate 2 roi data
  for (d1 in 1:dim1 ){
    for (d2 in 1:dim2 ){
      tempcount <- tempcount + 1
      
      temp_v <- mvrnorm(1, rep(0,roi), SigCov)
      roi1[d1,d2,1] = temp_v[1] #+ beta_spat1[tempcount]
      roi2[d1,d2,1] = temp_v[2] #+ beta_spat1[tempcount]
      
      for (t1 in 2:alltime){
        temp_v <- mvrnorm(1, rep(0,roi),SigCov)
        roi1[d1,d2,t1] = phi * roi1[d1,d2,t1-1] + temp_v[1] #+ beta_spat1[tempcount]
        roi2[d1,d2,t1] = phi * roi2[d1,d2,t1-1] + temp_v[2]#+ beta_spat1[tempcount]
      }
    }
  }
  roidata <- list(roi1,roi2)
  #roidata <- addbeta(roidata, X,beta1, beta2)
  
  sizeroi <- dim(roidata[[1]])
  roidata1 <- roidata[[1]] 
  roidata2 <- roidata[[2]]
  
  tempcount <- 0
  
  if (block==TRUE){

    for (d1 in 1:sizeroi[1]){
      for (d2 in 1:sizeroi[2]){
        tempcount <- tempcount + 1
        
        roidata1[d1,d2, ] <- roidata1[d1,d2, ] + 
          as.matrix(X) %*% (beta1 + beta_spat1[tempcount])
        roidata2[d1,d2, ] <- roidata2[d1,d2, ] + 
          as.matrix(X) %*% (beta2 ++ beta_spat2[tempcount])
      }
    }
    
      }
  
  
  roidata <- list(roidata1,roidata2)
  return(roidata)
}

# function to generate identical error data
gen2roi_error_same <- function(input, X=X, block=TRUE){
  
  dim1 <- input$N.dim1
  dim2 <- input$N.dim2
  alltime <- input$N.time
  roi <- 2 
  SigCov  <- matrix( c(1,input$correlation, input$correlation, 1), ncol=2 )*input$randomsigma
  phi <- input$phi
  beta1 <- c(1+input$effect_size, 1)
  beta2 <- c(2,2)
  
  roi1 <- roi2 <- array(NA, c(dim1,dim2, alltime))
  
  temp1 <- temp2 <- rep(NA, alltime)
  
  temp1[1] <- rnorm(1,0,1)
  temp2[1] <- rnorm(1,0,1)
  
  
  for (t1 in 2:alltime){
    temp_v <- rnorm(1, rep(0,roi),SigCov)
    temp1[t1] = phi * temp1[t1-1] + temp_v[1]
    temp2[t1] = phi * temp1[t1-1] + temp_v[2]
  }
  
  
  # generate 2 roi data
  for (d1 in 1:dim1 ){
    for (d2 in 1:dim2 ){
      roi1[d1,d2,] = temp1
      roi2[d1,d2,] = temp2
    }
  }
  roidata <- list(roi1,roi2)
  if (block==TRUE){
    roidata <- addbeta(roidata, X,beta1, beta2)
  }
  return(roidata)
}


# function to do spatial smoothing
spatial_smoothing <- function(data, input){
  GauSigma <- input$GauSigma
  for (i in 1:length(data)){
    roidata <- data[[i]]
    sizeroi <- dim(roidata)
    for ( t1 in 1: sizeroi[3] ){
      roidata[,,t1]  <- applyFilter(x = roidata[,,t1] ,
                                    kernel = convKernel(sigma = GauSigma, k = 'gaussian'))
    }
    data[[i]] <- roidata
  }
  return(data)
}

# function to normalize vector
scalar1 <- function(x) {x / sqrt(sum(x^2))}

# function to normalize data
normalseries <- function(data){
  for (i in 1:length(data)){
    roidata <- data[[i]]
    sizeroi <- dim(roidata)
    for ( d1 in 1: sizeroi[1] ){
      for (d2 in 1:sizeroi[2]){
        roidata[d1,d2,] <- scalar1(roidata[d1,d2,])
      }
    }
    data[[i]] <- roidata
  }
  return(data)
}

# function to do mean approach
block_mean <- function(data, X=X, con=con){
  
  roicoef <- rep(NA,length(data))
  for (i in 1:length(data)){
    roidata <- data[[i]]
    sizeroi <- dim(roidata)
    
    meandata <- rep(NA, sizeroi[3])  
    
    for (t1 in 1: sizeroi[3] ){
      meandata[t1] <- mean(roidata[,,t1])
    }
    alldata <- data.frame(y = meandata, x1 = X[,1], x2 = X[,2])
    m1 <- gls( y  ~ x1 + x2 -1 , data=alldata, correlation=corARMA(p=1,q=0))
    roicoef[i] <- m1$coefficients %*% con
  }
  return(roicoef)
}

# function to do double wavelet
block_dw <- function(data, input, X=X, con=con, level=1){
  
  waveP <- input$waveP
  waveT <- input$waveT
  
  roicoef <- rep(NA,length(data))
  
  sizeroi <- dim(data[[1]])
  
  tempdw2d <- mydwt.2d(matrix(0,sizeroi[1],sizeroi[2] ), waveP, level)
  size2dLL <- dim(tempdw2d$LL1) # low-low
  
  tempdw1d <- mydwt(rep(0, sizeroi[3]), waveP, level)
  size1dLL <- length(tempdw1d$s1) # low
  
  dw_roidata_2d <- array(NA, c(size2dLL[1],size2dLL[2], sizeroi[3]))
  #dw_roidata_last <- array(NA, c(size2dLL[1],size2dLL[2],size1dLL))
  
  coefdw <- matrix(NA, size2dLL[1],size2dLL[2])
  
  #finalcoef <- 
    
  s1 <- mydwt(X[,1], waveP, level)
  s2 <- mydwt(X[,2], waveP, level)
  
  for (i in 1:length(data)){
    roidata <- data[[i]]
    
    for ( t1 in 1: sizeroi[3] ){
      tempdw  <- mydwt.2d(roidata[,,t1], waveP, level )
      dw_roidata_2d[,,t1] <- tempdw$LL1
    }
    
    for (d1 in 1:size2dLL[1]){
      for (d2 in 1: size2dLL[2]){
        tempdw1d <- mydwt(dw_roidata_2d[d1,d2,] , waveP, level)
        
        tempdata <- data.frame(y=tempdw1d$s1, x1 = s1$s1, x2 = s2$s1)
        
        x <- as.matrix(tempdata[,2:3])
        
        coeffi <- t(solve(t(x) %*% x ) %*% t(x) %*% tempdata[,1])
        #m1 <- lm( y ~ x1 + x2 -1 , data = tempdata)
        coefdw[d1, d2] <- coeffi %*% con
      }
    }
    
    
    roicoef[i] <- mean(coefdw)
  }
  return(roicoef)
}

# function to calculate pdist
naive_pdist <- function(A) {
  # A: matrix with obersvation vectors 
  #         (nrow = number of observations)
  #
  # B: matrix with another set of vectors
  #          (e.g. cluster centers)
  result = matrix(ncol=nrow(A), nrow=nrow(A))
  for (i in 1:nrow(A))
    for (j in 1:nrow(A))
      result[i,j] = sqrt(sum( (A[i,] - A[j,])^2 ))
    result
}

# function to add effect size to data roi 1
addbeta <- function(data, X,  beta1=c(2, 1), beta2=c(2,2)){
  
  sizeroi <- dim(data[[1]])
  roidata1 <- data[[1]] 
  roidata2 <- data[[2]]
  
  for (d1 in 1:sizeroi[1]){
    for (d2 in 1:sizeroi[2]){
      roidata1[d1,d2, ] <- roidata1[d1,d2, ] + as.matrix(X) %*% beta1
      roidata2[d1,d2, ] <- roidata2[d1,d2, ] + as.matrix(X) %*% beta2
    }
  }
  data <- list(roidata1,roidata2)
  return(data)
}

# generate independent voxel first
hrf <- read.csv("hrf.csv")
X <- hrf[1:128,]

mywave.filter <- function (name) 
{
  select.sym8 <- function() {
    L <- 16
    g <- c(-0.003382416, -0.0005421323, 0.0316950878, 0.0076074873, 
           -0.1432942384, -0.0612733591, 0.4813596513, 0.7771857517, 
           0.3644418948, -0.0519458381, -0.0272190299, 0.0491371797, 
           0.003808752, -0.0149522583, -0.0003029205, 0.0018899503)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.db3 <- function() {
    L <- 6
    g <- c(0.0352262919, -0.0854412739, -0.13501102, 0.4598775021, 
           0.8068915093, 0.332670553)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.haar <- function() {
    L <- 2
    g <- c(0.707106781186547, 0.707106781186547)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.d4 <- function() {
    L <- 4
    g <- c(0.482962913144534, 0.836516303737808, 0.224143868042013, 
           -0.12940952255126)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.mb4 <- function() {
    L <- 4
    g <- c(0.4801755, 0.8372545, 0.2269312, -0.1301477)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.bs3.1 <- function() {
    L <- 4
    g <- c(0.1767767, 0.5303301, 0.5303301, 0.1767767)
    h <- qmf(g)
    gd <- c(0.3535534, 1.06066, -1.06066, -0.3535534)
    hd <- qmf(g)
    return(list(length = L, hpf = h, lpf = g, dhpf = hd, 
                dlpf = gd))
  }
  select.w4 <- function() {
    L <- 4
    g <- c(-1, 3, 3, -1)/8
    h <- c(-1, 3, -3, 1)/8
    return(list(length = L, hpf = h, lpf = g))
  }
  select.fk4 <- function() {
    L <- 4
    g <- c(0.653927555569765, 0.753272492839487, 0.0531792287790598, 
           -0.0461657148152177)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.d6 <- function() {
    L <- 6
    g <- c(0.332670552950083, 0.806891509311093, 0.459877502118491, 
           -0.135011020010255, -0.0854412738820267, 0.0352262918857096)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.fk6 <- function() {
    L <- 6
    g <- c(0.42791503242231, 0.812919643136907, 0.356369511070187, 
           -0.146438681272577, -0.0771777574069701, 0.0406258144232379)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.d8 <- function() {
    L <- 8
    g <- c(0.230377813307443, 0.714846570548406, 0.630880767935879, 
           -0.0279837694166834, -0.187034811717913, 0.0308413818353661, 
           0.0328830116666778, -0.0105974017850021)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.fk8 <- function() {
    L <- 8
    g <- c(0.3492381118638, 0.782683620384065, 0.475265135079471, 
           -0.0996833284505732, -0.15997809743403, 0.0431066681065162, 
           0.0425816316775818, -0.0190001788537359)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.la8 <- function() {
    L <- 8
    g <- c(-0.0757657147893567, -0.0296355276459604, 0.497618667632563, 
           0.803738751805386, 0.297857795605605, -0.0992195435769564, 
           -0.0126039672622638, 0.0322231006040782)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.mb8 <- function() {
    L <- 8
    g <- rev(c(-0.1673619, 0.01847751, 0.5725771, 0.7351331, 
               0.2947855, -0.1108673, 0.007106015, 0.06436345))
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.bl14 <- function() {
    L <- 14
    g <- c(0.0120154192834842, 0.0172133762994439, -0.0649080035533744, 
           -0.064131289818917, 0.360218460898555, 0.781921593296555, 
           0.483610915693782, -0.0568044768822707, -0.101010920866413, 
           0.0447423494687405, 0.0204642075778225, -0.0181266051311065, 
           -0.0032832978473081, 0.0022918339541009)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.fk14 <- function() {
    L <- 14
    g <- c(0.260371769291396, 0.686891477239599, 0.611554653959511, 
           0.0514216541421191, -0.245613928162192, -0.0485753390858553, 
           0.124282560921513, 0.0222267396224631, -0.0639973730391417, 
           -0.00507437254997285, 0.029779711590379, -0.00329747915270872, 
           -0.00927061337444824, 0.00351410097043596)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.d16 <- function() {
    L <- 16
    g <- c(0.0544158422431049, 0.312871590914303, 0.67563073629729, 
           0.585354683654191, -0.0158291052563816, -0.28401554296157, 
           0.0004724845739124, 0.128747426620484, -0.0173693010018083, 
           -0.0440882539307952, 0.0139810279173995, 0.0087460940474061, 
           -0.0048703529934518, -0.000391740373377, 0.0006754494064506, 
           -0.0001174767841248)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.la16 <- function() {
    L <- 16
    g <- c(-0.0033824159513594, -0.0005421323316355, 0.0316950878103452, 
           0.0076074873252848, -0.143294238351054, -0.0612733590679088, 
           0.481359651259201, 0.777185751699748, 0.364441894835956, 
           -0.0519458381078751, -0.0272190299168137, 0.0491371796734768, 
           0.0038087520140601, -0.0149522583367926, -0.0003029205145516, 
           0.0018899503329007)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.mb16 <- function() {
    L <- 16
    g <- rev(c(-0.0130277, 0.02173677, 0.1136116, -0.0577657, 
               -0.2278359, 0.1188725, 0.6349228, 0.6701646, 0.2345342, 
               -0.05656657, -0.01987986, 0.05474628, -0.02483876, 
               -0.04984698, 0.009620427, 0.005765899))
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.la20 <- function() {
    L <- 20
    g <- c(0.000770159809103, 9.56326707837e-05, -0.0086412992759401, 
           -0.0014653825833465, 0.0459272392237649, 0.0116098939129724, 
           -0.159494278857531, -0.0708805358108615, 0.471690666842659, 
           0.769510037014339, 0.383826761225382, -0.0355367403054689, 
           -0.0319900568281631, 0.049994972079156, 0.0057649120455518, 
           -0.020354939803946, -0.000804358934537, 0.0045931735836703, 
           5.7036084339e-05, -0.0004593294205481)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.bl20 <- function() {
    L <- 20
    g <- c(0.0008625782242896, 0.0007154205305517, -0.0070567640909701, 
           0.0005956827305406, 0.0496861265075979, 0.0262403647054251, 
           -0.121552106157816, -0.0150192395413644, 0.513709872833405, 
           0.766954836501085, 0.340216013511079, -0.0878787107378667, 
           -0.0670899071680668, 0.0338423550064691, -0.0008687519578684, 
           -0.0230054612862905, -0.0011404297773324, 0.0050716491945793, 
           0.0003401492622332, -0.0004101159165852)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.fk22 <- function() {
    L <- 22
    g <- c(0.193896107759957, 0.589452190929428, 0.670084962942026, 
           0.21562984913477, -0.228028855771577, -0.164465715268843, 
           0.11154914372207, 0.110155264934066, -0.0660845167937792, 
           -0.0718416819231261, 0.0435423676255571, 0.0447752121844098, 
           -0.0297428807492741, -0.0259708730890212, 0.020284486066678, 
           0.0129642494110898, -0.0128859905624436, -0.00483843263644019, 
           0.00717380316527169, 0.00036128556221949, -0.00267699163858104, 
           0.000880577368638464)
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  select.mb24 <- function() {
    L <- 24
    g <- rev(c(-2.132706e-05, 0.0004745736, 0.0007456041, 
               -0.004879053, -0.001482995, 0.04199576, -0.002658282, 
               -0.006559513, 0.1019512, 0.1689456, 0.1243531, 0.1949147, 
               0.4581101, 0.6176385, 0.2556731, -0.3091111, -0.3622424, 
               -0.004575448, 0.1479342, 0.01027154, -0.01644859, 
               -0.002062335, 0.001193006, 5.361301e-05))
    h <- qmf(g)
    return(list(length = L, hpf = h, lpf = g))
  }
  switch(name, haar = select.haar(), d4 = select.d4(), mb4 = select.mb4(), 
         w4 = select.w4(), bs3.1 = select.bs3.1(), fk4 = select.fk4(), 
         d6 = select.d6(), fk6 = select.fk6(), d8 = select.d8(), 
         fk8 = select.fk8(), la8 = select.la8(), mb8 = select.mb8(), 
         bl14 = select.bl14(), fk14 = select.fk14(), d16 = select.d16(), 
         la16 = select.la16(), mb16 = select.mb16(), la20 = select.la20(), 
         bl20 = select.bl20(), fk22 = select.fk22(), mb24 = select.mb24(), 
         sym8 = select.sym8(), db3 = select.db3(),
         stop("Invalid selection for wave.filter"))
}


mydwt <- function (x, wf = "la8", n.levels = 4, boundary = "periodic") 
{
  switch(boundary, reflection = x <- c(x, rev(x)), periodic = invisible(), 
         stop("Invalid boundary rule in dwt"))
  N <- length(x)
  J <- n.levels
  if (N/2^J != trunc(N/2^J)) 
    stop("Sample size is not divisible by 2^J")
  if (2^J > N) 
    stop("wavelet transform exceeds sample size in dwt")
  dict <- mywave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"
  y <- vector("list", J + 1)
  names(y) <- c(paste("d", 1:J, sep = ""), paste("s", J, sep = ""))
  for (j in 1:J) {
    W <- V <- numeric(N/2^j)
    out <- .C("dwt", as.double(x), as.integer(N/2^(j - 1)), 
              L, h, g, W = as.double(W), V = as.double(V), PACKAGE = "waveslim")[6:7]
    y[[j]] <- out$W
    x <- out$V
  }
  y[[J + 1]] <- x
  class(y) <- "dwt"
  attr(y, "wavelet") <- wf
  attr(y, "boundary") <- boundary
  return(y)
}


mydwt.2d <- function (x, wf, J = 4, boundary = "periodic") 
{
  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"
  dict <- mywave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"
  z <- matrix(0, m/2, n/2)
  storage.mode(z) <- "double"
  x.wt <- vector("list", 3 * J + 1)
  x.names <- NULL
  for (j in 1:J) {
    out <- .C("two_D_dwt", Image = as.double(x), Rows = m, 
              Cols = n, filter.length = L, hpf = h, lpf = g, LL = z, 
              LH = z, HL = z, HH = z, PACKAGE = "waveslim")[7:10]
    if (j < J) {
      index <- (3 * j - 2):(3 * j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, 
                                   j, sep = ""))
      x <- out[[1]]
      m <- dim(x)[1]
      storage.mode(m) <- "integer"
      n <- dim(x)[2]
      storage.mode(n) <- "integer"
      z <- matrix(0, m/2, n/2)
      storage.mode(z) <- "double"
    }
    else {
      index <- (3 * j):(3 * (j + 1)) - 2
      x.wt[index] <- out[c(2:4, 1)]
      x.names <- c(x.names, sapply(names(out)[c(2:4, 1)], 
                                   paste, j, sep = ""))
    }
  }
  names(x.wt) <- x.names
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  attr(x.wt, "class") <- "dwt.2d"
  x.wt
}



