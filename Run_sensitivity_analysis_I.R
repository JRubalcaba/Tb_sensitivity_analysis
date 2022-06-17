################################# Sensitivity analyses Code 1 ################################# 
# R code from: "Predicting functional trait variation across climatic gradients in ectotherms"
# Juan G. Rubalcaba, Sidney F. Gouveia, Fabricio Villalobos, Miguel Á. Olalla-Tárraga, Jennifer Sunday.

# This code generates the main dataframe containing the results of the sensitivity analysis, i.e.
# body temperature and thermal performance as a function of different values of body mass,
# skin absorbance to shortwave radiation, preferred temperature, thermoregulatory ability, 
# and critical thermal limits (CTmax and CTmin)

setwd("") # set working directory

load("Sources/xy.values.RData") # List of cells and lon/lat values
load("Sources/cells.RData")

################################# LOAD FUNCTIONS and DATA ################################# 
require(NicheMapR)
require(Rcpp)

# Function to compute behavioral thermoregulation (Rubalcaba et al. 2019 Am Nat) 
sourceCpp("Sources/Thermoregulation_Rubalcaba_etal_2019_AmNat.cpp")

# Functions of the heat transfer model (Supplementary Material - transient heat model)
compute_theta <- function(Ta, Tg, v, RH, M, A){
  l = sqrt(A)
  
  nu = -1.1555e-14*(Ta+273)^3 + 9.5728e-11*(Ta+273)^2 + 3.7604e-08*(Ta+273) - 3.4484e-06
  kf = 1.5207e-11*(Ta+273)^3 - 4.8574e-08*(Ta+273)^2 + 1.0184e-04*(Ta+273) - 3.9333e-04
  
  Re = l * v / nu   
  c = 0.37
  n = 0.7
  Nu = 2 + c * Re^n
  hc = Nu * kf / l
  
  ksub = 0.027  
  ts = 0.025 * (0.001 * M / (3.1416 * 1000))^0.2 
  hg = ksub / ts
  
  RH_prop = RH * 0.01
  ps_a = exp(77.3450 + 0.0057 * (Ta+273) - 7235 / (Ta+273)) / (Ta+273)^8.2
  rho_a = RH_prop * 2.2 * ps_a / (Ta+273)
  
  ps_b = exp(77.3450 + 0.0057 * (Ta+273) - 7235 / (Ta+273)) / (Ta+273)^8.2  # changed Tb by Ta to match analythical model
  rho_b = 2.2 * ps_b / (Ta+273)
  
  Ra = 4 * 0.95 * 5.670367e-8 * (Ta+273)^3
  Rg = 4 * 0.95 * 5.670367e-8 * (Tg+273)^3
  
  theta = (1-0.6) * (Ra + hc) + 0.6 * (Rg + hg)
  
  return(theta)
} 
compute_Te <- function(S, Ta, Tg, v, RH, M, A, a){
  l <- sqrt(A)
  
  # heat transfer
  nu = -1.1555e-14*(Ta+273)^3 + 9.5728e-11*(Ta+273)^2 + 3.7604e-08*(Ta+273) - 3.4484e-06
  kf = 1.5207e-11*(Ta+273)^3 - 4.8574e-08*(Ta+273)^2 + 1.0184e-04*(Ta+273) - 3.9333e-04
  Re = l * v / nu
  c = 0.37
  n = 0.7
  Nu = 2 + c * Re^n
  hc = Nu * kf / l
  ksub = 0.027  
  ts = 0.025 * (0.001 * M / (3.1416 * 1000))^0.2 
  hg = ksub / ts
  Ra = 4 * 0.95 * 5.670367e-8 * (Ta+273)^3
  Rg = 4 * 0.95 * 5.670367e-8 * (Tg+273)^3
  
  # Constant j
  Ad = A * (1-0.6) 
  Ag = A * 0.6
  C = 3.7 * M
  
  j = Ad / C * (a * S + (Ra + hc) * Ta) + Ag / C * Tg * (Rg + hg)

  # constant theta
  theta <- A / C * ((1-0.6) * (Ra + hc) + 0.6 * (Rg + hg))

  # equilibrium
  Te <- j / theta
  return(Te)
} 

# Proportion of direct solar radiation through the canopy (Online Methods)
RAD_canopy <- function(z, RAD, LAI){
  omega_p <- function(ac=0.8, z, LAI){
    Kbc <- sqrt((1 + tan(z*pi/180)^2) / 2.00132)
    omega_p <- exp(-sqrt(ac) * Kbc * LAI)
    return(omega_p)
  }
  # Solar radiation reaching below the canopy
  RAD_p <- RAD * omega_p(0.8, z, LAI)
  return(RAD_p)
}

# Thermal performance function (Deutsch et al. 2008 PNAS)
perf_function <- function(temp, Topt, CTmax, CTmin, Wmax=1){
  sigma = (Topt - CTmin) / 4
  y1 <- exp(-((temp - Topt)/(2*sigma))^2)
  y2 <- 1 - ((temp - Topt) / (Topt - CTmax))^2
  y1[temp > Topt] <- y2[temp > Topt]
  y1[y1<0] <- 0
  return(y1)
}

# Leaf Area Index and shade datasets for shade levels calculations (see Online Methods)
load("Sources/LAI_matrix.RData")
load("Sources/shade_matrix.RData")

################################# TRAIT VALUES #################################
mean_values <- data.frame(M_mean = 50, M_sd = 25,   # log10
                          a_mean = 0.8, a_sd=0.05,
                          Topt_mean = 35, Topt_sd = 4.34,
                          lambda_mean = 1.5, lambda_sd = 0.4)
  
# set parameters for numerical integration (see Rubalcaba et al. 2019 Am Nat for details)

Tmin = -50  # Min and max for integration 
Tmax = 100
m = 100   # number of meshpoints
h = (Tmax - Tmin)/m
meshpts = Tmin + ((1:m) - 1/2) * h # vector of meshpoints
q_i <- c(1,1)

################################# RUN SENSITIVITY ANALYSIS ALGORITHM #################################

require(lhs)

parameters = 6 # Mass, absorbance, Tpref, Lambda (thermreogulatory ability), CTmax, CTmin
iterations = 20

# Set input trait distributions for sensitivity analysis
set.seed(111112)
A <- optimumLHS(n=2*iterations, k=parameters)
B <- array(NA, dim=c(2*iterations,parameters))
B[,1] <- qlnorm(A[,1], meanlog = 3.5, sdlog = 0.8) # Mass
B[,2] <- qnorm(A[,2], mean = 0.8, sd = 0.07)       # Absorbance
B[,3] <- qnorm(A[,3], mean = 35, sd = 4.3)         # Tpref
B[,4] <- qlnorm(A[,4], meanlog = 0.02, sdlog = 0.4)   # Lambda
B[,5] <- qnorm(A[,5], mean = 41.78236, sd = 3.260548) # CTmax
B[,6] <- qnorm(A[,6], mean = 7.639902, sd = 3.743959) # CTmin

curve(dlnorm(x, meanlog = 3.5, sdlog = 0.8), 0, 200, xlab="Body mass (g)", ylab="P density", cex.axis=1.1, cex.lab=1.5)
curve(dnorm(x, mean = 0.8, sd = 0.07), 0.5, 1, xlab="Absorbance (-)", ylab="P density", cex.axis=1.1, cex.lab=1.5)
curve(dnorm(x, mean = 35, sd = 4.3), 20, 50, xlab="T pref (degC)", ylab="P density", cex.axis=1.1, cex.lab=1.5)
curve(dlnorm(x, meanlog = 0.02, sdlog = 0.4), 0, 3, xlab="Lambda (-)", ylab="P density", cex.axis=1.1, cex.lab=1.5)
curve(dnorm(x, mean = 41.78236, sd = 3.260548), 20, 60, xlab="CTmax (degC)", ylab="P density", cex.axis=1.1, cex.lab=1.5)
curve(dnorm(x, mean = 7.639902, sd = 3.743959), -10, 20, xlab="CTmin (degC)", ylab="P density", cex.axis=1.1, cex.lab=1.5)

# Build table for sensitivity analysis
sensit_matrix <- data.frame("x"=NA,"y"=NA,"cell"=NA,"iter"=NA,"trait"=rep(1:7,iterations),"M"=NA,"A"=NA,"height"=NA,"a"=NA,"Topt"=NA,"lambda"=NA, "CTmax"=NA, "CTmin"=NA,"Te_sun"=NA,"Te_shade"=NA,"Tb_mean"=NA,"deviation"=NA,"ThermalOpp"=NA, "IntegPerformance"=NA) 
for(iter in 1:iterations){
  i <- 7 * iter - 6
  j <- 2 * iter - 1
  sensit_matrix$iter[i:(i+6)] <- iter
  
  ########## 1. Set initial values
  sensit_matrix$M[i] <- B[j,1]
  sensit_matrix$A[i] <- 0.0314 * pi * (sensit_matrix$M[i]/1000)^(2/3)
  sensit_matrix$height[i] <- (1 + 0.005 * sensit_matrix$M[i]) * 1e-2
  sensit_matrix$a[i] <- B[j,2]
  sensit_matrix$Topt[i] <- B[j,3]
  sensit_matrix$lambda[i] <- B[j,4]
  sensit_matrix$CTmax[i] <- B[j,5]
  sensit_matrix$CTmin[i] <- B[j,6]

  ########## 2. Change M
  sensit_matrix$M[i+1] <- B[j+1,1]
  sensit_matrix$A[i+1] <- 0.0314 * pi * (sensit_matrix$M[i+1]/1000)^(2/3)
  sensit_matrix$height[i+1] <- (1 + 0.005 * sensit_matrix$M[i+1]) * 1e-2
  sensit_matrix$a[i+1] <- sensit_matrix$a[i]
  sensit_matrix$Topt[i+1] <- sensit_matrix$Topt[i]
  sensit_matrix$lambda[i+1] <- sensit_matrix$lambda[i]
  sensit_matrix$CTmax[i+1] <- sensit_matrix$CTmax[i]
  sensit_matrix$CTmin[i+1] <- sensit_matrix$CTmin[i]
  
  ########## 3. Change a
  sensit_matrix$M[i+2] <- sensit_matrix$M[i+1]
  sensit_matrix$A[i+2] <- sensit_matrix$A[i+1]
  sensit_matrix$height[i+2] <- sensit_matrix$height[i+1]
  sensit_matrix$a[i+2] <- B[j+1,2]
  sensit_matrix$Topt[i+2] <- sensit_matrix$Topt[i+1]
  sensit_matrix$lambda[i+2] <- sensit_matrix$lambda[i+1]
  sensit_matrix$CTmax[i+2] <- sensit_matrix$CTmax[i+1]
  sensit_matrix$CTmin[i+2] <- sensit_matrix$CTmin[i+1]
  
  ########## 4. Change Topt
  sensit_matrix$M[i+3] <- sensit_matrix$M[i+2]
  sensit_matrix$A[i+3] <- sensit_matrix$A[i+2]
  sensit_matrix$height[i+3] <- sensit_matrix$height[i+2]
  sensit_matrix$a[i+3] <- sensit_matrix$a[i+2]
  sensit_matrix$Topt[i+3] <- B[j+1,3]
  sensit_matrix$lambda[i+3] <- sensit_matrix$lambda[i+2]
  sensit_matrix$CTmax[i+3] <- sensit_matrix$CTmax[i+2]
  sensit_matrix$CTmin[i+3] <- sensit_matrix$CTmin[i+2]
  
  ########## 5. Change lambda
  sensit_matrix$M[i+4] <- sensit_matrix$M[i+3]
  sensit_matrix$A[i+4] <- sensit_matrix$A[i+3]
  sensit_matrix$height[i+4] <- sensit_matrix$height[i+3]
  sensit_matrix$a[i+4] <- sensit_matrix$a[i+3]
  sensit_matrix$Topt[i+4] <- sensit_matrix$Topt[i+3] 
  sensit_matrix$lambda[i+4] <- B[j+1,4]
  sensit_matrix$CTmax[i+4] <- sensit_matrix$CTmax[i+3]
  sensit_matrix$CTmin[i+4] <- sensit_matrix$CTmin[i+3]
  
  ########## 6. Change CTmax
  sensit_matrix$M[i+5] <- sensit_matrix$M[i+4]
  sensit_matrix$A[i+5] <- sensit_matrix$A[i+4]
  sensit_matrix$height[i+5] <- sensit_matrix$height[i+4]
  sensit_matrix$a[i+5] <- sensit_matrix$a[i+4]
  sensit_matrix$Topt[i+5] <- sensit_matrix$Topt[i+4] 
  sensit_matrix$lambda[i+5] <- sensit_matrix$lambda[i+4]
  sensit_matrix$CTmax[i+5] <- B[j+1,5]
  sensit_matrix$CTmin[i+5] <- sensit_matrix$CTmin[i+4]
  
  ########## 7. Change CTmin
  sensit_matrix$M[i+6] <- sensit_matrix$M[i+5]
  sensit_matrix$A[i+6] <- sensit_matrix$A[i+5]
  sensit_matrix$height[i+6] <- sensit_matrix$height[i+5]
  sensit_matrix$a[i+6] <- sensit_matrix$a[i+5]
  sensit_matrix$Topt[i+6] <- sensit_matrix$Topt[i+5] 
  sensit_matrix$lambda[i+6] <- sensit_matrix$lambda[i+5]
  sensit_matrix$CTmax[i+6] <- sensit_matrix$CTmax[i+5]
  sensit_matrix$CTmin[i+6] <- B[j+1,6]
}

pairs(sensit_matrix[,6:13]) # we expect no correlations between Mass, abs, Tpref and lambda, CTmax, and CTmin

# Run sensitivity analysis

for(cell in 1:length(cells)){
  # Subset months of activity:
  if(xy.values[cell,2] > 45) months_activ <- 4:8 # temperate-boreal
  if(xy.values[cell,2] > 24 && xy.values[cell,2] < 44) months_activ <- 3:9 # subtropical
  if(xy.values[cell,2] > -23 && xy.values[cell,2] < 23) months_activ <- 1:12 # equatorial
  if(xy.values[cell,2] > -44 && xy.values[cell,2] < -24) months_activ <- c(10,11,12,1,2,3,4) # subtropical (south)
  if(xy.values[cell,2] < -45) months_activ <- c(11,12,1,2,3) # austral
  
  for(iter in 1:iterations){
    i <- 7 * iter - 6
    
    ########## 1. Initial values ##########
    height <- sensit_matrix$height[i]
    
    shade_level <- mean(shade_matrix[cell,months_activ])
    lai_level <- mean(LAI_matrix[cell,months_activ])
    
    tryCatch(
      micro <- micro_global(run.gads=0, loc = xy.values[cell,], minshade = 0, maxshade = shade_level, Usrhyt = height),
      error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    ## Full sun
    metout_sun <- as.data.frame(micro$metout)
    soil_sun <- as.data.frame(micro$soil)
    
    ## Shade
    metout_shade <- as.data.frame(micro$shadmet)
    soil_shade <- as.data.frame(micro$shadsoil)

    months <- unique(metout_sun$DOY)
    metout_sun <- metout_sun[metout_sun$DOY %in% months[months_activ],]
    soil_sun <- soil_sun[soil_sun$DOY %in% months[months_activ],]
    metout_shade <- metout_shade[metout_shade$DOY %in% months[months_activ],]
    soil_shade <- soil_shade[soil_shade$DOY %in% months[months_activ],]
    
    # subset daytime
    daytime <- which(subset(metout_sun)$ZEN < 90)
    metout_sun <- metout_sun[daytime,]
    soil_sun <- soil_sun[daytime,]
    metout_shade <- metout_shade[daytime,]
    soil_shade <- soil_shade[daytime,]

    ## Operative temperatures
    TIME <- subset(metout_sun)$TIME
    MONTH <- subset(metout_sun)$DOY
    RAD_sun  <- subset(metout_sun)$SOLR
    AIRT_sun  <- subset(metout_sun)$TALOC
    SOILT_sun  <- subset(soil_sun)$D0cm
    RH_sun  <- subset(metout_sun)$RHLOC
    V_sun  <- subset(metout_sun)$VLOC
    AIRT_shade <- subset(metout_shade)$TALOC
    SOILT_shade <- subset(soil_shade)$D0cm
    RH_shade <- subset(metout_shade)$RHLOC
    V_shade <- subset(metout_shade)$VLOC
    RAD_0 <- subset(metout_shade)$SOLR
    ZEN <- subset(metout_shade)$ZEN
    RAD_shade <- RAD_canopy(z=ZEN, RAD=RAD_0, LAI=lai_level)
    
    M <- sensit_matrix$M[i]
    A <- sensit_matrix$A[i]
    a <- sensit_matrix$a[i]
    
    Te_sun <- Te_shade <- numeric(length(TIME))
    for(t in 1:length(Te_sun)){
      Te_sun[t] <- compute_Te(S=RAD_sun[t], Ta=AIRT_sun[t], Tg=SOILT_sun[t], v=V_sun[t], RH=RH_sun[t], M=M, A=A, a=a)
      Te_shade[t] <- compute_Te(S=RAD_shade[t], Ta=AIRT_sun[t], Tg=SOILT_shade[t], v=V_shade[t], RH=RH_shade[t], M=M, A=A, a=a)
    }

    hc <- compute_theta(Ta=AIRT_shade, Tg=SOILT_shade, v=V_shade, RH=RH_shade, M=M, A=A)
    
    ## Tb distributions
    lambda <- sensit_matrix$lambda[i]
    Topt <- sensit_matrix$Topt[i]

    mpar <- c(n = 2, hc = hc[1], lambda = lambda, C = 3.6*M, A = A, Topt = Topt, sigma = 1, delta_T = 60) # input data to compute Tb distributions
    tempdist <- array(NA, dim=c(m,length(TIME)))
    for(hour in 1:length(TIME)){
      Te <- c(Te_shade[hour],Te_sun[hour]) # operative temperatures sun and shade
      mpar[2] <- hc[hour]
      P <- array(0, dim = c(m,m)) # matrix to store result
      P <- kernel_P(meshpts, meshpts, Te, q_i, mpar) # kernel P (probability of transition from one temperature to another)
      v <- GetEigenvector(P) # get the eigenvector of the kernel P (this is the Tb distribution)
      tempdist[,hour] <- v/sum(v) # normalize Tb distribution (see Rubalcaba et al. 2019 Am Nat for details)
    }
    Tb_dist <- rowMeans(tempdist)

    ## Extract values
    threshold = 10
    sensit_matrix$Te_sun[i] <- mean(Te_sun)
    sensit_matrix$Te_shade[i] <- mean(Te_shade)
    sensit_matrix$Tb_mean[i] <- sum(Tb_dist * meshpts)
    sensit_matrix$deviation[i] <- sum(Tb_dist * abs(meshpts - Topt))
    sensit_matrix$ThermalOpp[i] <- sum(Tb_dist[which(meshpts > (Topt-threshold) & meshpts < (Topt+threshold))])
    
    CTmax <- sensit_matrix$CTmax[i]
    CTmin <- sensit_matrix$CTmin[i]
    perf <- perf_function(meshpts, Topt=Topt, CTmax = CTmax, CTmin = CTmin)
    sensit_matrix$IntegPerformance[i] <- sum(perf * Tb_dist) * h
    
    ########## 2. Change Mass ##########
    height <- sensit_matrix$height[i+1] # note that animal's height changes as a function of body mass and thus the microclimate model needs to run again for this new height value

    tryCatch(
      micro <- micro_global(run.gads=0, loc = xy.values[cell,], minshade = 0, maxshade = shade_level, Usrhyt = height),
      error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    ## Full sun
    metout_sun <- as.data.frame(micro$metout)
    soil_sun <- as.data.frame(micro$soil)
    
    ## Shade
    metout_shade <- as.data.frame(micro$shadmet)
    soil_shade <- as.data.frame(micro$shadsoil)
    
    # Subset months of activity:
    months <- unique(metout_sun$DOY)
    metout_sun <- metout_sun[metout_sun$DOY %in% months[months_activ],]
    soil_sun <- soil_sun[soil_sun$DOY %in% months[months_activ],]
    metout_shade <- metout_shade[metout_shade$DOY %in% months[months_activ],]
    soil_shade <- soil_shade[soil_shade$DOY %in% months[months_activ],]
    
    # subset daytime
    daytime <- which(subset(metout_sun)$ZEN < 90)
    metout_sun <- metout_sun[daytime,]
    soil_sun <- soil_sun[daytime,]
    metout_shade <- metout_shade[daytime,]
    soil_shade <- soil_shade[daytime,]
    
    ## Operative temperatures
    TIME <- subset(metout_sun)$TIME
    MONTH <- subset(metout_sun)$DOY
    RAD_sun  <- subset(metout_sun)$SOLR
    AIRT_sun  <- subset(metout_sun)$TALOC
    SOILT_sun  <- subset(soil_sun)$D0cm
    RH_sun  <- subset(metout_sun)$RHLOC
    V_sun  <- subset(metout_sun)$VLOC
    AIRT_shade <- subset(metout_shade)$TALOC
    SOILT_shade <- subset(soil_shade)$D0cm
    RH_shade <- subset(metout_shade)$RHLOC
    V_shade <- subset(metout_shade)$VLOC
    RAD_0 <- subset(metout_shade)$SOLR
    ZEN <- subset(metout_shade)$ZEN
    RAD_shade <- RAD_canopy(z=ZEN, RAD=RAD_0, LAI=lai_level)
    
    M <- sensit_matrix$M[i+1]
    A <- sensit_matrix$A[i+1]
    a <- sensit_matrix$a[i+1]
    
    Te_sun <- Te_shade <- numeric(length(TIME))
    for(t in 1:length(Te_sun)){
      Te_sun[t] <- compute_Te(S=RAD_sun[t], Ta=AIRT_sun[t], Tg=SOILT_sun[t], v=V_sun[t], RH=RH_sun[t], M=M, A=A, a=a)
      Te_shade[t] <- compute_Te(S=RAD_shade[t], Ta=AIRT_sun[t], Tg=SOILT_shade[t], v=V_shade[t], RH=RH_shade[t], M=M, A=A, a=a)
    }
    hc <- compute_theta(Ta=AIRT_shade, Tg=SOILT_shade, v=V_shade, RH=RH_shade, M=M, A=A)
    
    ## Tb distributions
    lambda <- sensit_matrix$lambda[i+1]
    Topt <- sensit_matrix$Topt[i+1]
    
    mpar <- c(n = 2, hc = hc[1], lambda = lambda, C = 3.6*M, A = A, Topt = Topt, sigma = 1, delta_T = 60) 
    tempdist <- array(NA, dim=c(m,length(TIME)))
    for(hour in 1:length(TIME)){
      Te <- c(Te_shade[hour],Te_sun[hour])
      mpar[2] <- hc[hour]
      P <- array(0, dim = c(m,m)) 
      P <- kernel_P(meshpts, meshpts, Te, q_i, mpar) 
      v <- GetEigenvector(P) 
      tempdist[,hour] <- v/sum(v) 
    }
    Tb_dist <- rowMeans(tempdist)
    
    ## Extract values
    sensit_matrix$Te_sun[i+1] <- mean(Te_sun)
    sensit_matrix$Te_shade[i+1] <- mean(Te_shade)
    sensit_matrix$Tb_mean[i+1] <- sum(Tb_dist * meshpts)
    sensit_matrix$deviation[i+1] <- sum(Tb_dist * abs(meshpts - Topt))
    sensit_matrix$ThermalOpp[i+1] <- sum(Tb_dist[which(meshpts > (Topt-threshold) & meshpts < (Topt+threshold))])
    
    CTmax <- sensit_matrix$CTmax[i+1]
    CTmin <- sensit_matrix$CTmin[i+1]
    perf <- perf_function(meshpts, Topt=Topt, CTmax = CTmax, CTmin = CTmin)
    sensit_matrix$IntegPerformance[i+1] <- sum(perf * Tb_dist) * h
    
    ########## 3. Change Absorbance ##########

    M <- sensit_matrix$M[i+2]
    A <- sensit_matrix$A[i+2]
    a <- sensit_matrix$a[i+2]
    
    Te_sun <- Te_shade <- numeric(length(TIME))
    for(t in 1:length(Te_sun)){
      Te_sun[t] <- compute_Te(S=RAD_sun[t], Ta=AIRT_sun[t], Tg=SOILT_sun[t], v=V_sun[t], RH=RH_sun[t], M=M, A=A, a=a)
      Te_shade[t] <- compute_Te(S=RAD_shade[t], Ta=AIRT_sun[t], Tg=SOILT_shade[t], v=V_shade[t], RH=RH_shade[t], M=M, A=A, a=a)
    }
    hc <- compute_theta(Ta=AIRT_shade, Tg=SOILT_shade, v=V_shade, RH=RH_shade, M=M, A=A)
    
    ## Tb distributions
    lambda <- sensit_matrix$lambda[i+2]
    Topt <- sensit_matrix$Topt[i+2]
    
    mpar <- c(n = 2, hc = hc[1], lambda = lambda, C = 3.6*M, A = A, Topt = Topt, sigma = 1, delta_T = 60) 
    tempdist <- array(NA, dim=c(m,length(TIME)))
    for(hour in 1:length(TIME)){
      Te <- c(Te_shade[hour],Te_sun[hour]) 
      mpar[2] <- hc[hour]
      P <- array(0, dim = c(m,m)) 
      P <- kernel_P(meshpts, meshpts, Te, q_i, mpar) 
      v <- GetEigenvector(P) 
      tempdist[,hour] <- v/sum(v) 
    }
    Tb_dist <- rowMeans(tempdist)
    
    ## Extract values
    sensit_matrix$Te_sun[i+2] <- mean(Te_sun)
    sensit_matrix$Te_shade[i+2] <- mean(Te_shade)
    sensit_matrix$Tb_mean[i+2] <- sum(Tb_dist * meshpts)
    sensit_matrix$deviation[i+2] <- sum(Tb_dist * abs(meshpts - Topt))
    sensit_matrix$ThermalOpp[i+2] <- sum(Tb_dist[which(meshpts > (Topt-threshold) & meshpts < (Topt+threshold))])
    
    CTmax <- sensit_matrix$CTmax[i+2]
    CTmin <- sensit_matrix$CTmin[i+2]
    perf <- perf_function(meshpts, Topt=Topt, CTmax = CTmax, CTmin = CTmin)
    sensit_matrix$IntegPerformance[i+2] <- sum(perf * Tb_dist) * h
    
    ########## 4. Change Topt ##########

    ## Tb distributions
    lambda <- sensit_matrix$lambda[i+3]
    Topt <- sensit_matrix$Topt[i+3]
    
    mpar <- c(n = 2, hc = hc[1], lambda = lambda, C = 3.6*M, A = A, Topt = Topt, sigma = 1, delta_T = 60) 
    tempdist <- array(NA, dim=c(m,length(TIME)))
    for(hour in 1:length(TIME)){
      Te <- c(Te_shade[hour],Te_sun[hour]) 
      mpar[2] <- hc[hour]
      P <- array(0, dim = c(m,m)) 
      P <- kernel_P(meshpts, meshpts, Te, q_i, mpar) 
      v <- GetEigenvector(P) 
      tempdist[,hour] <- v/sum(v) 
    }
    Tb_dist <- rowMeans(tempdist)
    
    ## Extract values
    sensit_matrix$Te_sun[i+3] <- mean(Te_sun)
    sensit_matrix$Te_shade[i+3] <- mean(Te_shade)
    sensit_matrix$Tb_mean[i+3] <- sum(Tb_dist * meshpts)
    sensit_matrix$deviation[i+3] <- sum(Tb_dist * abs(meshpts - Topt))
    sensit_matrix$ThermalOpp[i+3] <- sum(Tb_dist[which(meshpts > (Topt-threshold) & meshpts < (Topt+threshold))])
    
    CTmax <- sensit_matrix$CTmax[i+3]
    CTmin <- sensit_matrix$CTmin[i+3]
    perf <- perf_function(meshpts, Topt=Topt, CTmax = CTmax, CTmin = CTmin)
    sensit_matrix$IntegPerformance[i+3] <- sum(perf * Tb_dist) * h
    
    ########## 5. Change Lambda ##########
    
    ## Tb distributions
    lambda <- sensit_matrix$lambda[i+4]
    Topt <- sensit_matrix$Topt[i+4]
    
    mpar <- c(n = 2, hc = hc[1], lambda = lambda, C = 3.6*M, A = A, Topt = Topt, sigma = 1, delta_T = 60) 
    tempdist <- array(NA, dim=c(m,length(TIME)))
    for(hour in 1:length(TIME)){
      Te <- c(Te_shade[hour],Te_sun[hour]) 
      mpar[2] <- hc[hour]
      P <- array(0, dim = c(m,m)) 
      P <- kernel_P(meshpts, meshpts, Te, q_i, mpar) 
      v <- GetEigenvector(P) 
      tempdist[,hour] <- v/sum(v) 
    }
    Tb_dist <- rowMeans(tempdist)
    
    ## Extract values
    sensit_matrix$Te_sun[i+4] <- mean(Te_sun)
    sensit_matrix$Te_shade[i+4] <- mean(Te_shade)
    sensit_matrix$Tb_mean[i+4] <- sum(Tb_dist * meshpts)
    sensit_matrix$deviation[i+4] <- sum(Tb_dist * abs(meshpts - Topt))
    sensit_matrix$ThermalOpp[i+4] <- sum(Tb_dist[which(meshpts > (Topt-threshold) & meshpts < (Topt+threshold))])

    CTmax <- sensit_matrix$CTmax[i+4]
    CTmin <- sensit_matrix$CTmin[i+4]
    perf <- perf_function(meshpts, Topt=Topt, CTmax = CTmax, CTmin = CTmin)
    sensit_matrix$IntegPerformance[i+4] <- sum(perf * Tb_dist) * h
    
    ########## 6. Change CTmax ##########

    ## Extract values
    sensit_matrix$Te_sun[i+5] <- mean(Te_sun)
    sensit_matrix$Te_shade[i+5] <- mean(Te_shade)
    sensit_matrix$Tb_mean[i+5] <- sum(Tb_dist * meshpts)
    sensit_matrix$deviation[i+5] <- sum(Tb_dist * abs(meshpts - Topt))
    sensit_matrix$ThermalOpp[i+5] <- sum(Tb_dist[which(meshpts > (Topt-threshold) & meshpts < (Topt+threshold))])
    
    Topt <- sensit_matrix$Topt[i+5]
    CTmax <- sensit_matrix$CTmax[i+5]
    CTmin <- sensit_matrix$CTmin[i+5]
    perf <- perf_function(meshpts, Topt=Topt, CTmax = CTmax, CTmin = CTmin)
    sensit_matrix$IntegPerformance[i+5] <- sum(perf * Tb_dist) * h
    
    ########## 7. Change CTmin ##########
    
    ## Extract values
    sensit_matrix$Te_sun[i+6] <- mean(Te_sun)
    sensit_matrix$Te_shade[i+6] <- mean(Te_shade)
    sensit_matrix$Tb_mean[i+6] <- sum(Tb_dist * meshpts)
    sensit_matrix$deviation[i+6] <- sum(Tb_dist * abs(meshpts - Topt))
    sensit_matrix$ThermalOpp[i+6] <- sum(Tb_dist[which(meshpts > (Topt-threshold) & meshpts < (Topt+threshold))])
    
    Topt <- sensit_matrix$Topt[i+6]
    CTmax <- sensit_matrix$CTmax[i+6]
    CTmin <- sensit_matrix$CTmin[i+6]
    perf <- perf_function(meshpts, Topt=Topt, CTmax = CTmax, CTmin = CTmin)
    sensit_matrix$IntegPerformance[i+6] <- sum(perf * Tb_dist) * h
    
    ####
    print(paste((i+6)/140, "of cell", cell))
  }

  sensit_matrix[,1] <- xy.values[cell,1]
  sensit_matrix[,2] <- xy.values[cell,2]
  sensit_matrix[,3] <- cells[cell]
  
  save(sensit_matrix, file=paste0("sensit_matrix_",cells[cell],".RData"))
}
 

