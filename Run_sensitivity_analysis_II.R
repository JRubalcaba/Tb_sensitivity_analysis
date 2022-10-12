################################# Sensitivity analyses Code 2 ################################# 
# R code from: "Predicting functional trait variation across climatic gradients in ectotherms"
# Juan G. Rubalcaba, Sidney F. Gouveia, Fabricio Villalobos, Miguel Á. Olalla-Tárraga, Jennifer Sunday.

# This code uses the dataframe generated in Code 1 to compute the sensitivity of thermal performance 
# to variations in body mass,skin absorbance to shortwave radiation, preferred temperature, thermoregulatory ability, 
# and critical thermal limits (CTmax and CTmin).

setwd("") # set working directory

load("Sources/xy.values.RData") # List of cells and lon/lat values
load("Sources/cells.RData")

files <- list.files(pattern = "sensit_matrix_") # run Code 1 first

mat_betas <- data.frame("cells"=rep(NA,length(files)), "beta_M"=NA, "beta_a"=NA, "beta_Topt"=NA, "beta_lambda"=NA, "beta_CTmax"=NA, "beta_CTmin"=NA)
mat_var <- data.frame("cells"=rep(NA,length(files)), "var_M"=NA, "var_a"=NA, "var_Topt"=NA, "var_lambda"=NA, "var_CTmax"=NA, "var_CTmin"=NA)
mat_means <- data.frame("cells"=rep(NA,length(files)),"meanPerf"=NA,"multi_traits"=NA, "multi_microclim"=NA, "meanTb"=NA,"meanTeSun"=NA,"meanTeShade"=NA)
for(i in 1:length(files)){
  load(files[i])
  sensit_matrix_total <- sensit_matrix
  
  # Change in body mass
  M_change <- sensit_matrix_total$M[which(sensit_matrix_total$trait==2)] - sensit_matrix_total$M[which(sensit_matrix_total$trait==1)]
  dev_change_M <- sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==2)] - sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==1)]
  mod_M <- lm((dev_change_M) ~ (M_change))
  
  # Change in absorbance
  a_change <- sensit_matrix_total$a[which(sensit_matrix_total$trait==3)] - sensit_matrix_total$a[which(sensit_matrix_total$trait==2)]
  dev_change_a <- sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==3)] - sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==2)]
  mod_a <- lm((dev_change_a) ~ (a_change))
  
  # Change in preferred temperature
  Topt_change <- sensit_matrix_total$Topt[which(sensit_matrix_total$trait==4)] - sensit_matrix_total$Topt[which(sensit_matrix_total$trait==3)]
  dev_change_Topt <- sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==4)] - sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==3)]
  mod_Topt <- lm((dev_change_Topt) ~ (Topt_change))
  
  # Change in thermoregulatory ability
  lambda_change <- sensit_matrix_total$lambda[which(sensit_matrix_total$trait==5)] - sensit_matrix_total$lambda[which(sensit_matrix_total$trait==4)]
  dev_change_lambda <- sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==5)] - sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==4)]
  mod_lambda <- lm((dev_change_lambda) ~ (lambda_change))
  
  # Change in CTmax
  CTmax_change <- sensit_matrix_total$CTmax[which(sensit_matrix_total$trait==6)] - sensit_matrix_total$CTmax[which(sensit_matrix_total$trait==5)]
  dev_change_CTmax <- sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==6)] - sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==5)]
  mod_CTmax <- lm((dev_change_CTmax) ~ (CTmax_change))
  
  # Change in CTmin
  CTmin_change <- sensit_matrix_total$CTmin[which(sensit_matrix_total$trait==7)] - sensit_matrix_total$CTmin[which(sensit_matrix_total$trait==6)]
  dev_change_CTmin <- sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==7)] - sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==6)]
  mod_CTmin <- lm((dev_change_CTmin) ~ (CTmin_change))
  
  # Check computed relationships
  # plot(dev_change_M ~ (M_change), ylab="Change in performance", xlab="Change in mass", pch=20)
  # abline(mod_M)
  # plot(dev_change_a ~ (a_change), ylab="Change in performance", xlab="Change in absorbance", pch=20)
  # abline(mod_a)
  # plot(dev_change_Topt ~ (Topt_change), ylab="Change in performance", xlab="Change in Tpref", pch=20)
  # abline(mod_Topt)
  # plot(dev_change_lambda ~ (lambda_change), ylab="Change in performance", xlab="Change in lambda", pch=20)
  # abline(mod_lambda)
  # plot(dev_change_CTmax ~ (CTmax_change), ylab="Change in performance", xlab="Change in CTmax", pch=20)
  # abline(mod_CTmax)
  # plot(dev_change_CTmin ~ (CTmin_change), ylab="Change in performance", xlab="Change in CTmin", pch=20)
  # abline(mod_CTmin)
  
  totalTb_change_traits <- sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==7)] - sensit_matrix_total$IntegPerformance[which(sensit_matrix_total$trait==1)]

  # Extract sensitivity coefficients: slope of the relationship between change in performance vs change in trait value
  mat_betas$beta_M[i] <- coef(mod_M)[2]
  mat_betas$beta_a[i] <- coef(mod_a)[2]
  mat_betas$beta_Topt[i] <- coef(mod_Topt)[2]
  mat_betas$beta_lambda[i] <- coef(mod_lambda)[2]
  mat_betas$beta_CTmax[i] <- coef(mod_CTmax)[2]
  mat_betas$beta_CTmin[i] <- coef(mod_CTmin)[2]
  
  # Extract sensitivity coefficients in absolute value
  abs_mod_M <- lm(abs(dev_change_M) ~ abs(M_change))
  abs_mod_a <- lm(abs(dev_change_a) ~ abs(a_change))
  abs_mod_Topt <- lm(abs(dev_change_Topt) ~ abs(Topt_change))
  abs_mod_lambda <- lm(abs(dev_change_lambda) ~abs (lambda_change))
  abs_mod_CTmax <- lm(abs(dev_change_CTmax) ~ abs(CTmax_change))
  abs_mod_CTmin <- lm(abs(dev_change_CTmin) ~ abs(CTmin_change))
  
  mat_var$var_M[i] <- coef(abs_mod_M)[2]
  mat_var$var_a[i] <- coef(abs_mod_a)[2]
  mat_var$var_Topt[i] <- coef(abs_mod_Topt)[2]
  mat_var$var_lambda[i] <- coef(abs_mod_lambda)[2]
  mat_var$var_CTmax[i] <- coef(abs_mod_CTmax)[2]
  mat_var$var_CTmin[i] <- coef(abs_mod_CTmin)[2]

  # Other results
  mat_means$meanTb[i] <- mean(sensit_matrix_total$Tb_mean)
  mat_means$multi_traits[i] <- mean(abs(totalTb_change_traits))
  mat_means$meanTeSun[i] <- mean(sensit_matrix_total$Te_sun)
  mat_means$meanTeShade[i] <- mean(sensit_matrix_total$Te_shade)
  mat_means$meanPerf[i] <- mean(sensit_matrix_total$IntegPerformance)
  
  mat_betas$cells[i] <- sensit_matrix_total$cell[1]
  mat_var$cells[i] <- sensit_matrix_total$cell[1]
  mat_means$cells[i] <- sensit_matrix_total$cell[1]
  
  print(100*i/length(files))
}

# Generate maps of sensitivity values for each trait

betaM_data <- data.frame(xy.values, betaM_data=mat_betas$beta_M)
betaM_map <- rasterFromXYZ(betaM_data)

betaa_data <- data.frame(xy.values, betaa_data=mat_betas$beta_a)
betaa_map <- rasterFromXYZ(betaa_data)

betaTopt_data <- data.frame(xy.values, betaTopt_data=mat_betas$beta_Topt)
betaTopt_map <- rasterFromXYZ(betaTopt_data)

betalambda_data <- data.frame(xy.values, betalambda_data=mat_betas$beta_lambda)
betalambda_map <- rasterFromXYZ(betalambda_data)

betaCTmax_data <- data.frame(xy.values, betaCTmax_data=mat_betas$beta_CTmax)
betaCTmax_map <- rasterFromXYZ(betaCTmax_data)

betaCTmin_data <- data.frame(xy.values, betaCTmin_data=mat_betas$beta_CTmin)
betaCTmin_map <- rasterFromXYZ(betaCTmin_data)

varM_data <- data.frame(xy.values, varM_data=mat_var$var_M)
varM_map <- rasterFromXYZ(varM_data)

vara_data <- data.frame(xy.values, vara_data=mat_var$var_a)
vara_map <- rasterFromXYZ(vara_data)

varTopt_data <- data.frame(xy.values, varTopt_data=mat_var$var_Topt)
varTopt_map <- rasterFromXYZ(varTopt_data)

varlambda_data <- data.frame(xy.values, varlambda_data=mat_var$var_lambda)
varlambda_map <- rasterFromXYZ(varlambda_data)

varCTmax_data <- data.frame(xy.values,varCTmax_data= mat_var$var_CTmax)
varCTmax_map <- rasterFromXYZ(varCTmax_data)

varCTmin_data <- data.frame(xy.values, varCTmin_data=mat_var$var_CTmin)
varCTmin_map <- rasterFromXYZ(varCTmin_data)

df <- data.frame(xy.values, mat_means$meanPerf)
meanPerf_map <- rasterFromXYZ(df)

df <- data.frame(xy.values, mat_means$meanTb)
meanTb_map <- rasterFromXYZ(df)

df <- data.frame(xy.values, mat_means$meanTeSun)
meanTeSun_map <- rasterFromXYZ(df)

df <- data.frame(xy.values, mat_means$meanTeShade)
meanTeShade_map <- rasterFromXYZ(df)

df <- data.frame(xy.values, mat_means$multi_traits)
meanmultitraits_map <- rasterFromXYZ(df)

plot(betaM_map)
plot(betaa_map)
plot(betaTopt_map)
plot(betalambda_map)
plot(betaCTmax_map)
plot(betaCTmin_map)
