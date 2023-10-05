###
# Project: Methods Comparison
# Data:    Abrolhos BOSS & BRUV fish
# Task:    Modelling fish abundance w/ FSSGAM
# Author:  Charlotte (Claude)
# Date:    July 2023
##

rm(list=ls())

# libraries----
library(devtools)
# devtools::install_github("beckyfisher/FSSgam_package") # Run once
library(FSSgam)
library(tidyverse)
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(doSNOW)
# devtools::install_github("UWAMEGFisheries/GlobalArchive") # Run once
library(GlobalArchive)
library(ggplot2)
library(corrr)
library(gstat)

#### SET DIRECTORIES AND READ IN DATA ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")

study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

##### READ IN FORMATTED DATA ####
setwd(data_dir)

dat <-  readRDS(paste0(name, sep="_", "dat_cti.rds")) %>% 
  glimpse()

# Set the predictors for modeling
pred.vars <- c("depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended", "sd.relief") 

# Check the correlations between predictor variables
summary(dat[,pred.vars])

correlate(dat[,pred.vars], use = "complete.obs") %>%  
  gather(-term, key = "colname", value = "cor") %>% 
  dplyr::filter(abs(cor) > 0.5)

# Check for transformations
par(mfrow = c(3, 2))
for (i in pred.vars) {
  x <- dat[ , i]
  x = as.numeric(unlist(x))
  hist((x))
  plot((x), main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x + 1))
  plot(log(x + 1))
}

resp.vars <- "CTI"

# Run the full subset model selection
savedir <- out_dir
use.dat <- as.data.frame(dat) 
factor.vars <- c("method") # Method as a factors with 2 levels
out.all     <- list()
var.imp     <- list()

# Loop through the FSS function for each response variable----
for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  #use.dat <- as.data.frame(dat[which(dat$Maturity2 == resp.vars[i]), ])
  Model1  <- gam(CTI ~ s(depth, k=3, bs='cr'),
                 family = gaussian(),  data = use.dat)
  
  model.set <- generate.model.set(use.dat = use.dat,
                                  test.fit = Model1,
                                  pred.vars.cont = pred.vars,
                                  pred.vars.fact = factor.vars,
                                  linear.vars = "depth",
                                  k = 3,
                                  factor.smooth.interactions = F
  )
  out.list <- fit.model.set(model.set,
                            max.models = 600,
                            parallel = T)
  names(out.list)
  
  out.list$failed.models # examine the list of failed models
  mod.table <- out.list$mod.data.out  # look at the model selection table
  mod.table <- mod.table[order(mod.table$AICc), ]
  mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
  out.i   <- mod.table[which(mod.table$delta.AICc <= 2), ]
  out.all <- c(out.all,list(out.i))
  var.imp <- c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2
  
  # plot the best models
  for(m in 1:nrow(out.i)){
    best.model.name = as.character(out.i$modname[m])
    png(file = paste(savedir, paste(name, m, resp.vars[i], "mod_fits.png", sep = "_"), sep = "/"))
    if(best.model.name != "null"){
      par(mfrow = c(3, 1), mar = c(9, 4, 3, 1))
      best.model = out.list$success.models[[best.model.name]]
      plot(best.model,all.terms = T, pages = 1, residuals = T, pch = 16)
      mtext(side = 2, text = resp.vars[i], outer = F)}  
    dev.off()
  }
}
plot(best.model$residuals)
best.model$residuals
# Save model fits and importance scores
names(out.all) <- resp.vars
names(var.imp) <- resp.vars
all.mod.fits   <- do.call("rbind",out.all)
all.var.imp    <- do.call("rbind",var.imp)
write.csv(all.mod.fits[ , -2], file = paste(savedir, paste(name, "CTI","all.mod.fits.csv", sep = "_"), sep = "/"))
write.csv(all.var.imp, file = paste(savedir, paste(name, "CTI","all.var.imp.csv", sep = "_"), sep = "/"))

#* Check top models make sense, particularly for 95% 0s ####

Best1 <- gam(CTI~s(biog, k = 3, bs = "cr") + s(detrended, k=3, bs="cr") + method,
             family=gaussian(), data=use.dat)

summary(Best1)
gam.check(Best1)

