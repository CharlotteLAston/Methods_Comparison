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

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")

study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

##### READ IN FORMATTED DATA ####
setwd(data_dir)

dat <-  readRDS(paste0(name, sep="_", "dat_length.rds")) %>% 
  glimpse()

# Set the predictors for modeling
pred.vars <- c("depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") 

# Format data

dat.response <- dat %>% 
  filter(status %in% "No-take") %>% 
  group_by(method, sample, scientific, Maturity) %>% 
  summarise(Abundance = length(Maturity)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Maturity", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
  mutate(greater_mat_125 = ifelse(is.na(greater_mat_125), 0, greater_mat_125),
         greater_mat_less_125 = ifelse(is.na(greater_mat_less_125), 0, greater_mat_less_125),
         greater_50_less_mat = ifelse(is.na(greater_50_less_mat), 0, greater_50_less_mat)) %>% 
  pivot_longer(cols= c(greater_mat_125, greater_mat_less_125, greater_50_less_mat), names_to="Maturity", values_to="Abundance") %>% 
  mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
  mutate(Maturity = fct_recode(Maturity, "< 50 Length Maturity" = "less_50", ">50 Maturity but < Maturity"="greater_50_less_mat",
                               "> Length Maturity but < 1.25x Maturity" = "greater_mat_less_125",
                               "> 1.25x Maturity"="greater_mat_125")) %>% 
  group_by(method, sample, Maturity) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  glimpse()

dat.preds <- dat %>% 
  dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") %>% 
  distinct()

dat <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample))


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

# Check to make sure Response vector has not more than 80% zeros
unique.vars <- unique(as.character(dat$Maturity))

resp.vars <- character()
for(i in 1:length(unique.vars)){
  temp.dat <- dat[which(dat$Maturity == unique.vars[i]), ]
  if(length(which(temp.dat$Abundance == 0)) / nrow(temp.dat) < 0.9){ # Change here
    resp.vars <- c(resp.vars, unique.vars[i])}
}
resp.vars <- unique(as.character(dat$Maturity))

# Run the full subset model selection
savedir <- out_dir
use.dat <- as.data.frame(dat) 
factor.vars <- c("method") # Method as a factors with 2 levels
out.all     <- list()
var.imp     <- list()

# Loop through the FSS function for each response variable----
for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  use.dat <- as.data.frame(dat[which(dat$Maturity == resp.vars[i]), ])
  Model1  <- gam(Abundance ~ s(depth, k=3, bs='cr'),
                 family = tw(),  data = use.dat)
  
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

# Save model fits and importance scores
names(out.all) <- resp.vars
names(var.imp) <- resp.vars
all.mod.fits   <- do.call("rbind",out.all)
all.var.imp    <- do.call("rbind",var.imp)
write.csv(all.mod.fits[ , -2], file = paste(savedir, paste(name, "mat_groups", "all.mod.fits.csv", sep = "_"), sep = "/"))
write.csv(all.var.imp, file = paste(savedir, paste(name, "mat_groups","all.var.imp.csv", sep = "_"), sep = "/"))

#* Check top models make sense, particularly for 95% 0s ####

use.dat <- dat %>% 
  filter(Maturity %in% c("> 1.25x Maturity"))

Best1 <- gam(Abundance~s(tpi, k = 3, bs = "cr") + method,
             family=tw(), data=use.dat)

summary(Best1)
gam.check(Best1)

use.dat <- dat %>% 
  filter(Maturity %in% c("> Length Maturity but < 1.25x Maturity"))

Best1 <- gam(Abundance~s(macroalgae, k = 3, bs = "cr") + method,
             family=tw(), data=use.dat)

summary(Best1)
gam.check(Best1)

use.dat <- dat%>% 
  filter(Maturity %in% c(">50 Maturity but < Maturity"))

Best1 <- gam(Abundance~ s(macroalgae, k = 3, bs = "cr") + method,
             family=tw(), data=use.dat)

summary(Best1)
gam.check(Best1)


#### MODELS OF ALL LENGTH GROUPS THAN BY SPECIES ####

# Check to make sure Response vector has not more than 80% zeros
dat <-  readRDS(paste0(name, sep="_", "dat_length.rds")) %>% 
  glimpse()

# Set the predictors for modeling
pred.vars <- c("depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") 

# Format data

dat.response <- dat %>% 
  filter(status %in% "No-take") %>% 
  group_by(method, sample, scientific, Maturity) %>% 
  summarise(Abundance = length(Maturity)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Maturity", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
  mutate(greater_mat_125 = ifelse(is.na(greater_mat_125), 0, greater_mat_125),
         greater_mat_less_125 = ifelse(is.na(greater_mat_less_125), 0, greater_mat_less_125),
         greater_50_less_mat = ifelse(is.na(greater_50_less_mat), 0, greater_50_less_mat)) %>% 
  pivot_longer(cols= c(greater_mat_125, greater_mat_less_125, greater_50_less_mat), names_to="Maturity", values_to="Abundance") %>% 
  mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
  mutate(Maturity = fct_recode(Maturity, "< 50 Length Maturity" = "less_50", ">50 Maturity but < Maturity"="greater_50_less_mat",
                               "> Length Maturity but < 1.25x Maturity" = "greater_mat_less_125",
                               "> 1.25x Maturity"="greater_mat_125")) %>% 
  glimpse()

dat.preds <- dat %>% 
  dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") %>% 
  distinct()

dat.species <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample)) %>% 
  mutate(Mat.Species = paste0(Maturity, sep="_", scientific))


unique.vars <- unique(as.character(dat.species$Mat.Species))

resp.vars <- character()
for(i in 1:length(unique.vars)){
  temp.dat <- dat.species[which(dat.species$Mat.Species == unique.vars[i]), ]
  if(length(which(temp.dat$Abundance == 0)) / nrow(temp.dat) < 0.95){ # Change here
    resp.vars <- c(resp.vars, unique.vars[i])}
}
resp.vars   

# Run the full subset model selection
savedir <- out_dir
use.dat <- as.data.frame(dat.species) 
factor.vars <- c("method") # Method as a factors with 2 levels
out.all     <- list()
var.imp     <- list()

# Loop through the FSS function for each response variable----
for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  use.dat <- as.data.frame(dat.species[which(dat.species$Mat.Species == resp.vars[i]), ])
  Model1  <- gam(Abundance ~ s(depth, k=3, bs='cr') + method,
                 family = tw(),  data = use.dat)
  
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

# Save model fits and importance scores
names(out.all) <- resp.vars
names(var.imp) <- resp.vars
all.mod.fits   <- do.call("rbind",out.all)
all.var.imp    <- do.call("rbind",var.imp)
write.csv(all.mod.fits[ , -2], file = paste(savedir, paste(name, "mat_groups", "by_species", "all.mod.fits.csv", sep = "_"), sep = "/"))
write.csv(all.var.imp, file = paste(savedir, paste(name, "mat_groups", "by_species", "all.var.imp.csv", sep = "_"), sep = "/"))

#* Check top models make sense, particularly for 95% 0s ####

use.dat <- dat.species %>% 
  filter(Mat.Species %in% c("> 1.25x Maturity_Labridae Choerodon rubescens"))

Best1 <- gam(Abundance~s(tpi, k = 3, bs = "cr") + depth + method,
             family=tw(), data=use.dat)

summary(Best1)
gam.check(Best1)

use.dat <- dat.species %>% 
  filter(Mat.Species %in% c("> 1.25x Maturity_Lethrinidae Lethrinus miniatus"))

Best1 <- gam(Abundance~s(detrended, k = 3, bs = "cr") + s(roughness, k = 3, bs = "cr") + method,
             family=tw(), data=use.dat)

summary(Best1)
gam.check(Best1)

use.dat <- dat.species %>% 
  filter(Mat.Species %in% c("> Length Maturity but < 1.25x Maturity_Lethrinidae Lethrinus miniatus"))

Best1 <- gam(Abundance~ s(biog, k = 3, bs = "cr") + s(detrended, k = 3, bs = "cr") + s(macroalgae, k = 3, bs = "cr") + method,
             family=tw(), data=use.dat)

summary(Best1)
gam.check(Best1)

use.dat <- dat.species%>% 
  filter(Mat.Species %in% c(">50 Maturity but < Maturity_Lethrinidae Lethrinus miniatus"))

Best1 <- gam(Abundance~ s(macroalgae, k = 3, bs = "cr") + method,
             family=tw(), data=use.dat)

summary(Best1)
gam.check(Best1)

use.dat <- dat.species %>% 
  filter(Mat.Species %in% c("> 1.25x Maturity_Sparidae Chrysophrys auratus"))

Best1 <- gam(Abundance~ s(macroalgae, k = 3, bs = "cr") + method,
             family=tw(), data=use.dat)

summary(Best1)
gam.check(Best1)

use.dat <- dat.species%>% 
  filter(Mat.Species %in% c("> Length Maturity but < 1.25x Maturity_Sparidae Chrysophrys auratus"))

Best1 <- gam(Abundance~ s(biog, k = 3, bs = "cr") + s(macroalgae, k = 3, bs = "cr") + method,
             family=tw(), data=use.dat)

summary(Best1)
gam.check(Best1)


#### CHECKING FOR WHICH METHOD HAS MORE ZEROES ####

unique.vars <- unique(as.character(dat$Maturity))
methods <- c("BRUV", "BOSS")

Method.0 <- array(0, dim=c(3,3))
Method.0[,1] <- unique.vars
colnames(Method.0) <- c("Method", "BRUV", "BOSS")
  
for(M in 1:2){
  
  dat.temp <- dat %>% 
    filter(method %in% c(methods[M]))
  
  for(i in 1:length(unique.vars)){
    temp.dat <- dat.temp %>% 
      filter(Maturity %in% c(unique.vars[i]))
    temp.0 <- (length(which(temp.dat$Abundance == 0)) / nrow(temp.dat))*100
    
    Method.0[i,M+1] <- temp.0
  }
}
  
  

