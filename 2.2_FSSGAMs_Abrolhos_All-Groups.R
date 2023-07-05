###
# Project: Parks - Abrolhos
# Data:    BOSS & BRUV fish, habitat
# Task:    Modelling fish lengths w/ FSSGAM
# author:  Claude, Brooke, Kingsley
# date:    Nov-Dec 2021
##

# Part 1-FSS modeling
rm(list=ls())

## librarys
# detach("package:plyr", unload=TRUE)#will error - don't worry
library(tidyr)
library(dplyr)
options(dplyr.width = Inf) #enables head() to display all coloums
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(gplots)
library(RColorBrewer)
library(doParallel) #this can removed?
library(doSNOW)
library(gamm4)
library(RCurl) #needed to download data from GitHub
library(FSSgam)
library(GlobalArchive)
library(ggplot2)

## set study name
study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

#### SET DIRECTORIES AND READ IN DATA ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")

setwd(data_dir)
dat <- readRDS("dat.length.rds")%>%
  glimpse()

abrolhos.indicator <- c("Sparidae Chrysophrys auratus", "Labridae Choerodon rubescens", "Lethrinidae Lethrinus miniatus",
                        "Glaucosomatidae Glaucosoma herbraicum", "Berycidae Centroberyx gerrardi")



#### ABUNDANCE INDICATOR FISH >MATURE LENGTH ####
dat1 <- dat %>% 
  mutate(length = ifelse(is.na(length), 0, length)) %>% 
  mutate(number = ifelse(scientific %in% abrolhos.indicator, number, 0)) %>% 
  mutate(length.mat = ifelse(scientific %in% "Sparidae Chrysophrys auratus", 262,
                             ifelse(scientific %in% "Labridae Choerodon rubescens", 270,
                                    ifelse(scientific %in% "Lethrinidae Lethrinus miniatus", 361,
                                           ifelse(scientific %in% "Glaucosomatidae Glaucosoma herbraicum", 301,
                                                  ifelse(scientific %in% "Berycidae Centroberyx gerrardi", 250, 0)))))) %>% 
  mutate(number = ifelse(length >= length.mat, number, 0))%>% 
  group_by(sample) %>% 
  mutate(abundance = sum(number)) %>% 
  ungroup() %>% 
  mutate(sample = as.factor(sample)) %>% 
  distinct(sample, .keep_all=T) %>% 
  dplyr::select("sample", "method", "depth", "macroalgae", "biog", "mean.relief","tpi","roughness","detrended", "abundance")
  

# Set the predictors for modeling
pred.vars <- c("depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") 

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
unique.vars <- unique(as.character(dat$scientific))

resp.vars <- character()
for(i in 1:length(unique.vars)){
  temp.dat <- dat[which(dat$scientific == unique.vars[i]), ]
  if(length(which(temp.dat$maxn == 0)) / nrow(temp.dat) < 0.8){ # Change here
    resp.vars <- c(resp.vars, unique.vars[i])}
}
resp.vars   

# Run the full subset model selection
savedir <- "output/"
use.dat <- as.data.frame(dat) 
factor.vars <- c("status") # Status as a factors with 2 levels
out.all     <- list()
var.imp     <- list()

# Loop through the FSS function for each Taxa----
for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  use.dat <- as.data.frame(dat[which(dat$scientific == resp.vars[i]), ])
  Model1  <- gam(maxn ~ s(depth, k = 3, bs='cr'),
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
write.csv(all.mod.fits[ , -2], file = paste(savedir, paste(name, "all.mod.fits.csv", sep = "_"), sep = "/"))
write.csv(all.var.imp, file = paste(savedir, paste(name, "all.var.imp.csv", sep = "_"), sep = "/"))

##### ABUNDANCE INDICATOR SPECIES >MATURE LENGTH*2 ####
dat2 <- dat %>% 
  mutate(length = ifelse(is.na(length), 0, length)) %>% 
  mutate(number = ifelse(scientific %in% abrolhos.indicator, number, 0)) %>% 
  mutate(length.mat = ifelse(scientific %in% "Sparidae Chrysophrys auratus", 262,
                             ifelse(scientific %in% "Labridae Choerodon rubescens", 270,
                                    ifelse(scientific %in% "Lethrinidae Lethrinus miniatus", 361,
                                           ifelse(scientific %in% "Glaucosomatidae Glaucosoma herbraicum", 301,
                                                  ifelse(scientific %in% "Berycidae Centroberyx gerrardi", 250, 0)))))) %>% 
  filter(length>=(length.mat*2)) %>% 
  group_by(sample) %>% 
  mutate(abundance = sum(number)) %>% 
  ungroup() %>% 
  distinct(sample, .keep_all=T) %>% 
  dplyr::select("sample", "method", "depth", "macroalgae", "biog", "mean.relief","tpi","roughness","detrended", "abundance")



# # Re-set the predictors for modeling
pred.vars = c("depth", 
              "macroalgae", 
              "biog", 
              "mean.relief",
              "tpi",
              "roughness",
              "detrended")



#* Run the full subset model selection ####
savedir <- "fssgam_output"

use.dat=as.data.frame(dat2) 
str(use.dat)

name<- paste(study,"indicator_length_mat_x2",sep="_")

factor.vars=c("method")# Method as a Factor with two levels
out.all=list()
var.imp=list()

### Can't use random effect because there aren't enough samples 
### Can't do for mat *2 because we didn't see any of them on the BRUVs for indicator species 
use.dat <- as.data.frame(dat2) 
use.dat$method <- as.factor(use.dat$method)
use.dat$sample <- as.factor(use.dat$sample)
Model1=gam(abundance~s(depth,k=3,bs='cr'),
           family=tw(),  data=use.dat)

model.set=generate.model.set(use.dat=use.dat,
                             test.fit=Model1,
                             pred.vars.cont=pred.vars,
                             pred.vars.fact=factor.vars,
                             factor.smooth.interactions = NA,
                             # smooth.smooth.interactions = c("depth", "biog"),
                             k=3
                             #null.terms="s(sample,bs='re')"
)
out.list=fit.model.set(model.set,
                       max.models=600,
                       parallel=T)
names(out.list)

out.list$failed.models # examine the list of failed models
mod.table=out.list$mod.data.out  # look at the model selection table
mod.table=mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
out.i=mod.table[which(mod.table$delta.AICc<=5),]
out.all=c(out.all,list(out.i))
# var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2

# plot the best models
setwd(fig_dir)
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[m])
  
  png(file=paste(name,"mod_fits.png",sep="_"))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model=out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    #mtext(side=2,text=resp.vars[i],outer=F)
  }  
  dev.off()
}


#* Model fits and importance ####
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
setwd(out_dir)
write.csv(all.mod.fits[ , -5], file = paste(name, "all.mod.fits.csv", sep = "_"))
write.csv(all.var.imp, file = paste(name, "all.var.imp.csv", sep = "_"))

#### OVERALL ABUNDANCE OF INDICATOR SPECIES ####
dat3 <- dat %>% 
  mutate(number = ifelse(scientific %in% abrolhos.indicator, number, 0)) %>% 
  mutate(length.mat = ifelse(scientific %in% "Sparidae Chrysophrys auratus", 262,
                             ifelse(scientific %in% "Labridae Choerodon rubescens", 270,
                                    ifelse(scientific %in% "Lethrinidae Lethrinus miniatus", 361,
                                           ifelse(scientific %in% "Glaucosomatidae Glaucosoma herbraicum", 301,
                                                  ifelse(scientific %in% "Berycidae Centroberyx gerrardi", 250, 0)))))) %>% 
  group_by(sample) %>% 
  mutate(abundance = sum(number)) %>% 
  ungroup() %>% 
  distinct(sample, .keep_all=T) %>% 
  dplyr::select("sample", "method", "depth", "macroalgae", "biog", "mean.relief","tpi","roughness","detrended", "abundance")



# # Re-set the predictors for modeling
pred.vars = c("depth", 
              "macroalgae", 
              "biog", 
              "mean.relief",
              "tpi",
              "roughness",
              "detrended")


# changed to 90% - smaller than legal size included

#* Run the full subset model selection ####
savedir <- "fssgam_output"

use.dat=as.data.frame(dat3) 
str(use.dat)

name<- paste(study,"indicator_abundance",sep="_")

factor.vars=c("method")# Method as a Factor with two levels
out.all=list()
var.imp=list()

### Can't use random effect because there aren't enough samples 
### Can't do for mat *2 because we didn't see any of them on the BRUVs for indicator species 
use.dat$method <- as.factor(use.dat$method)
use.dat$sample <- as.factor(use.dat$sample)
Model1=gam(abundance~s(depth,k=3,bs='cr'),
           family=tw(),  data=use.dat)

model.set=generate.model.set(use.dat=use.dat,
                             test.fit=Model1,
                             pred.vars.cont=pred.vars,
                             pred.vars.fact=factor.vars,
                             factor.smooth.interactions = NA,
                             # smooth.smooth.interactions = c("depth", "biog"),
                             k=3
                             #null.terms="s(sample,bs='re')"
)
out.list=fit.model.set(model.set,
                       max.models=600,
                       parallel=T)
names(out.list)

out.list$failed.models # examine the list of failed models
mod.table=out.list$mod.data.out  # look at the model selection table
mod.table=mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
out.i=mod.table[which(mod.table$delta.AICc<=5),]
out.all=c(out.all,list(out.i))
# var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2

# plot the best models
setwd(fig_dir)
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[m])
  
  png(file=paste(name,"mod_fits.png",sep="_"))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model=out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    #mtext(side=2,text=resp.vars[i],outer=F)
  }  
  dev.off()
}


#* Model fits and importance ####
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
setwd(out_dir)
write.csv(all.mod.fits[ , -5], file = paste(name, "all.mod.fits.csv", sep = "_"))
write.csv(all.var.imp, file = paste(name, "all.var.imp.csv", sep = "_"))

##### ABUNDANCE INDICATOR SPECIES <MATURE LENGTH ####

dat4 <- dat %>% 
  mutate(length = ifelse(is.na(length), 0, length)) %>% 
  mutate(number = ifelse(scientific %in% abrolhos.indicator, number, 0)) %>% 
  mutate(length.mat = ifelse(scientific %in% "Sparidae Chrysophrys auratus", 262,
                             ifelse(scientific %in% "Labridae Choerodon rubescens", 270,
                                    ifelse(scientific %in% "Lethrinidae Lethrinus miniatus", 361,
                                           ifelse(scientific %in% "Glaucosomatidae Glaucosoma herbraicum", 301,
                                                  ifelse(scientific %in% "Berycidae Centroberyx gerrardi", 250, 0)))))) %>% 
  mutate(number = ifelse(length > length.mat, number, 0))%>% 
  group_by(sample) %>% 
  mutate(abundance = sum(number)) %>% 
  ungroup() %>% 
  mutate(sample = as.factor(sample)) %>% 
  distinct(sample, .keep_all=T) %>% 
  dplyr::select("sample", "method", "depth", "macroalgae", "biog", "mean.relief", "sd.relief","tpi","roughness","detrended", "abundance")



# Re-set the predictors for modeling
pred.vars = c("depth", 
              "macroalgae", 
              "biog", 
              "mean.relief",
              "sd.relief",
              "tpi",
              "roughness",
              "detrended")


# changed to 90% - smaller than legal size included

# Run the full subset model selection
savedir <- "fssgam_output"

use.dat=as.data.frame(dat4) 
str(use.dat)

name<- paste(study,"indicator_length_less_mat",sep="_")

factor.vars=c("method")# Method as a Factor with two levels
out.all=list()
var.imp=list()

## Can't use random effect because there aren't enough samples
use.dat$method <- as.factor(use.dat$method)
use.dat$sample <- as.factor(use.dat$sample)
Model1=gam(abundance~s(depth,k=3,bs='cr'),
           family=tw(),  data=use.dat)

model.set=generate.model.set(use.dat=use.dat,
                             test.fit=Model1,
                             pred.vars.cont=pred.vars,
                             pred.vars.fact=factor.vars,
                             factor.smooth.interactions = NA,
                             # smooth.smooth.interactions = c("depth", "biog"),
                             k=3
                             #null.terms="s(sample,bs='re')"
)
out.list=fit.model.set(model.set,
                       max.models=600,
                       parallel=T)
names(out.list)

out.list$failed.models # examine the list of failed models
mod.table=out.list$mod.data.out  # look at the model selection table
mod.table=mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
out.i=mod.table[which(mod.table$delta.AICc<=5),]
out.all=c(out.all,list(out.i))
# var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2

# plot the best models
setwd(fig_dir)
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[m])
  
  png(file=paste(name,"mod_fits.png",sep="_"))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model=out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    #mtext(side=2,text=resp.vars[i],outer=F)
  }  
  dev.off()
}


# Model fits and importance
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
setwd(out_dir)
write.csv(all.mod.fits[ , -5], file = paste(name, "all.mod.fits.csv", sep = "_"))
write.csv(all.var.imp, file = paste(name, "all.var.imp.csv", sep = "_"))


