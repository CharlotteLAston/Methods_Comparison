rm(list=ls())

library(dplyr)
library(tidyr)
# library(gridExtra)
# library(grid)
library(GlobalArchive)
library(stringr)
library(ggplot2)
library(gamm4)
# library(ggmap)
# library(rgdal)
# library(raster)
# library(png)
# library(cowplot)
# library(tibble)
library(purrr)

# set theme
# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend.background = element_rect(fill="white"),
    legend.background = element_blank(),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=10),
    legend.title = element_blank(),
    legend.position = c(0.2, 0.8),
    text=element_text(size=10),
    strip.text.y = element_text(size = 10,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=10),
    axis.title.y=element_text(vjust=0.6, angle=90, size=10),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

## Set working directory----
working.dir <- getwd()
setwd(working.dir)
#OR Set manually once

# Load and join the raw data ----
# Maxn
dat.maxn <- readRDS("data/Tidy/dat.maxn.rds")%>%
  dplyr::filter(location%in%"NPZ6")%>%
  dplyr::rename(number = maxn)%>%
  glimpse()

# Lengths
dat.length <- readRDS("data/Tidy/dat.length.rds") %>%
  dplyr::filter(location%in%"NPZ6") %>%
  glimpse()

# Stack both
dat <- bind_rows(dat.maxn,dat.length)

# Model outputs from fss-gam script - filtered to the top model for each taxa
# This needs a bit of tidying in the original export
topmodels <- readRDS("data/Tidy/all.mod.fits.RDS") %>%
  dplyr::mutate(rownames = row.names(.)) %>% # Need to fix this in the original export 
  dplyr::mutate(taxa = c("total.abundance", "total.abundance",
                         "total.abundance", "species.richness")) %>% # Need to fix this in the original exports - row names not correct
  tidyr::separate(col = modname, into = c("term1", "term2", "term3"), sep = "\\+", remove = F) %>%
  dplyr::mutate(no.terms = rowSums(!is.na(dplyr::select(., starts_with("term"))))) %>%
  dplyr::group_by(taxa) %>% # Per response
  dplyr::arrange(AICc, no.terms) %>% # Sort it by AICc and the number of terms
  slice(1) %>% # Slice off the top model
  ungroup() %>% # Just to be safe
  glimpse()
  
topmodels$formula <- gsub("\"","\'", topmodels$formula) # Swap the " for ' - otherwise errors when pasting in
topmodels$formula <- substr(topmodels$formula, 1, nchar(topmodels$formula) - 9) # Remove the fixed effect for method

resp.vars <- unique(topmodels$taxa) # Number of taxa from top-models fss-gam export
factor.vars <- "" # Turn this on if you have factor variables in the model

# Make the dataframes to predict data
for (i in 1:length(resp.vars)) { # Start of loop to format the data - loop through each taxa
  
  # Create table of model terms only
  modterms <- topmodels %>% 
    pivot_longer(cols = starts_with("term"), names_to = "NA", values_to = "terms") %>%
    dplyr::select(taxa, terms) %>%
    dplyr::filter(!is.na(terms)) %>%
    group_by(taxa) %>%
    dplyr::mutate(term.no = row_number()) %>%
    ungroup() %>%
    glimpse()
  
  # Filtered model information
  modterms_f <- modterms %>% 
    dplyr::filter(taxa %in% resp.vars[i]) %>%
    glimpse()
  
  # List of covariates in the model
  covars <- unique(modterms_f$terms) # Just the covariates of that model - same as var.names?
  
  # Raw covariate data
  covar_f <- dat %>%
    dplyr::select(c(all_of(covars))) # Only columns in the original data that match the covariates of the model
  
  # Raw data per taxa
  dat_f <- dat %>%
    dplyr::filter(scientific %in% resp.vars[i])
  
  # Set the model - from the top model dataframe loaded from fssgam output
  mod <- gam(eval(parse(text=(paste0("number~", topmodels$formula[i])))), family = "tw",data = dat_f)
  
  n.vars <- ncol(covar_f) # Number of variables from the subset data
  var.names <- colnames(covar_f) # Names of the covariates
  
  for (i in 1:n.vars) { # Start of the loop to predict the data
    temp1 <- data.frame(sequ = seq(min(covar_f[i]), max(covar_f[i]), length = 20),
                        var = var.names[i])
    
    means <- tibble::rownames_to_column(data.frame(means = colMeans(covar_f)), "var") %>%
      spread(var, means)
    means <- do.call("rbind", replicate(20, means, simplify = FALSE))
    
    temp1 <- temp1 %>%
      cbind(means) %>%
      dplyr::select(-all_of(as.character(unique(temp1$var))))
    
    names(temp1)[names(temp1) == names(temp1[1])] <- as.character(unique(temp1[2]))
    
    fits <- predict.gam(mod, newdata = temp1, type = 'response', se.fit=T)
    
    temp2 = temp1 %>% 
      data.frame (fits) %>%
      group_by_at(1) %>% # Only change here
      dplyr::summarise(number = mean(fit),se.fit = mean(se.fit)) %>%
      dplyr::mutate(se.pos = number + se.fit,
                    se.neg = number - se.fit) %>%
      ungroup()
    assign(as.character(paste(unique(modterms_f$taxa), names(temp1[1]), sep = ".")), temp2)
    
    # Make the plots now
    if (i == 1) {
      # Make the plot including title
      ggmod <- ggplot() +
        ylab("")+
        xlab(as.character(names(temp1[1])))+
        geom_point(data = dat_f,aes_string(x = var.names[i], y = "number"),
                   alpha = 0.2, size = 1,show.legend = F) +
        geom_line(data = eval(parse(text = paste(unique(modterms_f$taxa), var.names[i], sep = "."))),
                  aes_string(x = var.names[i],y = "number"),alpha = 0.5) +
        geom_line(data = eval(parse(text = paste(unique(modterms_f$taxa), var.names[i], sep = "."))),
                  aes_string(x = var.names[i],y = "se.pos"),
                  linetype = "dashed",alpha = 0.5) +
        geom_line(data = eval(parse(text = paste(unique(modterms_f$taxa), var.names[i], sep = "."))),
                  aes_string(x = var.names[i],y = "se.neg"),
                  linetype = "dashed",alpha = 0.5) +
        theme_classic() +
        Theme1 +
        ggtitle(as.character(unique(modterms_f$taxa))) +
        theme(plot.title = element_text(hjust = 0))
      assign(paste("ggmod", as.character(paste(unique(modterms_f$taxa), names(temp1[1]), sep = ".")), sep = "."), ggmod)
    }
    else {
      # Make the plot without title - for second and third terms in the model
      ggmod <- ggplot() +
        ylab("")+
        xlab(as.character(names(temp1[1])))+
        geom_point(data = dat_f,aes_string(x = var.names[i], y = "number"),
                   alpha = 0.2, size = 1,show.legend = F) +
        geom_line(data = eval(parse(text = paste(unique(modterms_f$taxa), var.names[i], sep = "."))),
                  aes_string(x = var.names[i],y = "number"),alpha = 0.5) +
        geom_line(data = eval(parse(text = paste(unique(modterms_f$taxa), var.names[i], sep = "."))),
                  aes_string(x = var.names[i],y = "se.pos"),
                  linetype = "dashed",alpha = 0.5) +
        geom_line(data = eval(parse(text = paste(unique(modterms_f$taxa), var.names[i], sep = "."))),
                  aes_string(x = var.names[i],y = "se.neg"),
                  linetype = "dashed",alpha = 0.5) +
        theme_classic() +
        Theme1
      assign(paste("ggmod", as.character(paste(unique(modterms_f$taxa), names(temp1[1]), sep = ".")), sep = "."), ggmod)
    }
  } # End of the Predict loop

} # End of the loop to format data 
