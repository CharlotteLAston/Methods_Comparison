rm(list=ls()) # Clear memory

## Load Libraries ----
# To connect to GlobalArchive
library(devtools)
install_github("UWAMEGFisheries/GlobalArchive") #to check for updates
library(GlobalArchive)
# To connect to GitHub
library(RCurl)
library(R.utils)
# To tidy data
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
# to connect to googlesheets
library(googlesheets4)

## Set your working directory ----
working.dir<- dirname(rstudioapi::getActiveDocumentContext()$path)

## Save these directory names to use later----
full.data <- paste(working.dir,"full_data",sep="/")
tidy.dir <- paste(working.dir,"tidy_data",sep="/")
download.dir <- paste(working.dir,"EM Export",sep="/")
staging.dir <-  paste(working.dir,"staging",sep="/")

setwd(working.dir)

# ABROLHOS ####

### Set Study Name ----
# Change this to suit your study name. This will also be the prefix on your final saved files.
study<-"2021-05_Abrolhos_BOSS" 

### Metadata ----
metadata <- ga.list.files("_metadata.csv") %>% # list all files ending in "_metadata.csv"
  purrr::map_df(~ga.read.files_em.csv(.)) %>% # combine into dataframe
  dplyr::select(campaignid,sample,latitude,longitude,date,time,location,status,site,depth,observer,successful.count,successful.length) %>% # This line ONLY keep the 15 columns listed. Remove or turn this line off to keep all columns (Turn off with a # at the front).
  mutate(campaignid = ifelse(str_detect(campaignid, "_metadata.csv"), str_replace(campaignid, "_metadata.csv", ""), campaignid)) %>% 
  dplyr::filter(campaignid%in%c("2021-05_Abrolhos_BOSS"))%>%
  glimpse()

unique(metadata$campaignid) # check the number of campaigns in metadata, and the campaign name

setwd(staging.dir)
write.csv(metadata,paste(study,"metadata.csv",sep="_"),row.names = FALSE)

### Combine Points and Count files into maxn ----
points.files <-ga.list.files("_Points.txt") # list all files ending in "Lengths.txt"
points.files$lines<-sapply(points.files,countLines) # Count lines in files (to avoid empty files breaking the script)
points<-as.data.frame(points.files)%>%
  dplyr::mutate(campaign=row.names(.))%>%
  filter(lines>1)%>% # filter out all empty text files
  dplyr::select(campaign)%>%
  as_vector(.)%>% # remove all empty files
  purrr::map_df(~ga.read.files_em.txt(.))%>%
  dplyr::filter(campaignid%in%c("2021-05_Abrolhos_BOSS"))

maxn<-points%>%
  dplyr::select(-c(sample)) %>%
  dplyr::rename(sample = period) %>%
  dplyr::group_by(campaignid,sample,filename,periodtime,frame,family,genus,species)%>%
  dplyr::mutate(number=as.numeric(number))%>%
  dplyr::summarise(maxn=sum(number))%>%
  dplyr::group_by(campaignid,sample,family,genus,species)%>%
  dplyr::slice(which.max(maxn))%>%
  dplyr::ungroup()%>%
  dplyr::filter(!is.na(maxn))%>%
  dplyr::select(-frame)%>%
  tidyr::replace_na(list(maxn=0))%>%
  dplyr::mutate(maxn=as.numeric(maxn))%>%
  dplyr::filter(maxn>0)%>%
  dplyr::inner_join(metadata)%>%
  dplyr::filter(successful.count=="Y")%>%
  dplyr::filter(maxn>0)%>%
  glimpse()

### Save MaxN file ----
setwd(staging.dir)
write.csv(maxn,paste(study,"maxn.csv",sep="_"),row.names = FALSE)

unique(maxn$sample)


### Combine Length, Lengths and 3D point files into length3dpoints----
length3dpoints<-ga.create.em.length3dpoints()%>%
  dplyr::select(-c(time,comment))%>% # take time out as there is also a time column in the metadata
  dplyr::select(-c(sample)) %>%
  dplyr::rename(sample=period)%>%
  dplyr::inner_join(metadata)%>%
  dplyr::filter(successful.length=="Y")%>%
  glimpse()

### Save length files ----
setwd(staging.dir)
write.csv(length3dpoints,paste(study,"length3dpoints.csv",sep="_"),row.names = FALSE)

# NINGALOO PT CLOATES ####
### Set Study Name ----
# Change this to suit your study name. This will also be the prefix on your final saved files.
study<-"Ningaloo_PtCloates_BOSS" 

### Metadata ----
setwd(download.dir)
metadata <- ga.list.files("_Metadata.csv") %>% # list all files ending in "_metadata.csv"
  purrr::map_df(~ga.read.files_em.csv(.)) %>% # combine into dataframe
  mutate(sample = ifelse(is.na(sample), sample...2, sample)) %>% 
  dplyr::select(campaignid,sample,latitude,longitude,date,time,location,status,site,depth,successful.count,successful.length) %>% # This line ONLY keep the 15 columns listed. Remove or turn this line off to keep all columns (Turn off with a # at the front).
  dplyr::mutate(campaignid = str_replace(campaignid, "_metadata.csv", "")) %>% 
  dplyr::filter(campaignid%in%c("2022-05_PtCloates_Naked-BOSS","2022-05_PtCloates_Squid-BOSS","2021-08_PtCloates_Flasher-BOSS", "2021-08_PtCloates_Squid-BOSS")) %>%
  dplyr::mutate(campaignid = str_replace(campaignid, "Squid-", "")) %>%
  dplyr::mutate(campaignid = str_replace(campaignid, "Naked-", "")) %>%
  dplyr::mutate(campaignid = str_replace(campaignid, "Flasher-", "")) %>%
  filter(successful.count=="Yes"|successful.count=="Y") %>% 
  filter(!location %in% c("Deep Yardie"))

unique(metadata$campaignid) # check the number of campaigns in metadata, and the campaign name

setwd(staging.dir)
write.csv(metadata,paste(study,"metadata.csv",sep="_"),row.names = FALSE)

### Combine Points and Count files into maxn ----
points.files <-ga.list.files("_Points.txt") # list all files ending in "Lengths.txt"
points.files$lines<-sapply(points.files,countLines) # Count lines in files (to avoid empty files breaking the script)
points<-as.data.frame(points.files)%>%
  dplyr::mutate(campaign=row.names(.))%>%
  filter(lines>1)%>% # filter out all empty text files
  dplyr::select(campaign)%>%
  as_vector(.)%>% # remove all empty files
  purrr::map_df(~ga.read.files_em.txt(.))%>%
  dplyr::mutate(campaignid = str_replace(campaignid, "Squid-", "")) %>%
  dplyr::mutate(campaignid = str_replace(campaignid, "Naked-", "")) %>%
  dplyr::mutate(campaignid = str_replace(campaignid, "Flasher-", "")) %>%
  dplyr::filter(campaignid%in%c("2022-05_PtCloates_BOSS", "2021-08_PtCloates_BOSS"))

maxn<-points%>%
  dplyr::select(-c(sample)) %>%
  dplyr::rename(sample = period) %>%
  dplyr::group_by(campaignid,sample,filename,periodtime,frame,family,genus,species) %>%
  dplyr::mutate(number=as.numeric(number)) %>%
  dplyr::summarise(maxn=sum(number)) %>%
  dplyr::group_by(campaignid,sample,family,genus,species) %>%
  dplyr::slice(which.max(maxn)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(maxn)) %>%
  dplyr::select(-frame) %>%
  tidyr::replace_na(list(maxn=0)) %>%
  dplyr::mutate(maxn=as.numeric(maxn)) %>%
  dplyr::filter(maxn>0) %>%
  dplyr::inner_join(metadata) %>%
  dplyr::filter(successful.count=="Yes"|successful.count=="Y") %>%
  dplyr::filter(maxn>0)%>%
  glimpse()

test <- metadata %>%
  anti_join(maxn)%>%
  glimpse()

### Save MaxN file ----
setwd(staging.dir)
write.csv(maxn,paste(study,"maxn.csv",sep="_"),row.names = FALSE)


### Combine Length, Lengths and 3D point files into length3dpoints----
length3dpoints<-ga.create.em.length3dpoints()%>%
  dplyr::select(-c(time,comment))%>% # take time out as there is also a time column in the metadata
  dplyr::select(-c(sample)) %>%
  dplyr::rename(sample=period)%>%
  dplyr::mutate(campaignid = str_replace(campaignid, "Squid-", "")) %>%
  dplyr::mutate(campaignid = str_replace(campaignid, "Naked-", "")) %>%
  dplyr::mutate(campaignid = str_replace(campaignid, "Flasher-", "")) %>%
  dplyr::inner_join(metadata)%>%
  dplyr::filter(successful.length=="Yes")%>%
  glimpse()

### Save length files ----
setwd(staging.dir)
write.csv(length3dpoints,paste(study,"length3dpoints.csv",sep="_"),row.names = FALSE)

