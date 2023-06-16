
###
# Project: Parks - Ningaloo Post-Survey
# Data:    BRUVS, BOSS
# Task:    Getting data from labsheets
# author:  Kingsley Griffin
# date:    July 2021
##

library(googlesheets4)

# get sampling data

bruvd <- as.data.frame(read_sheet("https://docs.google.com/spreadsheets/d/1ZfW-XJKP0BmY2UXPNquTxnO5-iHnG9Kw3UuJbALCcrs/edit#gid=814068592",
                                  sheet = "2021-05_Abrolhos_stereo-BRUVs"))
bossd <- as.data.frame(read_sheet("https://docs.google.com/spreadsheets/d/1ZfW-XJKP0BmY2UXPNquTxnO5-iHnG9Kw3UuJbALCcrs/edit#gid=814068592",
                                  sheet = "2021-05_Abrolhos_BOSS"))

saveRDS(bossd, 'data/2105_abrolhos_boss.rds')
saveRDS(bruvd, 'data/2105_abrolhos_bruv.rds')
