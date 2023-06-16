library(tidyverse)
library(GlobalArchive)

boss <- read.delim("data/raw/TM Export/2021-05_Abrolhos_BOSS_Dot Point Measurements.txt",
                   header = T, skip = 4, stringsAsFactors = FALSE) %>% # read in the file
  dplyr::select(-c(FieldOfView, Fine)) %>%
  dplyr::mutate(CAAB_CODE = as.character(CAAB_CODE)) %>%
  dplyr::mutate(CAAB_CODE = ifelse(Broad %in% c("Unknown", "Open Water"), "00000001",
                                   ifelse(Type %in% "Ecklonia radiata", "54079009", 
                                          ifelse(Broad %in% "Invertebrate Complex", "20000000", CAAB_CODE)))) %>%
  dplyr::mutate(Broad = ifelse(CAAB_CODE %in% c("", NA, "NA"), "Unknown", Broad),
                CAAB_CODE = ifelse(CAAB_CODE %in% c("", NA, "NA"), "00000001", CAAB_CODE)) %>%
  glimpse() # preview

bruv <- read.delim("data/raw/TM Export/2021-05_Abrolhos_stereo-BRUVs_Forwards_Dot Point Measurements.txt",
                   header = T, skip = 4, stringsAsFactors = FALSE) %>% # read in the file
  dplyr::select(-c(FieldOfView, Fine)) %>%
  dplyr::mutate(CAAB_CODE = as.character(CAAB_CODE)) %>%
  dplyr::mutate(CAAB_CODE = ifelse(Broad %in% c("Unknown", "Open Water"), "00000001",
                                   ifelse(Type %in% "Ecklonia radiata", "54079009", 
                                          ifelse(Broad %in% "Invertebrate Complex", "20000000", CAAB_CODE)))) %>%
  dplyr::mutate(Broad = ifelse(CAAB_CODE %in% c("", NA, "NA"), "Unknown", Broad),
                CAAB_CODE = ifelse(CAAB_CODE %in% c("", NA, "NA"), "00000001", CAAB_CODE)) %>%
  glimpse() # preview

write.table(boss, "data/staging/CHECKEMTEST_2021-05_Abrolhos_BOSS_Dot Point Measurements.txt",
            sep = "\t", row.names = F, col.names = T, quote = F, na = "")

write.table(bruv, "data/staging/CHECKEMTEST_2021-05_Abrolhos_stereo-BRUVs_Forwards_Dot Point Measurements.txt",
            sep = "\t", row.names = F, col.names = T, quote = F, na = "")

boss.r <- read.delim("data/raw/TM Export/2021-05_Abrolhos_BOSS_Relief_Dot Point Measurements.txt",
                     header = T, skip = 4, stringsAsFactors = FALSE) %>%
  glimpse()

bruv.r
