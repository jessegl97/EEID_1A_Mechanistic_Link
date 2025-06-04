#Data Cleaning Script EEID21 Antibody Analysis
#Requires: EEID2021_master_data_20220503.csv, EEID_1a_ELISA.csv
#Produces: reinfection_response.csv, reinfection_response.rds"

master <- read.csv("EEID2021_master_data_20220503.csv")
ab <- read.csv("EEID_1a_ELISA.csv")
m.ab <- left_join(master, ab, by = "bird_ID", keep=TRUE)

m.ab[m.ab == "f "] <- "f" #remove spaces

m.ab <- m.ab %>%
  dplyr::select(- date.y, -bird_ID.y, -band_number.y, -dpi, -(X:X.13)) #remove extra columns

m.ab <- m.ab %>%
  filter(experiment_location == "vt") #include only birds at VT

#rename columns
names(m.ab)[names(m.ab) == "date.x"] <- "date"
names(m.ab)[names(m.ab) == "band_number.x"] <- "band_number"
names(m.ab)[names(m.ab) == "bird_ID.x"] <- "bird_ID"
names(m.ab)[names(m.ab) == "Avg.OD"] <- "elisa_od"
names(m.ab)[names(m.ab) == "CV"] <- "elisa_cv"
names(m.ab)[names(m.ab) == "dppi"] <- "dpi"

m.ab$sex <- as.factor(m.ab$sex)
m.ab$sex <- factor(m.ab$sex, levels = c("f", "m"), labels = c("Female", "Male"))
m.ab$primary_treatment <- as.factor(m.ab$primary_treatment)
levels(m.ab$primary_treatment) <- c("High", "Low", "Sham")
m.ab$primary_treatment <- factor(m.ab$primary_treatment, levels = c("Sham", "Low", "High"))
m.ab$l_eye_score <- as.numeric(m.ab$l_eye_score) #Converts chr na to numeric NA on days where eyescore was not recorded
m.ab$r_eye_score <- as.numeric(m.ab$r_eye_score) #Converts chr na to numeric NA on days where eyescore was not recorded
#total eye score = sum of l and r eye score
m.ab$tes <- m.ab$l_eye_score + m.ab$r_eye_score
m.ab$band_number <- as.factor(m.ab$band_number)
m.ab$secondary_dose <- as.numeric(m.ab$secondary_dose)

#some sample days were split between two days. 
#One day experimental birds were measured (dpi -8 and dpi 14), the next day infected birds were measured (dpi -7, dpi15)
#this code makes a new column dpi.new that combines dpi 14 and 15
m.ab$dpi.new <- ifelse(m.ab$dpi %in% c(14, 15), "14", as.integer(m.ab$dpi))

#i then use dpi.new to replace the dpi column to combine dpi -7 and -8
m.ab$dpi <- ifelse(m.ab$dpi %in% c(-7, -8), "-8", as.integer(m.ab$dpi.new))

m.ab$dpi.new <- NULL #delete dpi.new column

m.ab$dpi <- as.numeric(m.ab$dpi)

#make baseline eyescore 0
m.ab <- m.ab %>%
  mutate(tes = ifelse(dpi == -8 & is.na(tes), 0, tes))

#install.packages("gtsummary")
library(gtsummary)

#Initial Sample Sizes
m.ab %>%
  filter(dpi == -8)%>%
  dplyr::select(primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=primary_treatment
  )%>%
  modify_header(
    label ~ "**All Birds Before Removal**"
  )


#Infection Thresholds
m.ab$threshold_cutoff = 50
m.ab$seropos_cutoff = 0.061
m.ab$pathology_cutoff = 0

#generate infection column: If quantity > 50, inf = 1, else, inf=0
m.ab <- m.ab %>%
   arrange(band_number, dpi) %>%

#   # Replace NA in quantity with 0 conservatively
   mutate(quantity_clean = coalesce(quantity, 0)) %>%
#   
#   # Flag individual days meeting quantity criteria
   mutate(day_inf = if_else(quantity_clean > 50, 1, 0)) %>%
#   
   group_by(band_number) %>%
#   
#   # Individual bird infection status based on any day meeting criteria
   mutate(
    inf = if_else(any(day_inf == 1), 1, 0),
    inf_pri = if_else(any(day_inf == 1 & dpi < 42), 1, 0),
    inf_sec = if_else(any(day_inf == 1 & dpi > 42), 1, 0)
  ) %>%

  ungroup()

m.ab %>%
  filter(dpi == -8) %>%
  dplyr::select(primary_treatment, secondary_dose, inf, inf_pri, inf_sec)%>%
  tbl_summary(
    by=primary_treatment
  )%>%
  modify_header(
    label ~ "**All Birds Before Removal**"
  )

#omit birds (n=6)
#omit 2505 from data set for analysis as it was positive at quarantine (n=1)
m.ab <- m.ab %>%
  filter(band_number != 2505)

m.ab %>%
  filter(dpi == -8) %>%
  dplyr::select(primary_treatment, secondary_dose, inf_pri, inf_sec)%>%
  tbl_summary(
    by=primary_treatment
  )%>%
  modify_header(
    label ~ "**All Birds After Seropos -8 Removed**"
  )

#omit 2452 (1.99E2 @ DPPI 7) and 2562 (1.68E2 DPPI -8) bc positive qpcr and they are shams (n=2)
m.ab <- m.ab %>% filter(!band_number %in% c(2452, 2562)) 

#omit sham birds with eye scores (0.5) (2419, 2434) (n=2)
m.ab <- m.ab %>% filter(!band_number %in% c(2419, 2434))

#omit sham bird with path load (2514) (67.6 @ DPPI 41) (n=1)
m.ab <- m.ab %>% filter(band_number != 2514)

m.ab %>%
  filter(dpi == -8) %>%
  dplyr::select(primary_treatment, secondary_dose, inf_pri, inf_sec)%>%
  tbl_summary(
    by=primary_treatment
  )%>%
  modify_header(
    label ~ "**Birds After Bad Shams and Seropos -8 Removed**"
  )

##Some birds did not fully recover after primary infection. These birds were included in primary analysis, but will be excluded
#from secondary analysis. I have identified these birds below, but *will not remove them until secondary analysis*.

#check birds that were not recovered by secondary infection: 2274, 2469, 2520, 2494 (2514 not recovered but omitted above)
not_recovered <- m.ab %>% 
  filter(band_number %in% c(2274, 2469, 2520, 2494))%>%
  dplyr::select(band_number, dpi, primary_treatment, secondary_dose, quantity, tes)

not_recovered %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**Birds Not Recovered DPI 41**"
  )

#Starting Sample Sizes
m.ab %>%
  filter(dpi == -8)%>%
  dplyr::select(primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=primary_treatment
  )%>%
  modify_header(
    label ~ "**Starting Sample Size**"
  )


#Sample Sizes Secondary Pathology (removing 4 unrecovered birds)
m.ab %>% filter(!band_number %in% c(unique(not_recovered$band_number))) %>%
  filter(dpi == -8)%>%
  dplyr::select(primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=primary_treatment
  )%>%
  modify_header(
    label ~ "**Sample Size Secondary**"
  )

head(m.ab)
tail(m.ab)

#save final dataframe
#saveRDS(m.ab, file ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Final Dataframes/reinfection_response.rds")
#write.csv(m.ab, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Final Dataframes/reinfection_response.csv", row.names=FALSE)
