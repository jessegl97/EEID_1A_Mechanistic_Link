#Data Cleaning Script EEID21 Antibody Analysis

master <- read.csv("EEID2021_master_data_20220503.csv")
ab <- read.csv("EEID_1a_ELISA.csv")
m.ab <- left_join(master, ab, by = "bird_ID", keep=TRUE)

m.ab[m.ab == "f "] <- "f"

m.ab <- m.ab %>%
  dplyr::select(- date.y, -bird_ID.y, -band_number.y, -dpi, -(X:X.13)) #remove extra columns

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
m.ab$l_eye_score <- as.numeric(m.ab$l_eye_score)
m.ab$r_eye_score <- as.numeric(m.ab$r_eye_score)
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

m.ab <- m.ab %>%
  filter(experiment_location == "vt") #include only birds at VT

#set color scheme primary treatment
pri_colors <- c("#619CFF", "#00BA38", "#F8766D")
#set color scheme secondary treatment
sec_colors <- c("#8C754B", "#77AB59", "#59A5D8", "#9F77D9", "#FA8072")
sec_default <- c( "#00BA38", "#619CFF", "#F564E3", "#BA00FF", "#F8766D")
def_colors <- c("#F8766D", "#00BA38", "#619CFF", "#F564E3", "#BA00FF", "#7A76A6")

#omit 2505 from data set for analysis as it was positive at quarantine
m.ab <- m.ab %>%
  filter(band_number != 2505)