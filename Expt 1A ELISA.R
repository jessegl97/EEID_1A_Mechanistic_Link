#################################
#EEID 2021 Expt 1A ELISA Analysis
#Jesse Garrett-Larsen
#Winter 2024

#Analysis of IgY antibody levels taken from plasma of birds infected with varying doses of MG (VA94)
#################################

rm(list=ls())

library(tidyverse)
library(dplyr)
#install.packages("glmmTMB")
library(glmmTMB)
library(effects)
library(AICcmodavg)

#read in data
ab_old <- read.csv("expt1a_elisa_raw.csv")
master <- read.csv("EEID2021_master_data_20220503.csv")
ab <- read.csv("EEID_1a_ELISA.csv")
#merge
#m.ab <- merge(master, ab, by = "bird_ID", all=TRUE)
m.ab <- left_join(master, ab, by = "bird_ID", keep=TRUE)
####Data Formatting####
#replace all "f " with f
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

m.ab$l_eye_score <- as.numeric(m.ab$l_eye_score)
m.ab$r_eye_score <- as.numeric(m.ab$r_eye_score)

#generate threshold cutoffs for seropositivity and path load
m.ab$threshold_cutoff = 50
m.ab$seropos_cutoff = 0.061


#generate infected column - if copy number > 50, infected.
#for each timepoint only - can become infected or uninfected at next timepoint
m.ab <- m.ab %>%
  mutate(infected = ifelse(quantity>threshold_cutoff, 1, 0))

#generate infection data; if path load > 50 copies = infected
#1 = bird was successfully infected at any point during the priming phase of the experiment
#coalesce function replaces any NA with 0 here. 
#To be conservative, I am always assuming there is no path load unless it's measured with qPCR
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(ever_infected_qpcr = ifelse(any(coalesce(dpi, 0) < 42 & coalesce(infected, 0) == 1), 1, 0)) %>%
  ungroup()

#generate infection data; if path load > 50 copies any time after secondary infection -> 1
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate (infected_sec = ifelse(any(coalesce(dpi, 0) > 42 & coalesce(infected, 0) == 1), 1, 0)) %>%
  ungroup()

#check
ggplot(m.ab%>% filter(dpi >42), aes(x=secondary_treatment, y=log10_quantity1))+
  geom_jitter(width=0.25)+
  geom_hline(yintercept = log10(50))+
  facet_wrap(~infected_sec ~ dpi, ncol=4)

#looks good
ggplot(m.ab %>%filter(dpi > 41), aes(x=dpi, y=log10_quantity1))+
  geom_jitter()#+
  facet_wrap(~infected_d56)

#total eye score = sum of l and r eye score
m.ab$tes <- m.ab$l_eye_score + m.ab$r_eye_score

ggplot(m.ab %>% filter(dpi < 42), aes(x=tes, y=tes, color=band_number))+
  geom_jitter(height=0.1)+
  facet_wrap(~ever_infected_qpcr)

#who has an eye score but no path load? 2451
m.ab%>%
  filter(dpi < 42 & ever_infected_qpcr == 0 & tes == 2.5)

m.ab%>%
  filter(band_number == 2451)%>%
  dplyr::select(dpi, quantity, l_eye_score, r_eye_score, elisa_od, primary_treatment, secondary_treatment)
#weird that there was no path load dpi 7

#generate seropositivity data; if elisa OD > 0.061 = seropositive
#1 = bird seroconverted at any point during the priming phase of the experiment
m.ab <- m.ab %>%
  mutate(seropos = ifelse(elisa_od>seropos_cutoff, 1, 0))

m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(ever_seropos = ifelse(any(coalesce(dpi, 0) < 42 & coalesce(seropos, 0) == 1), 1, 0)) %>%
  ungroup()

unique(m.ab$seropos)

#check seropositivity cutoff
ggplot(m.ab, aes(x=primary_treatment, y=elisa_od, color=primary_treatment))+
  geom_jitter()+
  geom_hline(yintercept=0.061)+
  facet_wrap(~seropos)


#check seropositivity against qPCR cutoff
ggplot(m.ab %>% filter(dpi < 42), aes(x=primary_treatment, y=elisa_od, color=primary_treatment))+
  geom_jitter()+
  geom_hline(yintercept=0.061)+
  facet_wrap(~ever_infected_qpcr)

serop_neg <- m.ab %>%
  group_by(band_number)%>%
  filter(elisa_od > 0.061 & quantity > 50  & dpi == 15)
serop_neg
unique(serop_neg$band_number)

#check for seropositive birds from quarantine from dataset
m.ab %>%
  filter(dpi == -8) %>%
  filter(elisa_od >= 0.061)%>%
  dplyr::select(band_number, elisa_od, primary_treatment, secondary_treatment, infected, seropos, dpi, ELISA.Run.Date)


#check for seropositive sham birds on experimental days that have not been inoculated
m.ab %>% 
  filter(dpi <= 41)%>%
  filter(primary_treatment == "sham")%>%
  filter(elisa_od>0.061)


#omit 2505 from data set for analysis as it was positive at quarantine
m.ab <- m.ab %>%
  filter(band_number != 2505)

#which birds are seropositive but qPCR negative?
m.ab%>%
  filter(ever_infected_qpcr == 0 & elisa_od > 0.061 & dpi < 42)%>%
  dplyr::select(band_number, dpi, tes, elisa_od, quantity, primary_treatment)

#2398, 2451, 2470, 2515
a<-m.ab%>%
  filter(band_number %in% c(2398, 2451, 2470, 2515))%>%
  dplyr::select(band_number, dpi, tes, elisa_od, quantity, primary_treatment)
a

#could just be very resistant birds who develop no eye score or path load

#some sample days were split between two days. 
#One day experimental birds were measured (dpi -8 and dpi 14), the next day infected birds were measured (dpi -7, dpi15)
#this code makes a new column dpi.new that combines dpi 14 and 15
m.ab$dpi.new <- ifelse(m.ab$dpi %in% c(14, 15), "14", as.integer(m.ab$dpi))
unique(m.ab$dpi.new)
#i then use dpi.new to replace the dpi column to combine dpi -7 and -8
m.ab$dpi <- ifelse(m.ab$dpi %in% c(-7, -8), "-8", as.integer(m.ab$dpi.new))
unique(m.ab$dpi)
m.ab$dpi.new <- NULL #delete dpi.new column

m.ab %>%
  filter(dpi == "NA")
m.ab$dpi <- as.numeric(m.ab$dpi)
any(is.na(m.ab$dpi))
m.ab <- m.ab %>%
  filter(experiment_location == "vt") #include only birds at VT


# Check for missing combinations of band_number and dpi
missing_combinations <- expand.grid(band_number = unique(m.ab$band_number), dpi = unique(m.ab$dpi))
missing_combinations <- anti_join(missing_combinations, m.ab, by = c("band_number", "dpi"))

# Display missing combinations
print(missing_combinations)
#2375 DPI 41 did not have enough plasma for ELISA

#probably could remove for this analysis

#general overview of antibody responses
ggplot(m.ab, aes(x=dpi, y=elisa_od, color=primary_treatment))+
  geom_point()+
  geom_hline(yintercept = 0.061)+
  geom_smooth(aes(x=dpi, y=elisa_od, group=band_number), alpha=0.5, size=0.5)+
  facet_wrap(~primary_treatment)
unique(m.ab$dpi)

#####Side Quest: Low dose seroconversion + outcomes
m.ab.low <- m.ab %>%
  filter(primary_treatment == "low")

ggplot(m.ab.low, aes(x=dpi, y=log10_quantity1, color=as.factor(ever_seropos)))+
  geom_point()+
  geom_smooth(aes(x=dpi, y=log10_quantity1, group=band_number), alpha=0.5, size=0.75, se=FALSE)+
  ylim(c(0,6))+
  facet_wrap(~secondary_dose, ncol=5)

ggplot(m.ab, aes(x=dpi, y=log10_quantity1, color=as.factor(ever_seropos)))+
  geom_point()+
  geom_line(aes(x=dpi, y=log10_quantity1, group=band_number), alpha=0.5, size=0.75)+
  ylim(c(0,6))+
  facet_wrap(~primary_dose~secondary_dose, ncol=5)

#were birds that seroconverted at low doses better protected upon secondary infection?
glm.a <- glm(infected_sec ~ ever_seropos, data=m.ab, family=binomial())
summary(glm.a)
plot(allEffects(glm.a))
simulateResiduals(glm.a, plot=T)
####Does inoculation dose predict antibodies?####
####A: Does inoculation dose predict antibody levels across primary infection?####

#data frame with only primary infection (no PID 56)
p.ab <- m.ab %>%
  filter(dpi <=41 & elisa_od >0) #elisa_od > 0 to filter out samples that need to be re-run

#are birds with higher elisa_od also the birds with higher elisa_cvs? Not really 
ggplot(p.ab, aes(x=elisa_cv, y=elisa_od))+
  geom_point()+
  geom_hline(yintercept = 0.061)

#make comparisons to control in model output
p.ab$primary_treatment <- factor(p.ab$primary_treatment, levels = c("sham", "low", "high"))
#gamma distribution because antibody data looks exponential
#gamma distribution; band_number as random effect b/c longitudinal

#confirm all elisa_ods were the same before infection
lm0 <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi <0), family=Gamma())
summary(lm0)
plot(allEffects(lm0))
#yes

#Does primary treatment predict antibody levels across primary infection?
#I use primary_treatment here (categorical) to look at differences between low and high dose
lm1 <- glmmTMB(elisa_od~primary_treatment + (1|band_number), data=p.ab, family=Gamma())
summary(lm1)
library(DHARMa)
l1r <- simulateResiduals(lm1)
plot(l1r)

plot(allEffects(lm1)) #note inverse link function flips y axis


#model comparison 
p1 <- glm(elisa_od~primary_treatment + (1|band_number), data=p.ab, family=Gamma())
p2<- glm(elisa_od~primary_treatment + sex + (1|band_number), data=p.ab, family=Gamma())
p3<- glm(elisa_od~primary_treatment * sex + (1|band_number), data=p.ab, family=Gamma())
p4 <- glm(elisa_od~1, + (1|band_number), data=p.ab, family=Gamma())
p5 <- glm(elisa_od~primary_treatment + mass + (1|band_number), data=p.ab, family=Gamma())
p6<- glm(elisa_od~primary_dose + (1|band_number), data=p.ab, family=Gamma())

#AICc
aictab(cand.set=list(p1, p2, p3, p4, p5, p6), modnames=c("p1","p2", "p3", "p4", "p5", "p6"))

#p1 best model
library(multcomp)

#hypothesis testing
# Specify the contrast for the comparison
contrast_matrix <- matrix(c(0, 1, -1), nrow = 1)

# Perform the contrast using glht
contrast_test <- glht(lm1, contrast_matrix)

# Summarize the results
summary(contrast_test)
#high primary dose had significantly higher antibody levels than low primary dose (p < 0.0001)

#antibodies across all days primary infection
ggplot(data=p.ab, aes(x=as.factor(primary_dose), y= elisa_od, shape=primary_treatment))+
  geom_jitter(width=0.1)+
  geom_hline(yintercept = 0.061, linetype="dashed")+
  stat_summary(aes(group=primary_treatment), fun=mean, geom="point", alpha=1, size=10, shape="-", color=
                 "red")+
  stat_summary(aes(group=primary_treatment, shape = primary_treatment), fun.y=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               geom= "errorbar", size=0.5, width=0.25)+
  labs(title = "ELISA OD Primary Infection", y="ELISA OD", x= "Primary Treatment", shape="Primary Treatment")



##### (1) Does inoculation dose predict antibody levels on day 14? ####

#glm b/c only using one dpi so no need for bird_id as a fixed effect
lm2 <- glm(elisa_od~primary_treatment, data=p.ab %>% filter(dpi == 14), family=Gamma())
hist(p.ab$elisa_od)
summary(lm2)


plot(allEffects(lm2))
#Primary treatment predicts antibody levels on day 14/15
hist(resid(lm2))
qqnorm(resid(lm2))
qqline(resid(lm2))
shapiro.test(resid(lm2))


#graph showing antibody levels on dpi 14/15 by primary_treatment
ggplot(data=p.ab %>% filter(dpi== 14), aes(x=as.factor(primary_dose), y= elisa_od, shape=primary_treatment))+
  geom_jitter(width=0.1)+
  geom_hline(yintercept = 0.061, linetype="dashed")+
  stat_summary(aes(group=primary_treatment), fun=mean, geom="point", alpha=1, size=10, shape="-")+
  stat_summary(aes(group=primary_treatment, shape = primary_treatment), fun.y=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               geom= "errorbar", size=0.5, width=0.25)+
  labs(title = "ELISA OD DPI 14", y="ELISA OD", x= "Primary Treatment", shape="Primary Treatment")


#model comparison
pr1 <- glm(elisa_od~primary_treatment, data=p.ab %>% filter(dpi == 14), family=Gamma())
pr2<- glm(elisa_od~primary_treatment + sex , data=p.ab %>% filter(dpi == 14), family=Gamma())
pr2.5<- glm(elisa_od~primary_treatment * sex , data=p.ab %>% filter(dpi == 14), family=Gamma())
pr3 <- glm(elisa_od~primary_treatment + room , data=p.ab %>% filter(dpi == 14), family=Gamma())
pr4 <- glm(elisa_od~primary_treatment * room , data=p.ab %>% filter(dpi == 14), family=Gamma())
pr5 <- glm(elisa_od~1, data=p.ab %>% filter(dpi == 14), family=Gamma())
pr6 <- glm(elisa_od~primary_treatment + mass, data=p.ab %>% filter(dpi == 14), family=Gamma())
pr7<- glm(elisa_od~primary_dose, data=p.ab %>% filter(dpi == 14), family=Gamma())

#AICc
aictab(cand.set=list(pr1, pr2, pr2.5, pr3, pr4, pr5, pr6, pr7), modnames=c("pr1","pr2", "pr2.5", "pr3", "pr4", "pr5", "pr6", "pr7"))

#Model selection based on AICc:
#       K    AICc Delta_AICc AICcWt Cum.Wt     LL
#pr1    4 -637.45       0.00   0.54   0.54 322.89
#pr2    5 -635.41       2.03   0.19   0.73 322.95
#pr3   10 -635.20       2.25   0.18   0.91 328.55
#pr2.5  7 -633.88       3.57   0.09   1.00 324.41
#pr7    3 -622.45      15.00   0.00   1.00 314.32
#pr4   22 -609.07      28.37   0.00   1.00 331.40
#pr5    2 -553.23      84.22   0.00   1.00 278.66
#pr6   49 -545.16      92.29   0.00   1.00 353.40

#pr1 is best model

pr1
summary(pr1)
plot(allEffects(pr1))

#Primary treatment predicts antibody levels on day 14/15

#                       Estimate Std. Error t value  Pr(>|t|)    
#(Intercept)             22.372      1.406   15.91  < 2e-16 ***
#primary_treatmentlow    -5.949      1.608   -3.70 0.000322 ***
#primary_treatmenthigh  -11.742      1.492   -7.87  1.5e-12 ***
  
summary(pr2) #sex does not influence

#Primary dose predicts antibody levels on day 14 post inoculation (n=127, Estimate = -11.74, SE=1.49)
t1 <- m.ab%>%  filter(dpi == 14) %>%
  group_by( primary_treatment)%>%
  summarize(count = n())
t1
primary_treatment_14_name_n = c(
  'sham' = 'Sham (n= 51)',
  'low' = 'Low (n= 51)',
  'high' = 'High (n=53)'
  )

#####model predictions: primary_treatment####
dat.new=expand.grid(primary_treatment=unique(p.ab$primary_treatment))#new grid to put predictions into
dat.new$yhat = predict(pr1, type="response", newdata=dat.new, re.form=NA) #predicted values based off pr1 model
head(dat.new)

#plot predicted values over raw data
pid14.pred <- ggplot(data=p.ab %>% filter(dpi == 14), aes(x=primary_treatment, y=elisa_od, shape=primary_treatment))+
  geom_jitter(size=2, width=0.1)+
  geom_point(data=dat.new, aes(x=primary_treatment, y=yhat, shape=primary_treatment), size=10, shape="-", color="brown")+#model predictions
  scale_shape_manual(values = c(16, 17, 15), labels = primary_treatment_14_name_n)+
  labs(title="MG Antibodies DPI 14", x="Primary Treatment", y="ELISA OD", shape="Primary Treatment")+
  scale_x_discrete(labels = primary_treatment_14_name_n)

pid14.pred

######primary_dose - Primary Treatment treats the inoculation dose as a factor whereas primary_dose maintains the numerical values####

#model comparison
prd1 <- glm(elisa_od~primary_dose, data=p.ab %>% filter(dpi == 14), family=Gamma())
prd2<- glm(elisa_od~primary_dose + sex , data=p.ab %>% filter(dpi == 14), family=Gamma())
prd2.5<- glm(elisa_od~primary_dose * sex , data=p.ab %>% filter(dpi == 14), family=Gamma())
prd3 <- glm(elisa_od~primary_dose + room , data=p.ab %>% filter(dpi == 14), family=Gamma())
prd4 <- glm(elisa_od~primary_dose * room , data=p.ab %>% filter(dpi == 14), family=Gamma())
prd5 <- glm(elisa_od~1, data=p.ab %>% filter(dpi == 14), family=Gamma())
prd6 <- glm(elisa_od~primary_dose + mass, data=p.ab %>% filter(dpi == 14), family=Gamma())

#AICc
aictab(cand.set=list(prd1, prd2, prd2.5, prd3, prd4, prd5, prd6), modnames=c("prd1", "prd2", "prd2.5", "prd3", "prd4", "prd5", "prd6"))

#Model selection based on AICc:
  
#  K    AICc Delta_AICc AICcWt Cum.Wt     LL
#prd1    3 -622.45       0.00   0.51   0.51 314.32
#prd2.5  5 -620.44       2.01   0.19   0.69 315.47
#prd2    4 -620.43       2.02   0.18   0.88 314.38
#prd3    9 -619.64       2.81   0.12   1.00 319.59
#prd4   15 -606.24      16.20   0.00   1.00 320.28
#prd5    2 -553.23      69.22   0.00   1.00 278.66
#prd6   48 -529.45      93.00   0.00   1.00 342.88

#prd1 best model
summary(prd1)
plot(allEffects(prd1))

#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    1.833e+01  7.764e-01  23.614  < 2e-16 ***
# primary_dose  -2.574e-04  3.194e-05  -8.058 5.29e-13 ***

#Primary dose predicts antibody levels on day 14/15 post inoculation (n=127, Estimate = -2.57e-4, SE=3.19e-5)
length(which(p.ab$dpi == 14))

#####model predictions: primary_dose####
dat.new=expand.grid(primary_dose=unique(p.ab$primary_dose))#new grid to put predictions into
dat.new$yhat = predict(prd1, type="response", newdata=dat.new, re.form=NA) #predicted values based off prd1 model
head(dat.new)


#continuous
dat.new2 <- data.frame(primary_dose = seq(min(p.ab$primary_dose), max(p.ab$primary_dose), length.out = 10))
dat.new2$yhat <- predict(prd1, type = "response", interval = "confidence", level = 0.95, newdata = dat.new2)
head(dat.new2)

#plot predicted values over raw data
pid14.dose.pred <- ggplot(data=p.ab %>% filter(dpi == 14), aes(x=as.factor(primary_dose), y=elisa_od, shape=as.factor(primary_dose)))+
  geom_jitter(size=2, width=0.1)+
  geom_point(data=dat.new, aes(x=as.factor(primary_dose), y=yhat, shape=primary_dose), size=10, shape="-", color="brown")+#model predictions
  labs(x="Primary Dose", y="ELISA OD", shape="Primary Dose")

pid14.dose.pred

#plot continuous predicted values over raw data
pid14.dose.pred.c <- ggplot(data=p.ab %>% filter(dpi == 14), aes(x=primary_dose, y=elisa_od))+
  geom_jitter(size=2, width=0.1)+
  geom_line(data=dat.new2, aes(x=primary_dose, y=yhat), color="brown")+#model predictions
  labs(x="Primary Dose", y="ELISA OD", shape="Primary Dose")

pid14.dose.pred.c

##### (2) Does inoculation dose predict antibody levels on day 41?####
#gamma distribution because antibody data looks exponential
#glm b/c only using one dpi so no need for bird_id as a fixed effect
lm3 <- glm(elisa_od~primary_treatment, data=p.ab %>% filter(dpi==41), family=Gamma())
summary(lm3)
lm3
simulateResiduals(lm3, plot=T)

#Conditional model:
#                          Estimate   Std. Error z value Pr(>|z|)    
#(Intercept)            22.2469     0.5685  39.130  < 2e-16 ***
#primary_treatmentlow   -1.4552     0.7782  -1.870   0.0634 .  
#primary_treatmenthigh  -4.9142     0.7130  -6.893 1.42e-10 ***

#Birds that received a high primary dose still had elevated antibody levels on day 41 post infection (n=153, Estimate = -4.91, SE = 0.71).
lm3
plot(allEffects(lm3))
hist(resid(lm3))
l3r <- simulateResiduals(lm3)
plot(l3r)
qqnorm(resid(lm3))
qqline(resid(lm3))

#model comparison
p1 <- glm(elisa_od~primary_treatment, data=p.ab %>% filter(dpi==41), family=Gamma())
p2<- glm(elisa_od~primary_treatment + sex , data=p.ab %>% filter(dpi==41), family=Gamma())
p3<- glm(elisa_od~primary_treatment * sex , data=p.ab %>% filter(dpi==41), family=Gamma())
p4 <- glm(elisa_od~primary_treatment + room , data=p.ab %>% filter(dpi==41), family=Gamma())
p5 <- glm(elisa_od~primary_treatment * room , data=p.ab %>% filter(dpi==41), family=Gamma())
p6 <- glm(elisa_od~primary_treatment + room + quantity, data=p.ab %>% filter(dpi==41), family=Gamma())
p7 <- glm(elisa_od~primary_treatment + quantity , data=p.ab %>% filter(dpi==41), family=Gamma())
p8 <- glm(elisa_od~1, data=p.ab %>% filter(dpi==41), family=Gamma())
p9 <- glm(elisa_od~primary_treatment + mass, data= p.ab %>% filter(dpi==41), family=Gamma())
p10<- glm(elisa_od~primary_dose, data=p.ab %>% filter(dpi==41), family=Gamma())

#AICc
aictab(cand.set=list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10), modnames=c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10"))
#p4 is the best model
#Model selection based on AICc:
  
#  K     AICc Delta_AICc AICcWt Cum.Wt     LL
#p4  10 -1049.98       0.00   0.69   0.69 535.76
#p6  11 -1048.18       1.80   0.28   0.97 536.03
#p5  22 -1043.77       6.21   0.03   1.00 547.78
#p2   5 -1037.41      12.57   0.00   1.00 523.91
#p1   4 -1036.88      13.10   0.00   1.00 522.57
#p10  3 -1035.19      14.79   0.00   1.00 520.68
#p7   5 -1034.75      15.23   0.00   1.00 522.58
#p3   7 -1033.73      16.25   0.00   1.00 524.25
#p8   2  -984.45      65.53   0.00   1.00 494.27
#p9  50  -966.53      83.46   0.00   1.00 558.26

summary(p4)
plot(allEffects(p4))

#table w/ number of infected birds by room
t1 <- m.ab%>%
  filter(dpi == 41) %>%
  group_by(room, primary_treatment)%>%
  summarize(count = n())
print(t1)
#There are more birds in room a which means that the average elisa_od is going to be higher compared to the other rooms
#This should explain why there is an effect of room in the models 

#do not include room in models
#AICc
aictab(cand.set=list(p1, p2, p3, p7, p8, p9, p10), modnames=c("p1", "p2", "p3", "p7", "p8", "p9", "p10"))
#p2 is the best model

#Model selection based on AICc:
  
#     K     AICc Delta_AICc AICcWt Cum.Wt     LL
#p2   5 -1037.41       0.00   0.40   0.40 523.91
#p1   4 -1036.88       0.53   0.30   0.70 522.57
#p10  3 -1035.19       2.21   0.13   0.83 520.68
#p7   5 -1034.75       2.66   0.10   0.94 522.58
#p3   7 -1033.73       3.67   0.06   1.00 524.25
#p8   2  -984.45      52.95   0.00   1.00 494.27
#p9  50  -966.53      70.88   0.00   1.00 558.26

p2
summary(p2)
plot(allEffects(p2))

#                           Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)            16.9198    0.5110   33.111     < 2e-16 ***
#  primary_treatmentlow    3.4652     0.6806   5.091     1.06e-06 ***
#  primary_treatmentsham   4.9199     0.7098   6.932     1.17e-10 ***
#  sexm                    0.8300     0.5742   1.445     0.15    


#facet_wrap labels
sex_names <- c(
              'm' = "Male",
              'f' = "Female"
              )
#Antibodies on DPI 41
g.ab <- ggplot(data = m.ab %>% filter(dpi == 41), aes(x = fct_rev(primary_treatment), y = elisa_od, shape = primary_treatment)) +
  geom_jitter(size=1.5, width = 0.25) +
  stat_summary(aes(group=primary_treatment), fun=mean, geom="point", alpha=1, size=10, shape="-", color="red")+
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype='dashed') +
  labs(x="Inoculation Dose", y="MG Antibodies OD", shape = "Primary Treatment")+
  facet_wrap(~sex, labeller=as_labeller(sex_names))
g.ab

####model predictions: primary_treatment####

#with sex
dat.new=expand.grid(primary_treatment=unique(p.ab$primary_treatment),
                    sex=unique(p.ab$sex))#new grid to put predictions into
dat.new$yhat = predict(p2, type="response", newdata=dat.new, re.form=NA) #predicted values based off glm.phago model
head(dat.new)

#plot predicted values over raw data
pid41.sex.pred <- ggplot(data = m.ab %>% filter(dpi == 41), aes(x = fct_rev(primary_treatment), y = elisa_od, shape = sex)) +
  geom_jitter(size=1.5, width = 0.25) +
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype='dashed') +
  ylab("MG Antibodies OD") +
  xlab("Inoculation Dose")+
  geom_point(data=dat.new, aes(x=primary_treatment, y=yhat, shape=sex), size=10, shape="-", color="brown4")+#model predictions
  facet_wrap(~sex, labeller=as_labeller(sex_names))


pid41.sex.pred

#without sex
dat.newer=expand.grid(primary_treatment=unique(p.ab$primary_treatment))#new grid to put predictions into
dat.newer$yhat = predict(p1, type="response", newdata=dat.newer, re.form=NA) #predicted values based off p1.5 model
head(dat.newer)

#plot predicted values over raw data
pid41.pred <- ggplot(data = m.ab %>% filter(dpi == 41), aes(x = fct_rev(primary_treatment), y = elisa_od, shape = primary_treatment)) +
  geom_jitter(size=1.5, width = 0.25) +
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype='dashed') +
  labs(x="Inoculation Dose", y= "MG Antibodies OD", shape = "Primary Treatment")+
  geom_point(data=dat.newer, aes(x=primary_treatment, y=yhat, shape=primary_treatment), size=10, shape="-", color="brown4")#model predictions

pid41.pred

#graph showing antibody levels on dpi 14/15 by primary_treatment with mean, error, and model predictions
ggplot(data=p.ab %>% filter(dpi== 41), aes(x=as.factor(primary_treatment), y= elisa_od, shape=primary_treatment))+
  geom_jitter(width=0.1)+
  geom_hline(yintercept = 0.061, linetype="dashed")+
  stat_summary(aes(group=primary_treatment), fun=mean, geom="point", alpha=1, size=3, shape="x", color="red")+
  stat_summary(aes(group=primary_treatment, shape = primary_treatment), fun.y=mean, color="red",
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               geom= "errorbar", size=0.5, width=0.25)+
  geom_point(data=dat.newer, aes(x=as.factor(primary_treatment), y=yhat, shape=primary_treatment), size=10, shape="-", color="brown4")+#model predictions
  labs(title = "Elisa OD DPI 41", y="ELISA OD", x= "Primary Treatment", shape="Primary Treatment")

######primary_dose - Primary Treatment treats the inoculation dose as a factor whereas primary_dose maintains the numerical values####

#model comparison
pr1 <- glm(elisa_od~primary_dose, data=p.ab %>% filter(dpi==41), family=Gamma())
pr2<- glm(elisa_od~primary_dose + sex , data=p.ab %>% filter(dpi==41), family=Gamma())
pr3<- glm(elisa_od~primary_dose * sex , data=p.ab %>% filter(dpi==41), family=Gamma())
pr4 <- glm(elisa_od~primary_dose + quantity , data=p.ab %>% filter(dpi==41), family=Gamma())
pr5 <- glm(elisa_od~1, data=p.ab %>% filter(dpi==41), family=Gamma())
pr6 <- glm(elisa_od~primary_dose + mass, data= p.ab %>% filter(dpi==41), family=Gamma())


#AICc
aictab(cand.set=list(pr1, pr2, pr3, pr4, pr5, pr6), modnames=c("pr1", "pr2", "pr3", "pr4", "pr5", "pr6"))
#Model selection based on AICc:
  
#  K     AICc Delta_AICc AICcWt Cum.Wt     LL
#pr2  4 -1035.69       0.00   0.41   0.41 521.98
#pr1  3 -1035.19       0.49   0.32   0.74 520.68
#pr3  5 -1033.64       2.04   0.15   0.89 522.03
#pr4  4 -1033.10       2.59   0.11   1.00 520.69
#pr5  2  -984.45      51.24   0.00   1.00 494.27
#pr6 49  -968.11      67.58   0.00   1.00 556.84

#pr2 best model
summary(pr2)
plot(allEffects(pr2))

#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)  21.1491420  0.4849462  43.611  < 2e-16 ***
#  primary_dose -0.0001413  0.0000197  -7.173 3.13e-11 ***
#  sexm          0.8305667  0.5811476   1.429    0.155    

#Primary dose predicts antibody levels on day 41 post inoculation when controlling for sex (n=153, Estimate = 18.3, SE=0.776)

#without sex
summary(pr1)
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  21.5563342  0.3982652   54.13  < 2e-16 ***
#  primary_dose -0.0001411  0.0000198   -7.13 3.88e-11 ***

#Primary dose predicts antibody levels on day 41 post inoculation (n=153, Estimate = 21.56, SE = 0.398)
#####model predictions: primary_dose####
dat.new=expand.grid(primary_dose=unique(p.ab$primary_dose),
                    sex=unique(p.ab$sex))#new grid to put predictions into
dat.new$yhat = predict(pr2, type="response", newdata=dat.new, re.form=NA) #predicted values based off pr2
head(dat.new)

#plot predicted values over raw data
pid41.dose.sex.pred <- ggplot(data=p.ab %>% filter(dpi == 14), aes(x=as.factor(primary_dose), y=elisa_od, shape=as.factor(primary_dose)))+
  geom_jitter(size=2, width=0.1)+
  geom_point(data=dat.new, aes(x=as.factor(primary_dose), y=yhat, shape=primary_dose), size=10, shape="-", color="brown")+#model predictions
  labs(x="Primary Dose", y="ELISA OD", shape="Primary Dose")+
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype='dashed') +
  facet_wrap(~sex, labeller=as_labeller(sex_names))

pid41.dose.sex.pred

#continuous model predictions
dat.newer2 <- data.frame(primary_dose = seq(min(p.ab$primary_dose), max(p.ab$primary_dose), length.out = 10),
                         sex=unique(p.ab$sex))
dat.newer2$yhat <- predict(pr2, type = "response", interval = "confidence", level = 0.95, newdata = dat.newer2)
head(dat.newer2)

pid41.dose.sex.pred.c <- ggplot(data=p.ab %>% filter(dpi == 14), aes(x=primary_dose, y=elisa_od))+
  geom_jitter(size=2, width=0.1)+
  geom_line(data=dat.newer2, aes(x=primary_dose, y=yhat, color=sex))+#model predictions
  labs(x="Primary Dose", y="ELISA OD", shape="Primary Dose")+
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype='dashed')# +
  facet_wrap(~sex, labeller=as_labeller(sex_names))

pid41.dose.sex.pred.c

#no sex
dat.newer3 <- data.frame(primary_dose = seq(min(p.ab$primary_dose), max(p.ab$primary_dose), length.out = 100))
dat.newer3$yhat <- predict(pr1, type = "response", interval = "confidence", level = 0.95, newdata = dat.newer3)
head(dat.newer3)

pid41.dose.pred.c <- ggplot(data=p.ab %>% filter(dpi == 14), aes(x=primary_dose, y=elisa_od))+
  geom_jitter(size=2, width=0.1)+
  geom_line(data=dat.newer3, aes(x=primary_dose, y=yhat), color="brown")+#model predictions
  labs(title = "ELISA OD PID 41", x="Primary Dose", y="ELISA OD", shape="Primary Dose")+
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype='dashed')

pid41.dose.pred.c

#plot predicted values over raw data
pid41.dose.pred <- ggplot(data=p.ab %>% filter(dpi == 14), aes(x=as.factor(primary_dose), y=elisa_od, shape=as.factor(primary_dose)))+
  geom_jitter(size=2, width=0.1)+
  geom_point(data=dat.new, aes(x=as.factor(primary_dose), y=yhat, shape=primary_dose), size=10, shape="-", color="brown")+#model predictions
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype='dashed') +
  labs(x="Primary Dose", y="ELISA OD", shape="Primary Dose")

pid41.dose.pred


#### (3) Do primary and secondary inoculation dose predict antibody levels on day 56?####
table(m.ab$primary_treatment)
table(m.ab$secondary_treatment)
table(m.ab$sex)
table(m.ab$dpi)
m.ab %>% 
  filter(dpi==56)%>%
  duplicated(select = band_number)

# Filter the data for dpi == 56
filtered_data <- m.ab %>% filter(dpi == 14)

# Identify band_numbers occurring more than once
band_numbers_more_than_once <- filtered_data %>%
  group_by(band_number) %>%
  filter(n() > 1) %>%
  pull(band_number)

# Display the list of band_numbers
band_numbers_more_than_once

summary(m.ab$elisa_od)
table(m.ab$secondary_treatment)
m.ab$secondary_dose <- as.integer(m.ab$secondary_dose)
table(m.ab$secondary_dose, m.ab$dpi)
table(m.ab$secondary_dose, m.ab$dpsi)
m.ab$primary_treatment <- factor(m.ab$primary_treatment, levels = c("sham", "low", "high"))

#use categorical primary treatment, but keep numerical secondary dose because the secondary dose is of interest
#new df with only dpsi 14
s.ab <- m.ab %>%
  filter(dpsi > 0 & elisa_od !=0)
s.ab
m.ab$dpsi <- as.numeric(m.ab$dpsi)
table(m.ab$dpsi)

lm4 <- glm(elisa_od~primary_treatment * secondary_dose, data=s.ab, family=Gamma())
summary(lm4)
plot(allEffects(lm4))
plot(resid(lm4))
qqnorm(resid(lm4))
qqline(resid(lm4))
anova(lm4, "III")
residlm4 <- simulateResiduals(lm4)
plot(residlm4)


#Inoculation with a high dose results in a more robust antibody response at lower secondary doses.

#Model Comparison
s1 <- glm(elisa_od~primary_treatment*secondary_dose, data=s.ab, family=Gamma())
s2 <- glm(elisa_od~primary_treatment+secondary_dose, data=s.ab, family=Gamma())
s3 <- glm(elisa_od~primary_treatment*secondary_dose + sex, data=s.ab, family=Gamma())
s4 <- glm(elisa_od~primary_treatment+secondary_dose+sex, data=s.ab, family=Gamma())
s5 <- glm(elisa_od~1, data=s.ab, family=Gamma())

#AICc
aictab(cand.set=list(s1, s2, s3, s4, s5), modnames=c("s1", "s2", "s3", "s4", "s5"))

#Model selection based on AICc:
  
#  K    AICc Delta_AICc AICcWt Cum.Wt     LL
#s1 7 -830.07       0.00   0.51   0.51 422.46
#s3 8 -830.01       0.06   0.49   1.00 423.56
#s4 6 -810.01      20.07   0.00   1.00 411.32
#s2 5 -809.90      20.17   0.00   1.00 410.18
#s5 2 -781.93      48.14   0.00   1.00 393.01

#model 1 best
summary(s1)
plot(allEffects(s1))
summary(s3)
plot(allEffects(s3))


ggplot(m.ab %>% filter(elisa_cv !=0 & dpi == 14 | dpi == 41), aes(x=primary_treatment,
                                   y=elisa_od, color = primary_treatment, shape = primary_treatment))+
  geom_jitter(size=1.5, width = 0.25)+
  stat_summary(fun.y=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
               #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
               geom= "errorbar", size=0.5, width=.25)+
  labs(title = "Antibody OD DPI 14 + 41", x="Primary Treatment", y="ELISA OD", shape="Primary Treatment", color = "Primary Treatment")+
  stat_summary(aes(group=primary_treatment), color="black", fun=mean, shape = "-", size=2)+
  geom_hline(yintercept=0.061, linetype="dashed")+
  scale_colour_manual(values=c( "ivory4", "blue", "red"))+
  #scale_color_manual(values=c( "white", "gray75", "gray60", "gray50", "gray20"),
  #                   labels = c("0", "30", "100", "300", "7000"))#+
  #facet_wrap(~fct_rev(primary_treatment), ncol=5)+
  facet_wrap(~dpi)+
  theme_bw()

#Secondary infection by primary infection  
ggplot(s.ab, aes(x=primary_treatment,
                                   y=elisa_od, color = as.factor(secondary_dose), shape = as.factor(secondary_dose)))+
  geom_jitter(size=1.5, width = 0.25)+
  stat_summary(fun.y=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
               #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
               geom= "errorbar", size=0.2, width=.25)+
  labs(title = "Antibody OD DPI 56", x="Secondary Dose", y="ELISA OD", shape="Secondary Dose", color = "Secondary Dose")+
  #stat_summary(aes(group=primary_treatment), color="black", fun=mean, shape = "-", size=2)+
  geom_hline(yintercept=0.061, linetype="dashed")+
  #scale_color_manual(values=c( "white", "gray75", "gray60", "gray50", "gray20"),
  #                   labels = c("0", "30", "100", "300", "7000"))#+
  #facet_wrap(~fct_rev(primary_treatment), ncol=5)+
  theme_bw()



s1 <- glm(elisa_od~primary_treatment*secondary_dose, data=m.ab %>% filter(dpi == 56), family=Gamma())

dat.new=expand.grid(secondary_dose=seq(min(p.ab$secondary_dose), max(p.ab$secondary_dose), length.out = 100),
                    primary_treatment=unique(m.ab$primary_treatment),
                    elisa_od=unique(m.ab$elisa_od),
                    sex=unique(m.ab$sex),
                    band_number=unique(m.ab$band_number))#new grid to put predictions into
dat.new$yhat = predict(s1, type="response", newdata=dat.new, re.form=NA) #predicted values
head(dat.new)


#plot showing predicted
ggplot(m.ab %>% filter(dpi == 56), aes(y=elisa_od, x=secondary_dose))+
  geom_point(aes(y=elisa_od, x=secondary_dose, shape = as.factor(secondary_dose), color=primary_treatment))+
  stat_summary(aes(group=secondary_treatment), fun=mean, shape = 4, color="brown4")+
  geom_line(data=dat.new, aes(x=secondary_dose, y=yhat, color=primary_treatment, group = primary_treatment))+
  labs(title = "Primary Treatment x Secondary Dose Model Predictions", 
       x="Secondary Dose", y="ELISA OD Day 56", color = "Primary Treatment", shape = "Secondary Dose")+
  #geom_hline(yintercept=0.061, linetype="dashed")+
  scale_x_continuous(breaks = unique(m.ab$secondary_dose), labels = unique(m.ab$secondary_dose))#+
  facet_wrap(~fct_rev(primary_treatment))

####Do antibody levels on day 14 predict infection status on Day 56?####

    
#make a dataframe with elisa_od from dpi 14 for every band_number and infection status from dpi 56 for every band_number
m.ab1<-m.ab %>%
  dplyr::select(dpi, quantity, elisa_od, band_number, infected, ever_infected_qpcr, 
                primary_treatment, secondary_treatment, primary_dose, secondary_dose, sex) #make new df with just relevant variables

#subset data for dpi 14 and 56
sub.ab <- m.ab1 %>%
  filter(dpi %in% c("14", "56"))
head(sub.ab)
colnames(sub.ab)
sub.ab$dpi <- as.numeric(sub.ab$dpi)
sub.ab$infection_56 <- ifelse(dpi == 56 & !is.na(infected) & infected == 1, 1, 0)
  
elisa14 <- m.ab1 %>%
  filter(dpi == 14)%>% #dpi 14 for elisa
  dplyr:: select(elisa_od, band_number, dpi, infected, ever_infected_qpcr, primary_treatment,
                 secondary_treatment, primary_dose, secondary_dose, sex) #df with only elisa od, band number, dpi from 14

str(m.ab1)
quant56 <- m.ab1 %>%
  filter(dpi == 56)%>% #dpi 56 only for path load
  dplyr:: select(quantity, band_number, dpi, infected, ever_infected, 
                 primary_treatment, secondary_treatment, primary_dose, secondary_dose, sex) #df with only quantity, band number, dpi from d56
    

#add column to elisa14 df "from quant56 df row quantity, 
#take all data where the band number from elisa14 match band number from quant56
#returns new column with missing points as NA
elisa14$quantity.d56 <- quant56$quantity[match(elisa14$band_number,quant56$band_number)]
elisa14$dpi.check <- quant56$dpi[match(elisa14$band_number,quant56$band_number)] #check band numbers to see if match works
elisa14$dpi.check.no <- quant56$band_number[match(elisa14$band_number,quant56$band_number)] #check band numbers to see if match works
elisa14$infected.d56 <-quant56$infected[match(elisa14$band_number,quant56$band_number)] #add variable specifying infection d56
elisa14
mbirds <- m.ab %>%
  filter(dpi == 56 | dpi == 14) %>%
  select(dpi, band_number, elisa_od, quantity)


mbirds %>%
  group_by(dpi)%>%
  count(dpi)

elisa14 %>%
  count(elisa_od, quantity.d56[quantity.d56 == "NA"])




# Get all possible band_number values
all_band_numbers <- m.ab %>% pull(band_number) %>% unique()

# Create a data frame with all combinations of dpi and band_number
all_combinations <- expand.grid(dpi = unique(elisa14$dpi), band_number = all_band_numbers)

# Merge with the actual data to identify missing band_numbers
missing_band_numbers <- all_combinations %>%
  anti_join(elisa14, by = c("dpi", "band_number"))

# Display the missing band_numbers
print(missing_band_numbers)

to_rerun <- ab %>%
  filter(Avg.OD == 0)
table(to_rerun$dpi)

#do antibody levels on day 14 predict infection status on day 56?
elisa14 <- elisa14%>%
  drop_na(quantity.d56)

lm6 <- glm(infected.d56~elisa_od,data=elisa14, family=nbinom1())
summary(lm6)
plot(allEffects(lm6))
hist(resid(lm6))



#model selection
i1 <- glm(infected.d56~elisa_od ,data=elisa14, family=nbinom1())
i2 <- glm(infected.d56~elisa_od + primary_treatment ,data=elisa14, family=nbinom1())
i3 <- glm(infected.d56~elisa_od * primary_treatment,data=elisa14, family=nbinom1())
i4 <- glm(infected.d56~elisa_od + sex ,data=elisa14, family=nbinom1())
i5 <- glm(infected.d56~elisa_od + secondary_treatment,data=elisa14, family=nbinom1())
i6 <- glm(infected.d56~elisa_od + primary_treatment * secondary_treatment,data=elisa14, family=nbinom1())
i7 <- glm(infected.d56~elisa_od + primary_treatment + secondary_treatment ,data=elisa14, family=nbinom1())
i8 <- glm(infected.d56~1 ,data=elisa14, family=nbinom1())

summary(i1)
summary(i7)
plot(allEffects(i7))

#why isn't this working??
aictab(cand.set=list(i1, i2, i3, i4, i5, i6, i7, i8), modnames=c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8"))


#infd56 vs ab14
dat.new=expand.grid(primary=unique(elisa41$secondary_treatment),
                    secondary_dose=unique(elisa41$secondary_dose),
                    infected.d56=unique(elisa41$infected.d56),
                    elisa_od=unique(elisa41$elisa_od),
                    band_number=unique(elisa41$band_number))#new grid to put predictions into
dat.new$yhat = predict(lm3, type="response", newdata=dat.new, re.form=NA) #predicted values based off glm.phago model
head(dat.new)

#plot showing predicted
ggplot(elisa41, aes(y=infected.d56, x=elisa_od, color=as.factor(primary_treatment)))+
  geom_point(shape=1, size = 5)+
  stat_smooth(method = "glm",
              method.args = list(family="binomial"), se=FALSE,
              fullrange = TRUE)+
  labs(x="ELISA OD Day 41", y="Infection Status Day 56")+
  geom_vline(xintercept=0.061, linetype="dashed")+
  xlim(c(0.04, 0.2))+
  facet_wrap(~secondary_dose)


####Do antibody levels on day 41 predict infection status on day 56?####
m.ab1<-m.ab %>%
  dplyr::select(dpi, quantity, elisa_od, band_number, infected, ever_infected, 
                primary_treatment, secondary_treatment, primary_dose, secondary_dose, sex) #make new df with just relevant variables


elisa41 <- m.ab %>%
  filter(dpi == 41)%>% #dpi 41 for elisa
  dplyr:: select(elisa_od, band_number, dpi, infected, primary_treatment, secondary_treatment, primary_dose, secondary_dose) #df with only elisa od, band number, dpi from 41

quant56 <- m.ab %>%
  filter(dpi == 56 & band_number != 2375)%>% #dpi 56 only for path load; drop 2375 because there was no elisa for it on dpi 41
  select(quantity, band_number, dpi, infected, ever_infected, 
                 primary_treatment, secondary_treatment, primary_dose, secondary_dose, sex) #df with only quantity, band number, dpi from d56

summary(quant56)
summary(elisa41)

mvt <- master %>%
  filter(experiment_location == "vt")

table(mvt$dppi)
table(m.ab$dpi)
#add column to elisa41 df "from quant56 df row quantity, 
#take all data where the band number from elisa41 match band number from quant56
#returns new column with missing points as NA
elisa41$quantity_dpi_56 <- quant56$quantity[match(elisa41$band_number,quant56$band_number)]
elisa41$dpi.check <- quant56$dpi[match(elisa41$band_number,quant56$band_number)] #check band numbers to see if match works
elisa41$dpi.check.no <- quant56$band_number[match(elisa41$band_number,quant56$band_number)] #check band numbers to see if match works
#elisa41$infected.d56 <-quant56$infected[match(elisa41$band_number,quant56$band_number)] #add variable specifying infection d56

elisa41$threshold_cutoff = 50


elisa41 <- elisa41 %>%
  mutate(infected.d56 = ifelse(quantity_dpi_56>threshold_cutoff, 1, 0))

# Get all possible band_number values
all_band_numbers <- m.ab %>% pull(band_number) %>% unique()


# Create a data frame with all combinations of dpi and band_number
all_combinations <- expand.grid(dpi = unique(elisa41$dpi), band_number = all_band_numbers)


# Merge with the actual data to identify missing band_numbers
missing_band_numbers <- all_combinations %>%
  anti_join(elisa41, by = c("dpi", "band_number"))

# Display the missing band_numbers
print(missing_band_numbers)
#2375 missing - didn't have enough plasma on pid 41
summary(elisa41)


# Get all possible band_number values
all_band_numbers <- m.ab %>% pull(band_number) %>% unique()

# Create a data frame with all combinations of dpi and band_number
all_combinations <- expand.grid(dpi = unique(quant56$dpi), band_number = all_band_numbers)

# Merge with the actual data to identify missing band_numbers
missing_band_numbers56 <- all_combinations %>%
  anti_join(quant56, by = c("dpi", "band_number"))

# Display the missing band_numbers
print(missing_band_numbers56)

#2435 2462, 2496 to be re-run
#2375 omitted b/c there wasn't plasma for an elisa on 41

pr<-m.ab %>%
  filter(dpi == 56) %>%
  select(quantity, band_number, elisa_od, Rerun.)
pr
count(pr)
ab <- m.ab %>%
  filter(dpi == 41)%>%
  select(elisa_od, band_number, primary_treatment)
ab
ab %>%
  group_by(primary_treatment)

elisa41 %>%
  drop_na(quantity)

#do antibody levels on d41 predict infection on day 56, controlling for secondary inoculation?
lm3 <- glm(infected.d56~elisa_od + secondary_dose + (1|band_number), data=elisa41, family="binomial")
summary(lm3)
plot(allEffects(lm3))

#does primary treatment predict infection on day 56, controlling for secondary inoculation?
lm3.5 <- glm(infected.d56~primary_treatment + secondary_dose + (1|band_number), data=elisa41, family="binomial")
summary(lm3.5)
plot(allEffects(lm3.5))

#do antibody levels on d14 predict infection on day 56?
lm4 <- glm(infected.d56~elisa_od + secondary_dose + (1|band_number), data=elisa41, family="binomial")
summary(lm4)
plot(allEffects(lm4))

lm3 <- glm(infected.d56~elisa_od + secondary_dose, data=elisa41, family="binomial")
summary(lm3)
#infd56 vs ab41
dat.new=expand.grid(secondary_treatment=unique(elisa41$secondary_treatment),
                    secondary_dose=unique(elisa41$secondary_dose),
                    primary_treatment=unique(elisa41$primary_treatment),
                    infected.d56=unique(elisa41$infected.d56),
                    elisa_od=unique(elisa41$elisa_od),
                    band_number=unique(elisa41$band_number))#new grid to put predictions into
dat.new$yhat = predict(lm3, type="response", newdata=dat.new, re.form=NA) #predicted values based off glm.phago model
head(dat.new)

#plot showing predicted values - Antibodies on dpi 41 predict infection probability on dpsi 14
ggplot(elisa41 %>% filter(infected.d56 != "NA"), aes(y=infected.d56, x=elisa_od))+
  geom_point(aes(color=fct_rev(as.factor(infected.d56))), shape=1, size = 3)+
  geom_line(data = dat.new, aes(x = elisa_od, y = yhat)) +
  #stat_smooth(method = "glm",
  #            method.args = list(family="binomial"), se=FALSE,
  #            fullrange = TRUE)+
  labs(x="ELISA OD Day 41", y="Infection Status Day 56", color = "Infection Status Day 56")+
  geom_vline(xintercept=0.061, linetype="dashed")+
  xlim(c(0.04, 0.13))+
  facet_wrap(~primary_treatment~secondary_dose, ncol=5)

#primary + secondary
lm3.5 <- glm(infected.d56~elisa_od + primary_dose * secondary_dose, data=elisa41, family="binomial")
summary(lm3.5)
plot(allEffects(lm3.5))
plot(resid(lm3.5))
qqnorm(resid(lm3.5))
qqline(resid(lm3.5))

lm3.75 <-glm(infected.d56~elisa_od, data=elisa41, family="binomial")
summary(lm3.75)

dat.new=expand.grid(secondary_treatment=unique(elisa41$secondary_treatment),
                    secondary_dose=unique(elisa41$secondary_dose),
                    primary_treatment=unique(elisa41$primary_treatment),
                    primary_dose=unique(elisa41$primary_dose),
                    infected.d56=unique(elisa41$infected.d56),
                    elisa_od=unique(elisa41$elisa_od),
                    band_number=unique(elisa41$band_number))#new grid to put predictions into
dat.new$yhat = predict(lm3.5, type="response", newdata=dat.new, re.form=NA) #predicted values based off glm.phago model
head(dat.new)

ggplot(elisa41%>% filter(infected.d56 != "NA"), aes(y=infected.d56, x=elisa_od))+
  geom_point((aes(color=fct_rev(as.factor(infected.d56)))), shape=1, size = 3)+
  geom_point(data = dat.new, aes(x = elisa_od, y = yhat, color=as.factor(secondary_dose))) +
  labs(x="ELISA OD Day 41", y="Infection Status Day 56", color = "Infection Status Day 56")+
  geom_vline(xintercept=0.061, linetype="dashed", alpha=0.3)+
  xlim(c(0.04, 0.13))+
  facet_wrap(~secondary_treatment)

ggplot(elisa41, aes(x=elisa_od, y=infected.d56, color = as.factor(primary_dose)))+
  geom_point()+
  geom_line(data = dat.new, aes(x = elisa_od, y = yhat, color = primary_treatment)) 
unique(m.ab$dpi)
#### (6) Does prior exposure induce heterogeneity in antibody levels?####
p.ab <- p.ab %>%
  filter(dpi !=56)
#subset antibodies by groups to look at their variance
#subset by primary treatment (number doesn't matter) with 1|band_number as random effect
#variance of random effect should be highest in highest heterogeneity groups
library(lme4)
#across all primary infection days, which groups are most heterogeneous?
glm.all <- glmer(elisa_od~1 +(1|band_number), data=p.ab, family=Gamma())
glm.sham <- glmer(elisa_od~1 +(1|band_number), data=subset(p.ab, primary_treatment == "sham"), family=Gamma())
glm.low <- glmer(elisa_od~1 +(1|band_number), data=subset(p.ab, primary_treatment == "low"), family=Gamma())
glm.high <- glmer(elisa_od~1 +(1|band_number), data=subset(p.ab, primary_treatment == "high"), family=Gamma())

#glm.all.nest<-glm.all <- glmer(elisa_od~1 +(1|band_number) +(1|primary_treatment), data=p.ab, family=Gamma())
#summary(glm.all.nest)
summary(glm.all); summary(glm.sham); summary(glm.low); summary(glm.high)
#options - use SD to get a CI around each four variance estimates OR
#add up the log lik from the 4models and calculate AIC

#only days 14 and 41
#glm.all: V = 11.5496; AIC = -2343.8; SD = 3.39885; LogLik = 1174.9
#glm.sham: V = 0.284649; AIC = -1099; SD = 0.534; LogLik = 552.5
#glm.low:  V = 8.65716; AIC = -864.6; SD = 2.9423; LogLik = 435.3
#glm.high: V = 5.8528; AIC = -707.7; SD = 2.4192; LogLik = 356.8

#days -8, 14, and 41
#glm.all: V = 10.89235; AIC = -2543.4; SD = 3.3004; LogLik = 1274.7
#glm.sham: V = 0.267613; AIC = -1245.6; SD = 0.51731; LogLik = 625.8
#glm.low:  V = 7.58692; AIC = -928.5; SD = 2.7544; LogLik = 467.3
#glm.high: V = 5.2232; AIC = -740.9; SD = 2.2854; LogLik = 373.5
AIC(glm.sham,glm.low, glm.high, glm.all)

##Days 14 and 41 only
#quick graph showing differences in variance by group
#var <- c(11.5496, 0.284649, 8.65716, 5.8528) #variance from table above
#ci <- c(3.39885, 0.534, 2.9423, 2.4192) #Standard Deviation from table above
#g <- c("All", "Sham", "Low", "High") #Groups

##Days -8, 14, and 41 
var <- c(10.89235, 0.267613, 7.58692, 5.2232) #variance from table above
ci <- c(3.3004, 0.51731, 2.7544, 2.2854) #Standard Deviation from table above
g <- c("All", "Sham", "Low", "High") #Groups

#new df with variance and error for error bars (min/max)
variability <- data.frame(g,var, ci)
variability$min <- variability$var - variability$ci #new column for lowerbound error
variability$max <- variability$var + variability$ci #new column for upperbound error

g.var <- ggplot(variability %>% filter(g != "All"), aes(x=fct_rev(g), y=var, color=g))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=min, ymax=max), width=0.5)+ #add errorbars representing +/- 1 SD
  labs(title = "Primary Infection Antibody Variance", x="Primary Treatment", y="Variance", color="Treatment")+
  scale_color_manual(values=c("gray15", "gray60", "gray85"))+
  theme_minimal()

g.var


m.ab%>%
  group_by(primary_treatment)%>%
  summarise(count=n_distinct(band_number))

m.ab%>%
  group_by(primary_treatment, dpi)%>%
  summarise(count=n_distinct(band_number))

####Residual Models####
unique(p.ab$dpi)
unique(p.ab$elisa_od)
p.ab$dpi.new <- p.ab$dpi + 9
#p.ab$band_number <- as.character(p.ab$band_number)
glm <- glmmTMB(elisa_od ~ (band_number|dpi), data=p.ab, family=Gamma())
glm.a <- glmmTMB(elisa_od ~ (dpi.new|band_number), data=p.ab, family=Gamma())
glm.b <- glmmTMB(elisa_od ~ (1|band_number) + (1|dpi), data=p.ab, family=nbinom2)
summary(glm)

p.ab$resid <- resid(glm)
hist(p.ab$resid)
lm1 <- lm(resid ~ primary_treatment, data=p.ab)

anova(lm1)
summary(lm1)

ggplot(p.ab, aes(x=primary_treatment, y=resid))+
  geom_jitter()+
  theme_bw()

#specify intercept as 0.045
ab.1441 <- p.ab %>%
  filter(dpi == 14 | dpi == 41)

glm.c <- glm(elisa_od ~ 1  + (band_number|dpi), data= ab.1441, family=Gamma())
ab.1441$resid <- resid(glm.c)

lm.a <- lm(resid~primary_treatment, data= ab.1441)

ggplot(ab.1441, aes(x=primary_treatment, y=resid))+
  geom_jitter(width=0.25)

anova(lm.a)
summary(lm.a)
#for only pid 41
ab.41 <- p.ab %>%
  filter(dpi == 41)

glm41 <- glmmTMB(elisa_od ~ 1, data=ab.41, family=Gamma())
ab.41$resid <- resid(glm41)

lm2 <- lm(resid ~ primary_treatment, data=ab.41)

anova(lm2)
summary(lm2)

ggplot(ab.41, aes(x=primary_treatment, y=resid))+
  geom_point()+
  theme_bw()

####Histograms####
#new df with labels for graphs - treatment + n
treat_names_n <- c(
  'high' = "High (n=53)",
  'low' = "Low (n=51)",
  'sham' = "Sham (n=50)"
)
#new df with labels for graphs - treatment only
treat_names <- c(
  'high' = "High",
  'low' = "Low",
  'sham' = "Sham"
)

#histogram showing distribution of antibodies across all dpi

#Add variable indicating the elisa_od mean for each primary_treatment
p.ab <- p.ab %>%
  group_by(primary_treatment)%>%
  mutate(groupmean = mean(elisa_od))

#Histogram showing distribution of elisa_od across all sample days by primary_treatment
ggplot(p.ab%>% filter(dpi != -8), aes(x=elisa_od, fill=primary_treatment))+
  #geom_density(color="black", alpha=0.5)+
  geom_histogram(binwidth = 0.005, position="identity", alpha=2, color="black", size=0.35)+
  geom_vline(xintercept=0.061, linetype="dashed", alpha=0.5)+
  xlim(c(0.04, 0.19))+
  geom_vline(data = p.ab %>% filter(dpi != -8), aes(xintercept = groupmean, color = "red"), linetype = "dotted",alpha = 1, show.legend=FALSE) +
  labs(y= "Count", x= "ELISA OD", fill="Primary Treatment", color= "Average OD")+
  scale_fill_manual(values=c( "white", "gray","black"), labels=treat_names)+
  scale_color_manual(values = c("red"), guide = guide_legend(title= NULL))+
  facet_wrap(~primary_treatment~dpi, ncol=2)+
  theme_bw()

#Add variable indicating the elisa_od mean for each primary_treatment
p.ab <- p.ab %>%
  group_by(primary_treatment, dpi)%>%
  mutate(groupmean.dpi = mean(elisa_od))

#Histogram showing 
ggplot(p.ab, aes(x=elisa_od, fill=primary_treatment))+
  #geom_density(color="black", alpha=0.5)+
  geom_histogram(binwidth = 0.005, position="identity", alpha=2, color="black", size=0.35)+
  geom_vline(xintercept=0.061, linetype="dashed", alpha=0.5)+
  xlim(c(0.04, 0.19))+
  geom_vline(data = p.ab, aes(xintercept = groupmean.dpi, color = "red"), linetype = "dotted",alpha = 1, show.legend=FALSE) +
  labs(y= "Count", x= "ELISA OD", fill="Primary Treatment", color= "Average OD")+
  scale_fill_manual(values=c("black", "gray", "white"), labels=treat_names)+
  scale_color_manual(values = c("red"), guide = guide_legend(title= NULL))+
  facet_wrap(~primary_treatment~dpi, ncol=2)+
  theme_bw()

#density plot showing distribution of antibodies across all dpi 
ggplot(p.ab, aes(x=elisa_od, fill=fct_rev(primary_treatment)))+
  geom_density(color="black", alpha=0.5, aes(fill=fct_rev(primary_treatment)))+ #fct_rev(primary_treatment) reverses order of overlay
  geom_vline(xintercept=0.061, linetype="dashed", alpha=0.5)+
  xlim(c(0.04, 0.19))+
  labs(y= "Count", x= "ELISA OD", fill="Primary Treatment")+
  #scale_fill_manual(values=c("white", "gray", "black"))+
  facet_wrap(~dpi, nrow=3)+
  theme_bw()

#density plot overlayed with histogram showing distribution of antibodies across all dpi
ggplot(p.ab, aes(x=elisa_od, fill=primary_treatment))+
  geom_density(color="black", alpha=0.5, aes(fill=primary_treatment))+ #fct_rev(primary_treatment) reverses order of overlay
  geom_histogram(aes(x=elisa_od, fill=fct_rev(primary_treatment)), binwidth = 0.001, position="identity", alpha=2, color="black", size=0.35)+
  geom_vline(xintercept=0.061, linetype="dashed", alpha=0.5)+
  xlim(c(0.04, 0.19))+
  labs(y= "Count", x= "ELISA OD", fill="Primary Treatment")#+
  scale_fill_manual(values=c("black", "gray", "white"))+
  facet_grid(c("primary_treatment", "dpi"),
             labeller = "label_both")+ #labeller = as_labeller(treat_names_n)
  theme_bw()


#On pid 41, which groups are most heterogeneous?
p.ab.41 <- p.ab %>%
  filter(dpi == 41)

#This is not working - I need to find a way to calculate variance for one day only with no repeated measures.
#nested random effect?
##This use of models doesn't really make sense to me - I want to look at between individual variability instead of within 
    #this works well for within individual variability, but not between.
glm.41.all <- glmer(elisa_od~1 + (1|dpi)+(1|primary_treatment), data=p.ab, family=Gamma())
glm.41.all <- glmmTMB(elisa_od~1 + (1|dpi), data=p.ab, family=Gamma())
glm.41.sham <- glmmTMB(elisa_od~1 +(1|dpi), data=subset(p.ab, primary_treatment == "sham"), family=Gamma())
glm.41.low <- glmmTMB(elisa_od~1 +(1|dpi), data=subset(p.ab, primary_treatment == "low"), family=Gamma())
glm.41.high <- glmmTMB(elisa_od~1 +(1|dpi), data=subset(p.ab, primary_treatment == "high"), family=Gamma())

#subset by dpi have primary_treatment as random effect


summary(glm.41.all); summary(glm.41.sham); summary(glm.41.low); summary(glm.41.high)
#options - use SD to get a CI around each four variance estimates OR
#add up the log lik from the 4models and calculate AIC

#41.all: V = ; AIC = ; SD = ; LogLik = 
#41.sham: V = ; AIC = ; SD = ; LogLik = 
#41.low:  V = ; AIC = ; SD = ; LogLik = 
#41.high: V = ; AIC = ; SD = ; LogLik = 

AIC(glm.41.sham,glm.41.low, glm.41.high, glm.41.all)

#density plot overlayed with histogram showing distribution of antibodies across all dpi

#####calculate cv by hand####
#calculate CV
m.cv.all <- p.ab %>% 
  group_by(primary_treatment)%>%
  summarise(bird_cv = sd(elisa_od)/mean(elisa_od), #this is variation within groups 
            bird_sd = sd(elisa_od))

#generate error bar values
m.cv.all$max <- m.cv.all$bird_cv + m.cv.all$bird_sd

m.cv.all$min <- m.cv.all$bird_cv - m.cv.all$bird_sd

#Graph of CV calculated from dpi -8 to 41 with error bars +/- 1 SD
cv.all.primary <- ggplot(m.cv.all, aes(x=primary_treatment, y=bird_cv, color=primary_treatment))+
  geom_col(aes(fill=primary_treatment), color="black", size=0.5)+
  geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.5)+
  labs(x="Primary Treatment", y="CV", fill="Treatment")+
  scale_fill_manual(values=c( "gray90", "gray50", "gray25"))+
  theme_minimal()

cv.all.primary

#CV of each bird over each day to compare to variance
m.cv.ind.birds <- p.ab %>% 
  group_by(band_number, primary_treatment)%>%
  summarise(bird_cv = sd(elisa_od)/mean(elisa_od), 
            bird_sd = sd(elisa_od))%>%
  ungroup()
m.cv.ind.birds

#generate error bar values
m.cv.ind.birds$max <- m.cv.ind.birds$bird_cv + m.cv.ind.birds$bird_sd

m.cv.ind.birds$min <- m.cv.ind.birds$bird_cv - m.cv.ind.birds$bird_sd

#Graph of CV calculated from dpi -8 to 41 with error bars +/- 1 SD
cv.ind.birds.primary <- ggplot(m.cv.ind.birds, aes(x=fct_rev(primary_treatment), y=bird_cv, color=primary_treatment))+
  geom_col(aes(fill=primary_treatment), color="black", size=0.5)+
  geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.5)+
  labs(x="Primary Treatment", y="CV", fill="Treatment")+
  scale_fill_manual(values=c( "gray25", "gray60", "gray90"),
                    labels = c("High", "Low", "Sham"))+
  theme_minimal()

cv.ind.birds.primary

#Calculate CV each day
m.cv <- p.ab %>% 
  group_by(dpi, primary_treatment)%>%
  summarise(bird_cv = sd(elisa_od)/mean(elisa_od),
            bird_sd = sd(elisa_od))
m.cv

#generate error bar values
m.cv$max <- m.cv$bird_cv + m.cv$bird_sd
m.cv$min <- m.cv$bird_cv - m.cv$bird_sd

cv.primary <- ggplot(m.cv %>% filter(dpi != -8), aes(x=fct_rev(primary_treatment), y=bird_cv, color=primary_treatment))+
  geom_col(aes(fill=primary_treatment), color="black", size=0.5)+
  geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.5)+
  labs(x="Primary Treatment", y="CV", fill="Treatment")+
  scale_fill_manual(values=c(  "gray90", "gray60", "gray25"),
                    labels = c("High", "Low", "Sham"))+
  facet_wrap(~dpi, ncol=5)+ 
  theme_bw()

cv.primary

#Calculate CV each day including secondary
s.cv <- m.ab %>% 
  group_by(dpi, primary_treatment, secondary_dose)%>%
  summarise(bird_cv = sd(elisa_od)/mean(elisa_od),
            bird_sd = sd(elisa_od))
s.cv

#generate error bar values
s.cv$max <- s.cv$bird_cv + s.cv$bird_sd
s.cv$min <- s.cv$bird_cv - s.cv$bird_sd

#Graph of CV calculated from dpi -8 to 56
cv.secondary <- ggplot(s.cv%>%filter(dpi==56), aes(x=as.factor(secondary_dose), y=bird_cv, color=primary_treatment))+
  geom_col(aes(fill=as.factor(secondary_dose)), color="black", size=0.5)+
  geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.5)+
  labs(x="Primary Treatment", y="CV", fill="Treatment")+
  scale_fill_manual(values=c( "white", "gray75", "gray60", "gray50", "gray20"),
                    labels = c("0", "30", "100", "300", "7000"))+
  facet_wrap(~fct_rev(primary_treatment), ncol=5)+
  theme_bw()

cv.secondary

#Graph of elisa_od dpi 56
elisa.secondary <- ggplot(m.ab%>%filter(dpi==56), aes(x=as.factor(secondary_dose), y=elisa_od, color=as.factor(secondary_dose), shape = as.factor(secondary_dose)))+
  geom_jitter(aes(color=as.factor(secondary_dose)), size=1.5)+
  #geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.5)+
  labs(x="Primary Treatment", y="ELISA OD", fill="Treatment")+
  stat_summary(aes(group=secondary_treatment), color="black", fun=mean, shape = "-", size=3)+
  geom_hline(yintercept=0.061)+
  scale_color_manual(values=c( "white", "gray75", "gray60", "gray50", "gray20"),
                    labels = c("0", "30", "100", "300", "7000"))+
  facet_wrap(~fct_rev(primary_treatment), labeller = as_labeller(treat_names), ncol=5)+
  theme_gray()

elisa.secondary
cv.secondary + elisa.secondary

######include secondary infection now####
s.ab <- m.ab %>%
  filter(elisa_od !=0 & dpi == 56)


#histograms of pid 56
ggplot(s.ab, aes(x=elisa_od, fill=primary_treatment))+
  #geom_density(color="black", alpha=0.5, aes(fill=primary_treatment))+ #fct_rev(primary_treatment) reverses order of overlay
  geom_histogram(aes(x=elisa_od, fill=primary_treatment), binwidth = 0.0075, position="identity", alpha=2, color="black", size=0.35)+
  geom_vline(xintercept=0.061, linetype="dashed", alpha=0.5)+
  xlim(c(0.04, 0.19))+
  labs(y= "Count", x= "ELISA OD", fill="Primary Treatment")+
  facet_grid(c("primary_treatment", "secondary_dose"),
             labeller = "label_both")+ #labeller = as_labeller(treat_names_n)
  theme_bw()

#density plots of pid 56
ggplot(s.ab, aes(x=elisa_od, fill=primary_treatment))+
  geom_density(color="black", alpha=0.5, aes(fill=primary_treatment))+ #fct_rev(primary_treatment) reverses order of overlay
  geom_histogram(aes(x=elisa_od, fill=primary_treatment), binwidth = 0.0025, position="identity", alpha=2, color="black", size=0.35)+
  geom_vline(xintercept=0.061, linetype="dashed", alpha=0.5)+
  xlim(c(0.04, 0.19))+
  labs(y= "Count", x= "ELISA OD", fill="Primary Treatment")+
  facet_grid(c("primary_treatment", "secondary_dose"),
             labeller = "label_both")+ #labeller = as_labeller(treat_names_n)
  theme_bw()

#All antibodies including dpi 56
g.ab <- ggplot(data = p.ab, aes(x = primary_treatment, y = elisa_od, color = primary_treatment)) +
  geom_jitter(size=0.5) +
  #scale_colour_manual(values=c( "red", "blue", "ivory4"))+
  #geom_line(aes(group = as.factor(band_number.x)), size = 0.5, alpha = 0.75) +
  stat_summary(aes(group=primary_treatment), fun=mean, geom="point", alpha=1, size=5, shape="+", color="black")+
  #stat_summary(aes(group=primary_treatment, color = primary_treatment), fun.y=mean,
  #             fun.min = function(x) mean(x)-sd(x),
  #             fun.max = function(x) mean(x)+sd(x),
               #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
               #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
  #             geom= "errorbar", size=0.75)+
  
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75) +
  ylab("MG Antibodies OD") +
  xlab("Inoculation Dose")
g.ab

#antibodies by primary treatment > DPI -8 to 41
g.ab.primary <- ggplot(data = m.ab %>% filter(elisa_od != 0), aes(x = primary_treatment, y = elisa_od, shape = primary_treatment)) +
  geom_jitter(size=1, width=0.25) +
  stat_summary(aes(group=primary_treatment), fun=mean, geom="point", alpha=1, size=5, shape="+", color="brown4")+
  #stat_summary(aes(group=primary_treatment, color = primary_treatment), fun.y=mean,
  #             fun.min = function(x) mean(x)-sd(x),
  #             fun.max = function(x) mean(x)+sd(x),
  #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
  #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
  #             geom= "errorbar", size=0.75)+
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75) +
  ylab("MG Antibodies OD") +
  xlab("Inoculation Dose")+
  facet_wrap(~dpi)
g.ab.primary
#continuous look at antibodies
ggplot(data = m.ab %>% filter(elisa_od != 0), aes(x = dpi, y = elisa_od, color = primary_treatment)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.5) +
  scale_colour_manual(values=c( "ivory4", "blue", "red"))+
  geom_line(aes(group = as.factor(band_number)), size = 0.5, alpha = 0.25) +
  stat_summary(aes(group=primary_treatment), fun=mean, geom="point", alpha=1, size=3)+
  stat_summary(aes(group=primary_treatment), fun=mean, geom="line", alpha=1, size=1)+
  #stat_summary(aes(group=primary_treatment, color = primary_treatment), fun.y=mean,
  #             fun.min = function(x) mean(x)-sd(x),
  #             fun.max = function(x) mean(x)+sd(x),
               #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
               #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
  #             geom= "errorbar", size=0.25, width = 0.25, alpha=1)+
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype="dashed") +
  labs(x="Days Post Inoculation", y="MG Antibodies OD", color= "Primary Treatment")+
  theme_minimal()

g.ab.all

 #new data frame with just dppi 56 elisa od and treatments

sid14 <- m.ab%>%
  select(dpi, bird_ID, band_number.x, elisa_od, primary_treatment, primary_dose, secondary_treatment, secondary_dose)
sid14


g.ab.sid14 <- ggplot(data = sid14 %>% filter(dpi=="56"), aes(x = primary_treatment, y = elisa_od, color = primary_treatment)) +
  #geom_point(size = 0.5) +
  geom_jitter(size=0.5) +
  scale_colour_manual(values=c( "red", "blue", "ivory4"))+
  stat_summary(aes(group=primary_treatment), fun=mean, geom="point", alpha=1, size=3, shape=3)+
  stat_summary(aes(group=primary_treatment, color = primary_treatment), fun.y=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
               #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
               geom= "errorbar", size=0.75)+
  
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75) +
  ylab("MG Antibodies OD PID 56 (SID 14)") +
  xlab("Inoculation Dose")+
  facet_wrap(~secondary_dose)

g.ab.sid14

