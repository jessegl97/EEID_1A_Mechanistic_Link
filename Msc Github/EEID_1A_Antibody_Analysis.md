EEID 1A Antibody Analysis
================
Jesse Garrett-Larsen
2024-01-26

Heyo here we go: Load packages

``` r
rm(list=ls())
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.3     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(dplyr)
#install.packages("glmmTMB")
library(glmmTMB)
```

    ## Warning in checkMatrixPackageVersion(): Package version inconsistency detected.
    ## TMB was built with Matrix version 1.6.0
    ## Current Matrix version is 1.5.4.1
    ## Please re-install 'TMB' from source using install.packages('TMB', type = 'source') or ask CRAN for a binary version of 'TMB' matching CRAN's 'Matrix' package

``` r
library(effects)
```

    ## Loading required package: carData
    ## lattice theme set by effectsTheme()
    ## See ?effectsTheme for details.

``` r
library(AICcmodavg)
#install.packages("emmeans")
library(emmeans)
library(DHARMa)
```

    ## This is DHARMa 0.4.6. For overview type '?DHARMa'. For recent changes, type news(package = 'DHARMa')

``` r
master <- read.csv("EEID2021_master_data_20220503.csv")
ab <- read.csv("EEID_1a_ELISA.csv")
m.ab <- left_join(master, ab, by = "bird_ID", keep=TRUE)
```

Format Data

``` r
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

m.ab$primary_treatment <- as.factor(m.ab$primary_treatment)
levels(m.ab$primary_treatment) <- c("High", "Low", "Sham")
m.ab$primary_treatment <- factor(m.ab$primary_treatment, levels = c("Sham", "Low", "High"))
m.ab$l_eye_score <- as.numeric(m.ab$l_eye_score)
```

    ## Warning: NAs introduced by coercion

``` r
m.ab$r_eye_score <- as.numeric(m.ab$r_eye_score)
```

    ## Warning: NAs introduced by coercion

Generate Threshold Cutoffs for Pathogen Load and Antibody Levels

``` r
m.ab$threshold_cutoff = 50
m.ab$seropos_cutoff = 0.061
```

Generate new infected and seroconversion columns based off of cutoffs

``` r
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
  mutate(infected_prim = ifelse(any(coalesce(dpi, 0) < 42 & coalesce(infected, 0) == 1), 1, 0)) %>%
  ungroup()

#generate infection data; if path load > 50 copies any time after secondary infection -> 1
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate (infected_sec = ifelse(any(coalesce(dpi, 0) > 42 & coalesce(infected, 0) == 1), 1, 0)) %>%
  ungroup()

#generate seropositivity data; if elisa OD > 0.061 = seropositive
#1 = bird seroconverted at any point during the priming phase of the experiment
m.ab <- m.ab %>%
  mutate(seropos = ifelse(elisa_od>seropos_cutoff, 1, 0))

m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(ever_seropos = ifelse(any(coalesce(dpi, 0) < 42 & coalesce(seropos, 0) == 1), 1, 0)) %>%
  ungroup()
```

``` r
#total eye score = sum of l and r eye score
m.ab$tes <- m.ab$l_eye_score + m.ab$r_eye_score

ggplot(m.ab %>% filter(dpi < 42), aes(x=tes, y=tes, color=band_number))+
  geom_jitter(height=0.1)+
  facet_wrap(~infected_prim)
```

    ## Warning: Removed 156 rows containing missing values (`geom_point()`).

![](EEID_1A_Antibody_Analysis_files/figure-gfm/Eye%20Score%20Formatting-1.png)<!-- -->

``` r
#2415 has eye score but no path load - super resistant?
```

``` r
#check for seropositive birds from quarantine from dataset
m.ab %>%
  filter(dpi == -8) %>%
  filter(elisa_od >= 0.061)%>%
  dplyr::select(band_number, elisa_od, primary_treatment, secondary_treatment, infected, seropos, dpi, ELISA.Run.Date)
```

    ## # A tibble: 1 × 8
    ##   band_number elisa_od primary_treatment secondary_treatment infected seropos
    ##         <int>    <dbl> <fct>                           <int>    <dbl>   <dbl>
    ## 1        2505    0.071 Low                                 3        0       1
    ## # ℹ 2 more variables: dpi <int>, ELISA.Run.Date <chr>

``` r
#omit 2505 from data set for analysis as it was positive at quarantine
m.ab <- m.ab %>%
  filter(band_number != 2505)

#which birds are seropositive but qPCR negative?
m.ab%>%
  filter(infected_prim == 0 & elisa_od > 0.061 & dpi < 42)%>%
  dplyr::select(band_number, dpi, tes, elisa_od, quantity, primary_treatment)
```

    ## # A tibble: 4 × 6
    ##   band_number   dpi   tes elisa_od quantity primary_treatment
    ##         <int> <int> <dbl>    <dbl>    <dbl> <fct>            
    ## 1        2398    14   0      0.062       NA Low              
    ## 2        2451    14   2.5    0.09        NA High             
    ## 3        2470    14   0      0.071       NA High             
    ## 4        2515    14   1.5    0.073       NA Low

``` r
#some sample days were split between two days. 
#One day experimental birds were measured (dpi -8 and dpi 14), the next day infected birds were measured (dpi -7, dpi15)
#this code makes a new column dpi.new that combines dpi 14 and 15
m.ab$dpi.new <- ifelse(m.ab$dpi %in% c(14, 15), "14", as.integer(m.ab$dpi))

#i then use dpi.new to replace the dpi column to combine dpi -7 and -8
m.ab$dpi <- ifelse(m.ab$dpi %in% c(-7, -8), "-8", as.integer(m.ab$dpi.new))

m.ab$dpi.new <- NULL #delete dpi.new column

m.ab$dpi <- as.numeric(m.ab$dpi)

m.ab <- m.ab %>%
  filter(experiment_location == "vt") #include only birds at VT
```

General overview of antibody responses

``` r
ggplot(m.ab, aes(x=dpi, y=elisa_od, color=fct_rev(primary_treatment)))+
  geom_point()+
  geom_hline(yintercept = 0.061)+
  geom_smooth(aes(x=dpi, y=elisa_od, group=band_number), alpha=0.5, size=0.5)+
  #scale_color_manual(values = c("red", "green4","blue"))+
  labs(x="Days Post Infection", y="ELISA OD", color="Primary Treatment")+
  facet_wrap(~primary_treatment)
```

    ## `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

![](EEID_1A_Antibody_Analysis_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
\####Does Inoculation Dose Predict Antibody Levels?####

``` r
#data frame with only primary infection (no PID 56)
p.ab <- m.ab %>%
  filter(dpi <=41)
```

I use a gamma distribution for the antibody data With models spanning
multiple days, band_number is included as a random effect

Do I use a log link function or inverse link function? Compare AIC and
logLik Skipping ahead here with models, but showing that *inverse link
works better*

``` r
#pre-infection
lm0a <- glm(elisa_od~primary_treatment, data=p.ab %>% filter(dpi <0), family=Gamma())
lm0b <- glm(elisa_od~primary_treatment, data=p.ab %>% filter(dpi <0), family=Gamma(log))

aictab(cand.set=list(lm0a, lm0b), modnames=c("inv","log"))
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##     K     AICc Delta_AICc AICcWt Cum.Wt     LL
    ## log 4 -1357.07          0    0.5    0.5 682.67
    ## inv 4 -1357.07          0    0.5    1.0 682.67

``` r
#no difference between the two

#Across all infection
ps1 <- glmmTMB(elisa_od~primary_treatment + (1|band_number), data=p.ab, family=Gamma())
ps2 <- glmmTMB(elisa_od~primary_treatment + (1|band_number), data=p.ab, family=Gamma(log))

aictab(cand.set=list(ps1, ps2), modnames=c("inv","log"))
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##     K     AICc Delta_AICc AICcWt Cum.Wt      LL
    ## inv 5 -2548.74       0.00   0.97   0.97 1279.44
    ## log 5 -2541.45       7.29   0.03   1.00 1275.79

``` r
#inverse link function is slightly better

#DPI 14
p1a <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi ==14), family=Gamma())
p2a <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi ==14), family=Gamma(log))

aictab(cand.set=list(p1a, p2a), modnames=c("inv","log"))
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##     K    AICc Delta_AICc AICcWt Cum.Wt     LL
    ## log 4 -719.09          0    0.5    0.5 363.69
    ## inv 4 -719.09          0    0.5    1.0 363.69

``` r
#no difference 

#DPI 41
p1b <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi ==41), family=Gamma())
p2b <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi ==41), family=Gamma(log))

aictab(cand.set=list(p1b, p2b), modnames=c("inv","log"))
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##     K     AICc Delta_AICc AICcWt Cum.Wt     LL
    ## log 4 -1044.25          0    0.5    0.5 526.26
    ## inv 4 -1044.25          0    0.5    1.0 526.26

``` r
#no difference
```

Baseline antibody levels were not significantly different from each
other before inoculation.

emmeans shows estimated marginal means for variables and allows for
post-hoc pairwise comparisons using the Tukey-Kramer method to minimize
type 1 error.

``` r
lm0 <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi <0), family=Gamma())
summary(lm0)
```

    ##  Family: Gamma  ( inverse )
    ## Formula:          elisa_od ~ primary_treatment
    ## Data: p.ab %>% filter(dpi < 0)
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  -1357.3  -1345.2    682.7  -1365.3      151 
    ## 
    ## 
    ## Dispersion estimate for Gamma family (sigma^2): 0.00405 
    ## 
    ## Conditional model:
    ##                       Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)           21.54626    0.19209  112.17   <2e-16 ***
    ## primary_treatmentLow  -0.18053    0.27052   -0.67    0.505    
    ## primary_treatmentHigh  0.02476    0.26923    0.09    0.927    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(allEffects(lm0))
```

![](EEID_1A_Antibody_Analysis_files/figure-gfm/Baseline%20antibody%20levels%20were%20not%20significantly%20different-1.png)<!-- -->

``` r
emm_r <- emmeans(lm0, ~primary_treatment)
pairs(emm_r)
```

    ##  contrast    estimate    SE  df z.ratio p.value
    ##  Sham - Low    0.1805 0.271 Inf   0.667  0.7825
    ##  Sham - High  -0.0248 0.269 Inf  -0.092  0.9953
    ##  Low - High   -0.2053 0.268 Inf  -0.766  0.7241
    ## 
    ## Note: contrasts are still on the inverse scale 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

Using primary_treatment (categorical), I show that birds inoculated with
MG produce antibodies in a dose-dependent manner.

Across all days primary infection, inoculation dose

``` r
lm1 <- glmmTMB(elisa_od~primary_treatment + (1|band_number), data=p.ab, family=Gamma())
summary(lm1)
```

    ##  Family: Gamma  ( inverse )
    ## Formula:          elisa_od ~ primary_treatment + (1 | band_number)
    ## Data: p.ab
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  -2548.9  -2528.3   1279.4  -2558.9      444 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups      Name        Variance Std.Dev.
    ##  band_number (Intercept) 2.571    1.604   
    ## Number of obs: 449, groups:  band_number, 155
    ## 
    ## Dispersion estimate for Gamma family (sigma^2): 0.0631 
    ## 
    ## Conditional model:
    ##                       Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)            22.1212     0.5214   42.43  < 2e-16 ***
    ## primary_treatmentLow   -2.6771     0.6910   -3.87 0.000107 ***
    ## primary_treatmentHigh  -6.8396     0.6414  -10.66  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
l1r <- simulateResiduals(lm1)
plot(l1r)
```

![](EEID_1A_Antibody_Analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#model comparison 
p1 <- glmmTMB(elisa_od~primary_treatment + (1|band_number), data=p.ab, family=Gamma())
p2<- glmmTMB(elisa_od~primary_treatment + sex + (1|band_number), data=p.ab, family=Gamma())
p3<- glmmTMB(elisa_od~primary_treatment * sex + (1|band_number), data=p.ab, family=Gamma())
p4 <- glmmTMB(elisa_od~1 + (1|band_number), data=p.ab, family=Gamma())
p5 <- glmmTMB(elisa_od~primary_treatment + mass + (1|band_number), data=p.ab, family=Gamma())

aictab(cand.set=list(p1, p2, p3, p4, p5), modnames=c("p1","p2", "p3", "p4", "p5"))
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##     K     AICc Delta_AICc AICcWt Cum.Wt      LL
    ## p1  5 -2548.74       0.00   0.55   0.55 1279.44
    ## p2  6 -2547.86       0.88   0.35   0.90 1280.03
    ## p3  8 -2545.39       3.36   0.10   1.00 1280.86
    ## p5 66 -2512.91      35.83   0.00   1.00 1334.03
    ## p4  3 -2458.05      90.69   0.00   1.00 1232.05

``` r
emm_results <- emmeans(lm1, ~ primary_treatment, scale="response")
pairs(emm_results)
```

    ##  contrast    estimate    SE  df z.ratio p.value
    ##  Sham - Low      2.68 0.691 Inf   3.874  0.0003
    ##  Sham - High     6.84 0.641 Inf  10.663  <.0001
    ##  Low - High      4.16 0.587 Inf   7.089  <.0001
    ## 
    ## Note: contrasts are still on the inverse scale 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

``` r
#antibodies across all days primary infection
ggplot(data=p.ab, aes(x=dpi, y= elisa_od, shape=primary_treatment,
                      color=fct_rev(primary_treatment)))+
  geom_jitter(width=5, size=1.5)+
  geom_line(aes(x=dpi, y=elisa_od, group = as.factor(band_number)))+
  geom_hline(yintercept = 0.061, linetype="dashed")+
  stat_summary(aes(group=primary_treatment), fun=mean, geom="line")+
  stat_summary(aes(group=primary_treatment, shape = primary_treatment), fun.y=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               geom= "errorbar", size=.75, width=2)+
  labs(title = "Antibody Levels Primary Infection", y="Antibody Levels", 
       x= "Primary Treatment", shape="Primary Treatment", color="Primary Treatment")
```

    ## Warning: The `fun.y` argument of `stat_summary()` is deprecated as of ggplot2 3.3.0.
    ## ℹ Please use the `fun` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning in stat_summary(aes(group = primary_treatment, shape =
    ## primary_treatment), : Ignoring unknown aesthetics: shape

    ## Warning: Removed 637 rows containing non-finite values (`stat_summary()`).
    ## Removed 637 rows containing non-finite values (`stat_summary()`).

    ## Warning: Removed 637 rows containing missing values (`geom_point()`).

    ## Warning: Removed 4 rows containing missing values (`geom_line()`).

![](EEID_1A_Antibody_Analysis_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
#Final Model
lm2 <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi == 14), family=Gamma())
summary(lm2)
```

    ##  Family: Gamma  ( inverse )
    ## Formula:          elisa_od ~ primary_treatment
    ## Data: p.ab %>% filter(dpi == 14)
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   -719.4   -707.6    363.7   -727.4      136 
    ## 
    ## 
    ## Dispersion estimate for Gamma family (sigma^2): 0.078 
    ## 
    ## Conditional model:
    ##                       Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)             22.424      1.030  21.776  < 2e-16 ***
    ## primary_treatmentLow    -5.928      1.219  -4.864 1.15e-06 ***
    ## primary_treatmentHigh  -11.879      1.106 -10.736  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emm_results <- emmeans(lm2, ~ primary_treatment, scale="response")
emm_results
```

    ##  primary_treatment emmean    SE  df asymp.LCL asymp.UCL
    ##  Sham                22.4 1.030 Inf     20.41      24.4
    ##  Low                 16.5 0.652 Inf     15.22      17.8
    ##  High                10.5 0.405 Inf      9.75      11.3
    ## 
    ## Results are given on the inverse (not the response) scale. 
    ## Confidence level used: 0.95

``` r
pairs(emm_results)
```

    ##  contrast    estimate    SE  df z.ratio p.value
    ##  Sham - Low      5.93 1.219 Inf   4.864  <.0001
    ##  Sham - High    11.88 1.106 Inf  10.736  <.0001
    ##  Low - High      5.95 0.767 Inf   7.758  <.0001
    ## 
    ## Note: contrasts are still on the inverse scale 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

``` r
plot(allEffects(lm2))
```

![](EEID_1A_Antibody_Analysis_files/figure-gfm/Antibody%20levels%20on%20day%2014%20post%20infection%20were%20sig%20different%20between%20groups-1.png)<!-- -->

``` r
#model comparison
pr1 <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi == 14), family=Gamma())
pr2<- glmmTMB(elisa_od~primary_treatment + sex , data=p.ab %>% filter(dpi == 14), family=Gamma())
pr3<- glmmTMB(elisa_od~primary_treatment * sex , data=p.ab %>% filter(dpi == 14), family=Gamma())
pr4 <- glmmTMB(elisa_od~primary_treatment * room , data=p.ab %>% filter(dpi == 14), family=Gamma())
pr5 <- glmmTMB(elisa_od~1, data=p.ab %>% filter(dpi == 14), family=Gamma())
pr6 <- glmmTMB(elisa_od~primary_treatment + mass, data=p.ab %>% filter(dpi == 14), family=Gamma())
pr7 <- glmmTMB(elisa_od~primary_treatment + room , data=p.ab %>% filter(dpi == 14), family=Gamma())

#AICc
aictab(cand.set=list(pr1, pr2, pr3, pr4, pr5, pr6, pr7), modnames=c("pr1","pr2", "pr3", "pr4", "pr5", "pr6", "pr7"))
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##      K    AICc Delta_AICc AICcWt Cum.Wt     LL
    ## pr1  4 -719.09       0.00   0.39   0.39 363.69
    ## pr7 10 -718.99       0.09   0.37   0.75 370.35
    ## pr2  5 -717.36       1.73   0.16   0.92 363.90
    ## pr3  7 -715.99       3.09   0.08   1.00 365.42
    ## pr4 22 -694.97      24.12   0.00   1.00 373.81
    ## pr6 50 -636.94      82.15   0.00   1.00 397.12
    ## pr5  2 -612.99     106.09   0.00   1.00 308.54

``` r
#There were 15 blood samples that were lost before they were run and are not included in this graph
missing_samples <- p.ab %>%
  filter(dpi %in% c(-8, 14, 41) & is.na(elisa_od)) %>%
  dplyr::select(band_number, bird_ID, dpi, elisa_od)
print(missing_samples)
```

    ## # A tibble: 16 × 4
    ##    band_number bird_ID   dpi elisa_od
    ##          <int> <chr>   <dbl>    <dbl>
    ##  1        2496 14_2496    14       NA
    ##  2        2385 15_2385    14       NA
    ##  3        2395 15_2395    14       NA
    ##  4        2407 15_2407    14       NA
    ##  5        2417 15_2417    14       NA
    ##  6        2439 15_2439    14       NA
    ##  7        2450 15_2450    14       NA
    ##  8        2460 15_2460    14       NA
    ##  9        2489 15_2489    14       NA
    ## 10        2498 15_2498    14       NA
    ## 11        2511 15_2511    14       NA
    ## 12        2537 15_2537    14       NA
    ## 13        2546 15_2546    14       NA
    ## 14        2563 15_2563    14       NA
    ## 15        2565 15_2565    14       NA
    ## 16        2375 41_2375    41       NA

``` r
#Model predictions
dat.new=expand.grid(primary_treatment=unique(p.ab$primary_treatment))#new grid to put predictions into
dat.new$yhat = predict(lm2, type="response", newdata=dat.new, re.form=NA) #predicted values based off lm2 model
head(dat.new)
```

    ##   primary_treatment       yhat
    ## 1              Sham 0.04459459
    ## 2              High 0.09483017
    ## 3               Low 0.06061996

``` r
#graph showing antibody levels on day 14 post-infection with model predictions
pid14.pred <- ggplot(data=p.ab %>% filter(dpi == 14), aes(x=primary_treatment, y=elisa_od, shape=primary_treatment))+
  geom_jitter(size=2, width=0.25)+
  geom_point(data=dat.new, aes(x=primary_treatment, y=yhat, shape=primary_treatment), size=10, shape="-", color="brown")+#model predictions
  scale_shape_manual(values = c(16, 17, 15))+
  geom_hline(yintercept=0.061, linetype="dashed", alpha=0.5)+
  labs(title="MG Antibodies Day 14", x="Primary Treatment", y="Antibody Levels", shape="Primary Treatment")

pid14.pred
```

    ## Warning: Removed 15 rows containing missing values (`geom_point()`).

![](EEID_1A_Antibody_Analysis_files/figure-gfm/Antibody%20levels%20on%20day%2014%20post%20infection%20were%20sig%20different%20between%20groups-2.png)<!-- -->
I also looked at primary_dose, a continuous variable to test whether
inoculation dose affects antibody response.

``` r
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
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##         K    AICc Delta_AICc AICcWt Cum.Wt     LL
    ## prd1    3 -699.47       0.00   0.39   0.39 352.82
    ## prd3    9 -698.76       0.71   0.28   0.67 359.07
    ## prd2.5  5 -697.74       1.72   0.17   0.83 354.10
    ## prd2    4 -697.74       1.73   0.17   1.00 353.02
    ## prd4   15 -685.96      13.51   0.00   1.00 359.92
    ## prd6   49 -617.53      81.94   0.00   1.00 384.98
    ## prd5    2 -612.94      86.53   0.00   1.00 308.51

``` r
summary(prd1)
```

    ## 
    ## Call:
    ## glm(formula = elisa_od ~ primary_dose, family = Gamma(), data = p.ab %>% 
    ##     filter(dpi == 14))
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.874e+01  7.291e-01  25.699  < 2e-16 ***
    ## primary_dose -2.737e-04  3.001e-05  -9.119 8.04e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.1295845)
    ## 
    ##     Null deviance: 23.968  on 139  degrees of freedom
    ## Residual deviance: 12.896  on 138  degrees of freedom
    ##   (15 observations deleted due to missingness)
    ## AIC: -699.64
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
#continuous model predictions
dat.new2 <- data.frame(primary_dose = seq(min(p.ab$primary_dose), max(p.ab$primary_dose), length.out = 10))
dat.new2$yhat <- predict(prd1, type = "response", interval = "confidence", level = 0.95, newdata = dat.new2)

#plot continuous predicted values over raw data
pid14.dose.pred.c <- ggplot(data=p.ab %>% filter(dpi == 14), aes(x=primary_dose, y=elisa_od))+
  geom_jitter(size=2, width=0.5, shape = 1)+
  geom_line(data=dat.new2, aes(x=primary_dose, y=yhat), color="brown")+#model predictions
  labs(x="Primary Dose", y="Antibody Levels", shape="Primary Dose")

pid14.dose.pred.c
```

    ## Warning: Removed 15 rows containing missing values (`geom_point()`).

![](EEID_1A_Antibody_Analysis_files/figure-gfm/Primary%20Dose%20Continuous%20Day%2014-1.png)<!-- -->
Alternatively, I can subset the primary_dose into facors based on
primary_dose. This maintains the continuous nature of the data, but
allows for pairwise comparisons using emmeans.

``` r
# Create a new factor variable based on primary_dose
p.ab <- p.ab %>% mutate(primary_dose_group = cut(primary_dose, breaks = c(-Inf, 100, 2000, Inf), labels = c("Sham", "Low", "High")))

# Fit glm with  new factor variable
prdc1 <- glm(elisa_od ~ primary_dose_group, data = p.ab %>% filter(dpi == 14), family = Gamma())

# Calculate estimated marginal means (EMMs)
emm_results <- emmeans(prdc1, specs = "primary_dose_group")

# Conduct pairwise comparisons between the groups
pairs(emm_results)
```

    ##  contrast    estimate    SE  df t.ratio p.value
    ##  Sham - Low      5.93 1.393 137   4.257  0.0001
    ##  Sham - High    11.88 1.264 137   9.396  <.0001
    ##  Low - High      5.95 0.877 137   6.790  <.0001
    ## 
    ## Note: contrasts are still on the inverse scale 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

Antibody levels return to baseline in the low (p = 0.149), but not high
dose (p \< 0.001) groups by day 41.

``` r
lm3 <- glm(elisa_od~primary_treatment, data=p.ab %>% filter(dpi==41), family=Gamma())
summary(lm3)
```

    ## 
    ## Call:
    ## glm(formula = elisa_od ~ primary_treatment, family = Gamma(), 
    ##     data = p.ab %>% filter(dpi == 41))
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            22.2319     0.5617  39.577  < 2e-16 ***
    ## primary_treatmentLow   -1.4505     0.7725  -1.878   0.0624 .  
    ## primary_treatmentHigh  -4.8947     0.7073  -6.921  1.2e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.03256064)
    ## 
    ##     Null deviance: 5.7018  on 153  degrees of freedom
    ## Residual deviance: 3.9511  on 151  degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## AIC: -1044.5
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
emm_results<- emmeans(lm3, ~ primary_treatment)
pairs(emm_results)
```

    ##  contrast    estimate    SE  df t.ratio p.value
    ##  Sham - Low      1.45 0.773 151   1.878  0.1487
    ##  Sham - High     4.89 0.707 151   6.921  <.0001
    ##  Low - High      3.44 0.683 151   5.046  <.0001
    ## 
    ## Note: contrasts are still on the inverse scale 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

``` r
#model comparison
p1 <- glm(elisa_od~primary_treatment, data=p.ab %>% filter(dpi==41), family=Gamma())
p2<- glm(elisa_od~primary_treatment + sex , data=p.ab %>% filter(dpi==41), family=Gamma())
p3<- glm(elisa_od~primary_treatment * sex , data=p.ab %>% filter(dpi==41), family=Gamma())
p4 <- glm(elisa_od~primary_treatment + quantity , data=p.ab %>% filter(dpi==41), family=Gamma())
p5 <- glm(elisa_od~1, data=p.ab %>% filter(dpi==41), family=Gamma())
p6 <- glm(elisa_od~primary_treatment + mass, data= p.ab %>% filter(dpi==41), family=Gamma())

#AICc
aictab(cand.set=list(p1, p2, p3, p4, p5, p6), modnames=c("p1", "p2", "p3", "p4", "p5", "p6"))
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##     K     AICc Delta_AICc AICcWt Cum.Wt     LL
    ## p2  5 -1044.86       0.00   0.46   0.46 527.63
    ## p1  4 -1044.25       0.60   0.34   0.81 526.26
    ## p4  5 -1042.12       2.74   0.12   0.93 526.26
    ## p3  7 -1041.21       3.65   0.07   1.00 527.99
    ## p5  2  -991.66      53.19   0.00   1.00 497.87
    ## p6 50  -975.02      69.84   0.00   1.00 562.27

``` r
summary(p2)
```

    ## 
    ## Call:
    ## glm(formula = elisa_od ~ primary_treatment + sex, family = Gamma(), 
    ##     data = p.ab %>% filter(dpi == 41))
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            21.8125     0.6256  34.869  < 2e-16 ***
    ## primary_treatmentLow   -1.4418     0.7689  -1.875   0.0627 .  
    ## primary_treatmentHigh  -4.8921     0.7039  -6.950 1.04e-10 ***
    ## sexm                    0.8383     0.5715   1.467   0.1445    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.0322655)
    ## 
    ##     Null deviance: 5.7018  on 153  degrees of freedom
    ## Residual deviance: 3.8817  on 150  degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## AIC: -1045.3
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
#p2 is the best model
#but sex adds complexity and is not significant (0.1445) so remove

#p1 is second best - use
summary(p1)
```

    ## 
    ## Call:
    ## glm(formula = elisa_od ~ primary_treatment, family = Gamma(), 
    ##     data = p.ab %>% filter(dpi == 41))
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            22.2319     0.5617  39.577  < 2e-16 ***
    ## primary_treatmentLow   -1.4505     0.7725  -1.878   0.0624 .  
    ## primary_treatmentHigh  -4.8947     0.7073  -6.921  1.2e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.03256064)
    ## 
    ##     Null deviance: 5.7018  on 153  degrees of freedom
    ## Residual deviance: 3.9511  on 151  degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## AIC: -1044.5
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
#model predictions
dat.newer=expand.grid(primary_treatment=unique(p.ab$primary_treatment))#new grid to put predictions into
dat.newer$yhat = predict(lm3, type="response", newdata=dat.newer, re.form=NA) #predicted values based off lm3  model
head(dat.newer)
```

    ##   primary_treatment       yhat
    ## 1              Sham 0.04498039
    ## 2              High 0.05767925
    ## 3               Low 0.04812000

``` r
#plot predicted values over raw data
pid41.pred <- ggplot(data = m.ab %>% filter(dpi == 41), aes(x = primary_treatment, y = elisa_od, shape = primary_treatment)) +
  geom_jitter(size=1.5, width = 0.25) +
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype='dashed') +
  labs(x="Inoculation Dose", y= "MG Antibodies OD", shape = "Primary Treatment")+
  geom_point(data=dat.newer, aes(x=primary_treatment, y=yhat, shape=primary_treatment), size=10, shape="-", color="brown4")#model predictions

pid41.pred
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](EEID_1A_Antibody_Analysis_files/figure-gfm/Antibody%20levels%20on%20day%2041%20post%20infection%20were%20sig%20different%20between%20groups-1.png)<!-- -->

Antibody levels return to baseline in the low (p = 0.149), but not high
dose (p \< 0.001) groups by day 41.

``` r
#model comparison
pr1 <- glm(elisa_od~primary_dose, data=p.ab %>% filter(dpi==41), family=Gamma())
pr2<- glm(elisa_od~primary_dose + sex , data=p.ab %>% filter(dpi==41), family=Gamma())
pr3<- glm(elisa_od~primary_dose * sex , data=p.ab %>% filter(dpi==41), family=Gamma())
pr4 <- glm(elisa_od~primary_dose + quantity , data=p.ab %>% filter(dpi==41), family=Gamma())
pr5 <- glm(elisa_od~1, data=p.ab %>% filter(dpi==41), family=Gamma())
pr6 <- glm(elisa_od~primary_dose + mass, data= p.ab %>% filter(dpi==41), family=Gamma())


#AICc
aictab(cand.set=list(pr1, pr2, pr3, pr4, pr5, pr6), modnames=c("pr1", "pr2", "pr3", "pr4", "pr5", "pr6"))
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##      K     AICc Delta_AICc AICcWt Cum.Wt     LL
    ## pr2  4 -1043.15       0.00   0.42   0.42 525.71
    ## pr1  3 -1042.54       0.62   0.31   0.74 524.35
    ## pr3  5 -1041.13       2.02   0.15   0.89 525.77
    ## pr4  4 -1040.44       2.72   0.11   1.00 524.35
    ## pr5  2  -991.66      51.49   0.00   1.00 497.87
    ## pr6 49  -976.76      66.39   0.00   1.00 560.94

``` r
#pr2 best model
#but sex adds complexity and is not significant (0.146) so remove

summary(pr1)
```

    ## 
    ## Call:
    ## glm(formula = elisa_od ~ primary_dose, family = Gamma(), data = p.ab %>% 
    ##     filter(dpi == 41))
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   2.155e+01  3.956e-01  54.471  < 2e-16 ***
    ## primary_dose -1.408e-04  1.972e-05  -7.138 3.64e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.03335949)
    ## 
    ##     Null deviance: 5.7018  on 153  degrees of freedom
    ## Residual deviance: 4.0499  on 152  degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## AIC: -1042.7
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
#continuous model predictions
dat.newer2 <- data.frame(primary_dose = seq(min(p.ab$primary_dose), max(p.ab$primary_dose), length.out = 100))
dat.newer2$yhat <- predict(pr1, type = "response", interval = "confidence", level = 0.95, newdata = dat.newer2)
head(dat.newer2)
```

    ##   primary_dose       yhat
    ## 1       0.0000 0.04640312
    ## 2     303.0303 0.04649517
    ## 3     606.0606 0.04658758
    ## 4     909.0909 0.04668036
    ## 5    1212.1212 0.04677351
    ## 6    1515.1515 0.04686703

``` r
pid41.dose.pred.c <- ggplot(data=p.ab %>% filter(dpi == 41), aes(x=primary_dose, y=elisa_od))+
  geom_jitter(size=2, width=0.1)+
  geom_line(data=dat.newer2, aes(x=primary_dose, y=yhat), color="brown")+#model predictions
  labs(title = "ELISA OD PID 41", x="Primary Dose", y="ELISA OD", shape="Primary Dose")+
  geom_hline(yintercept = 0.061, color = "black", alpha = 0.75, linetype='dashed')

pid41.dose.pred.c
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](EEID_1A_Antibody_Analysis_files/figure-gfm/Primary%20Dose%20Continuous%20Day%2041-1.png)<!-- -->

``` r
#pairwise comparisons
p.ab <- p.ab %>% mutate(primary_dose_group = cut(primary_dose, breaks = c(-Inf, 100, 2000, Inf), labels = c("Sham", "Low", "High")))

# Fit glm with  new factor variable
pr1 <- glm(elisa_od~primary_dose_group, data=p.ab %>% filter(dpi==41), family=Gamma())

# Calculate estimated marginal means (EMMs)
emm_results <- emmeans(pr1, specs = "primary_dose_group")

# Conduct pairwise comparisons between the groups
pairs(emm_results)
```

    ##  contrast    estimate    SE  df t.ratio p.value
    ##  Sham - Low      1.45 0.773 151   1.878  0.1487
    ##  Sham - High     4.89 0.707 151   6.921  <.0001
    ##  Low - High      3.44 0.683 151   5.046  <.0001
    ## 
    ## Note: contrasts are still on the inverse scale 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

*Summary: *
