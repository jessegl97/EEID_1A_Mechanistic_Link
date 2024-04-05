###Extra EEID 1A Antibody Code
```{r Proportion Seroconverted  vs Infected day 14}
#what proportion of birds seroconverted by day 14?
library(gt)

library(gtExtras)
#proportion seropositive
m.abs<- m.ab %>%
  filter(dpi %in% c(-8, 7, 14, 41))%>%
  group_by(dpi, primary_treatment)%>%
  summarise(prop_seropos = sum(seropos == 1, na.rm=TRUE) / n(),
            n_seropos = sum(seropos == 1, na.rm=TRUE),
            n_seroneg = sum(seropos == 0, na.rm=TRUE),
            prop_inf = sum(inf == 1, na.rm=TRUE)/ n(),
            n_inf = sum(inf == 1, na.rm=TRUE),
            n_not_inf = sum(inf == 0),
            seropos_prim = sum(seropos_prim == 1, na.rm=TRUE),
            inf_prim = sum(inf_prim == 1, na.rm=TRUE),
            total_birds =n())
m.abs

table_gt <- m.abs %>%
  gt() %>%
  fmt_number(
    columns=vars(prop_seropos, prop_inf),
    decimals = 2
  )%>%
  tab_header(
    title = "Summary Table",
    subtitle = "Summary of Proportion Seropositive and Infected by Treatment Group"
  ) %>%
  cols_label(
    dpi = "Day Post Infection",
    primary_treatment = "Primary Treatment",
    prop_seropos = "Proportion Seropositive",
    n_seropos = "Number Seropositive",
    n_seroneg = "Number Seronegative",
    prop_inf = "Proportion Infected",
    n_inf = "Number Infected",
    n_not_inf = "Number Not Infected",
    seropos_prim = "Seropositive Primary",
    inf_prim = "Infected Primary",
    total_birds = "Total Birds",
  )

# Display the table
table_gt


#proportion seropositive
ggplot(m.abs, aes(x=dpi, y=prop_seropos, color=primary_treatment, shape=primary_treatment))+
  geom_point(size = 1)+
  scale_color_manual(values = pri_colors)+
  ylim(c(0,1))+
  geom_text(aes(label = round(prop_seropos, 2)), vjust = -0.25, size = 3) +
  labs(title= "Proportion Seropositive", x = "Day Post Infection", y = "Proportion Seropositive", color = "Primary Treatment", shape = "Primary Treatment")#+
guides(shape = FALSE)

m.ab.p<- m.ab %>%
  filter(dpi %in% c(-8, 7, 14, 41)) %>%
  group_by(dpi, primary_treatment)%>%
  summarise(prop_seropos = sum(seropos == 1, na.rm=TRUE) / n(), na.rm=TRUE,
            n_seropos = sum(seropos == 1, na.rm=TRUE),
            n_seroneg = sum(seropos == 0, na.rm=TRUE),
            prop_inf = sum(inf == 1, na.rm=TRUE)/ n(), na.rm=TRUE,
            n_inf = sum(inf == 1, na.rm=TRUE),
            n_not_inf = sum(inf == 0, na.rm=TRUE),
            seropos_prim = sum(seropos_prim == 1, na.rm=TRUE),
            inf_prim = sum(inf_prim == 1, na.rm=TRUE),
            total_birds =n(), na.rm=TRUE)
m.ab.p

#proportion infected primary
ggplot(m.ab.p, aes(x=dpi, y=prop_inf, color=primary_treatment, shape=primary_treatment))+
  geom_point(size = 1)+
  scale_color_manual(values = pri_colors)+
  ylim(c(0,1))+
  #geom_text(aes(label = paste(n_inf, "/", total_birds)), vjust = -0.2, size = 3) +
  geom_text(aes(label = round(prop_inf, 2)), vjust = -0.2, size = 3) +
  labs(title = "Proportion Infected Primary", x = "Day Post Infection", y = "Proportion Infected", color = "Primary Treatment", shape = "Primary Treatment")

#proportion seropositive primary
ggplot(m.abs %>% filter(dpi %in% c(-8, 14, 41)), aes(x=dpi, y=prop_seropos, color=primary_treatment, shape=primary_treatment))+
  geom_point(size = 1)+
  scale_color_manual(values = pri_colors)+
  ylim(c(0,1))+
  geom_text(aes(label = round(prop_seropos, 2)), vjust = -0.2, size = 3) +
  labs(title = "Proportion Seropositive", x = "Day Post Infection", y = "Proportion Seropositive", color = "Primary Treatment", shape = "Primary Treatment")


m.ab.inf<- m.ab %>%
  filter(dpi %in% c(-8, 7, 14, 41, 46, 49, 56, 63)) %>%
  group_by(dpi, primary_treatment, secondary_dose)%>%
  summarise(prop_seropos = sum(seropos == 1, na.rm=TRUE) / n(), na.rm=TRUE,
            n_seropos = sum(seropos == 1, na.rm=TRUE),
            n_seroneg = sum(seropos == 0, na.rm=TRUE),
            prop_inf = sum(inf == 1, na.rm=TRUE)/ n(), na.rm=TRUE,
            n_inf = sum(inf == 1, na.rm=TRUE),
            n_not_inf = sum(inf == 0, na.rm=TRUE),
            seropos_prim = sum(seropos_prim == 1, na.rm=TRUE),
            inf_prim = sum(inf_prim == 1, na.rm=TRUE),
            total_birds =n(), na.rm=TRUE)

#proportion infected
ggplot(m.ab.inf %>% filter(dpi > 41), aes(x=dpi, y=prop_inf, color=as.factor(secondary_dose), shape=primary_treatment, groups=secondary_dose))+
  geom_point(size = 1)+
  scale_color_manual(values = sec_colors)+
  geom_path(aes(linetype=fct_rev(primary_treatment)))+
  ylim(c(0,1))+
  geom_text(aes(label = round(prop_inf, 2)), vjust = -0.2, size = 3) +
  labs(title = "Proportion Infected", x = "Day Post Infection", y = "Proportion Infected", color = "Primary Treatment", shape = "Primary Treatment")+
  facet_grid(~primary_treatment)

#infected secondary
m.ab.sec<- m.ab %>%
  filter(dpi %in% c(46, 49, 56, 63)) %>%
  group_by(primary_treatment, secondary_dose)%>%
  summarise(prop_seropos = sum(seropos == 1, na.rm=TRUE) / n(), na.rm=TRUE,
            n_seropos = sum(seropos == 1, na.rm=TRUE),
            n_seroneg = sum(seropos == 0, na.rm=TRUE),
            prop_inf = sum(inf == 1, na.rm=TRUE)/ n(), na.rm=TRUE,
            n_inf = sum(inf == 1, na.rm=TRUE),
            n_not_inf = sum(inf == 0, na.rm=TRUE),
            seropos_prim = sum(seropos_prim == 1, na.rm=TRUE),
            inf_prim = sum(inf_prim == 1, na.rm=TRUE),
            total_birds =n(), na.rm=TRUE)

#proportion infected secondary
ggplot(m.ab.sec, aes(x=log(secondary_dose), y=prop_inf, color=as.factor(secondary_dose), shape=primary_treatment, groups=secondary_dose))+
  geom_point(size = 1)+
  scale_color_manual(values = sec_colors)+
  ylim(c(0,1))+
  geom_text(aes(label = round(prop_inf, 2)), vjust = -0.2, size = 3) +
  labs(title = "Proportion Infected", x = "log10(secondary dose)", y = "Proportion Infected", color = "Primary Treatment", shape = "Primary Treatment")+
  facet_grid(~primary_treatment)
```
```{r If birds seroconverted, were they protected from secondary pathology?}
#does seropositivity during primary infection predict whether or not a bird will develop pathology upon secondary infection?
lm.sp1 <- glm(inf_sec ~ elisa_od + primary_treatment*secondary_dose, data=m.ab %>% filter(dpi == 14), family=binomial(link= "logit"))
summary(lm.sp1)
plot(allEffects(lm.sp1))
car::Anova(lm.sp1, type="III")
emm_results<- emmeans(lm.sp1, ~ primary_treatment*secondary_dose)
pairs(emm_results)
simulateResiduals(lm.sp1, plot=T)
dat.new=expand.grid(secondary_dose=unique(m.ab$secondary_dose),
                    primary_treatment=unique(m.ab$primary_treatment),
                    diseased_sec=unique(m.ab$diseased_sec),
                    seropos_prim=unique(m.ab$seropos_prim),
                    elisa_od = unique(m.ab$elisa_od))#new grid to put predictions into
dat.new$yhat = predict(lm.sp1, type="response", newdata=dat.new, re.form=NA) #predicted values
head(dat.new)

ggplot(m.ab, aes(y=diseased_sec, x=seropos_prim))+
  geom_jitter(size=1, height=0, width=0.5, aes(shape=primary_treatment, color=primary_treatment))+
  geom_point(data = dat.new, aes(x =seropos_prim, y = yhat, shape = primary_treatment, color = primary_treatment)) +
  geom_path(data = dat.new, aes(x = seropos_prim, y = yhat, color = primary_treatment, linetype = primary_treatment))+
  labs(x="Seroconversion Status During Primary Infection", y="Disease During Secondary Infection", color = "Primary Treatment", shape = "Primary Treatment", linetype = "Primary Treatment")+
  scale_x_continuous(breaks = c(0, 1), labels = c("N", "Y"))+
  scale_y_continuous(breaks = c(0,1), labels = c("0", "1"))+
  scale_color_manual(values = c(pri_colors))+
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))+
  facet_grid(~secondary_dose)




#lm.sp2 <- glm(tes ~ elisa_od + primary_treatment*secondary_dose, data=m.ab, family=binomial(link= "logit"))


#install.packages("vcd")
library(vcd)
#mosaic plot
mosaic(~primary_treatment + secondary_dose + seropos_prim + diseased_sec, data=m.ab,
       highlighting = "primary_treatment", highlighting_fill = pri_colors)

mosaic(~primary_treatment + secondary_dose + seropos_prim + diseased_sec, data=m.ab,
       highlighting = "secondary_dose", highlighting_fill = sec_colors)

mosaic(~primary_treatment + secondary_dose + inf_prim + diseased_sec, data=m.ab,
       highlighting = "primary_treatment", highlighting_fill = pri_colors)

```
```{r infection during primary vs pathology secondary}
#does being infected during primary infeciton predict whether a bird will develop pathology upon secondary infection?
m.ab.inf <- m.ab %>%
  filter(primary_treatment != "Sham")
lm.inf <- glm(diseased_sec ~ inf_prim + primary_treatment*secondary_dose, data=m.ab.inf, family=binomial(link="logit"))
summary(lm.inf)
resid <- simulateResiduals(lm.inf)
plot(resid)
plot(allEffects(lm.inf))

dat.new=expand.grid(secondary_dose=unique(m.ab.inf$secondary_dose),
                    primary_treatment=unique(m.ab.inf$primary_treatment),
                    diseased_sec=unique(m.ab.inf$diseased_sec),
                    inf_prim=unique(m.ab.inf$inf_prim))#new grid to put predictions into
dat.new$yhat = predict(lm.inf, type="response", newdata=dat.new, re.form=NA) #predicted values
head(dat.new)




ggplot(m.ab.inf, aes(y=diseased_sec, x=inf_prim))+
  geom_jitter(size=1, height=0, width=0.5, aes(shape=primary_treatment))+
  geom_point(data = dat.new, aes(x =inf_prim, y = yhat, color=as.factor(secondary_dose), shape = primary_treatment)) +
  geom_path(data = dat.new, aes(x = inf_prim, y = yhat, color=as.factor(secondary_dose), linetype = primary_treatment))+
  labs(x="Infection Status During Primary Infection", y="Disease During Secondary Infection", color = "Secondary Dose", shape = "Primary Treatment", linetype = "Primary Treatment")+
  scale_x_continuous(breaks = c(0, 1), labels = c("N", "Y"))+
  scale_y_continuous(breaks = c(0,1), labels = c("0", "1"))+
  scale_linetype_manual(values = c("dashed", "solid"))+
  scale_shape_manual(values = c(17, 15))+
  facet_grid(~secondary_dose)

#with sham
lm.dis <- glm(diseased_sec ~ inf_prim + primary_treatment*secondary_dose, data=m.ab, family=binomial(link="logit"))
summary(lm.dis)
resid <- simulateResiduals(lm.dis)
plot(resid)

dat.new=expand.grid(secondary_dose=unique(m.ab$secondary_dose),
                    primary_treatment=unique(m.ab$primary_treatment),
                    diseased_sec=unique(m.ab$diseased_sec),
                    inf_prim=unique(m.ab$inf_prim))#new grid to put predictions into
dat.new$yhat = predict(lm.dis, type="response", newdata=dat.new, re.form=NA) #predicted values
head(dat.new)

m.ab.sec<- m.ab %>%
  filter(dpi %in% c(46, 49, 56, 63)) %>%
  group_by(dpi, primary_treatment, secondary_dose, inf_prim)%>%
  summarise(prop_seropos_sec = sum(seropos == 1, na.rm=TRUE) / n(), na.rm=TRUE,
            n_seropos_sec = sum(seropos == 1, na.rm=TRUE),
            n_seroneg_sec = sum(seropos == 0, na.rm=TRUE),
            prop_inf_sec = sum(inf == 1, na.rm=TRUE)/ n(), na.rm=TRUE,
            n_inf_sec = sum(inf == 1, na.rm=TRUE),
            n_not_inf_sec = sum(inf == 0, na.rm=TRUE),
            seropos_sec = sum(seropos_sec == 1, na.rm=TRUE),
            inf_sec = sum(inf_sec == 1, na.rm=TRUE),
            dis_sec = sum(diseased_sec == 1, na.rm =TRUE),
            total_birds =n(), na.rm=TRUE)

ggplot(m.ab, aes(y=diseased_sec, x=inf_prim))+
  geom_jitter(size=1, height=0, width=0.5, aes(shape=primary_treatment))+
  geom_point(data = dat.new, aes(x =inf_prim, y = yhat, color=as.factor(secondary_dose), shape = primary_treatment)) +
  geom_path(data = dat.new, aes(x = inf_prim, y = yhat, color=as.factor(secondary_dose), linetype = primary_treatment))+
  labs(x="Infection Status During Primary Infection", y="Disease During Secondary Infection", color = "Secondary Dose", shape = "Primary Treatment", linetype = "Primary Treatment")+
  scale_x_continuous(breaks = c(0, 1), labels = c("N", "Y"))+
  scale_y_continuous(breaks = c(0,1), labels = c("0", "1"))+
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))+
  facet_grid(~secondary_dose)

ggplot(m.ab, aes(y=diseased_sec, x=inf_prim))+
  geom_jitter(size=1, height=0, width=0.4, aes(shape=primary_treatment, color = primary_treatment))+
  geom_point(data = dat.new, aes(x =inf_prim, y = yhat, shape = primary_treatment, color = primary_treatment)) +
  geom_path(data = dat.new, aes(x = inf_prim, y = yhat, color = primary_treatment, linetype = primary_treatment))+
  geom_point(data=m.ab.sec, aes(x=inf_prim, y=dis_sec/total_birds, color=primary_treatment), shape = 45, size=10)+
  labs(x="Infection Status During Primary Infection", y="Disease During Secondary Infection", color = "Primary Treatment", shape = "Primary Treatment", linetype = "Primary Treatment")+
  scale_x_continuous(breaks = c(0, 1), labels = c("N", "Y"))+
  scale_y_continuous(breaks = c(0,1), labels = c("0", "1"))+
  scale_color_manual(values = c(pri_colors))+
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))+
  facet_grid(~secondary_dose)
```
```{r infection during primary vs infection secondary}
#infected secondary vs infected primary
lm.dis <- glm(inf_sec ~ inf_prim + primary_treatment*secondary_dose, data=m.ab, family=binomial(link="logit"))
summary(lm.dis)
resid <- simulateResiduals(lm.dis)
plot(resid)

dat.new=expand.grid(secondary_dose=unique(m.ab$secondary_dose),
                    primary_treatment=unique(m.ab$primary_treatment),
                    inf_sec=unique(m.ab$inf_sec),
                    inf_prim=unique(m.ab$inf_prim))#new grid to put predictions into
dat.new$yhat = predict(lm.dis, type="response", newdata=dat.new, re.form=NA) #predicted values
head(dat.new)

ggplot(m.ab, aes(y=diseased_sec, x=inf_prim))+
  geom_jitter(size=1, height=0, width=0.5, aes(shape=primary_treatment))+
  geom_point(data = dat.new, aes(x =inf_prim, y = yhat, color=as.factor(secondary_dose), shape = primary_treatment)) +
  geom_path(data = dat.new, aes(x = inf_prim, y = yhat, color=as.factor(secondary_dose), linetype = primary_treatment))+
  labs(x="Infection Status During Primary Infection", y="Disease During Secondary Infection", color = "Secondary Dose", shape = "Primary Treatment", linetype = "Primary Treatment")+
  scale_x_continuous(breaks = c(0, 1), labels = c("N", "Y"))+
  scale_y_continuous(breaks = c(0,1), labels = c("0", "1"))+
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))+
  facet_grid(~secondary_dose)

ggplot(m.ab, aes(y=diseased_sec, x=inf_prim))+
  geom_jitter(size=1, height=0, width=0.4, aes(shape=primary_treatment, color = primary_treatment))+
  geom_point(data = dat.new, aes(x =inf_prim, y = yhat, shape = primary_treatment, color = primary_treatment)) +
  geom_path(data = dat.new, aes(x = inf_prim, y = yhat, color = primary_treatment, linetype = primary_treatment))+
  labs(x="Infection Status During Primary Infection", y="Disease During Secondary Infection", color = "Primary Treatment", shape = "Primary Treatment", linetype = "Primary Treatment")+
  scale_x_continuous(breaks = c(0, 1), labels = c("N", "Y"))+
  scale_y_continuous(breaks = c(0,1), labels = c("0", "1"))+
  scale_color_manual(values = c(pri_colors))+
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))+
  facet_grid(~secondary_dose)

ggplot(m.ab, aes(x=inf_prim, y=inf_sec, color=primary_treatment))+
  #geom_jitter(width=0.4, height =0.05, alpha = 0.5)+
  geom_count(position=position_jitter(height = 0, width = 0.3))+
  geom_point(data = dat.new, aes(x =inf_prim, y = yhat, shape = primary_treatment, color = primary_treatment)) +
  geom_path(data = dat.new, aes(x = inf_prim, y = yhat, color = primary_treatment, linetype = primary_treatment))+
  scale_x_continuous(breaks = c(0, 1), labels = c("N", "Y"))+
  scale_y_continuous(breaks = c(0,1), labels = c("0", "1"))+
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))+
  scale_color_manual(values = c(pri_colors))+
  labs(x="Infected During Primary Infection", y="Infection Status During Secondary Infection", color = "Primary Treatment", shape = "Primary Treatment", linetype = "Primary Treatment")+
  facet_wrap(~secondary_dose, nrow=1)

ggplot(m.ab, aes(x=inf_prim, y=inf_sec, color=primary_treatment))+
  geom_jitter(width=0.4, height =0, alpha = 0.5)+
  geom_point(data = dat.new, aes(x =inf_prim, y = yhat, shape = primary_treatment, color = primary_treatment)) +
  geom_path(data = dat.new, aes(x = inf_prim, y = yhat, color = primary_treatment, linetype = primary_treatment))+
  scale_x_continuous(breaks = c(0, 1), labels = c("N", "Y"))+
  scale_y_continuous(breaks = c(0,1), labels = c("0", "1"))+
  scale_color_manual(values = c(pri_colors))+
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))+
  labs(x="Infected During Primary Infection", y="Infection Status During Secondary Infection", color = "Primary Treatment", shape = "Primary Treatment", linetype = "Primary Treatment")+
  facet_wrap(~secondary_dose, nrow=1)
```
```{r Are antibody levels on day 14 predictive of susceptibility upon secondary challenge}
wibird
#remove birds that weren't recovered by secondary inoculation
wibird <- wibird %>%
  filter(!(band_number %in% c(2274, 2514, 2469, 2520, 2494)))
wibird$log10.sec_dose <- round(log10(wibird$secondary_dose+0.001), digits = 5)
glm.ab14 <- glm(inf_sec ~ elisa_od_14 + log10.sec_dose, data=wibird, family=binomial())
hist(resid(glm.ab14))

summary(glm.ab14)
car::Anova(glm.ab14, type = 3)
simulateResiduals(glm.ab14, plot=T)

library(epitools)
dat.new=expand.grid(log10.sec_dose=unique(wibird$log10.sec_dose),
                    inf_sec=unique(wibird$inf_sec),
                    elisa_od_14=unique(wibird$elisa_od_14))#new grid to put predictions into
dat.new$yhat=predict(glm.ab14, type="response", newdata = dat.new)
#prediction intervals
preds = predict(glm.ab14, type = "link", newdata = dat.new, se.fit =T)
#bind se's and fitted points
dat.new = cbind(dat.new, preds)
#inverse link function
ilink <- family(glm.ab14)$linkinv
#back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

ggplot(wibird, aes(x=(elisa_od_14), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_point(alpha=0.5)+
  geom_line(data=dat.new, aes(x=(elisa_od_14), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_14, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.1) + 
  scale_color_manual(values = c(sec_colors))+
  scale_fill_manual(values = c(sec_colors))+
  labs(x="Antibody Levels Day 14", y= "Susceptibility (Secondary Infection 0|1)", color="log10(Secondary Dose)", fill ="log10(Secondary Dose)")+
  facet_wrap(~log10.sec_dose, nrow=1)+
  scale_x_continuous(name = "Antibody Levels Day 14",
                     labels = scales::number_format(accuracy = 0.05))+
  scale_x_reverse()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```
```{r eval=FALSE, include=FALSE}
#does being infected during primary infection predict whether a bird will develop pathology upon secondary infection?
m.ab.inf <- m.ab1# %>%
#  filter(primary_treatment != "Sham")
lm.inf <- glm(inf_sec ~ inf_prim + seropos_prim + primary_treatment*secondary_dose, data=m.ab.inf, family=binomial(link="logit"))
summary(lm.inf)
plot(allEffects(lm.inf))

dat.new=expand.grid(secondary_dose=unique(m.ab.inf$secondary_dose),
                    primary_treatment=unique(m.ab.inf$primary_treatment),
                    diseased_sec=unique(m.ab.inf$diseased_sec),
                    inf_prim=unique(m.ab.inf$inf_prim))#new grid to put predictions into
dat.new$yhat = predict(lm.inf, type="response", newdata=dat.new, re.form=NA) #predicted values
head(dat.new)

ggplot(m.ab.inf, aes(y=diseased_sec, x=inf_prim))+
  geom_jitter(size=1, height=0, width=0.5, aes(shape=primary_treatment))+
  geom_point(data = dat.new, aes(x =inf_prim, y = yhat, color=as.factor(secondary_dose), shape = fct_rev(primary_treatment))) +
  geom_path(data = dat.new, aes(x = inf_prim, y = yhat, color=as.factor(secondary_dose), linetype = fct_rev(primary_treatment)))+
  labs(x="Infection Status During Primary Infection", y="Relative Likelihood of Disease During Secondary Infection", color = "Secondary Dose", shape = "Primary Treatment", linetype = "Primary Treatment")+
  scale_x_continuous(breaks = c(0, 1), labels = c("N", "Y"))+
  scale_y_continuous(breaks = c(0,1), labels = c("N", "Y"))+
  facet_grid(~secondary_dose)
```
```{r Variance across DPI}
####Residual Models####
p.ab <- p.ab %>%
  filter(dpi <= 41 )%>%
  drop_na(band_number, dpi, elisa_od)



#Levene's Test for Equality of Variances
library(car)
leveneTest(elisa_od ~ primary_treatment, data=p.ab)
##Variance is not equal##

#glm <- glmmTMB(elisa_od ~ 1 + (band_number|dpi), data=p.ab, family=Gamma())
glm.a <- glmmTMB(elisa_od ~ 1 + (1|dpi), data=p.ab, family=Gamma())

summary(glm.a)
resid(glm.a)
hist(resid(glm.a))
p.ab$resid <- resid(glm.a)
hist(p.ab$resid)

ggplot(p.ab, aes(x=dpi, y=elisa_od, color=fct_rev(primary_treatment)))+
  geom_point()

ggplot(p.ab, aes(x=dpi, y=resid, color=fct_rev(primary_treatment)))+
  geom_point()+
  labs(x="Days Post Infection", y="Residuals", color="Primary Treatment")+
  facet_wrap(~primary_treatment)

lm1 <- lm(abs(resid) ~ primary_treatment, data=p.ab)
summary(lm1)
```
```{r eval=FALSE, include=FALSE}
ab.41 <- p.ab %>%
  filter(dpi == 41)

glm.all <- glmmTMB(elisa_od~1 + (1|band_number), data=ab.41, family=Gamma())
glm.sham <- glmmTMB(elisa_od~1, data=subset(ab.41, primary_treatment == "Sham"), family=Gamma())
glm.low <- glmmTMB(elisa_od~1, data=subset(ab.41, primary_treatment == "Low"), family=Gamma())
glm.high <- glmmTMB(elisa_od~1, data=subset(ab.41, primary_treatment == "High"), family=Gamma())

hist(abs(resid(glm.sham)))
hist(abs(resid(glm.low)))
hist(abs(resid(glm.high)))

summary(glm.sham)
ab.1441$resid <- resid(glm.all)
lm.a <- lm(abs(resid)~primary_treatment, data= ab.1441)
summary(lm.a)
ggplot(ab.1441, aes(x=primary_treatment, y=resid))+
  geom_jitter(width=0.25)


ggplot(ab.1441, aes(x=primary_treatment, y=abs(resid)))+
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
```
