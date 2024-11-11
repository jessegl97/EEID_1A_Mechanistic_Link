####Eye Score Across Priming Infection####
p.abe <- p.ab
p.abe$band_number <- as.factor(p.abe$band_number)
p.abe$dpi.f <- as.factor(p.abe$dpi)

ggplot(p.abe, aes(x=dpi, y=tes, color=fct_rev(primary_treatment)))+
  geom_jitter(height=0.1, width=0.5, alpha=0.5)+
  stat_summary(aes(x=dpi, y=tes, group=primary_treatment, fill=fct_rev(primary_treatment)),
               fun.y=mean, geom="point", size=2, shape=25)+
  stat_summary(aes(x=dpi, y=tes, group=primary_treatment),
             fun.y=mean, geom="line", size=1)

pv.ep <- p.abe %>%
  dplyr::select(dpi, tes, primary_treatment)

pv.ep$tes <- pv.ep$tes+0.001

pv.ep <- pv.ep %>%
  group_by(dpi, primary_treatment)%>%
  drop_na(tes)%>%
  dplyr::summarise(mean_tes = mean(tes),
                   bird_cv = sd(tes)/mean(tes),
                   bird_sd = sd(tes),
                   bird_pv = PV(tes),
                   elisa_od = tes)
pv.ep

p.abepv<- left_join(p.abe, pv.ep)

ggplot(p.abe, aes(x=dpi, y=tes, color=fct_rev(primary_treatment)))+
  geom_jitter(height=0.1, width=0.5, alpha=0.5)+
  # stat_summary(aes(x=dpi, y=tes, group=primary_treatment, fill=fct_rev(primary_treatment)),
  #             fun.y=mean, geom="point", size=2, shape=25)+
  # stat_summary(aes(x=dpi, y=tes, group=primary_treatment),
  #             fun.y=mean, geom="line", size=1)+
  geom_point(data=pv.ep, aes(x=dpi, y=bird_pv, color=fct_rev(primary_treatment)), shape=15, size=3)+
  geom_line(data=pv.ep, aes(x=dpi, y=bird_pv, color=fct_rev(primary_treatment)), linetype="dashed")


p.abe$tes <- p.abe$tes+0.001
e.max <- p.abe %>%
  dplyr::select(band_number, tes, primary_treatment, dpi)%>%
  group_by(band_number)%>%
  mutate(max_tes = max(tes))

e.max <- e.max %>%
  filter(dpi==14)

lm.tes.max <- glmmTMB(max_tes~ primary_treatment, data=e.max, family=Gamma())
summary(lm.tes.max)

mod.pred <- lm.tes.max
#Model predictions
dat.new=expand.grid(primary_treatment=unique(e.max$primary_treatment),
                    max_tes = unique(e.max$max_tes))
#sex = unique(p.abt$sex))#new grid to put predictions into
dat.new$yhat = predict(mod.pred, type="response", newdata=dat.new, re.form=NA)
preds = predict(mod.pred, type="link", newdata=dat.new, se.fit=TRUE, re.form=NA)
dat.new = cbind(dat.new, preds)

ilink <- family(mod.pred)$linkinv
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

e.cv.m <- e.max %>%
  group_by(primary_treatment)%>%
  dplyr::reframe(
    group_cv = sd(tes)/mean(tes),
    group_sd = sd(tes),
    group_pv = PV(tes),
    max_tes = max_tes)
e.max.v <- left_join(e.max, e.cv.m)



ggplot(e.cv.m)+
  geom_crossbar(data=dat.new, aes(ymin=Lower, ymax=Upper, x=primary_treatment, y=yhat, color=fct_rev(primary_treatment),
                                  fill=fct_rev(primary_treatment)), alpha=0.01, width=0.5, stroke=0.01)+
  geom_jitter(aes(x=primary_treatment, y=max_tes, color=primary_treatment),
             alpha=0.5, size=4, height=0, width=0.1)+
  geom_point(data=e.cv.m, aes(y=group_pv*10, x=primary_treatment), shape=25, size=4)+
  geom_point(data=e.cv.m, aes(y=group_pv*10, x=primary_treatment), shape=25, size=4, alpha=0.01)+
  scale_color_manual(values=pri_colors)+
  labs(x="Primary Treatment", y="Max TES", color="Primary Treatment", fill="Primary Treatment")+
  scale_y_continuous(
    sec.axi=sec_axis(~.*1, name="PV*10")
  )+
  theme(
    axis.title.y = element_text(color = "black", size=11),
    axis.text.y = element_text(color="black"),
    axis.text.y.right = element_text(color="black", face="bold"),
    axis.title.y.right = element_text(color = "black", size=11, face = "bold")
  )

p.abe$tes1 <- p.abe$tes+0.001

ggplot(p.abe, aes(x=tes, color=fct_rev(primary_treatment), fill=fct_rev(primary_treatment)))+
         geom_histogram()+
  facet_wrap(~primary_treatment)

lm1e <- glmmTMB(tes ~ primary_treatment*as.factor(dpi) + (1|band_number),
               data=p.abe %>% filter(dpi %in% c(-8, 14, 41)),
               ziformula = ~ primary_treatment,
               family = poisson)

summary(lm1e)
simulateResiduals(lm1e, plot=T)

lm1test <- glmmTMB(tes ~ primary_treatment*as.factor(dpi)+ (1|band_number), data=p.abe)
resids <- simulateResiduals(lm1test)
testZeroInflation(resids)

mod.pred <- lm1e
#Model predictions
dat.new=expand.grid(primary_treatment=unique(p.abe$primary_treatment),
                    tes = unique(p.abe$tes),
                    dpi = unique(p.abe$dpi))
#generate yhat values
dat.new$yhat = predict(mod.pred, type="response", newdata=dat.new, re.form=NA)
#generate SEM
preds = predict(mod.pred, type="link", newdata=dat.new, se.fit=TRUE, re.form=NA)
dat.new = cbind(dat.new, preds)

ilink <- family(mod.pred)$linkinv
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

pv.ep1 <- p.abe %>%
  dplyr::select(dpi, tes, primary_treatment)

pv.ep1$tes <- pv.ep1$tes+0.001

pv.ep1 <- pv.ep1 %>%
  group_by(dpi, primary_treatment)%>%
  drop_na(tes)%>%
  dplyr::summarise(mean_tes = mean(tes),
                   bird_cv = sd(tes)/mean(tes),
                   bird_sd = sd(tes),
                   bird_pv = PV(tes),
                   elisa_od = tes)
pv.ep1


ggplot(p.abe, aes(x=dpi))+
  geom_jitter(aes(x=as.numeric(dpi), y=tes, color=fct_rev(primary_treatment)),
              alpha=0.5, size=2, height=0, width=01)+
  geom_point(data=dat.new, aes(y=yhat, x=as.numeric(dpi), color=primary_treatment), shape=16, size=3)+
  geom_point(data=dat.new, aes(y=yhat, x=as.numeric(dpi), color=primary_treatment), shape=1, size=3, color="black")+
  geom_path(data=dat.new, aes(y=yhat, x=as.numeric(dpi), color=primary_treatment))+
  geom_errorbar(data=dat.new, aes(y=yhat, x=as.numeric(dpi), ymax=Upper, ymin=Lower), width=0.25)+
  geom_point(data=pv.ep1, aes(y=bird_pv*10, x=as.numeric(dpi), fill=fct_rev(primary_treatment)), shape=25, size=4)+
  geom_point(data=pv.ep1, aes(y=bird_pv*10, x=as.numeric(dpi)), color="black", shape=25, size=4)+
  geom_path(data=pv.ep1, aes(y=bird_pv*10, x=as.numeric(dpi), color=primary_treatment), linetype="dotdash")+
  scale_color_manual(values=pri_colors)+
  scale_fill_manual(values=pri_colors)+
  labs(x="Days Post Infection", y="Total Eye Score", color="Primary Treatment", fill="Primary Treatment")+
  scale_y_continuous(
    sec.axi=sec_axis(~.*1, name="PV*10")
  )+
  theme(
    axis.title.y = element_text(color = "black", size=11),
    axis.text.y = element_text(color="black"),
    axis.text.y.right = element_text(color="black", face="bold"),
    axis.title.y.right = element_text(color = "black", size=11, face = "bold")
  )
