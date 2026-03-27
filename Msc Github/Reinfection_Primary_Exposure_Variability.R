#Reinfection Primary Exposure Variability Code
##directly cut and pasted out of Reinfection_Final_Code.R line ~956

####Calculate Variability in Total Eye Score Primary####
#Using all individuals regardless of infection status
e.cv <- p.aba %>% 
  group_by(dpi, primary_treatment) %>%
  reframe(
    tes.new = tes.new,
    mean_tes = mean(tes.new),
    band_number = band_number,
    bird_cv = calculate_cv(tes.new),
    bird_sd = sd(tes.new),
    bird_pv = calculate_pv(tes.new),
    bird_v2 = calculate_v2(tes.new),
    dpi.f = dpi.f,
    tes.new = tes.new,
    
    # Bootstrap for confidence intervals
    cv_bootstrap = list(replicate(n_boot, calculate_cv(sample(tes.new, replace = TRUE)))),
    v2_bootstrap = list(replicate(n_boot, calculate_v2(sample(tes.new, replace = TRUE)))),
    pv_bootstrap = list(replicate(n_boot, calculate_pv(sample(tes.new, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals
  mutate(
    cv_lower_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.025)),
    cv_upper_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.975)),
    v2_lower_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.025)),
    v2_upper_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.975)),
    pv_lower_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.025)),
    pv_upper_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.975))
  ) %>%
  # Remove bootstrap columns
  select(-cv_bootstrap, -v2_bootstrap, -pv_bootstrap) %>%
  ungroup()

#add to df and format names
p.aba.e <- left_join(p.aba, e.cv, by=c("dpi.f", "primary_treatment", "band_number"))
p.aba.e$tes.new <- p.aba.e$tes.new.x
p.aba.e <- p.aba.e %>%
  select(-tes.new.y, -tes.new.x)

#just days samples were taken
p.aba.e <- p.aba.e %>%
  filter(dpi.f %in% c(-8, 7, 14, 21, 28, 35, 41))

#summary table
summary_tibble.e <- e.cv %>%
  group_by(dpi.f, primary_treatment) %>%
  summarize(
    CV = mean(bird_cv, na.rm = TRUE),    # Calculate the mean of bird_cv for each group
    SD = mean(bird_sd, na.rm = TRUE),    # Calculate the mean of bird_sd for each group
    V2 = mean(bird_v2, na.rm = TRUE),    # Calculate the mean of bird_v2 for each group
    PV = mean(bird_pv, na.rm = TRUE),    # Calculate the mean of bird_pv for each group
    tes.new = mean(tes.new, na.rm = TRUE),# Calculate the average elisa_od for each group
    lower_ci_pv = mean(pv_lower_ci),
    upper_ci_pv = mean(pv_upper_ci),
    lower_ci_cv = mean(cv_lower_ci),
    upper_ci_cv = mean(cv_upper_ci),
    lower_ci_v2 = mean(v2_lower_ci),
    upper_ci_v2 = mean(v2_upper_ci),
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  ) %>%
  as_tibble()
print(summary_tibble.e)
#write.csv(summary_tibble.e, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/ELISA Variability Primary/tes_var_not_removed.csv", row.names=FALSE)

#visualize variability in eye score
ggplot(summary_tibble.e, aes(x=dpi.f, color=primary_treatment))+
  geom_point(aes(x=dpi.f, y=PV, color=primary_treatment, shape="PV"), size=2) +
  geom_point(aes(x=dpi.f, y=V2, color=primary_treatment, shape="V2"), size=2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=dpi.f, ymin = lower_ci_pv, ymax = upper_ci_pv, y=PV, color=primary_treatment), width=0.2, lty="dashed")+
  geom_errorbar(aes(x=dpi.f, ymin = lower_ci_v2, ymax = upper_ci_v2, y=V2, color=primary_treatment), position=position_dodge(width=0.5) , width=0.2)+
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 17, "V2" = 18),  # Assign shapes
                     labels = c("PV" = "PV", "V2" = "CV^2")) +
  scale_color_manual(values=(pri_colors))+
  labs(x="Days Post Inoculation", y="Variability", shape="Variability Metric", color="Primary Treatment")

####Calculate Variability in Raw Pathogen Load Primary####
####Raw Pathogen Load###
#Using all individuals regardless of infection status use p.aba
#Use p.abi for only infected
q.cv <- p.aba%>% 
  group_by(dpi, primary_treatment) %>%
  drop_na(quantity1) %>%
  reframe(
    bird_cv = calculate_cv(quantity1),
    bird_pv = calculate_pv(quantity1),
    bird_v2 = calculate_v2(quantity1),
    dpi.f = dpi.f,
    quantity1 = quantity1,
    inf_pri = inf_pri,
    band_number = band_number,
    
    # Bootstrap for confidence intervals
    cv_bootstrap = list(replicate(n_boot, calculate_cv(sample(quantity1, replace = TRUE)))),
    v2_bootstrap = list(replicate(n_boot, calculate_v2(sample(quantity1, replace = TRUE)))),
    pv_bootstrap = list(replicate(n_boot, calculate_pv(sample(quantity1, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals
  mutate(
    cv_lower_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.025)),
    cv_upper_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.975)),
    v2_lower_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.025)),
    v2_upper_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.975)),
    pv_lower_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.025)),
    pv_upper_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.975))
  ) %>%
  # Remove bootstrap columns
  select(-cv_bootstrap, -v2_bootstrap, -pv_bootstrap) %>%
  ungroup()

#add to df and format names
p.aba.q <- left_join(p.aba, q.cv, by=c("dpi.f", "primary_treatment", "band_number"))
p.aba.q$quantity1 <- p.aba.q$quantity1.x
p.aba.q <- p.aba.q %>%
  select(-quantity1.y, -quantity1.x)

#just days samples were taken
p.aba.q <- p.aba.q %>%
  filter(dpi.f %in% c(-8, 7, 41))

#summary table
summary_tibble.q <- q.cv %>%
  group_by(dpi.f, primary_treatment) %>%
  summarize(
    CV = mean(bird_cv, na.rm = TRUE),    # Calculate the mean of bird_cv for each group
    V2 = mean(bird_v2, nar.rm = TRUE),    # Calculate the mean of bird_v2 for each group
    PV = mean(bird_pv, na.rm = TRUE),    # Calculate the mean of bird_pv for each group
    quantity1 = mean(quantity1, na.rm = TRUE),# Calculate the average elisa_od for each group
    lower_ci_pv = mean(pv_lower_ci),
    upper_ci_pv = mean(pv_upper_ci),
    lower_ci_cv = mean(cv_lower_ci),
    upper_ci_cv = mean(cv_upper_ci),
    lower_ci_v2 = mean(v2_lower_ci),
    upper_ci_v2 = mean(v2_upper_ci),
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  ) %>%
  as_tibble()
print(summary_tibble.q)
#write.csv(summary_tibble.q, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/ELISA Variability Primary/quant_var_not_removed.csv", row.names=FALSE)

#visualize variability in pathogen load
ggplot(summary_tibble.q, aes(x=dpi.f, color=primary_treatment))+
  geom_point(aes(x=dpi.f, y=PV, color=primary_treatment, shape="PV"), size=2) +
  geom_point(aes(x=dpi.f, y=V2, color=primary_treatment, shape="V2"), size=2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=dpi.f, ymin = lower_ci_pv, ymax = upper_ci_pv, y=PV, color=primary_treatment), width=0.2, lty="dashed")+
  geom_errorbar(aes(x=dpi.f, ymin = lower_ci_v2, ymax = upper_ci_v2, y=V2, color=primary_treatment), position=position_dodge(width=0.5) , width=0.2)+
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 17, "V2" = 18),  # Assign shapes
                     labels = c("PV" = "PV", "V2" = "CV^2")) +
  scale_color_manual(values=(pri_colors))+
  labs(x="Days Post Inoculation", y="Variability", shape="Variability Metric", color="Primary Treatment")

ggplot(p.abi, aes(x=dpi.f, y=quantity1, color=primary_treatment))+
  geom_jitter(width=0.25, height=0)+
  scale_color_manual(values=(pri_colors))+
  scale_y_log10()




####Log10 Pathogen Load####
#log10 variability
p.aba$quantity1 <- p.aba$quantity1+0.0001
p.abi$quantity1 <- p.abi$quantity1+0.0001
q.cv.10 <- p.aba %>% 
  group_by(dpi, primary_treatment) %>%
  drop_na(quantity1) %>%
  reframe(
    bird_cv = calculate_cv(log10(quantity1)),
    bird_pv = calculate_pv(log10(quantity1)),
    bird_v2 = calculate_v2(log10(quantity1)),
    dpi.f = dpi.f,
    quantity1 = quantity1,
    band_number = band_number,
    
    # Bootstrap for confidence intervals
    cv_bootstrap = list(replicate(n_boot, calculate_cv(sample(log10(quantity1), replace = TRUE)))),
    v2_bootstrap = list(replicate(n_boot, calculate_v2(sample(log10(quantity1), replace = TRUE)))),
    pv_bootstrap = list(replicate(n_boot, calculate_pv(sample(log10(quantity1), replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals
  mutate(
    cv_lower_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  # Remove bootstrap columns
  select(-cv_bootstrap, -v2_bootstrap, -pv_bootstrap) %>%
  ungroup()

#add to df and format names
p.aba.q10 <- left_join(p.aba, q.cv.10, by=c("dpi.f", "primary_treatment", "band_number"))
p.aba.q10$quantity1 <- p.aba.q10$quantity1.x
p.aba.q10 <- p.aba.q10 %>%
  select(-quantity1.y, -quantity1.x)

#just days samples were taken
p.aba.q10 <- p.aba.q10 %>%
  filter(dpi.f %in% c(-8, 7, 41))

#summary table
summary_tibble.q10 <- q.cv.10 %>%
  group_by(dpi.f, primary_treatment) %>%
  summarize(
    CV = mean(bird_cv, na.rm = TRUE),    # Calculate the mean of bird_cv for each group
    V2 = mean(bird_v2, nar.rm = TRUE),    # Calculate the mean of bird_v2 for each group
    PV = mean(bird_pv, na.rm = TRUE),    # Calculate the mean of bird_pv for each group
    quantity1 = mean(quantity1, na.rm = TRUE),# Calculate the average elisa_od for each group
    lower_ci_pv = mean(pv_lower_ci),
    upper_ci_pv = mean(pv_upper_ci),
    lower_ci_cv = mean(cv_lower_ci),
    upper_ci_cv = mean(cv_upper_ci),
    lower_ci_v2 = mean(v2_lower_ci),
    upper_ci_v2 = mean(v2_upper_ci),
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  ) %>%
  as_tibble()
print(summary_tibble.q10)
#write.csv(summary_tibble.q10, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/ELISA Variability Primary/log10quant_var_not_removed.csv", row.names=FALSE)

#visualize variability in pathogen load
ggplot(summary_tibble.q10, aes(x=dpi.f, color=primary_treatment))+
  geom_point(aes(x=dpi.f, y=PV, color=primary_treatment, shape="PV"), size=2) +
  #geom_point(aes(x=dpi.f, y=V2, color=primary_treatment, shape="V2"), size=2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(x=dpi.f, ymin = lower_ci_pv, ymax = upper_ci_pv, y=PV, color=primary_treatment), width=0.2, lty="dashed")+
  #geom_errorbar(aes(x=dpi.f, ymin = lower_ci_v2, ymax = upper_ci_v2, y=V2, color=primary_treatment), position=position_dodge(width=0.5) , width=0.2)+
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 17, "V2" = 18),  # Assign shapes
                     labels = c("PV" = "PV", "V2" = "CV^2")) +
  scale_color_manual(values=(pri_colors))+
  labs(x="Days Post Inoculation", y="Variability", shape="Variability Metric", color="Primary Treatment")

ggplot(p.aba, aes(x=dpi.f, y=quantity1, color=primary_treatment))+
  geom_hline(yintercept = 51)+
  geom_jitter(width=0.25, height=0)+
  scale_color_manual(values=(pri_colors))+
  scale_y_log10()

#Compare
ggplot(summary_tibble.q10)+
  geom_point(aes(x=dpi.f, y=PV, color=primary_treatment, shape="PV Log"), position = position_dodge(width=0.3))+
  geom_errorbar(aes(x=dpi.f, ymin = lower_ci_pv, ymax = upper_ci_pv, y=PV, color=primary_treatment),
                position = position_dodge(width=0.3), width=0.1, lty="solid")+
  geom_point(data=summary_tibble.q, aes(x=dpi.f, y=PV, color=primary_treatment, shape="PV Raw"), size=3, position=position_dodge(width=0.5))+
  geom_errorbar(data=summary_tibble.q, aes(x=dpi.f, ymin = lower_ci_pv, ymax = upper_ci_pv, y=PV, color=primary_treatment),
                width=0.2, lty="dashed", position=position_dodge(width=0.5))+
  scale_color_manual(values=(pri_colors))+
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV Log" = 17, "PV Raw" = 23),  # Assign shapes
                     labels = c("PV Log" = "PV Log", "PV Raw" = "PV Raw")) +
  labs(title="Log10(Quantity) vs Raw Quantity Variability", y="Variability", x="DPI")

####Variability in Max Eyescore Primary####

summary_tibble.e
summary_tibble.q


#Primary challenge variability; eye score
#Using all individuals regardless of infection status
e.cv.s <- p.aba %>% 
  group_by(dpi, primary_treatment) %>%
  reframe(
    tes.new = tes.new,
    mean_tes = mean(tes.new),
    band_number = band_number,
    bird_cv = calculate_cv(tes.new),
    bird_sd = sd(tes.new),
    bird_pv = calculate_pv(tes.new),
    bird_v2 = calculate_v2(tes.new),
    dpi.f = dpi.f,
    tes.new = tes.new,
    
    # Bootstrap for confidence intervals
    cv_bootstrap = list(replicate(n_boot, calculate_cv(sample(tes.new, replace = TRUE)))),
    v2_bootstrap = list(replicate(n_boot, calculate_v2(sample(tes.new, replace = TRUE)))),
    pv_bootstrap = list(replicate(n_boot, calculate_pv(sample(tes.new, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals
  mutate(
    cv_lower_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.025)),
    cv_upper_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.975)),
    v2_lower_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.025)),
    v2_upper_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.975)),
    pv_lower_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.025)),
    pv_upper_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.975))
  ) %>%
  # Remove bootstrap columns
  select(-cv_bootstrap, -v2_bootstrap, -pv_bootstrap) %>%
  ungroup()

#add to df and format names
p.aba.e <- left_join(p.aba, e.cv, by=c("dpi.f", "primary_treatment", "band_number"))
p.aba.e$tes.new <- p.aba.e$tes.new.x
p.aba.e <- p.aba.e %>%
  select(-tes.new.y, -tes.new.x)

#just days samples were taken
p.aba.e <- p.aba.e %>%
  filter(dpi.f %in% c(-8, 7, 14, 21, 28, 35, 41))

#summary table
summary_tibble.e.s <- e.cv %>%
  group_by(dpi.f, primary_treatment) %>%
  summarize(
    CV = mean(bird_cv, na.rm = TRUE),    # Calculate the mean of bird_cv for each group
    SD = mean(bird_sd, na.rm = TRUE),    # Calculate the mean of bird_sd for each group
    V2 = mean(bird_v2, na.rm = TRUE),    # Calculate the mean of bird_v2 for each group
    PV = mean(bird_pv, na.rm = TRUE),    # Calculate the mean of bird_pv for each group
    tes.new = mean(tes.new, na.rm = TRUE),# Calculate the average elisa_od for each group
    lower_ci_pv = mean(pv_lower_ci),
    upper_ci_pv = mean(pv_upper_ci),
    lower_ci_cv = mean(cv_lower_ci),
    upper_ci_cv = mean(cv_upper_ci),
    lower_ci_v2 = mean(v2_lower_ci),
    upper_ci_v2 = mean(v2_upper_ci),
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  ) %>%
  as_tibble()
print(summary_tibble.e.s)

#####Max pathogen load and tes of infected birds only Primary####
#Calculate max eye score and pathogen load
#Calculate max pathogen load primary
max.quant.pri <- p.aba %>%
  group_by(band_number, primary_treatment) %>%
  reframe(max_quantity = max(quantity1, na.rm = TRUE),
          inf_pri = max(inf_pri))

#How many infected primary
max.quant.pri %>%
  dplyr::select(primary_treatment, inf_pri)%>%
  group_by(primary_treatment, inf_pri) %>%
  tbl_summary(
    by = inf_pri
  )%>%
  modify_header(
    label ~ "**Inf Pri By Primary Dose**",
  )

#Calculate max eye score primary
max.tes.pri <- p.aba %>%
  group_by(band_number, primary_treatment)%>%
  reframe(max_tes = max(tes.new, na.rm=TRUE),
          inf_pri = max(inf_pri))

max.pri<-left_join(max.quant.pri, max.tes.pri, by="band_number")
max.pri$inf_pri <- max.pri$inf_pri.x
max.pri$primary_treatment <- max.pri$primary_treatment.x
max.pri <- max.pri %>%
  dplyr::select(-inf_pri.x, -inf_pri.y, -primary_treatment.x, -primary_treatment.y)


#variability in primary infected only
max.p.v <- max.pri %>% 
  filter(inf_pri == 1)%>%
  group_by(primary_treatment) %>%
  reframe(
    # Metrics for max_tes
    max_tes = max_tes,
    mean_tes = mean(max_tes-0.01),
    se_tes = SE(max_tes-0.01),
    inf_pri = inf_pri,
    bird_cv_tes = calculate_cv(max_tes),
    bird_sd_tes = sd(max_tes),
    bird_pv_tes = calculate_pv(max_tes),
    bird_v2_tes = calculate_v2(max_tes),
    
    # Metrics for max_quantity
    max_quantity = max_quantity,
    mean_quantity = mean(max_quantity),
    se_quantity = SE(max_quantity -1),
    bird_cv_quantity = calculate_cv(max_quantity),
    bird_sd_quantity = sd(max_quantity),
    bird_pv_quantity = calculate_pv(max_quantity),
    bird_v2_quantity = calculate_v2(max_quantity),
    
    #Other info
    band_number = band_number,
    primary_treatment = primary_treatment,
    
    # Bootstrap for confidence intervals (max_tes)
    cv_bootstrap_tes = list(replicate(n_boot, calculate_cv(sample(max_tes, replace = TRUE)))),
    v2_bootstrap_tes = list(replicate(n_boot, calculate_v2(sample(max_tes, replace = TRUE)))),
    pv_bootstrap_tes = list(replicate(n_boot, calculate_pv(sample(max_tes, replace = TRUE)))),
    
    # Bootstrap for confidence intervals (max_quantity)
    cv_bootstrap_quantity = list(replicate(n_boot, calculate_cv(sample(max_quantity, replace = TRUE)))),
    v2_bootstrap_quantity = list(replicate(n_boot, calculate_v2(sample(max_quantity, replace = TRUE)))),
    pv_bootstrap_quantity = list(replicate(n_boot, calculate_pv(sample(max_quantity, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals for both max_tes and max_quantity
  mutate(
    # max_tes CIs
    cv_lower_ci_tes = map_dbl(cv_bootstrap_tes, ~ quantile(.x, 0.025)),
    cv_upper_ci_tes = map_dbl(cv_bootstrap_tes, ~ quantile(.x, 0.975)),
    v2_lower_ci_tes = map_dbl(v2_bootstrap_tes, ~ quantile(.x, 0.025)),
    v2_upper_ci_tes = map_dbl(v2_bootstrap_tes, ~ quantile(.x, 0.975)),
    pv_lower_ci_tes = map_dbl(pv_bootstrap_tes, ~ quantile(.x, 0.025)),
    pv_upper_ci_tes = map_dbl(pv_bootstrap_tes, ~ quantile(.x, 0.975)),
    
    # max_quantity CIs
    cv_lower_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.025)),
    cv_upper_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.975)),
    v2_lower_ci_quantity = map_dbl(v2_bootstrap_quantity, ~ quantile(.x, 0.025)),
    v2_upper_ci_quantity = map_dbl(v2_bootstrap_quantity, ~ quantile(.x, 0.975)),
    pv_lower_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.025)),
    pv_upper_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.975))
  ) %>%
  # Remove bootstrap columns
  # select(-cv_bootstrap_tes, -v2_bootstrap_tes, -pv_bootstrap_tes,
  #        -cv_bootstrap_quantity, -v2_bootstrap_quantity, -pv_bootstrap_quantity) %>%
  ungroup()

# Add to df and format names
max.pri.v <- left_join(max.pri, max.p.v, by = c("primary_treatment", "band_number"))

# Summary table for both max_tes and max_quantity
summary_tibble.p.v <- max.p.v %>%
  group_by( primary_treatment) %>%
  summarize(
    # Summary statistics for max_tes
    CV_tes = mean(bird_cv_tes, na.rm = TRUE),
    SD_tes = mean(bird_sd_tes, na.rm = TRUE),
    V2_tes = mean(bird_v2_tes, na.rm = TRUE),
    PV_tes = mean(bird_pv_tes, na.rm = TRUE),
    SE_tes = mean(se_tes, na.rm = TRUE),
    
    lower_ci_pv_tes = mean(pv_lower_ci_tes),
    upper_ci_pv_tes = mean(pv_upper_ci_tes),
    lower_ci_cv_tes = mean(cv_lower_ci_tes),
    upper_ci_cv_tes = mean(cv_upper_ci_tes),
    lower_ci_v2_tes = mean(v2_lower_ci_tes),
    upper_ci_v2_tes = mean(v2_upper_ci_tes),
    
    # Summary statistics for max_quantity
    CV_quantity = mean(bird_cv_quantity, na.rm = TRUE),
    SD_quantity = mean(bird_sd_quantity, na.rm = TRUE),
    V2_quantity = mean(bird_v2_quantity, na.rm = TRUE),
    PV_quantity = mean(bird_pv_quantity, na.rm = TRUE),
    SE_quantity = mean(se_quantity, na.rm = TRUE),
    
    lower_ci_pv_quantity = mean(pv_lower_ci_quantity),
    upper_ci_pv_quantity = mean(pv_upper_ci_quantity),
    lower_ci_cv_quantity = mean(cv_lower_ci_quantity),
    upper_ci_cv_quantity = mean(cv_upper_ci_quantity),
    lower_ci_v2_quantity = mean(v2_lower_ci_quantity),
    upper_ci_v2_quantity = mean(v2_upper_ci_quantity),
    
    # Number of individuals
    n_individuals = n_distinct(band_number),
    n_inf = n_distinct(band_number[inf_pri == 1]),
    n_uninf = n_distinct(band_number[inf_pri == 0]),
    mean_tes = mean(mean_tes),
    mean_quantity = mean(mean_quantity)
  ) %>%
  as_tibble()

print(summary_tibble.p.v)

#Primary Values 
p.cv <- summary_tibble.p.v[, c(1, 2, 6, 9, 10, 27, 13, 19, 20, 28, 17)]
p.pv <- summary_tibble.p.v[, c(1, 5, 7, 8, 27, 6, 16, 18, 19, 28, 17)]
p.var <- summary_tibble.p.v %>% 
  dplyr::select(primary_treatment, CV_tes, PV_tes, mean_tes, SE_tes, CV_quantity, PV_quantity, mean_quantity, SE_quantity, n_inf)
