# 90 day mortality 
all_patients <- read_xlsx("all_pt.xlsx")

all_patients$stage <- factor(all_patients$stage, levels = c(0, 1, 2, 3))
all_patients$EthnicGroup2 <- factor(all_patients$EthnicGroup2, levels = c("White", "Asian", "Black", "Other", "Unknown"))
all_patients$quintile_study <- factor (all_patients$quintile_study, levels = c( "1", "2", "3", "4", "5"), labels = c("1", "2", "3", "4", "5"))


surv_90 <- all_patients %>% 
  mutate (censor = ifelse(died90new == 0, 90, NA)) %>%
  mutate (futime = coalesce(censor, died_daysnew)) 


# adjusted for age and sex 

age_sex_90 <- coxph(Surv(futime, died90new) ~ stage + age_div + Sex, data = surv_90)

summary (age_sex_90)

stage_90 <- with(surv_90, 
                 data.frame(stage = factor(c(0, 1, 2, 3)), 
                            age_div = c(rep(median(age_div), 4)),
                            Sex = c(rep("Male", 4))))

fit <- survfit(age_sex_90, newdata = stage_90)  

age_sex_90_surv_plot <- ggsurvplot(fit, data= surv_90, legend.title="",conf.int = TRUE, legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), ggtheme = theme_minimal(), break.time.by = 15)

age_sex_90_surv_plot

png("uni_90_surv.png", units="cm", width=11.4, height=10.0, res=400)
age_sex_90_surv_plot 
dev.off()

# number at risk 

age_sex_90_strat <- coxph(Surv(futime, died90new) ~ strata(stage) + age_div + Sex, data = surv_90)

fit <- survfit(age_sex_90_strat, newdata = stage_90)

risk_tab <- ggrisktable(fit, surv_90, break.time.by = 15)
age_sex_90_rt <- risk.table(risk_tab, 0, 90, 15) %>% mutate (
  no_AKI = paste (no_AKI, "(", no_AKI[1] - no_AKI, ")", sep = ""), 
  stage_1 = paste (stage_1, "(", stage_1[1] - stage_1, ")", sep = ""), 
  stage_2 = paste (stage_2, "(", stage_2[1] - stage_2, ")", sep = ""), 
  stage_3 = paste (stage_3, "(", stage_3[1] - stage_3, ")", sep = ""), 
)

uni_90_rt <- age_sex_90_rt %>% gather ( key = "stage", value = "n", no_AKI:stage_3) %>% 
  spread(key = timepoint, value = n)

write_xlsx(uni_90_rt, "uni_90_rt.xlsx")

# make log forest 

age_sex_x <- c("stage1" = "AKI Stage 1", "stage2" = "AKI Stage 2", "stage3" = "AKI Stage 3",
               "age_div" = "Age", "SexMale" = "Male sex")

age_sex_90_round <- HR_extractor(age_sex_90)

age_sex_90_round$predictor <- factor(age_sex_90_round$predictor, levels = rev(c("SexMale", "age_div","stage1", "stage2", "stage3")))

age_sex_90_fp <- make_log_forest (age_sex_90_round, breaks = 1, xlabs = age_sex_x, xcoords = 0.75,
                               a = 1, b  = 3.2, c = 7, alpha = 0.3, barwidth = 0.2)

png("uni_90_fp.png", units="cm", width=17.1, height=12.1, res=400)
age_sex_90_fp 
dev.off()


### all covariates 


surv_90 <- rename(surv_90, "IMD" = quintile_study)
surv_90$IMD <- factor (surv_90$IMD, levels = c(1:5))

add_cov <- c("EthnicGroup2","IMD", "smoking", "obes", "Diabetes", "HTN", "ckd_gfr")

all_cov_90 <- coxph(Surv(futime, died90new) ~ stage + age_div + Sex + EthnicGroup2 + IMD + smoking + obes + Diabetes + HTN + ckd_gfr, data = surv_90)

summary(all_cov_90)

stage_df <- with(surv_90, 
                 data.frame(stage = factor(c(0, 1, 2, 3)), 
                            AdmissionAge = rep(median(age_div), na.rm = TRUE, 4),
                            Sex = c(rep("Male", 4)),
                            EthnicGroup2 = c(rep("White", 4)), 
                            IMD = factor(c(rep(1, 4))), 
                            smoking = c(rep(0, 4)), 
                            obes  = c(rep(1, 4)), 
                            Diabetes = c(rep(0, 4)), 
                            HTN = c(rep(0, 4)), 
                            ckd_gfr = c(rep(0, 4))))

fit <- survfit(all_cov_90, newdata = stage_df)     

all_cov_90_surv_plot <- ggsurvplot(fit, data= surv_90, legend.title="", conf.int = TRUE, legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), ggtheme = theme_minimal(), break.time.by = 15)

all_cov_90_surv_plot

#at risk table 
all_cov_90_strat <- coxph(Surv(futime, died90new) ~ strata(stage) + AdmissionAge + Sex + EthnicGroup2 + IMD + smoking + obes + Diabetes + HTN + ckd_gfr, data = surv_90)
fit <- survfit(all_cov_90_strat, newdata = stage_df) 
risk_tab <- ggrisktable(fit, surv_90, break.time.by = 15)
all_cov_90_rt <- risk.table(risk_tab, 0, 90, 15)
write.csv(all_cov_90_rt %>% mutate (deaths_no_AKI = no_AKI[1] - no_AKI, deaths_1 = stage_1[1] - stage_1, deaths_2 = stage_2[1] - stage_2, deaths_3 = stage_3[1] - stage_3) %>% select (timepoint, no_AKI, deaths_no_AKI, stage_1, deaths_1, stage_2, deaths_2, stage_3, deaths_3), "all_cov_90_rt.csv")



# forest plot 


covariate_names <- c(stage = "Stage" , AdmissionAge = "Age", Sex = "Sex", Ethnicity = "Ethnicity", IMD = "IMD", smoking = "Smoking", obes = "Obesity", Diabetes = "Diabetes", HTN = "HTN", ckd_gfr = "CKD")

surv_90 %>% rename("Ethnicity" = EthnicGroup2, "Stage" = stage)  %>% analyse_multivariate(vars(futime, died90new), covariates = vars(Stage, AdmissionAge, Sex, Ethnicity, IMD, smoking, obes, Diabetes, HTN, ckd_gfr), 
                                                                                          covariate_name_dict = covariate_names) -> all_cov_result_90
print(all_cov_result_90)

all_cov_90_fp <- forest_plot(all_cov_result_90, factor_labeller = covariate_names, labels_displayed = c("factor"), factor_id_sep = ": ")

all_cov_90_fp

# kaplan meier

kpm <- survfit(Surv(futime, died90new) ~ stage, data = surv_90)

kpm_90 <- ggsurvplot(kpm, data = surv_90, conf.int = TRUE, legend.title = " ",legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), ggtheme = theme_minimal(), break.time.by = 15)

png("kpm_90.png", units="cm", width=11.4, height=10.0, res=400)
kpm_90
dev.off()
