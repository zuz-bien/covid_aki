# Sensitivity analysis - patients with baseline SCr only

all_patients <- read_xlsx("all_pt.xlsx")

all_pt <- all_patients %>% filter (is.na(cr_base) == FALSE)

all_pt$EthnicGroup2 <- factor (all_pt$EthnicGroup2, levels = c("White", "Asian", "Black", "Other", "Unknown"))
all_pt$quintile_study <- factor (all_pt$quintile_study, levels = c( "1", "2", "3", "4", "5"), labels = c("1", "2", "3", "4", "5"))
all_pt$stage <- factor(all_pt$stage, levels = c(0, 1, 2, 3), labels = c(0, 1, 2, 3)) 

median(all_pt$AdmissionAge)

# TABLE ONE 

myVars <- c("AdmissionAge", "Sex", "EthnicGroup2", "quintile_study", "smoking", "rend", 
            "obes","ihd", "ami", "chf", "Diabetes", "HTN","pvd", "cevd", "copd", "ld","dementia", "canc",
            "cci_bin", "rfscat", "hfrscat", "base_gfr","low_gfr_base","last_gfr", "low_gfr_final", 
            "crp", "crp_bin","icu", "icumv", "rrt", "iculos", "icuorgan", "hosp_days", 
            "diednew", "died_daysnew", "died30new", "died90new","make90", "dis_place")

catVars <- c("Sex", "EthnicGroup2", "quintile_study", "smoking",
             "bmicat2", "obes", "rend","ihd", "ami", "chf", 
             "Diabetes", "HTN","pvd", "cevd", "copd", "ld","dementia", "canc", 
             "cci_bin", "rfscat", "hfrscat", "crp_bin","low_gfr_base", "low_gfr_final",
             "icu", "icumv", "rrt", "icuorgan", "make90",
             "dis_place", "diednew", "died30new", "died90new")

tab <- CreateTableOne(vars = myVars, strata = "stage", data = all_pt, factorVars = catVars)

write.csv(print(tab, nonnormal = c("AdmissionAge", "baseline_SCr_mgdl", "last_SCr_mgdl", "base_gfr", "last_gfr", "died_daysnew", "hosp_days", "iculos", "crp"), formatOptions = list(big.mark = ",")), "tableone_baseline.csv")

write_xlsx(as.data.frame(sapply (all_pt %>% select (myVars), function (x) length(x) - sum(is.na(x)))) %>% rownames_to_column() %>% rename (var = 1, n = 2), "n_tab_baseline.xlsx")

# MAKE 90 TABLE

make90 <- all_pt %>% filter(make90 == 1)

make90$survunrec <- case_when (make90$died90new == 0 & make90$unrecovery == 1 ~ 1)
make90$survunrec <- replace_na(make90$survunrec, 0)
make90$survrrt <- case_when (make90$died90new == 0 & make90$rrt == 1 ~ 1)
make90$survrrt <- replace_na(make90$survrrt, 0)

makeVars <- c("survrrt", "survunrec")
catVars<- c("survrrt", "survunrec")

tab_make90sens <- CreateTableOne(vars = makeVars, factorVars = catVars, strata = "stage", data = make90)
write.csv(print(tab_make90sens, formatOptions = list(big.mark = ",")), "tab_make90sense.csv")

# UNIVARIATE 30

surv_df <- all_pt %>% 
  mutate (censor = ifelse(died30new == 0, 30, NA)) %>%
  mutate (futime = coalesce(censor, died_daysnew))

# adjusted for age and sex 

age_sex_adjusted <- coxph(Surv(futime, died30new) ~ stage + age_div + Sex, data = surv_df)

stage_df <- with(surv_df, 
                 data.frame(stage = factor(c(0, 1, 2, 3)), 
                            age_div = c(rep(median(age_div), 4)),
                            Sex = c(rep("Male", 4))))

fit <- survfit(age_sex_adjusted, newdata = stage_df)  

summary(age_sex_adjusted)

age_sex_adjusted_surv_plot <- ggsurvplot(fit, data= surv_df, legend.title="", conf.int = TRUE, legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), ggtheme = theme_minimal(), break.time.by = 5)

age_sex_adjusted_surv_plot


png("uni_30_baseline_cox.png", units="cm", width=11.4, height=10.0, res=400)
age_sex_adjusted_surv_plot
dev.off()

# forest plot

age_sex_x <- c("stage1" = "AKI Stage 1", "stage2" = "AKI Stage 2", "stage3" = "AKI Stage 3",
               "age_div" = "Age", "SexMale" = "Male sex")

age_sex_round <- HR_extractor(age_sex_adjusted)

age_sex_round$predictor <- factor(age_sex_round$predictor, levels = rev(c("SexMale", "age_div","stage1", "stage2", "stage3")))

age_sex_fp <- make_log_forest (age_sex_round, breaks = 1, xlabs = age_sex_x, xcoords = 0.75,
                               a = 1, b  = 3.2, c = 7, alpha = 0.3, barwidth = 0.2)

png("uni_30_baseline_fp.png", units="cm", width=17.1, height=12.1, res=400)
age_sex_fp
dev.off()



# MULTIVARIATE 30

all_covariates <- coxph(Surv(futime, died30new) ~ stage + age_div + Sex + EthnicGroup2 + quintile_study + smoking + obes + Diabetes + HTN + low_gfr_baseline, data = surv_df)

stage_df <- with(surv_df, 
                 data.frame(stage = factor(c(0, 1, 2, 3)), 
                            age_div = rep(median(age_div), 4), 
                            Sex = c(rep("Male", 4)),
                            EthnicGroup2 = c(rep("White", 4)), 
                            quintile_study = factor(c(rep(1, 4))), 
                            smoking = c(rep(0, 4)), 
                            obes  = c(rep(0, 4)), 
                            Diabetes = c(rep(0, 4)), 
                            HTN = c(rep(0, 4)), 
                            low_gfr_baseline = c(rep(0, 4))))

fit <- survfit(all_covariates, newdata = stage_df)     

all_covariates_surv_plot <- ggsurvplot(fit, data= all_patients,legend.title="", conf.int = TRUE, legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), ggtheme = theme_minimal(), ylim = c(0, 1))

all_covariates_surv_plot

png("multi_30_baseline_surv.png", units="cm", width=11.4, height=10.0, res=400)
all_covariates_surv_plot 
dev.off()

summary(all_covariates)

# Multivariate 30 forest plot

all_cov_x <- c("stage1" = "AKI Stage 1", "stage2" = "AKI Stage 2", "stage3" = "AKI Stage 3",
               "EthnicGroup2Asian" = "Ethnicity: Asian", "EthnicGroup2Black" = "Ethnicity: Black",
               "EthnicGroup2Other" = "Ethnicity: Other", "EthnicGroup2Unknown" = "Ethnicity: Unknown",
               "SexMale" = "Male sex", "smoking" = "Smoking", "obes" = "Obesity", "low_gfr_baseline" = "CKD",
               "age_div" = "Age",
               # "age_quint1" = "Age 16-47", "age_quint2" = "Age 48-59", "age_quint4" = "Age 60-70", "age_quint5" = "Age 71+",
               "quintile_study5" = "IMD 5", "quintile_study4" = "IMD 4", "quintile_study3" = "IMD 3", "quintile_study2" = "IMD 2") 

all_covariates_round <- HR_extractor(all_covariates)

all_covariates_round$predictor <- factor(all_covariates_round$predictor, levels <- 
                                           rev(c("stage1", "stage2", "stage3" ,
                                                 "SexMale", "age_div",
                                                 "EthnicGroup2Asian", "EthnicGroup2Black",  "EthnicGroup2Other", "EthnicGroup2Unknown", 
                                                 "quintile_study5" ,"quintile_study4" , "quintile_study3", "quintile_study2" , 
                                                 "smoking" , "obes", "low_gfr_baseline", "HTN", "Diabetes"))
)                                             

multi_30_fp <- make_log_forest (all_covariates_round, breaks = 1, title = " ", xlabs = all_cov_x, a = 2, b = 6, c = 15)

png("multi_30_baseline_fp.png", units="cm", width=17.1, height=17.1, res=400)
multi_30_fp 
dev.off()


# CREATE SURV 90 DATAFRAME (90 day mortality)


surv_90 <- all_pt %>% 
  mutate (censor = ifelse(died90new == 0, 90, NA)) %>%
  mutate (futime = coalesce(censor, died_daysnew)) 



# UNIVARIATE 90 

#Cox 

age_sex_90 <- coxph(Surv(futime, died90new) ~ stage + age_div + Sex, data = surv_90)

summary (age_sex_90)

stage_90 <- with(surv_90, 
                 data.frame(stage = factor(c(0, 1, 2, 3)), 
                            age_div = c(rep(median(age_div), 4)),
                            Sex = c(rep("Male", 4))))

fit <- survfit(age_sex_90, newdata = stage_90)  

age_sex_90_surv_plot <- ggsurvplot(fit, data= surv_90, legend.title="", conf.int = TRUE, legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), ggtheme = theme_minimal(), break.time.by = 15)

age_sex_90_surv_plot

png("uni_90_baseline_surv.png", units="cm", width=11.4, height=10.0, res=400)
age_sex_90_surv_plot 
dev.off()

# forest plot


age_sex_x <- c("stage1" = "AKI Stage 1", "stage2" = "AKI Stage 2", "stage3" = "AKI Stage 3",
               "age_div" = "Age", "SexMale" = "Male sex")

age_sex_90_round <- HR_extractor(age_sex_90)

age_sex_90_round$predictor <- factor(age_sex_90_round$predictor, levels = rev(c("SexMale", "age_div","stage1", "stage2", "stage3")))

age_sex_90_fp <- make_log_forest (age_sex_90_round, breaks = 1, xlabs = age_sex_x, xcoords = 0.75,
                                  a = 1, b  = 3.2, c = 7, alpha = 0.3, barwidth = 0.2)

png("uni_90_baseline_fp.png", units="cm", width=17.1, height=12.1, res=400)
age_sex_90_fp 
dev.off()


# MULTIVARIATE 90

add_cov <- c("EthnicGroup2","IMD", "smoking", "obes", "Diabetes", "HTN", "low_gfr_base")

all_cov_90 <- coxph(Surv(futime, died90new) ~ stage + age_div + Sex + EthnicGroup2 + quintile_study + smoking + obes + Diabetes + HTN + low_gfr_baseline, data = surv_90)

summary(all_cov_90)

stage_df <- with(surv_90, 
                 data.frame(stage = factor(c(0, 1, 2, 3)), 
                            age_div = rep(median(age_div), 4),
                            Sex = c(rep("Male", 4)),
                            EthnicGroup2 = c(rep("White", 4)), 
                            quintile_study = factor(c(rep(1, 4))), 
                            smoking = c(rep(0, 4)), 
                            obes  = c(rep(1, 4)), 
                            Diabetes = c(rep(0, 4)), 
                            HTN = c(rep(0, 4)), 
                            low_gfr_baseline = c(rep(0, 4))))

fit <- survfit(all_cov_90, newdata = stage_df)     

all_cov_90_surv_plot <- ggsurvplot(fit, data= surv_90, legend.title="", conf.int = TRUE, legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), ggtheme = theme_minimal(), break.time.by = 15)

all_cov_90_surv_plot

png("multi_90_baseline_surv.png", units="cm", width=11.4, height=10.0, res=400)
all_cov_90_surv_plot
dev.off()

# forest plot 


all_cov_x <- c("stage1" = "AKI Stage 1", "stage2" = "AKI Stage 2", "stage3" = "AKI Stage 3",
               "EthnicGroup2Asian" = "Ethnicity: Asian", "EthnicGroup2Black" = "Ethnicity: Black",
               "EthnicGroup2Other" = "Ethnicity: Other", "EthnicGroup2Unknown" = "Ethnicity: Unknown",
               "SexMale" = "Male sex", "smoking" = "Smoking", "obes" = "Obesity", "low_gfr_baseline" = "CKD",
               "age_div" = "Age",
               # "age_quint1" = "Age 16-47", "age_quint2" = "Age 48-59", "age_quint4" = "Age 60-70", "age_quint5" = "Age 71+",
               "quintile_study5" = "IMD 5", "quintile_study4" = "IMD 4", "quintile_study3" = "IMD 3", "quintile_study2" = "IMD 2") 

all_covariates_round <- HR_extractor(all_cov_90)

all_covariates_round$predictor <- factor(all_covariates_round$predictor, levels <- 
                                           rev(c("stage1", "stage2", "stage3" ,
                                             "SexMale", "age_div",
                                             "EthnicGroup2Asian", "EthnicGroup2Black",  "EthnicGroup2Other", "EthnicGroup2Unknown", 
                                                 "quintile_study5" ,"quintile_study4" , "quintile_study3", "quintile_study2" , 
                                                 "smoking" , "obes", "low_gfr_baseline", "HTN", "Diabetes"))
)                                             

multi_90_fp <- make_log_forest (all_covariates_round, breaks = 1, title = " ", xlabs = all_cov_x, a = 2, b = 6, c = 15)

png("multi_90_baseline_fp.png", units="cm", width=17.1, height=17.1, res=400)
multi_90_fp 
dev.off()


## AKI RISK FACTORS 

glm.baseline <- glm(aki_bin ~  
                  age_div +  
                  Sex +
                  EthnicGroup2 + 
                  quintile_study + 
                  smoking + 
                  # bmicat2 +
                  obes +
                  # rend + 
                    low_gfr_baseline +
                  chf + 
                  Diabetes + 
                  HTN + 
                  cevd + 
                  crp_bin
                ,na.action = na.omit,
                data = all_pt, family = binomial)



OR_baseline <- OR_extractor(glm.baseline)

xaxis <- c("age_div" = "Age", "SexMale" = "Male sex", 
           "EthnicGroup2Black" = "Ethnicity: Black", "EthnicGroup2Unknown" = "Ethnicity: Unknown",
           "EthnicGroup2Other" = "Ethnicity: Other", "EthnicGroup2Asian" = "Ethnicity: Asian", 
           "quintile_study2" = "IMD 2", "quintile_study3" = "IMD 3", "quintile_study4" = "IMD 4",
           "quintile_study5" = "IMD 5", "smoking" = "Smoking", "obes" = "Obesity", "low_gfr_baseline" = "CKD", 
           "chf" = "CHF",  "Diabetes" = "Diabetes", "HTN" = "HTN", "cevd" = "CEVD", "crp_bin" = "CRP > 145")

OR_baseline$predictor <- factor(OR_baseline$predictor, levels <- 
                                        rev(c("SexMale", "age_div",
                                              "EthnicGroup2Asian", "EthnicGroup2Black",  "EthnicGroup2Other", "EthnicGroup2Unknown", 
                                              "quintile_study5" ,"quintile_study4" , "quintile_study3", "quintile_study2" , 
                                              "smoking" , "obes", "low_gfr_baseline", "HTN", "Diabetes","cevd", "chf", "crp_bin")))




aki_rf_baseline <- make_log_forest(OR_baseline, breaks =1, xlabs = xaxis, title = "", a = 1, b = 5, c = 12, metric = "OR")

png("aki_rf_baseline.png", units="cm", width=17.1, height=17.1, res=400)
aki_rf_baseline
dev.off()


glm.baseline2 <- lrm(aki_bin ~  
                      age_div +  
                      Sex +
                      EthnicGroup2 + 
                      quintile_study + 
                      smoking + 
                      # bmicat2 +
                      obes +
                      # rend + 
                       low_gfr_baseline +
                      chf + 
                      Diabetes + 
                      HTN + 
                      cevd + 
                      crp_bin, data = all_pt)

glm.baseline2


