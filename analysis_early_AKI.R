setwd("/Users/Zuzanna_Bien/Data/Files")

all_patients <- read_xlsx("all_patients_early_AKI.xlsx")

# find SCr closest to day 90 for eGFR calculations (as discussed with John, he wanted 
# the value higher than but closest to day 90 if possible; and if not, then value lower
# than and closest to D90 (hence why this is slightly convoluted))

last_SCr_after_90 <- bloods %>% filter (ID %in% all_patients$ID) %>%
  filter (result_day >= -1) %>% 
  mutate(ClinicalEventResult = as.numeric(ClinicalEventResult)) %>%
  na.omit() %>%
  arrange (ID, result_day)%>% 
  filter (result_day >= 90) %>% 
  group_by(ID) %>% filter (abs(result_day - 90) == min(abs(result_day - 90))) %>%
  filter (ClinicalEventResult == max(ClinicalEventResult)) %>% slice(1L)
  
last_SCr_before_90 <-  bloods %>% filter (ID %in% all_patients$ID, !ID %in% last_SCr_after_90$ID) %>%
  filter (result_day >= -1) %>% 
  mutate(ClinicalEventResult = as.numeric(ClinicalEventResult)) %>%
  na.omit() %>%
  arrange (ID, result_day)%>% 
  group_by(ID) %>% filter (abs(result_day - 90) == min(abs(result_day - 90))) %>%
  filter (ClinicalEventResult == max(ClinicalEventResult)) %>% slice(1L)

last_SCr <- bind_rows(last_SCr_before_90, last_SCr_after_90)

missing_SCr_ID <- setdiff(patients_1873$ID, last_SCr$ID)

# here is the final number of patients who have had any SCr results recorded between day -1 and
# the end of follow up period, this is where number 1855 comes from 
all_patients_1855 <- merge(all_patients, last_SCr %>% select (ID, last_SCr = ClinicalEventResult))

# dataframe construction

all_pt <- all_patients_1855   %>% mutate(baseline_SCr = coalesce(cr_base, cr_imp), 
                                  crp_bin = ifelse(crp > 145, 1, 0),
                                  aki_bin = ifelse(stage == 0, 0, 1),
                                  cci_bin = case_when (
                                    cci ==0 ~ "0",
                                    cci %in% c(1, 2) ~"1-2",
                                    cci %in% c(3, 4) ~"3-4",
                                    cci >= 5 ~"5+"), 
                                  age_div = AdmissionAge/10,
                                  cr_base_mgdl = cr_base/88.4,
                                  baseline_SCr_mgdl = baseline_SCr/88.4, 
                                  last_SCr_mgdl = last_SCr/ 88.4, 
                                  black = ifelse (EthnicGroup2 == "Black", TRUE, FALSE))

#from now on, "baseline SCr/ gfr" will be used to refer to both imputed and measured values
#whereas "base SCr/gfr " is just looking at measured values 

#calculate baseline egfr for all measured and imputed values of baseline SCr as the function can't
#handle NAs

all_pt$baseline_gfr <- mapply(ckdepicre,age=all_pt$AdmissionAge,sex=all_pt$Sex=="Male",cre=all_pt$baseline_SCr_mgdl, race=all_pt$black) 

#replace the base eGFR value with NA for patients who had imputed baseline SCr
all_pt$base_gfr <- mapply(ckdepicre,age=all_pt$AdmissionAge,sex=all_pt$Sex=="Male",cre=all_pt$baseline_SCr_mgdl, race=all_pt$black) 
all_pt$base_gfr <- ifelse(is.na(all_pt$cr_base) == TRUE, NA, all_pt$base_gfr)

#calculate last eGFR
all_pt$last_gfr <-  mapply(ckdepicre,age=all_pt$AdmissionAge,sex=all_pt$Sex=="Male",cre=all_pt$last_SCr_mgdl, race=all_pt$black)

#'unrecovery' = patients whose last eGFR was <0.7 baseline 
all_pt$unrecovery <- ifelse(all_pt$last_gfr/all_pt$gfr < 0.7, 1, 0)

all_pt$make90 <- with(all_pt, 
                      case_when (
                        unrecovery == 1 ~ 1,
                        rrt == 1 ~ 1, 
                        died90new == 1 ~ 1
                      ))

all_pt$make90krt <- with(all_pt, case_when(make90 == 1 & died90new == 0 & rrt == 1 ~ 1))
all_pt$make90unrecovery <- with(all_pt, case_when(make90 == 1 & died90new == 0 & unrecovery == 1 ~ 1))
all_pt$delta <- all_pt$last_SCr_mgdl - all_pt$baseline_SCr_mgdl

#base = only for people with cr values measures
#baseline = for those with measured AND imputed scrs 
all_pat <- all_pt %>% mutate (make90 = ifelse(is.na(make90) == TRUE, 0, 1),
                              low_gfr_final = ifelse (last_gfr <60, 1, 0), 
                              low_gfr_baseline = ifelse(baseline_gfr < 60, 1, 0),
                              low_gfr_base = ifelse(base_gfr < 60, 1, 0))

write_xlsx(all_pat, "all_pt.xlsx")

# TABLE ONE

all_pt <- read_xlsx("all_pt.xlsx")

myVars <- c("AdmissionAge", "Sex", "EthnicGroup2", "quintile_study", "smoking", "rend", 
            "obes","ihd", "ami", "chf", "Diabetes", "HTN","pvd", "cevd", "copd", "ld","dementia", "canc",
             "cci_bin", "rfscat", "hfrscat", "base_gfr","low_gfr_base","last_gfr", "low_gfr_final", 
            "crp", "crp_bin","icu", "icumv", "rrt", "iculos", "icuorgan", "hosp_days", 
             "diednew", "died_daysnew", "died30new", "make90", "died90new", "dis_place")

catVars <- c("Sex", "EthnicGroup2", "quintile_study", "smoking",
             "bmicat2", "obes", "rend","ihd", "ami", "chf", 
             "Diabetes", "HTN","pvd", "cevd", "copd", "ld","dementia", "canc", 
             "cci_bin", "rfscat", "hfrscat", "crp_bin","low_gfr_base", "low_gfr_final",
             "icu", "icumv", "rrt", "icuorgan", "make90", 
             "dis_place", "diednew", "died30new", "died90new")

tab_bin <- CreateTableOne(vars = myVars, strata = "aki_bin", data = all_pt, factorVars = catVars, addOverall = TRUE)

write.csv(print(tab_bin, nonnormal = c("AdmissionAge", "baseline_SCr_mgdl", "last_SCr_mgdl", "base_gfr", "last_gfr", "died_daysnew", "hosp_days", "iculos", "crp"), formatOptions = list(big.mark = ",")), "tableone_bin.csv")

write_xlsx(as.data.frame(sapply (all_pt %>% select (myVars), function (x) length(x) - sum(is.na(x)))) %>% rownames_to_column() %>% rename (var = 1, n = 2), "n_tab_bin.xlsx")

# TABLE ONE - FOUR GROUPS 

tab <- CreateTableOne(vars = myVars, strata = "stage", data = all_pt, factorVars = catVars)

write.csv(print(tab, nonnormal = c("AdmissionAge", "baseline_SCr_mgdl", "last_SCr_mgdl", "base_gfr", "last_gfr", "died_daysnew", "hosp_days", "iculos", "crp"), formatOptions = list(big.mark = ",")), "tableone.csv")

# TABLE ONE - MAKE 90

make90 <- all_pt %>% filter(make90 == 1)

make90$survunrec <- case_when (make90$died90new == 0 & make90$unrecovery == 1 ~ 1)
make90$survunrec <- replace_na(make90$survunrec, 0)
make90$survrrt <- case_when (make90$died90new == 0 & make90$rrt == 1 ~ 1)
make90$survrrt <- replace_na(make90$survrrt, 0)

makeVars <- c("survrrt", "survunrec")
catVars<- c("survrrt", "survunrec")

tab_make90krt <- CreateTableOne(vars = makeVars, factorVars = catVars, strata = "aki_bin", data = make90, addOverall = TRUE)
write.csv(print(tab_make90krt, formatOptions = list(big.mark = ",")), "tab_make90krt.csv")

tab_make90stage <- CreateTableOne(vars = makeVars, factorVars = catVars, strata = "stage", data = make90)
write.csv(print(tab_make90stage, formatOptions = list(big.mark = ",")), "tab_make90stage.csv")



# SURVIVAL ANALYSIS 

## dataframe prep 
all_pt <- read_xlsx("all_pt.xlsx")

all_pt$EthnicGroup2 <- factor (all_pt$EthnicGroup2, levels = c("White", "Asian", "Black", "Other", "Unknown"))
all_pt$quintile_study <- factor (all_pt$quintile_study, levels = c( "1", "2", "3", "4", "5"), labels = c("1", "2", "3", "4", "5"))
all_pt$stage <- factor(all_pt$stage, levels = c(0, 1, 2, 3), labels = c(0, 1, 2, 3)) 
all_pt$aki_bin <- factor(all_pt$aki_bin, levels = c(0, 1))


surv_df <- all_pt %>% 
  mutate (censor = ifelse(died30new == 0, 30, NA)) %>%
  mutate (futime = coalesce(censor, died_daysnew))

# Age and sex adjusted Cox plot

age_sex_adjusted <- coxph(Surv(futime, died30new) ~ stage + age_div + Sex, data = surv_df)

stage_df <- with(surv_df, 
                 data.frame(stage = factor(c(0, 1, 2, 3)), 
                            age_div = c(rep(median(age_div), 4)),
                            Sex = c(rep("Male", 4))))

fit <- survfit(age_sex_adjusted, newdata = stage_df)  

summary(age_sex_adjusted)

age_sex_adjusted_surv_plot <- ggsurvplot(fit, data= surv_df, conf.int = TRUE, legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), legend.title="", ggtheme = theme_minimal(), break.time.by = 5)

age_sex_adjusted_surv_plot


png("uni_30_cox.png", units="cm", width=11.4, height=10.0, res=400)
age_sex_adjusted_surv_plot
dev.off()


# Age and sex-adjusted forest plot
age_sex_x <- c("stage1" = "AKI Stage 1", "stage2" = "AKI Stage 2", "stage3" = "AKI Stage 3",
"age_div" = "Age", "SexMale" = "Male sex")

age_sex_round <- HR_extractor(age_sex_adjusted)

age_sex_round$predictor <- factor(age_sex_round$predictor, levels = rev(c("SexMale", "age_div","stage1", "stage2", "stage3")))
                                                 
age_sex_fp <- make_log_forest (age_sex_round, breaks = 1, xlabs = age_sex_x, xcoords = 0.75,
                               a = 1, b  = 3.2, c = 7, alpha = 0.3, barwidth = 0.2)

png("uni_30_fp.png", units="cm", width=17.1, height=12.1, res=400)
age_sex_fp
dev.off()


## risk table 

risk.table <- function(df, t0, tmax, breaks) {
  a <- as.data.frame (matrix (df$data$n.risk, ncol = 4), row.names = c(seq(t0, tmax, by = breaks))) %>% rownames_to_column()
  colnames (a) <- c("timepoint", "no_AKI", "stage_1", "stage_2", "stage_3")
  a
}


age_sex_adjusted_strat <- coxph(Surv(futime, died30new) ~ strata(stage) + age_div + Sex, data = surv_df)

fit <- survfit(age_sex_adjusted_strat, newdata = stage_df)  

fit$n.event

risk_tab <- ggrisktable(fit, surv_df, break.time.by = 5)

uni_30 <- risk.table(risk_tab, 0, 31, 5) %>% mutate (
  no_AKI = paste (no_AKI, "(", no_AKI[1] - no_AKI, ")", sep = ""), 
  stage_1 = paste (stage_1, "(", stage_1[1] - stage_1, ")", sep = ""), 
  stage_2 = paste (stage_2, "(", stage_2[1] - stage_2, ")", sep = ""), 
  stage_3 = paste (stage_3, "(", stage_3[1] - stage_3, ")", sep = ""), 
) 

uni_30_rt <- uni_30 %>% gather ( key = "stage", value = "n", no_AKI:stage_3) %>% 
  spread(key = timepoint, value = n) 

#careful - patients who died on the last day are misclassified (use this table to inspect)
#might need adjusting manually
data.frame(fit$n.risk, fit$n.censor, fit$n.event)

write_xlsx(uni_30_rt, "uni_30_rt.xlsx")


## Ethnicity x AKI interaction 

interaction_simple <- coxph(Surv(futime, died30new) ~ aki_bin*EthnicGroup2 + age_div + Sex, data = surv_df)

anova(interaction_simple)

## Multivariate analysis  

# Kaplan meier plot

surv_object <- Surv(time= surv_df$futime, event = surv_df$died30new) 
kpm <- survfit(Surv(futime, died30new) ~ stage, data = surv_df)

multi_30_kpm <- ggsurvplot(kpm, data = surv_df, pval = FALSE)

png("multi_30_kpm.png", units="cm", width=11.4, height=10.0, res=400)
multi_30_kpm
dev.off()

            
# Multivariate Cox plot 

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

all_covariates_surv_plot <- ggsurvplot(fit, data= all_patients, legend.title="", conf.int = TRUE, legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), ggtheme = theme_minimal(), ylim = c(0, 1))

all_covariates_surv_plot

png("multi_30_surv.png", units="cm", width=11.4, height=10.0, res=400)
all_covariates_surv_plot 
dev.off()

summary(all_covariates)

# Side-by-side Cox and Kaplan-Meier plot

require("survminer")
splots <- list()
splots[[2]] <- ggsurvplot(fit, data= surv_df, conf.int = TRUE,legend = "top",  legend.labs = c(" ", "  ", "   ", "    "),  ggtheme = theme_minimal())
splots[[1]] <- ggsurvplot(kpm, data = surv_df, conf.int = TRUE, legend.title = " ",legend.labs = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), ggtheme = theme_minimal())

# Arrange multiple ggsurvplots and print the output
p <- arrange_ggsurvplots(splots, print = TRUE,
                    ncol = 2, nrow = 1 )


png("multi_30_compare.png", units="cm", width=17.1, height=12.0, res=400)
p
dev.off()

# Multivariate forest plot 

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

png("multi_30_fp.png", units="cm", width=17.1, height=17.1, res=400)
multi_30_fp 
dev.off()



## Figures S2A-E

setwd("~/Data/Files")

all_pt <- read_xlsx("all_pt.xlsx")

all_pt$stage <- factor(all_pt$stage, levels = c(0, 1, 2, 3))

all_pt_ID <- all_pt$ID

#only patients with bloods during admission period
admission_SCr <- bloods %>% filter (ID %in% all_pt_ID, result_day >= -1, ClinicalEvent == "Creatinine Serum") %>% 
  mutate(ClinicalEventResult = as.numeric (ClinicalEventResult)) %>% na.omit() %>%
  arrange(ID, result_day) %>% group_by(ID) %>% slice_head() %>% select(ID, first_SCr =ClinicalEventResult)

a <- data.frame()

for (i in all_pt$ID) {
  admission <- all_pt %>% filter (ID == i) %>% select (hosp_days) %>% unlist (use.names = FALSE)
  admission_bloods <- bloods %>% filter (ID == i, result_day >= -1, result_day <= admission)
  a <- rbind(a, admission_bloods)
}

#patients who did not have any blood tests during their admission 
no_bloods_during_adm <- setdiff(all_pt$ID, a$ID)

#find peak value during admission 
peak_SCr <- a %>% 
  group_by(ID) %>% mutate(ClinicalEventResult = as.numeric (ClinicalEventResult)) %>% 
  na.omit() %>% summarise (peak_SCr = max(ClinicalEventResult))


creat <- merge(admission_SCr, peak_SCr, by = "ID", all = TRUE)

#find last SCr value on admission (i.e. the discharge value)
discharge_SCr <- a %>% 
  group_by(ID) %>% mutate(ClinicalEventResult = as.numeric (ClinicalEventResult)) %>% 
  na.omit() %>% arrange(ID, result_day) %>% slice_tail() %>% select(ID, discharge_SCr =ClinicalEventResult)


creats <- merge(creat, discharge_SCr, by = "ID", all = TRUE)

## final SCr (for people who had bloods at least 30 days post-discharge)

postdc <-data.frame()


for (i in all_pt$ID) {
  admission <- all_pt %>% filter (ID == i) %>% select (hosp_days) %>% unlist (use.names = FALSE)
  post_admission_bloods <- bloods %>% filter (ID == i, result_day >= admission + 30)
  postdc <- rbind(postdc, post_admission_bloods)
}

postdc$d90 <- (as.numeric(postdc$result_day-90)^2)^0.5

# df of patients with results at least 30 days post-discharge (721 patients)                                  
D90_SCr <- postdc %>% mutate (ClinicalEventResult = as.numeric(ClinicalEventResult)) %>% na.omit() %>% group_by (ID) %>% arrange (ID, d90) %>% 
  slice_head() %>% select (ID, D90_SCr = ClinicalEventResult)   

D90_SCrID <- D90_SCr$ID
                                            
final_creats <- merge(creats, D90_SCr, by = "ID", all = TRUE)                     

# added final scr 
scr <- merge (all_pt %>% select (ID, stage, cr_base), final_creats, by = "ID")

scr1 <- scr %>% pivot_longer (cols = cr_base:D90_SCr) %>% mutate (name = factor(name, levels = c("cr_base", "first_SCr", "peak_SCr", "discharge_SCr", "D90_SCr")), valuemgdl = value/88.4)

#hospital survivors only (n = 1361)

hosp_survivors <- all_pt %>% filter (diednew == 0 | died_daysnew > hosp_days) %>% select(ID) %>% unlist(use.names = FALSE)

scr_survivors <- scr %>% filter (ID %in% hosp_survivors)

writexl::write_xlsx(scr_survivors, "scr_survivors.xlsx")

# plot S2A
scr_survivors_long <- scr1 %>% filter (ID %in% hosp_survivors) %>% mutate (stage = factor(stage, levels = c(seq(0, 4))))

boxplot <- ggboxplot(scr_survivors_long %>% filter(name != "D90_SCr") %>% na.omit(), x = "stage", y = "valuemgdl", color = "name", outlier.shape = NA) + 
  scale_color_discrete(name = NULL, breaks= c("cr_base", "first_SCr", "peak_SCr", "discharge_SCr"), 
                       labels = c("Baseline", "First", "Peak", "Last")) +
  scale_x_discrete(breaks  = c(seq(0, 3)), labels = c("No AKI", seq(1:3))) +
  ylab( "Serum Creatinine [mg/dl]")

png("boxplot.png", units="cm", width=17.1, height=14, res=400)
boxplot
dev.off()

#S2A log scale version 

boxplot_log <- ggboxplot(scr_survivors_long, x = "stage", y = "valuemgdl", color = "name") + 
  scale_color_discrete(name = NULL, breaks= c("cr_base", "first_SCr", "peak_SCr", "discharge_SCr"), 
                       labels = c("Baseline", "First", "Peak", "Last")) +
  scale_x_discrete(breaks  = c(seq(0, 3)), labels = c("No AKI", seq(1:3))) +
  scale_y_log10(name = "Serum Creatinine [mg/dl]")

png("boxplot_log.png", units="cm", width=17.1, height=14, res=400)
boxplot_log
dev.off()

## Plot S2B (those with scr results post discharge)


boxplot_2 <- ggboxplot(scr_survivors_long %>% filter(ID %in% D90_SCrID), x = "stage", y = "valuemgdl", color = "name", outlier.shape = NA) + 
 ylab( "Serum Creatinine [mg/dl]") +
  scale_color_discrete(name = NULL, breaks= c("cr_base", "first_SCr", "peak_SCr", "discharge_SCr", "D90_SCr"), 
                       labels = c("Baseline", "First", "Peak", "Last", "D90")) +
  scale_x_discrete(breaks  = c(seq(0, 3)), labels = c("No AKI", seq(1:3))) 

png("boxplot_2.png", units="cm", width=17.1, height=14, res=400)
boxplot_2
dev.off()


## plot S2C

scr_survivors <- read_xlsx("scr_survivors.xlsx")


creats_graph1 <- ggplot (data = scr_survivors, aes (x = cr_base/88.4, y = discharge_SCr/88.4, colour = stage)) + 
  geom_point(alpha = 0.8) + 
  scale_colour_manual(labels = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), values = c("grey", "#7846B4", "#FF9933", "#FF3399")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Baseline SCr [mg/dl]", y = "Last SCr before discharge [mg/dl]") + 
  theme_minimal() + 
  theme(legend.title = element_blank(), 
        legend.position = "top") +
  coord_cartesian(xlim = c(0.3, 10), ylim = c(0.3, 10)) +
  geom_abline(slope = 1, intercept = 0, colour = "darkgrey", linetype = "dotted")

creats_graph1

png("baseline_vs_discharge_SCr.png", units="cm", width=17.1, height=14, res=400)
creats_graph1
dev.off()

#Plot S2D by AKI yes/no
scr_survivors$aki_bin <- as.factor(ifelse(scr_survivors$stage == 0, 0, 1))

creats_graph2 <- ggplot (data = scr_survivors %>% filter(ID %in% D90_SCrID), aes (x = cr_base/88.4, y = D90_SCr/88.4, colour = aki_bin)) + 
  geom_point(alpha = 0.8) + 
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Baseline SCr [mg/dl]", y = "Day 90 SCr [mg/dl]") + 
  theme_minimal() + 
  theme(legend.title = element_blank(), 
        legend.position = "top") +
  scale_colour_manual(labels = c("No AKI", "Early AKI"), values = c("grey", "#7846B4")) +
  #coord_cartesian(xlim = c(0.03, 3.3), ylim = c(0.03, 3.3)) +
  geom_abline(slope = 1, intercept = 0, colour = "darkgrey", linetype = "dotted")

creats_graph2
  
png("baseline_vs_D90_SCr.png", units="cm", width=17.1, height=14, res=400)
creats_graph2
dev.off()

#Plot S2D by AKI stage 

creats_graph3 <- ggplot (data = scr_survivors %>% filter(ID %in% D90_SCrID), aes (x = cr_base/88.4, y = D90_SCr/88.4, colour = stage)) + 
  geom_point(alpha = 0.8) + 
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Baseline SCr [mg/dl]", y = "Day 90 SCr [mg/dl]") + 
  theme_minimal() + 
  theme(legend.title = element_blank(), 
        legend.position = "top") +
  scale_colour_manual(labels = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), values = c("grey", "#7846B4", "#FF9933", "#FF3399")) +
  #coord_cartesian(xlim = c(0.03, 3.3), ylim = c(0.03, 3.3)) +
  geom_abline(slope = 1, intercept = 0, colour = "darkgrey", linetype = "dotted")


png("baseline_vs_D90_SCr_stages.png", units="cm", width=17.1, height=14, res=400)
creats_graph3
dev.off()

## Statistical tests 

# data preparation 
library(rstatix)


# Plot S2C

# Baseline vs discharge 
stats <- scr_survivors %>% mutate(across(-c(ID, stage, aki_bin), function (x){x/88.4})) %>% group_by (stage) %>% summarise (base_median = median(cr_base, na.rm = TRUE), base_q1 = quantile(cr_base, 1/4, na.rm = TRUE), base_q3 = quantile (cr_base, 3/4, na.rm = TRUE), 
  last_median = median(discharge_SCr, na.rm = TRUE), last_q1 = quantile(discharge_SCr, 1/4, na.rm = TRUE), last_q3 = quantile (discharge_SCr, 3/4, na.rm = TRUE)) %>% 
  mutate(across(-stage, function(x) {format(round(x, digits = 2), nsmall = 2)}))

# no AKI

st0 <- scr_survivors %>% filter (stage == 0)
st0res <- wilcox.test(st0$cr_base, st0$discharge_SCr, paired = TRUE)
st0res$p.value

# AKI stage 1
st1 <- scr_survivors %>% filter (stage == 1)
st1res <- wilcox.test(st1$cr_base, st1$discharge_SCr, paired = TRUE)
st1res$p.value

# AKI stage 2
st2 <- scr_survivors %>% filter (stage == 2)
st2res <- wilcox.test(st2$cr_base, st2$discharge_SCr, paired = TRUE)
st2res$p.value

#AKI stage 3
st3 <- scr_survivors %>% filter (stage == 3)
st3res <- wilcox.test(st3$cr_base, st3$discharge_SCr, paired = TRUE)
st3res$p.value

# Plot S2D

stats2 <- scr_survivors %>% filter (ID %in% D90_SCrID) %>% 
  mutate(across(-c(ID, stage, aki_bin), function (x){x/88.4})) %>% group_by (aki_bin) %>% 
  summarise (base_median = median(cr_base, na.rm = TRUE), base_q1 = quantile(cr_base, 1/4, na.rm = TRUE), base_q3 = quantile (cr_base, 3/4, na.rm = TRUE), 
  final_median = median(D90_SCr, na.rm = TRUE), final_q1 = quantile(D90_SCr, 1/4, na.rm = TRUE), final_q3 = quantile (D90_SCr, 3/4, na.rm = TRUE)) %>% 
  mutate(across(-aki_bin, function(x) {format(round(x, digits = 2), nsmall = 2)}))

aki0 <- scr_survivors %>% filter (ID %in% D90_SCrID, stage == 0)
aki0res <- wilcox.test(aki0$cr_base, aki0$D90_SCr, paired = TRUE)
aki0res$p.value

aki1 <- scr_survivors %>% filter (ID %in% D90_SCrID, stage != 0)
aki1res <- wilcox.test(aki1$cr_base, aki1$D90_SCr, paired = TRUE)
aki1res$p.value


# stages
stats3 <- scr_survivors %>% filter (ID %in% D90_SCrID) %>%  mutate(across(-c(ID, stage, aki_bin), function (x){x/88.4})) %>% 
  group_by (stage) %>% 
  summarise (base_median = median(cr_base, na.rm = TRUE), base_q1 = quantile(cr_base, 1/4, na.rm = TRUE), base_q3 = quantile (cr_base, 3/4, na.rm = TRUE), 
             final_median = median(D90_SCr, na.rm = TRUE), final_q1 = quantile(D90_SCr, 1/4, na.rm = TRUE), final_q3 = quantile (D90_SCr, 3/4, na.rm = TRUE)) %>% 
  mutate(across(-stage, function(x) {format(round(x, digits = 2), nsmall = 2)}))

# AKI stage 0
st0 <- scr_survivors %>% filter (ID %in% D90_SCrID, stage == 0)
st0res <- wilcox.test(st0$cr_base, st0$D90_SCr, paired = TRUE)
st0res$p.value

# AKI stage 1
st1 <- scr_survivors %>% filter (ID %in% D90_SCrID,stage == 1)
st1res <- wilcox.test(st1$cr_base, st1$D90_SCr, paired = TRUE)
st1res$p.value

# AKI stage 2
st2 <- scr_survivors %>% filter (ID %in% D90_SCrID,stage == 2)
st2res <- wilcox.test(st2$cr_base, st2$D90_SCr, paired = TRUE)
st2res$p.value

#AKI stage 3
st3 <- scr_survivors %>% filter (ID %in% D90_SCrID,stage == 3)
st3res <- wilcox.test(st3$cr_base, st3$D90_SCr, paired = TRUE)
st3res$p.value




# Unused plot: discharge:admission creatinine vs length of stay 

scr_survivors_days <- merge (scr_survivors, all_pt %>% select (ID, hosp_days), by = "ID") %>% mutate (ratio = discharge_SCr/first_SCr)

creats_los <- ggplot (data = scr_survivors_days, aes (x = hosp_days, y = ratio, colour = stage)) + 
  geom_point(alpha = 0.8) + 
  scale_colour_manual(labels = c("No AKI", "Stage 1", "Stage 2", "Stage 3"), values = c("grey", "#7846B4", "#FF9933", "#FF3399")) +
  labs(x = "Length of stay [days]", y = "Discharge:admission SCr ratio") + 
  theme_minimal() + 
  theme(legend.title = element_blank(), 
        legend.position = "top")+ 
  coord_cartesian(ylim = c(-1, 10)) +
  geom_hline(yintercept = 1, colour = "darkgrey", linetype = "dotted") + 
  scale_y_continuous(breaks = c(seq(-1, 9, by = 2)))

creats_los
    
png("creats_los.png", units="cm", width=17.1, height=14, res=400)
creats_los
dev.off()

