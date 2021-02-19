setwd("/Users/Zuzanna_Bien/Data/Files")

### EARLY LIST - FIRST 7 DAYS OF ADMISSION 
all_pt <- read_xlsx("all_pt.xlsx")
dead7 <- all_pt %>% filter (died_daysnew <= 7)

# exclude everyone who died before D7

trajectories_pt <- all_pt %>% filter (!died_daysnew <= 7 | is.na(died_daysnew) == TRUE )
trajectories_bloods <- bloods %>% filter (ID %in% trajectories_pt$ID)


early_aki <- trajectories_pt %>% filter (aki_bin == 1)
early_AKI_ID <- early_aki$ID
no_aki <- trajectories_pt %>% filter (stage == 0) 
no_AKI_ID <- no_aki$ID


early_t <- trajectories_bloods %>% 
  filter (ClinicalEvent == "Creatinine Serum" & result_day >= -1 & result_day <= 7 ) %>% 
  arrange(ID, result_day) %>% mutate (ClinicalEventResult = as.numeric (ClinicalEventResult)) %>% na.omit ()

# DETECTION OF EARLY -PERSISTENT AKI BASED ON THE ABSOLUTE CRITERION 

early_t_list <- split(early_t, f = early_t$ID) 

t <- data.frame()
v <- vector()

for (a in 1:length(early_t_list)) {
  
  b <- early_t_list[[a]]
  
  for (i in nrow(b):2) {
    for (x in 1:(i-1)) {
      v <- c(b$ID[1], 
             b$ClinicalEventResult[i],
             b$ClinicalEventResult[x],
             b$result_day[i],
             b$result_day[x],
             b$result_time[i],
             b$result_time[x])
      t <- rbind(t, v, stringsAsFactors = FALSE)
    }
  }
}


colnames (t) <- c("ID","SCr_max", "SCr_min", "D_max", "D_min", "T_max", "T_min")
cols.num <- c("SCr_max", "SCr_min", "D_max", "D_min")
t[cols.num] <- sapply(t[cols.num],as.numeric)


# create a vector of IDs of patients with persistent AKI based on the absolute criterion 

absolute_Scr <- t %>% na.omit() %>% mutate (SCr_diff = SCr_max - SCr_min, D_diff = D_max - D_min, T_diff = T_max >= T_min)

# patients who met the absolute criterion in the last 48h of the first 7day period
ab_pers_7 <- absolute_Scr %>% filter (SCr_diff > 26 & D_diff <=2 & D_min %in% c(6, 7)) %>% arrange(ID) %>% select(ID, SCr_min, SCr_max, SCr_diff, D_min, D_max, D_diff, T_min, T_max, T_diff)

ab_pers_ID <- unique(ab_pers_7$ID)

length(ab_pers_ID)


# EARLY PERSISTENT AKI - DETECTION BASED ON THE RELATIVE CRITERION

peak_scr  <- early_t %>% na.omit () %>% filter (result_day %in% c(6:7)) %>%  group_by (ID) %>% 
  summarise (peak_scr_d7 = max(ClinicalEventResult))

pers_rel_ID <- merge(trajectories_pt %>% select (ID, baseline_SCr), peak_scr, by = "ID") %>%
  mutate (scr_change = peak_scr_d7/baseline_SCr) %>%
  mutate (pers_rel_aki = ifelse (scr_change > 1.5, 1, 0)) %>% 
  filter (pers_rel_aki == 1) %>% select (ID) %>% unlist(use.names = FALSE)

## Combine all patients with persistent AKI with those who received RRT (they will be classified as persistent AKI)

ab_rel_pers_ID <- union(ab_pers_ID, pers_rel_ID)
rrt_ID <- early_aki %>% filter (rrt == 1) %>% select (ID)  %>% unlist(use.names = FALSE)
pers_ID <- as.numeric(union(ab_rel_pers_ID, rrt_ID))

length(pers_ID)

# RELAPSED AKI: BASED ON ABSOLUTE CRITERION 

recov <- early_aki %>% filter (!ID %in% pers_ID)
recovered <- trajectories_bloods %>% filter (ClinicalEvent == "Creatinine Serum" & ID %in% recov$ID & result_day >= 7) %>%
  arrange (ID, result_day) %>% mutate (ClinicalEventResult = as.numeric(ClinicalEventResult)) %>% na.omit()


recovered_list <- split(recovered, f = recovered$ID) 

n <- data.frame()
v <- vector()

for (a in 1:length(recovered_list)) {
  
  b <- recovered_list[[a]]
  
  for (i in nrow(b):2) {
    for (x in 1:(i-1)) {
      v <- c(b$ID[1], 
             b$ClinicalEventResult[i],
             b$ClinicalEventResult[x],
             b$result_day[i],
             b$result_day[x],
             b$result_time[i],
             b$result_time[x])
      n <- rbind(n, v, stringsAsFactors = FALSE)
    }
  }
}


colnames (n) <- c("ID","SCr_max", "SCr_min", "D_max", "D_min", "T_max", "T_min")
cols.num <- c("SCr_max", "SCr_min", "D_max", "D_min")
n[cols.num] <- sapply(n[cols.num],as.numeric)

write_xlsx(n, "n.xlsx")

relapsed_ab_Scr <- n %>% mutate (SCr_diff = SCr_max - SCr_min, D_diff = D_max - D_min, T_diff = T_max >= T_min)

relapsed_aki_ab_SCr <- relapsed_ab_Scr %>% filter (SCr_diff > 26 & D_diff <=2 & D_diff >=0) %>% filter (!(D_diff == 2 & T_diff == TRUE))%>% arrange(ID) %>% select(ID, SCr_min, SCr_max, SCr_diff, D_min, D_max, D_diff, T_min, T_max, T_diff)

relapsed_ab_ID <- unique(relapsed_aki_ab_SCr$ID)


### RELAPSED AKI: DETECTION BASED ON RELATIVE CRITERION 

peak_scr2  <- recovered %>% filter (result_day > 7) %>% group_by (ID) %>% summarise (peak_scr = max(ClinicalEventResult))

relapsed_rel_ID <- merge(trajectories_pt %>% select (ID, baseline_SCr), peak_scr2, by = "ID") %>%
  mutate (scr_change = peak_scr/baseline_SCr) %>%
  mutate (pers_rel_aki = ifelse (scr_change > 1.5, 1, 0)) %>% 
  filter (pers_rel_aki == 1) %>% select (ID) %>% unlist(use.names = FALSE)


#all patients with relapsed ID - based on relative and absolute criterion 
relapsed_ID <- union(relapsed_rel_ID, relapsed_ab_ID)

length(relapsed_ID)

## recovered no relapse 

non_recovered_ID <- union(pers_ID, relapsed_ID)
recovered_ID <- setdiff(early_aki$ID, non_recovered_ID)
length(recovered_ID)

subgroups <- trajectories_pt %>% 
  mutate (fate = case_when (
    ID %in% relapsed_ID ~ "early - relapsed", 
    ID %in% pers_ID ~ "early - persistent", 
    ID %in%recovered_ID ~ "early - recovered", 
    ID %in% no_AKI_ID ~ "no AKI"))

write_xlsx(subgroups, "subgroups.xlsx")
  
# LATE AKI 

late <- trajectories_bloods %>% 
  filter (ID %in% no_AKI_ID & ClinicalEvent == "Creatinine Serum" & result_day > 7) %>%
  arrange(ID, result_day) %>% mutate (ClinicalEventResult = as.numeric (ClinicalEventResult)) %>% na.omit ()

# LATE AKKI: RELATIVE CRITERION 
peak_scr_late  <- late %>% group_by (ID) %>% summarise (peak_scr_late = max(ClinicalEventResult))

late_rel_ID <- merge(trajectories_pt %>% select (ID, baseline_SCr), peak_scr_late, by = "ID") %>%
  mutate (scr_change = peak_scr_late/baseline_SCr) %>%
  mutate (late_rel_aki = ifelse (scr_change > 1.5, 1, 0)) %>% 
  filter (late_rel_aki == 1) %>% select (ID) %>% unlist(use.names = FALSE)

length(late_rel_ID)

# LATE AKI: ABSOLUTE CRITERION 

late_list <- split(late, f = late$ID) 

o <- data.frame()
v <- vector()

for (a in 1:length(late_list)) {
  
  b <- late_list[[a]]
  
  for (i in nrow(b):2) {
    for (x in 1:(i-1)) {
      v <- c(b$ID[1], 
             b$ClinicalEventResult[i],
             b$ClinicalEventResult[x],
             b$result_day[i],
             b$result_day[x],
             b$result_time[i],
             b$result_time[x])
      o <- rbind(o, v, stringsAsFactors = FALSE)
    }
  }
}

colnames (o) <- c("ID","SCr_max", "SCr_min", "D_max", "D_min", "T_max", "T_min")
cols.num <- c("SCr_max", "SCr_min", "D_max", "D_min")
o[cols.num] <- sapply(o[cols.num],as.numeric)


write_xlsx(o, "o.xlsx")

late_absolute_SCr <- o %>% na.omit() %>% mutate (D_diff = D_max - D_min, T_diff = T_max >= T_min, SCr_diff = SCr_max - SCr_min)

late_ab_SCr <- late_absolute_SCr %>% filter (SCr_diff > 26 & D_diff <=2) %>% filter (!(D_diff == 2 & T_diff == TRUE)) %>% arrange(ID) %>% select(ID, SCr_min, SCr_max, SCr_diff, D_min, D_max, D_diff, T_min, T_max, T_diff)

write_xlsx(late_ab_SCr, "late_absolute_SCr.xlxs")

# find unique cases 

late_ab_SCr <- read_xlsx("late_absolute_SCr.xlxs")

ab_late_ID <- unique (late_ab_SCr$ID)

length(ab_late_ID)

## combine the two criteria  

late_rrt_ID <- trajectories_pt %>% filter (aki_bin == 0, rrt == 1) %>% select (ID) %>% unlist (use.names = FALSE)

late_ID <- union(late_rel_ID, ab_late_ID) 
all_late_ID <- union(late_ID, late_rrt_ID)

length(all_late_ID)

## combined with the earlier subgroup labels

subgroups <- read_xlsx("subgroups.xlsx")
  
fates <- subgroups %>% mutate (fate = ifelse (ID %in% all_late_ID, "late", fate))

write_xlsx(fates, "fates.xlsx")

#### SUBGROUP ANALYSIS 

## Table 1

final_fate <- read_xlsx("fates.xlsx")

final_fate$fate <- factor (final_fate$fate, levels = c("no AKI", "early - recovered", "early - relapsed", "early - persistent", "late"))
final_fate$stage <- factor(final_fate$stage, levels = c(0, 1, 2, 3), labels = c(0, 1, 2, 3))
final_fate$quintile_study <- factor (final_fate$quintile_study, levels = c(seq(1,5)))
final_fate$EthnicGroup2 <- factor (final_fate$EthnicGroup2, levels = c("White", "Asian", "Black", "Other", "Unknown"))



final_fates <- merge(final_fate, all_pt %>% 
                       select (ID, baseline_SCr, last_SCr, baseline_SCr_mgdl, last_SCr_mgdl, base_gfr, last_gfr, low_gfr_base, low_gfr_final, unrecovery, make90), by = "ID") %>% 
  select (-"baseline_SCr_mgdl.x", -"last_SCr_mgdl.x", -"base_gfr.x", -"last_gfr.x", - "make90.x" , -"low_gfr_final.x", -"low_gfr_base.x" , -"last_SCr.x", -"baseline_SCr.x", -"unrecovery.x" ) %>%
  rename ("baseline_SCr" = "baseline_SCr.y", "last_SCr" = "last_SCr.y","baseline_SCr_mgdl" = "baseline_SCr_mgdl.y", "last_SCr_mgdl"="last_SCr_mgdl.y",
          "base_gfr"="base_gfr.y", "last_gfr"="last_gfr.y", "low_gfr_base" ="low_gfr_base.y", "low_gfr_final"="low_gfr_final.y", "unrecovery"="unrecovery.y" , "make90" ="make90.y"  )
  
sum(is.na(final_fates$base_gfr) == FALSE)

write_xlsx(final_fates, "final_fates.xlsx")

# TABLE ONE

final_fates <- read_xlsx("final_fates.xlsx")
final_fates$fate <- factor (final_fates$fate, levels = c("no AKI", "early - recovered", "early - relapsed", "early - persistent", "late"))

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

tab_sens <- CreateTableOne(vars = myVars, strata = "fate", data = final_fates, factorVars = catVars)

write.csv(print(tab_sens, nonnormal = c("AdmissionAge", "baseline_SCr_mgdl", "last_SCr_mgdl", "base_gfr", "last_gfr", "died_daysnew", "hosp_days", "iculos", "crp"), formatOptions = list(big.mark = ",")), "tableone_traj.csv")


write_xlsx(as.data.frame(sapply (final_fates %>% select (myVars), function (x) length(x) - sum(is.na(x)))) %>% rownames_to_column() %>% rename (var = 1, n = 2), "n_tab_traj.xlsx")

# TABLE ONE: MAKE 90
make90 <- final_fates %>% filter(make90 == 1)

make90$survunrec <- case_when (make90$died90new == 0 & make90$unrecovery == 1 ~ 1)
make90$survunrec <- replace_na(make90$survunrec, 0)
make90$survrrt <- case_when (make90$died90new == 0 & make90$rrt == 1 ~ 1)
make90$survrrt <- replace_na(make90$survrrt, 0)

makeVars <- c("survrrt", "survunrec")
catVars<- c("survrrt", "survunrec")

tab_make90traj <- CreateTableOne(vars = makeVars, factorVars = catVars, strata = "fate", data = make90, addOverall = TRUE)
write.csv(print(tab_make90traj, formatOptions = list(big.mark = ",")), "tab_make90traj.csv")


##figure S6

final_fate$fate <- factor (final_fate$fate, levels = rev(c("no AKI", "early - recovered", "early - relapsed", "early - persistent", "late")))

ethnicity_trajectory_barchart <- ggplot(data = final_fate, aes(x = EthnicGroup2, fill = fate)) + 
  geom_bar(position = "fill", alpha = 0.7) +
  theme_minimal() + 
  coord_flip() + 
  theme (legend.title = element_blank(), 
         legend.text = element_text(size = 6),
         axis.text.x = element_text(size = 6),
         axis.text.y = element_text(size = 6),
         legend.position = "top",
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         plot.margin = margin(10, 40, 10, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(breaks = c("no AKI", "early - recovered", "early - relapsed", "early - persistent", "late"), 
                    labels = c("No AKI", "Early: Recovered", "Early: Relapsed", "Early: Persistent", "Late"),
                    values = c( "lightgrey", "#F8766D", "#00C19F", "#DB72FB","#619CFF")) + xlab("Ethnicity") 

ethnicity_trajectory_barchart    

png("ethnicity_trajectory_barchart_horizontal.png", units="cm", width=11.4, height=10.0, res=400)
ethnicity_trajectory_barchart
dev.off()


# Logistic regression for survival to 90 days 

final_fates_baseline_scr <- merge(final_fates, all_pt%>% select (ID, low_gfr_baseline), by = "ID")

final_fates_baseline_scr$EthnicGroup2 <- factor (final_fates_baseline_scr$EthnicGroup2, levels = c("White", "Asian", "Black", "Other", "Unknown"))

glm.trajectories <- glm(died90new ~  
                  fate +
                  age_div +  
                  Sex +
                  EthnicGroup2 + 
                  quintile_study + 
                  smoking + 
                  obes +
                  low_gfr_baseline +
                  chf + 
                  Diabetes + 
                  HTN + 
                  cevd + 
                  crp_bin, data = final_fates_baseline_scr, na.action = na.omit, family = binomial)


summary(glm.trajectories)

OR_trajectories <- OR_extractor(glm.trajectories)


xaxis <- c("age_div" = "Age", "SexMale" = "Male sex", 
           "fateearly - recovered" = "Early AKI: Recovered", "fateearly - relapsed" = "Early AKI: Relapsed", 
           "fateearly - persistent" = "Early AKI: Persistent", "fatelate" = "Late AKI",
           "EthnicGroup2Black" = "Ethnicity: Black", "EthnicGroup2Unknown" = "Ethnicity: Unknown",
           "EthnicGroup2Other" = "Ethnicity: Other", "EthnicGroup2Asian" = "Ethnicity: Asian", 
           "quintile_study2" = "IMD 2", "quintile_study3" = "IMD 3", "quintile_study4" = "IMD 4",
           "quintile_study5" = "IMD 5", "smoking" = "Smoking", "obes" = "Obesity", "low_gfr_baseline" = "CKD", 
           "chf" = "CHF",  "Diabetes" = "Diabetes", "HTN" = "HTN", "cevd" = "CEVD", "crp_bin" = "CRP > 145")

OR_trajectories$predictor <- factor(OR_trajectories$predictor, levels <- 
                                        rev(c("fateearly - recovered", "fateearly - persistent", "fateearly - relapsed",  "fatelate",
                                              "SexMale", "age_div",
                                              "EthnicGroup2Asian", "EthnicGroup2Black",  "EthnicGroup2Other", "EthnicGroup2Unknown", 
                                              "quintile_study5" ,"quintile_study4" , "quintile_study3", "quintile_study2" , 
                                              "smoking" , "obes", "low_gfr_baseline", "HTN", "Diabetes","cevd", "chf", "crp_bin")))




trajectories_fp_glm <- make_log_forest(OR_trajectories, breaks =5, xlabs = xaxis, title = "", a = 10, b = 50, c =150, metric = "OR")


png("glm_multi_90_trajectories.png", units="cm", width=17.1, height=17.1, res=400)
trajectories_fp_glm
dev.off()


# obs and event s

glm.trajectories <- lrm(died90new ~  
                          fate +
                          age_div +  
                          Sex +
                          EthnicGroup2 + 
                          quintile_study + 
                          smoking + 
                          obes +
                          low_gfr_baseline +
                          chf + 
                          Diabetes + 
                          HTN + 
                          cevd + 
                          crp_bin, data = final_fates_baseline_scr)
glm.trajectories

## BOXPLOTS


final_fates <- read_xlsx("final_fates.xlsx")
final_fates$fate <- factor(final_fates$fate, levels = c("no AKI", "early - recovered", "early - relapsed", "early - persistent", "late"))

scr_survivors <- read_xlsx("scr_survivors.xlsx")

scr_survivors_traj <- merge (scr_survivors, final_fates %>% select (ID, fate),  by = "ID")

scr_survivors_traj_long <- scr_survivors_traj %>% select (-stage) %>% pivot_longer (cols = cr_base:D90_SCr) %>% mutate (name = factor(name, levels = c("cr_base", "first_SCr", "peak_SCr", "discharge_SCr", "D90_SCr")), valuemgdl = value/88.4)

# Figure S2E
boxplot <- ggboxplot(scr_survivors_traj_long %>% filter (name != "D90_SCr"), x = "fate", y = "valuemgdl", color = "name", outlier.shape = NA) + 
  scale_color_discrete(name = NULL, breaks= c("cr_base", "first_SCr", "peak_SCr", "discharge_SCr"), 
                       labels = c("Baseline", "First", "Peak", "Last")) +
  scale_x_discrete(breaks  = c("no AKI", "early - recovered", "early - relapsed", "early - persistent", "late"), labels = c("No AKI", "Early: Recovered", "Early: Relapsed", "Early: Persistent", "Late")) +
  ylab( "Serum Creatinine [mg/dl]") +
  theme(axis.text.x = element_text(angle = 15, hjust = 0.5, vjust = 0.5)) + xlab (NULL) 

png("boxplot_traj1.png", units="cm", width=17.1, height=14, res=400)
boxplot
dev.off()

#Figure S2F
boxplot_traj2 <- ggboxplot(scr_survivors_traj_long %>% filter (ID %in%D90_SCrID), x = "fate", y = "valuemgdl", color = "name", outlier.shape = NA) + 
 # scale_y_log10(name = "Serum Creatinine [mg/dl]") +
  scale_color_discrete(name = NULL, breaks= c("cr_base", "first_SCr", "peak_SCr", "discharge_SCr", "D90_SCr"), 
                       labels = c("Baseline", "First", "Peak", "Last", "D90")) +
  scale_x_discrete(breaks  = c("no AKI", "early - recovered", "early - relapsed", "early - persistent", "late"), labels = c("No AKI", "Early: Recovered", "Early: Relapsed", "Early: Persistent", "Late")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 0.5, vjust = 0.5)) + xlab (NULL) + ylab("Serum Creatinine [mg/dl]")

boxplot_traj2 

png("boxplot_traj2.png", units="cm", width=17.1, height=14, res=400)
boxplot_traj2
dev.off()

#S2F log scale
boxplot_traj2log <- ggboxplot(creatinine1, x = "fate", y = "valuemgdl", color = "name") + 
  scale_y_log10(name = "Serum Creatinine [mg/dl]") +
  scale_color_discrete(name = NULL, breaks= c("cr_base", "first_SCr", "peak_SCr", "last_SCr", "last_post_discharge"), 
                       labels = c("Baseline", "First", "Peak", "Last", "D90")) +
  scale_x_discrete(breaks  = c("no AKI", "early - recovered", "early - relapsed", "early - persistent", "late"), labels = c("No AKI", "Early: Recovered", "Early: Relapsed", "Early: Persistent", "Late")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 0.5, vjust = 0.5)) + xlab (NULL) 

boxplot_traj2log 

png("boxplot_traj2log.png", units="cm", width=17.1, height=14, res=400)
boxplot_traj2log
dev.off()
