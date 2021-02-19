all_pt <- read_xlsx("all_pt.xlsx")


all_pt$EthnicGroup2 <- factor (all_pt$EthnicGroup2, levels = c("White", "Asian", "Black", "Other", "Unknown"))
all_pt$quintile_study <- factor (all_pt$quintile_study, levels = c( "1", "2", "3", "4", "5"), labels = c("1", "2", "3", "4", "5"))
all_pt$stage <- factor(all_pt$stage, levels = c(0, 1, 2, 3))


## ICU admission 
glm.icu <- glm(icu ~ stage +
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

summary(glm.icu)

odds_ratios_icu <- OR_extractor(glm.icu)

x_axis <- c("stage1" = "AKI Stage 1", "stage2" = "AKI Stage 2", "stage3" = "AKI Stage 3", 
            "smoking" = "Smoking", "SexMale" = "Male sex", 
            "quintile_study5" = "IMD 5", "quintile_study4" = "IMD 4", 
            "quintile_study3" = "IMD 3", "quintile_study2" = "IMD 2", 
            "obes" = "Obesity", "HTN1" = "HTN", "crp_bin" = "CRP > 145", 
            "EthnicGroup2Other" = "Ethnicity: Other", "EthnicGroup2Black" = "Ethnicity: Black", 
            "EthnicGroup2Asian" = "Ethnicity: Asian", "EthnicGroup2Unknown" = "Ethnicity: Unknown","Diabetes1" = "Diabetes", "low_gfr_baseline" = "CKD",
            "chf" = "CHF", "cevd" = "CEVD", "age_div" = "Age")

odds_ratios_icu$predictor <- factor(odds_ratios_icu$predictor, levels <- 
                                        rev(c("stage1", "stage2", "stage3", "SexMale", "age_div",
                                              "EthnicGroup2Asian", "EthnicGroup2Black",  "EthnicGroup2Other", "EthnicGroup2Unknown", 
                                              "quintile_study5" ,"quintile_study4" , "quintile_study3", "quintile_study2" , 
                                              "smoking" , "obes", "low_gfr_baseline", "HTN", "Diabetes","cevd", "chf", "crp_bin")))

icu <- make_log_forest (odds_ratios_icu, breaks = 10, title = " ", xlabs = x_axis, a = 40, b = 350, c = 2000, metric = "OR")

icu_fp <- icu + scale_y_log10(breaks = c(1, 5, 10, 20, 40), labels = c(1, 5, 10, 20, 40))


png("icu_fp.png", units="cm", width=17.1, height=17.1, res=400)
icu_fp 
dev.off()



# n obs and nevents

dd <- datadist (all_pt)
options (datadist = "dd")


glm.icu2 <- lrm(icu ~ stage +
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
                ,data = all_pt)

glm.icu2

# icumv


glm.icumv <- glm(icumv ~ stage +
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

summary(glm.icumv)



odds_ratios_icumv <- OR_extractor(glm.icumv)

odds_ratios_icumv$predictor <- factor(odds_ratios_icumv$predictor, levels <- 
                                      rev(c("stage1", "stage2", "stage3", "SexMale", "age_div",
                                            "EthnicGroup2Asian", "EthnicGroup2Black",  "EthnicGroup2Other", "EthnicGroup2Unknown", 
                                            "quintile_study5" ,"quintile_study4" , "quintile_study3", "quintile_study2" , 
                                            "smoking" , "obes", "low_gfr_baseline", "HTN", "Diabetes","cevd", "chf", "crp_bin")))

icumv <- make_log_forest (odds_ratios_icumv, breaks = 100, title = " ", xlabs = x_axis, a = 500, b = 7000, c = 70000, metric = "OR")

icumv_fp <- icumv + scale_y_log10(breaks = c(1, 10, 50, 200), labels = c(1, 10, 50, 200))

png("icumv_fp.png", units="cm", width=17.1, height=17.1, res=400)
icumv_fp 
dev.off()

# n obs and events 

dd <- datadist (all_pt)
options (datadist = "dd")
glm.icumv2 <- lrm(icumv ~ stage +
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

glm.icumv2



## rrt 


glm.rrt <- glm(rrt ~ stage + 
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
                   hlh_bin, na.action = na.omit, 
                 data = all_pt, family = binomial)

summary(glm.rrt)


odds_ratios_rrt <- OR_extractor(glm.rrt)

odds_ratios_rrt$predictor <- factor(odds_ratios_rrt$predictor, levels <- 
                                        rev(c("stage1", "stage2", "stage3", "SexMale", "age_div",
                                              "EthnicGroup2Asian", "EthnicGroup2Black",  "EthnicGroup2Other", "EthnicGroup2Unknown", 
                                              "quintile_study5" ,"quintile_study4" , "quintile_study3", "quintile_study2" , 
                                              "smoking" , "obes", "low_gfr_baseline", "HTN", "Diabetes","cevd", "chf", "hlh_bin")))

rrt <- make_log_forest (odds_ratios_rrt, breaks = 100, title = " ", xlabs = x_axis, a = 400, b = 5000, c = 60000, metric = "OR")


rrt_fp <- rrt + scale_y_log10(breaks = c(1, 10, 50, 200), labels = c(1, 10, 50, 200))


png("rrt_fp.png", units="cm", width=17.1, height=17.1, res=400)
rrt_fp 
dev.off()

# n obs and events 
glm.rrt <- lrm(rrt ~ stage + 
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
                 hlh_bin
               ,data = all_pt)

glm.rrt

## MAKE 90 

glm.make90 <- glm(make90 ~ stage + 
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
                 crp_bin, na.action = na.omit, 
                data = all_pt, family = binomial)


odds_ratios_make90 <- OR_extractor(glm.make90)

odds_ratios_make90$predictor <- factor(odds_ratios_make90$predictor, levels <- 
                                      rev(c("stage1", "stage2", "stage3", "SexMale", "age_div",
                                            "EthnicGroup2Asian", "EthnicGroup2Black",  "EthnicGroup2Other", "EthnicGroup2Unknown", 
                                            "quintile_study5" ,"quintile_study4" , "quintile_study3", "quintile_study2" , 
                                            "smoking" , "obes", "low_gfr_baseline", "HTN", "Diabetes","cevd", "chf", "crp_bin")))

make90 <- make_log_forest (odds_ratios_make90, breaks = 5, title = " ", xlabs = x_axis, a = 40, b = 250, c = 1000, metric = "OR")


make90_fp <- make90 + scale_y_log10(breaks = c(1, 5, 10, 20, 40), labels = c(1, 5, 10, 20, 40))


png("make90_fp.png", units="cm", width=17.1, height=17.1, res=400)
make90_fp 
dev.off()


# n obs and events
glm.make90.2 <- lrm(make90 ~ stage + 
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
                    crp_bin,data = all_pt)

glm.make90.2
