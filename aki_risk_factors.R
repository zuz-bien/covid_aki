aki_fates <- read_xlsx("all_pt.xlsx")

# glm model 

aki_fates$EthnicGroup2 <- factor (aki_fates$EthnicGroup2, levels = c("White", "Asian", "Black", "Other", "Unknown"))
aki_fates$quintile_study <- factor (aki_fates$quintile_study, levels = c( "1", "2", "3", "4", "5"), labels = c("1", "2", "3", "4", "5"))

glm.fit1 <- glm(aki_bin ~  
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
                 data = aki_fates, family = binomial)



dd <- datadist (aki_fates)
options (datadist = "dd")


glm.fit2 <- lrm(aki_bin ~  
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
                  crp_bin, data = aki_fates)


#round values

odds_ratios_round <- OR_extractor(glm.fit1)

xaxis <- c("age_div" = "Age", "SexMale" = "Male sex", 
           "EthnicGroup2Black" = "Ethnicity: Black", "EthnicGroup2Unknown" = "Ethnicity: Unknown",
           "EthnicGroup2Other" = "Ethnicity: Other", "EthnicGroup2Asian" = "Ethnicity: Asian", 
           "quintile_study2" = "IMD 2", "quintile_study3" = "IMD 3", "quintile_study4" = "IMD 4",
           "quintile_study5" = "IMD 5", "smoking" = "Smoking", "obes" = "Obesity", "low_gfr_baseline" = "CKD", 
           "chf" = "CHF",  "Diabetes" = "Diabetes", "HTN" = "HTN", "cevd" = "CEVD", "crp_bin" = "CRP > 145")

odds_ratios_round$predictor <- factor(odds_ratios_round$predictor, levels <- 
                                           rev(c("SexMale", "age_div",
                                                 "EthnicGroup2Asian", "EthnicGroup2Black",  "EthnicGroup2Other", "EthnicGroup2Unknown", 
                                                 "quintile_study5" ,"quintile_study4" , "quintile_study3", "quintile_study2" , 
                                                 "smoking" , "obes", "low_gfr_baseline", "HTN", "Diabetes","cevd", "chf", "crp_bin")))


aki_rf <- make_log_forest(odds_ratios_round, breaks =1, xlabs = xaxis, title = "", a = 1, b = 5, c = 12, metric = "OR")

png("aki_rf.png", units="cm", width=17.1, height=17.1, res=400)
aki_rf 
dev.off()

