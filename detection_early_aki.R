# DETECTION OF EARLY AKI 

#Set working directory and read in the files 

setwd("/Users/Zuzanna_Bien/Data/Files")

blood <- read.csv("covidaki_bloods_01dec20.csv", stringsAsFactors = FALSE)
patients<- read.csv("covidaki_pt_04feb21.csv", stringsAsFactors = FALSE)

# only use creatinine values 
bloods <- blood %>% filter (ClinicalEvent == "Creatinine Serum")

# filter out patients with ESRF and those who did not have any creatinine values recorded on admission
ESRF_ID <- patients %>% filter (esrf == 1) %>% select (ID) %>% unlist (use.names = FALSE)
no_SCr <- setdiff(blood$ID, bloods$ID)
excluded_pt <- union(no_SCr, ESRF_ID)

#data frame after excluding the above groups 
#note, more patients will later be excluded due to missingness; 1873 is not the final value
patients_1873 <- patients %>% filter (!ID %in%excluded_pt)

### DETECTION OF EARLY AKI 

#filter bloods to only include results between day -1 and 7
early <- bloods %>% 
  filter (result_day >= -1 & result_day <= 7 ) %>% 
  filter (!ID %in% excluded_pt) %>% 
  arrange(ID, result_day) 

early$ClinicalEventResult <- as.numeric (early$ClinicalEventResult)

# DETECTION OF AKI BASED ON THE >26 in 48h CRITERION 

early_list <- split(early, f = early$ID) 

#initialise empty vector and data frame
m <- data.frame()
v <- vector()

# loop through the results to create a new dataframe which pulls each creatinine value with
# date and time stamp for each row of results for each individual patient
# and compares it to the values from the row prior (arranged chronologically)

for (a in 1:length(early_list)) {
  
  b <- early_list[[a]]

for (i in nrow(b):2) {
  for (x in 1:(i-1)) {
    v <- c(b$ID[1], 
           b$ClinicalEventResult[i],
           b$ClinicalEventResult[x],
           b$result_day[i],
           b$result_day[x],
           b$result_time[i],
           b$result_time[x])
  m <- rbind(m, v, stringsAsFactors = FALSE)
  }
}
}

colnames (m) <- c("ID","SCr_max", "SCr_min", "D_max", "D_min", "T_max", "T_min")
cols.num <- c("SCr_max", "SCr_min", "D_max", "D_min")
m[cols.num] <- sapply(m[cols.num],as.numeric)

# m = data frame with timestamped creatinine values 
write_xlsx(m, "m.xlsx")

# SCr_max here means creatinine value taken at the later time. SCr_diff is the difference
# in creatinine between later and earlier samples
# Similarly D_max is the day stamp of the later value and T_max is the timestamp of the later value
# Because of the way the times are formatted, you can compare them using >= (i.e. 13:40 >= 11:15 will equal TRUE)
absolute_Scr <- m %>% na.omit() %>% mutate (SCr_diff = SCr_max - SCr_min, D_diff = D_max - D_min, T_diff = T_max >= T_min)

# Now filter for rows where creatinine difference is > 26 and the day difference is either
# >2, or it is equal to 2 but the later value was taken at a later time in the day. I.e. creatinine difference
# of >26 between values from Day 1 at 11:00 and Day 3 at 11:30 will count as AKI, but not if the 
# Day 1 value was at 11:30 and D3 at 11:00 (as this is <48h) 
aki_absolute_SCr <- absolute_Scr %>% filter (SCr_diff > 26 & D_diff <=2) %>% filter (!(D_diff == 2 & T_diff == TRUE)) %>% arrange(as.numeric(ID)) %>% select(ID, SCr_min, SCr_max, SCr_diff, D_min, D_max, D_diff, T_min, T_max, T_diff)

write_xlsx(aki_absolute_SCr, "aki_absolute_SCr.xlsx")

aki_absolute_SCr <- read_xlsx("aki_absolute_SCr.xlsx")

### DETECTION OF AKI BASED ON THE RELATIVE CRITERION (SCr > 1.5 times baseline)

#pull out baseline cretinine values 
baseline <- patients_1873 %>% select(ID, cr_base, cr_imp) %>% transmute(ID = ID, baseline_SCr = coalesce(cr_base, cr_imp)) 

#pull out the highest creatinine value on admission 
peak <- early %>% group_by(ID) %>% na.omit () %>% summarise (peak_SCr = max(ClinicalEventResult))

#store IDs of patients with the highest creatinine >353
highscr <- peak %>% filter (peak_SCr > 353) %>% select (ID) %>% unlist (use.names = FALSE)

#stage according to peak value divided by baseline
base_peak <- merge(baseline, peak, by = "ID")

relative_SCr <- base_peak %>% 
  mutate (relative_SCr_change = peak_SCr / baseline_SCr)  %>%
  mutate (stage = cut(relative_SCr_change, breaks =  c(-Inf, 1.5, 2, 3, Inf),  include.lowest=TRUE, labels = c(0, 1, 2, 3)))

relative_SCr$stage <- as.character(relative_SCr$stage)
aki_relative_SCr <- relative_SCr %>% filter(stage != "0")

write_xlsx(aki_relative_SCr, "aki_relative_SCr.xlsx")

## merge the two 

aki_absolute_SCr <- read_xlsx("aki_absolute_SCr.xlsx")
aki_relative_SCr <- read_xlsx("aki_relative_SCr.xlsx")

#IDs of patients with early AKI according to the absolute criterion
ID_absolute <- unique(aki_absolute_SCr$ID)

#IDs of patients with early AKI according to the relative criterion
ID_relative <- unique(aki_relative_SCr$ID)

# setdiff - rows that appear in the first dataset, but not the second
# IDs of patients who only qualified based on the absolute criterion
absolute_only <- as.numeric(setdiff(ID_absolute, ID_relative))
# IDs of patients who only qualified based on the relative criterion
relative_only <- setdiff(ID_relative, ID_absolute)

#all cases of AKI
all_AKI_IDs <- union(ID_absolute, ID_relative)

# create a dataframe of all AKI patients 

aki_absolute_SCr$ID <- as.numeric(aki_absolute_SCr$ID)
aki_absolute_SCr$stage <- factor(1)

absolute <- data.frame (ID = c(unique(aki_absolute_SCr$ID)),stage_ab = factor(1))
relative <- aki_relative_SCr %>% select (ID, stage_rel = stage)

# assign stage to patients with early AKI
aki <- merge (absolute, relative, by = "ID", all = TRUE)  %>% mutate (stage = coalesce (stage_rel, stage_ab))

# merge the dataframe with all patient data with the AKI stages
aki_patients <- merge (patients_1873, aki %>% select (ID, stage), all = TRUE) %>% mutate (stage = ifelse(is.na(stage) == TRUE, "0", stage)) 

# reclassify all patients who ended up having RRT as stage 3
rrt_aki_patients <- aki_patients %>% filter (!stage == "0") %>% mutate (s_rrt = ifelse(rrt == 1, "3", NA)) %>%
  mutate(stage_rrt = coalesce(s_rrt, stage))

all_patients_rrt <- merge (aki_patients, rrt_aki_patients %>% select(ID, stage_rrt), by = "ID", all = TRUE)

all_patients <- all_patients_rrt %>% mutate (stage = coalesce (stage_rrt, stage))

# reclassify all patients with SCr > 353 as stage 3 
highscrpt <- all_patients %>% mutate(stage = ifelse (ID %in% highscr == TRUE, 3, stage))

write_xlsx(highscrpt, "all_patients_early_AKI.xlsx")








