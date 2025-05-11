###PROJECT COLDIA: Analyses for article: Childhood acute leukemia and type 1 diabetes in children:
#a nationwide case-control study###

#Load necessary packages
library(lubridate)
library(dplyr)
library(survival)

#Set working directory
setwd("SET WORKING DIRECTORY HERE")

#--- Leukemia datasets ---

#Import leukemia datasets
leuk_info <- read.csv("./Leukemia_Cases_Cancer_Registry.csv",sep = ";") #leukemia cases
leuk_all <- read.csv("./Leukemia_Cases_Controls_DVV.csv",sep = ";") #leukemia cases and their age- and sex-matched controls

#Exclude lymphoma and chronic leukemia cases based on ICD-10 codes from leuk_info
leuk_info_remove_icd <- subset(leuk_info, iarccrgtools_icd10 %in% c("C835", "C837", "C911", "C914", "C921", "C933", "C929", "C927", "C919", "C918"))
remove_icd_ID <- leuk_info_remove_icd$FID
leuk_info <- subset(leuk_info, !(FID %in% remove_icd_ID))

#Rename columns for consistency
names(leuk_all)[names(leuk_all) == "FID.tutkimushenkil?."] <- "group"  
names(leuk_all)[names(leuk_all) == "FID.verrokki"] <- "ID0" 
names(leuk_all)[names(leuk_all) == "Syntym?.p?.iv?."] <- "bd"
names(leuk_all)[names(leuk_all) == "Sukupuoli"] <- "sex"

#Add binary indicator for leukemia cases (1 = case, 0 = controls) 
case_ids = unique(leuk_all$group)
leuk_all$DV = ifelse(test = leuk_all$ID0 %in% case_ids, 1, 0)

#Recode group variable into consecutive integers
previous_group = leuk_all$group[1]
i = 1
for (j in 2:length(leuk_all$group))
{
  if(leuk_all$group[j] == previous_group){
    leuk_all$group[j] = i
  }else{
    previous_group = leuk_all$group[j]
    i = i + 1
    leuk_all$group[j] = i
  }
}
leuk_all$group[1] = 1

#Exclude matched controls corresponding to previously excluded cases
leuk_all_remove_icd <- subset(leuk_all, ID0 %in% remove_icd_ID)
leuk_all_remove_icd_groups <- leuk_all_remove_icd$group
leuk_all = subset(leuk_all, !(group %in% leuk_all_remove_icd_groups))

#Merge leukemia case-control dataset with Cancer Registry data
leuk_all_info = merge(leuk_all, leuk_info, all.x = TRUE, by.x = "ID0", by.y = "FID")
leuk_all_info = leuk_all_info[order(-leuk_all_info$DV),]
leuk_all_info = leuk_all_info[order(leuk_all_info$group),]

#Convert selected columns to appropriate data types
leuk_all_info$dg_age = as.numeric(sub(",",".",leuk_all_info$dg_age))
leuk_all_info$bd = ymd(leuk_all_info$bd)
leuk_all_info$group = as.integer(leuk_all_info$group)
leuk_all_info$dg_date = ymd(leuk_all_info$dg_date)

#Remove unnecessary columns
leuk_all_info$Verrokin.nro = NULL
leuk_all_info$sex.y = NULL

#Calculate reference diagnosis date for controls basd on case age at diagnosis
case_group = leuk_all_info$group[1]
case_dg_age = leuk_all_info$dg_age[1]
for (j in 2:length(leuk_all_info$group))
{
  if(leuk_all_info$group[j] == case_group){
    leuk_all_info$dg_date[j] = leuk_all_info$bd[j] + days(floor(case_dg_age*365.25))
  }else{
    previous_group = leuk_all_info$group[j]
    case_group = leuk_all_info$group[j]
    case_dg_age = leuk_all_info$dg_age[j]
  }
}

#Rename diagnosis date column to reference date
names(leuk_all_info)[names(leuk_all_info) == "dg_date"] <- "ref.date"

#Remove duplicate patient records, keeping the earlies diagnosis (excluding relapse)
leuk_all_info_ordered = leuk_all_info[order(leuk_all_info$ref.date, decreasing = FALSE),]
duplicated <- leuk_all_info_ordered[duplicated(leuk_all_info_ordered$ID0),]
leuk_all_info = leuk_all_info_ordered[!duplicated(leuk_all_info_ordered$ID0, fromLast = FALSE),]
table(leuk_all_info$DV) #1626 leukemia cases and 4877 controls

#Rename sex variable
names(leuk_all_info)[names(leuk_all_info) == "sex.x"] <- "sex"

#Round age at diagnosis to three decimal places
leuk_all_info$dg_age <- round(leuk_all_info$dg_age,digits = 3)

#Create subgroup datasets for AML and AML
##Define ICD-10 and morphology codes for AML
AML_icd = c("C920", "C930", "C923", "C924", "C925", "C940", "C942", "C928")
AML_morpho = c(9840, 9861, 9866, 9867, 9873, 9874, 9891, 9910)

#Subset AML cases
leuk_AML_case_ <- leuk_all_info[leuk_all_info$iarccrgtools_icd1 %in% AML_icd
                                | leuk_all_info$morpho %in% AML_morpho ,]

groups_AML_case_ <- unique(leuk_AML_case_$group)
AML_ <- leuk_all_info[leuk_all_info$group %in% groups_AML_case_,]

#Convert selected columns to appropriate data types
as.numeric(leuk_all_info$bd)
as.numeric(leuk_all_info$ref.date)

# Calculate age at reference date
leuk_all_info$age <- as.double(difftime(leuk_all_info$ref.date, leuk_all_info$bd, units="days"))/365.25
summary(leuk_all_info$age) #max 17.993

#AML age-specific subgroups
AML01_ <-AML_[AML_$age <1,]
AML19_ <- AML_[AML_$age >=1 & AML_$age < 10,]
AML1017_ <- AML_[AML_$age >=10 & AML_$age < 18,]

#AML sex-specific subgroups
AML_female_ <- AML_[AML_$sex==2,] #2 = female
AML_male_ <- AML_[AML_$sex==1,] #1 = male

#Define ICD-10 and morphology codes for ALL
ALL_icd = c("C910")
ALL_morpho = c(9811, 9812, 9816, 9820, 9835, 9836, 9837)

# Subset ALL cases
leuk_ALL_case_ <- leuk_all_info[leuk_all_info$iarccrgtools_icd1 %in% ALL_icd
                                | leuk_all_info$morpho %in% ALL_morpho ,]
groups_ALL_case_ <- unique(leuk_ALL_case_$group)
ALL_ <- leuk_all_info[leuk_all_info$group %in% groups_ALL_case_,]

#ALL age-specific subgroups
ALL_1.5_5_ <- ALL_[ALL_$age >=1.5 & ALL_$age <6,]
ALL01_ <-ALL_[ALL_$age <1,]
ALL19_ <- ALL_[ALL_$age >=1 & ALL_$age < 10,]
ALL1017_ <- ALL_[ALL_$age >=10 & ALL_$age < 18,]

#ALL sex-specific subgroups
ALL_female_ <- ALL_[ALL_$sex==2,]
ALL_male_ <- ALL_[ALL_$sex==1,]

#Subset remaining acute leukemia cases (non-AML and non-ALL)
other_types <- leuk_all_info %>%
  anti_join(ALL_, by = "ID0") %>%
  anti_join(AML_, by = "ID0")

#--- Diabetes datasets ---

#Import diabetes datasets
dm_info <- read.csv("./Diabetes_Cases_Medical_Reimbursements.csv",sep = ";")
dm_all <- read.csv("./Diabetes_Cases_Controls_DVV.csv",sep = ";")

#Rename columns for consistency
names(dm_all)[names(dm_all) == "FID_Tutkimushenkil?.n_henkil?.tunnus"] <- "group"  
names(dm_all)[names(dm_all) == "FID_verrokki"] <- "ID0" 
names(dm_all)[names(dm_all) == "Syntym?.p?.iv?."] <- "bd"
names(dm_all)[names(dm_all) == "Sukupuoli"] <- "sex"

#Add binary indicator for diabetes cases
case_ids_dm = unique(dm_all$group)
dm_all$DV = ifelse(test = dm_all$ID0 %in% case_ids_dm, 1, 0)

# Restrict to medical reimbursement codes associated with diabetes
new <- dm_info[dm_info$KORVAUSOIKEUS_KOODI %in% c("103", "171", "371", "177", "382"),]
#Dataset new contains reimbursement codes 103, 171, 371, 177 and 382 that refer to insulin medication

#Identify unique case IDs
new_id <- unique(new$FID) #16 967 unique ids

#Sort reimbursement data by approval date
new_ordered <- new %>% arrange(new$KORVAUSOIKEUS_ALPV)

#Retain first reimbursement event for each individual
final_dm_info <- new_ordered %>% distinct(new_ordered$FID, .keep_all = TRUE)

#Recode group variable into consecutive integers
previous_group = dm_all$group[1]
i = 1
for (j in 2:length(dm_all$group))
{
  if(dm_all$group[j] == previous_group){
    dm_all$group[j] = i
  }else{
    previous_group = dm_all$group[j]
    i = i + 1
    dm_all$group[j] = i
  }
}
dm_all$group[1] = 1

#Merge diabetes case information with case-control data but retain dataset only for cases
dm_case_info <- merge(final_dm_info, dm_all, all.x=TRUE, by.x = "FID", by.y ="ID0")
#dm_case_info contains only info of cases, not controls

# Rename diagnosis date column
names(dm_case_info)[names(dm_case_info) == "KORVAUSOIKEUS_ALPV"] <- "dg_date"

#Convert columns to appropriate data types
dm_case_info$bd = ymd(dm_case_info$bd)
dm_case_info$group = as.integer(dm_case_info$group)
dm_case_info$dg_date = ymd(dm_case_info$dg_date)

#Calculate age at diabetes diagnosis
as.numeric(dm_case_info$bd)
as.numeric(dm_case_info$dg_date)
dm_case_info$dg_age <- as.double(difftime(dm_case_info$dg_date, dm_case_info$bd, units="days"))/365.25

summary(dm_case_info$dg_age) #max 38 years

#Exclude cases diagnosed over 18 years of age and their matched controls
remove <- subset(dm_case_info, dg_age > 18) #16 cases over 18 years old
remove_ids <- remove$FID
dm_case_info <- subset(dm_case_info, !(FID %in% remove_ids))

#Rename diagnosis date column to reference date
names(dm_case_info)[names(dm_case_info) == "dg_date"] <- "ref.date"

# Identify and merge matched controls
cases_groups <- dm_case_info$group
dm_controls <- dm_all[dm_all$group %in% cases_groups & dm_all$DV==0,]
dm_cases <- dm_all[dm_all$group %in% cases_groups & dm_all$DV==1,]

#Combine the correct cases and controls
dm_case_controls <- rbind(dm_cases, dm_controls)

#Merge full diabetes case-control dataset
dm_all_info <- merge(dm_case_controls, dm_case_info, by.x = "ID0", by.y = "FID", all.x=TRUE )

#Rename columns for consistency
names(dm_all_info)[names(dm_all_info) == "group.x"] <- "group"
names(dm_all_info)[names(dm_all_info) == "Verrokin.nro.x"] <- "Verrokin.nro"
names(dm_all_info)[names(dm_all_info) == "bd.x"] <- "bd"
names(dm_all_info)[names(dm_all_info) == "sex.x"] <- "sex"
names(dm_all_info)[names(dm_all_info) == "DV.x"] <- "DV"

#Convert birth date and reference date columns to date format
as.numeric(dm_all_info$bd)
as.numeric(dm_all_info$dg_date)
names(dm_all_info)[names(dm_all_info) == "ref.date"] <- "dg_date"
dm_all_info$bd = ymd(dm_all_info$bd)
dm_all_info$dg_date = ymd(dm_all_info$dg_date)

#Sort dataset first by dependent variable (DV: case/control status) and then by matched group identifier
dm_all_info = dm_all_info[order(-dm_all_info$DV),]
dm_all_info = dm_all_info[order(dm_all_info$group),]

#Assign reference diagnosis dates for controls, based on the diagnosis age of the matched case within each group
case_dg_age = dm_all_info$dg_age[1]
for (j in 2:length(dm_all_info$group))
{
  if(dm_all_info$group[j] == case_group){
    dm_all_info$dg_date[j] = dm_all_info$bd[j] + days(floor(case_dg_age*365.25))
  }else{
    previous_group = dm_all_info$group[j]
    case_group = dm_all_info$group[j]
    case_dg_age = dm_all_info$dg_age[j]
  }
}

#Rename 'dg_date' column back to 'ref.date' to indicate the reference date for cases and controls
names(dm_all_info)[names(dm_all_info) == "dg_date"] <- "ref.date"

#Calculate age at reference date (in years) for both cases and controls
dm_all_info$age <- as.double(difftime(dm_all_info$ref.date, dm_all_info$bd, units="days"))/365.25



#Merge leukemia dataset with diabetes dataset by patient identifier (ID0)
leukemia_DM <- leuk_all_info %>%
  left_join(dm_all_info, by = "ID0")

#Rename duplicated columns for clarity: leukemia variables (.x) and diabetes variables (.y)
names(leukemia_DM)[names(leukemia_DM) == "group.x"] <- "group_leuk"  
names(leukemia_DM)[names(leukemia_DM) == "bd.x"] <- "bd_leuk" 
names(leukemia_DM)[names(leukemia_DM) == "sex.x"] <- "sex_leuk"
names(leukemia_DM)[names(leukemia_DM) == "DV.x"] <- "DV_leuk"
names(leukemia_DM)[names(leukemia_DM) == "dg_age.x"] <- "dg_age_leuk"

names(leukemia_DM)[names(leukemia_DM) == "group.y"] <- "group_DM"  
names(leukemia_DM)[names(leukemia_DM) == "bd.y"] <- "bd_DM" 
names(leukemia_DM)[names(leukemia_DM) == "sex.y"] <- "sex_DM"
names(leukemia_DM)[names(leukemia_DM) == "DV.y"] <- "DV_DM"
names(leukemia_DM)[names(leukemia_DM) == "dg_age.y"] <- "dg_age_DM"

#Recode diabetes status: 0 = no type 1 diabetes, 1 = has type 1 diabetes
leukemia_DM$DV_DM <- ifelse(is.na(leukemia_DM$DV_DM), 0 , leukemia_DM$DV_DM)

#Identify individuals with both leukemia and type 1 diabetes
both <- leukemia_DM[leukemia_DM$DV_leuk ==1 & leukemia_DM$DV_DM == 1,] #22 leukemia cases have type 1 diabetes

#--- Removal of individuals with Down syndrome ---

#Import Down syndrome registries
down_8792 <- read.csv("./Down_Finnish_Institute_for_Health_and_Welfare.csv",sep = ";")
down_8719 <-  read.csv("./Down_The_Register_of_Congenital_Malformations.csv",sep = ";")

#Harmonize ID column name across datasets
names(down_8792)[names(down_8792) == "FID"] <- "ID0"  
names(down_8719)[names(down_8719) == "FID"] <- "ID0"

#Mark Down syndrome cases with indicator variable
down_8719$down <- 1
down_8792$down <- 1

#Remove duplicates within each Down syndrome dataset
down_8792 <- distinct(down_8792, ID0, .keep_all =  TRUE)
down_8719 <- distinct(down_8719, ID0, .keep_all =  TRUE)

#Remove unnecessary columns
down_8719$ICD9 <- NULL
down_8719$ICD10 <- NULL
down_8792 <- down_8792[, c("ID0", "down")]

#Combine the two Down syndrome datasets and remove any duplicates
down <- rbind(down_8719, down_8792)
down <- distinct(down, ID0, .keep_all =  TRUE)

#Merge Down syndrome information into leukemia_DM
leukemia_DM <- merge(leukemia_DM,down, by="ID0", all.x=TRUE)
leukemia_DM <- distinct(leukemia_DM, ID0, .keep_all =  TRUE)

#Assess the number of Down syndrome cases among leukemia cases and controls
leukemia_DM_down <- leukemia_DM[which(leukemia_DM$down==1),]
table(leukemia_DM_down$DV_leuk) #63 leukemia cases, 4 controls

#Remove individuals with Down syndrome
leukemia_DM <- leukemia_DM[-which(leukemia_DM$down==1),]
table(leukemia_DM$DV_leuk)

#Remove controls matched to cases with Down syndrome
leukemia_DM_down_cases <- leukemia_DM_down[leukemia_DM_down$DV_leuk==1,]
group_down <- leukemia_DM_down_cases$group_leuk
leukemia_DM <- leukemia_DM[!leukemia_DM$group_leuk %in% group_down,]

#--- Removal of individuals with a history of pancreatitis ---

#Import five datasets covering different time periods
data1 <- read.csv("./panc1.csv",sep = ";")
data2 <- read.csv("/panc2.csv", sep = ";")
data3 <-  read.csv("/panc3.csv",sep = ";")
data4 <-  read.csv("/panc4.csv",sep = ";")
data5 <-  read.csv("/panc5.csv",sep = ";")

#Identify acute pancreatitis cases based on ICD codes in each dataset

#Dataset 1: ICD-10 acute pancreatitis codes
values1 <- c("K85", "K85.3#", "K85.8", "K85.9") #ICD-10 codes indicating acute pancreatitis
panc1 <- data1[data1$ICD10 %in% values1,]

#Dataset 2: Additional ICD-10 codes
columns2 <- c(data2$KOODI)
values2 <- c("K85", "K850", "K853", "K858", "K859") #ICD-10 codes indicating acute pancreatitis
panc2 <- data2[data2$KOODI %in% values2,]

#Dataset 3: ICD-9 codes; no pancreatitis found
columns3 <- c(data3$DG1,data3$DG2, data3$DG3, data3$DG4)
values3 <- c("5770")
panc3 <- data3[columns3 %in% values3,]

#Dataset 4: ICD-9 codes; no pancreatitis found
columns4 <- c(data4$PDG, data4$SDG1, data4$SDG2, data4$SDG3)
values4 <- c("5770E")
panc4 <- data4[columns4 %in% values4,]

#Dataset 5: ICD-9 acute pancreatitis codes
columns5 <- c(data5$PDGO, data5$PDGE, data5$SDG1O, data5$SDG1E,
              data5$SDG2O, data5$SDG2E, data5$PDG, data5$SDG1, data5$SDG2 )
values5 <- c("5770A")
panc5 <- data5[data5$PDG %in% values5,]

#Standardize column names and remove irrelevant columns
names(panc1)[names(panc1) == "FID"] <- "ID0"  
names(panc1)[names(panc1) == "KAYNTI_ALKOI"] <- "dg_date"
names(panc1)[names(panc1) == "ICD10"] <- "dg_code"
panc1$JARJESTYS <- NULL

names(panc2)[names(panc2) == "FID"] <- "ID0"  
names(panc2)[names(panc2) == "TUPVA"] <- "dg_date"
names(panc2)[names(panc2) == "KOODI"] <- "dg_code"
panc2$KENTTA <- NULL
panc2$N <- NULL

names(panc5)[names(panc5) == "FID"] <- "ID0"  
names(panc5)[names(panc5) == "TUPVA"] <- "dg_date"
names(panc5)[names(panc5) == "PDG"] <- "dg_code"
panc5 <- panc5[, c("ID0", "dg_date", "dg_code")]

#Combine pancreatitis datasets (panc3 and panc4 had no relevant cases)
panc <- rbind(panc1, panc2, panc5)

#Mark pancreatitis cases with indicator variable
panc$panc <- 1

#Rename diagnosis date for pancreatitis
names(panc)[names(panc) == "dg_date"] <- "panc_date"

#Extract only date part (remove time) and convert to Date format
panc$panc_date <- sub("[0-9]{2}:[0-9]{2}$","",panc$panc_date)
panc$panc_date <- as.Date(panc$panc_date, format = "%d.%m.%Y")

#Order and de-duplicate pancreatitis dataset based on earliest diagnosis
panc_order = panc[order(format(panc$panc_date, "%Y")),]
panc = panc_order[!duplicated(panc_order$ID0, fromLast = FALSE),]

#Merge pancreatitis status into leukemia_DM
leukemia_DM <- merge(leukemia_DM,panc, by="ID0", all.x=TRUE)

#Assess the number of pancreatitis cases among leukemia cases and controls
leukemia_DM_panc <- leukemia_DM[which(leukemia_DM$panc==1),]
table(leukemia_DM_panc$DV_leuk) #30 leukemia cases, 2 controls

#Remove individuals with pancreatitis
leukemia_DM <- leukemia_DM[-which(leukemia_DM$panc==1),]

#Remove controls matched to cases with pancreatitis
leukemia_DM_panc_cases <- leukemia_DM_panc[leukemia_DM_panc$DV_leuk==1,]
group_panc <- leukemia_DM_panc_cases$group_leuk
leukemia_DM<- leukemia_DM[!leukemia_DM$group_leuk %in% group_panc,]

#--- Conditional logistic regression model ---

#Fit conditional logistic regression model for leukemia outcome
res_leuk <- clogit(DV_leuk ~ DV_DM + strata(leukemia_DM$group_leuk), data = leukemia_DM)
summary(res_leuk)

#Define a function to extract and format ORs and 95% CIs from a logistic model
calculate_conf <- function(logmodel, name) {
  conf = round(exp(confint(logmodel)),2)
  summary = summary(logmodel)
  odds_ratio = round(summary$coefficients[2], 2)
  temp_result = c(name, as.character(odds_ratio), as.character(conf[1]), as.character(conf[2]))
  return(temp_result)
}

#Initialize an empty dataframe to store results
column_names = c("name", "OR", "lower_95", "upper_95")
conf_results = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results) <- column_names

# Create stratified datasets based on sex
# 1 = male
leuk_male <- leukemia_DM[leukemia_DM$sex_leuk==1,]

# 2 = female
leuk_female <- leukemia_DM[leukemia_DM$sex_leuk==2,]

#Calculate age at reference date
as.numeric(leukemia_DM$bd_leuk)
as.numeric(leukemia_DM$ref.date.x)
leukemia_DM$age <- as.double(difftime(leukemia_DM$ref.date.x, leukemia_DM$bd_leuk, units="days"))/365.25

#Create age group datasets
#0-0.99 years old
leuk_0_1 <- leukemia_DM[leukemia_DM$age <1,]
#1-9.99 years old
leuk_1_9 <- leukemia_DM[leukemia_DM$age >=1 & leukemia_DM$age < 10,]
#10-17.99 years old
leuk_1017 <- leukemia_DM[leukemia_DM$age >=10 & leukemia_DM$age < 18,]

#Subset AML cases based on ICD and morphology codes
AML_icd = c("C920", "C930", "C923", "C924", "C925", "C940", "C942", "C928")
AML_morpho = c(9840, 9861, 9866, 9867, 9873, 9874, 9891, 9910)
leuk_AML_case <- leukemia_DM[leukemia_DM$iarccrgtools_icd1 %in% AML_icd
                             | leukemia_DM$morpho %in% AML_morpho ,]
groups_AML_case <- unique(leuk_AML_case$group_leuk)
AML <- leukemia_DM[leukemia_DM$group_leuk %in% groups_AML_case,]

#Stratify AML cases by age and sex
AML01 <-AML[AML$age <1,]
AML19 <- AML[AML$age >=1 & AML$age < 10,]
AML1017 <- AML[AML$age >=10 & AML$age < 18,]
AML_female <- AML[AML$sex_leuk==2,]
AML_male <- AML[AML$sex_leuk==1,]

#Subset ALL cases based on ICD and morphology codes
ALL_icd = c("C910")
ALL_morpho = c(9811, 9812, 9816, 9820, 9835, 9836, 9837)
leuk_ALL_case <- leukemia_DM[leukemia_DM$iarccrgtools_icd1 %in% ALL_icd
                             | leukemia_DM$morpho %in% ALL_morpho ,]
groups_ALL_case <- unique(leuk_ALL_case$group_leuk)
ALL <- leukemia_DM[leukemia_DM$group_leuk %in% groups_ALL_case, ]

#Stratify ALL cases by age and sex
ALL_1.5_5 <- ALL[ALL$age >=1.5 & ALL$age <6,]
ALL01 <-ALL[ALL$age <1,]
ALL19 <- ALL[ALL$age >=1 & ALL$age < 10,]
ALL1017 <- ALL[ALL$age >=10 & ALL$age < 18,]
#ALL sexes
ALL_female <- ALL[ALL$sex_leuk==2,]
ALL_male <- ALL[ALL$sex_leuk==1,]

#Create list of datasets and corresponding names
datasets = list(leukemia_DM, leuk_male, leuk_female, leuk_0_1, leuk_1_9, leuk_1017, AML, AML01, AML19, AML1017, AML_female, AML_male, ALL, ALL_1.5_5, ALL01, ALL19, ALL1017, ALL_female, ALL_male)
names = list("all", "male", "female", "0-0.99", "1-9.99", "10-17.99", "AML","AML 0-0.99", "AML 1-9.99", "AML 10-17.99","AML_female", "AML_male", "ALL", "ALL 1.5-5.99", "ALL 0-0.99", "ALL 1-9.99", "ALL 10-17.99", "ALL_female", "ALL_male")

#Fit models across all datasets
for(i in 1:length(datasets)){
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk <- clogit(DV_leuk ~ DV_DM + strata(dataset$group_leuk), data = dataset)
  conf_results = rbind(conf_results, calculate_conf(res_leuk, name))
  
  #Model assuming diabetes precedes leukemia
  dataset$var <- dataset$DV_DM
  dataset$var[which(dataset$DV_DM == 1 & 
                      (dataset$ref.date.x < dataset$ref.date.y))] <- 0
  res_leuk_order <- clogit(DV_leuk ~ var + strata(dataset$group_leuk), data = dataset)
  conf_results = rbind(conf_results, calculate_conf(res_leuk_order, paste0(name, "_dm_prec")))
  
  # Model assuming leukemia precedes diabetes
  dataset$var2 <- dataset$DV_DM
  dataset$var2[which(dataset$DV_DM == 1 & 
                       (dataset$ref.date.x > dataset$ref.date.y))] <- 0
  res_leuk_order2 <- clogit(DV_leuk ~ var2 + strata(dataset$group_leuk), data = dataset)
  conf_results = rbind(conf_results, calculate_conf(res_leuk_order2, paste0(name, "_leuk_prec")))
}
colnames(conf_results) <- column_names


#Check distribution of ALL cases by age group and sex
ALL_cases <- ALL[ALL$DV_leuk==1,]
age1 <- ALL01[ALL01$DV_leuk==1,] #46
age2 <- ALL19[ALL19$DV_leuk==1,]
age3 <- ALL1017[ALL1017$DV_leuk==1,]
ALL_female <- ALL_cases[ALL_cases$sex_leuk==2,]
ALL_male <- ALL_cases[ALL_cases$sex_leuk==1,]

#Check distribution of AML cases by age group and sex
age4 <- leuk_AML_case[leuk_AML_case$dg_age_leuk <1,] 
age5 <- leuk_AML_case[leuk_AML_case$dg_age_leuk >=1 & leuk_AML_case$dg_age_leuk < 10,]
age6 <- leuk_AML_case[leuk_AML_case$dg_age_leuk >=10 & leuk_AML_case$dg_age_leuk < 18,]
AML_female <- leuk_AML_case[leuk_AML_case$sex_leuk==2,]
AML_male <- leuk_AML_case[leuk_AML_case$sex_leuk==1,]

#Test interaction: does the OR of oldest ALL age group differ significantly from younger age groups?

#Model without interaction term
model1 <- clogit(DV_leuk ~ DV_DM + strata(ALL$group_leuk), data = ALL)

#Model with interaction between diabetes and age group
ALL$age_group <- cut(ALL$age,
                     breaks = c(0, 10,18),
                     labels = c("1", "2"),
                     include.lowest = TRUE,
                     right = FALSE)
model2 <- clogit(DV_leuk ~ DV_DM * as.factor(age_group) + strata(ALL$group_leuk), data = ALL)

#Likelihood ratio test
anova(model1, model2, test = "Chisq") #p=0.02

# Test interaction: does diabetes effect differ by sex in ALL cases?

#Model without interaction term
model1 <- clogit(DV_leuk ~ DV_DM + strata(ALL$group_leuk), data = ALL)

#Model with interaction between diabetes and sex
model2 <- clogit(DV_leuk ~ DV_DM * as.factor(sex_leuk) + strata(ALL$group_leuk), data = ALL)

# Likelihood ratio test
anova(model1, model2, test = "Chisq") #p=0.72

#--- Multivariable analysis based on information from the Medical Birth Register---

#Load supplementary dataset containing birth register information
synre_l <- read.csv("C:/Users/JuliaVentela-b45/Documents/Analyses/data/vanhat_dm/THL/FD_2022_1096_THL_THL2022_1096_synre.csv",sep = ";")

# Rename the identifier column for consistency
names(synre_l)[names(synre_l) == "LAPSI_FID"] <- "ID0"

#Merge birth register information with the primary leukemia dataset
leukemia_DM_synre <- merge(leukemia_DM,synre_l, by="ID0", all.x=TRUE)

# Newborn size analysis

#Identify cases with missing size at birth (SGAC) data
SGAC_na <- leukemia_DM_synre[is.na(leukemia_DM_synre$SGAC),] #1,123
table(SGAC_na$DV_leuk) #Missing information concerns 266 cases and 857 controls

#Identify entries where SGAC is coded as missing (value = 9)
sgac_missing <- leukemia_DM_synre[leukemia_DM_synre$SGAC==9,]
sgac_missing_ <- subset(sgac_missing, !sgac_missing$ID0 == "NA.NA")
table(sgac_missing_$DV_leuk) #8 cases and 18 controls

#Assess birth dates among entries with missing SGAC
before_1987 <- sum(SGAC_na$bd_leuk < as.Date("1987-01-01")) #867
after_1986 <- sum(SGAC_na$bd_leuk >= as.Date("1987-01-01")) #256

#Recode SGAC variable: 1 and 2 = not LGA (0), 3 = LGA (1)
leukemia_DM_synre$SGAC[leukemia_DM_synre$SGAC=="1"] <- "0"
leukemia_DM_synre$SGAC[leukemia_DM_synre$SGAC=="2"] <- "0"
leukemia_DM_synre$SGAC[leukemia_DM_synre$SGAC=="3"] <- "1"
leukemia_DM_synre$SGAC[leukemia_DM_synre$SGAC==9] <- NA

#Calculate number of LGA cases and controls
cases_LGA <- leukemia_DM_synre[leukemia_DM_synre$DV_leuk==1 &
                                 leukemia_DM_synre$SGAC==1,]
cases_LGA <- subset(cases_LGA, !cases_LGA$ID0 == "NA.NA") #63 cases


controls_LGA <- leukemia_DM_synre[leukemia_DM_synre$DV_leuk==0 &
                                 leukemia_DM_synre$SGAC==1,]
controls_LGA <- subset(controls_LGA, !controls_LGA$ID0 == "NA.NA") #126 controls

# Maternal smoking analysis

#Identify missing data for maternal smoking (TUPAKOINTITUNNUS)
smoking_na <- leukemia_DM_synre[is.na(leukemia_DM_synre$TUPAKOINTITUNNUS),] #1,123
table(smoking_na$DV_leuk) #Missing information concerns 266 cases and 857 controls

#Identify smoking status entries coded as missing (value = 9)
smoking_missing <- leukemia_DM_synre[leukemia_DM_synre$TUPAKOINTITUNNUS==9,]
smoking_missing_ <- subset(smoking_missing, !smoking_missing$ID0 == "NA.NA")
table(smoking_missing_$DV_leuk) #30 cases and 87 controls

# Recode maternal smoking status: 
# 2,3,4 = smoked (1); 1 = non-smoking (0); 9 = missing (NA)
leukemia_DM_synre$TUPAKOINTITUNNUS[leukemia_DM_synre$TUPAKOINTITUNNUS==1] <- 0
leukemia_DM_synre$TUPAKOINTITUNNUS[leukemia_DM_synre$TUPAKOINTITUNNUS==2] <- 1
leukemia_DM_synre$TUPAKOINTITUNNUS[leukemia_DM_synre$TUPAKOINTITUNNUS==3] <- 1
leukemia_DM_synre$TUPAKOINTITUNNUS[leukemia_DM_synre$TUPAKOINTITUNNUS==4] <- 1
leukemia_DM_synre$TUPAKOINTITUNNUS[leukemia_DM_synre$TUPAKOINTITUNNUS==9] <- NA

#Calculate number of cases and controls with maternal smoking
cases_smoking <- leukemia_DM_synre[leukemia_DM_synre$DV_leuk==1 &
                                     leukemia_DM_synre$TUPAKOINTITUNNUS==1,]
cases_smoking <- subset(cases_smoking, !cases_smoking$ID0 == "NA.NA") #215 cases

controls_smoking <- leukemia_DM_synre[leukemia_DM_synre$DV_leuk==0 &
                                     leukemia_DM_synre$TUPAKOINTITUNNUS==1,]
controls_smoking <- subset(controls_smoking, !controls_smoking$ID0 == "NA.NA") #547 controls

#Mode of delivery analysis

#Identify missing data for mode of delivery (SYNNYTYSTAPATUNNUS)
synnytys_na <- leukemia_DM_synre[is.na(leukemia_DM_synre$SYNNYTYSTAPATUNNUS),] #1,123
table(synnytys_na$DV_leuk) #Missing information concerns 266 cases and 857 controls

#Identify delivery mode entries coded as missing (value = 9)
birth_missing <- leukemia_DM_synre[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==9,]
birth_missing_ <- subset(birth_missing, !birth_missing$ID0 == "NA.NA")
table(birth_missing_$DV_leuk) #4 cases ja 8 controls

#Recode mode of delivery: 5 = pre-labour caesarean (1), others = 0
leukemia_DM_synre$SYNNYTYSTAPATUNNUS[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==1] <- 0
leukemia_DM_synre$SYNNYTYSTAPATUNNUS[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==2] <- 0
leukemia_DM_synre$SYNNYTYSTAPATUNNUS[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==3] <- 0
leukemia_DM_synre$SYNNYTYSTAPATUNNUS[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==4] <- 0
leukemia_DM_synre$SYNNYTYSTAPATUNNUS[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==5] <- 1
leukemia_DM_synre$SYNNYTYSTAPATUNNUS[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==6] <- 0
leukemia_DM_synre$SYNNYTYSTAPATUNNUS[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==7] <- 0
leukemia_DM_synre$SYNNYTYSTAPATUNNUS[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==8] <- 0
leukemia_DM_synre$SYNNYTYSTAPATUNNUS[leukemia_DM_synre$SYNNYTYSTAPATUNNUS==9] <- NA

#Calculate number of cases and controls with prelabour caesarean section
cases_csection <- leukemia_DM_synre[leukemia_DM_synre$DV_leuk==1 &
                                 leukemia_DM_synre$SYNNYTYSTAPATUNNUS==1,]
cases_section <- subset(cases_csection, !cases_csection$ID0 == "NA.NA") #84

controls_csection <- leukemia_DM_synre[leukemia_DM_synre$DV_leuk==0 &
                                      leukemia_DM_synre$SYNNYTYSTAPATUNNUS==1,]
controls_section <- subset(controls_csection, !controls_csection$ID0 == "NA.NA") #237

# Univariate conditional logistic regression models
#Smoking status adjusted model
res_leuk_adj_smoking <- clogit(DV_leuk ~ DV_DM +  as.character(TUPAKOINTITUNNUS) + strata(leukemia_DM$group_leuk), data = leukemia_DM_synre)
summary(res_leuk_adj_smoking) #OR = 1.2 (95% CI 1.0-1.4)

#Mode of delivery adjusted model
res_leuk_adj_delivery <- clogit(DV_leuk ~  as.character(SYNNYTYSTAPATUNNUS) + strata(leukemia_DM$group_leuk), data = leukemia_DM_synre)
summary(res_leuk_adj_delivery) #OR = 1.1 (95% CI 0.8-1.4)
#Non-significant, exclude from the multivariable analysis

# Newborn size (LGA) adjusted model
res_leuk_adj_LGA <- clogit(DV_leuk ~  as.character(SGAC) + strata(leukemia_DM$group_leuk), data = leukemia_DM_synre)
summary(res_leuk_adj_LGA) #OR = 1.5 (95% CI 1.1-2.0)

#Multivariable model: newborn size and maternal smoking

#Build multivariable model
res_leuk_adj <- clogit(DV_leuk ~  DV_DM + as.character(SGAC) + as.character(TUPAKOINTITUNNUS)
                       + strata(leukemia_DM_synre$group_leuk), data = leukemia_DM_synre)
summary(res_leuk_adj) #OR = 1.7 (95% CI 0.9-3.0)

#Remove cases with missing covariate data
leukemia_DM_synre <- leukemia_DM_synre[!is.na(leukemia_DM_synre$TUPAKOINTITUNNUS),]
leukemia_DM_synre <- leukemia_DM_synre[!is.na(leukemia_DM_synre$SGAC),]

#Function to extract odds ratios and confidence intervals
calculate_conf_2 <- function(logmodel, name) {
  conf = round(exp(confint(logmodel)),2)[1,]
  summary = summary(logmodel)
  odds_ratio = round(exp(summary$coefficients[1]), 2)
  temp_result = c(name, as.character(odds_ratio), as.character(conf[1]), as.character(conf[2]))
  return(temp_result)
}

#Initialize results dataframe
column_names = c("name", "OR", "lower_95", "upper_95")
conf_results_multi = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_multi) <- column_names

#Create stratified datasets based on sex
#1 = male
leuk_male <- leukemia_DM_synre[leukemia_DM_synre$sex_leuk==1,]

#2 = female
leuk_female <- leukemia_DM_synre[leukemia_DM_synre$sex_leuk==2,]

#Calculate age at reference date
as.numeric(leukemia_DM_synre$bd_leuk)
as.numeric(leukemia_DM_synre$ref.date.x)
leukemia_DM_synre$age <- as.double(difftime(leukemia_DM_synre$ref.date.x, leukemia_DM_synre$bd_leuk, units="days"))/365.25

# Create age group datasets
#0-0.99 years old
leuk_0_1 <- leukemia_DM_synre[leukemia_DM_synre$age <1,]

#1-9.99 years old
leuk_1_9 <- leukemia_DM_synre[leukemia_DM_synre$age >=1 & leukemia_DM_synre$age < 10,]

#10-17.99 years old
leuk_1017 <- leukemia_DM_synre[leukemia_DM_synre$age >=10 & leukemia_DM_synre$age < 18,]

#Subset AML cases based on ICD and morphology codes
AML_icd = c("C920", "C930", "C923", "C924", "C925", "C940", "C942", "C928")
AML_morpho = c(9840, 9861, 9866, 9867, 9873, 9874, 9891, 9910)
leuk_AML_case <- leukemia_DM_synre[leukemia_DM_synre$iarccrgtools_icd1 %in% AML_icd
                                   | leukemia_DM_synre$morpho %in% AML_morpho ,]
groups_AML_case <- unique(leuk_AML_case$group_leuk)
AML <- leukemia_DM_synre[leukemia_DM_synre$group_leuk %in% groups_AML_case,]

#Stratify AML cases by age and sex
AML01 <-AML[AML$age <1,]
AML19 <- AML[AML$age >=1 & AML$age < 10,]
AML1017 <- AML[AML$age >=10 & AML$age < 18,]
AML_female <- AML[AML$sex_leuk==2,]
AML_male <- AML[AML$sex_leuk==1,]

#Subset ALL cases based on ICD and morphology codes
ALL_icd = c("C910")
ALL_morpho = c(9811, 9812, 9816, 9820, 9835, 9836, 9837)
leuk_ALL_case <- leukemia_DM_synre[leukemia_DM_synre$iarccrgtools_icd1 %in% ALL_icd
                                   | leukemia_DM_synre$morpho %in% ALL_morpho ,]

groups_ALL_case <- unique(leuk_ALL_case$group_leuk)
ALL <- leukemia_DM_synre[leukemia_DM_synre$group_leuk %in% groups_ALL_case, ]

#Stratify ALL cases by age and sex
ALL_1.5_5 <- ALL[ALL$age >=1.5 & ALL$age <6,]
ALL01 <-ALL[ALL$age <1,]
ALL19 <- ALL[ALL$age >=1 & ALL$age < 10,]
ALL1017 <- ALL[ALL$age >=10 & ALL$age < 18,]
ALL_female <- ALL[ALL$sex_leuk==2,]
ALL_male <- ALL[ALL$sex_leuk==1,]

#Create list of datasets and corresponding names
datasets = list(leukemia_DM_synre, leuk_male, leuk_female, leuk_0_1, leuk_1_9, leuk_1017, AML, AML01, AML19, AML1017, AML_female, ALL_male, ALL, ALL_1.5_5, ALL01, ALL19, ALL1017, ALL_female, ALL_male)
names = list("all", "male", "female", "0-0.99", "1-9.99", "10-17.99", "AML","AML 0-0.99", "AML 1-9.99", "AML 10-17.99","AML_female", "AML_male", "ALL", "ALL 1.5-5.99", "ALL 0-0.99", "ALL 1-9.99", "ALL 10-17.99", "ALL_female", "ALL_male")

#Fit models across all datasets
for(i in 1:length(datasets)){
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk <- clogit(DV_leuk ~ DV_DM + as.character(SGAC) + as.character(TUPAKOINTITUNNUS) + strata(dataset$group_leuk), data = dataset)
  conf_results_multi = rbind(conf_results_multi, calculate_conf_2(res_leuk, name))
  
  #Model assuming diabetes precedes leukemia
  dataset$var <- dataset$DV_DM
  dataset$var[which(dataset$DV_DM == 1 & 
                      (dataset$ref.date.x < dataset$ref.date.y))] <- 0
  res_leuk_order <- clogit(DV_leuk ~ var + as.character(SGAC) + as.character(TUPAKOINTITUNNUS)+ strata(dataset$group_leuk), data = dataset)
  conf_results_multi = rbind(conf_results_multi, calculate_conf_2(res_leuk_order, paste0(name, "_dm_prec")))
  
  #Model assuming leukemia precedes diabetes
  dataset$var2 <- dataset$DV_DM
  dataset$var2[which(dataset$DV_DM == 1 & 
                       (dataset$ref.date.x > dataset$ref.date.y))] <- 0
  
  res_leuk_order2 <- clogit(DV_leuk ~ var2 + as.character(SGAC) + as.character(TUPAKOINTITUNNUS)+ strata(dataset$group_leuk), data = dataset)
  conf_results_multi = rbind(conf_results_multi, calculate_conf_2(res_leuk_order2, paste0(name, "_leuk_prec")))
}
colnames(conf_results_multi) <- column_names