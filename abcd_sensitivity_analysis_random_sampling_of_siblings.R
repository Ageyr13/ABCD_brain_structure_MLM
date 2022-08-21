#Title: Sensitivity Analysis - Randomly selecting 1 sib per family - ABCD Release 4.0
#Date: 08/19/2022

#clear workspace
rm(list=ls())

##Load libraries, read files, data management

#load packages
library(psych)              # use for descriptive statistics
library(ggplot2)            # use for plots
library(lme4)               # fits mixed models
library(lmerTest)           # provides t-tests for fixed effects
library(dplyr)              # load dplyr package
library(plyr)               # load plyr package
library(patchwork)          # use for combining figures

#set working directory
setwd("V:/ABCD/Data/Release4.0/Analysis/MLM/updated")

#import wide dataset
pfactor_wide <- read.csv("abcd_brain_structure_factor_scores_wide.csv",header=TRUE)
#import long dataset
pfactor <- read.csv("abcd_brain_structure_factor_scores_long.csv",header=TRUE)

#recode all missing values from -99 to NA
pfactor_wide[pfactor_wide == -99] <- NA
pfactor[pfactor == -99] <- NA

#recode wave from values of 1,2,3 to 0,1,2
pfactor$wave <- pfactor$Wave1 - 1

#create new variable renaming nums to id
pfactor$id <- pfactor$nums
pfactor_wide$id <- pfactor_wide$ï..nums

#create new variables renaming TICs of no interest
pfactor$age <- pfactor$interview_age_baseline
pfactor$sex <- pfactor$sex_coded

pfactor_wide$age <- pfactor_wide$interview_age_baseline
pfactor_wide$sex <- pfactor_wide$sex_coded

#rescale and rename global brain variables
pfactor$wb_cort_vol <- pfactor$smri_vol_cdk_total/10000
pfactor$meanwb_ct <- pfactor$smri_thick_cdk_mean
pfactor$wb_cort_area <- pfactor$smri_area_cdk_total/10000
pfactor$subcort_vol <- pfactor$smri_vol_scs_subcorticalgv/1000

pfactor_wide$wb_cort_vol <- pfactor_wide$smri_vol_cdk_total/10000
pfactor_wide$meanwb_ct <- pfactor_wide$smri_thick_cdk_mean
pfactor_wide$wb_cort_area <- pfactor_wide$smri_area_cdk_total/10000
pfactor_wide$subcort_vol <- pfactor_wide$smri_vol_scs_subcorticalgv/1000

#create dataframe of siblings only in the long format
pfactor_long_sibs <- pfactor[which(pfactor$siblings!=0 & pfactor$cbcl_demo_inclusion_use!=0 & pfactor$smri_inclusion_use!=0 &
                                     pfactor$wave!="NA" & pfactor$p_ho!="NA" & pfactor$smri_area_cdk_banksstslh!="NA"),]
#create dataframe of singletons only in the long format
pfactor_long_nosibs <- pfactor[which(pfactor$siblings!=1 & pfactor$cbcl_demo_inclusion_use!=0 & pfactor$smri_inclusion_use!=0 &
                                       pfactor$wave!="NA" & pfactor$p_ho!="NA" & pfactor$smri_area_cdk_banksstslh!="NA"),]
#create dataframe of siblings only in the wide format
pfactor_wide_sibs <- pfactor_wide[which(pfactor_wide$siblings!=0 & pfactor_wide$cbcl_demo_inclusion_use!=0 & pfactor_wide$smri_inclusion_use!=0),]

#randomly generate 100 columns of numbers 1 to 3745 (n=3745 siblings)
for (i in 1:100){
  sibs[i] <- sample(1:3745, 3745, replace=F)
}
write.csv(sibs, "100_random_samples_of_sibs.csv") #write csv file

#combine 100 columns of random numbers with the siblings only wide dataframe
sibs <- read.csv("100_random_samples_of_sibs.csv",header=TRUE)
pfactor_wide_sibs <- cbind(pfactor_wide_sibs, sibs)
pfactor_wide_sibs$V1 <- pfactor_wide_sibs$matrix.unlist.sibs...nrow...length.sibs...byrow...TRUE.
 
#create 100 dataframes selecting only the sibling with the maximum number within their family
V1 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V1==max(V1))

#reduce dataframe by selecting only id and random samples columns
V1_reduced <- subset(V1, select=c("id", "V1"))

#merge into long dataset to get data ready for analysis
pfactor_long_sibs_V1 <- merge(pfactor_long_sibs, V1_reduced)
pfactor_long_nosibs[, 'V1'] <- NA
pfactor_long_sibs_nosibs_V1 <- rbind(pfactor_long_sibs_V1, pfactor_long_nosibs) # data is now ready for analysis

V2 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V2==max(V2))
V2_reduced <- subset(V2, select=c("id", "V2"))

pfactor_long_sibs_V2 <- merge(pfactor_long_sibs, V2_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V2 = V1)
pfactor_long_sibs_nosibs_V2 <- rbind(pfactor_long_sibs_V2, pfactor_long_nosibs)

V3 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V3==max(V3))
V3_reduced <- subset(V3, select=c("id", "V3"))

pfactor_long_sibs_V3 <- merge(pfactor_long_sibs, V3_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V3 = V2)
pfactor_long_sibs_nosibs_V3 <- rbind(pfactor_long_sibs_V3, pfactor_long_nosibs)

V4 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V4==max(V4))
V4_reduced <- subset(V4, select=c("id", "V4"))

pfactor_long_sibs_V4 <- merge(pfactor_long_sibs, V4_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V4 = V3)
pfactor_long_sibs_nosibs_V4 <- rbind(pfactor_long_sibs_V4, pfactor_long_nosibs)

V5 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V5==max(V5))
V5_reduced <- subset(V5, select=c("id", "V5"))

pfactor_long_sibs_V5 <- merge(pfactor_long_sibs, V5_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V5 = V4)
pfactor_long_sibs_nosibs_V5 <- rbind(pfactor_long_sibs_V5, pfactor_long_nosibs)

V6 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V6==max(V6))
V6_reduced <- subset(V6, select=c("id", "V6"))

pfactor_long_sibs_V6 <- merge(pfactor_long_sibs, V6_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V6 = V5)
pfactor_long_sibs_nosibs_V6 <- rbind(pfactor_long_sibs_V6, pfactor_long_nosibs)

V7 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V7==max(V7))
V7_reduced <- subset(V7, select=c("id", "V7"))

pfactor_long_sibs_V7 <- merge(pfactor_long_sibs, V7_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V7 = V6)
pfactor_long_sibs_nosibs_V7 <- rbind(pfactor_long_sibs_V7, pfactor_long_nosibs)

V8 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V8==max(V8))
V8_reduced <- subset(V8, select=c("id", "V8"))

pfactor_long_sibs_V8 <- merge(pfactor_long_sibs, V8_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V8 = V7)
pfactor_long_sibs_nosibs_V8 <- rbind(pfactor_long_sibs_V8, pfactor_long_nosibs)

V9 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V9==max(V9))
V9_reduced <- subset(V9, select=c("id", "V9"))

pfactor_long_sibs_V9 <- merge(pfactor_long_sibs, V9_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V9 = V8)
pfactor_long_sibs_nosibs_V9 <- rbind(pfactor_long_sibs_V9, pfactor_long_nosibs)

V10 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V10==max(V10))
V10_reduced <- subset(V10, select=c("id", "V10"))

pfactor_long_sibs_V10 <- merge(pfactor_long_sibs, V10_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V10 = V9)
pfactor_long_sibs_nosibs_V10 <- rbind(pfactor_long_sibs_V10, pfactor_long_nosibs)

V11 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V11==max(V11))
V11_reduced <- subset(V11, select=c("id", "V11"))

pfactor_long_sibs_V11 <- merge(pfactor_long_sibs, V11_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V11 = V10)
pfactor_long_sibs_nosibs_V11 <- rbind(pfactor_long_sibs_V11, pfactor_long_nosibs)

V12 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V12==max(V12))
V12_reduced <- subset(V12, select=c("id", "V12"))

pfactor_long_sibs_V12 <- merge(pfactor_long_sibs, V12_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V12 = V11)
pfactor_long_sibs_nosibs_V12 <- rbind(pfactor_long_sibs_V12, pfactor_long_nosibs)

V13 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V13==max(V13))
V13_reduced <- subset(V13, select=c("id", "V13"))

pfactor_long_sibs_V13 <- merge(pfactor_long_sibs, V13_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V13 = V12)
pfactor_long_sibs_nosibs_V13 <- rbind(pfactor_long_sibs_V13, pfactor_long_nosibs)

V14 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V14==max(V14))
V14_reduced <- subset(V14, select=c("id", "V14"))

pfactor_long_sibs_V14 <- merge(pfactor_long_sibs, V14_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V14 = V13)
pfactor_long_sibs_nosibs_V14 <- rbind(pfactor_long_sibs_V14, pfactor_long_nosibs)

V15 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V15==max(V15))
V15_reduced <- subset(V15, select=c("id", "V15"))

pfactor_long_sibs_V15 <- merge(pfactor_long_sibs, V15_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V15 = V14)
pfactor_long_sibs_nosibs_V15 <- rbind(pfactor_long_sibs_V15, pfactor_long_nosibs)

V16 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V16==max(V16))
V16_reduced <- subset(V16, select=c("id", "V16"))

pfactor_long_sibs_V16 <- merge(pfactor_long_sibs, V16_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V16 = V15)
pfactor_long_sibs_nosibs_V16 <- rbind(pfactor_long_sibs_V16, pfactor_long_nosibs)

V17 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V17==max(V17))
V17_reduced <- subset(V17, select=c("id", "V17"))

pfactor_long_sibs_V17 <- merge(pfactor_long_sibs, V17_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V17 = V16)
pfactor_long_sibs_nosibs_V17 <- rbind(pfactor_long_sibs_V17, pfactor_long_nosibs)

V18 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V18==max(V18))
V18_reduced <- subset(V18, select=c("id", "V18"))

pfactor_long_sibs_V18 <- merge(pfactor_long_sibs, V18_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V18 = V17)
pfactor_long_sibs_nosibs_V18 <- rbind(pfactor_long_sibs_V18, pfactor_long_nosibs)

V19 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V19==max(V19))
V19_reduced <- subset(V19, select=c("id", "V19"))

pfactor_long_sibs_V19 <- merge(pfactor_long_sibs, V19_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V19 = V18)
pfactor_long_sibs_nosibs_V19 <- rbind(pfactor_long_sibs_V19, pfactor_long_nosibs)

V20 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V20==max(V20))
V20_reduced <- subset(V20, select=c("id", "V20"))

pfactor_long_sibs_V20 <- merge(pfactor_long_sibs, V20_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V20 = V19)
pfactor_long_sibs_nosibs_V20 <- rbind(pfactor_long_sibs_V20, pfactor_long_nosibs)

V21 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V21==max(V21))
V21_reduced <- subset(V21, select=c("id", "V21"))

pfactor_long_sibs_V21 <- merge(pfactor_long_sibs, V21_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V21 = V20)
pfactor_long_sibs_nosibs_V21 <- rbind(pfactor_long_sibs_V21, pfactor_long_nosibs)

V22 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V22==max(V22))
V22_reduced <- subset(V22, select=c("id", "V22"))

pfactor_long_sibs_V22 <- merge(pfactor_long_sibs, V22_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V22 = V21)
pfactor_long_sibs_nosibs_V22 <- rbind(pfactor_long_sibs_V22, pfactor_long_nosibs)

V23 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V23==max(V23))
V23_reduced <- subset(V23, select=c("id", "V23"))

pfactor_long_sibs_V23 <- merge(pfactor_long_sibs, V23_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V23 = V22)
pfactor_long_sibs_nosibs_V23 <- rbind(pfactor_long_sibs_V23, pfactor_long_nosibs)

V24 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V24==max(V24))
V24_reduced <- subset(V24, select=c("id", "V24"))

pfactor_long_sibs_V24 <- merge(pfactor_long_sibs, V24_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V24 = V23)
pfactor_long_sibs_nosibs_V24 <- rbind(pfactor_long_sibs_V24, pfactor_long_nosibs)

V25 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V25==max(V25))
V25_reduced <- subset(V25, select=c("id", "V25"))

pfactor_long_sibs_V25 <- merge(pfactor_long_sibs, V25_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V25 = V24)
pfactor_long_sibs_nosibs_V25 <- rbind(pfactor_long_sibs_V25, pfactor_long_nosibs)

V26 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V26==max(V26))
V26_reduced <- subset(V26, select=c("id", "V26"))

pfactor_long_sibs_V26 <- merge(pfactor_long_sibs, V26_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V26 = V25)
pfactor_long_sibs_nosibs_V26 <- rbind(pfactor_long_sibs_V26, pfactor_long_nosibs)

V27 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V27==max(V27))
V27_reduced <- subset(V27, select=c("id", "V27"))

pfactor_long_sibs_V27 <- merge(pfactor_long_sibs, V27_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V27 = V26)
pfactor_long_sibs_nosibs_V27 <- rbind(pfactor_long_sibs_V27, pfactor_long_nosibs)

V28 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V28==max(V28))
V28_reduced <- subset(V28, select=c("id", "V28"))

pfactor_long_sibs_V28 <- merge(pfactor_long_sibs, V28_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V28 = V27)
pfactor_long_sibs_nosibs_V28 <- rbind(pfactor_long_sibs_V28, pfactor_long_nosibs)

V29 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V29==max(V29))
V29_reduced <- subset(V29, select=c("id", "V29"))

pfactor_long_sibs_V29 <- merge(pfactor_long_sibs, V29_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V29 = V28)
pfactor_long_sibs_nosibs_V29 <- rbind(pfactor_long_sibs_V29, pfactor_long_nosibs)

V30 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V30==max(V30))
V30_reduced <- subset(V30, select=c("id", "V30"))

pfactor_long_sibs_V30 <- merge(pfactor_long_sibs, V30_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V30 = V29)
pfactor_long_sibs_nosibs_V30 <- rbind(pfactor_long_sibs_V30, pfactor_long_nosibs)

V31 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V31==max(V31))
V31_reduced <- subset(V31, select=c("id", "V31"))

pfactor_long_sibs_V31 <- merge(pfactor_long_sibs, V31_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V31 = V30)
pfactor_long_sibs_nosibs_V31 <- rbind(pfactor_long_sibs_V31, pfactor_long_nosibs)

V32 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V32==max(V32))
V32_reduced <- subset(V32, select=c("id", "V32"))

pfactor_long_sibs_V32 <- merge(pfactor_long_sibs, V32_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V32 = V31)
pfactor_long_sibs_nosibs_V32 <- rbind(pfactor_long_sibs_V32, pfactor_long_nosibs)

V33 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V33==max(V33))
V33_reduced <- subset(V33, select=c("id", "V33"))

pfactor_long_sibs_V33 <- merge(pfactor_long_sibs, V33_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V33 = V32)
pfactor_long_sibs_nosibs_V33 <- rbind(pfactor_long_sibs_V33, pfactor_long_nosibs)

V34 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V34==max(V34))
V34_reduced <- subset(V34, select=c("id", "V34"))

pfactor_long_sibs_V34 <- merge(pfactor_long_sibs, V34_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V34 = V33)
pfactor_long_sibs_nosibs_V34 <- rbind(pfactor_long_sibs_V34, pfactor_long_nosibs)

V35 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V35==max(V35))
V35_reduced <- subset(V35, select=c("id", "V35"))

pfactor_long_sibs_V35 <- merge(pfactor_long_sibs, V35_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V35 = V34)
pfactor_long_sibs_nosibs_V35 <- rbind(pfactor_long_sibs_V35, pfactor_long_nosibs)

V36 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V36==max(V36))
V36_reduced <- subset(V36, select=c("id", "V36"))

pfactor_long_sibs_V36 <- merge(pfactor_long_sibs, V36_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V36 = V35)
pfactor_long_sibs_nosibs_V36 <- rbind(pfactor_long_sibs_V36, pfactor_long_nosibs)

V37 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V37==max(V37))
V37_reduced <- subset(V37, select=c("id", "V37"))

pfactor_long_sibs_V37 <- merge(pfactor_long_sibs, V37_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V37 = V36)
pfactor_long_sibs_nosibs_V37 <- rbind(pfactor_long_sibs_V37, pfactor_long_nosibs)

V38 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V38==max(V38))
V38_reduced <- subset(V38, select=c("id", "V38"))

pfactor_long_sibs_V38 <- merge(pfactor_long_sibs, V38_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V38 = V37)
pfactor_long_sibs_nosibs_V38 <- rbind(pfactor_long_sibs_V38, pfactor_long_nosibs)

V39 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V39==max(V39))
V39_reduced <- subset(V39, select=c("id", "V39"))

pfactor_long_sibs_V39 <- merge(pfactor_long_sibs, V39_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V39 = V38)
pfactor_long_sibs_nosibs_V39 <- rbind(pfactor_long_sibs_V39, pfactor_long_nosibs)

V40 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V40==max(V40))
V40_reduced <- subset(V40, select=c("id", "V40"))

pfactor_long_sibs_V40 <- merge(pfactor_long_sibs, V40_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V40 = V39)
pfactor_long_sibs_nosibs_V40 <- rbind(pfactor_long_sibs_V40, pfactor_long_nosibs)

V41 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V41==max(V41))
V41_reduced <- subset(V41, select=c("id", "V41"))

pfactor_long_sibs_V41 <- merge(pfactor_long_sibs, V41_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V41 = V40)
pfactor_long_sibs_nosibs_V41 <- rbind(pfactor_long_sibs_V41, pfactor_long_nosibs)

V42 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V42==max(V42))
V42_reduced <- subset(V42, select=c("id", "V42"))

pfactor_long_sibs_V42 <- merge(pfactor_long_sibs, V42_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V42 = V41)
pfactor_long_sibs_nosibs_V42 <- rbind(pfactor_long_sibs_V42, pfactor_long_nosibs)

V43 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V43==max(V43))
V43_reduced <- subset(V43, select=c("id", "V43"))

pfactor_long_sibs_V43 <- merge(pfactor_long_sibs, V43_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V43 = V42)
pfactor_long_sibs_nosibs_V43 <- rbind(pfactor_long_sibs_V43, pfactor_long_nosibs)

V44 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V44==max(V44))
V44_reduced <- subset(V44, select=c("id", "V44"))

pfactor_long_sibs_V44 <- merge(pfactor_long_sibs, V44_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V44 = V43)
pfactor_long_sibs_nosibs_V44 <- rbind(pfactor_long_sibs_V44, pfactor_long_nosibs)

V45 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V45==max(V45))
V45_reduced <- subset(V45, select=c("id", "V45"))

pfactor_long_sibs_V45 <- merge(pfactor_long_sibs, V45_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V45 = V44)
pfactor_long_sibs_nosibs_V45 <- rbind(pfactor_long_sibs_V45, pfactor_long_nosibs)

V46 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V46==max(V46))
V46_reduced <- subset(V46, select=c("id", "V46"))

pfactor_long_sibs_V46 <- merge(pfactor_long_sibs, V46_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V46 = V45)
pfactor_long_sibs_nosibs_V46 <- rbind(pfactor_long_sibs_V46, pfactor_long_nosibs)

V47 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V47==max(V47))
V47_reduced <- subset(V47, select=c("id", "V47"))

pfactor_long_sibs_V47 <- merge(pfactor_long_sibs, V47_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V47 = V46)
pfactor_long_sibs_nosibs_V47 <- rbind(pfactor_long_sibs_V47, pfactor_long_nosibs)

V48 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V48==max(V48))
V48_reduced <- subset(V48, select=c("id", "V48"))

pfactor_long_sibs_V48 <- merge(pfactor_long_sibs, V48_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V48 = V47)
pfactor_long_sibs_nosibs_V48 <- rbind(pfactor_long_sibs_V48, pfactor_long_nosibs)

V49 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V49==max(V49))
V49_reduced <- subset(V49, select=c("id", "V49"))

pfactor_long_sibs_V49 <- merge(pfactor_long_sibs, V49_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V49 = V48)
pfactor_long_sibs_nosibs_V49 <- rbind(pfactor_long_sibs_V49, pfactor_long_nosibs)

V50 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V50==max(V50))
V50_reduced <- subset(V50, select=c("id", "V50"))

pfactor_long_sibs_V50 <- merge(pfactor_long_sibs, V50_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V50 = V49)
pfactor_long_sibs_nosibs_V50 <- rbind(pfactor_long_sibs_V50, pfactor_long_nosibs)

V51 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V51==max(V51))
V51_reduced <- subset(V51, select=c("id", "V51"))

pfactor_long_sibs_V51 <- merge(pfactor_long_sibs, V51_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V51 = V50)
pfactor_long_sibs_nosibs_V51 <- rbind(pfactor_long_sibs_V51, pfactor_long_nosibs)

V52 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V52==max(V52))
V52_reduced <- subset(V52, select=c("id", "V52"))

pfactor_long_sibs_V52 <- merge(pfactor_long_sibs, V52_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V52 = V51)
pfactor_long_sibs_nosibs_V52 <- rbind(pfactor_long_sibs_V52, pfactor_long_nosibs)

V53 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V53==max(V53))
V53_reduced <- subset(V53, select=c("id", "V53"))

pfactor_long_sibs_V53 <- merge(pfactor_long_sibs, V53_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V53 = V52)
pfactor_long_sibs_nosibs_V53 <- rbind(pfactor_long_sibs_V53, pfactor_long_nosibs)

V54 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V54==max(V54))
V54_reduced <- subset(V54, select=c("id", "V54"))

pfactor_long_sibs_V54 <- merge(pfactor_long_sibs, V54_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V54 = V53)
pfactor_long_sibs_nosibs_V54 <- rbind(pfactor_long_sibs_V54, pfactor_long_nosibs)

V55 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V55==max(V55))
V55_reduced <- subset(V55, select=c("id", "V55"))

pfactor_long_sibs_V55 <- merge(pfactor_long_sibs, V55_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V55 = V54)
pfactor_long_sibs_nosibs_V55 <- rbind(pfactor_long_sibs_V55, pfactor_long_nosibs)

V56 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V56==max(V56))
V56_reduced <- subset(V56, select=c("id", "V56"))

pfactor_long_sibs_V56 <- merge(pfactor_long_sibs, V56_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V56 = V55)
pfactor_long_sibs_nosibs_V56 <- rbind(pfactor_long_sibs_V56, pfactor_long_nosibs)

V57 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V57==max(V57))
V57_reduced <- subset(V57, select=c("id", "V57"))

pfactor_long_sibs_V57 <- merge(pfactor_long_sibs, V57_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V57 = V56)
pfactor_long_sibs_nosibs_V57 <- rbind(pfactor_long_sibs_V57, pfactor_long_nosibs)

V58 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V58==max(V58))
V58_reduced <- subset(V58, select=c("id", "V58"))

pfactor_long_sibs_V58 <- merge(pfactor_long_sibs, V58_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V58 = V57)
pfactor_long_sibs_nosibs_V58 <- rbind(pfactor_long_sibs_V58, pfactor_long_nosibs)

V59 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V59==max(V59))
V59_reduced <- subset(V59, select=c("id", "V59"))

pfactor_long_sibs_V59 <- merge(pfactor_long_sibs, V59_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V59 = V58)
pfactor_long_sibs_nosibs_V59 <- rbind(pfactor_long_sibs_V59, pfactor_long_nosibs)

V60 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V60==max(V60))
V60_reduced <- subset(V60, select=c("id", "V60"))

pfactor_long_sibs_V60 <- merge(pfactor_long_sibs, V60_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V60 = V59)
pfactor_long_sibs_nosibs_V60 <- rbind(pfactor_long_sibs_V60, pfactor_long_nosibs)

V61 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V61==max(V61))
V61_reduced <- subset(V61, select=c("id", "V61"))

pfactor_long_sibs_V61 <- merge(pfactor_long_sibs, V61_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V61 = V60)
pfactor_long_sibs_nosibs_V61 <- rbind(pfactor_long_sibs_V61, pfactor_long_nosibs)

V62 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V62==max(V62))
V62_reduced <- subset(V62, select=c("id", "V62"))

pfactor_long_sibs_V62 <- merge(pfactor_long_sibs, V62_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V62 = V61)
pfactor_long_sibs_nosibs_V62 <- rbind(pfactor_long_sibs_V62, pfactor_long_nosibs)

V63 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V63==max(V63))
V63_reduced <- subset(V63, select=c("id", "V63"))

pfactor_long_sibs_V63 <- merge(pfactor_long_sibs, V63_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V63 = V62)
pfactor_long_sibs_nosibs_V63 <- rbind(pfactor_long_sibs_V63, pfactor_long_nosibs)

V64 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V64==max(V64))
V64_reduced <- subset(V64, select=c("id", "V64"))

pfactor_long_sibs_V64 <- merge(pfactor_long_sibs, V64_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V64 = V63)
pfactor_long_sibs_nosibs_V64 <- rbind(pfactor_long_sibs_V64, pfactor_long_nosibs)

V65 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V65==max(V65))
V65_reduced <- subset(V65, select=c("id", "V65"))

pfactor_long_sibs_V65 <- merge(pfactor_long_sibs, V65_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V65 = V64)
pfactor_long_sibs_nosibs_V65 <- rbind(pfactor_long_sibs_V65, pfactor_long_nosibs)

V66 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V66==max(V66))
V66_reduced <- subset(V66, select=c("id", "V66"))

pfactor_long_sibs_V66 <- merge(pfactor_long_sibs, V66_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V66 = V65)
pfactor_long_sibs_nosibs_V66 <- rbind(pfactor_long_sibs_V66, pfactor_long_nosibs)

V67 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V67==max(V67))
V67_reduced <- subset(V67, select=c("id", "V67"))

pfactor_long_sibs_V67 <- merge(pfactor_long_sibs, V67_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V67 = V66)
pfactor_long_sibs_nosibs_V67 <- rbind(pfactor_long_sibs_V67, pfactor_long_nosibs)

V68 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V68==max(V68))
V68_reduced <- subset(V68, select=c("id", "V68"))

pfactor_long_sibs_V68 <- merge(pfactor_long_sibs, V68_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V68 = V67)
pfactor_long_sibs_nosibs_V68 <- rbind(pfactor_long_sibs_V68, pfactor_long_nosibs)

V69 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V69==max(V69))
V69_reduced <- subset(V69, select=c("id", "V69"))

pfactor_long_sibs_V69 <- merge(pfactor_long_sibs, V69_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V69 = V68)
pfactor_long_sibs_nosibs_V69 <- rbind(pfactor_long_sibs_V69, pfactor_long_nosibs)

V70 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V70==max(V70))
V70_reduced <- subset(V70, select=c("id", "V70"))

pfactor_long_sibs_V70 <- merge(pfactor_long_sibs, V70_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V70 = V69)
pfactor_long_sibs_nosibs_V70 <- rbind(pfactor_long_sibs_V70, pfactor_long_nosibs)

V71 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V71==max(V71))
V71_reduced <- subset(V71, select=c("id", "V71"))

pfactor_long_sibs_V71 <- merge(pfactor_long_sibs, V71_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V71 = V70)
pfactor_long_sibs_nosibs_V71 <- rbind(pfactor_long_sibs_V71, pfactor_long_nosibs)

V72 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V72==max(V72))
V72_reduced <- subset(V72, select=c("id", "V72"))

pfactor_long_sibs_V72 <- merge(pfactor_long_sibs, V72_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V72 = V71)
pfactor_long_sibs_nosibs_V72 <- rbind(pfactor_long_sibs_V72, pfactor_long_nosibs)

V73 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V73==max(V73))
V73_reduced <- subset(V73, select=c("id", "V73"))

pfactor_long_sibs_V73 <- merge(pfactor_long_sibs, V73_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V73 = V72)
pfactor_long_sibs_nosibs_V73 <- rbind(pfactor_long_sibs_V73, pfactor_long_nosibs)

V74 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V74==max(V74))
V74_reduced <- subset(V74, select=c("id", "V74"))

pfactor_long_sibs_V74 <- merge(pfactor_long_sibs, V74_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V74 = V73)
pfactor_long_sibs_nosibs_V74 <- rbind(pfactor_long_sibs_V74, pfactor_long_nosibs)

V75 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V75==max(V75))
V75_reduced <- subset(V75, select=c("id", "V75"))

pfactor_long_sibs_V75 <- merge(pfactor_long_sibs, V75_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V75 = V74)
pfactor_long_sibs_nosibs_V75 <- rbind(pfactor_long_sibs_V75, pfactor_long_nosibs)

V76 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V76==max(V76))
V76_reduced <- subset(V76, select=c("id", "V76"))

pfactor_long_sibs_V76 <- merge(pfactor_long_sibs, V76_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V76 = V75)
pfactor_long_sibs_nosibs_V76 <- rbind(pfactor_long_sibs_V76, pfactor_long_nosibs)

V77 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V77==max(V77))
V77_reduced <- subset(V77, select=c("id", "V77"))

pfactor_long_sibs_V77 <- merge(pfactor_long_sibs, V77_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V77 = V76)
pfactor_long_sibs_nosibs_V77 <- rbind(pfactor_long_sibs_V77, pfactor_long_nosibs)

V78 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V78==max(V78))
V78_reduced <- subset(V78, select=c("id", "V78"))

pfactor_long_sibs_V78 <- merge(pfactor_long_sibs, V78_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V78 = V77)
pfactor_long_sibs_nosibs_V78 <- rbind(pfactor_long_sibs_V78, pfactor_long_nosibs)

V79 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V79==max(V79))
V79_reduced <- subset(V79, select=c("id", "V79"))

pfactor_long_sibs_V79 <- merge(pfactor_long_sibs, V79_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V79 = V78)
pfactor_long_sibs_nosibs_V79 <- rbind(pfactor_long_sibs_V79, pfactor_long_nosibs)

V80 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V80==max(V80))
V80_reduced <- subset(V80, select=c("id", "V80"))

pfactor_long_sibs_V80 <- merge(pfactor_long_sibs, V80_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V80 = V79)
pfactor_long_sibs_nosibs_V80 <- rbind(pfactor_long_sibs_V80, pfactor_long_nosibs)

V81 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V81==max(V81))
V81_reduced <- subset(V81, select=c("id", "V81"))

pfactor_long_sibs_V81 <- merge(pfactor_long_sibs, V81_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V81 = V80)
pfactor_long_sibs_nosibs_V81 <- rbind(pfactor_long_sibs_V81, pfactor_long_nosibs)

V82 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V82==max(V82))
V82_reduced <- subset(V82, select=c("id", "V82"))

pfactor_long_sibs_V82 <- merge(pfactor_long_sibs, V82_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V82 = V81)
pfactor_long_sibs_nosibs_V82 <- rbind(pfactor_long_sibs_V82, pfactor_long_nosibs)

V83 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V83==max(V83))
V83_reduced <- subset(V83, select=c("id", "V83"))

pfactor_long_sibs_V83 <- merge(pfactor_long_sibs, V83_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V83 = V82)
pfactor_long_sibs_nosibs_V83 <- rbind(pfactor_long_sibs_V83, pfactor_long_nosibs)

V84 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V84==max(V84))
V84_reduced <- subset(V84, select=c("id", "V84"))

pfactor_long_sibs_V84 <- merge(pfactor_long_sibs, V84_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V84 = V83)
pfactor_long_sibs_nosibs_V84 <- rbind(pfactor_long_sibs_V84, pfactor_long_nosibs)

V85 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V85==max(V85))
V85_reduced <- subset(V85, select=c("id", "V85"))

pfactor_long_sibs_V85 <- merge(pfactor_long_sibs, V85_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V85 = V84)
pfactor_long_sibs_nosibs_V85 <- rbind(pfactor_long_sibs_V85, pfactor_long_nosibs)

V86 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V86==max(V86))
V86_reduced <- subset(V86, select=c("id", "V86"))

pfactor_long_sibs_V86 <- merge(pfactor_long_sibs, V86_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V86 = V85)
pfactor_long_sibs_nosibs_V86 <- rbind(pfactor_long_sibs_V86, pfactor_long_nosibs)

V87 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V87==max(V87))
V87_reduced <- subset(V87, select=c("id", "V87"))

pfactor_long_sibs_V87 <- merge(pfactor_long_sibs, V87_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V87 = V86)
pfactor_long_sibs_nosibs_V87 <- rbind(pfactor_long_sibs_V87, pfactor_long_nosibs)

V88 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V88==max(V88))
V88_reduced <- subset(V88, select=c("id", "V88"))

pfactor_long_sibs_V88 <- merge(pfactor_long_sibs, V88_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V88 = V87)
pfactor_long_sibs_nosibs_V88 <- rbind(pfactor_long_sibs_V88, pfactor_long_nosibs)

V89 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V89==max(V89))
V89_reduced <- subset(V89, select=c("id", "V89"))

pfactor_long_sibs_V89 <- merge(pfactor_long_sibs, V89_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V89 = V88)
pfactor_long_sibs_nosibs_V89 <- rbind(pfactor_long_sibs_V89, pfactor_long_nosibs)

V90 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V90==max(V90))
V90_reduced <- subset(V90, select=c("id", "V90"))

pfactor_long_sibs_V90 <- merge(pfactor_long_sibs, V90_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V90 = V89)
pfactor_long_sibs_nosibs_V90 <- rbind(pfactor_long_sibs_V90, pfactor_long_nosibs)

V91 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V91==max(V91))
V91_reduced <- subset(V91, select=c("id", "V91"))

pfactor_long_sibs_V91 <- merge(pfactor_long_sibs, V91_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V91 = V90)
pfactor_long_sibs_nosibs_V91 <- rbind(pfactor_long_sibs_V91, pfactor_long_nosibs)

V92 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V92==max(V92))
V92_reduced <- subset(V92, select=c("id", "V92"))

pfactor_long_sibs_V92 <- merge(pfactor_long_sibs, V92_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V92 = V91)
pfactor_long_sibs_nosibs_V92 <- rbind(pfactor_long_sibs_V92, pfactor_long_nosibs)

V93 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V93==max(V93))
V93_reduced <- subset(V93, select=c("id", "V93"))

pfactor_long_sibs_V93 <- merge(pfactor_long_sibs, V93_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V93 = V92)
pfactor_long_sibs_nosibs_V93 <- rbind(pfactor_long_sibs_V93, pfactor_long_nosibs)

V94 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V94==max(V94))
V94_reduced <- subset(V94, select=c("id", "V94"))

pfactor_long_sibs_V94 <- merge(pfactor_long_sibs, V94_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V94 = V93)
pfactor_long_sibs_nosibs_V94 <- rbind(pfactor_long_sibs_V94, pfactor_long_nosibs)

V95 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V95==max(V95))
V95_reduced <- subset(V95, select=c("id", "V95"))

pfactor_long_sibs_V95 <- merge(pfactor_long_sibs, V95_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V95 = V94)
pfactor_long_sibs_nosibs_V95 <- rbind(pfactor_long_sibs_V95, pfactor_long_nosibs)

V96 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V96==max(V96))
V96_reduced <- subset(V96, select=c("id", "V96"))

pfactor_long_sibs_V96 <- merge(pfactor_long_sibs, V96_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V96 = V95)
pfactor_long_sibs_nosibs_V96 <- rbind(pfactor_long_sibs_V96, pfactor_long_nosibs)

V97 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V97==max(V97))
V97_reduced <- subset(V97, select=c("id", "V97"))

pfactor_long_sibs_V97 <- merge(pfactor_long_sibs, V97_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V97 = V96)
pfactor_long_sibs_nosibs_V97 <- rbind(pfactor_long_sibs_V97, pfactor_long_nosibs)

V98 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V98==max(V98))
V98_reduced <- subset(V98, select=c("id", "V98"))

pfactor_long_sibs_V98 <- merge(pfactor_long_sibs, V98_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V98 = V97)
pfactor_long_sibs_nosibs_V98 <- rbind(pfactor_long_sibs_V98, pfactor_long_nosibs)

V99 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V99==max(V99))
V99_reduced <- subset(V99, select=c("id", "V99"))

pfactor_long_sibs_V99 <- merge(pfactor_long_sibs, V99_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V99 = V98)
pfactor_long_sibs_nosibs_V99 <- rbind(pfactor_long_sibs_V99, pfactor_long_nosibs)

V100 <- pfactor_wide_sibs %>% 
  group_by(rel_family_id) %>%
  filter(V100==max(V100))
V100_reduced <- subset(V100, select=c("id", "V100"))

pfactor_long_sibs_V100 <- merge(pfactor_long_sibs, V100_reduced)
pfactor_long_nosibs <-  dplyr::rename(pfactor_long_nosibs, V100 = V99)
pfactor_long_sibs_nosibs_V100 <- rbind(pfactor_long_sibs_V100, pfactor_long_nosibs)


##Global structure analysis
#get data ready for analysis

#Group-mean center global brain structure predictors by site 
groupmeans_V1 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V1, mean)
names(groupmeans_V1) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V1$site_meanct_cent <- groupmeans_V1$site_meanct - mean(groupmeans_V1$site_meanct)
groupmeans_V1$site_totalvol_cent <- groupmeans_V1$site_totalvol - mean(groupmeans_V1$site_totalvol)
groupmeans_V1$site_totalarea_cent <- groupmeans_V1$site_totalarea - mean(groupmeans_V1$site_totalarea)
groupmeans_V1$site_subcort_vol_cent <- groupmeans_V1$site_subcort_vol - mean(groupmeans_V1$site_subcort_vol)

groupmeans_V2 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V2, mean)
names(groupmeans_V2) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V2$site_meanct_cent <- groupmeans_V2$site_meanct - mean(groupmeans_V2$site_meanct)
groupmeans_V2$site_totalvol_cent <- groupmeans_V2$site_totalvol - mean(groupmeans_V2$site_totalvol)
groupmeans_V2$site_totalarea_cent <- groupmeans_V2$site_totalarea - mean(groupmeans_V2$site_totalarea)
groupmeans_V2$site_subcort_vol_cent <- groupmeans_V2$site_subcort_vol - mean(groupmeans_V2$site_subcort_vol)

groupmeans_V3 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V3, mean)
names(groupmeans_V3) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V3$site_meanct_cent <- groupmeans_V3$site_meanct - mean(groupmeans_V3$site_meanct)
groupmeans_V3$site_totalvol_cent <- groupmeans_V3$site_totalvol - mean(groupmeans_V3$site_totalvol)
groupmeans_V3$site_totalarea_cent <- groupmeans_V3$site_totalarea - mean(groupmeans_V3$site_totalarea)
groupmeans_V3$site_subcort_vol_cent <- groupmeans_V3$site_subcort_vol - mean(groupmeans_V3$site_subcort_vol)

groupmeans_V4 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V4, mean)
names(groupmeans_V4) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V4$site_meanct_cent <- groupmeans_V4$site_meanct - mean(groupmeans_V4$site_meanct)
groupmeans_V4$site_totalvol_cent <- groupmeans_V4$site_totalvol - mean(groupmeans_V4$site_totalvol)
groupmeans_V4$site_totalarea_cent <- groupmeans_V4$site_totalarea - mean(groupmeans_V4$site_totalarea)
groupmeans_V4$site_subcort_vol_cent <- groupmeans_V4$site_subcort_vol - mean(groupmeans_V4$site_subcort_vol)

groupmeans_V5 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V5, mean)
names(groupmeans_V5) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V5$site_meanct_cent <- groupmeans_V5$site_meanct - mean(groupmeans_V5$site_meanct)
groupmeans_V5$site_totalvol_cent <- groupmeans_V5$site_totalvol - mean(groupmeans_V5$site_totalvol)
groupmeans_V5$site_totalarea_cent <- groupmeans_V5$site_totalarea - mean(groupmeans_V5$site_totalarea)
groupmeans_V5$site_subcort_vol_cent <- groupmeans_V5$site_subcort_vol - mean(groupmeans_V5$site_subcort_vol)

groupmeans_V6 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V6, mean)
names(groupmeans_V6) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V6$site_meanct_cent <- groupmeans_V6$site_meanct - mean(groupmeans_V6$site_meanct)
groupmeans_V6$site_totalvol_cent <- groupmeans_V6$site_totalvol - mean(groupmeans_V6$site_totalvol)
groupmeans_V6$site_totalarea_cent <- groupmeans_V6$site_totalarea - mean(groupmeans_V6$site_totalarea)
groupmeans_V6$site_subcort_vol_cent <- groupmeans_V6$site_subcort_vol - mean(groupmeans_V6$site_subcort_vol)

groupmeans_V7 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V7, mean)
names(groupmeans_V7) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V7$site_meanct_cent <- groupmeans_V7$site_meanct - mean(groupmeans_V7$site_meanct)
groupmeans_V7$site_totalvol_cent <- groupmeans_V7$site_totalvol - mean(groupmeans_V7$site_totalvol)
groupmeans_V7$site_totalarea_cent <- groupmeans_V7$site_totalarea - mean(groupmeans_V7$site_totalarea)
groupmeans_V7$site_subcort_vol_cent <- groupmeans_V7$site_subcort_vol - mean(groupmeans_V7$site_subcort_vol)

groupmeans_V8 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V8, mean)
names(groupmeans_V8) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V8$site_meanct_cent <- groupmeans_V8$site_meanct - mean(groupmeans_V8$site_meanct)
groupmeans_V8$site_totalvol_cent <- groupmeans_V8$site_totalvol - mean(groupmeans_V8$site_totalvol)
groupmeans_V8$site_totalarea_cent <- groupmeans_V8$site_totalarea - mean(groupmeans_V8$site_totalarea)
groupmeans_V8$site_subcort_vol_cent <- groupmeans_V8$site_subcort_vol - mean(groupmeans_V8$site_subcort_vol)

groupmeans_V9 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V9, mean)
names(groupmeans_V9) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V9$site_meanct_cent <- groupmeans_V9$site_meanct - mean(groupmeans_V9$site_meanct)
groupmeans_V9$site_totalvol_cent <- groupmeans_V9$site_totalvol - mean(groupmeans_V9$site_totalvol)
groupmeans_V9$site_totalarea_cent <- groupmeans_V9$site_totalarea - mean(groupmeans_V9$site_totalarea)
groupmeans_V9$site_subcort_vol_cent <- groupmeans_V9$site_subcort_vol - mean(groupmeans_V9$site_subcort_vol)

groupmeans_V10 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V10, mean)
names(groupmeans_V10) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V10$site_meanct_cent <- groupmeans_V10$site_meanct - mean(groupmeans_V10$site_meanct)
groupmeans_V10$site_totalvol_cent <- groupmeans_V10$site_totalvol - mean(groupmeans_V10$site_totalvol)
groupmeans_V10$site_totalarea_cent <- groupmeans_V10$site_totalarea - mean(groupmeans_V10$site_totalarea)
groupmeans_V10$site_subcort_vol_cent <- groupmeans_V10$site_subcort_vol - mean(groupmeans_V10$site_subcort_vol)

groupmeans_V11 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V11, mean)
names(groupmeans_V11) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V11$site_meanct_cent <- groupmeans_V11$site_meanct - mean(groupmeans_V11$site_meanct)
groupmeans_V11$site_totalvol_cent <- groupmeans_V11$site_totalvol - mean(groupmeans_V11$site_totalvol)
groupmeans_V11$site_totalarea_cent <- groupmeans_V11$site_totalarea - mean(groupmeans_V11$site_totalarea)
groupmeans_V11$site_subcort_vol_cent <- groupmeans_V11$site_subcort_vol - mean(groupmeans_V11$site_subcort_vol)

groupmeans_V12 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V12, mean)
names(groupmeans_V12) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V12$site_meanct_cent <- groupmeans_V12$site_meanct - mean(groupmeans_V12$site_meanct)
groupmeans_V12$site_totalvol_cent <- groupmeans_V12$site_totalvol - mean(groupmeans_V12$site_totalvol)
groupmeans_V12$site_totalarea_cent <- groupmeans_V12$site_totalarea - mean(groupmeans_V12$site_totalarea)
groupmeans_V12$site_subcort_vol_cent <- groupmeans_V12$site_subcort_vol - mean(groupmeans_V12$site_subcort_vol)

groupmeans_V13 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V13, mean)
names(groupmeans_V13) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V13$site_meanct_cent <- groupmeans_V13$site_meanct - mean(groupmeans_V13$site_meanct)
groupmeans_V13$site_totalvol_cent <- groupmeans_V13$site_totalvol - mean(groupmeans_V13$site_totalvol)
groupmeans_V13$site_totalarea_cent <- groupmeans_V13$site_totalarea - mean(groupmeans_V13$site_totalarea)
groupmeans_V13$site_subcort_vol_cent <- groupmeans_V13$site_subcort_vol - mean(groupmeans_V13$site_subcort_vol)

groupmeans_V14 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V14, mean)
names(groupmeans_V14) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V14$site_meanct_cent <- groupmeans_V14$site_meanct - mean(groupmeans_V14$site_meanct)
groupmeans_V14$site_totalvol_cent <- groupmeans_V14$site_totalvol - mean(groupmeans_V14$site_totalvol)
groupmeans_V14$site_totalarea_cent <- groupmeans_V14$site_totalarea - mean(groupmeans_V14$site_totalarea)
groupmeans_V14$site_subcort_vol_cent <- groupmeans_V14$site_subcort_vol - mean(groupmeans_V14$site_subcort_vol)

groupmeans_V15 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V15, mean)
names(groupmeans_V15) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V15$site_meanct_cent <- groupmeans_V15$site_meanct - mean(groupmeans_V15$site_meanct)
groupmeans_V15$site_totalvol_cent <- groupmeans_V15$site_totalvol - mean(groupmeans_V15$site_totalvol)
groupmeans_V15$site_totalarea_cent <- groupmeans_V15$site_totalarea - mean(groupmeans_V15$site_totalarea)
groupmeans_V15$site_subcort_vol_cent <- groupmeans_V15$site_subcort_vol - mean(groupmeans_V15$site_subcort_vol)

groupmeans_V16 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V16, mean)
names(groupmeans_V16) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V16$site_meanct_cent <- groupmeans_V16$site_meanct - mean(groupmeans_V16$site_meanct)
groupmeans_V16$site_totalvol_cent <- groupmeans_V16$site_totalvol - mean(groupmeans_V16$site_totalvol)
groupmeans_V16$site_totalarea_cent <- groupmeans_V16$site_totalarea - mean(groupmeans_V16$site_totalarea)
groupmeans_V16$site_subcort_vol_cent <- groupmeans_V16$site_subcort_vol - mean(groupmeans_V16$site_subcort_vol)

groupmeans_V17 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V17, mean)
names(groupmeans_V17) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V17$site_meanct_cent <- groupmeans_V17$site_meanct - mean(groupmeans_V17$site_meanct)
groupmeans_V17$site_totalvol_cent <- groupmeans_V17$site_totalvol - mean(groupmeans_V17$site_totalvol)
groupmeans_V17$site_totalarea_cent <- groupmeans_V17$site_totalarea - mean(groupmeans_V17$site_totalarea)
groupmeans_V17$site_subcort_vol_cent <- groupmeans_V17$site_subcort_vol - mean(groupmeans_V17$site_subcort_vol)

groupmeans_V18 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V18, mean)
names(groupmeans_V18) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V18$site_meanct_cent <- groupmeans_V18$site_meanct - mean(groupmeans_V18$site_meanct)
groupmeans_V18$site_totalvol_cent <- groupmeans_V18$site_totalvol - mean(groupmeans_V18$site_totalvol)
groupmeans_V18$site_totalarea_cent <- groupmeans_V18$site_totalarea - mean(groupmeans_V18$site_totalarea)
groupmeans_V18$site_subcort_vol_cent <- groupmeans_V18$site_subcort_vol - mean(groupmeans_V18$site_subcort_vol)

groupmeans_V19 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V19, mean)
names(groupmeans_V19) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V19$site_meanct_cent <- groupmeans_V19$site_meanct - mean(groupmeans_V19$site_meanct)
groupmeans_V19$site_totalvol_cent <- groupmeans_V19$site_totalvol - mean(groupmeans_V19$site_totalvol)
groupmeans_V19$site_totalarea_cent <- groupmeans_V19$site_totalarea - mean(groupmeans_V19$site_totalarea)
groupmeans_V19$site_subcort_vol_cent <- groupmeans_V19$site_subcort_vol - mean(groupmeans_V19$site_subcort_vol)

groupmeans_V20 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V20, mean)
names(groupmeans_V20) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V20$site_meanct_cent <- groupmeans_V20$site_meanct - mean(groupmeans_V20$site_meanct)
groupmeans_V20$site_totalvol_cent <- groupmeans_V20$site_totalvol - mean(groupmeans_V20$site_totalvol)
groupmeans_V20$site_totalarea_cent <- groupmeans_V20$site_totalarea - mean(groupmeans_V20$site_totalarea)
groupmeans_V20$site_subcort_vol_cent <- groupmeans_V20$site_subcort_vol - mean(groupmeans_V20$site_subcort_vol)

groupmeans_V21 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V21, mean)
names(groupmeans_V21) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V21$site_meanct_cent <- groupmeans_V21$site_meanct - mean(groupmeans_V21$site_meanct)
groupmeans_V21$site_totalvol_cent <- groupmeans_V21$site_totalvol - mean(groupmeans_V21$site_totalvol)
groupmeans_V21$site_totalarea_cent <- groupmeans_V21$site_totalarea - mean(groupmeans_V21$site_totalarea)
groupmeans_V21$site_subcort_vol_cent <- groupmeans_V21$site_subcort_vol - mean(groupmeans_V21$site_subcort_vol)

groupmeans_V22 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V22, mean)
names(groupmeans_V22) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V22$site_meanct_cent <- groupmeans_V22$site_meanct - mean(groupmeans_V22$site_meanct)
groupmeans_V22$site_totalvol_cent <- groupmeans_V22$site_totalvol - mean(groupmeans_V22$site_totalvol)
groupmeans_V22$site_totalarea_cent <- groupmeans_V22$site_totalarea - mean(groupmeans_V22$site_totalarea)
groupmeans_V22$site_subcort_vol_cent <- groupmeans_V22$site_subcort_vol - mean(groupmeans_V22$site_subcort_vol)

groupmeans_V23 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V23, mean)
names(groupmeans_V23) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V23$site_meanct_cent <- groupmeans_V23$site_meanct - mean(groupmeans_V23$site_meanct)
groupmeans_V23$site_totalvol_cent <- groupmeans_V23$site_totalvol - mean(groupmeans_V23$site_totalvol)
groupmeans_V23$site_totalarea_cent <- groupmeans_V23$site_totalarea - mean(groupmeans_V23$site_totalarea)
groupmeans_V23$site_subcort_vol_cent <- groupmeans_V23$site_subcort_vol - mean(groupmeans_V23$site_subcort_vol)

groupmeans_V24 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V24, mean)
names(groupmeans_V24) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V24$site_meanct_cent <- groupmeans_V24$site_meanct - mean(groupmeans_V24$site_meanct)
groupmeans_V24$site_totalvol_cent <- groupmeans_V24$site_totalvol - mean(groupmeans_V24$site_totalvol)
groupmeans_V24$site_totalarea_cent <- groupmeans_V24$site_totalarea - mean(groupmeans_V24$site_totalarea)
groupmeans_V24$site_subcort_vol_cent <- groupmeans_V24$site_subcort_vol - mean(groupmeans_V24$site_subcort_vol)

groupmeans_V25 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V25, mean)
names(groupmeans_V25) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V25$site_meanct_cent <- groupmeans_V25$site_meanct - mean(groupmeans_V25$site_meanct)
groupmeans_V25$site_totalvol_cent <- groupmeans_V25$site_totalvol - mean(groupmeans_V25$site_totalvol)
groupmeans_V25$site_totalarea_cent <- groupmeans_V25$site_totalarea - mean(groupmeans_V25$site_totalarea)
groupmeans_V25$site_subcort_vol_cent <- groupmeans_V25$site_subcort_vol - mean(groupmeans_V25$site_subcort_vol)

groupmeans_V26 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V26, mean)
names(groupmeans_V26) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V26$site_meanct_cent <- groupmeans_V26$site_meanct - mean(groupmeans_V26$site_meanct)
groupmeans_V26$site_totalvol_cent <- groupmeans_V26$site_totalvol - mean(groupmeans_V26$site_totalvol)
groupmeans_V26$site_totalarea_cent <- groupmeans_V26$site_totalarea - mean(groupmeans_V26$site_totalarea)
groupmeans_V26$site_subcort_vol_cent <- groupmeans_V26$site_subcort_vol - mean(groupmeans_V26$site_subcort_vol)

groupmeans_V27 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V27, mean)
names(groupmeans_V27) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V27$site_meanct_cent <- groupmeans_V27$site_meanct - mean(groupmeans_V27$site_meanct)
groupmeans_V27$site_totalvol_cent <- groupmeans_V27$site_totalvol - mean(groupmeans_V27$site_totalvol)
groupmeans_V27$site_totalarea_cent <- groupmeans_V27$site_totalarea - mean(groupmeans_V27$site_totalarea)
groupmeans_V27$site_subcort_vol_cent <- groupmeans_V27$site_subcort_vol - mean(groupmeans_V27$site_subcort_vol)

groupmeans_V28 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V28, mean)
names(groupmeans_V28) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V28$site_meanct_cent <- groupmeans_V28$site_meanct - mean(groupmeans_V28$site_meanct)
groupmeans_V28$site_totalvol_cent <- groupmeans_V28$site_totalvol - mean(groupmeans_V28$site_totalvol)
groupmeans_V28$site_totalarea_cent <- groupmeans_V28$site_totalarea - mean(groupmeans_V28$site_totalarea)
groupmeans_V28$site_subcort_vol_cent <- groupmeans_V28$site_subcort_vol - mean(groupmeans_V28$site_subcort_vol)

groupmeans_V29 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V29, mean)
names(groupmeans_V29) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V29$site_meanct_cent <- groupmeans_V29$site_meanct - mean(groupmeans_V29$site_meanct)
groupmeans_V29$site_totalvol_cent <- groupmeans_V29$site_totalvol - mean(groupmeans_V29$site_totalvol)
groupmeans_V29$site_totalarea_cent <- groupmeans_V29$site_totalarea - mean(groupmeans_V29$site_totalarea)
groupmeans_V29$site_subcort_vol_cent <- groupmeans_V29$site_subcort_vol - mean(groupmeans_V29$site_subcort_vol)

groupmeans_V30 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V30, mean)
names(groupmeans_V30) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V30$site_meanct_cent <- groupmeans_V30$site_meanct - mean(groupmeans_V30$site_meanct)
groupmeans_V30$site_totalvol_cent <- groupmeans_V30$site_totalvol - mean(groupmeans_V30$site_totalvol)
groupmeans_V30$site_totalarea_cent <- groupmeans_V30$site_totalarea - mean(groupmeans_V30$site_totalarea)
groupmeans_V30$site_subcort_vol_cent <- groupmeans_V30$site_subcort_vol - mean(groupmeans_V30$site_subcort_vol)

groupmeans_V31 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V31, mean)
names(groupmeans_V31) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V31$site_meanct_cent <- groupmeans_V31$site_meanct - mean(groupmeans_V31$site_meanct)
groupmeans_V31$site_totalvol_cent <- groupmeans_V31$site_totalvol - mean(groupmeans_V31$site_totalvol)
groupmeans_V31$site_totalarea_cent <- groupmeans_V31$site_totalarea - mean(groupmeans_V31$site_totalarea)
groupmeans_V31$site_subcort_vol_cent <- groupmeans_V31$site_subcort_vol - mean(groupmeans_V31$site_subcort_vol)

groupmeans_V32 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V32, mean)
names(groupmeans_V32) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V32$site_meanct_cent <- groupmeans_V32$site_meanct - mean(groupmeans_V32$site_meanct)
groupmeans_V32$site_totalvol_cent <- groupmeans_V32$site_totalvol - mean(groupmeans_V32$site_totalvol)
groupmeans_V32$site_totalarea_cent <- groupmeans_V32$site_totalarea - mean(groupmeans_V32$site_totalarea)
groupmeans_V32$site_subcort_vol_cent <- groupmeans_V32$site_subcort_vol - mean(groupmeans_V32$site_subcort_vol)

groupmeans_V33 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V33, mean)
names(groupmeans_V33) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V33$site_meanct_cent <- groupmeans_V33$site_meanct - mean(groupmeans_V33$site_meanct)
groupmeans_V33$site_totalvol_cent <- groupmeans_V33$site_totalvol - mean(groupmeans_V33$site_totalvol)
groupmeans_V33$site_totalarea_cent <- groupmeans_V33$site_totalarea - mean(groupmeans_V33$site_totalarea)
groupmeans_V33$site_subcort_vol_cent <- groupmeans_V33$site_subcort_vol - mean(groupmeans_V33$site_subcort_vol)

groupmeans_V34 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V34, mean)
names(groupmeans_V34) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V34$site_meanct_cent <- groupmeans_V34$site_meanct - mean(groupmeans_V34$site_meanct)
groupmeans_V34$site_totalvol_cent <- groupmeans_V34$site_totalvol - mean(groupmeans_V34$site_totalvol)
groupmeans_V34$site_totalarea_cent <- groupmeans_V34$site_totalarea - mean(groupmeans_V34$site_totalarea)
groupmeans_V34$site_subcort_vol_cent <- groupmeans_V34$site_subcort_vol - mean(groupmeans_V34$site_subcort_vol)

groupmeans_V35 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V35, mean)
names(groupmeans_V35) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V35$site_meanct_cent <- groupmeans_V35$site_meanct - mean(groupmeans_V35$site_meanct)
groupmeans_V35$site_totalvol_cent <- groupmeans_V35$site_totalvol - mean(groupmeans_V35$site_totalvol)
groupmeans_V35$site_totalarea_cent <- groupmeans_V35$site_totalarea - mean(groupmeans_V35$site_totalarea)
groupmeans_V35$site_subcort_vol_cent <- groupmeans_V35$site_subcort_vol - mean(groupmeans_V35$site_subcort_vol)

groupmeans_V36 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V36, mean)
names(groupmeans_V36) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V36$site_meanct_cent <- groupmeans_V36$site_meanct - mean(groupmeans_V36$site_meanct)
groupmeans_V36$site_totalvol_cent <- groupmeans_V36$site_totalvol - mean(groupmeans_V36$site_totalvol)
groupmeans_V36$site_totalarea_cent <- groupmeans_V36$site_totalarea - mean(groupmeans_V36$site_totalarea)
groupmeans_V36$site_subcort_vol_cent <- groupmeans_V36$site_subcort_vol - mean(groupmeans_V36$site_subcort_vol)

groupmeans_V37 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V37, mean)
names(groupmeans_V37) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V37$site_meanct_cent <- groupmeans_V37$site_meanct - mean(groupmeans_V37$site_meanct)
groupmeans_V37$site_totalvol_cent <- groupmeans_V37$site_totalvol - mean(groupmeans_V37$site_totalvol)
groupmeans_V37$site_totalarea_cent <- groupmeans_V37$site_totalarea - mean(groupmeans_V37$site_totalarea)
groupmeans_V37$site_subcort_vol_cent <- groupmeans_V37$site_subcort_vol - mean(groupmeans_V37$site_subcort_vol)

groupmeans_V38 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V38, mean)
names(groupmeans_V38) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V38$site_meanct_cent <- groupmeans_V38$site_meanct - mean(groupmeans_V38$site_meanct)
groupmeans_V38$site_totalvol_cent <- groupmeans_V38$site_totalvol - mean(groupmeans_V38$site_totalvol)
groupmeans_V38$site_totalarea_cent <- groupmeans_V38$site_totalarea - mean(groupmeans_V38$site_totalarea)
groupmeans_V38$site_subcort_vol_cent <- groupmeans_V38$site_subcort_vol - mean(groupmeans_V38$site_subcort_vol)

groupmeans_V39 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V39, mean)
names(groupmeans_V39) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V39$site_meanct_cent <- groupmeans_V39$site_meanct - mean(groupmeans_V39$site_meanct)
groupmeans_V39$site_totalvol_cent <- groupmeans_V39$site_totalvol - mean(groupmeans_V39$site_totalvol)
groupmeans_V39$site_totalarea_cent <- groupmeans_V39$site_totalarea - mean(groupmeans_V39$site_totalarea)
groupmeans_V39$site_subcort_vol_cent <- groupmeans_V39$site_subcort_vol - mean(groupmeans_V39$site_subcort_vol)

groupmeans_V40 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V40, mean)
names(groupmeans_V40) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V40$site_meanct_cent <- groupmeans_V40$site_meanct - mean(groupmeans_V40$site_meanct)
groupmeans_V40$site_totalvol_cent <- groupmeans_V40$site_totalvol - mean(groupmeans_V40$site_totalvol)
groupmeans_V40$site_totalarea_cent <- groupmeans_V40$site_totalarea - mean(groupmeans_V40$site_totalarea)
groupmeans_V40$site_subcort_vol_cent <- groupmeans_V40$site_subcort_vol - mean(groupmeans_V40$site_subcort_vol)

groupmeans_V41 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V41, mean)
names(groupmeans_V41) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V41$site_meanct_cent <- groupmeans_V41$site_meanct - mean(groupmeans_V41$site_meanct)
groupmeans_V41$site_totalvol_cent <- groupmeans_V41$site_totalvol - mean(groupmeans_V41$site_totalvol)
groupmeans_V41$site_totalarea_cent <- groupmeans_V41$site_totalarea - mean(groupmeans_V41$site_totalarea)
groupmeans_V41$site_subcort_vol_cent <- groupmeans_V41$site_subcort_vol - mean(groupmeans_V41$site_subcort_vol)

groupmeans_V42 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V42, mean)
names(groupmeans_V42) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V42$site_meanct_cent <- groupmeans_V42$site_meanct - mean(groupmeans_V42$site_meanct)
groupmeans_V42$site_totalvol_cent <- groupmeans_V42$site_totalvol - mean(groupmeans_V42$site_totalvol)
groupmeans_V42$site_totalarea_cent <- groupmeans_V42$site_totalarea - mean(groupmeans_V42$site_totalarea)
groupmeans_V42$site_subcort_vol_cent <- groupmeans_V42$site_subcort_vol - mean(groupmeans_V42$site_subcort_vol)

groupmeans_V43 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V43, mean)
names(groupmeans_V43) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V43$site_meanct_cent <- groupmeans_V43$site_meanct - mean(groupmeans_V43$site_meanct)
groupmeans_V43$site_totalvol_cent <- groupmeans_V43$site_totalvol - mean(groupmeans_V43$site_totalvol)
groupmeans_V43$site_totalarea_cent <- groupmeans_V43$site_totalarea - mean(groupmeans_V43$site_totalarea)
groupmeans_V43$site_subcort_vol_cent <- groupmeans_V43$site_subcort_vol - mean(groupmeans_V43$site_subcort_vol)

groupmeans_V44 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V44, mean)
names(groupmeans_V44) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V44$site_meanct_cent <- groupmeans_V44$site_meanct - mean(groupmeans_V44$site_meanct)
groupmeans_V44$site_totalvol_cent <- groupmeans_V44$site_totalvol - mean(groupmeans_V44$site_totalvol)
groupmeans_V44$site_totalarea_cent <- groupmeans_V44$site_totalarea - mean(groupmeans_V44$site_totalarea)
groupmeans_V44$site_subcort_vol_cent <- groupmeans_V44$site_subcort_vol - mean(groupmeans_V44$site_subcort_vol)

groupmeans_V45 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V45, mean)
names(groupmeans_V45) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V45$site_meanct_cent <- groupmeans_V45$site_meanct - mean(groupmeans_V45$site_meanct)
groupmeans_V45$site_totalvol_cent <- groupmeans_V45$site_totalvol - mean(groupmeans_V45$site_totalvol)
groupmeans_V45$site_totalarea_cent <- groupmeans_V45$site_totalarea - mean(groupmeans_V45$site_totalarea)
groupmeans_V45$site_subcort_vol_cent <- groupmeans_V45$site_subcort_vol - mean(groupmeans_V45$site_subcort_vol)

groupmeans_V46 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V46, mean)
names(groupmeans_V46) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V46$site_meanct_cent <- groupmeans_V46$site_meanct - mean(groupmeans_V46$site_meanct)
groupmeans_V46$site_totalvol_cent <- groupmeans_V46$site_totalvol - mean(groupmeans_V46$site_totalvol)
groupmeans_V46$site_totalarea_cent <- groupmeans_V46$site_totalarea - mean(groupmeans_V46$site_totalarea)
groupmeans_V46$site_subcort_vol_cent <- groupmeans_V46$site_subcort_vol - mean(groupmeans_V46$site_subcort_vol)

groupmeans_V47 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V47, mean)
names(groupmeans_V47) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V47$site_meanct_cent <- groupmeans_V47$site_meanct - mean(groupmeans_V47$site_meanct)
groupmeans_V47$site_totalvol_cent <- groupmeans_V47$site_totalvol - mean(groupmeans_V47$site_totalvol)
groupmeans_V47$site_totalarea_cent <- groupmeans_V47$site_totalarea - mean(groupmeans_V47$site_totalarea)
groupmeans_V47$site_subcort_vol_cent <- groupmeans_V47$site_subcort_vol - mean(groupmeans_V47$site_subcort_vol)

groupmeans_V48 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V48, mean)
names(groupmeans_V48) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V48$site_meanct_cent <- groupmeans_V48$site_meanct - mean(groupmeans_V48$site_meanct)
groupmeans_V48$site_totalvol_cent <- groupmeans_V48$site_totalvol - mean(groupmeans_V48$site_totalvol)
groupmeans_V48$site_totalarea_cent <- groupmeans_V48$site_totalarea - mean(groupmeans_V48$site_totalarea)
groupmeans_V48$site_subcort_vol_cent <- groupmeans_V48$site_subcort_vol - mean(groupmeans_V48$site_subcort_vol)

groupmeans_V49 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V49, mean)
names(groupmeans_V49) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V49$site_meanct_cent <- groupmeans_V49$site_meanct - mean(groupmeans_V49$site_meanct)
groupmeans_V49$site_totalvol_cent <- groupmeans_V49$site_totalvol - mean(groupmeans_V49$site_totalvol)
groupmeans_V49$site_totalarea_cent <- groupmeans_V49$site_totalarea - mean(groupmeans_V49$site_totalarea)
groupmeans_V49$site_subcort_vol_cent <- groupmeans_V49$site_subcort_vol - mean(groupmeans_V49$site_subcort_vol)

groupmeans_V50 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V50, mean)
names(groupmeans_V50) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V50$site_meanct_cent <- groupmeans_V50$site_meanct - mean(groupmeans_V50$site_meanct)
groupmeans_V50$site_totalvol_cent <- groupmeans_V50$site_totalvol - mean(groupmeans_V50$site_totalvol)
groupmeans_V50$site_totalarea_cent <- groupmeans_V50$site_totalarea - mean(groupmeans_V50$site_totalarea)
groupmeans_V50$site_subcort_vol_cent <- groupmeans_V50$site_subcort_vol - mean(groupmeans_V50$site_subcort_vol)

groupmeans_V51 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V51, mean)
names(groupmeans_V51) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V51$site_meanct_cent <- groupmeans_V51$site_meanct - mean(groupmeans_V51$site_meanct)
groupmeans_V51$site_totalvol_cent <- groupmeans_V51$site_totalvol - mean(groupmeans_V51$site_totalvol)
groupmeans_V51$site_totalarea_cent <- groupmeans_V51$site_totalarea - mean(groupmeans_V51$site_totalarea)
groupmeans_V51$site_subcort_vol_cent <- groupmeans_V51$site_subcort_vol - mean(groupmeans_V51$site_subcort_vol)

groupmeans_V52 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V52, mean)
names(groupmeans_V52) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V52$site_meanct_cent <- groupmeans_V52$site_meanct - mean(groupmeans_V52$site_meanct)
groupmeans_V52$site_totalvol_cent <- groupmeans_V52$site_totalvol - mean(groupmeans_V52$site_totalvol)
groupmeans_V52$site_totalarea_cent <- groupmeans_V52$site_totalarea - mean(groupmeans_V52$site_totalarea)
groupmeans_V52$site_subcort_vol_cent <- groupmeans_V52$site_subcort_vol - mean(groupmeans_V52$site_subcort_vol)

groupmeans_V53 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V53, mean)
names(groupmeans_V53) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V53$site_meanct_cent <- groupmeans_V53$site_meanct - mean(groupmeans_V53$site_meanct)
groupmeans_V53$site_totalvol_cent <- groupmeans_V53$site_totalvol - mean(groupmeans_V53$site_totalvol)
groupmeans_V53$site_totalarea_cent <- groupmeans_V53$site_totalarea - mean(groupmeans_V53$site_totalarea)
groupmeans_V53$site_subcort_vol_cent <- groupmeans_V53$site_subcort_vol - mean(groupmeans_V53$site_subcort_vol)

groupmeans_V54 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V54, mean)
names(groupmeans_V54) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V54$site_meanct_cent <- groupmeans_V54$site_meanct - mean(groupmeans_V54$site_meanct)
groupmeans_V54$site_totalvol_cent <- groupmeans_V54$site_totalvol - mean(groupmeans_V54$site_totalvol)
groupmeans_V54$site_totalarea_cent <- groupmeans_V54$site_totalarea - mean(groupmeans_V54$site_totalarea)
groupmeans_V54$site_subcort_vol_cent <- groupmeans_V54$site_subcort_vol - mean(groupmeans_V54$site_subcort_vol)

groupmeans_V55 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V55, mean)
names(groupmeans_V55) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V55$site_meanct_cent <- groupmeans_V55$site_meanct - mean(groupmeans_V55$site_meanct)
groupmeans_V55$site_totalvol_cent <- groupmeans_V55$site_totalvol - mean(groupmeans_V55$site_totalvol)
groupmeans_V55$site_totalarea_cent <- groupmeans_V55$site_totalarea - mean(groupmeans_V55$site_totalarea)
groupmeans_V55$site_subcort_vol_cent <- groupmeans_V55$site_subcort_vol - mean(groupmeans_V55$site_subcort_vol)

groupmeans_V56 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V56, mean)
names(groupmeans_V56) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V56$site_meanct_cent <- groupmeans_V56$site_meanct - mean(groupmeans_V56$site_meanct)
groupmeans_V56$site_totalvol_cent <- groupmeans_V56$site_totalvol - mean(groupmeans_V56$site_totalvol)
groupmeans_V56$site_totalarea_cent <- groupmeans_V56$site_totalarea - mean(groupmeans_V56$site_totalarea)
groupmeans_V56$site_subcort_vol_cent <- groupmeans_V56$site_subcort_vol - mean(groupmeans_V56$site_subcort_vol)

groupmeans_V57 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V57, mean)
names(groupmeans_V57) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V57$site_meanct_cent <- groupmeans_V57$site_meanct - mean(groupmeans_V57$site_meanct)
groupmeans_V57$site_totalvol_cent <- groupmeans_V57$site_totalvol - mean(groupmeans_V57$site_totalvol)
groupmeans_V57$site_totalarea_cent <- groupmeans_V57$site_totalarea - mean(groupmeans_V57$site_totalarea)
groupmeans_V57$site_subcort_vol_cent <- groupmeans_V57$site_subcort_vol - mean(groupmeans_V57$site_subcort_vol)

groupmeans_V58 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V58, mean)
names(groupmeans_V58) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V58$site_meanct_cent <- groupmeans_V58$site_meanct - mean(groupmeans_V58$site_meanct)
groupmeans_V58$site_totalvol_cent <- groupmeans_V58$site_totalvol - mean(groupmeans_V58$site_totalvol)
groupmeans_V58$site_totalarea_cent <- groupmeans_V58$site_totalarea - mean(groupmeans_V58$site_totalarea)
groupmeans_V58$site_subcort_vol_cent <- groupmeans_V58$site_subcort_vol - mean(groupmeans_V58$site_subcort_vol)

groupmeans_V59 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V59, mean)
names(groupmeans_V59) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V59$site_meanct_cent <- groupmeans_V59$site_meanct - mean(groupmeans_V59$site_meanct)
groupmeans_V59$site_totalvol_cent <- groupmeans_V59$site_totalvol - mean(groupmeans_V59$site_totalvol)
groupmeans_V59$site_totalarea_cent <- groupmeans_V59$site_totalarea - mean(groupmeans_V59$site_totalarea)
groupmeans_V59$site_subcort_vol_cent <- groupmeans_V59$site_subcort_vol - mean(groupmeans_V59$site_subcort_vol)

groupmeans_V60 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V60, mean)
names(groupmeans_V60) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V60$site_meanct_cent <- groupmeans_V60$site_meanct - mean(groupmeans_V60$site_meanct)
groupmeans_V60$site_totalvol_cent <- groupmeans_V60$site_totalvol - mean(groupmeans_V60$site_totalvol)
groupmeans_V60$site_totalarea_cent <- groupmeans_V60$site_totalarea - mean(groupmeans_V60$site_totalarea)
groupmeans_V60$site_subcort_vol_cent <- groupmeans_V60$site_subcort_vol - mean(groupmeans_V60$site_subcort_vol)

groupmeans_V61 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V61, mean)
names(groupmeans_V61) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V61$site_meanct_cent <- groupmeans_V61$site_meanct - mean(groupmeans_V61$site_meanct)
groupmeans_V61$site_totalvol_cent <- groupmeans_V61$site_totalvol - mean(groupmeans_V61$site_totalvol)
groupmeans_V61$site_totalarea_cent <- groupmeans_V61$site_totalarea - mean(groupmeans_V61$site_totalarea)
groupmeans_V61$site_subcort_vol_cent <- groupmeans_V61$site_subcort_vol - mean(groupmeans_V61$site_subcort_vol)

groupmeans_V62 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V62, mean)
names(groupmeans_V62) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V62$site_meanct_cent <- groupmeans_V62$site_meanct - mean(groupmeans_V62$site_meanct)
groupmeans_V62$site_totalvol_cent <- groupmeans_V62$site_totalvol - mean(groupmeans_V62$site_totalvol)
groupmeans_V62$site_totalarea_cent <- groupmeans_V62$site_totalarea - mean(groupmeans_V62$site_totalarea)
groupmeans_V62$site_subcort_vol_cent <- groupmeans_V62$site_subcort_vol - mean(groupmeans_V62$site_subcort_vol)

groupmeans_V63 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V63, mean)
names(groupmeans_V63) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V63$site_meanct_cent <- groupmeans_V63$site_meanct - mean(groupmeans_V63$site_meanct)
groupmeans_V63$site_totalvol_cent <- groupmeans_V63$site_totalvol - mean(groupmeans_V63$site_totalvol)
groupmeans_V63$site_totalarea_cent <- groupmeans_V63$site_totalarea - mean(groupmeans_V63$site_totalarea)
groupmeans_V63$site_subcort_vol_cent <- groupmeans_V63$site_subcort_vol - mean(groupmeans_V63$site_subcort_vol)

groupmeans_V64 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V64, mean)
names(groupmeans_V64) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V64$site_meanct_cent <- groupmeans_V64$site_meanct - mean(groupmeans_V64$site_meanct)
groupmeans_V64$site_totalvol_cent <- groupmeans_V64$site_totalvol - mean(groupmeans_V64$site_totalvol)
groupmeans_V64$site_totalarea_cent <- groupmeans_V64$site_totalarea - mean(groupmeans_V64$site_totalarea)
groupmeans_V64$site_subcort_vol_cent <- groupmeans_V64$site_subcort_vol - mean(groupmeans_V64$site_subcort_vol)

groupmeans_V65 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V65, mean)
names(groupmeans_V65) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V65$site_meanct_cent <- groupmeans_V65$site_meanct - mean(groupmeans_V65$site_meanct)
groupmeans_V65$site_totalvol_cent <- groupmeans_V65$site_totalvol - mean(groupmeans_V65$site_totalvol)
groupmeans_V65$site_totalarea_cent <- groupmeans_V65$site_totalarea - mean(groupmeans_V65$site_totalarea)
groupmeans_V65$site_subcort_vol_cent <- groupmeans_V65$site_subcort_vol - mean(groupmeans_V65$site_subcort_vol)

groupmeans_V66 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V66, mean)
names(groupmeans_V66) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V66$site_meanct_cent <- groupmeans_V66$site_meanct - mean(groupmeans_V66$site_meanct)
groupmeans_V66$site_totalvol_cent <- groupmeans_V66$site_totalvol - mean(groupmeans_V66$site_totalvol)
groupmeans_V66$site_totalarea_cent <- groupmeans_V66$site_totalarea - mean(groupmeans_V66$site_totalarea)
groupmeans_V66$site_subcort_vol_cent <- groupmeans_V66$site_subcort_vol - mean(groupmeans_V66$site_subcort_vol)

groupmeans_V67 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V67, mean)
names(groupmeans_V67) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V67$site_meanct_cent <- groupmeans_V67$site_meanct - mean(groupmeans_V67$site_meanct)
groupmeans_V67$site_totalvol_cent <- groupmeans_V67$site_totalvol - mean(groupmeans_V67$site_totalvol)
groupmeans_V67$site_totalarea_cent <- groupmeans_V67$site_totalarea - mean(groupmeans_V67$site_totalarea)
groupmeans_V67$site_subcort_vol_cent <- groupmeans_V67$site_subcort_vol - mean(groupmeans_V67$site_subcort_vol)

groupmeans_V68 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V68, mean)
names(groupmeans_V68) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V68$site_meanct_cent <- groupmeans_V68$site_meanct - mean(groupmeans_V68$site_meanct)
groupmeans_V68$site_totalvol_cent <- groupmeans_V68$site_totalvol - mean(groupmeans_V68$site_totalvol)
groupmeans_V68$site_totalarea_cent <- groupmeans_V68$site_totalarea - mean(groupmeans_V68$site_totalarea)
groupmeans_V68$site_subcort_vol_cent <- groupmeans_V68$site_subcort_vol - mean(groupmeans_V68$site_subcort_vol)

groupmeans_V69 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V69, mean)
names(groupmeans_V69) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V69$site_meanct_cent <- groupmeans_V69$site_meanct - mean(groupmeans_V69$site_meanct)
groupmeans_V69$site_totalvol_cent <- groupmeans_V69$site_totalvol - mean(groupmeans_V69$site_totalvol)
groupmeans_V69$site_totalarea_cent <- groupmeans_V69$site_totalarea - mean(groupmeans_V69$site_totalarea)
groupmeans_V69$site_subcort_vol_cent <- groupmeans_V69$site_subcort_vol - mean(groupmeans_V69$site_subcort_vol)

groupmeans_V70 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V70, mean)
names(groupmeans_V70) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V70$site_meanct_cent <- groupmeans_V70$site_meanct - mean(groupmeans_V70$site_meanct)
groupmeans_V70$site_totalvol_cent <- groupmeans_V70$site_totalvol - mean(groupmeans_V70$site_totalvol)
groupmeans_V70$site_totalarea_cent <- groupmeans_V70$site_totalarea - mean(groupmeans_V70$site_totalarea)
groupmeans_V70$site_subcort_vol_cent <- groupmeans_V70$site_subcort_vol - mean(groupmeans_V70$site_subcort_vol)

groupmeans_V71 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V71, mean)
names(groupmeans_V71) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V71$site_meanct_cent <- groupmeans_V71$site_meanct - mean(groupmeans_V71$site_meanct)
groupmeans_V71$site_totalvol_cent <- groupmeans_V71$site_totalvol - mean(groupmeans_V71$site_totalvol)
groupmeans_V71$site_totalarea_cent <- groupmeans_V71$site_totalarea - mean(groupmeans_V71$site_totalarea)
groupmeans_V71$site_subcort_vol_cent <- groupmeans_V71$site_subcort_vol - mean(groupmeans_V71$site_subcort_vol)

groupmeans_V72 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V72, mean)
names(groupmeans_V72) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V72$site_meanct_cent <- groupmeans_V72$site_meanct - mean(groupmeans_V72$site_meanct)
groupmeans_V72$site_totalvol_cent <- groupmeans_V72$site_totalvol - mean(groupmeans_V72$site_totalvol)
groupmeans_V72$site_totalarea_cent <- groupmeans_V72$site_totalarea - mean(groupmeans_V72$site_totalarea)
groupmeans_V72$site_subcort_vol_cent <- groupmeans_V72$site_subcort_vol - mean(groupmeans_V72$site_subcort_vol)

groupmeans_V73 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V73, mean)
names(groupmeans_V73) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V73$site_meanct_cent <- groupmeans_V73$site_meanct - mean(groupmeans_V73$site_meanct)
groupmeans_V73$site_totalvol_cent <- groupmeans_V73$site_totalvol - mean(groupmeans_V73$site_totalvol)
groupmeans_V73$site_totalarea_cent <- groupmeans_V73$site_totalarea - mean(groupmeans_V73$site_totalarea)
groupmeans_V73$site_subcort_vol_cent <- groupmeans_V73$site_subcort_vol - mean(groupmeans_V73$site_subcort_vol)

groupmeans_V74 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V74, mean)
names(groupmeans_V74) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V74$site_meanct_cent <- groupmeans_V74$site_meanct - mean(groupmeans_V74$site_meanct)
groupmeans_V74$site_totalvol_cent <- groupmeans_V74$site_totalvol - mean(groupmeans_V74$site_totalvol)
groupmeans_V74$site_totalarea_cent <- groupmeans_V74$site_totalarea - mean(groupmeans_V74$site_totalarea)
groupmeans_V74$site_subcort_vol_cent <- groupmeans_V74$site_subcort_vol - mean(groupmeans_V74$site_subcort_vol)

groupmeans_V75 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V75, mean)
names(groupmeans_V75) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V75$site_meanct_cent <- groupmeans_V75$site_meanct - mean(groupmeans_V75$site_meanct)
groupmeans_V75$site_totalvol_cent <- groupmeans_V75$site_totalvol - mean(groupmeans_V75$site_totalvol)
groupmeans_V75$site_totalarea_cent <- groupmeans_V75$site_totalarea - mean(groupmeans_V75$site_totalarea)
groupmeans_V75$site_subcort_vol_cent <- groupmeans_V75$site_subcort_vol - mean(groupmeans_V75$site_subcort_vol)

groupmeans_V76 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V76, mean)
names(groupmeans_V76) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V76$site_meanct_cent <- groupmeans_V76$site_meanct - mean(groupmeans_V76$site_meanct)
groupmeans_V76$site_totalvol_cent <- groupmeans_V76$site_totalvol - mean(groupmeans_V76$site_totalvol)
groupmeans_V76$site_totalarea_cent <- groupmeans_V76$site_totalarea - mean(groupmeans_V76$site_totalarea)
groupmeans_V76$site_subcort_vol_cent <- groupmeans_V76$site_subcort_vol - mean(groupmeans_V76$site_subcort_vol)

groupmeans_V77 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V77, mean)
names(groupmeans_V77) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V77$site_meanct_cent <- groupmeans_V77$site_meanct - mean(groupmeans_V77$site_meanct)
groupmeans_V77$site_totalvol_cent <- groupmeans_V77$site_totalvol - mean(groupmeans_V77$site_totalvol)
groupmeans_V77$site_totalarea_cent <- groupmeans_V77$site_totalarea - mean(groupmeans_V77$site_totalarea)
groupmeans_V77$site_subcort_vol_cent <- groupmeans_V77$site_subcort_vol - mean(groupmeans_V77$site_subcort_vol)

groupmeans_V78 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V78, mean)
names(groupmeans_V78) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V78$site_meanct_cent <- groupmeans_V78$site_meanct - mean(groupmeans_V78$site_meanct)
groupmeans_V78$site_totalvol_cent <- groupmeans_V78$site_totalvol - mean(groupmeans_V78$site_totalvol)
groupmeans_V78$site_totalarea_cent <- groupmeans_V78$site_totalarea - mean(groupmeans_V78$site_totalarea)
groupmeans_V78$site_subcort_vol_cent <- groupmeans_V78$site_subcort_vol - mean(groupmeans_V78$site_subcort_vol)

groupmeans_V79 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V79, mean)
names(groupmeans_V79) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V79$site_meanct_cent <- groupmeans_V79$site_meanct - mean(groupmeans_V79$site_meanct)
groupmeans_V79$site_totalvol_cent <- groupmeans_V79$site_totalvol - mean(groupmeans_V79$site_totalvol)
groupmeans_V79$site_totalarea_cent <- groupmeans_V79$site_totalarea - mean(groupmeans_V79$site_totalarea)
groupmeans_V79$site_subcort_vol_cent <- groupmeans_V79$site_subcort_vol - mean(groupmeans_V79$site_subcort_vol)

groupmeans_V80 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V80, mean)
names(groupmeans_V80) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V80$site_meanct_cent <- groupmeans_V80$site_meanct - mean(groupmeans_V80$site_meanct)
groupmeans_V80$site_totalvol_cent <- groupmeans_V80$site_totalvol - mean(groupmeans_V80$site_totalvol)
groupmeans_V80$site_totalarea_cent <- groupmeans_V80$site_totalarea - mean(groupmeans_V80$site_totalarea)
groupmeans_V80$site_subcort_vol_cent <- groupmeans_V80$site_subcort_vol - mean(groupmeans_V80$site_subcort_vol)

groupmeans_V81 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V81, mean)
names(groupmeans_V81) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V81$site_meanct_cent <- groupmeans_V81$site_meanct - mean(groupmeans_V81$site_meanct)
groupmeans_V81$site_totalvol_cent <- groupmeans_V81$site_totalvol - mean(groupmeans_V81$site_totalvol)
groupmeans_V81$site_totalarea_cent <- groupmeans_V81$site_totalarea - mean(groupmeans_V81$site_totalarea)
groupmeans_V81$site_subcort_vol_cent <- groupmeans_V81$site_subcort_vol - mean(groupmeans_V81$site_subcort_vol)

groupmeans_V82 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V82, mean)
names(groupmeans_V82) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V82$site_meanct_cent <- groupmeans_V82$site_meanct - mean(groupmeans_V82$site_meanct)
groupmeans_V82$site_totalvol_cent <- groupmeans_V82$site_totalvol - mean(groupmeans_V82$site_totalvol)
groupmeans_V82$site_totalarea_cent <- groupmeans_V82$site_totalarea - mean(groupmeans_V82$site_totalarea)
groupmeans_V82$site_subcort_vol_cent <- groupmeans_V82$site_subcort_vol - mean(groupmeans_V82$site_subcort_vol)

groupmeans_V83 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V83, mean)
names(groupmeans_V83) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V83$site_meanct_cent <- groupmeans_V83$site_meanct - mean(groupmeans_V83$site_meanct)
groupmeans_V83$site_totalvol_cent <- groupmeans_V83$site_totalvol - mean(groupmeans_V83$site_totalvol)
groupmeans_V83$site_totalarea_cent <- groupmeans_V83$site_totalarea - mean(groupmeans_V83$site_totalarea)
groupmeans_V83$site_subcort_vol_cent <- groupmeans_V83$site_subcort_vol - mean(groupmeans_V83$site_subcort_vol)

groupmeans_V84 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V84, mean)
names(groupmeans_V84) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V84$site_meanct_cent <- groupmeans_V84$site_meanct - mean(groupmeans_V84$site_meanct)
groupmeans_V84$site_totalvol_cent <- groupmeans_V84$site_totalvol - mean(groupmeans_V84$site_totalvol)
groupmeans_V84$site_totalarea_cent <- groupmeans_V84$site_totalarea - mean(groupmeans_V84$site_totalarea)
groupmeans_V84$site_subcort_vol_cent <- groupmeans_V84$site_subcort_vol - mean(groupmeans_V84$site_subcort_vol)

groupmeans_V85 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V85, mean)
names(groupmeans_V85) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V85$site_meanct_cent <- groupmeans_V85$site_meanct - mean(groupmeans_V85$site_meanct)
groupmeans_V85$site_totalvol_cent <- groupmeans_V85$site_totalvol - mean(groupmeans_V85$site_totalvol)
groupmeans_V85$site_totalarea_cent <- groupmeans_V85$site_totalarea - mean(groupmeans_V85$site_totalarea)
groupmeans_V85$site_subcort_vol_cent <- groupmeans_V85$site_subcort_vol - mean(groupmeans_V85$site_subcort_vol)

groupmeans_V86 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V86, mean)
names(groupmeans_V86) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V86$site_meanct_cent <- groupmeans_V86$site_meanct - mean(groupmeans_V86$site_meanct)
groupmeans_V86$site_totalvol_cent <- groupmeans_V86$site_totalvol - mean(groupmeans_V86$site_totalvol)
groupmeans_V86$site_totalarea_cent <- groupmeans_V86$site_totalarea - mean(groupmeans_V86$site_totalarea)
groupmeans_V86$site_subcort_vol_cent <- groupmeans_V86$site_subcort_vol - mean(groupmeans_V86$site_subcort_vol)

groupmeans_V87 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V87, mean)
names(groupmeans_V87) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V87$site_meanct_cent <- groupmeans_V87$site_meanct - mean(groupmeans_V87$site_meanct)
groupmeans_V87$site_totalvol_cent <- groupmeans_V87$site_totalvol - mean(groupmeans_V87$site_totalvol)
groupmeans_V87$site_totalarea_cent <- groupmeans_V87$site_totalarea - mean(groupmeans_V87$site_totalarea)
groupmeans_V87$site_subcort_vol_cent <- groupmeans_V87$site_subcort_vol - mean(groupmeans_V87$site_subcort_vol)

groupmeans_V88 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V88, mean)
names(groupmeans_V88) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V88$site_meanct_cent <- groupmeans_V88$site_meanct - mean(groupmeans_V88$site_meanct)
groupmeans_V88$site_totalvol_cent <- groupmeans_V88$site_totalvol - mean(groupmeans_V88$site_totalvol)
groupmeans_V88$site_totalarea_cent <- groupmeans_V88$site_totalarea - mean(groupmeans_V88$site_totalarea)
groupmeans_V88$site_subcort_vol_cent <- groupmeans_V88$site_subcort_vol - mean(groupmeans_V88$site_subcort_vol)

groupmeans_V89 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V89, mean)
names(groupmeans_V89) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V89$site_meanct_cent <- groupmeans_V89$site_meanct - mean(groupmeans_V89$site_meanct)
groupmeans_V89$site_totalvol_cent <- groupmeans_V89$site_totalvol - mean(groupmeans_V89$site_totalvol)
groupmeans_V89$site_totalarea_cent <- groupmeans_V89$site_totalarea - mean(groupmeans_V89$site_totalarea)
groupmeans_V89$site_subcort_vol_cent <- groupmeans_V89$site_subcort_vol - mean(groupmeans_V89$site_subcort_vol)

groupmeans_V90 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V90, mean)
names(groupmeans_V90) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V90$site_meanct_cent <- groupmeans_V90$site_meanct - mean(groupmeans_V90$site_meanct)
groupmeans_V90$site_totalvol_cent <- groupmeans_V90$site_totalvol - mean(groupmeans_V90$site_totalvol)
groupmeans_V90$site_totalarea_cent <- groupmeans_V90$site_totalarea - mean(groupmeans_V90$site_totalarea)
groupmeans_V90$site_subcort_vol_cent <- groupmeans_V90$site_subcort_vol - mean(groupmeans_V90$site_subcort_vol)

groupmeans_V91 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V91, mean)
names(groupmeans_V91) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V91$site_meanct_cent <- groupmeans_V91$site_meanct - mean(groupmeans_V91$site_meanct)
groupmeans_V91$site_totalvol_cent <- groupmeans_V91$site_totalvol - mean(groupmeans_V91$site_totalvol)
groupmeans_V91$site_totalarea_cent <- groupmeans_V91$site_totalarea - mean(groupmeans_V91$site_totalarea)
groupmeans_V91$site_subcort_vol_cent <- groupmeans_V91$site_subcort_vol - mean(groupmeans_V91$site_subcort_vol)

groupmeans_V92 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V92, mean)
names(groupmeans_V92) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V92$site_meanct_cent <- groupmeans_V92$site_meanct - mean(groupmeans_V92$site_meanct)
groupmeans_V92$site_totalvol_cent <- groupmeans_V92$site_totalvol - mean(groupmeans_V92$site_totalvol)
groupmeans_V92$site_totalarea_cent <- groupmeans_V92$site_totalarea - mean(groupmeans_V92$site_totalarea)
groupmeans_V92$site_subcort_vol_cent <- groupmeans_V92$site_subcort_vol - mean(groupmeans_V92$site_subcort_vol)

groupmeans_V93 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V93, mean)
names(groupmeans_V93) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V93$site_meanct_cent <- groupmeans_V93$site_meanct - mean(groupmeans_V93$site_meanct)
groupmeans_V93$site_totalvol_cent <- groupmeans_V93$site_totalvol - mean(groupmeans_V93$site_totalvol)
groupmeans_V93$site_totalarea_cent <- groupmeans_V93$site_totalarea - mean(groupmeans_V93$site_totalarea)
groupmeans_V93$site_subcort_vol_cent <- groupmeans_V93$site_subcort_vol - mean(groupmeans_V93$site_subcort_vol)

groupmeans_V94 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V94, mean)
names(groupmeans_V94) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V94$site_meanct_cent <- groupmeans_V94$site_meanct - mean(groupmeans_V94$site_meanct)
groupmeans_V94$site_totalvol_cent <- groupmeans_V94$site_totalvol - mean(groupmeans_V94$site_totalvol)
groupmeans_V94$site_totalarea_cent <- groupmeans_V94$site_totalarea - mean(groupmeans_V94$site_totalarea)
groupmeans_V94$site_subcort_vol_cent <- groupmeans_V94$site_subcort_vol - mean(groupmeans_V94$site_subcort_vol)

groupmeans_V95 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V95, mean)
names(groupmeans_V95) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V95$site_meanct_cent <- groupmeans_V95$site_meanct - mean(groupmeans_V95$site_meanct)
groupmeans_V95$site_totalvol_cent <- groupmeans_V95$site_totalvol - mean(groupmeans_V95$site_totalvol)
groupmeans_V95$site_totalarea_cent <- groupmeans_V95$site_totalarea - mean(groupmeans_V95$site_totalarea)
groupmeans_V95$site_subcort_vol_cent <- groupmeans_V95$site_subcort_vol - mean(groupmeans_V95$site_subcort_vol)

groupmeans_V96 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V96, mean)
names(groupmeans_V96) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V96$site_meanct_cent <- groupmeans_V96$site_meanct - mean(groupmeans_V96$site_meanct)
groupmeans_V96$site_totalvol_cent <- groupmeans_V96$site_totalvol - mean(groupmeans_V96$site_totalvol)
groupmeans_V96$site_totalarea_cent <- groupmeans_V96$site_totalarea - mean(groupmeans_V96$site_totalarea)
groupmeans_V96$site_subcort_vol_cent <- groupmeans_V96$site_subcort_vol - mean(groupmeans_V96$site_subcort_vol)

groupmeans_V97 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V97, mean)
names(groupmeans_V97) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V97$site_meanct_cent <- groupmeans_V97$site_meanct - mean(groupmeans_V97$site_meanct)
groupmeans_V97$site_totalvol_cent <- groupmeans_V97$site_totalvol - mean(groupmeans_V97$site_totalvol)
groupmeans_V97$site_totalarea_cent <- groupmeans_V97$site_totalarea - mean(groupmeans_V97$site_totalarea)
groupmeans_V97$site_subcort_vol_cent <- groupmeans_V97$site_subcort_vol - mean(groupmeans_V97$site_subcort_vol)

groupmeans_V98 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V98, mean)
names(groupmeans_V98) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V98$site_meanct_cent <- groupmeans_V98$site_meanct - mean(groupmeans_V98$site_meanct)
groupmeans_V98$site_totalvol_cent <- groupmeans_V98$site_totalvol - mean(groupmeans_V98$site_totalvol)
groupmeans_V98$site_totalarea_cent <- groupmeans_V98$site_totalarea - mean(groupmeans_V98$site_totalarea)
groupmeans_V98$site_subcort_vol_cent <- groupmeans_V98$site_subcort_vol - mean(groupmeans_V98$site_subcort_vol)

groupmeans_V99 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V99, mean)
names(groupmeans_V99) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V99$site_meanct_cent <- groupmeans_V99$site_meanct - mean(groupmeans_V99$site_meanct)
groupmeans_V99$site_totalvol_cent <- groupmeans_V99$site_totalvol - mean(groupmeans_V99$site_totalvol)
groupmeans_V99$site_totalarea_cent <- groupmeans_V99$site_totalarea - mean(groupmeans_V99$site_totalarea)
groupmeans_V99$site_subcort_vol_cent <- groupmeans_V99$site_subcort_vol - mean(groupmeans_V99$site_subcort_vol)

groupmeans_V100 <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_long_sibs_nosibs_V100, mean)
names(groupmeans_V100) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans_V100$site_meanct_cent <- groupmeans_V100$site_meanct - mean(groupmeans_V100$site_meanct)
groupmeans_V100$site_totalvol_cent <- groupmeans_V100$site_totalvol - mean(groupmeans_V100$site_totalvol)
groupmeans_V100$site_totalarea_cent <- groupmeans_V100$site_totalarea - mean(groupmeans_V100$site_totalarea)
groupmeans_V100$site_subcort_vol_cent <- groupmeans_V100$site_subcort_vol - mean(groupmeans_V100$site_subcort_vol)
groupmeans_V100$site_subcort_vol_cent <- groupmeans_V100$site_subcort_vol - mean(groupmeans_V100$site_subcort_vol)

#merge centered variables (groupmeans) with each dataset
pfactor_long_sibs_nosibs_V1 <- join(pfactor_long_sibs_nosibs_V1, groupmeans_V1, by="site_id")
pfactor_long_sibs_nosibs_V2 <- join(pfactor_long_sibs_nosibs_V2, groupmeans_V2, by="site_id")
pfactor_long_sibs_nosibs_V3 <- join(pfactor_long_sibs_nosibs_V3, groupmeans_V3, by="site_id")
pfactor_long_sibs_nosibs_V4 <- join(pfactor_long_sibs_nosibs_V4, groupmeans_V4, by="site_id")
pfactor_long_sibs_nosibs_V5 <- join(pfactor_long_sibs_nosibs_V5, groupmeans_V5, by="site_id")
pfactor_long_sibs_nosibs_V6 <- join(pfactor_long_sibs_nosibs_V6, groupmeans_V6, by="site_id")
pfactor_long_sibs_nosibs_V7 <- join(pfactor_long_sibs_nosibs_V7, groupmeans_V7, by="site_id")
pfactor_long_sibs_nosibs_V8 <- join(pfactor_long_sibs_nosibs_V8, groupmeans_V8, by="site_id")
pfactor_long_sibs_nosibs_V9 <- join(pfactor_long_sibs_nosibs_V9, groupmeans_V9, by="site_id")
pfactor_long_sibs_nosibs_V10 <- join(pfactor_long_sibs_nosibs_V10, groupmeans_V10, by="site_id")
pfactor_long_sibs_nosibs_V11 <- join(pfactor_long_sibs_nosibs_V11, groupmeans_V11, by="site_id")
pfactor_long_sibs_nosibs_V12 <- join(pfactor_long_sibs_nosibs_V12, groupmeans_V12, by="site_id")
pfactor_long_sibs_nosibs_V13 <- join(pfactor_long_sibs_nosibs_V13, groupmeans_V13, by="site_id")
pfactor_long_sibs_nosibs_V14 <- join(pfactor_long_sibs_nosibs_V14, groupmeans_V14, by="site_id")
pfactor_long_sibs_nosibs_V15 <- join(pfactor_long_sibs_nosibs_V15, groupmeans_V15, by="site_id")
pfactor_long_sibs_nosibs_V16 <- join(pfactor_long_sibs_nosibs_V16, groupmeans_V16, by="site_id")
pfactor_long_sibs_nosibs_V17 <- join(pfactor_long_sibs_nosibs_V17, groupmeans_V17, by="site_id")
pfactor_long_sibs_nosibs_V18 <- join(pfactor_long_sibs_nosibs_V18, groupmeans_V18, by="site_id")
pfactor_long_sibs_nosibs_V19 <- join(pfactor_long_sibs_nosibs_V19, groupmeans_V19, by="site_id")
pfactor_long_sibs_nosibs_V20 <- join(pfactor_long_sibs_nosibs_V20, groupmeans_V20, by="site_id")
pfactor_long_sibs_nosibs_V21 <- join(pfactor_long_sibs_nosibs_V21, groupmeans_V21, by="site_id")
pfactor_long_sibs_nosibs_V22 <- join(pfactor_long_sibs_nosibs_V22, groupmeans_V22, by="site_id")
pfactor_long_sibs_nosibs_V23 <- join(pfactor_long_sibs_nosibs_V23, groupmeans_V23, by="site_id")
pfactor_long_sibs_nosibs_V24 <- join(pfactor_long_sibs_nosibs_V24, groupmeans_V24, by="site_id")
pfactor_long_sibs_nosibs_V25 <- join(pfactor_long_sibs_nosibs_V25, groupmeans_V25, by="site_id")
pfactor_long_sibs_nosibs_V26 <- join(pfactor_long_sibs_nosibs_V26, groupmeans_V26, by="site_id")
pfactor_long_sibs_nosibs_V27 <- join(pfactor_long_sibs_nosibs_V27, groupmeans_V27, by="site_id")
pfactor_long_sibs_nosibs_V28 <- join(pfactor_long_sibs_nosibs_V28, groupmeans_V28, by="site_id")
pfactor_long_sibs_nosibs_V29 <- join(pfactor_long_sibs_nosibs_V29, groupmeans_V29, by="site_id")
pfactor_long_sibs_nosibs_V30 <- join(pfactor_long_sibs_nosibs_V30, groupmeans_V30, by="site_id")
pfactor_long_sibs_nosibs_V31 <- join(pfactor_long_sibs_nosibs_V31, groupmeans_V31, by="site_id")
pfactor_long_sibs_nosibs_V32 <- join(pfactor_long_sibs_nosibs_V32, groupmeans_V32, by="site_id")
pfactor_long_sibs_nosibs_V33 <- join(pfactor_long_sibs_nosibs_V33, groupmeans_V33, by="site_id")
pfactor_long_sibs_nosibs_V34 <- join(pfactor_long_sibs_nosibs_V34, groupmeans_V34, by="site_id")
pfactor_long_sibs_nosibs_V35 <- join(pfactor_long_sibs_nosibs_V35, groupmeans_V35, by="site_id")
pfactor_long_sibs_nosibs_V36 <- join(pfactor_long_sibs_nosibs_V36, groupmeans_V36, by="site_id")
pfactor_long_sibs_nosibs_V37 <- join(pfactor_long_sibs_nosibs_V37, groupmeans_V37, by="site_id")
pfactor_long_sibs_nosibs_V38 <- join(pfactor_long_sibs_nosibs_V38, groupmeans_V38, by="site_id")
pfactor_long_sibs_nosibs_V39 <- join(pfactor_long_sibs_nosibs_V39, groupmeans_V39, by="site_id")
pfactor_long_sibs_nosibs_V40 <- join(pfactor_long_sibs_nosibs_V40, groupmeans_V40, by="site_id")
pfactor_long_sibs_nosibs_V41 <- join(pfactor_long_sibs_nosibs_V41, groupmeans_V41, by="site_id")
pfactor_long_sibs_nosibs_V42 <- join(pfactor_long_sibs_nosibs_V42, groupmeans_V42, by="site_id")
pfactor_long_sibs_nosibs_V43 <- join(pfactor_long_sibs_nosibs_V43, groupmeans_V43, by="site_id")
pfactor_long_sibs_nosibs_V44 <- join(pfactor_long_sibs_nosibs_V44, groupmeans_V44, by="site_id")
pfactor_long_sibs_nosibs_V45 <- join(pfactor_long_sibs_nosibs_V45, groupmeans_V45, by="site_id")
pfactor_long_sibs_nosibs_V46 <- join(pfactor_long_sibs_nosibs_V46, groupmeans_V46, by="site_id")
pfactor_long_sibs_nosibs_V47 <- join(pfactor_long_sibs_nosibs_V47, groupmeans_V47, by="site_id")
pfactor_long_sibs_nosibs_V48 <- join(pfactor_long_sibs_nosibs_V48, groupmeans_V48, by="site_id")
pfactor_long_sibs_nosibs_V49 <- join(pfactor_long_sibs_nosibs_V49, groupmeans_V49, by="site_id")
pfactor_long_sibs_nosibs_V50 <- join(pfactor_long_sibs_nosibs_V50, groupmeans_V50, by="site_id")
pfactor_long_sibs_nosibs_V51 <- join(pfactor_long_sibs_nosibs_V51, groupmeans_V51, by="site_id")
pfactor_long_sibs_nosibs_V52 <- join(pfactor_long_sibs_nosibs_V52, groupmeans_V52, by="site_id")
pfactor_long_sibs_nosibs_V53 <- join(pfactor_long_sibs_nosibs_V53, groupmeans_V53, by="site_id")
pfactor_long_sibs_nosibs_V54 <- join(pfactor_long_sibs_nosibs_V54, groupmeans_V54, by="site_id")
pfactor_long_sibs_nosibs_V55 <- join(pfactor_long_sibs_nosibs_V55, groupmeans_V55, by="site_id")
pfactor_long_sibs_nosibs_V56 <- join(pfactor_long_sibs_nosibs_V56, groupmeans_V56, by="site_id")
pfactor_long_sibs_nosibs_V57 <- join(pfactor_long_sibs_nosibs_V57, groupmeans_V57, by="site_id")
pfactor_long_sibs_nosibs_V58 <- join(pfactor_long_sibs_nosibs_V58, groupmeans_V58, by="site_id")
pfactor_long_sibs_nosibs_V59 <- join(pfactor_long_sibs_nosibs_V59, groupmeans_V59, by="site_id")
pfactor_long_sibs_nosibs_V60 <- join(pfactor_long_sibs_nosibs_V60, groupmeans_V60, by="site_id")
pfactor_long_sibs_nosibs_V61 <- join(pfactor_long_sibs_nosibs_V61, groupmeans_V61, by="site_id")
pfactor_long_sibs_nosibs_V62 <- join(pfactor_long_sibs_nosibs_V62, groupmeans_V62, by="site_id")
pfactor_long_sibs_nosibs_V63 <- join(pfactor_long_sibs_nosibs_V63, groupmeans_V63, by="site_id")
pfactor_long_sibs_nosibs_V64 <- join(pfactor_long_sibs_nosibs_V64, groupmeans_V64, by="site_id")
pfactor_long_sibs_nosibs_V65 <- join(pfactor_long_sibs_nosibs_V65, groupmeans_V65, by="site_id")
pfactor_long_sibs_nosibs_V66 <- join(pfactor_long_sibs_nosibs_V66, groupmeans_V66, by="site_id")
pfactor_long_sibs_nosibs_V67 <- join(pfactor_long_sibs_nosibs_V67, groupmeans_V67, by="site_id")
pfactor_long_sibs_nosibs_V68 <- join(pfactor_long_sibs_nosibs_V68, groupmeans_V68, by="site_id")
pfactor_long_sibs_nosibs_V69 <- join(pfactor_long_sibs_nosibs_V69, groupmeans_V69, by="site_id")
pfactor_long_sibs_nosibs_V70 <- join(pfactor_long_sibs_nosibs_V70, groupmeans_V70, by="site_id")
pfactor_long_sibs_nosibs_V71 <- join(pfactor_long_sibs_nosibs_V71, groupmeans_V71, by="site_id")
pfactor_long_sibs_nosibs_V72 <- join(pfactor_long_sibs_nosibs_V72, groupmeans_V72, by="site_id")
pfactor_long_sibs_nosibs_V73 <- join(pfactor_long_sibs_nosibs_V73, groupmeans_V73, by="site_id")
pfactor_long_sibs_nosibs_V74 <- join(pfactor_long_sibs_nosibs_V74, groupmeans_V74, by="site_id")
pfactor_long_sibs_nosibs_V75 <- join(pfactor_long_sibs_nosibs_V75, groupmeans_V75, by="site_id")
pfactor_long_sibs_nosibs_V76 <- join(pfactor_long_sibs_nosibs_V76, groupmeans_V76, by="site_id")
pfactor_long_sibs_nosibs_V77 <- join(pfactor_long_sibs_nosibs_V77, groupmeans_V77, by="site_id")
pfactor_long_sibs_nosibs_V78 <- join(pfactor_long_sibs_nosibs_V78, groupmeans_V78, by="site_id")
pfactor_long_sibs_nosibs_V79 <- join(pfactor_long_sibs_nosibs_V79, groupmeans_V79, by="site_id")
pfactor_long_sibs_nosibs_V80 <- join(pfactor_long_sibs_nosibs_V80, groupmeans_V80, by="site_id")
pfactor_long_sibs_nosibs_V81 <- join(pfactor_long_sibs_nosibs_V81, groupmeans_V81, by="site_id")
pfactor_long_sibs_nosibs_V82 <- join(pfactor_long_sibs_nosibs_V82, groupmeans_V82, by="site_id")
pfactor_long_sibs_nosibs_V83 <- join(pfactor_long_sibs_nosibs_V83, groupmeans_V83, by="site_id")
pfactor_long_sibs_nosibs_V84 <- join(pfactor_long_sibs_nosibs_V84, groupmeans_V84, by="site_id")
pfactor_long_sibs_nosibs_V85 <- join(pfactor_long_sibs_nosibs_V85, groupmeans_V85, by="site_id")
pfactor_long_sibs_nosibs_V86 <- join(pfactor_long_sibs_nosibs_V86, groupmeans_V86, by="site_id")
pfactor_long_sibs_nosibs_V87 <- join(pfactor_long_sibs_nosibs_V87, groupmeans_V87, by="site_id")
pfactor_long_sibs_nosibs_V88 <- join(pfactor_long_sibs_nosibs_V88, groupmeans_V88, by="site_id")
pfactor_long_sibs_nosibs_V89 <- join(pfactor_long_sibs_nosibs_V89, groupmeans_V89, by="site_id")
pfactor_long_sibs_nosibs_V90 <- join(pfactor_long_sibs_nosibs_V90, groupmeans_V90, by="site_id")
pfactor_long_sibs_nosibs_V91 <- join(pfactor_long_sibs_nosibs_V91, groupmeans_V91, by="site_id")
pfactor_long_sibs_nosibs_V92 <- join(pfactor_long_sibs_nosibs_V92, groupmeans_V92, by="site_id")
pfactor_long_sibs_nosibs_V93 <- join(pfactor_long_sibs_nosibs_V93, groupmeans_V93, by="site_id")
pfactor_long_sibs_nosibs_V94 <- join(pfactor_long_sibs_nosibs_V94, groupmeans_V94, by="site_id")
pfactor_long_sibs_nosibs_V95 <- join(pfactor_long_sibs_nosibs_V95, groupmeans_V95, by="site_id")
pfactor_long_sibs_nosibs_V96 <- join(pfactor_long_sibs_nosibs_V96, groupmeans_V96, by="site_id")
pfactor_long_sibs_nosibs_V97 <- join(pfactor_long_sibs_nosibs_V97, groupmeans_V97, by="site_id")
pfactor_long_sibs_nosibs_V98 <- join(pfactor_long_sibs_nosibs_V98, groupmeans_V98, by="site_id")
pfactor_long_sibs_nosibs_V99 <- join(pfactor_long_sibs_nosibs_V99, groupmeans_V99, by="site_id")
pfactor_long_sibs_nosibs_V100 <- join(pfactor_long_sibs_nosibs_V100, groupmeans_V100, by="site_id")

#Center global brain structure variables by subject for each dataframe
pfactor_long_sibs_nosibs_V1$subjmeanct <- pfactor_long_sibs_nosibs_V1$meanwb_ct - pfactor_long_sibs_nosibs_V1$site_meanct
pfactor_long_sibs_nosibs_V1$subjtotalvol <- pfactor_long_sibs_nosibs_V1$wb_cort_vol - pfactor_long_sibs_nosibs_V1$site_totalvol
pfactor_long_sibs_nosibs_V1$subjtotalarea <- pfactor_long_sibs_nosibs_V1$wb_cort_area - pfactor_long_sibs_nosibs_V1$site_totalarea
pfactor_long_sibs_nosibs_V1$subjsubcort_vol <- pfactor_long_sibs_nosibs_V1$subcort_vol - pfactor_long_sibs_nosibs_V1$site_subcort_vol

pfactor_long_sibs_nosibs_V2$subjmeanct <- pfactor_long_sibs_nosibs_V2$meanwb_ct - pfactor_long_sibs_nosibs_V2$site_meanct
pfactor_long_sibs_nosibs_V2$subjtotalvol <- pfactor_long_sibs_nosibs_V2$wb_cort_vol - pfactor_long_sibs_nosibs_V2$site_totalvol
pfactor_long_sibs_nosibs_V2$subjtotalarea <- pfactor_long_sibs_nosibs_V2$wb_cort_area - pfactor_long_sibs_nosibs_V2$site_totalarea
pfactor_long_sibs_nosibs_V2$subjsubcort_vol <- pfactor_long_sibs_nosibs_V2$subcort_vol - pfactor_long_sibs_nosibs_V2$site_subcort_vol

pfactor_long_sibs_nosibs_V3$subjmeanct <- pfactor_long_sibs_nosibs_V3$meanwb_ct - pfactor_long_sibs_nosibs_V3$site_meanct
pfactor_long_sibs_nosibs_V3$subjtotalvol <- pfactor_long_sibs_nosibs_V3$wb_cort_vol - pfactor_long_sibs_nosibs_V3$site_totalvol
pfactor_long_sibs_nosibs_V3$subjtotalarea <- pfactor_long_sibs_nosibs_V3$wb_cort_area - pfactor_long_sibs_nosibs_V3$site_totalarea
pfactor_long_sibs_nosibs_V3$subjsubcort_vol <- pfactor_long_sibs_nosibs_V3$subcort_vol - pfactor_long_sibs_nosibs_V3$site_subcort_vol

pfactor_long_sibs_nosibs_V4$subjmeanct <- pfactor_long_sibs_nosibs_V4$meanwb_ct - pfactor_long_sibs_nosibs_V4$site_meanct
pfactor_long_sibs_nosibs_V4$subjtotalvol <- pfactor_long_sibs_nosibs_V4$wb_cort_vol - pfactor_long_sibs_nosibs_V4$site_totalvol
pfactor_long_sibs_nosibs_V4$subjtotalarea <- pfactor_long_sibs_nosibs_V4$wb_cort_area - pfactor_long_sibs_nosibs_V4$site_totalarea
pfactor_long_sibs_nosibs_V4$subjsubcort_vol <- pfactor_long_sibs_nosibs_V4$subcort_vol - pfactor_long_sibs_nosibs_V4$site_subcort_vol

pfactor_long_sibs_nosibs_V5$subjmeanct <- pfactor_long_sibs_nosibs_V5$meanwb_ct - pfactor_long_sibs_nosibs_V5$site_meanct
pfactor_long_sibs_nosibs_V5$subjtotalvol <- pfactor_long_sibs_nosibs_V5$wb_cort_vol - pfactor_long_sibs_nosibs_V5$site_totalvol
pfactor_long_sibs_nosibs_V5$subjtotalarea <- pfactor_long_sibs_nosibs_V5$wb_cort_area - pfactor_long_sibs_nosibs_V5$site_totalarea
pfactor_long_sibs_nosibs_V5$subjsubcort_vol <- pfactor_long_sibs_nosibs_V5$subcort_vol - pfactor_long_sibs_nosibs_V5$site_subcort_vol

pfactor_long_sibs_nosibs_V6$subjmeanct <- pfactor_long_sibs_nosibs_V6$meanwb_ct - pfactor_long_sibs_nosibs_V6$site_meanct
pfactor_long_sibs_nosibs_V6$subjtotalvol <- pfactor_long_sibs_nosibs_V6$wb_cort_vol - pfactor_long_sibs_nosibs_V6$site_totalvol
pfactor_long_sibs_nosibs_V6$subjtotalarea <- pfactor_long_sibs_nosibs_V6$wb_cort_area - pfactor_long_sibs_nosibs_V6$site_totalarea
pfactor_long_sibs_nosibs_V6$subjsubcort_vol <- pfactor_long_sibs_nosibs_V6$subcort_vol - pfactor_long_sibs_nosibs_V6$site_subcort_vol

pfactor_long_sibs_nosibs_V7$subjmeanct <- pfactor_long_sibs_nosibs_V7$meanwb_ct - pfactor_long_sibs_nosibs_V7$site_meanct
pfactor_long_sibs_nosibs_V7$subjtotalvol <- pfactor_long_sibs_nosibs_V7$wb_cort_vol - pfactor_long_sibs_nosibs_V7$site_totalvol
pfactor_long_sibs_nosibs_V7$subjtotalarea <- pfactor_long_sibs_nosibs_V7$wb_cort_area - pfactor_long_sibs_nosibs_V7$site_totalarea
pfactor_long_sibs_nosibs_V7$subjsubcort_vol <- pfactor_long_sibs_nosibs_V7$subcort_vol - pfactor_long_sibs_nosibs_V7$site_subcort_vol

pfactor_long_sibs_nosibs_V8$subjmeanct <- pfactor_long_sibs_nosibs_V8$meanwb_ct - pfactor_long_sibs_nosibs_V8$site_meanct
pfactor_long_sibs_nosibs_V8$subjtotalvol <- pfactor_long_sibs_nosibs_V8$wb_cort_vol - pfactor_long_sibs_nosibs_V8$site_totalvol
pfactor_long_sibs_nosibs_V8$subjtotalarea <- pfactor_long_sibs_nosibs_V8$wb_cort_area - pfactor_long_sibs_nosibs_V8$site_totalarea
pfactor_long_sibs_nosibs_V8$subjsubcort_vol <- pfactor_long_sibs_nosibs_V8$subcort_vol - pfactor_long_sibs_nosibs_V8$site_subcort_vol

pfactor_long_sibs_nosibs_V9$subjmeanct <- pfactor_long_sibs_nosibs_V9$meanwb_ct - pfactor_long_sibs_nosibs_V9$site_meanct
pfactor_long_sibs_nosibs_V9$subjtotalvol <- pfactor_long_sibs_nosibs_V9$wb_cort_vol - pfactor_long_sibs_nosibs_V9$site_totalvol
pfactor_long_sibs_nosibs_V9$subjtotalarea <- pfactor_long_sibs_nosibs_V9$wb_cort_area - pfactor_long_sibs_nosibs_V9$site_totalarea
pfactor_long_sibs_nosibs_V9$subjsubcort_vol <- pfactor_long_sibs_nosibs_V9$subcort_vol - pfactor_long_sibs_nosibs_V9$site_subcort_vol

pfactor_long_sibs_nosibs_V10$subjmeanct <- pfactor_long_sibs_nosibs_V10$meanwb_ct - pfactor_long_sibs_nosibs_V10$site_meanct
pfactor_long_sibs_nosibs_V10$subjtotalvol <- pfactor_long_sibs_nosibs_V10$wb_cort_vol - pfactor_long_sibs_nosibs_V10$site_totalvol
pfactor_long_sibs_nosibs_V10$subjtotalarea <- pfactor_long_sibs_nosibs_V10$wb_cort_area - pfactor_long_sibs_nosibs_V10$site_totalarea
pfactor_long_sibs_nosibs_V10$subjsubcort_vol <- pfactor_long_sibs_nosibs_V10$subcort_vol - pfactor_long_sibs_nosibs_V10$site_subcort_vol

pfactor_long_sibs_nosibs_V11$subjmeanct <- pfactor_long_sibs_nosibs_V11$meanwb_ct - pfactor_long_sibs_nosibs_V11$site_meanct
pfactor_long_sibs_nosibs_V11$subjtotalvol <- pfactor_long_sibs_nosibs_V11$wb_cort_vol - pfactor_long_sibs_nosibs_V11$site_totalvol
pfactor_long_sibs_nosibs_V11$subjtotalarea <- pfactor_long_sibs_nosibs_V11$wb_cort_area - pfactor_long_sibs_nosibs_V11$site_totalarea
pfactor_long_sibs_nosibs_V11$subjsubcort_vol <- pfactor_long_sibs_nosibs_V11$subcort_vol - pfactor_long_sibs_nosibs_V11$site_subcort_vol

pfactor_long_sibs_nosibs_V12$subjmeanct <- pfactor_long_sibs_nosibs_V12$meanwb_ct - pfactor_long_sibs_nosibs_V12$site_meanct
pfactor_long_sibs_nosibs_V12$subjtotalvol <- pfactor_long_sibs_nosibs_V12$wb_cort_vol - pfactor_long_sibs_nosibs_V12$site_totalvol
pfactor_long_sibs_nosibs_V12$subjtotalarea <- pfactor_long_sibs_nosibs_V12$wb_cort_area - pfactor_long_sibs_nosibs_V12$site_totalarea
pfactor_long_sibs_nosibs_V12$subjsubcort_vol <- pfactor_long_sibs_nosibs_V12$subcort_vol - pfactor_long_sibs_nosibs_V12$site_subcort_vol

pfactor_long_sibs_nosibs_V13$subjmeanct <- pfactor_long_sibs_nosibs_V13$meanwb_ct - pfactor_long_sibs_nosibs_V13$site_meanct
pfactor_long_sibs_nosibs_V13$subjtotalvol <- pfactor_long_sibs_nosibs_V13$wb_cort_vol - pfactor_long_sibs_nosibs_V13$site_totalvol
pfactor_long_sibs_nosibs_V13$subjtotalarea <- pfactor_long_sibs_nosibs_V13$wb_cort_area - pfactor_long_sibs_nosibs_V13$site_totalarea
pfactor_long_sibs_nosibs_V13$subjsubcort_vol <- pfactor_long_sibs_nosibs_V13$subcort_vol - pfactor_long_sibs_nosibs_V13$site_subcort_vol

pfactor_long_sibs_nosibs_V14$subjmeanct <- pfactor_long_sibs_nosibs_V14$meanwb_ct - pfactor_long_sibs_nosibs_V14$site_meanct
pfactor_long_sibs_nosibs_V14$subjtotalvol <- pfactor_long_sibs_nosibs_V14$wb_cort_vol - pfactor_long_sibs_nosibs_V14$site_totalvol
pfactor_long_sibs_nosibs_V14$subjtotalarea <- pfactor_long_sibs_nosibs_V14$wb_cort_area - pfactor_long_sibs_nosibs_V14$site_totalarea
pfactor_long_sibs_nosibs_V14$subjsubcort_vol <- pfactor_long_sibs_nosibs_V14$subcort_vol - pfactor_long_sibs_nosibs_V14$site_subcort_vol

pfactor_long_sibs_nosibs_V15$subjmeanct <- pfactor_long_sibs_nosibs_V15$meanwb_ct - pfactor_long_sibs_nosibs_V15$site_meanct
pfactor_long_sibs_nosibs_V15$subjtotalvol <- pfactor_long_sibs_nosibs_V15$wb_cort_vol - pfactor_long_sibs_nosibs_V15$site_totalvol
pfactor_long_sibs_nosibs_V15$subjtotalarea <- pfactor_long_sibs_nosibs_V15$wb_cort_area - pfactor_long_sibs_nosibs_V15$site_totalarea
pfactor_long_sibs_nosibs_V15$subjsubcort_vol <- pfactor_long_sibs_nosibs_V15$subcort_vol - pfactor_long_sibs_nosibs_V15$site_subcort_vol

pfactor_long_sibs_nosibs_V16$subjmeanct <- pfactor_long_sibs_nosibs_V16$meanwb_ct - pfactor_long_sibs_nosibs_V16$site_meanct
pfactor_long_sibs_nosibs_V16$subjtotalvol <- pfactor_long_sibs_nosibs_V16$wb_cort_vol - pfactor_long_sibs_nosibs_V16$site_totalvol
pfactor_long_sibs_nosibs_V16$subjtotalarea <- pfactor_long_sibs_nosibs_V16$wb_cort_area - pfactor_long_sibs_nosibs_V16$site_totalarea
pfactor_long_sibs_nosibs_V16$subjsubcort_vol <- pfactor_long_sibs_nosibs_V16$subcort_vol - pfactor_long_sibs_nosibs_V16$site_subcort_vol

pfactor_long_sibs_nosibs_V17$subjmeanct <- pfactor_long_sibs_nosibs_V17$meanwb_ct - pfactor_long_sibs_nosibs_V17$site_meanct
pfactor_long_sibs_nosibs_V17$subjtotalvol <- pfactor_long_sibs_nosibs_V17$wb_cort_vol - pfactor_long_sibs_nosibs_V17$site_totalvol
pfactor_long_sibs_nosibs_V17$subjtotalarea <- pfactor_long_sibs_nosibs_V17$wb_cort_area - pfactor_long_sibs_nosibs_V17$site_totalarea
pfactor_long_sibs_nosibs_V17$subjsubcort_vol <- pfactor_long_sibs_nosibs_V17$subcort_vol - pfactor_long_sibs_nosibs_V17$site_subcort_vol

pfactor_long_sibs_nosibs_V18$subjmeanct <- pfactor_long_sibs_nosibs_V18$meanwb_ct - pfactor_long_sibs_nosibs_V18$site_meanct
pfactor_long_sibs_nosibs_V18$subjtotalvol <- pfactor_long_sibs_nosibs_V18$wb_cort_vol - pfactor_long_sibs_nosibs_V18$site_totalvol
pfactor_long_sibs_nosibs_V18$subjtotalarea <- pfactor_long_sibs_nosibs_V18$wb_cort_area - pfactor_long_sibs_nosibs_V18$site_totalarea
pfactor_long_sibs_nosibs_V18$subjsubcort_vol <- pfactor_long_sibs_nosibs_V18$subcort_vol - pfactor_long_sibs_nosibs_V18$site_subcort_vol

pfactor_long_sibs_nosibs_V19$subjmeanct <- pfactor_long_sibs_nosibs_V19$meanwb_ct - pfactor_long_sibs_nosibs_V19$site_meanct
pfactor_long_sibs_nosibs_V19$subjtotalvol <- pfactor_long_sibs_nosibs_V19$wb_cort_vol - pfactor_long_sibs_nosibs_V19$site_totalvol
pfactor_long_sibs_nosibs_V19$subjtotalarea <- pfactor_long_sibs_nosibs_V19$wb_cort_area - pfactor_long_sibs_nosibs_V19$site_totalarea
pfactor_long_sibs_nosibs_V19$subjsubcort_vol <- pfactor_long_sibs_nosibs_V19$subcort_vol - pfactor_long_sibs_nosibs_V19$site_subcort_vol

pfactor_long_sibs_nosibs_V20$subjmeanct <- pfactor_long_sibs_nosibs_V20$meanwb_ct - pfactor_long_sibs_nosibs_V20$site_meanct
pfactor_long_sibs_nosibs_V20$subjtotalvol <- pfactor_long_sibs_nosibs_V20$wb_cort_vol - pfactor_long_sibs_nosibs_V20$site_totalvol
pfactor_long_sibs_nosibs_V20$subjtotalarea <- pfactor_long_sibs_nosibs_V20$wb_cort_area - pfactor_long_sibs_nosibs_V20$site_totalarea
pfactor_long_sibs_nosibs_V20$subjsubcort_vol <- pfactor_long_sibs_nosibs_V20$subcort_vol - pfactor_long_sibs_nosibs_V20$site_subcort_vol

pfactor_long_sibs_nosibs_V21$subjmeanct <- pfactor_long_sibs_nosibs_V21$meanwb_ct - pfactor_long_sibs_nosibs_V21$site_meanct
pfactor_long_sibs_nosibs_V21$subjtotalvol <- pfactor_long_sibs_nosibs_V21$wb_cort_vol - pfactor_long_sibs_nosibs_V21$site_totalvol
pfactor_long_sibs_nosibs_V21$subjtotalarea <- pfactor_long_sibs_nosibs_V21$wb_cort_area - pfactor_long_sibs_nosibs_V21$site_totalarea
pfactor_long_sibs_nosibs_V21$subjsubcort_vol <- pfactor_long_sibs_nosibs_V21$subcort_vol - pfactor_long_sibs_nosibs_V21$site_subcort_vol

pfactor_long_sibs_nosibs_V22$subjmeanct <- pfactor_long_sibs_nosibs_V22$meanwb_ct - pfactor_long_sibs_nosibs_V22$site_meanct
pfactor_long_sibs_nosibs_V22$subjtotalvol <- pfactor_long_sibs_nosibs_V22$wb_cort_vol - pfactor_long_sibs_nosibs_V22$site_totalvol
pfactor_long_sibs_nosibs_V22$subjtotalarea <- pfactor_long_sibs_nosibs_V22$wb_cort_area - pfactor_long_sibs_nosibs_V22$site_totalarea
pfactor_long_sibs_nosibs_V22$subjsubcort_vol <- pfactor_long_sibs_nosibs_V22$subcort_vol - pfactor_long_sibs_nosibs_V22$site_subcort_vol

pfactor_long_sibs_nosibs_V23$subjmeanct <- pfactor_long_sibs_nosibs_V23$meanwb_ct - pfactor_long_sibs_nosibs_V23$site_meanct
pfactor_long_sibs_nosibs_V23$subjtotalvol <- pfactor_long_sibs_nosibs_V23$wb_cort_vol - pfactor_long_sibs_nosibs_V23$site_totalvol
pfactor_long_sibs_nosibs_V23$subjtotalarea <- pfactor_long_sibs_nosibs_V23$wb_cort_area - pfactor_long_sibs_nosibs_V23$site_totalarea
pfactor_long_sibs_nosibs_V23$subjsubcort_vol <- pfactor_long_sibs_nosibs_V23$subcort_vol - pfactor_long_sibs_nosibs_V23$site_subcort_vol

pfactor_long_sibs_nosibs_V24$subjmeanct <- pfactor_long_sibs_nosibs_V24$meanwb_ct - pfactor_long_sibs_nosibs_V24$site_meanct
pfactor_long_sibs_nosibs_V24$subjtotalvol <- pfactor_long_sibs_nosibs_V24$wb_cort_vol - pfactor_long_sibs_nosibs_V24$site_totalvol
pfactor_long_sibs_nosibs_V24$subjtotalarea <- pfactor_long_sibs_nosibs_V24$wb_cort_area - pfactor_long_sibs_nosibs_V24$site_totalarea
pfactor_long_sibs_nosibs_V24$subjsubcort_vol <- pfactor_long_sibs_nosibs_V24$subcort_vol - pfactor_long_sibs_nosibs_V24$site_subcort_vol

pfactor_long_sibs_nosibs_V25$subjmeanct <- pfactor_long_sibs_nosibs_V25$meanwb_ct - pfactor_long_sibs_nosibs_V25$site_meanct
pfactor_long_sibs_nosibs_V25$subjtotalvol <- pfactor_long_sibs_nosibs_V25$wb_cort_vol - pfactor_long_sibs_nosibs_V25$site_totalvol
pfactor_long_sibs_nosibs_V25$subjtotalarea <- pfactor_long_sibs_nosibs_V25$wb_cort_area - pfactor_long_sibs_nosibs_V25$site_totalarea
pfactor_long_sibs_nosibs_V25$subjsubcort_vol <- pfactor_long_sibs_nosibs_V25$subcort_vol - pfactor_long_sibs_nosibs_V25$site_subcort_vol

pfactor_long_sibs_nosibs_V26$subjmeanct <- pfactor_long_sibs_nosibs_V26$meanwb_ct - pfactor_long_sibs_nosibs_V26$site_meanct
pfactor_long_sibs_nosibs_V26$subjtotalvol <- pfactor_long_sibs_nosibs_V26$wb_cort_vol - pfactor_long_sibs_nosibs_V26$site_totalvol
pfactor_long_sibs_nosibs_V26$subjtotalarea <- pfactor_long_sibs_nosibs_V26$wb_cort_area - pfactor_long_sibs_nosibs_V26$site_totalarea
pfactor_long_sibs_nosibs_V26$subjsubcort_vol <- pfactor_long_sibs_nosibs_V26$subcort_vol - pfactor_long_sibs_nosibs_V26$site_subcort_vol

pfactor_long_sibs_nosibs_V27$subjmeanct <- pfactor_long_sibs_nosibs_V27$meanwb_ct - pfactor_long_sibs_nosibs_V27$site_meanct
pfactor_long_sibs_nosibs_V27$subjtotalvol <- pfactor_long_sibs_nosibs_V27$wb_cort_vol - pfactor_long_sibs_nosibs_V27$site_totalvol
pfactor_long_sibs_nosibs_V27$subjtotalarea <- pfactor_long_sibs_nosibs_V27$wb_cort_area - pfactor_long_sibs_nosibs_V27$site_totalarea
pfactor_long_sibs_nosibs_V27$subjsubcort_vol <- pfactor_long_sibs_nosibs_V27$subcort_vol - pfactor_long_sibs_nosibs_V27$site_subcort_vol

pfactor_long_sibs_nosibs_V28$subjmeanct <- pfactor_long_sibs_nosibs_V28$meanwb_ct - pfactor_long_sibs_nosibs_V28$site_meanct
pfactor_long_sibs_nosibs_V28$subjtotalvol <- pfactor_long_sibs_nosibs_V28$wb_cort_vol - pfactor_long_sibs_nosibs_V28$site_totalvol
pfactor_long_sibs_nosibs_V28$subjtotalarea <- pfactor_long_sibs_nosibs_V28$wb_cort_area - pfactor_long_sibs_nosibs_V28$site_totalarea
pfactor_long_sibs_nosibs_V28$subjsubcort_vol <- pfactor_long_sibs_nosibs_V28$subcort_vol - pfactor_long_sibs_nosibs_V28$site_subcort_vol

pfactor_long_sibs_nosibs_V29$subjmeanct <- pfactor_long_sibs_nosibs_V29$meanwb_ct - pfactor_long_sibs_nosibs_V29$site_meanct
pfactor_long_sibs_nosibs_V29$subjtotalvol <- pfactor_long_sibs_nosibs_V29$wb_cort_vol - pfactor_long_sibs_nosibs_V29$site_totalvol
pfactor_long_sibs_nosibs_V29$subjtotalarea <- pfactor_long_sibs_nosibs_V29$wb_cort_area - pfactor_long_sibs_nosibs_V29$site_totalarea
pfactor_long_sibs_nosibs_V29$subjsubcort_vol <- pfactor_long_sibs_nosibs_V29$subcort_vol - pfactor_long_sibs_nosibs_V29$site_subcort_vol

pfactor_long_sibs_nosibs_V30$subjmeanct <- pfactor_long_sibs_nosibs_V30$meanwb_ct - pfactor_long_sibs_nosibs_V30$site_meanct
pfactor_long_sibs_nosibs_V30$subjtotalvol <- pfactor_long_sibs_nosibs_V30$wb_cort_vol - pfactor_long_sibs_nosibs_V30$site_totalvol
pfactor_long_sibs_nosibs_V30$subjtotalarea <- pfactor_long_sibs_nosibs_V30$wb_cort_area - pfactor_long_sibs_nosibs_V30$site_totalarea
pfactor_long_sibs_nosibs_V30$subjsubcort_vol <- pfactor_long_sibs_nosibs_V30$subcort_vol - pfactor_long_sibs_nosibs_V30$site_subcort_vol

pfactor_long_sibs_nosibs_V31$subjmeanct <- pfactor_long_sibs_nosibs_V31$meanwb_ct - pfactor_long_sibs_nosibs_V31$site_meanct
pfactor_long_sibs_nosibs_V31$subjtotalvol <- pfactor_long_sibs_nosibs_V31$wb_cort_vol - pfactor_long_sibs_nosibs_V31$site_totalvol
pfactor_long_sibs_nosibs_V31$subjtotalarea <- pfactor_long_sibs_nosibs_V31$wb_cort_area - pfactor_long_sibs_nosibs_V31$site_totalarea
pfactor_long_sibs_nosibs_V31$subjsubcort_vol <- pfactor_long_sibs_nosibs_V31$subcort_vol - pfactor_long_sibs_nosibs_V31$site_subcort_vol

pfactor_long_sibs_nosibs_V32$subjmeanct <- pfactor_long_sibs_nosibs_V32$meanwb_ct - pfactor_long_sibs_nosibs_V32$site_meanct
pfactor_long_sibs_nosibs_V32$subjtotalvol <- pfactor_long_sibs_nosibs_V32$wb_cort_vol - pfactor_long_sibs_nosibs_V32$site_totalvol
pfactor_long_sibs_nosibs_V32$subjtotalarea <- pfactor_long_sibs_nosibs_V32$wb_cort_area - pfactor_long_sibs_nosibs_V32$site_totalarea
pfactor_long_sibs_nosibs_V32$subjsubcort_vol <- pfactor_long_sibs_nosibs_V32$subcort_vol - pfactor_long_sibs_nosibs_V32$site_subcort_vol

pfactor_long_sibs_nosibs_V33$subjmeanct <- pfactor_long_sibs_nosibs_V33$meanwb_ct - pfactor_long_sibs_nosibs_V33$site_meanct
pfactor_long_sibs_nosibs_V33$subjtotalvol <- pfactor_long_sibs_nosibs_V33$wb_cort_vol - pfactor_long_sibs_nosibs_V33$site_totalvol
pfactor_long_sibs_nosibs_V33$subjtotalarea <- pfactor_long_sibs_nosibs_V33$wb_cort_area - pfactor_long_sibs_nosibs_V33$site_totalarea
pfactor_long_sibs_nosibs_V33$subjsubcort_vol <- pfactor_long_sibs_nosibs_V33$subcort_vol - pfactor_long_sibs_nosibs_V33$site_subcort_vol

pfactor_long_sibs_nosibs_V34$subjmeanct <- pfactor_long_sibs_nosibs_V34$meanwb_ct - pfactor_long_sibs_nosibs_V34$site_meanct
pfactor_long_sibs_nosibs_V34$subjtotalvol <- pfactor_long_sibs_nosibs_V34$wb_cort_vol - pfactor_long_sibs_nosibs_V34$site_totalvol
pfactor_long_sibs_nosibs_V34$subjtotalarea <- pfactor_long_sibs_nosibs_V34$wb_cort_area - pfactor_long_sibs_nosibs_V34$site_totalarea
pfactor_long_sibs_nosibs_V34$subjsubcort_vol <- pfactor_long_sibs_nosibs_V34$subcort_vol - pfactor_long_sibs_nosibs_V34$site_subcort_vol

pfactor_long_sibs_nosibs_V35$subjmeanct <- pfactor_long_sibs_nosibs_V35$meanwb_ct - pfactor_long_sibs_nosibs_V35$site_meanct
pfactor_long_sibs_nosibs_V35$subjtotalvol <- pfactor_long_sibs_nosibs_V35$wb_cort_vol - pfactor_long_sibs_nosibs_V35$site_totalvol
pfactor_long_sibs_nosibs_V35$subjtotalarea <- pfactor_long_sibs_nosibs_V35$wb_cort_area - pfactor_long_sibs_nosibs_V35$site_totalarea
pfactor_long_sibs_nosibs_V35$subjsubcort_vol <- pfactor_long_sibs_nosibs_V35$subcort_vol - pfactor_long_sibs_nosibs_V35$site_subcort_vol

pfactor_long_sibs_nosibs_V36$subjmeanct <- pfactor_long_sibs_nosibs_V36$meanwb_ct - pfactor_long_sibs_nosibs_V36$site_meanct
pfactor_long_sibs_nosibs_V36$subjtotalvol <- pfactor_long_sibs_nosibs_V36$wb_cort_vol - pfactor_long_sibs_nosibs_V36$site_totalvol
pfactor_long_sibs_nosibs_V36$subjtotalarea <- pfactor_long_sibs_nosibs_V36$wb_cort_area - pfactor_long_sibs_nosibs_V36$site_totalarea
pfactor_long_sibs_nosibs_V36$subjsubcort_vol <- pfactor_long_sibs_nosibs_V36$subcort_vol - pfactor_long_sibs_nosibs_V36$site_subcort_vol

pfactor_long_sibs_nosibs_V37$subjmeanct <- pfactor_long_sibs_nosibs_V37$meanwb_ct - pfactor_long_sibs_nosibs_V37$site_meanct
pfactor_long_sibs_nosibs_V37$subjtotalvol <- pfactor_long_sibs_nosibs_V37$wb_cort_vol - pfactor_long_sibs_nosibs_V37$site_totalvol
pfactor_long_sibs_nosibs_V37$subjtotalarea <- pfactor_long_sibs_nosibs_V37$wb_cort_area - pfactor_long_sibs_nosibs_V37$site_totalarea
pfactor_long_sibs_nosibs_V37$subjsubcort_vol <- pfactor_long_sibs_nosibs_V37$subcort_vol - pfactor_long_sibs_nosibs_V37$site_subcort_vol

pfactor_long_sibs_nosibs_V38$subjmeanct <- pfactor_long_sibs_nosibs_V38$meanwb_ct - pfactor_long_sibs_nosibs_V38$site_meanct
pfactor_long_sibs_nosibs_V38$subjtotalvol <- pfactor_long_sibs_nosibs_V38$wb_cort_vol - pfactor_long_sibs_nosibs_V38$site_totalvol
pfactor_long_sibs_nosibs_V38$subjtotalarea <- pfactor_long_sibs_nosibs_V38$wb_cort_area - pfactor_long_sibs_nosibs_V38$site_totalarea
pfactor_long_sibs_nosibs_V38$subjsubcort_vol <- pfactor_long_sibs_nosibs_V38$subcort_vol - pfactor_long_sibs_nosibs_V38$site_subcort_vol

pfactor_long_sibs_nosibs_V39$subjmeanct <- pfactor_long_sibs_nosibs_V39$meanwb_ct - pfactor_long_sibs_nosibs_V39$site_meanct
pfactor_long_sibs_nosibs_V39$subjtotalvol <- pfactor_long_sibs_nosibs_V39$wb_cort_vol - pfactor_long_sibs_nosibs_V39$site_totalvol
pfactor_long_sibs_nosibs_V39$subjtotalarea <- pfactor_long_sibs_nosibs_V39$wb_cort_area - pfactor_long_sibs_nosibs_V39$site_totalarea
pfactor_long_sibs_nosibs_V39$subjsubcort_vol <- pfactor_long_sibs_nosibs_V39$subcort_vol - pfactor_long_sibs_nosibs_V39$site_subcort_vol

pfactor_long_sibs_nosibs_V40$subjmeanct <- pfactor_long_sibs_nosibs_V40$meanwb_ct - pfactor_long_sibs_nosibs_V40$site_meanct
pfactor_long_sibs_nosibs_V40$subjtotalvol <- pfactor_long_sibs_nosibs_V40$wb_cort_vol - pfactor_long_sibs_nosibs_V40$site_totalvol
pfactor_long_sibs_nosibs_V40$subjtotalarea <- pfactor_long_sibs_nosibs_V40$wb_cort_area - pfactor_long_sibs_nosibs_V40$site_totalarea
pfactor_long_sibs_nosibs_V40$subjsubcort_vol <- pfactor_long_sibs_nosibs_V40$subcort_vol - pfactor_long_sibs_nosibs_V40$site_subcort_vol

pfactor_long_sibs_nosibs_V41$subjmeanct <- pfactor_long_sibs_nosibs_V41$meanwb_ct - pfactor_long_sibs_nosibs_V41$site_meanct
pfactor_long_sibs_nosibs_V41$subjtotalvol <- pfactor_long_sibs_nosibs_V41$wb_cort_vol - pfactor_long_sibs_nosibs_V41$site_totalvol
pfactor_long_sibs_nosibs_V41$subjtotalarea <- pfactor_long_sibs_nosibs_V41$wb_cort_area - pfactor_long_sibs_nosibs_V41$site_totalarea
pfactor_long_sibs_nosibs_V41$subjsubcort_vol <- pfactor_long_sibs_nosibs_V41$subcort_vol - pfactor_long_sibs_nosibs_V41$site_subcort_vol

pfactor_long_sibs_nosibs_V42$subjmeanct <- pfactor_long_sibs_nosibs_V42$meanwb_ct - pfactor_long_sibs_nosibs_V42$site_meanct
pfactor_long_sibs_nosibs_V42$subjtotalvol <- pfactor_long_sibs_nosibs_V42$wb_cort_vol - pfactor_long_sibs_nosibs_V42$site_totalvol
pfactor_long_sibs_nosibs_V42$subjtotalarea <- pfactor_long_sibs_nosibs_V42$wb_cort_area - pfactor_long_sibs_nosibs_V42$site_totalarea
pfactor_long_sibs_nosibs_V42$subjsubcort_vol <- pfactor_long_sibs_nosibs_V42$subcort_vol - pfactor_long_sibs_nosibs_V42$site_subcort_vol

pfactor_long_sibs_nosibs_V43$subjmeanct <- pfactor_long_sibs_nosibs_V43$meanwb_ct - pfactor_long_sibs_nosibs_V43$site_meanct
pfactor_long_sibs_nosibs_V43$subjtotalvol <- pfactor_long_sibs_nosibs_V43$wb_cort_vol - pfactor_long_sibs_nosibs_V43$site_totalvol
pfactor_long_sibs_nosibs_V43$subjtotalarea <- pfactor_long_sibs_nosibs_V43$wb_cort_area - pfactor_long_sibs_nosibs_V43$site_totalarea
pfactor_long_sibs_nosibs_V43$subjsubcort_vol <- pfactor_long_sibs_nosibs_V43$subcort_vol - pfactor_long_sibs_nosibs_V43$site_subcort_vol

pfactor_long_sibs_nosibs_V44$subjmeanct <- pfactor_long_sibs_nosibs_V44$meanwb_ct - pfactor_long_sibs_nosibs_V44$site_meanct
pfactor_long_sibs_nosibs_V44$subjtotalvol <- pfactor_long_sibs_nosibs_V44$wb_cort_vol - pfactor_long_sibs_nosibs_V44$site_totalvol
pfactor_long_sibs_nosibs_V44$subjtotalarea <- pfactor_long_sibs_nosibs_V44$wb_cort_area - pfactor_long_sibs_nosibs_V44$site_totalarea
pfactor_long_sibs_nosibs_V44$subjsubcort_vol <- pfactor_long_sibs_nosibs_V44$subcort_vol - pfactor_long_sibs_nosibs_V44$site_subcort_vol

pfactor_long_sibs_nosibs_V45$subjmeanct <- pfactor_long_sibs_nosibs_V45$meanwb_ct - pfactor_long_sibs_nosibs_V45$site_meanct
pfactor_long_sibs_nosibs_V45$subjtotalvol <- pfactor_long_sibs_nosibs_V45$wb_cort_vol - pfactor_long_sibs_nosibs_V45$site_totalvol
pfactor_long_sibs_nosibs_V45$subjtotalarea <- pfactor_long_sibs_nosibs_V45$wb_cort_area - pfactor_long_sibs_nosibs_V45$site_totalarea
pfactor_long_sibs_nosibs_V45$subjsubcort_vol <- pfactor_long_sibs_nosibs_V45$subcort_vol - pfactor_long_sibs_nosibs_V45$site_subcort_vol

pfactor_long_sibs_nosibs_V46$subjmeanct <- pfactor_long_sibs_nosibs_V46$meanwb_ct - pfactor_long_sibs_nosibs_V46$site_meanct
pfactor_long_sibs_nosibs_V46$subjtotalvol <- pfactor_long_sibs_nosibs_V46$wb_cort_vol - pfactor_long_sibs_nosibs_V46$site_totalvol
pfactor_long_sibs_nosibs_V46$subjtotalarea <- pfactor_long_sibs_nosibs_V46$wb_cort_area - pfactor_long_sibs_nosibs_V46$site_totalarea
pfactor_long_sibs_nosibs_V46$subjsubcort_vol <- pfactor_long_sibs_nosibs_V46$subcort_vol - pfactor_long_sibs_nosibs_V46$site_subcort_vol

pfactor_long_sibs_nosibs_V47$subjmeanct <- pfactor_long_sibs_nosibs_V47$meanwb_ct - pfactor_long_sibs_nosibs_V47$site_meanct
pfactor_long_sibs_nosibs_V47$subjtotalvol <- pfactor_long_sibs_nosibs_V47$wb_cort_vol - pfactor_long_sibs_nosibs_V47$site_totalvol
pfactor_long_sibs_nosibs_V47$subjtotalarea <- pfactor_long_sibs_nosibs_V47$wb_cort_area - pfactor_long_sibs_nosibs_V47$site_totalarea
pfactor_long_sibs_nosibs_V47$subjsubcort_vol <- pfactor_long_sibs_nosibs_V47$subcort_vol - pfactor_long_sibs_nosibs_V47$site_subcort_vol

pfactor_long_sibs_nosibs_V48$subjmeanct <- pfactor_long_sibs_nosibs_V48$meanwb_ct - pfactor_long_sibs_nosibs_V48$site_meanct
pfactor_long_sibs_nosibs_V48$subjtotalvol <- pfactor_long_sibs_nosibs_V48$wb_cort_vol - pfactor_long_sibs_nosibs_V48$site_totalvol
pfactor_long_sibs_nosibs_V48$subjtotalarea <- pfactor_long_sibs_nosibs_V48$wb_cort_area - pfactor_long_sibs_nosibs_V48$site_totalarea
pfactor_long_sibs_nosibs_V48$subjsubcort_vol <- pfactor_long_sibs_nosibs_V48$subcort_vol - pfactor_long_sibs_nosibs_V48$site_subcort_vol

pfactor_long_sibs_nosibs_V49$subjmeanct <- pfactor_long_sibs_nosibs_V49$meanwb_ct - pfactor_long_sibs_nosibs_V49$site_meanct
pfactor_long_sibs_nosibs_V49$subjtotalvol <- pfactor_long_sibs_nosibs_V49$wb_cort_vol - pfactor_long_sibs_nosibs_V49$site_totalvol
pfactor_long_sibs_nosibs_V49$subjtotalarea <- pfactor_long_sibs_nosibs_V49$wb_cort_area - pfactor_long_sibs_nosibs_V49$site_totalarea
pfactor_long_sibs_nosibs_V49$subjsubcort_vol <- pfactor_long_sibs_nosibs_V49$subcort_vol - pfactor_long_sibs_nosibs_V49$site_subcort_vol

pfactor_long_sibs_nosibs_V50$subjmeanct <- pfactor_long_sibs_nosibs_V50$meanwb_ct - pfactor_long_sibs_nosibs_V50$site_meanct
pfactor_long_sibs_nosibs_V50$subjtotalvol <- pfactor_long_sibs_nosibs_V50$wb_cort_vol - pfactor_long_sibs_nosibs_V50$site_totalvol
pfactor_long_sibs_nosibs_V50$subjtotalarea <- pfactor_long_sibs_nosibs_V50$wb_cort_area - pfactor_long_sibs_nosibs_V50$site_totalarea
pfactor_long_sibs_nosibs_V50$subjsubcort_vol <- pfactor_long_sibs_nosibs_V50$subcort_vol - pfactor_long_sibs_nosibs_V50$site_subcort_vol

pfactor_long_sibs_nosibs_V51$subjmeanct <- pfactor_long_sibs_nosibs_V51$meanwb_ct - pfactor_long_sibs_nosibs_V51$site_meanct
pfactor_long_sibs_nosibs_V51$subjtotalvol <- pfactor_long_sibs_nosibs_V51$wb_cort_vol - pfactor_long_sibs_nosibs_V51$site_totalvol
pfactor_long_sibs_nosibs_V51$subjtotalarea <- pfactor_long_sibs_nosibs_V51$wb_cort_area - pfactor_long_sibs_nosibs_V51$site_totalarea
pfactor_long_sibs_nosibs_V51$subjsubcort_vol <- pfactor_long_sibs_nosibs_V51$subcort_vol - pfactor_long_sibs_nosibs_V51$site_subcort_vol

pfactor_long_sibs_nosibs_V52$subjmeanct <- pfactor_long_sibs_nosibs_V52$meanwb_ct - pfactor_long_sibs_nosibs_V52$site_meanct
pfactor_long_sibs_nosibs_V52$subjtotalvol <- pfactor_long_sibs_nosibs_V52$wb_cort_vol - pfactor_long_sibs_nosibs_V52$site_totalvol
pfactor_long_sibs_nosibs_V52$subjtotalarea <- pfactor_long_sibs_nosibs_V52$wb_cort_area - pfactor_long_sibs_nosibs_V52$site_totalarea
pfactor_long_sibs_nosibs_V52$subjsubcort_vol <- pfactor_long_sibs_nosibs_V52$subcort_vol - pfactor_long_sibs_nosibs_V52$site_subcort_vol

pfactor_long_sibs_nosibs_V53$subjmeanct <- pfactor_long_sibs_nosibs_V53$meanwb_ct - pfactor_long_sibs_nosibs_V53$site_meanct
pfactor_long_sibs_nosibs_V53$subjtotalvol <- pfactor_long_sibs_nosibs_V53$wb_cort_vol - pfactor_long_sibs_nosibs_V53$site_totalvol
pfactor_long_sibs_nosibs_V53$subjtotalarea <- pfactor_long_sibs_nosibs_V53$wb_cort_area - pfactor_long_sibs_nosibs_V53$site_totalarea
pfactor_long_sibs_nosibs_V53$subjsubcort_vol <- pfactor_long_sibs_nosibs_V53$subcort_vol - pfactor_long_sibs_nosibs_V53$site_subcort_vol

pfactor_long_sibs_nosibs_V54$subjmeanct <- pfactor_long_sibs_nosibs_V54$meanwb_ct - pfactor_long_sibs_nosibs_V54$site_meanct
pfactor_long_sibs_nosibs_V54$subjtotalvol <- pfactor_long_sibs_nosibs_V54$wb_cort_vol - pfactor_long_sibs_nosibs_V54$site_totalvol
pfactor_long_sibs_nosibs_V54$subjtotalarea <- pfactor_long_sibs_nosibs_V54$wb_cort_area - pfactor_long_sibs_nosibs_V54$site_totalarea
pfactor_long_sibs_nosibs_V54$subjsubcort_vol <- pfactor_long_sibs_nosibs_V54$subcort_vol - pfactor_long_sibs_nosibs_V54$site_subcort_vol

pfactor_long_sibs_nosibs_V55$subjmeanct <- pfactor_long_sibs_nosibs_V55$meanwb_ct - pfactor_long_sibs_nosibs_V55$site_meanct
pfactor_long_sibs_nosibs_V55$subjtotalvol <- pfactor_long_sibs_nosibs_V55$wb_cort_vol - pfactor_long_sibs_nosibs_V55$site_totalvol
pfactor_long_sibs_nosibs_V55$subjtotalarea <- pfactor_long_sibs_nosibs_V55$wb_cort_area - pfactor_long_sibs_nosibs_V55$site_totalarea
pfactor_long_sibs_nosibs_V55$subjsubcort_vol <- pfactor_long_sibs_nosibs_V55$subcort_vol - pfactor_long_sibs_nosibs_V55$site_subcort_vol

pfactor_long_sibs_nosibs_V56$subjmeanct <- pfactor_long_sibs_nosibs_V56$meanwb_ct - pfactor_long_sibs_nosibs_V56$site_meanct
pfactor_long_sibs_nosibs_V56$subjtotalvol <- pfactor_long_sibs_nosibs_V56$wb_cort_vol - pfactor_long_sibs_nosibs_V56$site_totalvol
pfactor_long_sibs_nosibs_V56$subjtotalarea <- pfactor_long_sibs_nosibs_V56$wb_cort_area - pfactor_long_sibs_nosibs_V56$site_totalarea
pfactor_long_sibs_nosibs_V56$subjsubcort_vol <- pfactor_long_sibs_nosibs_V56$subcort_vol - pfactor_long_sibs_nosibs_V56$site_subcort_vol

pfactor_long_sibs_nosibs_V57$subjmeanct <- pfactor_long_sibs_nosibs_V57$meanwb_ct - pfactor_long_sibs_nosibs_V57$site_meanct
pfactor_long_sibs_nosibs_V57$subjtotalvol <- pfactor_long_sibs_nosibs_V57$wb_cort_vol - pfactor_long_sibs_nosibs_V57$site_totalvol
pfactor_long_sibs_nosibs_V57$subjtotalarea <- pfactor_long_sibs_nosibs_V57$wb_cort_area - pfactor_long_sibs_nosibs_V57$site_totalarea
pfactor_long_sibs_nosibs_V57$subjsubcort_vol <- pfactor_long_sibs_nosibs_V57$subcort_vol - pfactor_long_sibs_nosibs_V57$site_subcort_vol

pfactor_long_sibs_nosibs_V58$subjmeanct <- pfactor_long_sibs_nosibs_V58$meanwb_ct - pfactor_long_sibs_nosibs_V58$site_meanct
pfactor_long_sibs_nosibs_V58$subjtotalvol <- pfactor_long_sibs_nosibs_V58$wb_cort_vol - pfactor_long_sibs_nosibs_V58$site_totalvol
pfactor_long_sibs_nosibs_V58$subjtotalarea <- pfactor_long_sibs_nosibs_V58$wb_cort_area - pfactor_long_sibs_nosibs_V58$site_totalarea
pfactor_long_sibs_nosibs_V58$subjsubcort_vol <- pfactor_long_sibs_nosibs_V58$subcort_vol - pfactor_long_sibs_nosibs_V58$site_subcort_vol

pfactor_long_sibs_nosibs_V59$subjmeanct <- pfactor_long_sibs_nosibs_V59$meanwb_ct - pfactor_long_sibs_nosibs_V59$site_meanct
pfactor_long_sibs_nosibs_V59$subjtotalvol <- pfactor_long_sibs_nosibs_V59$wb_cort_vol - pfactor_long_sibs_nosibs_V59$site_totalvol
pfactor_long_sibs_nosibs_V59$subjtotalarea <- pfactor_long_sibs_nosibs_V59$wb_cort_area - pfactor_long_sibs_nosibs_V59$site_totalarea
pfactor_long_sibs_nosibs_V59$subjsubcort_vol <- pfactor_long_sibs_nosibs_V59$subcort_vol - pfactor_long_sibs_nosibs_V59$site_subcort_vol

pfactor_long_sibs_nosibs_V60$subjmeanct <- pfactor_long_sibs_nosibs_V60$meanwb_ct - pfactor_long_sibs_nosibs_V60$site_meanct
pfactor_long_sibs_nosibs_V60$subjtotalvol <- pfactor_long_sibs_nosibs_V60$wb_cort_vol - pfactor_long_sibs_nosibs_V60$site_totalvol
pfactor_long_sibs_nosibs_V60$subjtotalarea <- pfactor_long_sibs_nosibs_V60$wb_cort_area - pfactor_long_sibs_nosibs_V60$site_totalarea
pfactor_long_sibs_nosibs_V60$subjsubcort_vol <- pfactor_long_sibs_nosibs_V60$subcort_vol - pfactor_long_sibs_nosibs_V60$site_subcort_vol

pfactor_long_sibs_nosibs_V61$subjmeanct <- pfactor_long_sibs_nosibs_V61$meanwb_ct - pfactor_long_sibs_nosibs_V61$site_meanct
pfactor_long_sibs_nosibs_V61$subjtotalvol <- pfactor_long_sibs_nosibs_V61$wb_cort_vol - pfactor_long_sibs_nosibs_V61$site_totalvol
pfactor_long_sibs_nosibs_V61$subjtotalarea <- pfactor_long_sibs_nosibs_V61$wb_cort_area - pfactor_long_sibs_nosibs_V61$site_totalarea
pfactor_long_sibs_nosibs_V61$subjsubcort_vol <- pfactor_long_sibs_nosibs_V61$subcort_vol - pfactor_long_sibs_nosibs_V61$site_subcort_vol

pfactor_long_sibs_nosibs_V62$subjmeanct <- pfactor_long_sibs_nosibs_V62$meanwb_ct - pfactor_long_sibs_nosibs_V62$site_meanct
pfactor_long_sibs_nosibs_V62$subjtotalvol <- pfactor_long_sibs_nosibs_V62$wb_cort_vol - pfactor_long_sibs_nosibs_V62$site_totalvol
pfactor_long_sibs_nosibs_V62$subjtotalarea <- pfactor_long_sibs_nosibs_V62$wb_cort_area - pfactor_long_sibs_nosibs_V62$site_totalarea
pfactor_long_sibs_nosibs_V62$subjsubcort_vol <- pfactor_long_sibs_nosibs_V62$subcort_vol - pfactor_long_sibs_nosibs_V62$site_subcort_vol

pfactor_long_sibs_nosibs_V63$subjmeanct <- pfactor_long_sibs_nosibs_V63$meanwb_ct - pfactor_long_sibs_nosibs_V63$site_meanct
pfactor_long_sibs_nosibs_V63$subjtotalvol <- pfactor_long_sibs_nosibs_V63$wb_cort_vol - pfactor_long_sibs_nosibs_V63$site_totalvol
pfactor_long_sibs_nosibs_V63$subjtotalarea <- pfactor_long_sibs_nosibs_V63$wb_cort_area - pfactor_long_sibs_nosibs_V63$site_totalarea
pfactor_long_sibs_nosibs_V63$subjsubcort_vol <- pfactor_long_sibs_nosibs_V63$subcort_vol - pfactor_long_sibs_nosibs_V63$site_subcort_vol

pfactor_long_sibs_nosibs_V64$subjmeanct <- pfactor_long_sibs_nosibs_V64$meanwb_ct - pfactor_long_sibs_nosibs_V64$site_meanct
pfactor_long_sibs_nosibs_V64$subjtotalvol <- pfactor_long_sibs_nosibs_V64$wb_cort_vol - pfactor_long_sibs_nosibs_V64$site_totalvol
pfactor_long_sibs_nosibs_V64$subjtotalarea <- pfactor_long_sibs_nosibs_V64$wb_cort_area - pfactor_long_sibs_nosibs_V64$site_totalarea
pfactor_long_sibs_nosibs_V64$subjsubcort_vol <- pfactor_long_sibs_nosibs_V64$subcort_vol - pfactor_long_sibs_nosibs_V64$site_subcort_vol

pfactor_long_sibs_nosibs_V65$subjmeanct <- pfactor_long_sibs_nosibs_V65$meanwb_ct - pfactor_long_sibs_nosibs_V65$site_meanct
pfactor_long_sibs_nosibs_V65$subjtotalvol <- pfactor_long_sibs_nosibs_V65$wb_cort_vol - pfactor_long_sibs_nosibs_V65$site_totalvol
pfactor_long_sibs_nosibs_V65$subjtotalarea <- pfactor_long_sibs_nosibs_V65$wb_cort_area - pfactor_long_sibs_nosibs_V65$site_totalarea
pfactor_long_sibs_nosibs_V65$subjsubcort_vol <- pfactor_long_sibs_nosibs_V65$subcort_vol - pfactor_long_sibs_nosibs_V65$site_subcort_vol

pfactor_long_sibs_nosibs_V66$subjmeanct <- pfactor_long_sibs_nosibs_V66$meanwb_ct - pfactor_long_sibs_nosibs_V66$site_meanct
pfactor_long_sibs_nosibs_V66$subjtotalvol <- pfactor_long_sibs_nosibs_V66$wb_cort_vol - pfactor_long_sibs_nosibs_V66$site_totalvol
pfactor_long_sibs_nosibs_V66$subjtotalarea <- pfactor_long_sibs_nosibs_V66$wb_cort_area - pfactor_long_sibs_nosibs_V66$site_totalarea
pfactor_long_sibs_nosibs_V66$subjsubcort_vol <- pfactor_long_sibs_nosibs_V66$subcort_vol - pfactor_long_sibs_nosibs_V66$site_subcort_vol

pfactor_long_sibs_nosibs_V67$subjmeanct <- pfactor_long_sibs_nosibs_V67$meanwb_ct - pfactor_long_sibs_nosibs_V67$site_meanct
pfactor_long_sibs_nosibs_V67$subjtotalvol <- pfactor_long_sibs_nosibs_V67$wb_cort_vol - pfactor_long_sibs_nosibs_V67$site_totalvol
pfactor_long_sibs_nosibs_V67$subjtotalarea <- pfactor_long_sibs_nosibs_V67$wb_cort_area - pfactor_long_sibs_nosibs_V67$site_totalarea
pfactor_long_sibs_nosibs_V67$subjsubcort_vol <- pfactor_long_sibs_nosibs_V67$subcort_vol - pfactor_long_sibs_nosibs_V67$site_subcort_vol

pfactor_long_sibs_nosibs_V68$subjmeanct <- pfactor_long_sibs_nosibs_V68$meanwb_ct - pfactor_long_sibs_nosibs_V68$site_meanct
pfactor_long_sibs_nosibs_V68$subjtotalvol <- pfactor_long_sibs_nosibs_V68$wb_cort_vol - pfactor_long_sibs_nosibs_V68$site_totalvol
pfactor_long_sibs_nosibs_V68$subjtotalarea <- pfactor_long_sibs_nosibs_V68$wb_cort_area - pfactor_long_sibs_nosibs_V68$site_totalarea
pfactor_long_sibs_nosibs_V68$subjsubcort_vol <- pfactor_long_sibs_nosibs_V68$subcort_vol - pfactor_long_sibs_nosibs_V68$site_subcort_vol

pfactor_long_sibs_nosibs_V69$subjmeanct <- pfactor_long_sibs_nosibs_V69$meanwb_ct - pfactor_long_sibs_nosibs_V69$site_meanct
pfactor_long_sibs_nosibs_V69$subjtotalvol <- pfactor_long_sibs_nosibs_V69$wb_cort_vol - pfactor_long_sibs_nosibs_V69$site_totalvol
pfactor_long_sibs_nosibs_V69$subjtotalarea <- pfactor_long_sibs_nosibs_V69$wb_cort_area - pfactor_long_sibs_nosibs_V69$site_totalarea
pfactor_long_sibs_nosibs_V69$subjsubcort_vol <- pfactor_long_sibs_nosibs_V69$subcort_vol - pfactor_long_sibs_nosibs_V69$site_subcort_vol

pfactor_long_sibs_nosibs_V70$subjmeanct <- pfactor_long_sibs_nosibs_V70$meanwb_ct - pfactor_long_sibs_nosibs_V70$site_meanct
pfactor_long_sibs_nosibs_V70$subjtotalvol <- pfactor_long_sibs_nosibs_V70$wb_cort_vol - pfactor_long_sibs_nosibs_V70$site_totalvol
pfactor_long_sibs_nosibs_V70$subjtotalarea <- pfactor_long_sibs_nosibs_V70$wb_cort_area - pfactor_long_sibs_nosibs_V70$site_totalarea
pfactor_long_sibs_nosibs_V70$subjsubcort_vol <- pfactor_long_sibs_nosibs_V70$subcort_vol - pfactor_long_sibs_nosibs_V70$site_subcort_vol

pfactor_long_sibs_nosibs_V71$subjmeanct <- pfactor_long_sibs_nosibs_V71$meanwb_ct - pfactor_long_sibs_nosibs_V71$site_meanct
pfactor_long_sibs_nosibs_V71$subjtotalvol <- pfactor_long_sibs_nosibs_V71$wb_cort_vol - pfactor_long_sibs_nosibs_V71$site_totalvol
pfactor_long_sibs_nosibs_V71$subjtotalarea <- pfactor_long_sibs_nosibs_V71$wb_cort_area - pfactor_long_sibs_nosibs_V71$site_totalarea
pfactor_long_sibs_nosibs_V71$subjsubcort_vol <- pfactor_long_sibs_nosibs_V71$subcort_vol - pfactor_long_sibs_nosibs_V71$site_subcort_vol

pfactor_long_sibs_nosibs_V72$subjmeanct <- pfactor_long_sibs_nosibs_V72$meanwb_ct - pfactor_long_sibs_nosibs_V72$site_meanct
pfactor_long_sibs_nosibs_V72$subjtotalvol <- pfactor_long_sibs_nosibs_V72$wb_cort_vol - pfactor_long_sibs_nosibs_V72$site_totalvol
pfactor_long_sibs_nosibs_V72$subjtotalarea <- pfactor_long_sibs_nosibs_V72$wb_cort_area - pfactor_long_sibs_nosibs_V72$site_totalarea
pfactor_long_sibs_nosibs_V72$subjsubcort_vol <- pfactor_long_sibs_nosibs_V72$subcort_vol - pfactor_long_sibs_nosibs_V72$site_subcort_vol

pfactor_long_sibs_nosibs_V73$subjmeanct <- pfactor_long_sibs_nosibs_V73$meanwb_ct - pfactor_long_sibs_nosibs_V73$site_meanct
pfactor_long_sibs_nosibs_V73$subjtotalvol <- pfactor_long_sibs_nosibs_V73$wb_cort_vol - pfactor_long_sibs_nosibs_V73$site_totalvol
pfactor_long_sibs_nosibs_V73$subjtotalarea <- pfactor_long_sibs_nosibs_V73$wb_cort_area - pfactor_long_sibs_nosibs_V73$site_totalarea
pfactor_long_sibs_nosibs_V73$subjsubcort_vol <- pfactor_long_sibs_nosibs_V73$subcort_vol - pfactor_long_sibs_nosibs_V73$site_subcort_vol

pfactor_long_sibs_nosibs_V74$subjmeanct <- pfactor_long_sibs_nosibs_V74$meanwb_ct - pfactor_long_sibs_nosibs_V74$site_meanct
pfactor_long_sibs_nosibs_V74$subjtotalvol <- pfactor_long_sibs_nosibs_V74$wb_cort_vol - pfactor_long_sibs_nosibs_V74$site_totalvol
pfactor_long_sibs_nosibs_V74$subjtotalarea <- pfactor_long_sibs_nosibs_V74$wb_cort_area - pfactor_long_sibs_nosibs_V74$site_totalarea
pfactor_long_sibs_nosibs_V74$subjsubcort_vol <- pfactor_long_sibs_nosibs_V74$subcort_vol - pfactor_long_sibs_nosibs_V74$site_subcort_vol

pfactor_long_sibs_nosibs_V75$subjmeanct <- pfactor_long_sibs_nosibs_V75$meanwb_ct - pfactor_long_sibs_nosibs_V75$site_meanct
pfactor_long_sibs_nosibs_V75$subjtotalvol <- pfactor_long_sibs_nosibs_V75$wb_cort_vol - pfactor_long_sibs_nosibs_V75$site_totalvol
pfactor_long_sibs_nosibs_V75$subjtotalarea <- pfactor_long_sibs_nosibs_V75$wb_cort_area - pfactor_long_sibs_nosibs_V75$site_totalarea
pfactor_long_sibs_nosibs_V75$subjsubcort_vol <- pfactor_long_sibs_nosibs_V75$subcort_vol - pfactor_long_sibs_nosibs_V75$site_subcort_vol

pfactor_long_sibs_nosibs_V76$subjmeanct <- pfactor_long_sibs_nosibs_V76$meanwb_ct - pfactor_long_sibs_nosibs_V76$site_meanct
pfactor_long_sibs_nosibs_V76$subjtotalvol <- pfactor_long_sibs_nosibs_V76$wb_cort_vol - pfactor_long_sibs_nosibs_V76$site_totalvol
pfactor_long_sibs_nosibs_V76$subjtotalarea <- pfactor_long_sibs_nosibs_V76$wb_cort_area - pfactor_long_sibs_nosibs_V76$site_totalarea
pfactor_long_sibs_nosibs_V76$subjsubcort_vol <- pfactor_long_sibs_nosibs_V76$subcort_vol - pfactor_long_sibs_nosibs_V76$site_subcort_vol

pfactor_long_sibs_nosibs_V77$subjmeanct <- pfactor_long_sibs_nosibs_V77$meanwb_ct - pfactor_long_sibs_nosibs_V77$site_meanct
pfactor_long_sibs_nosibs_V77$subjtotalvol <- pfactor_long_sibs_nosibs_V77$wb_cort_vol - pfactor_long_sibs_nosibs_V77$site_totalvol
pfactor_long_sibs_nosibs_V77$subjtotalarea <- pfactor_long_sibs_nosibs_V77$wb_cort_area - pfactor_long_sibs_nosibs_V77$site_totalarea
pfactor_long_sibs_nosibs_V77$subjsubcort_vol <- pfactor_long_sibs_nosibs_V77$subcort_vol - pfactor_long_sibs_nosibs_V77$site_subcort_vol

pfactor_long_sibs_nosibs_V78$subjmeanct <- pfactor_long_sibs_nosibs_V78$meanwb_ct - pfactor_long_sibs_nosibs_V78$site_meanct
pfactor_long_sibs_nosibs_V78$subjtotalvol <- pfactor_long_sibs_nosibs_V78$wb_cort_vol - pfactor_long_sibs_nosibs_V78$site_totalvol
pfactor_long_sibs_nosibs_V78$subjtotalarea <- pfactor_long_sibs_nosibs_V78$wb_cort_area - pfactor_long_sibs_nosibs_V78$site_totalarea
pfactor_long_sibs_nosibs_V78$subjsubcort_vol <- pfactor_long_sibs_nosibs_V78$subcort_vol - pfactor_long_sibs_nosibs_V78$site_subcort_vol

pfactor_long_sibs_nosibs_V79$subjmeanct <- pfactor_long_sibs_nosibs_V79$meanwb_ct - pfactor_long_sibs_nosibs_V79$site_meanct
pfactor_long_sibs_nosibs_V79$subjtotalvol <- pfactor_long_sibs_nosibs_V79$wb_cort_vol - pfactor_long_sibs_nosibs_V79$site_totalvol
pfactor_long_sibs_nosibs_V79$subjtotalarea <- pfactor_long_sibs_nosibs_V79$wb_cort_area - pfactor_long_sibs_nosibs_V79$site_totalarea
pfactor_long_sibs_nosibs_V79$subjsubcort_vol <- pfactor_long_sibs_nosibs_V79$subcort_vol - pfactor_long_sibs_nosibs_V79$site_subcort_vol

pfactor_long_sibs_nosibs_V80$subjmeanct <- pfactor_long_sibs_nosibs_V80$meanwb_ct - pfactor_long_sibs_nosibs_V80$site_meanct
pfactor_long_sibs_nosibs_V80$subjtotalvol <- pfactor_long_sibs_nosibs_V80$wb_cort_vol - pfactor_long_sibs_nosibs_V80$site_totalvol
pfactor_long_sibs_nosibs_V80$subjtotalarea <- pfactor_long_sibs_nosibs_V80$wb_cort_area - pfactor_long_sibs_nosibs_V80$site_totalarea
pfactor_long_sibs_nosibs_V80$subjsubcort_vol <- pfactor_long_sibs_nosibs_V80$subcort_vol - pfactor_long_sibs_nosibs_V80$site_subcort_vol

pfactor_long_sibs_nosibs_V81$subjmeanct <- pfactor_long_sibs_nosibs_V81$meanwb_ct - pfactor_long_sibs_nosibs_V81$site_meanct
pfactor_long_sibs_nosibs_V81$subjtotalvol <- pfactor_long_sibs_nosibs_V81$wb_cort_vol - pfactor_long_sibs_nosibs_V81$site_totalvol
pfactor_long_sibs_nosibs_V81$subjtotalarea <- pfactor_long_sibs_nosibs_V81$wb_cort_area - pfactor_long_sibs_nosibs_V81$site_totalarea
pfactor_long_sibs_nosibs_V81$subjsubcort_vol <- pfactor_long_sibs_nosibs_V81$subcort_vol - pfactor_long_sibs_nosibs_V81$site_subcort_vol

pfactor_long_sibs_nosibs_V82$subjmeanct <- pfactor_long_sibs_nosibs_V82$meanwb_ct - pfactor_long_sibs_nosibs_V82$site_meanct
pfactor_long_sibs_nosibs_V82$subjtotalvol <- pfactor_long_sibs_nosibs_V82$wb_cort_vol - pfactor_long_sibs_nosibs_V82$site_totalvol
pfactor_long_sibs_nosibs_V82$subjtotalarea <- pfactor_long_sibs_nosibs_V82$wb_cort_area - pfactor_long_sibs_nosibs_V82$site_totalarea
pfactor_long_sibs_nosibs_V82$subjsubcort_vol <- pfactor_long_sibs_nosibs_V82$subcort_vol - pfactor_long_sibs_nosibs_V82$site_subcort_vol

pfactor_long_sibs_nosibs_V83$subjmeanct <- pfactor_long_sibs_nosibs_V83$meanwb_ct - pfactor_long_sibs_nosibs_V83$site_meanct
pfactor_long_sibs_nosibs_V83$subjtotalvol <- pfactor_long_sibs_nosibs_V83$wb_cort_vol - pfactor_long_sibs_nosibs_V83$site_totalvol
pfactor_long_sibs_nosibs_V83$subjtotalarea <- pfactor_long_sibs_nosibs_V83$wb_cort_area - pfactor_long_sibs_nosibs_V83$site_totalarea
pfactor_long_sibs_nosibs_V83$subjsubcort_vol <- pfactor_long_sibs_nosibs_V83$subcort_vol - pfactor_long_sibs_nosibs_V83$site_subcort_vol

pfactor_long_sibs_nosibs_V84$subjmeanct <- pfactor_long_sibs_nosibs_V84$meanwb_ct - pfactor_long_sibs_nosibs_V84$site_meanct
pfactor_long_sibs_nosibs_V84$subjtotalvol <- pfactor_long_sibs_nosibs_V84$wb_cort_vol - pfactor_long_sibs_nosibs_V84$site_totalvol
pfactor_long_sibs_nosibs_V84$subjtotalarea <- pfactor_long_sibs_nosibs_V84$wb_cort_area - pfactor_long_sibs_nosibs_V84$site_totalarea
pfactor_long_sibs_nosibs_V84$subjsubcort_vol <- pfactor_long_sibs_nosibs_V84$subcort_vol - pfactor_long_sibs_nosibs_V84$site_subcort_vol

pfactor_long_sibs_nosibs_V85$subjmeanct <- pfactor_long_sibs_nosibs_V85$meanwb_ct - pfactor_long_sibs_nosibs_V85$site_meanct
pfactor_long_sibs_nosibs_V85$subjtotalvol <- pfactor_long_sibs_nosibs_V85$wb_cort_vol - pfactor_long_sibs_nosibs_V85$site_totalvol
pfactor_long_sibs_nosibs_V85$subjtotalarea <- pfactor_long_sibs_nosibs_V85$wb_cort_area - pfactor_long_sibs_nosibs_V85$site_totalarea
pfactor_long_sibs_nosibs_V85$subjsubcort_vol <- pfactor_long_sibs_nosibs_V85$subcort_vol - pfactor_long_sibs_nosibs_V85$site_subcort_vol

pfactor_long_sibs_nosibs_V86$subjmeanct <- pfactor_long_sibs_nosibs_V86$meanwb_ct - pfactor_long_sibs_nosibs_V86$site_meanct
pfactor_long_sibs_nosibs_V86$subjtotalvol <- pfactor_long_sibs_nosibs_V86$wb_cort_vol - pfactor_long_sibs_nosibs_V86$site_totalvol
pfactor_long_sibs_nosibs_V86$subjtotalarea <- pfactor_long_sibs_nosibs_V86$wb_cort_area - pfactor_long_sibs_nosibs_V86$site_totalarea
pfactor_long_sibs_nosibs_V86$subjsubcort_vol <- pfactor_long_sibs_nosibs_V86$subcort_vol - pfactor_long_sibs_nosibs_V86$site_subcort_vol

pfactor_long_sibs_nosibs_V87$subjmeanct <- pfactor_long_sibs_nosibs_V87$meanwb_ct - pfactor_long_sibs_nosibs_V87$site_meanct
pfactor_long_sibs_nosibs_V87$subjtotalvol <- pfactor_long_sibs_nosibs_V87$wb_cort_vol - pfactor_long_sibs_nosibs_V87$site_totalvol
pfactor_long_sibs_nosibs_V87$subjtotalarea <- pfactor_long_sibs_nosibs_V87$wb_cort_area - pfactor_long_sibs_nosibs_V87$site_totalarea
pfactor_long_sibs_nosibs_V87$subjsubcort_vol <- pfactor_long_sibs_nosibs_V87$subcort_vol - pfactor_long_sibs_nosibs_V87$site_subcort_vol

pfactor_long_sibs_nosibs_V88$subjmeanct <- pfactor_long_sibs_nosibs_V88$meanwb_ct - pfactor_long_sibs_nosibs_V88$site_meanct
pfactor_long_sibs_nosibs_V88$subjtotalvol <- pfactor_long_sibs_nosibs_V88$wb_cort_vol - pfactor_long_sibs_nosibs_V88$site_totalvol
pfactor_long_sibs_nosibs_V88$subjtotalarea <- pfactor_long_sibs_nosibs_V88$wb_cort_area - pfactor_long_sibs_nosibs_V88$site_totalarea
pfactor_long_sibs_nosibs_V88$subjsubcort_vol <- pfactor_long_sibs_nosibs_V88$subcort_vol - pfactor_long_sibs_nosibs_V88$site_subcort_vol

pfactor_long_sibs_nosibs_V89$subjmeanct <- pfactor_long_sibs_nosibs_V89$meanwb_ct - pfactor_long_sibs_nosibs_V89$site_meanct
pfactor_long_sibs_nosibs_V89$subjtotalvol <- pfactor_long_sibs_nosibs_V89$wb_cort_vol - pfactor_long_sibs_nosibs_V89$site_totalvol
pfactor_long_sibs_nosibs_V89$subjtotalarea <- pfactor_long_sibs_nosibs_V89$wb_cort_area - pfactor_long_sibs_nosibs_V89$site_totalarea
pfactor_long_sibs_nosibs_V89$subjsubcort_vol <- pfactor_long_sibs_nosibs_V89$subcort_vol - pfactor_long_sibs_nosibs_V89$site_subcort_vol

pfactor_long_sibs_nosibs_V90$subjmeanct <- pfactor_long_sibs_nosibs_V90$meanwb_ct - pfactor_long_sibs_nosibs_V90$site_meanct
pfactor_long_sibs_nosibs_V90$subjtotalvol <- pfactor_long_sibs_nosibs_V90$wb_cort_vol - pfactor_long_sibs_nosibs_V90$site_totalvol
pfactor_long_sibs_nosibs_V90$subjtotalarea <- pfactor_long_sibs_nosibs_V90$wb_cort_area - pfactor_long_sibs_nosibs_V90$site_totalarea
pfactor_long_sibs_nosibs_V90$subjsubcort_vol <- pfactor_long_sibs_nosibs_V90$subcort_vol - pfactor_long_sibs_nosibs_V90$site_subcort_vol

pfactor_long_sibs_nosibs_V91$subjmeanct <- pfactor_long_sibs_nosibs_V91$meanwb_ct - pfactor_long_sibs_nosibs_V91$site_meanct
pfactor_long_sibs_nosibs_V91$subjtotalvol <- pfactor_long_sibs_nosibs_V91$wb_cort_vol - pfactor_long_sibs_nosibs_V91$site_totalvol
pfactor_long_sibs_nosibs_V91$subjtotalarea <- pfactor_long_sibs_nosibs_V91$wb_cort_area - pfactor_long_sibs_nosibs_V91$site_totalarea
pfactor_long_sibs_nosibs_V91$subjsubcort_vol <- pfactor_long_sibs_nosibs_V91$subcort_vol - pfactor_long_sibs_nosibs_V91$site_subcort_vol

pfactor_long_sibs_nosibs_V92$subjmeanct <- pfactor_long_sibs_nosibs_V92$meanwb_ct - pfactor_long_sibs_nosibs_V92$site_meanct
pfactor_long_sibs_nosibs_V92$subjtotalvol <- pfactor_long_sibs_nosibs_V92$wb_cort_vol - pfactor_long_sibs_nosibs_V92$site_totalvol
pfactor_long_sibs_nosibs_V92$subjtotalarea <- pfactor_long_sibs_nosibs_V92$wb_cort_area - pfactor_long_sibs_nosibs_V92$site_totalarea
pfactor_long_sibs_nosibs_V92$subjsubcort_vol <- pfactor_long_sibs_nosibs_V92$subcort_vol - pfactor_long_sibs_nosibs_V92$site_subcort_vol

pfactor_long_sibs_nosibs_V93$subjmeanct <- pfactor_long_sibs_nosibs_V93$meanwb_ct - pfactor_long_sibs_nosibs_V93$site_meanct
pfactor_long_sibs_nosibs_V93$subjtotalvol <- pfactor_long_sibs_nosibs_V93$wb_cort_vol - pfactor_long_sibs_nosibs_V93$site_totalvol
pfactor_long_sibs_nosibs_V93$subjtotalarea <- pfactor_long_sibs_nosibs_V93$wb_cort_area - pfactor_long_sibs_nosibs_V93$site_totalarea
pfactor_long_sibs_nosibs_V93$subjsubcort_vol <- pfactor_long_sibs_nosibs_V93$subcort_vol - pfactor_long_sibs_nosibs_V93$site_subcort_vol

pfactor_long_sibs_nosibs_V94$subjmeanct <- pfactor_long_sibs_nosibs_V94$meanwb_ct - pfactor_long_sibs_nosibs_V94$site_meanct
pfactor_long_sibs_nosibs_V94$subjtotalvol <- pfactor_long_sibs_nosibs_V94$wb_cort_vol - pfactor_long_sibs_nosibs_V94$site_totalvol
pfactor_long_sibs_nosibs_V94$subjtotalarea <- pfactor_long_sibs_nosibs_V94$wb_cort_area - pfactor_long_sibs_nosibs_V94$site_totalarea
pfactor_long_sibs_nosibs_V94$subjsubcort_vol <- pfactor_long_sibs_nosibs_V94$subcort_vol - pfactor_long_sibs_nosibs_V94$site_subcort_vol

pfactor_long_sibs_nosibs_V95$subjmeanct <- pfactor_long_sibs_nosibs_V95$meanwb_ct - pfactor_long_sibs_nosibs_V95$site_meanct
pfactor_long_sibs_nosibs_V95$subjtotalvol <- pfactor_long_sibs_nosibs_V95$wb_cort_vol - pfactor_long_sibs_nosibs_V95$site_totalvol
pfactor_long_sibs_nosibs_V95$subjtotalarea <- pfactor_long_sibs_nosibs_V95$wb_cort_area - pfactor_long_sibs_nosibs_V95$site_totalarea
pfactor_long_sibs_nosibs_V95$subjsubcort_vol <- pfactor_long_sibs_nosibs_V95$subcort_vol - pfactor_long_sibs_nosibs_V95$site_subcort_vol

pfactor_long_sibs_nosibs_V96$subjmeanct <- pfactor_long_sibs_nosibs_V96$meanwb_ct - pfactor_long_sibs_nosibs_V96$site_meanct
pfactor_long_sibs_nosibs_V96$subjtotalvol <- pfactor_long_sibs_nosibs_V96$wb_cort_vol - pfactor_long_sibs_nosibs_V96$site_totalvol
pfactor_long_sibs_nosibs_V96$subjtotalarea <- pfactor_long_sibs_nosibs_V96$wb_cort_area - pfactor_long_sibs_nosibs_V96$site_totalarea
pfactor_long_sibs_nosibs_V96$subjsubcort_vol <- pfactor_long_sibs_nosibs_V96$subcort_vol - pfactor_long_sibs_nosibs_V96$site_subcort_vol

pfactor_long_sibs_nosibs_V97$subjmeanct <- pfactor_long_sibs_nosibs_V97$meanwb_ct - pfactor_long_sibs_nosibs_V97$site_meanct
pfactor_long_sibs_nosibs_V97$subjtotalvol <- pfactor_long_sibs_nosibs_V97$wb_cort_vol - pfactor_long_sibs_nosibs_V97$site_totalvol
pfactor_long_sibs_nosibs_V97$subjtotalarea <- pfactor_long_sibs_nosibs_V97$wb_cort_area - pfactor_long_sibs_nosibs_V97$site_totalarea
pfactor_long_sibs_nosibs_V97$subjsubcort_vol <- pfactor_long_sibs_nosibs_V97$subcort_vol - pfactor_long_sibs_nosibs_V97$site_subcort_vol

pfactor_long_sibs_nosibs_V98$subjmeanct <- pfactor_long_sibs_nosibs_V98$meanwb_ct - pfactor_long_sibs_nosibs_V98$site_meanct
pfactor_long_sibs_nosibs_V98$subjtotalvol <- pfactor_long_sibs_nosibs_V98$wb_cort_vol - pfactor_long_sibs_nosibs_V98$site_totalvol
pfactor_long_sibs_nosibs_V98$subjtotalarea <- pfactor_long_sibs_nosibs_V98$wb_cort_area - pfactor_long_sibs_nosibs_V98$site_totalarea
pfactor_long_sibs_nosibs_V98$subjsubcort_vol <- pfactor_long_sibs_nosibs_V98$subcort_vol - pfactor_long_sibs_nosibs_V98$site_subcort_vol

pfactor_long_sibs_nosibs_V99$subjmeanct <- pfactor_long_sibs_nosibs_V99$meanwb_ct - pfactor_long_sibs_nosibs_V99$site_meanct
pfactor_long_sibs_nosibs_V99$subjtotalvol <- pfactor_long_sibs_nosibs_V99$wb_cort_vol - pfactor_long_sibs_nosibs_V99$site_totalvol
pfactor_long_sibs_nosibs_V99$subjtotalarea <- pfactor_long_sibs_nosibs_V99$wb_cort_area - pfactor_long_sibs_nosibs_V99$site_totalarea
pfactor_long_sibs_nosibs_V99$subjsubcort_vol <- pfactor_long_sibs_nosibs_V99$subcort_vol - pfactor_long_sibs_nosibs_V99$site_subcort_vol

pfactor_long_sibs_nosibs_V100$subjmeanct <- pfactor_long_sibs_nosibs_V100$meanwb_ct - pfactor_long_sibs_nosibs_V100$site_meanct
pfactor_long_sibs_nosibs_V100$subjtotalvol <- pfactor_long_sibs_nosibs_V100$wb_cort_vol - pfactor_long_sibs_nosibs_V100$site_totalvol
pfactor_long_sibs_nosibs_V100$subjtotalarea <- pfactor_long_sibs_nosibs_V100$wb_cort_area - pfactor_long_sibs_nosibs_V100$site_totalarea
pfactor_long_sibs_nosibs_V100$subjsubcort_vol <- pfactor_long_sibs_nosibs_V100$subcort_vol - pfactor_long_sibs_nosibs_V100$site_subcort_vol

#loop of analysis run 100 times with each of the random samples for each of the psychopathology factor scores

#list of data frames to use in loop below
df <- list(pfactor_long_sibs_nosibs_V1,pfactor_long_sibs_nosibs_V2,pfactor_long_sibs_nosibs_V3,pfactor_long_sibs_nosibs_V4,
           pfactor_long_sibs_nosibs_V5,pfactor_long_sibs_nosibs_V6,pfactor_long_sibs_nosibs_V7,pfactor_long_sibs_nosibs_V8,
           pfactor_long_sibs_nosibs_V9,pfactor_long_sibs_nosibs_V10,pfactor_long_sibs_nosibs_V11,pfactor_long_sibs_nosibs_V12,
           pfactor_long_sibs_nosibs_V13,pfactor_long_sibs_nosibs_V14,pfactor_long_sibs_nosibs_V15,pfactor_long_sibs_nosibs_V16,
           pfactor_long_sibs_nosibs_V17,pfactor_long_sibs_nosibs_V18,pfactor_long_sibs_nosibs_V19,pfactor_long_sibs_nosibs_V20,
           pfactor_long_sibs_nosibs_V21,pfactor_long_sibs_nosibs_V22,pfactor_long_sibs_nosibs_V23,pfactor_long_sibs_nosibs_V24,
           pfactor_long_sibs_nosibs_V25,pfactor_long_sibs_nosibs_V26,pfactor_long_sibs_nosibs_V27,pfactor_long_sibs_nosibs_V28,
           pfactor_long_sibs_nosibs_V29,pfactor_long_sibs_nosibs_V30,pfactor_long_sibs_nosibs_V31,pfactor_long_sibs_nosibs_V32,
           pfactor_long_sibs_nosibs_V33,pfactor_long_sibs_nosibs_V34,pfactor_long_sibs_nosibs_V35,pfactor_long_sibs_nosibs_V36,
           pfactor_long_sibs_nosibs_V37,pfactor_long_sibs_nosibs_V38,pfactor_long_sibs_nosibs_V39,pfactor_long_sibs_nosibs_V40,
           pfactor_long_sibs_nosibs_V41,pfactor_long_sibs_nosibs_V42,pfactor_long_sibs_nosibs_V43,pfactor_long_sibs_nosibs_V44,
           pfactor_long_sibs_nosibs_V45,pfactor_long_sibs_nosibs_V46,pfactor_long_sibs_nosibs_V47,pfactor_long_sibs_nosibs_V48,
           pfactor_long_sibs_nosibs_V49,pfactor_long_sibs_nosibs_V50,pfactor_long_sibs_nosibs_V51,pfactor_long_sibs_nosibs_V52,
           pfactor_long_sibs_nosibs_V53,pfactor_long_sibs_nosibs_V54,pfactor_long_sibs_nosibs_V55,pfactor_long_sibs_nosibs_V56,
           pfactor_long_sibs_nosibs_V57,pfactor_long_sibs_nosibs_V58,pfactor_long_sibs_nosibs_V59,pfactor_long_sibs_nosibs_V60,
           pfactor_long_sibs_nosibs_V61,pfactor_long_sibs_nosibs_V62,pfactor_long_sibs_nosibs_V63,pfactor_long_sibs_nosibs_V64,
           pfactor_long_sibs_nosibs_V65,pfactor_long_sibs_nosibs_V66,pfactor_long_sibs_nosibs_V67,pfactor_long_sibs_nosibs_V68,
           pfactor_long_sibs_nosibs_V69,pfactor_long_sibs_nosibs_V70,pfactor_long_sibs_nosibs_V71,pfactor_long_sibs_nosibs_V72,
           pfactor_long_sibs_nosibs_V73,pfactor_long_sibs_nosibs_V74,pfactor_long_sibs_nosibs_V75,pfactor_long_sibs_nosibs_V76,
           pfactor_long_sibs_nosibs_V77,pfactor_long_sibs_nosibs_V78,pfactor_long_sibs_nosibs_V79,pfactor_long_sibs_nosibs_V80,
           pfactor_long_sibs_nosibs_V81,pfactor_long_sibs_nosibs_V82,pfactor_long_sibs_nosibs_V83,pfactor_long_sibs_nosibs_V84,
           pfactor_long_sibs_nosibs_V85,pfactor_long_sibs_nosibs_V86,pfactor_long_sibs_nosibs_V87,pfactor_long_sibs_nosibs_V88,
           pfactor_long_sibs_nosibs_V89,pfactor_long_sibs_nosibs_V90,pfactor_long_sibs_nosibs_V91,pfactor_long_sibs_nosibs_V92,
           pfactor_long_sibs_nosibs_V93,pfactor_long_sibs_nosibs_V94,pfactor_long_sibs_nosibs_V95,pfactor_long_sibs_nosibs_V96,
           pfactor_long_sibs_nosibs_V97,pfactor_long_sibs_nosibs_V98,pfactor_long_sibs_nosibs_V99,pfactor_long_sibs_nosibs_V100)


#conditional three-level growth model with global brain structure measures predicting the p factor scores from higher-order model
#looped across 100 dataframes of different random siblings included
#covariates: sex, age, scanner model dummies
#standardized betas  and p-values for intercepts (i.e., intb, intp) and slopes (i.e., slb, slp) were saved in a csv file
wholebrainst_p <- data.frame(x=c("CT_int","CT_sl","SA_int","SA_sl","cortvol_int","cortvol_sl","subcortvol_int","subcortvol_sl"))
for (i in 1:length(df)){ 
  CT_intb <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                            (1 + wave|site_id/id), #random intercept and slope for subject and site 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  CT_slb <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  CT_intp <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  CT_slp <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  SA_intb <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  SA_slb <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  SA_intp <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) +
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  SA_slp <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  cortvol_intb <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  cortvol_slb <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  cortvol_intp <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  cortvol_slp <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  subcortvol_intb <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  subcortvol_slb <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) +
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  subcortvol_intp <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  subcortvol_slp <- summary(lmer(scale(p_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  for (j in 2:101) #dataframe names
    newb <- c(CT_intb,CT_slb,SA_intb,SA_slb,cortvol_intb,cortvol_slb,subcortvol_intb,subcortvol_slb)
  wholebrainst_p <- cbind(wholebrainst_p,newb)
  names(wholebrainst_p)[ncol(wholebrainst_p)] <- names(sibs)[j]
  newp <- c(CT_intp,CT_slp,SA_intp,SA_slp,cortvol_intp,cortvol_slp,subcortvol_intp,subcortvol_slp)
  wholebrainst_p <- cbind(wholebrainst_p,newp)
  names(wholebrainst_p)[ncol(wholebrainst_p)] <- names(sibs)[j]
  
}
wholebrainst_p

names(wholebrainst_p) <- c("x","V1_beta","V1_pval","V2_beta","V2_pval","V3_beta","V3_pval","V4_beta","V4_pval",
                             "V5_beta","V5_pval","V6_beta","V6_pval","V7_beta","V7_pval","V8_beta","V8_pval",
                             "V9_beta","V9_pval","V10_beta","V10_pval","V11_beta","V11_pval","V12_beta","V12_pval",
                             "V13_beta","V13_pval","V14_beta","V14_pval","V15_beta","V15_pval","V16_beta","V16_pval",
                             "V17_beta","V17_pval","V18_beta","V18_pval","V19_beta","V19_pval","V20_beta","V20_pval",
                             "V21_beta","V21_pval","V22_beta","V22_pval","V23_beta","V23_pval","V24_beta","V24_pval",
                             "V25_beta","V25_pval","V26_beta","V26_pval","V27_beta","V27_pval","V28_beta","V28_pval",
                             "V29_beta","V29_pval","V30_beta","V30_pval","V31_beta","V31_pval","V32_beta","V32_pval",
                             "V33_beta","V33_pval","V34_beta","V34_pval","V35_beta","V35_pval","V36_beta","V36_pval",
                             "V37_beta","V37_pval","V38_beta","V38_pval","V39_beta","V39_pval","V40_beta","V40_pval",
                             "V41_beta","V41_pval","V42_beta","V42_pval","V43_beta","V43_pval","V44_beta","V44_pval",
                             "V45_beta","V45_pval","V46_beta","V46_pval","V47_beta","V47_pval","V48_beta","V48_pval",
                             "V49_beta","V49_pval","V50_beta","V50_pval","V51_beta","V51_pval","V52_beta","V52_pval",
                             "V53_beta","V53_pval","V54_beta","V54_pval","V55_beta","V55_pval","V56_beta","V56_pval",
                             "V57_beta","V57_pval","V58_beta","V58_pval","V59_beta","V59_pval","V60_beta","V60_pval",
                             "V61_beta","V61_pval","V62_beta","V62_pval","V63_beta","V63_pval","V64_beta","V64_pval",
                             "V65_beta","V65_pval","V66_beta","V66_pval","V67_beta","V67_pval","V68_beta","V68_pval",
                             "V69_beta","V69_pval","V70_beta","V70_pval","V71_beta","V71_pval","V72_beta","V72_pval",
                             "V73_beta","V73_pval","V74_beta","V74_pval","V75_beta","V75_pval","V76_beta","V76_pval",
                             "V77_beta","V77_pval","V78_beta","V78_pval","V79_beta","V79_pval","V80_beta","V80_pval",
                             "V81_beta","V81_pval","V82_beta","V82_pval","V83_beta","V83_pval","V84_beta","V84_pval",
                             "V85_beta","V85_pval","V86_beta","V86_pval","V87_beta","V87_pval","V88_beta","V88_pval",
                             "V89_beta","V89_pval","V90_beta","V90_pval","V91_beta","V91_pval","V92_beta","V92_pval",
                             "V93_beta","V93_pval","V94_beta","V94_pval","V95_beta","V95_pval","V96_beta","V96_pval",
                             "V97_beta","V97_pval","V98_beta","V98_pval","V99_beta","V99_pval","V100_beta","V100_pval")

write.csv(wholebrainst_p, "100_random_samples_output_p_final.csv") #write csv file

#conditional three-level growth model with global brain structure measures predicting the EXT factor scores from correlated factors model
#looped across 100 dataframes of different random siblings included
#covariates: sex, age, scanner model dummies
#standardized betas  and p-values for intercepts (i.e., intb, intp) and slopes (i.e., slb, slp) were saved in a csv file
wholebrainst_ext <- data.frame(x=c("CT_int","CT_sl","SA_int","SA_sl","cortvol_int","cortvol_sl","subcortvol_int","subcortvol_sl"))
for (i in 1:length(df)){ 
  CT_intb <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                            (1 + wave|site_id/id), #random intercept and slope for subject and site 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  CT_slb <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  CT_intp <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  CT_slp <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  SA_intb <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  SA_slb <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  SA_intp <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) +
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  SA_slp <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  cortvol_intb <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  cortvol_slb <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  cortvol_intp <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  cortvol_slp <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  subcortvol_intb <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  subcortvol_slb <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) +
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  subcortvol_intp <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  subcortvol_slp <- summary(lmer(scale(ext_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  for (j in 2:101) #dataframe names
    newb <- c(CT_intb,CT_slb,SA_intb,SA_slb,cortvol_intb,cortvol_slb,subcortvol_intb,subcortvol_slb)
  wholebrainst_ext <- cbind(wholebrainst_ext,newb)
  names(wholebrainst_ext)[ncol(wholebrainst_ext)] <- names(sibs)[j]
  newp <- c(CT_intp,CT_slp,SA_intp,SA_slp,cortvol_intp,cortvol_slp,subcortvol_intp,subcortvol_slp)
  wholebrainst_ext <- cbind(wholebrainst_ext,newp)
  names(wholebrainst_ext)[ncol(wholebrainst_ext)] <- names(sibs)[j]
  
}
wholebrainst_ext

names(wholebrainst_ext) <- c("x","V1_beta","V1_pval","V2_beta","V2_pval","V3_beta","V3_pval","V4_beta","V4_pval",
                             "V5_beta","V5_pval","V6_beta","V6_pval","V7_beta","V7_pval","V8_beta","V8_pval",
                             "V9_beta","V9_pval","V10_beta","V10_pval","V11_beta","V11_pval","V12_beta","V12_pval",
                             "V13_beta","V13_pval","V14_beta","V14_pval","V15_beta","V15_pval","V16_beta","V16_pval",
                             "V17_beta","V17_pval","V18_beta","V18_pval","V19_beta","V19_pval","V20_beta","V20_pval",
                             "V21_beta","V21_pval","V22_beta","V22_pval","V23_beta","V23_pval","V24_beta","V24_pval",
                             "V25_beta","V25_pval","V26_beta","V26_pval","V27_beta","V27_pval","V28_beta","V28_pval",
                             "V29_beta","V29_pval","V30_beta","V30_pval","V31_beta","V31_pval","V32_beta","V32_pval",
                             "V33_beta","V33_pval","V34_beta","V34_pval","V35_beta","V35_pval","V36_beta","V36_pval",
                             "V37_beta","V37_pval","V38_beta","V38_pval","V39_beta","V39_pval","V40_beta","V40_pval",
                             "V41_beta","V41_pval","V42_beta","V42_pval","V43_beta","V43_pval","V44_beta","V44_pval",
                             "V45_beta","V45_pval","V46_beta","V46_pval","V47_beta","V47_pval","V48_beta","V48_pval",
                             "V49_beta","V49_pval","V50_beta","V50_pval","V51_beta","V51_pval","V52_beta","V52_pval",
                             "V53_beta","V53_pval","V54_beta","V54_pval","V55_beta","V55_pval","V56_beta","V56_pval",
                             "V57_beta","V57_pval","V58_beta","V58_pval","V59_beta","V59_pval","V60_beta","V60_pval",
                             "V61_beta","V61_pval","V62_beta","V62_pval","V63_beta","V63_pval","V64_beta","V64_pval",
                             "V65_beta","V65_pval","V66_beta","V66_pval","V67_beta","V67_pval","V68_beta","V68_pval",
                             "V69_beta","V69_pval","V70_beta","V70_pval","V71_beta","V71_pval","V72_beta","V72_pval",
                             "V73_beta","V73_pval","V74_beta","V74_pval","V75_beta","V75_pval","V76_beta","V76_pval",
                             "V77_beta","V77_pval","V78_beta","V78_pval","V79_beta","V79_pval","V80_beta","V80_pval",
                             "V81_beta","V81_pval","V82_beta","V82_pval","V83_beta","V83_pval","V84_beta","V84_pval",
                             "V85_beta","V85_pval","V86_beta","V86_pval","V87_beta","V87_pval","V88_beta","V88_pval",
                             "V89_beta","V89_pval","V90_beta","V90_pval","V91_beta","V91_pval","V92_beta","V92_pval",
                             "V93_beta","V93_pval","V94_beta","V94_pval","V95_beta","V95_pval","V96_beta","V96_pval",
                             "V97_beta","V97_pval","V98_beta","V98_pval","V99_beta","V99_pval","V100_beta","V100_pval")

write.csv(wholebrainst_ext, "100_random_samples_output_ext_final.csv") #write csv file

#conditional three-level growth model with global brain structure measures predicting the INT factor scores from correlated factors model
#looped across 100 dataframes of different random siblings included
#covariates: sex, age, scanner model dummies
#standardized betas  and p-values for intercepts (i.e., intb, intp) and slopes (i.e., slb, slp) were saved in a csv file
wholebrainst_int <- data.frame(x=c("CT_int","CT_sl","SA_int","SA_sl","cortvol_int","cortvol_sl","subcortvol_int","subcortvol_sl"))
for (i in 1:length(df)){ 
  CT_intb <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                            (1 + wave|site_id/id), #random intercept and slope for subject and site 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  CT_slb <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  CT_intp <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  CT_slp <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  SA_intb <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  SA_slb <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  SA_intp <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) +
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  SA_slp <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  cortvol_intb <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  cortvol_slb <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  cortvol_intp <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  cortvol_slp <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  subcortvol_intb <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  subcortvol_slb <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) +
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  subcortvol_intp <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  subcortvol_slp <- summary(lmer(scale(int_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  for (j in 2:101) #dataframe names
    newb <- c(CT_intb,CT_slb,SA_intb,SA_slb,cortvol_intb,cortvol_slb,subcortvol_intb,subcortvol_slb)
  wholebrainst_int <- cbind(wholebrainst_int,newb)
  names(wholebrainst_int)[ncol(wholebrainst_int)] <- names(sibs)[j]
  newp <- c(CT_intp,CT_slp,SA_intp,SA_slp,cortvol_intp,cortvol_slp,subcortvol_intp,subcortvol_slp)
  wholebrainst_int <- cbind(wholebrainst_int,newp)
  names(wholebrainst_int)[ncol(wholebrainst_int)] <- names(sibs)[j]
  
}
wholebrainst_int

names(wholebrainst_int) <- c("x","V1_beta","V1_pval","V2_beta","V2_pval","V3_beta","V3_pval","V4_beta","V4_pval",
                             "V5_beta","V5_pval","V6_beta","V6_pval","V7_beta","V7_pval","V8_beta","V8_pval",
                             "V9_beta","V9_pval","V10_beta","V10_pval","V11_beta","V11_pval","V12_beta","V12_pval",
                             "V13_beta","V13_pval","V14_beta","V14_pval","V15_beta","V15_pval","V16_beta","V16_pval",
                             "V17_beta","V17_pval","V18_beta","V18_pval","V19_beta","V19_pval","V20_beta","V20_pval",
                             "V21_beta","V21_pval","V22_beta","V22_pval","V23_beta","V23_pval","V24_beta","V24_pval",
                             "V25_beta","V25_pval","V26_beta","V26_pval","V27_beta","V27_pval","V28_beta","V28_pval",
                             "V29_beta","V29_pval","V30_beta","V30_pval","V31_beta","V31_pval","V32_beta","V32_pval",
                             "V33_beta","V33_pval","V34_beta","V34_pval","V35_beta","V35_pval","V36_beta","V36_pval",
                             "V37_beta","V37_pval","V38_beta","V38_pval","V39_beta","V39_pval","V40_beta","V40_pval",
                             "V41_beta","V41_pval","V42_beta","V42_pval","V43_beta","V43_pval","V44_beta","V44_pval",
                             "V45_beta","V45_pval","V46_beta","V46_pval","V47_beta","V47_pval","V48_beta","V48_pval",
                             "V49_beta","V49_pval","V50_beta","V50_pval","V51_beta","V51_pval","V52_beta","V52_pval",
                             "V53_beta","V53_pval","V54_beta","V54_pval","V55_beta","V55_pval","V56_beta","V56_pval",
                             "V57_beta","V57_pval","V58_beta","V58_pval","V59_beta","V59_pval","V60_beta","V60_pval",
                             "V61_beta","V61_pval","V62_beta","V62_pval","V63_beta","V63_pval","V64_beta","V64_pval",
                             "V65_beta","V65_pval","V66_beta","V66_pval","V67_beta","V67_pval","V68_beta","V68_pval",
                             "V69_beta","V69_pval","V70_beta","V70_pval","V71_beta","V71_pval","V72_beta","V72_pval",
                             "V73_beta","V73_pval","V74_beta","V74_pval","V75_beta","V75_pval","V76_beta","V76_pval",
                             "V77_beta","V77_pval","V78_beta","V78_pval","V79_beta","V79_pval","V80_beta","V80_pval",
                             "V81_beta","V81_pval","V82_beta","V82_pval","V83_beta","V83_pval","V84_beta","V84_pval",
                             "V85_beta","V85_pval","V86_beta","V86_pval","V87_beta","V87_pval","V88_beta","V88_pval",
                             "V89_beta","V89_pval","V90_beta","V90_pval","V91_beta","V91_pval","V92_beta","V92_pval",
                             "V93_beta","V93_pval","V94_beta","V94_pval","V95_beta","V95_pval","V96_beta","V96_pval",
                             "V97_beta","V97_pval","V98_beta","V98_pval","V99_beta","V99_pval","V100_beta","V100_pval")

write.csv(wholebrainst_int, "100_random_samples_output_int_final.csv") #write csv file

#conditional three-level growth model with global brain structure measures predicting the ND factor scores from correlated factors model
#looped across 100 dataframes of different random siblings included
#covariates: sex, age, scanner model dummies
#standardized betas  and p-values for intercepts (i.e., intb, intp) and slopes (i.e., slb, slp) were saved in a csv file
wholebrainst_nd <- data.frame(x=c("CT_int","CT_sl","SA_int","SA_sl","cortvol_int","cortvol_sl","subcortvol_int","subcortvol_sl"))
for (i in 1:length(df)){ 
  CT_intb <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                            (1 + wave|site_id/id), #random intercept and slope for subject and site 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  CT_slb <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  CT_intp <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  CT_slp <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  SA_intb <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  SA_slb <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  SA_intp <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) +
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  SA_slp <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  cortvol_intb <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  cortvol_slb <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  cortvol_intp <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  cortvol_slp <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  subcortvol_intb <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  subcortvol_slb <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) +
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  subcortvol_intp <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  subcortvol_slp <- summary(lmer(scale(nd_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  for (j in 2:101) #dataframe names
    newb <- c(CT_intb,CT_slb,SA_intb,SA_slb,cortvol_intb,cortvol_slb,subcortvol_intb,subcortvol_slb)
  wholebrainst_nd <- cbind(wholebrainst_nd,newb)
  names(wholebrainst_nd)[ncol(wholebrainst_nd)] <- names(sibs)[j]
  newp <- c(CT_intp,CT_slp,SA_intp,SA_slp,cortvol_intp,cortvol_slp,subcortvol_intp,subcortvol_slp)
  wholebrainst_nd <- cbind(wholebrainst_nd,newp)
  names(wholebrainst_nd)[ncol(wholebrainst_nd)] <- names(sibs)[j]
  
}
wholebrainst_nd

names(wholebrainst_nd) <- c("x","V1_beta","V1_pval","V2_beta","V2_pval","V3_beta","V3_pval","V4_beta","V4_pval",
                             "V5_beta","V5_pval","V6_beta","V6_pval","V7_beta","V7_pval","V8_beta","V8_pval",
                             "V9_beta","V9_pval","V10_beta","V10_pval","V11_beta","V11_pval","V12_beta","V12_pval",
                             "V13_beta","V13_pval","V14_beta","V14_pval","V15_beta","V15_pval","V16_beta","V16_pval",
                             "V17_beta","V17_pval","V18_beta","V18_pval","V19_beta","V19_pval","V20_beta","V20_pval",
                             "V21_beta","V21_pval","V22_beta","V22_pval","V23_beta","V23_pval","V24_beta","V24_pval",
                             "V25_beta","V25_pval","V26_beta","V26_pval","V27_beta","V27_pval","V28_beta","V28_pval",
                             "V29_beta","V29_pval","V30_beta","V30_pval","V31_beta","V31_pval","V32_beta","V32_pval",
                             "V33_beta","V33_pval","V34_beta","V34_pval","V35_beta","V35_pval","V36_beta","V36_pval",
                             "V37_beta","V37_pval","V38_beta","V38_pval","V39_beta","V39_pval","V40_beta","V40_pval",
                             "V41_beta","V41_pval","V42_beta","V42_pval","V43_beta","V43_pval","V44_beta","V44_pval",
                             "V45_beta","V45_pval","V46_beta","V46_pval","V47_beta","V47_pval","V48_beta","V48_pval",
                             "V49_beta","V49_pval","V50_beta","V50_pval","V51_beta","V51_pval","V52_beta","V52_pval",
                             "V53_beta","V53_pval","V54_beta","V54_pval","V55_beta","V55_pval","V56_beta","V56_pval",
                             "V57_beta","V57_pval","V58_beta","V58_pval","V59_beta","V59_pval","V60_beta","V60_pval",
                             "V61_beta","V61_pval","V62_beta","V62_pval","V63_beta","V63_pval","V64_beta","V64_pval",
                             "V65_beta","V65_pval","V66_beta","V66_pval","V67_beta","V67_pval","V68_beta","V68_pval",
                             "V69_beta","V69_pval","V70_beta","V70_pval","V71_beta","V71_pval","V72_beta","V72_pval",
                             "V73_beta","V73_pval","V74_beta","V74_pval","V75_beta","V75_pval","V76_beta","V76_pval",
                             "V77_beta","V77_pval","V78_beta","V78_pval","V79_beta","V79_pval","V80_beta","V80_pval",
                             "V81_beta","V81_pval","V82_beta","V82_pval","V83_beta","V83_pval","V84_beta","V84_pval",
                             "V85_beta","V85_pval","V86_beta","V86_pval","V87_beta","V87_pval","V88_beta","V88_pval",
                             "V89_beta","V89_pval","V90_beta","V90_pval","V91_beta","V91_pval","V92_beta","V92_pval",
                             "V93_beta","V93_pval","V94_beta","V94_pval","V95_beta","V95_pval","V96_beta","V96_pval",
                             "V97_beta","V97_pval","V98_beta","V98_pval","V99_beta","V99_pval","V100_beta","V100_pval")

write.csv(wholebrainst_nd, "100_random_samples_output_nd_final.csv") #write csv file

#conditional three-level growth model with global brain structure measures predicting the SOMAT factor scores from correlated factors model
#looped across 100 dataframes of different random siblings included
#covariates: sex, age, scanner model dummies
#standardized betas  and p-values for intercepts (i.e., intb, intp) and slopes (i.e., slb, slp) were saved in a csv file
wholebrainst_som <- data.frame(x=c("CT_int","CT_sl","SA_int","SA_sl","cortvol_int","cortvol_sl","subcortvol_int","subcortvol_sl"))
for (i in 1:length(df)){ 
  CT_intb <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                            (1 + wave|site_id/id), #random intercept and slope for subject and site 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  CT_slb <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  CT_intp <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  CT_slp <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  SA_intb <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  SA_slb <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  SA_intp <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) +
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  SA_slp <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  cortvol_intb <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  cortvol_slb <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  cortvol_intp <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  cortvol_slp <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  subcortvol_intb <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  subcortvol_slb <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) +
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  subcortvol_intp <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  subcortvol_slp <- summary(lmer(scale(som_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  for (j in 2:101) #dataframe names
    newb <- c(CT_intb,CT_slb,SA_intb,SA_slb,cortvol_intb,cortvol_slb,subcortvol_intb,subcortvol_slb)
  wholebrainst_som <- cbind(wholebrainst_som,newb)
  names(wholebrainst_som)[ncol(wholebrainst_som)] <- names(sibs)[j]
  newp <- c(CT_intp,CT_slp,SA_intp,SA_slp,cortvol_intp,cortvol_slp,subcortvol_intp,subcortvol_slp)
  wholebrainst_som <- cbind(wholebrainst_som,newp)
  names(wholebrainst_som)[ncol(wholebrainst_som)] <- names(sibs)[j]
  
}
wholebrainst_som

names(wholebrainst_som) <- c("x","V1_beta","V1_pval","V2_beta","V2_pval","V3_beta","V3_pval","V4_beta","V4_pval",
                             "V5_beta","V5_pval","V6_beta","V6_pval","V7_beta","V7_pval","V8_beta","V8_pval",
                             "V9_beta","V9_pval","V10_beta","V10_pval","V11_beta","V11_pval","V12_beta","V12_pval",
                             "V13_beta","V13_pval","V14_beta","V14_pval","V15_beta","V15_pval","V16_beta","V16_pval",
                             "V17_beta","V17_pval","V18_beta","V18_pval","V19_beta","V19_pval","V20_beta","V20_pval",
                             "V21_beta","V21_pval","V22_beta","V22_pval","V23_beta","V23_pval","V24_beta","V24_pval",
                             "V25_beta","V25_pval","V26_beta","V26_pval","V27_beta","V27_pval","V28_beta","V28_pval",
                             "V29_beta","V29_pval","V30_beta","V30_pval","V31_beta","V31_pval","V32_beta","V32_pval",
                             "V33_beta","V33_pval","V34_beta","V34_pval","V35_beta","V35_pval","V36_beta","V36_pval",
                             "V37_beta","V37_pval","V38_beta","V38_pval","V39_beta","V39_pval","V40_beta","V40_pval",
                             "V41_beta","V41_pval","V42_beta","V42_pval","V43_beta","V43_pval","V44_beta","V44_pval",
                             "V45_beta","V45_pval","V46_beta","V46_pval","V47_beta","V47_pval","V48_beta","V48_pval",
                             "V49_beta","V49_pval","V50_beta","V50_pval","V51_beta","V51_pval","V52_beta","V52_pval",
                             "V53_beta","V53_pval","V54_beta","V54_pval","V55_beta","V55_pval","V56_beta","V56_pval",
                             "V57_beta","V57_pval","V58_beta","V58_pval","V59_beta","V59_pval","V60_beta","V60_pval",
                             "V61_beta","V61_pval","V62_beta","V62_pval","V63_beta","V63_pval","V64_beta","V64_pval",
                             "V65_beta","V65_pval","V66_beta","V66_pval","V67_beta","V67_pval","V68_beta","V68_pval",
                             "V69_beta","V69_pval","V70_beta","V70_pval","V71_beta","V71_pval","V72_beta","V72_pval",
                             "V73_beta","V73_pval","V74_beta","V74_pval","V75_beta","V75_pval","V76_beta","V76_pval",
                             "V77_beta","V77_pval","V78_beta","V78_pval","V79_beta","V79_pval","V80_beta","V80_pval",
                             "V81_beta","V81_pval","V82_beta","V82_pval","V83_beta","V83_pval","V84_beta","V84_pval",
                             "V85_beta","V85_pval","V86_beta","V86_pval","V87_beta","V87_pval","V88_beta","V88_pval",
                             "V89_beta","V89_pval","V90_beta","V90_pval","V91_beta","V91_pval","V92_beta","V92_pval",
                             "V93_beta","V93_pval","V94_beta","V94_pval","V95_beta","V95_pval","V96_beta","V96_pval",
                             "V97_beta","V97_pval","V98_beta","V98_pval","V99_beta","V99_pval","V100_beta","V100_pval")

write.csv(wholebrainst_som, "100_random_samples_output_som_final.csv") #write csv file

#conditional three-level growth model with global brain structure measures predicting the DETACH factor scores from correlated factors model
#looped across 100 dataframes of different random siblings included
#covariates: sex, age, scanner model dummies
#standardized betas  and p-values for intercepts (i.e., intb, intp) and slopes (i.e., slb, slp) were saved in a csv file
wholebrainst_det <- data.frame(x=c("CT_int","CT_sl","SA_int","SA_sl","cortvol_int","cortvol_sl","subcortvol_int","subcortvol_sl"))
for (i in 1:length(df)){ 
  CT_intb <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                            (1 + wave|site_id/id), #random intercept and slope for subject and site 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  CT_slb <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) +
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  CT_intp <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  CT_slp <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scale(subjmeanct) + site_meanct_cent*wave + scale(subjmeanct*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  SA_intb <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  SA_slb <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  SA_intp <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) +
                            (1 + wave|site_id/id), 
                          data = df[[i]], 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  SA_slp <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scale(subjtotalarea) + site_totalarea_cent*wave + scale(subjtotalarea*wave) + 
                           (1 + wave|site_id/id), 
                         data = df[[i]], 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  cortvol_intb <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  cortvol_slb <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  cortvol_intp <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                 (1 + wave|site_id/id), 
                               data = df[[i]], 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  cortvol_slp <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scale(subjtotalvol) + site_totalvol_cent*wave + scale(subjtotalvol*wave) + 
                                (1 + wave|site_id/id), 
                              data = df[[i]], 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  subcortvol_intb <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,1]
  subcortvol_slb <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) +
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,1]
  subcortvol_intp <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                    (1 + wave|site_id/id), 
                                  data = df[[i]], 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[10,5]
  subcortvol_slp <- summary(lmer(scale(det_ho) ~ wave + sex + age + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scale(subjsubcort_vol) + site_subcort_vol_cent*wave + scale(subjsubcort_vol*wave) + 
                                   (1 + wave|site_id/id), 
                                 data = df[[i]], 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[11,5]
  for (j in 2:101) #dataframe names
    newb <- c(CT_intb,CT_slb,SA_intb,SA_slb,cortvol_intb,cortvol_slb,subcortvol_intb,subcortvol_slb)
  wholebrainst_det <- cbind(wholebrainst_det,newb)
  names(wholebrainst_det)[ncol(wholebrainst_det)] <- names(sibs)[j]
  newp <- c(CT_intp,CT_slp,SA_intp,SA_slp,cortvol_intp,cortvol_slp,subcortvol_intp,subcortvol_slp)
  wholebrainst_det <- cbind(wholebrainst_det,newp)
  names(wholebrainst_det)[ncol(wholebrainst_det)] <- names(sibs)[j]
  
}
wholebrainst_det

names(wholebrainst_det) <- c("x","V1_beta","V1_pval","V2_beta","V2_pval","V3_beta","V3_pval","V4_beta","V4_pval",
                            "V5_beta","V5_pval","V6_beta","V6_pval","V7_beta","V7_pval","V8_beta","V8_pval",
                            "V9_beta","V9_pval","V10_beta","V10_pval","V11_beta","V11_pval","V12_beta","V12_pval",
                            "V13_beta","V13_pval","V14_beta","V14_pval","V15_beta","V15_pval","V16_beta","V16_pval",
                            "V17_beta","V17_pval","V18_beta","V18_pval","V19_beta","V19_pval","V20_beta","V20_pval",
                            "V21_beta","V21_pval","V22_beta","V22_pval","V23_beta","V23_pval","V24_beta","V24_pval",
                            "V25_beta","V25_pval","V26_beta","V26_pval","V27_beta","V27_pval","V28_beta","V28_pval",
                            "V29_beta","V29_pval","V30_beta","V30_pval","V31_beta","V31_pval","V32_beta","V32_pval",
                            "V33_beta","V33_pval","V34_beta","V34_pval","V35_beta","V35_pval","V36_beta","V36_pval",
                            "V37_beta","V37_pval","V38_beta","V38_pval","V39_beta","V39_pval","V40_beta","V40_pval",
                            "V41_beta","V41_pval","V42_beta","V42_pval","V43_beta","V43_pval","V44_beta","V44_pval",
                            "V45_beta","V45_pval","V46_beta","V46_pval","V47_beta","V47_pval","V48_beta","V48_pval",
                            "V49_beta","V49_pval","V50_beta","V50_pval","V51_beta","V51_pval","V52_beta","V52_pval",
                            "V53_beta","V53_pval","V54_beta","V54_pval","V55_beta","V55_pval","V56_beta","V56_pval",
                            "V57_beta","V57_pval","V58_beta","V58_pval","V59_beta","V59_pval","V60_beta","V60_pval",
                            "V61_beta","V61_pval","V62_beta","V62_pval","V63_beta","V63_pval","V64_beta","V64_pval",
                            "V65_beta","V65_pval","V66_beta","V66_pval","V67_beta","V67_pval","V68_beta","V68_pval",
                            "V69_beta","V69_pval","V70_beta","V70_pval","V71_beta","V71_pval","V72_beta","V72_pval",
                            "V73_beta","V73_pval","V74_beta","V74_pval","V75_beta","V75_pval","V76_beta","V76_pval",
                            "V77_beta","V77_pval","V78_beta","V78_pval","V79_beta","V79_pval","V80_beta","V80_pval",
                            "V81_beta","V81_pval","V82_beta","V82_pval","V83_beta","V83_pval","V84_beta","V84_pval",
                            "V85_beta","V85_pval","V86_beta","V86_pval","V87_beta","V87_pval","V88_beta","V88_pval",
                            "V89_beta","V89_pval","V90_beta","V90_pval","V91_beta","V91_pval","V92_beta","V92_pval",
                            "V93_beta","V93_pval","V94_beta","V94_pval","V95_beta","V95_pval","V96_beta","V96_pval",
                            "V97_beta","V97_pval","V98_beta","V98_pval","V99_beta","V99_pval","V100_beta","V100_pval")

write.csv(wholebrainst_det, "100_random_samples_output_det_final.csv") #write csv file

#Create histograms of distributions of standardized estimates (Figure S3)
random_samples <- read.csv("100_random_samples_output_ready_use.csv",header=TRUE)

p1 <- ggplot(random_samples, aes(x=p_ct_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(p_ct_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Intercept", y="Count")
p1

p2 <- ggplot(random_samples, aes(x=p_ct_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(p_ct_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Slope", y="Count")
p2

p3 <- ggplot(random_samples, aes(x=p_sa_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(p_sa_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Intercept", y="Count")
p3

p4 <- ggplot(random_samples, aes(x=p_sa_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(p_sa_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Slope", y="Count")
p4

p5 <- ggplot(random_samples, aes(x=p_cortvol_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(p_cortvol_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Intercept", y="Count")
p5

p6 <- ggplot(random_samples, aes(x=p_cortvol_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(p_cortvol_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Slope", y="Count")
p6

p7 <- ggplot(random_samples, aes(x=p_subcort_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(p_subcort_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Intercept", y="Count")
p7

p8 <- ggplot(random_samples, aes(x=p_subcort_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(p_subcort_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Slope", y="Count")
p8

patchwork1 <- (p5 | p3 | p7 | p1) /
  (p6 | p4 | p8 | p2)
patchwork1 + plot_annotation(
  title = 'Relations between Brain Structure and p Factor Scores',
)

ggsave("figureS3A.tiff", height=8, width=14, units="in", dpi=300)

p9 <- ggplot(random_samples, aes(x=ext_ct_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(ext_ct_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Intercept", y="Count")
p9

p10 <- ggplot(random_samples, aes(x=ext_ct_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(ext_ct_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Slope", y="Count")
p10

p11 <- ggplot(random_samples, aes(x=ext_sa_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(ext_sa_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Intercept", y="Count")
p11

p12 <- ggplot(random_samples, aes(x=ext_sa_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(ext_sa_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Slope", y="Count")
p12

p13 <- ggplot(random_samples, aes(x=ext_cortvol_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(ext_cortvol_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Intercept", y="Count")
p13

p14 <- ggplot(random_samples, aes(x=ext_cortvol_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(ext_cortvol_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Slope", y="Count")
p14

p15 <- ggplot(random_samples, aes(x=ext_subcort_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(ext_subcort_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Intercept", y="Count")
p15

p16 <- ggplot(random_samples, aes(x=ext_subcort_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(ext_subcort_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Slope", y="Count")
p16

patchwork2 <- (p13 | p11 | p15 | p9) /
  (p14 | p12 | p16 | p10)
patchwork2 + plot_annotation(
  title = 'Relations between Brain Structure and EXT Factor Scores',
)

ggsave("figureS3B.tiff", height=8, width=14, units="in", dpi=300)

p17 <- ggplot(random_samples, aes(x=int_ct_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(int_ct_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Intercept", y="Count")
p17

p18 <- ggplot(random_samples, aes(x=int_ct_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(int_ct_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Slope", y="Count")
p18

p19 <- ggplot(random_samples, aes(x=int_sa_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(int_sa_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Intercept", y="Count")
p19

p20 <- ggplot(random_samples, aes(x=int_sa_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(int_sa_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Slope", y="Count")
p20

p21 <- ggplot(random_samples, aes(x=int_cortvol_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(int_cortvol_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Intercept", y="Count")
p21

p22 <- ggplot(random_samples, aes(x=int_cortvol_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(int_cortvol_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Slope", y="Count")
p22

p23 <- ggplot(random_samples, aes(x=int_subcort_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(int_subcort_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Intercept", y="Count")
p23

p24 <- ggplot(random_samples, aes(x=int_subcort_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(int_subcort_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Slope", y="Count")
p24

patchwork3 <- (p21 | p19 | p23 | p17) /
  (p22 | p20 | p24 | p18)
patchwork3 + plot_annotation(
  title = 'Relations between Brain Structure and INT Factor Scores',
)

ggsave("figureS3C.tiff", height=8, width=14, units="in", dpi=300)

p25 <- ggplot(random_samples, aes(x=nd_ct_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(nd_ct_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Intercept", y="Count")
p25

p26 <- ggplot(random_samples, aes(x=nd_ct_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(nd_ct_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Slope", y="Count")
p26

p27 <- ggplot(random_samples, aes(x=nd_sa_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(nd_sa_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Intercept", y="Count")
p27

p28 <- ggplot(random_samples, aes(x=nd_sa_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(nd_sa_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Slope", y="Count")
p28

p29 <- ggplot(random_samples, aes(x=nd_cortvol_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(nd_cortvol_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Intercept", y="Count")
p29

p30 <- ggplot(random_samples, aes(x=nd_cortvol_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(nd_cortvol_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Slope", y="Count")
p30

p31 <- ggplot(random_samples, aes(x=nd_subcort_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(nd_subcort_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Intercept", y="Count")
p31

p32 <- ggplot(random_samples, aes(x=nd_subcort_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(nd_subcort_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Slope", y="Count")
p32

patchwork4 <- (p29 | p27 | p31 | p25) /
  (p30 | p28 | p32 | p26)
patchwork4 + plot_annotation(
  title = 'Relations between Brain Structure and ND Factor Scores',
)

ggsave("figureS3D.tiff", height=8, width=14, units="in", dpi=300)

p33 <- ggplot(random_samples, aes(x=som_ct_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(som_ct_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Intercept", y="Count")
p33

p34 <- ggplot(random_samples, aes(x=som_ct_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(som_ct_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Slope", y="Count")
p34

p35 <- ggplot(random_samples, aes(x=som_sa_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(som_sa_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Intercept", y="Count")
p35

p36 <- ggplot(random_samples, aes(x=som_sa_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(som_sa_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Slope", y="Count")
p36

p37 <- ggplot(random_samples, aes(x=som_cortvol_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(som_cortvol_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Intercept", y="Count")
p37

p38 <- ggplot(random_samples, aes(x=som_cortvol_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(som_cortvol_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Slope", y="Count")
p38

p39 <- ggplot(random_samples, aes(x=som_subcort_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(som_subcort_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Intercept", y="Count")
p39

p40 <- ggplot(random_samples, aes(x=som_subcort_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(som_subcort_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Slope", y="Count")
p40

patchwork5 <- (p37 | p35 | p39 | p33) /
  (p38 | p36 | p40 | p34)
patchwork5 + plot_annotation(
  title = 'Relations between Brain Structure and SOMAT Factor Scores',
)

ggsave("figureS3E.tiff", height=8, width=14, units="in", dpi=300)

p41 <- ggplot(random_samples, aes(x=det_ct_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(det_ct_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Intercept", y="Count")
p41

p42 <- ggplot(random_samples, aes(x=det_ct_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(det_ct_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Mean CT and Slope", y="Count")
p42

p43 <- ggplot(random_samples, aes(x=det_sa_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(det_sa_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Intercept", y="Count")
p43

p44 <- ggplot(random_samples, aes(x=det_sa_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(det_sa_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total SA and Slope", y="Count")
p44

p45 <- ggplot(random_samples, aes(x=det_cortvol_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(det_cortvol_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Intercept", y="Count")
p45

p46 <- ggplot(random_samples, aes(x=det_cortvol_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(det_cortvol_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Cortical Volume and Slope", y="Count")
p46

p47 <- ggplot(random_samples, aes(x=det_subcort_int_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(det_subcort_int_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Intercept", y="Count")
p47

p48 <- ggplot(random_samples, aes(x=det_subcort_sl_beta)) + geom_histogram(binwidth = 0.001, fill="white", color="black") + 
  geom_vline(aes(xintercept=mean(det_subcort_sl_beta)), color="blue", linetype="dashed", size=1) +
  labs(x="Total Subcortical Volume and Slope", y="Count")
p48

patchwork6 <- (p45 | p43 | p47 | p41) /
  (p46 | p44 | p48 | p42)
patchwork6 + plot_annotation(
  title = 'Relations between Brain Structure and DETACH Factor Scores',
)

ggsave("figureS3F.tiff", height=8, width=14, units="in", dpi=300)
dev.off()
