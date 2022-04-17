#Title: Longitudinal Three-Level Growth Models of Brain Structure Relations with Psychopathology Trajectories with ABCD Study Release 4.0
#Author: Adrienne L. Romer, PhD
#Date: 04/15/2022

#clear workspace
rm(list=ls())

##Load libraries, read files, data management

#load packages
library(psych)              # use for descriptive statistics
library(ggplot2)            # use for plots
library(lme4)               # fits mixed models
library(lmerTest)           # provides t-tests for fixed effects
library(nlme)               # fits mixed models with hetero residuals
library(dplyr)              # load dplyr package
library(performance)        #computes ICC
library(interactions)       # use for calculating simple intercepts and slopes
library(plyr)
library(ggseg)
library(ggseg3d)
library(forestplot)

#set working directory
setwd("V:/ABCD/Data/Release4.0/Analysis/MLM")

#import data
pfactor <- read.csv("abcd_p_factor_brain_longitudinal_unst_fs_long_rel4.csv",header=TRUE)

#recode all missing values from -99 to NA
pfactor[pfactor == -99] <- NA

#recode wave from values of 1,2,3 to 0,1,2
pfactor$wave <- pfactor$Index1 - 1

#remove missing values
pfactor_red <- pfactor[which(pfactor$wave!="NA" & pfactor$p!="NA" & pfactor$Hispanic_only!="NA"),]

#create new variable renaming nums to id
pfactor_red$id <- pfactor_red$nums

#create new variables renaming TICs of no interest
pfactor_red$age <- pfactor_red$interview_age_baseline
pfactor_red$sex <- pfactor_red$sex_coded
pfactor_red$black <- pfactor_red$Black_onerace
pfactor_red$asian <- pfactor_red$Asian_onerace
pfactor_red$hisp <- pfactor_red$Hispanic_only
pfactor_red$other <- pfactor_red$Other_multiracial

##Descriptive statistics, means, and histograms

#examine descriptives of factor scores to get sense of functional form of mean trajectory
describeBy(pfactor_red$p, group=pfactor_red$wave)
describeBy(pfactor_red$ext, group=pfactor_red$wave)
describeBy(pfactor_red$int, group=pfactor_red$wave)
describeBy(pfactor_red$nd, group=pfactor_red$wave)
describeBy(pfactor_red$som, group=pfactor_red$wave)
describeBy(pfactor_red$det, group=pfactor_red$wave)

#plot means of factor scores by wave
meansp <- aggregate(p~wave,pfactor_red,mean)
meansext <- aggregate(ext~wave,pfactor_red,mean)
meansint <- aggregate(int~wave,pfactor_red,mean)
meansnd <- aggregate(nd~wave,pfactor_red,mean)
meanssom <- aggregate(som~wave,pfactor_red,mean)
meansdet <- aggregate(det~wave,pfactor_red,mean)

ggplot(meansp,aes(wave, p)) +
  geom_line(inherit.aes=TRUE) +
  #geom_smooth(method='lm',se=FALSE,color="red") + 
  geom_point(color='black') +
  labs(x = "Wave  - 1") +
  labs(y = "Mean p Factor Scores")

ggplot(meansext,aes(wave, ext)) +
  geom_line(inherit.aes=TRUE) +
  #geom_smooth(method='lm',se=FALSE,color="red") + 
  geom_point(color='black') +
  labs(x = "Wave  - 1") +
  labs(y = "Mean EXT Factor Scores")

ggplot(meansint,aes(wave, int)) +
  geom_line(inherit.aes=TRUE) +
  #geom_smooth(method='lm',se=FALSE,color="red") + 
  geom_point(color='black') +
  labs(x = "Wave  - 1") +
  labs(y = "Mean INT Factor Scores")

ggplot(meansnd,aes(wave, nd)) +
  geom_line(inherit.aes=TRUE) +
  #geom_smooth(method='lm',se=FALSE,color="red") + 
  geom_point(color='black') +
  labs(x = "Wave  - 1") +
  labs(y = "Mean ND Factor Scores")

ggplot(meanssom,aes(wave, som)) +
  geom_line(inherit.aes=TRUE) +
  #geom_smooth(method='lm',se=FALSE,color="red") + 
  geom_point(color='black') +
  labs(x = "Wave  - 1") +
  labs(y = "Mean SOM Factor Scores")

ggplot(meansdet,aes(wave, det)) +
  geom_line(inherit.aes=TRUE) +
  #geom_smooth(method='lm',se=FALSE,color="red") + 
  geom_point(color='black') +
  labs(x = "Wave  - 1") +
  labs(y = "Mean DET Factor Scores")

#scale p, ext, int, nd, som, det factor scores
pfactor_red$scalep <- scale(pfactor_red$p)
pfactor_red$scaleext <- scale(pfactor_red$ext)
pfactor_red$scaleint <- scale(pfactor_red$int)
pfactor_red$scalend <- scale(pfactor_red$nd)
pfactor_red$scalesom <- scale(pfactor_red$som)
pfactor_red$scaledet <- scale(pfactor_red$det)

#Histograms of factor scores over wave

hist(pfactor_red$scalep,
     main="Marginal Histogram",
     xlab="p Factor Scores",
     xlim=c(-1, 8),
     ylim=c(0, 3.5),
     prob=TRUE,
     breaks=40)

hist(pfactor_red$scaleext,
     main="Marginal Histogram",
     xlab="EXT Factor Scores",
     xlim=c(-1, 8),
     ylim=c(0, 3.5),
     prob=TRUE,
     breaks=40)

hist(pfactor_red$scaleint,
     main="Marginal Histogram",
     xlab="INT Factor Scores",
     xlim=c(-1, 8),
     ylim=c(0, 3),
     prob=TRUE,
     breaks=40)

hist(pfactor_red$scalend,
     main="Marginal Histogram",
     xlab="ND Factor Scores",
     xlim=c(-1, 8),
     ylim=c(0, 3.5),
     prob=TRUE,
     breaks=40)

hist(pfactor_red$scalesom,
     main="Marginal Histogram",
     xlab="SOMAT Factor Scores",
     xlim=c(-1, 8),
     ylim=c(0, 3.5),
     prob=TRUE,
     breaks=40)

hist(pfactor_red$scaledet,
     main="Marginal Histogram",
     xlab="DETACH Factor Scores",
     xlim=c(-1, 8),
     ylim=c(0, 3.5),
     prob=TRUE,
     breaks=40)

##three level MLM of repeated measures of psychopathology factor scores nested within subject nested within site

#random effects ANOVA and intrasubject correlations of factor scores
re_ANOVAp_3level_p <- lmer(p ~ (1|site_id/id), data = pfactor_red)
summary(re_ANOVAp_3level_p)
icc(re_ANOVAp_3level_p)

re_ANOVAext_3level_ext <- lmer(ext ~ (1|site_id/id), data = pfactor_red)
summary(re_ANOVAext_3level_ext)
icc(re_ANOVAext_3level_ext)

re_ANOVAint_3level_int <- lmer(int ~ (1|site_id/id), data = pfactor_red)
summary(re_ANOVAint_3level_int)
icc(re_ANOVAint_3level_int)

re_ANOVAnd_3level_nd <- lmer(nd ~ (1|site_id/id), data = pfactor_red)
summary(re_ANOVAnd_3level_nd)
icc(re_ANOVAnd_3level_nd)

re_ANOVAsom_3level_som <- lmer(som ~ (1|site_id/id), data = pfactor_red)
summary(re_ANOVAsom_3level_som)
icc(re_ANOVAsom_3level_som)

re_ANOVAdet_3level_det <- lmer(det ~ (1|site_id/id), data = pfactor_red)
summary(re_ANOVAdet_3level_det)
icc(re_ANOVAdet_3level_det)


#Unconditional three-level linear growth model with random intercept and slope and homoscedastic residuals for each of the factor scores
uncond_growth_3level_model_p <- lmer(p ~ wave + (1 + wave|site_id/id), 
                                   data = pfactor_red, 
                                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(uncond_growth_3level_model_p)

uncond_growth_3level_model_ext <- lmer(ext ~ wave + (1 + wave|site_id/id), 
                                     data = pfactor_red, 
                                     control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(uncond_growth_3level_model_ext)

uncond_growth_3level_model_int <- lmer(int ~ wave + (1 + wave|site_id/id), 
                                       data = pfactor_red, 
                                       control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(uncond_growth_3level_model_int)

uncond_growth_3level_model_nd <- lmer(nd ~ wave + (1 + wave|site_id/id), 
                                       data = pfactor_red, 
                                       control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(uncond_growth_3level_model_nd)

uncond_growth_3level_model_som <- lmer(som ~ wave + (1 + wave|site_id/id), 
                                       data = pfactor_red, 
                                       control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(uncond_growth_3level_model_som)

uncond_growth_3level_model_det <- lmer(det ~ wave + (1 + wave|site_id/id), 
                                       data = pfactor_red, 
                                       control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(uncond_growth_3level_model_det)


##Global structure analysis

#rescale and rename global brain variables
pfactor_red$wb_cort_vol <- pfactor_red$smri_vol_cdk_total/10000
pfactor_red$meanwb_ct <- pfactor_red$smri_thick_cdk_mean
pfactor_red$wb_cort_area <- pfactor_red$smri_area_cdk_total/10000
pfactor_red$subcort_vol <- pfactor_red$smri_vol_scs_subcorticalgv/1000

#Conditional three-level linear growth model

#Group-mean center global brain structure predictors by site 
groupmeans <- aggregate(cbind(meanwb_ct,wb_cort_vol,wb_cort_area,subcort_vol)~site_id, pfactor_red, mean)
names(groupmeans) <- c("site_id","site_meanct","site_totalvol","site_totalarea","site_subcort_vol")
groupmeans$site_meanct_cent <- groupmeans$site_meanct - mean(groupmeans$site_meanct)
groupmeans$site_totalvol_cent <- groupmeans$site_totalvol - mean(groupmeans$site_totalvol)
groupmeans$site_totalarea_cent <- groupmeans$site_totalarea - mean(groupmeans$site_totalarea)
groupmeans$site_subcort_vol_cent <- groupmeans$site_subcort_vol - mean(groupmeans$site_subcort_vol)

#merge centered variables (groupmeans) with pfactor_red dataset
pfactor_red <- join(pfactor_red, groupmeans, by="site_id")

#Center global brain structure variables by subject 
pfactor_red$subjmeanct <- pfactor_red$meanwb_ct - pfactor_red$site_meanct
pfactor_red$subjtotalvol <- pfactor_red$wb_cort_vol - pfactor_red$site_totalvol
pfactor_red$subjtotalarea <- pfactor_red$wb_cort_area - pfactor_red$site_totalarea
pfactor_red$subjsubcort_vol <- pfactor_red$subcort_vol - pfactor_red$site_subcort_vol

#scale subjmeanct, subjmeanvol, subjmeanarea 
pfactor_red$scalesubjmeanct <- scale(pfactor_red$subjmeanct)
pfactor_red$scalesubjtotalvol <- scale(pfactor_red$subjtotalvol)
pfactor_red$scalesubjtotalarea <- scale(pfactor_red$subjtotalarea)
pfactor_red$scalesubjsubcort_vol <- scale(pfactor_red$subjsubcort_vol)
pfactor_red$scalesite_totalvol_cent <- scale(pfactor_red$site_totalvol_cent)
pfactor_red$scalesite_subcort_vol_cent <- scale(pfactor_red$site_subcort_vol_cent)
pfactor_red$scalesite_totalarea_cent <- scale(pfactor_red$site_totalarea_cent)
pfactor_red$scalesite_meanct_cent <- scale(pfactor_red$site_meanct_cent)

#conditional three-level growth model with global brain structure measures predicting the psychopathology factor scores (looped) 
#covariates: sex, age, race/ethnicity dummies, scanner model dummies
#standardized betas, standard errors, and p-values for intercepts (i.e., intb, intse, intp) and slopes (i.e., slb, slse, slp) were saved in a csv file

wholebrainst <- data.frame(x=c("CT_int","CT_sl","SA_int","SA_sl","cortvol_int","cortvol_sl","subcortvol_int","subcortvol_sl"))
for (i in 420:425){ # columns are factor scores
  CT_intb <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                             site_meanct_cent + scalesubjmeanct + site_meanct_cent*wave + scalesubjmeanct*wave +
                            (1 + wave|site_id/id), #random intercept and slope for subject and site 
                          data = pfactor_red, 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,1]
  CT_slb <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scalesubjmeanct + site_meanct_cent*wave + scalesubjmeanct*wave +
                           (1 + wave|site_id/id), 
                         data = pfactor_red, 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,1]
  CT_intse <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                             site_meanct_cent + scalesubjmeanct + site_meanct_cent*wave + scalesubjmeanct*wave + 
                             (1 + wave|site_id/id), 
                           data = pfactor_red, 
                           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,2]
  CT_slse <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scalesubjmeanct + site_meanct_cent*wave + scalesubjmeanct*wave + 
                            (1 + wave|site_id/id), 
                          data = pfactor_red, 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,2]
  CT_intp <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                            site_meanct_cent + scalesubjmeanct + site_meanct_cent*wave + scalesubjmeanct*wave + 
                            (1 + wave|site_id/id), 
                          data = pfactor_red, 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,5]
  CT_slp <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                           site_meanct_cent + scalesubjmeanct + site_meanct_cent*wave + scalesubjmeanct*wave + 
                           (1 + wave|site_id/id), 
                         data = pfactor_red, 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,5]
  SA_intb <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scalesubjtotalarea + site_totalarea_cent*wave + scalesubjtotalarea*wave + 
                            (1 + wave|site_id/id), 
                          data = pfactor_red, 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,1]
  SA_slb <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scalesubjtotalarea + site_totalarea_cent*wave + scalesubjtotalarea*wave + 
                           (1 + wave|site_id/id), 
                         data = pfactor_red, 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,1]
  SA_intse <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                             site_totalarea_cent + scalesubjtotalarea + site_totalarea_cent*wave + scalesubjtotalarea*wave +
                             (1 + wave|site_id/id), 
                           data = pfactor_red, 
                           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,2]
  SA_slse <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scalesubjtotalarea + site_totalarea_cent*wave + scalesubjtotalarea*wave +
                            (1 + wave|site_id/id), 
                          data = pfactor_red, 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,2]
  SA_intp <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                            site_totalarea_cent + scalesubjtotalarea + site_totalarea_cent*wave + scalesubjtotalarea*wave +
                            (1 + wave|site_id/id), 
                          data = pfactor_red, 
                          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,5]
  SA_slp <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                           site_totalarea_cent + scalesubjtotalarea + site_totalarea_cent*wave + scalesubjtotalarea*wave + 
                           (1 + wave|site_id/id), 
                         data = pfactor_red, 
                         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,5]
  cortvol_intb <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scalesubjtotalvol + site_totalvol_cent*wave + scalesubjtotalvol*wave + 
                                 (1 + wave|site_id/id), 
                               data = pfactor_red, 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,1]
  cortvol_slb <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scalesubjtotalvol + site_totalvol_cent*wave + scalesubjtotalvol*wave + 
                                (1 + wave|site_id/id), 
                              data = pfactor_red, 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,1]
  cortvol_intse <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                  site_totalvol_cent + scalesubjtotalvol + site_totalvol_cent*wave + scalesubjtotalvol*wave + 
                                  (1 + wave|site_id/id), 
                                data = pfactor_red, 
                                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,2]
  cortvol_slse <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scalesubjtotalvol + site_totalvol_cent*wave + scalesubjtotalvol*wave + 
                                 (1 + wave|site_id/id), 
                               data = pfactor_red, 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,2]
  cortvol_intp <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                 site_totalvol_cent + scalesubjtotalvol + site_totalvol_cent*wave + scalesubjtotalvol*wave + 
                                 (1 + wave|site_id/id), 
                               data = pfactor_red, 
                               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,5]
  cortvol_slp <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                site_totalvol_cent + scalesubjtotalvol + site_totalvol_cent*wave + scalesubjtotalvol*wave + 
                                (1 + wave|site_id/id), 
                              data = pfactor_red, 
                              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,5]
  subcortvol_intb <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scalesubjsubcort_vol + site_subcort_vol_cent*wave + scalesubjsubcort_vol*wave + 
                                    (1 + wave|site_id/id), 
                                  data = pfactor_red, 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,1]
  subcortvol_slb <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scalesubjsubcort_vol + site_subcort_vol_cent*wave + scalesubjsubcort_vol*wave +
                                   (1 + wave|site_id/id), 
                                 data = pfactor_red, 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,1]
  subcortvol_intse <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                     site_subcort_vol_cent + scalesubjsubcort_vol + site_subcort_vol_cent*wave + scalesubjsubcort_vol*wave +
                                     (1 + wave|site_id/id), 
                                   data = pfactor_red, 
                                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,2]
  subcortvol_slse <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scalesubjsubcort_vol + site_subcort_vol_cent*wave + scalesubjsubcort_vol*wave +
                                    (1 + wave|site_id/id), 
                                  data = pfactor_red, 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,2]
  subcortvol_intp <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                    site_subcort_vol_cent + scalesubjsubcort_vol + site_subcort_vol_cent*wave + scalesubjsubcort_vol*wave + 
                                    (1 + wave|site_id/id), 
                                  data = pfactor_red, 
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[14,5]
  subcortvol_slp <- summary(lmer(pfactor_red[[i]] ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +
                                   site_subcort_vol_cent + scalesubjsubcort_vol + site_subcort_vol_cent*wave + scalesubjsubcort_vol*wave + 
                                   (1 + wave|site_id/id), 
                                 data = pfactor_red, 
                                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))$coefficients[16,5]
  newb <- c(CT_intb,CT_slb,SA_intb,SA_slb,cortvol_intb,cortvol_slb,subcortvol_intb,subcortvol_slb)
  wholebrainst <- cbind(wholebrainst,newb)
  names(wholebrainst)[ncol(wholebrainst)] <- names(pfactor_red)[i]
  newse <- c(CT_intse,CT_slse,SA_intse,SA_slse,cortvol_intse,cortvol_slse,subcortvol_intse,subcortvol_slse)
  wholebrainst <- cbind(wholebrainst,newse)
  names(wholebrainst)[ncol(wholebrainst)] <- names(pfactor_red)[i]
  newp <- c(CT_intp,CT_slp,SA_intp,SA_slp,cortvol_intp,cortvol_slp,subcortvol_intp,subcortvol_slp)
  wholebrainst <- cbind(wholebrainst,newp)
  names(wholebrainst)[ncol(wholebrainst)] <- names(pfactor_red)[i]
}
wholebrainst

names(wholebrainst) <- c("x","p_beta","p_se","p_pval","ext_beta","ext_se","ext_pval","int_beta","int_se","int_pval",
                         "nd_beta","nd_se","nd_pval","som_beta","som_se","som_pval","det_beta","det_se","det_pval")

#calculate 95% CIs and create lower and upper bound variables 
wholebrainst$ci_lower_p <- wholebrainst$p_beta - 1.96*wholebrainst$p_se 
wholebrainst$ci_upper_p <- wholebrainst$p_beta + 1.96*wholebrainst$p_se 
wholebrainst$ci_lower_ext <- wholebrainst$ext_beta - 1.96*wholebrainst$ext_se 
wholebrainst$ci_upper_ext <- wholebrainst$ext_beta + 1.96*wholebrainst$ext_se 
wholebrainst$ci_lower_int <- wholebrainst$int_beta - 1.96*wholebrainst$int_se 
wholebrainst$ci_upper_int <- wholebrainst$int_beta + 1.96*wholebrainst$int_se 
wholebrainst$ci_lower_nd <- wholebrainst$nd_beta - 1.96*wholebrainst$nd_se 
wholebrainst$ci_upper_nd <- wholebrainst$nd_beta + 1.96*wholebrainst$nd_se 
wholebrainst$ci_lower_som <- wholebrainst$som_beta - 1.96*wholebrainst$som_se 
wholebrainst$ci_upper_som <- wholebrainst$som_beta + 1.96*wholebrainst$som_se
wholebrainst$ci_lower_det <- wholebrainst$det_beta - 1.96*wholebrainst$det_se 
wholebrainst$ci_upper_det <- wholebrainst$det_beta + 1.96*wholebrainst$det_se

write.csv(wholebrainst, "Global Structure MLM Analysis Output Standardized for Table.csv") #write csv file

#not looped 3level growth model

#Conditional three-level linear growth model: cross-level interactions, total volume predicting rates of change in p (switch out ext, int, nd, som, det for p)
model_cortvol <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma + 
                         scalesite_totalvol_cent + scalesubjtotalvol + scalesite_totalvol_cent*wave + scalesubjtotalvol*wave + 
                         (1 + wave|site_id/id), 
                       data = pfactor_red, 
                       control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(model_cortvol)
confint(model_cortvol)

model_subcort <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                             scalesite_subcort_vol_cent + scalesubjsubcort_vol + scalesite_subcort_vol_cent*wave + scalesubjsubcort_vol*wave +
                             (1 + wave|site_id/id), 
                           data = pfactor_red, 
                           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(model_subcort)
confint(model_subcort)

model_area <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                          scalesite_totalarea_cent + scalesubjtotalarea + scalesite_totalarea_cent*wave + scalesubjtotalarea*wave +
                          (1 + wave|site_id/id), 
                        data = pfactor_red, 
                        control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(model_area)
confint(model_area)

model_ct <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma + 
                        scalesite_meanct_cent + scalesubjmeanct + scalesite_meanct_cent*wave + scalesubjmeanct*wave + 
                        (1 + wave|site_id/id), 
                      data = pfactor_red, 
                      control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(model_ct)
confint(model_ct)

#Probing and plotting interaction of wave with subject-level mean CT for INT factor scores (Figure 1)
sim_slopes(model_ct, pred=wave, modx=scalesubjmeanct, modx.values=c(-1,0,1), 
           centered="none", cond.int=TRUE, johnson_neyman=FALSE, data = pfactor_red)

interact_plot(model_ct, pred = wave, modx = scalesubjmeanct, 
              modx.values=c(-1,0,1),
              modx.labels = c("1 SD Below Mean","Mean","1 SD Above Mean"),
              legend.main = "Mean Cortical Thickness",
              x.label = "Wave", y.label = "Internalizing Factor Scores")

##Parcel-wise analyses

#rescale and rename brain parcels
pfactor_red$banksstslh_vol <- pfactor_red$smri_vol_cdk_banksstslh/1000
pfactor_red$banksstsrh_vol <- pfactor_red$smri_vol_cdk_banksstsrh/1000
pfactor_red$cdacatelh_vol <- pfactor_red$smri_vol_cdk_cdacatelh/1000
pfactor_red$cdacaterh_vol <- pfactor_red$smri_vol_cdk_cdacaterh/1000
pfactor_red$cdmdfrlh_vol <- pfactor_red$smri_vol_cdk_cdmdfrlh/1000
pfactor_red$cdmdfrrh_vol <- pfactor_red$smri_vol_cdk_cdmdfrrh/1000
pfactor_red$cuneuslh_vol <- pfactor_red$smri_vol_cdk_cuneuslh/1000
pfactor_red$cuneusrh_vol <- pfactor_red$smri_vol_cdk_cuneusrh/1000
pfactor_red$ehinallh_vol <- pfactor_red$smri_vol_cdk_ehinallh/1000
pfactor_red$ehinalrh_vol <- pfactor_red$smri_vol_cdk_ehinalrh/1000
pfactor_red$fusiformlh_vol <- pfactor_red$smri_vol_cdk_fusiformlh/1000
pfactor_red$fusiformrh_vol <- pfactor_red$smri_vol_cdk_fusiformrh/1000
pfactor_red$ifpllh_vol <- pfactor_red$smri_vol_cdk_ifpllh/1000
pfactor_red$ifplrh_vol <- pfactor_red$smri_vol_cdk_ifplrh/1000
pfactor_red$iftmlh_vol <- pfactor_red$smri_vol_cdk_iftmlh/1000
pfactor_red$iftmrh_vol <- pfactor_red$smri_vol_cdk_iftmrh/1000
pfactor_red$ihcatelh_vol <- pfactor_red$smri_vol_cdk_ihcatelh/1000
pfactor_red$ihcaterh_vol <- pfactor_red$smri_vol_cdk_ihcaterh/1000
pfactor_red$locclh_vol <- pfactor_red$smri_vol_cdk_locclh/1000
pfactor_red$loccrh_vol <- pfactor_red$smri_vol_cdk_loccrh/1000
pfactor_red$lobfrlh_vol <- pfactor_red$smri_vol_cdk_lobfrlh/1000
pfactor_red$lobfrrh_vol <- pfactor_red$smri_vol_cdk_lobfrrh/1000
pfactor_red$linguallh_vol <- pfactor_red$smri_vol_cdk_linguallh/1000
pfactor_red$lingualrh_vol <- pfactor_red$smri_vol_cdk_lingualrh/1000
pfactor_red$mobfrlh_vol <- pfactor_red$smri_vol_cdk_mobfrlh/1000
pfactor_red$mobfrrh_vol <- pfactor_red$smri_vol_cdk_mobfrrh/1000
pfactor_red$mdtmlh_vol <- pfactor_red$smri_vol_cdk_mdtmlh/1000
pfactor_red$mdtmrh_vol <- pfactor_red$smri_vol_cdk_mdtmrh/1000
pfactor_red$parahpallh_vol <- pfactor_red$smri_vol_cdk_parahpallh/1000
pfactor_red$parahpalrh_vol <- pfactor_red$smri_vol_cdk_parahpalrh/1000
pfactor_red$paracnlh_vol <- pfactor_red$smri_vol_cdk_paracnlh/1000
pfactor_red$paracnrh_vol <- pfactor_red$smri_vol_cdk_paracnrh/1000
pfactor_red$parsopclh_vol <- pfactor_red$smri_vol_cdk_parsopclh/1000
pfactor_red$parsopcrh_vol <- pfactor_red$smri_vol_cdk_parsopcrh/1000
pfactor_red$parsobislh_vol <- pfactor_red$smri_vol_cdk_parsobislh/1000
pfactor_red$parsobisrh_vol <- pfactor_red$smri_vol_cdk_parsobisrh/1000
pfactor_red$parstgrislh_vol <- pfactor_red$smri_vol_cdk_parstgrislh/1000
pfactor_red$parstgrisrh_vol <- pfactor_red$smri_vol_cdk_parstgrisrh/1000
pfactor_red$pericclh_vol <- pfactor_red$smri_vol_cdk_pericclh/1000
pfactor_red$periccrh_vol <- pfactor_red$smri_vol_cdk_periccrh/1000
pfactor_red$postcnlh_vol <- pfactor_red$smri_vol_cdk_postcnlh/1000
pfactor_red$postcnrh_vol <- pfactor_red$smri_vol_cdk_postcnrh/1000
pfactor_red$ptcatelh_vol <- pfactor_red$smri_vol_cdk_ptcatelh/1000
pfactor_red$ptcaterh_vol <- pfactor_red$smri_vol_cdk_ptcaterh/1000
pfactor_red$precnlh_vol <- pfactor_red$smri_vol_cdk_precnlh/1000
pfactor_red$precnrh_vol <- pfactor_red$smri_vol_cdk_precnrh/1000
pfactor_red$pclh_vol <- pfactor_red$smri_vol_cdk_pclh/1000
pfactor_red$pcrh_vol <- pfactor_red$smri_vol_cdk_pcrh/1000
pfactor_red$rracatelh_vol <- pfactor_red$smri_vol_cdk_rracatelh/1000
pfactor_red$rracaterh_vol <- pfactor_red$smri_vol_cdk_rracaterh/1000
pfactor_red$rrmdfrlh_vol <- pfactor_red$smri_vol_cdk_rrmdfrlh/1000
pfactor_red$rrmdfrrh_vol <- pfactor_red$smri_vol_cdk_rrmdfrrh/1000
pfactor_red$sufrlh_vol <- pfactor_red$smri_vol_cdk_sufrlh/1000
pfactor_red$sufrrh_vol <- pfactor_red$smri_vol_cdk_sufrrh/1000
pfactor_red$supllh_vol <- pfactor_red$smri_vol_cdk_supllh/1000
pfactor_red$suplrh_vol <- pfactor_red$smri_vol_cdk_suplrh/1000
pfactor_red$sutmlh_vol <- pfactor_red$smri_vol_cdk_sutmlh/1000
pfactor_red$sutmrh_vol <- pfactor_red$smri_vol_cdk_sutmrh/1000
pfactor_red$smlh_vol <- pfactor_red$smri_vol_cdk_smlh/1000
pfactor_red$smrh_vol <- pfactor_red$smri_vol_cdk_smrh/1000
pfactor_red$frpolelh_vol <- pfactor_red$smri_vol_cdk_frpolelh/1000
pfactor_red$frpolerh_vol <- pfactor_red$smri_vol_cdk_frpolerh/1000
pfactor_red$tmpolelh_vol <- pfactor_red$smri_vol_cdk_tmpolelh/1000
pfactor_red$tmpolerh_vol <- pfactor_red$smri_vol_cdk_tmpolerh/1000
pfactor_red$trvtmlh_vol <- pfactor_red$smri_vol_cdk_trvtmlh/1000
pfactor_red$trvtmrh_vol <- pfactor_red$smri_vol_cdk_trvtmrh/1000
pfactor_red$insulalh_vol <- pfactor_red$smri_vol_cdk_insulalh/1000
pfactor_red$insularh_vol <- pfactor_red$smri_vol_cdk_insularh/1000

pfactor_red$crbcortexlh_vol <- pfactor_red$smri_vol_scs_crbcortexlh/1000
pfactor_red$tplh_vol <- pfactor_red$smri_vol_scs_tplh/1000
pfactor_red$caudatelh_vol <- pfactor_red$smri_vol_scs_caudatelh/1000
pfactor_red$putamenlh_vol <- pfactor_red$smri_vol_scs_putamenlh/1000
pfactor_red$pallidumlh_vol <- pfactor_red$smri_vol_scs_pallidumlh/1000
pfactor_red$bstem_vol <- pfactor_red$smri_vol_scs_bstem/1000
pfactor_red$hpuslh_vol <- pfactor_red$smri_vol_scs_hpuslh/1000
pfactor_red$amygdalalh_vol <- pfactor_red$smri_vol_scs_amygdalalh/1000
pfactor_red$aal_vol <- pfactor_red$smri_vol_scs_aal/1000
pfactor_red$vedclh_vol <- pfactor_red$smri_vol_scs_vedclh/1000
pfactor_red$crbcortexrh_vol <- pfactor_red$smri_vol_scs_crbcortexrh/1000
pfactor_red$tprh_vol <- pfactor_red$smri_vol_scs_tprh/1000
pfactor_red$caudaterh_vol <- pfactor_red$smri_vol_scs_caudaterh/1000
pfactor_red$putamenrh_vol <- pfactor_red$smri_vol_scs_putamenrh/1000
pfactor_red$pallidumrh_vol <- pfactor_red$smri_vol_scs_pallidumrh/1000
pfactor_red$hpusrh_vol  <- pfactor_red$smri_vol_scs_hpusrh/1000
pfactor_red$amygdalarh_vol <- pfactor_red$smri_vol_scs_amygdalarh/1000
pfactor_red$aar_vol <- pfactor_red$smri_vol_scs_aar/1000
pfactor_red$vedcrh_vol <- pfactor_red$smri_vol_scs_vedcrh/1000

pfactor_red$banksstslh_area <- pfactor_red$smri_area_cdk_banksstslh/1000
pfactor_red$banksstsrh_area <- pfactor_red$smri_area_cdk_banksstsrh/1000
pfactor_red$cdacatelh_area <- pfactor_red$smri_area_cdk_cdacatelh/1000
pfactor_red$cdacaterh_area <- pfactor_red$smri_area_cdk_cdacaterh/1000
pfactor_red$cdmdfrlh_area <- pfactor_red$smri_area_cdk_cdmdfrlh/1000
pfactor_red$cdmdfrrh_area <- pfactor_red$smri_area_cdk_cdmdfrrh/1000
pfactor_red$cuneuslh_area <- pfactor_red$smri_area_cdk_cuneuslh/1000
pfactor_red$cuneusrh_area <- pfactor_red$smri_area_cdk_cuneusrh/1000
pfactor_red$ehinallh_area <- pfactor_red$smri_area_cdk_ehinallh/1000
pfactor_red$ehinalrh_area <- pfactor_red$smri_area_cdk_ehinalrh/1000
pfactor_red$fusiformlh_area <- pfactor_red$smri_area_cdk_fusiformlh/1000
pfactor_red$fusiformrh_area <- pfactor_red$smri_area_cdk_fusiformrh/1000
pfactor_red$ifpllh_area <- pfactor_red$smri_area_cdk_ifpllh/1000
pfactor_red$ifplrh_area <- pfactor_red$smri_area_cdk_ifplrh/1000
pfactor_red$iftmlh_area <- pfactor_red$smri_area_cdk_iftmlh/1000
pfactor_red$iftmrh_area <- pfactor_red$smri_area_cdk_iftmrh/1000
pfactor_red$ihcatelh_area <- pfactor_red$smri_area_cdk_ihcatelh/1000
pfactor_red$ihcaterh_area <- pfactor_red$smri_area_cdk_ihcaterh/1000
pfactor_red$locclh_area <- pfactor_red$smri_area_cdk_locclh/1000
pfactor_red$loccrh_area <- pfactor_red$smri_area_cdk_loccrh/1000
pfactor_red$lobfrlh_area <- pfactor_red$smri_area_cdk_lobfrlh/1000
pfactor_red$lobfrrh_area <- pfactor_red$smri_area_cdk_lobfrrh/1000
pfactor_red$linguallh_area <- pfactor_red$smri_area_cdk_linguallh/1000
pfactor_red$lingualrh_area <- pfactor_red$smri_area_cdk_lingualrh/1000
pfactor_red$mobfrlh_area <- pfactor_red$smri_area_cdk_mobfrlh/1000
pfactor_red$mobfrrh_area <- pfactor_red$smri_area_cdk_mobfrrh/1000
pfactor_red$mdtmlh_area <- pfactor_red$smri_area_cdk_mdtmlh/1000
pfactor_red$mdtmrh_area <- pfactor_red$smri_area_cdk_mdtmrh/1000
pfactor_red$parahpallh_area <- pfactor_red$smri_area_cdk_parahpallh/1000
pfactor_red$parahpalrh_area <- pfactor_red$smri_area_cdk_parahpalrh/1000
pfactor_red$paracnlh_area <- pfactor_red$smri_area_cdk_paracnlh/1000
pfactor_red$paracnrh_area <- pfactor_red$smri_area_cdk_paracnrh/1000
pfactor_red$parsopclh_area <- pfactor_red$smri_area_cdk_parsopclh/1000
pfactor_red$parsopcrh_area <- pfactor_red$smri_area_cdk_parsopcrh/1000
pfactor_red$parsobislh_area <- pfactor_red$smri_area_cdk_parsobislh/1000
pfactor_red$parsobisrh_area <- pfactor_red$smri_area_cdk_parsobisrh/1000
pfactor_red$parstgrislh_area <- pfactor_red$smri_area_cdk_parstgrislh/1000
pfactor_red$parstgrisrh_area <- pfactor_red$smri_area_cdk_parstgrisrh/1000
pfactor_red$pericclh_area <- pfactor_red$smri_area_cdk_pericclh/1000
pfactor_red$periccrh_area <- pfactor_red$smri_area_cdk_periccrh/1000
pfactor_red$postcnlh_area <- pfactor_red$smri_area_cdk_postcnlh/1000
pfactor_red$postcnrh_area <- pfactor_red$smri_area_cdk_postcnrh/1000
pfactor_red$ptcatelh_area <- pfactor_red$smri_area_cdk_ptcatelh/1000
pfactor_red$ptcaterh_area <- pfactor_red$smri_area_cdk_ptcaterh/1000
pfactor_red$precnlh_area <- pfactor_red$smri_area_cdk_precnlh/1000
pfactor_red$precnrh_area <- pfactor_red$smri_area_cdk_precnrh/1000
pfactor_red$pclh_area <- pfactor_red$smri_area_cdk_pclh/1000
pfactor_red$pcrh_area <- pfactor_red$smri_area_cdk_pcrh/1000
pfactor_red$rracatelh_area <- pfactor_red$smri_area_cdk_rracatelh/1000
pfactor_red$rracaterh_area <- pfactor_red$smri_area_cdk_rracaterh/1000
pfactor_red$rrmdfrlh_area <- pfactor_red$smri_area_cdk_rrmdfrlh/1000
pfactor_red$rrmdfrrh_area <- pfactor_red$smri_area_cdk_rrmdfrrh/1000
pfactor_red$sufrlh_area <- pfactor_red$smri_area_cdk_sufrlh/1000
pfactor_red$sufrrh_area <- pfactor_red$smri_area_cdk_sufrrh/1000
pfactor_red$supllh_area <- pfactor_red$smri_area_cdk_supllh/1000
pfactor_red$suplrh_area <- pfactor_red$smri_area_cdk_suplrh/1000
pfactor_red$sutmlh_area <- pfactor_red$smri_area_cdk_sutmlh/1000
pfactor_red$sutmrh_area <- pfactor_red$smri_area_cdk_sutmrh/1000
pfactor_red$smlh_area <- pfactor_red$smri_area_cdk_smlh/1000
pfactor_red$smrh_area <- pfactor_red$smri_area_cdk_smrh/1000
pfactor_red$frpolelh_area <- pfactor_red$smri_area_cdk_frpolelh/1000
pfactor_red$frpolerh_area <- pfactor_red$smri_area_cdk_frpolerh/1000
pfactor_red$tmpolelh_area <- pfactor_red$smri_area_cdk_tmpolelh/1000
pfactor_red$tmpolerh_area <- pfactor_red$smri_area_cdk_tmpolerh/1000
pfactor_red$trvtmlh_area <- pfactor_red$smri_area_cdk_trvtmlh/1000
pfactor_red$trvtmrh_area <- pfactor_red$smri_area_cdk_trvtmrh/1000
pfactor_red$insulalh_area <- pfactor_red$smri_area_cdk_insulalh/1000
pfactor_red$insularh_area <- pfactor_red$smri_area_cdk_insularh/1000

pfactor_red$banksstslh_thick <- pfactor_red$smri_thick_cdk_banksstslh
pfactor_red$banksstsrh_thick <- pfactor_red$smri_thick_cdk_banksstsrh
pfactor_red$cdacatelh_thick <- pfactor_red$smri_thick_cdk_cdacatelh
pfactor_red$cdacaterh_thick <- pfactor_red$smri_thick_cdk_cdacaterh
pfactor_red$cdmdfrlh_thick <- pfactor_red$smri_thick_cdk_cdmdfrlh
pfactor_red$cdmdfrrh_thick <- pfactor_red$smri_thick_cdk_cdmdfrrh
pfactor_red$cuneuslh_thick <- pfactor_red$smri_thick_cdk_cuneuslh
pfactor_red$cuneusrh_thick <- pfactor_red$smri_thick_cdk_cuneusrh
pfactor_red$ehinallh_thick <- pfactor_red$smri_thick_cdk_ehinallh
pfactor_red$ehinalrh_thick <- pfactor_red$smri_thick_cdk_ehinalrh
pfactor_red$fusiformlh_thick <- pfactor_red$smri_thick_cdk_fusiformlh
pfactor_red$fusiformrh_thick <- pfactor_red$smri_thick_cdk_fusiformrh
pfactor_red$ifpllh_thick <- pfactor_red$smri_thick_cdk_ifpllh
pfactor_red$ifplrh_thick <- pfactor_red$smri_thick_cdk_ifplrh
pfactor_red$iftmlh_thick <- pfactor_red$smri_thick_cdk_iftmlh
pfactor_red$iftmrh_thick <- pfactor_red$smri_thick_cdk_iftmrh
pfactor_red$ihcatelh_thick <- pfactor_red$smri_thick_cdk_ihcatelh
pfactor_red$ihcaterh_thick <- pfactor_red$smri_thick_cdk_ihcaterh
pfactor_red$locclh_thick <- pfactor_red$smri_thick_cdk_locclh
pfactor_red$loccrh_thick <- pfactor_red$smri_thick_cdk_loccrh
pfactor_red$lobfrlh_thick <- pfactor_red$smri_thick_cdk_lobfrlh
pfactor_red$lobfrrh_thick <- pfactor_red$smri_thick_cdk_lobfrrh
pfactor_red$linguallh_thick <- pfactor_red$smri_thick_cdk_linguallh
pfactor_red$lingualrh_thick <- pfactor_red$smri_thick_cdk_lingualrh
pfactor_red$mobfrlh_thick <- pfactor_red$smri_thick_cdk_mobfrlh
pfactor_red$mobfrrh_thick <- pfactor_red$smri_thick_cdk_mobfrrh
pfactor_red$mdtmlh_thick <- pfactor_red$smri_thick_cdk_mdtmlh
pfactor_red$mdtmrh_thick <- pfactor_red$smri_thick_cdk_mdtmrh
pfactor_red$parahpallh_thick <- pfactor_red$smri_thick_cdk_parahpallh
pfactor_red$parahpalrh_thick <- pfactor_red$smri_thick_cdk_parahpalrh
pfactor_red$paracnlh_thick <- pfactor_red$smri_thick_cdk_paracnlh
pfactor_red$paracnrh_thick <- pfactor_red$smri_thick_cdk_paracnrh
pfactor_red$parsopclh_thick <- pfactor_red$smri_thick_cdk_parsopclh
pfactor_red$parsopcrh_thick <- pfactor_red$smri_thick_cdk_parsopcrh
pfactor_red$parsobislh_thick <- pfactor_red$smri_thick_cdk_parsobislh
pfactor_red$parsobisrh_thick <- pfactor_red$smri_thick_cdk_parsobisrh
pfactor_red$parstgrislh_thick <- pfactor_red$smri_thick_cdk_parstgrislh
pfactor_red$parstgrisrh_thick <- pfactor_red$smri_thick_cdk_parstgrisrh
pfactor_red$pericclh_thick <- pfactor_red$smri_thick_cdk_pericclh
pfactor_red$periccrh_thick <- pfactor_red$smri_thick_cdk_periccrh
pfactor_red$postcnlh_thick <- pfactor_red$smri_thick_cdk_postcnlh
pfactor_red$postcnrh_thick <- pfactor_red$smri_thick_cdk_postcnrh
pfactor_red$ptcatelh_thick <- pfactor_red$smri_thick_cdk_ptcatelh
pfactor_red$ptcaterh_thick <- pfactor_red$smri_thick_cdk_ptcaterh
pfactor_red$precnlh_thick <- pfactor_red$smri_thick_cdk_precnlh
pfactor_red$precnrh_thick <- pfactor_red$smri_thick_cdk_precnrh
pfactor_red$pclh_thick <- pfactor_red$smri_thick_cdk_pclh
pfactor_red$pcrh_thick <- pfactor_red$smri_thick_cdk_pcrh
pfactor_red$rracatelh_thick <- pfactor_red$smri_thick_cdk_rracatelh
pfactor_red$rracaterh_thick <- pfactor_red$smri_thick_cdk_rracaterh
pfactor_red$rrmdfrlh_thick <- pfactor_red$smri_thick_cdk_rrmdfrlh
pfactor_red$rrmdfrrh_thick <- pfactor_red$smri_thick_cdk_rrmdfrrh
pfactor_red$sufrlh_thick <- pfactor_red$smri_thick_cdk_sufrlh
pfactor_red$sufrrh_thick <- pfactor_red$smri_thick_cdk_sufrrh
pfactor_red$supllh_thick <- pfactor_red$smri_thick_cdk_supllh
pfactor_red$suplrh_thick <- pfactor_red$smri_thick_cdk_suplrh
pfactor_red$sutmlh_thick <- pfactor_red$smri_thick_cdk_sutmlh
pfactor_red$sutmrh_thick <- pfactor_red$smri_thick_cdk_sutmrh
pfactor_red$smlh_thick <- pfactor_red$smri_thick_cdk_smlh
pfactor_red$smrh_thick <- pfactor_red$smri_thick_cdk_smrh
pfactor_red$frpolelh_thick <- pfactor_red$smri_thick_cdk_frpolelh
pfactor_red$frpolerh_thick <- pfactor_red$smri_thick_cdk_frpolerh
pfactor_red$tmpolelh_thick <- pfactor_red$smri_thick_cdk_tmpolelh
pfactor_red$tmpolerh_thick <- pfactor_red$smri_thick_cdk_tmpolerh
pfactor_red$trvtmlh_thick <- pfactor_red$smri_thick_cdk_trvtmlh
pfactor_red$trvtmrh_thick <- pfactor_red$smri_thick_cdk_trvtmrh
pfactor_red$insulalh_thick <- pfactor_red$smri_thick_cdk_insulalh
pfactor_red$insularh_thick <- pfactor_red$smri_thick_cdk_insularh

#Centering parcel-wise brain volume predictors
groupmeans_vol <- aggregate(cbind(banksstslh_vol,cdacatelh_vol,cdmdfrlh_vol,cuneuslh_vol,ehinallh_vol,
                                  fusiformlh_vol,ifpllh_vol,iftmlh_vol,ihcatelh_vol,locclh_vol,
                                  lobfrlh_vol,linguallh_vol,mobfrlh_vol,mdtmlh_vol,parahpallh_vol,
                                  paracnlh_vol,parsopclh_vol,parsobislh_vol,parstgrislh_vol,
                                  pericclh_vol,postcnlh_vol,ptcatelh_vol,precnlh_vol,pclh_vol,rracatelh_vol,
                                  rrmdfrlh_vol,sufrlh_vol,supllh_vol,sutmlh_vol,smlh_vol,frpolelh_vol,
                                  tmpolelh_vol,trvtmlh_vol,insulalh_vol,banksstsrh_vol,cdacaterh_vol,
                                  cdmdfrrh_vol,cuneusrh_vol,ehinalrh_vol,fusiformrh_vol,ifplrh_vol,iftmrh_vol,
                                  ihcaterh_vol,loccrh_vol,lobfrrh_vol,lingualrh_vol,mobfrrh_vol,mdtmrh_vol,
                                  parahpalrh_vol,paracnrh_vol,parsopcrh_vol,parsobisrh_vol,parstgrisrh_vol,
                                  periccrh_vol,postcnrh_vol,ptcaterh_vol,precnrh_vol,pcrh_vol,rracaterh_vol,
                                  rrmdfrrh_vol,sufrrh_vol,suplrh_vol,sutmrh_vol,smrh_vol,frpolerh_vol,tmpolerh_vol,
                                  trvtmrh_vol,insularh_vol)~site_id, pfactor_red, mean)
names(groupmeans_vol) <- c("site_id","site_banksstslh_vol","site_cdacatelh_vol","site_cdmdfrlh_vol","site_cuneuslh_vol",
                           "site_ehinallh_vol","site_fusiformlh_vol","site_ifpllh_vol","site_iftmlh_vol","site_ihcatelh_vol",
                           "site_locclh_vol","site_lobfrlh_vol","site_linguallh_vol","site_mobfrlh_vol","site_mdtmlh_vol",
                           "site_parahpallh_vol","site_paracnlh_vol","site_parsopclh_vol","site_parsobislh_vol","site_parstgrislh_vol",
                           "site_pericclh_vol","site_postcnlh_vol","site_ptcatelh_vol","site_precnlh_vol","site_pclh_vol","site_rracatelh_vol",
                           "site_rrmdfrlh_vol","site_sufrlh_vol","site_supllh_vol","site_sutmlh_vol","site_smlh_vol","site_frpolelh_vol",
                           "site_tmpolelh_vol","site_trvtmlh_vol","site_insulalh_vol","site_banksstsrh_vol","site_cdacaterh_vol","site_cdmdfrrh_vol",
                           "site_cuneusrh_vol","site_ehinalrh_vol","site_fusiformrh_vol","site_ifplrh_vol","site_iftmrh_vol","site_ihcaterh_vol",
                           "site_loccrh_vol","site_lobfrrh_vol","site_lingualrh_vol","site_mobfrrh_vol","site_mdtmrh_vol",
                           "site_parahpalrh_vol","site_paracnrh_vol","site_parsopcrh_vol","site_parsobisrh_vol","site_parstgrisrh_vol",
                           "site_periccrh_vol","site_postcnrh_vol","site_ptcaterh_vol","site_precnrh_vol","site_pcrh_vol","site_rracaterh_vol",
                           "site_rrmdfrrh_vol","site_sufrrh_vol","site_suplrh_vol","site_sutmrh_vol","site_smrh_vol","site_frpolerh_vol",
                           "site_tmpolerh_vol","site_trvtmrh_vol","site_insularh_vol")
groupmeans_vol$site_banksstslh_vol_cent <- groupmeans_vol$site_banksstslh_vol - mean(groupmeans_vol$site_banksstslh_vol)
groupmeans_vol$site_cdacatelh_vol_cent <- groupmeans_vol$site_cdacatelh_vol - mean(groupmeans_vol$site_cdacatelh_vol)
groupmeans_vol$site_cdmdfrlh_vol_cent <- groupmeans_vol$site_cdmdfrlh_vol - mean(groupmeans_vol$site_cdmdfrlh_vol)
groupmeans_vol$site_cuneuslh_vol_cent <- groupmeans_vol$site_cuneuslh_vol - mean(groupmeans_vol$site_cuneuslh_vol)
groupmeans_vol$site_ehinallh_vol_cent <- groupmeans_vol$site_ehinallh_vol - mean(groupmeans_vol$site_ehinallh_vol)
groupmeans_vol$site_fusiformlh_vol_cent <- groupmeans_vol$site_fusiformlh_vol - mean(groupmeans_vol$site_fusiformlh_vol)
groupmeans_vol$site_ifpllh_vol_cent <- groupmeans_vol$site_ifpllh_vol - mean(groupmeans_vol$site_ifpllh_vol)
groupmeans_vol$site_iftmlh_vol_cent <- groupmeans_vol$site_iftmlh_vol - mean(groupmeans_vol$site_iftmlh_vol)
groupmeans_vol$site_ihcatelh_vol_cent <- groupmeans_vol$site_ihcatelh_vol - mean(groupmeans_vol$site_ihcatelh_vol)
groupmeans_vol$site_locclh_vol_cent <- groupmeans_vol$site_locclh_vol - mean(groupmeans_vol$site_locclh_vol)
groupmeans_vol$site_lobfrlh_vol_cent <- groupmeans_vol$site_lobfrlh_vol - mean(groupmeans_vol$site_lobfrlh_vol)
groupmeans_vol$site_linguallh_vol_cent <- groupmeans_vol$site_linguallh_vol - mean(groupmeans_vol$site_linguallh_vol)
groupmeans_vol$site_mobfrlh_vol_cent <- groupmeans_vol$site_mobfrlh_vol - mean(groupmeans_vol$site_mobfrlh_vol)
groupmeans_vol$site_mdtmlh_vol_cent <- groupmeans_vol$site_mdtmlh_vol - mean(groupmeans_vol$site_mdtmlh_vol)
groupmeans_vol$site_parahpallh_vol_cent <- groupmeans_vol$site_parahpallh_vol - mean(groupmeans_vol$site_parahpallh_vol)
groupmeans_vol$site_paracnlh_vol_cent <- groupmeans_vol$site_paracnlh_vol - mean(groupmeans_vol$site_paracnlh_vol)
groupmeans_vol$site_parsopclh_vol_cent <- groupmeans_vol$site_parsopclh_vol - mean(groupmeans_vol$site_parsopclh_vol)
groupmeans_vol$site_parsobislh_vol_cent <- groupmeans_vol$site_parsobislh_vol - mean(groupmeans_vol$site_parsobislh_vol)
groupmeans_vol$site_parstgrislh_vol_cent <- groupmeans_vol$site_parstgrislh_vol - mean(groupmeans_vol$site_parstgrislh_vol)
groupmeans_vol$site_pericclh_vol_cent <- groupmeans_vol$site_pericclh_vol - mean(groupmeans_vol$site_pericclh_vol)
groupmeans_vol$site_postcnlh_vol_cent <- groupmeans_vol$site_postcnlh_vol - mean(groupmeans_vol$site_postcnlh_vol)
groupmeans_vol$site_ptcatelh_vol_cent <- groupmeans_vol$site_ptcatelh_vol - mean(groupmeans_vol$site_ptcatelh_vol)
groupmeans_vol$site_precnlh_vol_cent <- groupmeans_vol$site_precnlh_vol - mean(groupmeans_vol$site_precnlh_vol)
groupmeans_vol$site_pclh_vol_cent <- groupmeans_vol$site_pclh_vol - mean(groupmeans_vol$site_pclh_vol)
groupmeans_vol$site_rracatelh_vol_cent <- groupmeans_vol$site_rracatelh_vol - mean(groupmeans_vol$site_rracatelh_vol)
groupmeans_vol$site_rrmdfrlh_vol_cent <- groupmeans_vol$site_rrmdfrlh_vol - mean(groupmeans_vol$site_rrmdfrlh_vol)
groupmeans_vol$site_sufrlh_vol_cent <- groupmeans_vol$site_sufrlh_vol - mean(groupmeans_vol$site_sufrlh_vol)
groupmeans_vol$site_supllh_vol_cent <- groupmeans_vol$site_supllh_vol - mean(groupmeans_vol$site_supllh_vol)
groupmeans_vol$site_sutmlh_vol_cent <- groupmeans_vol$site_sutmlh_vol - mean(groupmeans_vol$site_sutmlh_vol)
groupmeans_vol$site_smlh_vol_cent <- groupmeans_vol$site_smlh_vol - mean(groupmeans_vol$site_smlh_vol)
groupmeans_vol$site_frpolelh_vol_cent <- groupmeans_vol$site_frpolelh_vol - mean(groupmeans_vol$site_frpolelh_vol)
groupmeans_vol$site_tmpolelh_vol_cent <- groupmeans_vol$site_tmpolelh_vol - mean(groupmeans_vol$site_tmpolelh_vol)
groupmeans_vol$site_trvtmlh_vol_cent <- groupmeans_vol$site_trvtmlh_vol - mean(groupmeans_vol$site_trvtmlh_vol)
groupmeans_vol$site_insulalh_vol_cent <- groupmeans_vol$site_insulalh_vol - mean(groupmeans_vol$site_insulalh_vol)
groupmeans_vol$site_banksstsrh_vol_cent <- groupmeans_vol$site_banksstsrh_vol - mean(groupmeans_vol$site_banksstsrh_vol)
groupmeans_vol$site_cdacaterh_vol_cent <- groupmeans_vol$site_cdacaterh_vol - mean(groupmeans_vol$site_cdacaterh_vol)
groupmeans_vol$site_cdmdfrrh_vol_cent <- groupmeans_vol$site_cdmdfrrh_vol - mean(groupmeans_vol$site_cdmdfrrh_vol)
groupmeans_vol$site_cuneusrh_vol_cent <- groupmeans_vol$site_cuneusrh_vol - mean(groupmeans_vol$site_cuneusrh_vol)
groupmeans_vol$site_ehinalrh_vol_cent <- groupmeans_vol$site_ehinalrh_vol - mean(groupmeans_vol$site_ehinalrh_vol)
groupmeans_vol$site_fusiformrh_vol_cent <- groupmeans_vol$site_fusiformrh_vol - mean(groupmeans_vol$site_fusiformrh_vol)
groupmeans_vol$site_ifplrh_vol_cent <- groupmeans_vol$site_ifplrh_vol - mean(groupmeans_vol$site_ifplrh_vol)
groupmeans_vol$site_iftmrh_vol_cent <- groupmeans_vol$site_iftmrh_vol - mean(groupmeans_vol$site_iftmrh_vol)
groupmeans_vol$site_ihcaterh_vol_cent <- groupmeans_vol$site_ihcaterh_vol - mean(groupmeans_vol$site_ihcaterh_vol)
groupmeans_vol$site_loccrh_vol_cent <- groupmeans_vol$site_loccrh_vol - mean(groupmeans_vol$site_loccrh_vol)
groupmeans_vol$site_lobfrrh_vol_cent <- groupmeans_vol$site_lobfrrh_vol - mean(groupmeans_vol$site_lobfrrh_vol)
groupmeans_vol$site_lingualrh_vol_cent <- groupmeans_vol$site_lingualrh_vol - mean(groupmeans_vol$site_lingualrh_vol)
groupmeans_vol$site_mobfrrh_vol_cent <- groupmeans_vol$site_mobfrrh_vol - mean(groupmeans_vol$site_mobfrrh_vol)
groupmeans_vol$site_mdtmrh_vol_cent <- groupmeans_vol$site_mdtmrh_vol - mean(groupmeans_vol$site_mdtmrh_vol)
groupmeans_vol$site_parahpalrh_vol_cent <- groupmeans_vol$site_parahpalrh_vol - mean(groupmeans_vol$site_parahpalrh_vol)
groupmeans_vol$site_paracnrh_vol_cent <- groupmeans_vol$site_paracnrh_vol - mean(groupmeans_vol$site_paracnrh_vol)
groupmeans_vol$site_parsopcrh_vol_cent <- groupmeans_vol$site_parsopcrh_vol - mean(groupmeans_vol$site_parsopcrh_vol)
groupmeans_vol$site_parsobisrh_vol_cent <- groupmeans_vol$site_parsobisrh_vol - mean(groupmeans_vol$site_parsobisrh_vol)
groupmeans_vol$site_parstgrisrh_vol_cent <- groupmeans_vol$site_parstgrisrh_vol - mean(groupmeans_vol$site_parstgrisrh_vol)
groupmeans_vol$site_periccrh_vol_cent <- groupmeans_vol$site_periccrh_vol - mean(groupmeans_vol$site_periccrh_vol)
groupmeans_vol$site_postcnrh_vol_cent <- groupmeans_vol$site_postcnrh_vol - mean(groupmeans_vol$site_postcnrh_vol)
groupmeans_vol$site_ptcaterh_vol_cent <- groupmeans_vol$site_ptcaterh_vol - mean(groupmeans_vol$site_ptcaterh_vol)
groupmeans_vol$site_precnrh_vol_cent <- groupmeans_vol$site_precnrh_vol - mean(groupmeans_vol$site_precnrh_vol)
groupmeans_vol$site_pcrh_vol_cent <- groupmeans_vol$site_pcrh_vol - mean(groupmeans_vol$site_pcrh_vol)
groupmeans_vol$site_rracaterh_vol_cent <- groupmeans_vol$site_rracaterh_vol - mean(groupmeans_vol$site_rracaterh_vol)
groupmeans_vol$site_rrmdfrrh_vol_cent <- groupmeans_vol$site_rrmdfrrh_vol - mean(groupmeans_vol$site_rrmdfrrh_vol)
groupmeans_vol$site_sufrrh_vol_cent <- groupmeans_vol$site_sufrrh_vol - mean(groupmeans_vol$site_sufrrh_vol)
groupmeans_vol$site_suplrh_vol_cent <- groupmeans_vol$site_suplrh_vol - mean(groupmeans_vol$site_suplrh_vol)
groupmeans_vol$site_sutmrh_vol_cent <- groupmeans_vol$site_sutmrh_vol - mean(groupmeans_vol$site_sutmrh_vol)
groupmeans_vol$site_smrh_vol_cent <- groupmeans_vol$site_smrh_vol - mean(groupmeans_vol$site_smrh_vol)
groupmeans_vol$site_frpolerh_vol_cent <- groupmeans_vol$site_frpolerh_vol - mean(groupmeans_vol$site_frpolerh_vol)
groupmeans_vol$site_tmpolerh_vol_cent <- groupmeans_vol$site_tmpolerh_vol - mean(groupmeans_vol$site_tmpolerh_vol)
groupmeans_vol$site_trvtmrh_vol_cent <- groupmeans_vol$site_trvtmrh_vol - mean(groupmeans_vol$site_trvtmrh_vol)
groupmeans_vol$site_insularh_vol_cent <- groupmeans_vol$site_insularh_vol - mean(groupmeans_vol$site_insularh_vol)

#Centering parcel-wise brain thickness predictors
groupmeans_thick <- aggregate(cbind(banksstslh_thick,cdacatelh_thick,cdmdfrlh_thick,cuneuslh_thick,ehinallh_thick,
                                    fusiformlh_thick,ifpllh_thick,iftmlh_thick,ihcatelh_thick,locclh_thick,
                                    lobfrlh_thick,linguallh_thick,mobfrlh_thick,mdtmlh_thick,parahpallh_thick,
                                    paracnlh_thick,parsopclh_thick,parsobislh_thick,parstgrislh_thick,
                                    pericclh_thick,postcnlh_thick,ptcatelh_thick,precnlh_thick,pclh_thick,rracatelh_thick,
                                    rrmdfrlh_thick,sufrlh_thick,supllh_thick,sutmlh_thick,smlh_thick,frpolelh_thick,
                                    tmpolelh_thick,trvtmlh_thick,insulalh_thick,banksstsrh_thick,cdacaterh_thick,
                                    cdmdfrrh_thick,cuneusrh_thick,ehinalrh_thick,fusiformrh_thick,ifplrh_thick,iftmrh_thick,
                                    ihcaterh_thick,loccrh_thick,lobfrrh_thick,lingualrh_thick,mobfrrh_thick,mdtmrh_thick,
                                    parahpalrh_thick,paracnrh_thick,parsopcrh_thick,parsobisrh_thick,parstgrisrh_thick,
                                    periccrh_thick,postcnrh_thick,ptcaterh_thick,precnrh_thick,pcrh_thick,rracaterh_thick,
                                    rrmdfrrh_thick,sufrrh_thick,suplrh_thick,sutmrh_thick,smrh_thick,frpolerh_thick,tmpolerh_thick,
                                    trvtmrh_thick,insularh_thick)~site_id, pfactor_red, mean)
names(groupmeans_thick) <- c("site_id","site_banksstslh_thick","site_cdacatelh_thick","site_cdmdfrlh_thick","site_cuneuslh_thick",
                             "site_ehinallh_thick","site_fusiformlh_thick","site_ifpllh_thick","site_iftmlh_thick","site_ihcatelh_thick",
                             "site_locclh_thick","site_lobfrlh_thick","site_linguallh_thick","site_mobfrlh_thick","site_mdtmlh_thick",
                             "site_parahpallh_thick","site_paracnlh_thick","site_parsopclh_thick","site_parsobislh_thick","site_parstgrislh_thick",
                             "site_pericclh_thick","site_postcnlh_thick","site_ptcatelh_thick","site_precnlh_thick","site_pclh_thick","site_rracatelh_thick",
                             "site_rrmdfrlh_thick","site_sufrlh_thick","site_supllh_thick","site_sutmlh_thick","site_smlh_thick","site_frpolelh_thick",
                             "site_tmpolelh_thick","site_trvtmlh_thick","site_insulalh_thick","site_banksstsrh_thick","site_cdacaterh_thick","site_cdmdfrrh_thick",
                             "site_cuneusrh_thick","site_ehinalrh_thick","site_fusiformrh_thick","site_ifplrh_thick","site_iftmrh_thick","site_ihcaterh_thick",
                             "site_loccrh_thick","site_lobfrrh_thick","site_lingualrh_thick","site_mobfrrh_thick","site_mdtmrh_thick",
                             "site_parahpalrh_thick","site_paracnrh_thick","site_parsopcrh_thick","site_parsobisrh_thick","site_parstgrisrh_thick",
                             "site_periccrh_thick","site_postcnrh_thick","site_ptcaterh_thick","site_precnrh_thick","site_pcrh_thick","site_rracaterh_thick",
                             "site_rrmdfrrh_thick","site_sufrrh_thick","site_suplrh_thick","site_sutmrh_thick","site_smrh_thick","site_frpolerh_thick",
                             "site_tmpolerh_thick","site_trvtmrh_thick","site_insularh_thick")
groupmeans_thick$site_banksstslh_thick_cent <- groupmeans_thick$site_banksstslh_thick - mean(groupmeans_thick$site_banksstslh_thick)
groupmeans_thick$site_cdacatelh_thick_cent <- groupmeans_thick$site_cdacatelh_thick - mean(groupmeans_thick$site_cdacatelh_thick)
groupmeans_thick$site_cdmdfrlh_thick_cent <- groupmeans_thick$site_cdmdfrlh_thick - mean(groupmeans_thick$site_cdmdfrlh_thick)
groupmeans_thick$site_cuneuslh_thick_cent <- groupmeans_thick$site_cuneuslh_thick - mean(groupmeans_thick$site_cuneuslh_thick)
groupmeans_thick$site_ehinallh_thick_cent <- groupmeans_thick$site_ehinallh_thick - mean(groupmeans_thick$site_ehinallh_thick)
groupmeans_thick$site_fusiformlh_thick_cent <- groupmeans_thick$site_fusiformlh_thick - mean(groupmeans_thick$site_fusiformlh_thick)
groupmeans_thick$site_ifpllh_thick_cent <- groupmeans_thick$site_ifpllh_thick - mean(groupmeans_thick$site_ifpllh_thick)
groupmeans_thick$site_iftmlh_thick_cent <- groupmeans_thick$site_iftmlh_thick - mean(groupmeans_thick$site_iftmlh_thick)
groupmeans_thick$site_ihcatelh_thick_cent <- groupmeans_thick$site_ihcatelh_thick - mean(groupmeans_thick$site_ihcatelh_thick)
groupmeans_thick$site_locclh_thick_cent <- groupmeans_thick$site_locclh_thick - mean(groupmeans_thick$site_locclh_thick)
groupmeans_thick$site_lobfrlh_thick_cent <- groupmeans_thick$site_lobfrlh_thick - mean(groupmeans_thick$site_lobfrlh_thick)
groupmeans_thick$site_linguallh_thick_cent <- groupmeans_thick$site_linguallh_thick - mean(groupmeans_thick$site_linguallh_thick)
groupmeans_thick$site_mobfrlh_thick_cent <- groupmeans_thick$site_mobfrlh_thick - mean(groupmeans_thick$site_mobfrlh_thick)
groupmeans_thick$site_mdtmlh_thick_cent <- groupmeans_thick$site_mdtmlh_thick - mean(groupmeans_thick$site_mdtmlh_thick)
groupmeans_thick$site_parahpallh_thick_cent <- groupmeans_thick$site_parahpallh_thick - mean(groupmeans_thick$site_parahpallh_thick)
groupmeans_thick$site_paracnlh_thick_cent <- groupmeans_thick$site_paracnlh_thick - mean(groupmeans_thick$site_paracnlh_thick)
groupmeans_thick$site_parsopclh_thick_cent <- groupmeans_thick$site_parsopclh_thick - mean(groupmeans_thick$site_parsopclh_thick)
groupmeans_thick$site_parsobislh_thick_cent <- groupmeans_thick$site_parsobislh_thick - mean(groupmeans_thick$site_parsobislh_thick)
groupmeans_thick$site_parstgrislh_thick_cent <- groupmeans_thick$site_parstgrislh_thick - mean(groupmeans_thick$site_parstgrislh_thick)
groupmeans_thick$site_pericclh_thick_cent <- groupmeans_thick$site_pericclh_thick - mean(groupmeans_thick$site_pericclh_thick)
groupmeans_thick$site_postcnlh_thick_cent <- groupmeans_thick$site_postcnlh_thick - mean(groupmeans_thick$site_postcnlh_thick)
groupmeans_thick$site_ptcatelh_thick_cent <- groupmeans_thick$site_ptcatelh_thick - mean(groupmeans_thick$site_ptcatelh_thick)
groupmeans_thick$site_precnlh_thick_cent <- groupmeans_thick$site_precnlh_thick - mean(groupmeans_thick$site_precnlh_thick)
groupmeans_thick$site_pclh_thick_cent <- groupmeans_thick$site_pclh_thick - mean(groupmeans_thick$site_pclh_thick)
groupmeans_thick$site_rracatelh_thick_cent <- groupmeans_thick$site_rracatelh_thick - mean(groupmeans_thick$site_rracatelh_thick)
groupmeans_thick$site_rrmdfrlh_thick_cent <- groupmeans_thick$site_rrmdfrlh_thick - mean(groupmeans_thick$site_rrmdfrlh_thick)
groupmeans_thick$site_sufrlh_thick_cent <- groupmeans_thick$site_sufrlh_thick - mean(groupmeans_thick$site_sufrlh_thick)
groupmeans_thick$site_supllh_thick_cent <- groupmeans_thick$site_supllh_thick - mean(groupmeans_thick$site_supllh_thick)
groupmeans_thick$site_sutmlh_thick_cent <- groupmeans_thick$site_sutmlh_thick - mean(groupmeans_thick$site_sutmlh_thick)
groupmeans_thick$site_smlh_thick_cent <- groupmeans_thick$site_smlh_thick - mean(groupmeans_thick$site_smlh_thick)
groupmeans_thick$site_frpolelh_thick_cent <- groupmeans_thick$site_frpolelh_thick - mean(groupmeans_thick$site_frpolelh_thick)
groupmeans_thick$site_tmpolelh_thick_cent <- groupmeans_thick$site_tmpolelh_thick - mean(groupmeans_thick$site_tmpolelh_thick)
groupmeans_thick$site_trvtmlh_thick_cent <- groupmeans_thick$site_trvtmlh_thick - mean(groupmeans_thick$site_trvtmlh_thick)
groupmeans_thick$site_insulalh_thick_cent <- groupmeans_thick$site_insulalh_thick - mean(groupmeans_thick$site_insulalh_thick)
groupmeans_thick$site_banksstsrh_thick_cent <- groupmeans_thick$site_banksstsrh_thick - mean(groupmeans_thick$site_banksstsrh_thick)
groupmeans_thick$site_cdacaterh_thick_cent <- groupmeans_thick$site_cdacaterh_thick - mean(groupmeans_thick$site_cdacaterh_thick)
groupmeans_thick$site_cdmdfrrh_thick_cent <- groupmeans_thick$site_cdmdfrrh_thick - mean(groupmeans_thick$site_cdmdfrrh_thick)
groupmeans_thick$site_cuneusrh_thick_cent <- groupmeans_thick$site_cuneusrh_thick - mean(groupmeans_thick$site_cuneusrh_thick)
groupmeans_thick$site_ehinalrh_thick_cent <- groupmeans_thick$site_ehinalrh_thick - mean(groupmeans_thick$site_ehinalrh_thick)
groupmeans_thick$site_fusiformrh_thick_cent <- groupmeans_thick$site_fusiformrh_thick - mean(groupmeans_thick$site_fusiformrh_thick)
groupmeans_thick$site_ifplrh_thick_cent <- groupmeans_thick$site_ifplrh_thick - mean(groupmeans_thick$site_ifplrh_thick)
groupmeans_thick$site_iftmrh_thick_cent <- groupmeans_thick$site_iftmrh_thick - mean(groupmeans_thick$site_iftmrh_thick)
groupmeans_thick$site_ihcaterh_thick_cent <- groupmeans_thick$site_ihcaterh_thick - mean(groupmeans_thick$site_ihcaterh_thick)
groupmeans_thick$site_loccrh_thick_cent <- groupmeans_thick$site_loccrh_thick - mean(groupmeans_thick$site_loccrh_thick)
groupmeans_thick$site_lobfrrh_thick_cent <- groupmeans_thick$site_lobfrrh_thick - mean(groupmeans_thick$site_lobfrrh_thick)
groupmeans_thick$site_lingualrh_thick_cent <- groupmeans_thick$site_lingualrh_thick - mean(groupmeans_thick$site_lingualrh_thick)
groupmeans_thick$site_mobfrrh_thick_cent <- groupmeans_thick$site_mobfrrh_thick - mean(groupmeans_thick$site_mobfrrh_thick)
groupmeans_thick$site_mdtmrh_thick_cent <- groupmeans_thick$site_mdtmrh_thick - mean(groupmeans_thick$site_mdtmrh_thick)
groupmeans_thick$site_parahpalrh_thick_cent <- groupmeans_thick$site_parahpalrh_thick - mean(groupmeans_thick$site_parahpalrh_thick)
groupmeans_thick$site_paracnrh_thick_cent <- groupmeans_thick$site_paracnrh_thick - mean(groupmeans_thick$site_paracnrh_thick)
groupmeans_thick$site_parsopcrh_thick_cent <- groupmeans_thick$site_parsopcrh_thick - mean(groupmeans_thick$site_parsopcrh_thick)
groupmeans_thick$site_parsobisrh_thick_cent <- groupmeans_thick$site_parsobisrh_thick - mean(groupmeans_thick$site_parsobisrh_thick)
groupmeans_thick$site_parstgrisrh_thick_cent <- groupmeans_thick$site_parstgrisrh_thick - mean(groupmeans_thick$site_parstgrisrh_thick)
groupmeans_thick$site_periccrh_thick_cent <- groupmeans_thick$site_periccrh_thick - mean(groupmeans_thick$site_periccrh_thick)
groupmeans_thick$site_postcnrh_thick_cent <- groupmeans_thick$site_postcnrh_thick - mean(groupmeans_thick$site_postcnrh_thick)
groupmeans_thick$site_ptcaterh_thick_cent <- groupmeans_thick$site_ptcaterh_thick - mean(groupmeans_thick$site_ptcaterh_thick)
groupmeans_thick$site_precnrh_thick_cent <- groupmeans_thick$site_precnrh_thick - mean(groupmeans_thick$site_precnrh_thick)
groupmeans_thick$site_pcrh_thick_cent <- groupmeans_thick$site_pcrh_thick - mean(groupmeans_thick$site_pcrh_thick)
groupmeans_thick$site_rracaterh_thick_cent <- groupmeans_thick$site_rracaterh_thick - mean(groupmeans_thick$site_rracaterh_thick)
groupmeans_thick$site_rrmdfrrh_thick_cent <- groupmeans_thick$site_rrmdfrrh_thick - mean(groupmeans_thick$site_rrmdfrrh_thick)
groupmeans_thick$site_sufrrh_thick_cent <- groupmeans_thick$site_sufrrh_thick - mean(groupmeans_thick$site_sufrrh_thick)
groupmeans_thick$site_suplrh_thick_cent <- groupmeans_thick$site_suplrh_thick - mean(groupmeans_thick$site_suplrh_thick)
groupmeans_thick$site_sutmrh_thick_cent <- groupmeans_thick$site_sutmrh_thick - mean(groupmeans_thick$site_sutmrh_thick)
groupmeans_thick$site_smrh_thick_cent <- groupmeans_thick$site_smrh_thick - mean(groupmeans_thick$site_smrh_thick)
groupmeans_thick$site_frpolerh_thick_cent <- groupmeans_thick$site_frpolerh_thick - mean(groupmeans_thick$site_frpolerh_thick)
groupmeans_thick$site_tmpolerh_thick_cent <- groupmeans_thick$site_tmpolerh_thick - mean(groupmeans_thick$site_tmpolerh_thick)
groupmeans_thick$site_trvtmrh_thick_cent <- groupmeans_thick$site_trvtmrh_thick - mean(groupmeans_thick$site_trvtmrh_thick)
groupmeans_thick$site_insularh_thick_cent <- groupmeans_thick$site_insularh_thick - mean(groupmeans_thick$site_insularh_thick)

groupmeans_area <- aggregate(cbind(banksstslh_area,cdacatelh_area,cdmdfrlh_area,cuneuslh_area,ehinallh_area,
                                   fusiformlh_area,ifpllh_area,iftmlh_area,ihcatelh_area,locclh_area,
                                   lobfrlh_area,linguallh_area,mobfrlh_area,mdtmlh_area,parahpallh_area,
                                   paracnlh_area,parsopclh_area,parsobislh_area,parstgrislh_area,
                                   pericclh_area,postcnlh_area,ptcatelh_area,precnlh_area,pclh_area,rracatelh_area,
                                   rrmdfrlh_area,sufrlh_area,supllh_area,sutmlh_area,smlh_area,frpolelh_area,
                                   tmpolelh_area,trvtmlh_area,insulalh_area,banksstsrh_area,cdacaterh_area,
                                   cdmdfrrh_area,cuneusrh_area,ehinalrh_area,fusiformrh_area,ifplrh_area,iftmrh_area,
                                   ihcaterh_area,loccrh_area,lobfrrh_area,lingualrh_area,mobfrrh_area,mdtmrh_area,
                                   parahpalrh_area,paracnrh_area,parsopcrh_area,parsobisrh_area,parstgrisrh_area,
                                   periccrh_area,postcnrh_area,ptcaterh_area,precnrh_area,pcrh_area,rracaterh_area,
                                   rrmdfrrh_area,sufrrh_area,suplrh_area,sutmrh_area,smrh_area,frpolerh_area,tmpolerh_area,
                                   trvtmrh_area,insularh_area)~site_id, pfactor_red, mean)
names(groupmeans_area) <- c("site_id","site_banksstslh_area","site_cdacatelh_area","site_cdmdfrlh_area","site_cuneuslh_area",
                            "site_ehinallh_area","site_fusiformlh_area","site_ifpllh_area","site_iftmlh_area","site_ihcatelh_area",
                            "site_locclh_area","site_lobfrlh_area","site_linguallh_area","site_mobfrlh_area","site_mdtmlh_area",
                            "site_parahpallh_area","site_paracnlh_area","site_parsopclh_area","site_parsobislh_area","site_parstgrislh_area",
                            "site_pericclh_area","site_postcnlh_area","site_ptcatelh_area","site_precnlh_area","site_pclh_area","site_rracatelh_area",
                            "site_rrmdfrlh_area","site_sufrlh_area","site_supllh_area","site_sutmlh_area","site_smlh_area","site_frpolelh_area",
                            "site_tmpolelh_area","site_trvtmlh_area","site_insulalh_area","site_banksstsrh_area","site_cdacaterh_area","site_cdmdfrrh_area",
                            "site_cuneusrh_area","site_ehinalrh_area","site_fusiformrh_area","site_ifplrh_area","site_iftmrh_area","site_ihcaterh_area",
                            "site_loccrh_area","site_lobfrrh_area","site_lingualrh_area","site_mobfrrh_area","site_mdtmrh_area",
                            "site_parahpalrh_area","site_paracnrh_area","site_parsopcrh_area","site_parsobisrh_area","site_parstgrisrh_area",
                            "site_periccrh_area","site_postcnrh_area","site_ptcaterh_area","site_precnrh_area","site_pcrh_area","site_rracaterh_area",
                            "site_rrmdfrrh_area","site_sufrrh_area","site_suplrh_area","site_sutmrh_area","site_smrh_area","site_frpolerh_area",
                            "site_tmpolerh_area","site_trvtmrh_area","site_insularh_area")
groupmeans_area$site_banksstslh_area_cent <- groupmeans_area$site_banksstslh_area - mean(groupmeans_area$site_banksstslh_area)
groupmeans_area$site_cdacatelh_area_cent <- groupmeans_area$site_cdacatelh_area - mean(groupmeans_area$site_cdacatelh_area)
groupmeans_area$site_cdmdfrlh_area_cent <- groupmeans_area$site_cdmdfrlh_area - mean(groupmeans_area$site_cdmdfrlh_area)
groupmeans_area$site_cuneuslh_area_cent <- groupmeans_area$site_cuneuslh_area - mean(groupmeans_area$site_cuneuslh_area)
groupmeans_area$site_ehinallh_area_cent <- groupmeans_area$site_ehinallh_area - mean(groupmeans_area$site_ehinallh_area)
groupmeans_area$site_fusiformlh_area_cent <- groupmeans_area$site_fusiformlh_area - mean(groupmeans_area$site_fusiformlh_area)
groupmeans_area$site_ifpllh_area_cent <- groupmeans_area$site_ifpllh_area - mean(groupmeans_area$site_ifpllh_area)
groupmeans_area$site_iftmlh_area_cent <- groupmeans_area$site_iftmlh_area - mean(groupmeans_area$site_iftmlh_area)
groupmeans_area$site_ihcatelh_area_cent <- groupmeans_area$site_ihcatelh_area - mean(groupmeans_area$site_ihcatelh_area)
groupmeans_area$site_locclh_area_cent <- groupmeans_area$site_locclh_area - mean(groupmeans_area$site_locclh_area)
groupmeans_area$site_lobfrlh_area_cent <- groupmeans_area$site_lobfrlh_area - mean(groupmeans_area$site_lobfrlh_area)
groupmeans_area$site_linguallh_area_cent <- groupmeans_area$site_linguallh_area - mean(groupmeans_area$site_linguallh_area)
groupmeans_area$site_mobfrlh_area_cent <- groupmeans_area$site_mobfrlh_area - mean(groupmeans_area$site_mobfrlh_area)
groupmeans_area$site_mdtmlh_area_cent <- groupmeans_area$site_mdtmlh_area - mean(groupmeans_area$site_mdtmlh_area)
groupmeans_area$site_parahpallh_area_cent <- groupmeans_area$site_parahpallh_area - mean(groupmeans_area$site_parahpallh_area)
groupmeans_area$site_paracnlh_area_cent <- groupmeans_area$site_paracnlh_area - mean(groupmeans_area$site_paracnlh_area)
groupmeans_area$site_parsopclh_area_cent <- groupmeans_area$site_parsopclh_area - mean(groupmeans_area$site_parsopclh_area)
groupmeans_area$site_parsobislh_area_cent <- groupmeans_area$site_parsobislh_area - mean(groupmeans_area$site_parsobislh_area)
groupmeans_area$site_parstgrislh_area_cent <- groupmeans_area$site_parstgrislh_area - mean(groupmeans_area$site_parstgrislh_area)
groupmeans_area$site_pericclh_area_cent <- groupmeans_area$site_pericclh_area - mean(groupmeans_area$site_pericclh_area)
groupmeans_area$site_postcnlh_area_cent <- groupmeans_area$site_postcnlh_area - mean(groupmeans_area$site_postcnlh_area)
groupmeans_area$site_ptcatelh_area_cent <- groupmeans_area$site_ptcatelh_area - mean(groupmeans_area$site_ptcatelh_area)
groupmeans_area$site_precnlh_area_cent <- groupmeans_area$site_precnlh_area - mean(groupmeans_area$site_precnlh_area)
groupmeans_area$site_pclh_area_cent <- groupmeans_area$site_pclh_area - mean(groupmeans_area$site_pclh_area)
groupmeans_area$site_rracatelh_area_cent <- groupmeans_area$site_rracatelh_area - mean(groupmeans_area$site_rracatelh_area)
groupmeans_area$site_rrmdfrlh_area_cent <- groupmeans_area$site_rrmdfrlh_area - mean(groupmeans_area$site_rrmdfrlh_area)
groupmeans_area$site_sufrlh_area_cent <- groupmeans_area$site_sufrlh_area - mean(groupmeans_area$site_sufrlh_area)
groupmeans_area$site_supllh_area_cent <- groupmeans_area$site_supllh_area - mean(groupmeans_area$site_supllh_area)
groupmeans_area$site_sutmlh_area_cent <- groupmeans_area$site_sutmlh_area - mean(groupmeans_area$site_sutmlh_area)
groupmeans_area$site_smlh_area_cent <- groupmeans_area$site_smlh_area - mean(groupmeans_area$site_smlh_area)
groupmeans_area$site_frpolelh_area_cent <- groupmeans_area$site_frpolelh_area - mean(groupmeans_area$site_frpolelh_area)
groupmeans_area$site_tmpolelh_area_cent <- groupmeans_area$site_tmpolelh_area - mean(groupmeans_area$site_tmpolelh_area)
groupmeans_area$site_trvtmlh_area_cent <- groupmeans_area$site_trvtmlh_area - mean(groupmeans_area$site_trvtmlh_area)
groupmeans_area$site_insulalh_area_cent <- groupmeans_area$site_insulalh_area - mean(groupmeans_area$site_insulalh_area)
groupmeans_area$site_banksstsrh_area_cent <- groupmeans_area$site_banksstsrh_area - mean(groupmeans_area$site_banksstsrh_area)
groupmeans_area$site_cdacaterh_area_cent <- groupmeans_area$site_cdacaterh_area - mean(groupmeans_area$site_cdacaterh_area)
groupmeans_area$site_cdmdfrrh_area_cent <- groupmeans_area$site_cdmdfrrh_area - mean(groupmeans_area$site_cdmdfrrh_area)
groupmeans_area$site_cuneusrh_area_cent <- groupmeans_area$site_cuneusrh_area - mean(groupmeans_area$site_cuneusrh_area)
groupmeans_area$site_ehinalrh_area_cent <- groupmeans_area$site_ehinalrh_area - mean(groupmeans_area$site_ehinalrh_area)
groupmeans_area$site_fusiformrh_area_cent <- groupmeans_area$site_fusiformrh_area - mean(groupmeans_area$site_fusiformrh_area)
groupmeans_area$site_ifplrh_area_cent <- groupmeans_area$site_ifplrh_area - mean(groupmeans_area$site_ifplrh_area)
groupmeans_area$site_iftmrh_area_cent <- groupmeans_area$site_iftmrh_area - mean(groupmeans_area$site_iftmrh_area)
groupmeans_area$site_ihcaterh_area_cent <- groupmeans_area$site_ihcaterh_area - mean(groupmeans_area$site_ihcaterh_area)
groupmeans_area$site_loccrh_area_cent <- groupmeans_area$site_loccrh_area - mean(groupmeans_area$site_loccrh_area)
groupmeans_area$site_lobfrrh_area_cent <- groupmeans_area$site_lobfrrh_area - mean(groupmeans_area$site_lobfrrh_area)
groupmeans_area$site_lingualrh_area_cent <- groupmeans_area$site_lingualrh_area - mean(groupmeans_area$site_lingualrh_area)
groupmeans_area$site_mobfrrh_area_cent <- groupmeans_area$site_mobfrrh_area - mean(groupmeans_area$site_mobfrrh_area)
groupmeans_area$site_mdtmrh_area_cent <- groupmeans_area$site_mdtmrh_area - mean(groupmeans_area$site_mdtmrh_area)
groupmeans_area$site_parahpalrh_area_cent <- groupmeans_area$site_parahpalrh_area - mean(groupmeans_area$site_parahpalrh_area)
groupmeans_area$site_paracnrh_area_cent <- groupmeans_area$site_paracnrh_area - mean(groupmeans_area$site_paracnrh_area)
groupmeans_area$site_parsopcrh_area_cent <- groupmeans_area$site_parsopcrh_area - mean(groupmeans_area$site_parsopcrh_area)
groupmeans_area$site_parsobisrh_area_cent <- groupmeans_area$site_parsobisrh_area - mean(groupmeans_area$site_parsobisrh_area)
groupmeans_area$site_parstgrisrh_area_cent <- groupmeans_area$site_parstgrisrh_area - mean(groupmeans_area$site_parstgrisrh_area)
groupmeans_area$site_periccrh_area_cent <- groupmeans_area$site_periccrh_area - mean(groupmeans_area$site_periccrh_area)
groupmeans_area$site_postcnrh_area_cent <- groupmeans_area$site_postcnrh_area - mean(groupmeans_area$site_postcnrh_area)
groupmeans_area$site_ptcaterh_area_cent <- groupmeans_area$site_ptcaterh_area - mean(groupmeans_area$site_ptcaterh_area)
groupmeans_area$site_precnrh_area_cent <- groupmeans_area$site_precnrh_area - mean(groupmeans_area$site_precnrh_area)
groupmeans_area$site_pcrh_area_cent <- groupmeans_area$site_pcrh_area - mean(groupmeans_area$site_pcrh_area)
groupmeans_area$site_rracaterh_area_cent <- groupmeans_area$site_rracaterh_area - mean(groupmeans_area$site_rracaterh_area)
groupmeans_area$site_rrmdfrrh_area_cent <- groupmeans_area$site_rrmdfrrh_area - mean(groupmeans_area$site_rrmdfrrh_area)
groupmeans_area$site_sufrrh_area_cent <- groupmeans_area$site_sufrrh_area - mean(groupmeans_area$site_sufrrh_area)
groupmeans_area$site_suplrh_area_cent <- groupmeans_area$site_suplrh_area - mean(groupmeans_area$site_suplrh_area)
groupmeans_area$site_sutmrh_area_cent <- groupmeans_area$site_sutmrh_area - mean(groupmeans_area$site_sutmrh_area)
groupmeans_area$site_smrh_area_cent <- groupmeans_area$site_smrh_area - mean(groupmeans_area$site_smrh_area)
groupmeans_area$site_frpolerh_area_cent <- groupmeans_area$site_frpolerh_area - mean(groupmeans_area$site_frpolerh_area)
groupmeans_area$site_tmpolerh_area_cent <- groupmeans_area$site_tmpolerh_area - mean(groupmeans_area$site_tmpolerh_area)
groupmeans_area$site_trvtmrh_area_cent <- groupmeans_area$site_trvtmrh_area - mean(groupmeans_area$site_trvtmrh_area)
groupmeans_area$site_insularh_area_cent <- groupmeans_area$site_insularh_area - mean(groupmeans_area$site_insularh_area)

#Centering parcel-wise subcortical volume predictors
groupmeans_subcort <- aggregate(cbind(crbcortexlh_vol,tplh_vol,caudatelh_vol,putamenlh_vol,pallidumlh_vol,bstem_vol,
                                      hpuslh_vol,amygdalalh_vol,aal_vol,vedclh_vol,crbcortexrh_vol,tprh_vol,caudaterh_vol,
                                      putamenrh_vol,pallidumrh_vol,hpusrh_vol,amygdalarh_vol,aar_vol,vedcrh_vol)~site_id, pfactor_red, mean)
names(groupmeans_subcort) <- c("site_id","site_crbcortexlh_vol","site_tplh_vol","site_caudatelh_vol","site_putamenlh_vol",
                               "site_pallidumlh_vol","site_bstem_vol","site_hpuslh_vol","site_amygdalalh_vol","site_aal_vol",
                               "site_vedclh_vol","site_crbcortexrh_vol","site_tprh_vol","site_caudaterh_vol","site_putamenrh_vol",
                               "site_pallidumrh_vol","site_hpusrh_vol","site_amygdalarh_vol","site_aar_vol","site_vedcrh_vol")
groupmeans_subcort$site_crbcortexlh_vol_cent <- groupmeans_subcort$site_crbcortexlh_vol - mean(groupmeans_subcort$site_crbcortexlh_vol)
groupmeans_subcort$site_tplh_vol_cent <- groupmeans_subcort$site_tplh_vol - mean(groupmeans_subcort$site_tplh_vol)
groupmeans_subcort$site_caudatelh_vol_cent <- groupmeans_subcort$site_caudatelh_vol - mean(groupmeans_subcort$site_caudatelh_vol)
groupmeans_subcort$site_putamenlh_vol_cent <- groupmeans_subcort$site_putamenlh_vol - mean(groupmeans_subcort$site_putamenlh_vol)
groupmeans_subcort$site_pallidumlh_vol_cent <- groupmeans_subcort$site_pallidumlh_vol - mean(groupmeans_subcort$site_pallidumlh_vol)
groupmeans_subcort$site_bstem_vol_cent <- groupmeans_subcort$site_bstem_vol - mean(groupmeans_subcort$site_bstem_vol)
groupmeans_subcort$site_hpuslh_vol_cent <- groupmeans_subcort$site_hpuslh_vol - mean(groupmeans_subcort$site_hpuslh_vol)
groupmeans_subcort$site_amygdalalh_vol_cent <- groupmeans_subcort$site_amygdalalh_vol - mean(groupmeans_subcort$site_amygdalalh_vol)
groupmeans_subcort$site_aal_vol_cent <- groupmeans_subcort$site_aal_vol - mean(groupmeans_subcort$site_aal_vol)
groupmeans_subcort$site_vedclh_vol_cent <- groupmeans_subcort$site_vedclh_vol - mean(groupmeans_subcort$site_vedclh_vol)
groupmeans_subcort$site_crbcortexrh_vol_cent <- groupmeans_subcort$site_crbcortexrh_vol - mean(groupmeans_subcort$site_crbcortexrh_vol)
groupmeans_subcort$site_tprh_vol_cent <- groupmeans_subcort$site_tprh_vol - mean(groupmeans_subcort$site_tprh_vol)
groupmeans_subcort$site_caudaterh_vol_cent <- groupmeans_subcort$site_caudaterh_vol - mean(groupmeans_subcort$site_caudaterh_vol)
groupmeans_subcort$site_putamenrh_vol_cent <- groupmeans_subcort$site_putamenrh_vol - mean(groupmeans_subcort$site_putamenrh_vol)
groupmeans_subcort$site_pallidumrh_vol_cent <- groupmeans_subcort$site_pallidumrh_vol - mean(groupmeans_subcort$site_pallidumrh_vol)
groupmeans_subcort$site_hpusrh_vol_cent <- groupmeans_subcort$site_hpusrh_vol - mean(groupmeans_subcort$site_hpusrh_vol)
groupmeans_subcort$site_amygdalarh_vol_cent <- groupmeans_subcort$site_amygdalarh_vol - mean(groupmeans_subcort$site_amygdalarh_vol)
groupmeans_subcort$site_aar_vol_cent <- groupmeans_subcort$site_aar_vol - mean(groupmeans_subcort$site_aar_vol)
groupmeans_subcort$site_vedcrh_vol_cent <- groupmeans_subcort$site_vedcrh_vol - mean(groupmeans_subcort$site_vedcrh_vol)

#join the site centered cortical and subcortical volume and cortical thickness and area parcels
pfactor_red <- join(pfactor_red, groupmeans_vol, by="site_id")
pfactor_red <- join(pfactor_red, groupmeans_subcort, by="site_id")
pfactor_red <- join(pfactor_red, groupmeans_thick, by="site_id")
pfactor_red <- join(pfactor_red, groupmeans_area, by="site_id")

#Create subject-centered volume variables
pfactor_red$subj_banksstslh_vol <- pfactor_red$banksstslh_vol - pfactor_red$site_banksstslh_vol
pfactor_red$subj_cdacatelh_vol <- pfactor_red$cdacatelh_vol - pfactor_red$site_cdacatelh_vol
pfactor_red$subj_cdmdfrlh_vol <- pfactor_red$cdmdfrlh_vol - pfactor_red$site_cdmdfrlh_vol
pfactor_red$subj_cuneuslh_vol <- pfactor_red$cuneuslh_vol - pfactor_red$site_cuneuslh_vol
pfactor_red$subj_ehinallh_vol <- pfactor_red$ehinallh_vol - pfactor_red$site_ehinallh_vol
pfactor_red$subj_fusiformlh_vol <- pfactor_red$fusiformlh_vol - pfactor_red$site_fusiformlh_vol
pfactor_red$subj_ifpllh_vol <- pfactor_red$ifpllh_vol - pfactor_red$site_ifpllh_vol
pfactor_red$subj_iftmlh_vol <- pfactor_red$iftmlh_vol - pfactor_red$site_iftmlh_vol
pfactor_red$subj_ihcatelh_vol <- pfactor_red$ihcatelh_vol - pfactor_red$site_ihcatelh_vol
pfactor_red$subj_locclh_vol <- pfactor_red$locclh_vol - pfactor_red$site_locclh_vol
pfactor_red$subj_lobfrlh_vol <- pfactor_red$lobfrlh_vol - pfactor_red$site_lobfrlh_vol
pfactor_red$subj_linguallh_vol <- pfactor_red$linguallh_vol - pfactor_red$site_linguallh_vol
pfactor_red$subj_mobfrlh_vol <- pfactor_red$mobfrlh_vol - pfactor_red$site_mobfrlh_vol
pfactor_red$subj_mdtmlh_vol <- pfactor_red$mdtmlh_vol - pfactor_red$site_mdtmlh_vol
pfactor_red$subj_parahpallh_vol <- pfactor_red$parahpallh_vol - pfactor_red$site_parahpallh_vol
pfactor_red$subj_paracnlh_vol <- pfactor_red$paracnlh_vol - pfactor_red$site_paracnlh_vol
pfactor_red$subj_parsopclh_vol <- pfactor_red$parsopclh_vol - pfactor_red$site_parsopclh_vol
pfactor_red$subj_parsobislh_vol <- pfactor_red$parsobislh_vol - pfactor_red$site_parsobislh_vol
pfactor_red$subj_parstgrislh_vol <- pfactor_red$parstgrislh_vol - pfactor_red$site_parstgrislh_vol
pfactor_red$subj_pericclh_vol <- pfactor_red$pericclh_vol - pfactor_red$site_pericclh_vol
pfactor_red$subj_postcnlh_vol <- pfactor_red$postcnlh_vol - pfactor_red$site_postcnlh_vol
pfactor_red$subj_ptcatelh_vol <- pfactor_red$ptcatelh_vol - pfactor_red$site_ptcatelh_vol
pfactor_red$subj_precnlh_vol <- pfactor_red$precnlh_vol - pfactor_red$site_precnlh_vol
pfactor_red$subj_pclh_vol <- pfactor_red$pclh_vol - pfactor_red$site_pclh_vol
pfactor_red$subj_rracatelh_vol <- pfactor_red$rracatelh_vol - pfactor_red$site_rracatelh_vol
pfactor_red$subj_rrmdfrlh_vol <- pfactor_red$rrmdfrlh_vol - pfactor_red$site_rrmdfrlh_vol
pfactor_red$subj_sufrlh_vol <- pfactor_red$sufrlh_vol - pfactor_red$site_sufrlh_vol
pfactor_red$subj_supllh_vol <- pfactor_red$supllh_vol - pfactor_red$site_supllh_vol
pfactor_red$subj_sutmlh_vol <- pfactor_red$sutmlh_vol - pfactor_red$site_sutmlh_vol
pfactor_red$subj_smlh_vol <- pfactor_red$smlh_vol - pfactor_red$site_smlh_vol
pfactor_red$subj_frpolelh_vol <- pfactor_red$frpolelh_vol - pfactor_red$site_frpolelh_vol
pfactor_red$subj_tmpolelh_vol <- pfactor_red$tmpolelh_vol - pfactor_red$site_tmpolelh_vol
pfactor_red$subj_trvtmlh_vol <- pfactor_red$trvtmlh_vol - pfactor_red$site_trvtmlh_vol
pfactor_red$subj_insulalh_vol <- pfactor_red$insulalh_vol - pfactor_red$site_insulalh_vol
pfactor_red$subj_banksstsrh_vol <- pfactor_red$banksstsrh_vol - pfactor_red$site_banksstsrh_vol
pfactor_red$subj_cdacaterh_vol <- pfactor_red$cdacaterh_vol - pfactor_red$site_cdacaterh_vol
pfactor_red$subj_cdmdfrrh_vol <- pfactor_red$cdmdfrrh_vol - pfactor_red$site_cdmdfrrh_vol
pfactor_red$subj_cuneusrh_vol <- pfactor_red$cuneusrh_vol - pfactor_red$site_cuneusrh_vol
pfactor_red$subj_ehinalrh_vol <- pfactor_red$ehinalrh_vol - pfactor_red$site_ehinalrh_vol
pfactor_red$subj_fusiformrh_vol <- pfactor_red$fusiformrh_vol - pfactor_red$site_fusiformrh_vol
pfactor_red$subj_ifplrh_vol <- pfactor_red$ifplrh_vol - pfactor_red$site_ifplrh_vol
pfactor_red$subj_iftmrh_vol <- pfactor_red$iftmrh_vol - pfactor_red$site_iftmrh_vol
pfactor_red$subj_ihcaterh_vol <- pfactor_red$ihcaterh_vol - pfactor_red$site_ihcaterh_vol
pfactor_red$subj_loccrh_vol <- pfactor_red$loccrh_vol - pfactor_red$site_loccrh_vol
pfactor_red$subj_lobfrrh_vol <- pfactor_red$lobfrrh_vol - pfactor_red$site_lobfrrh_vol
pfactor_red$subj_lingualrh_vol <- pfactor_red$lingualrh_vol - pfactor_red$site_lingualrh_vol
pfactor_red$subj_mobfrrh_vol <- pfactor_red$mobfrrh_vol - pfactor_red$site_mobfrrh_vol
pfactor_red$subj_mdtmrh_vol <- pfactor_red$mdtmrh_vol - pfactor_red$site_mdtmrh_vol
pfactor_red$subj_parahpalrh_vol <- pfactor_red$parahpalrh_vol - pfactor_red$site_parahpalrh_vol
pfactor_red$subj_paracnrh_vol <- pfactor_red$paracnrh_vol - pfactor_red$site_paracnrh_vol
pfactor_red$subj_parsopcrh_vol <- pfactor_red$parsopcrh_vol - pfactor_red$site_parsopcrh_vol
pfactor_red$subj_parsobisrh_vol <- pfactor_red$parsobisrh_vol - pfactor_red$site_parsobisrh_vol
pfactor_red$subj_parstgrisrh_vol <- pfactor_red$parstgrisrh_vol - pfactor_red$site_parstgrisrh_vol
pfactor_red$subj_periccrh_vol <- pfactor_red$periccrh_vol - pfactor_red$site_periccrh_vol 
pfactor_red$subj_postcnrh_vol <- pfactor_red$postcnrh_vol - pfactor_red$site_postcnrh_vol
pfactor_red$subj_ptcaterh_vol <- pfactor_red$ptcaterh_vol - pfactor_red$site_ptcaterh_vol
pfactor_red$subj_precnrh_vol <- pfactor_red$precnrh_vol - pfactor_red$site_precnrh_vol
pfactor_red$subj_pcrh_vol <- pfactor_red$pcrh_vol - pfactor_red$site_pcrh_vol
pfactor_red$subj_rracaterh_vol <- pfactor_red$rracaterh_vol - pfactor_red$site_rracaterh_vol
pfactor_red$subj_rrmdfrrh_vol <- pfactor_red$rrmdfrrh_vol - pfactor_red$site_rrmdfrrh_vol
pfactor_red$subj_sufrrh_vol <- pfactor_red$sufrrh_vol - pfactor_red$site_sufrrh_vol
pfactor_red$subj_suplrh_vol <- pfactor_red$suplrh_vol - pfactor_red$site_suplrh_vol
pfactor_red$subj_sutmrh_vol <- pfactor_red$sutmrh_vol - pfactor_red$site_sutmrh_vol
pfactor_red$subj_smrh_vol <- pfactor_red$smrh_vol - pfactor_red$site_smrh_vol
pfactor_red$subj_frpolerh_vol <- pfactor_red$smrh_vol - pfactor_red$site_smrh_vol
pfactor_red$subj_tmpolerh_vol <- pfactor_red$tmpolerh_vol - pfactor_red$site_tmpolerh_vol
pfactor_red$subj_trvtmrh_vol <- pfactor_red$trvtmrh_vol - pfactor_red$site_trvtmrh_vol
pfactor_red$subj_insularh_vol <- pfactor_red$insularh_vol - pfactor_red$site_insularh_vol

#Create subject-centered thickness variables
pfactor_red$subj_banksstslh_thick <- pfactor_red$banksstslh_thick - pfactor_red$site_banksstslh_thick
pfactor_red$subj_cdacatelh_thick <- pfactor_red$cdacatelh_thick - pfactor_red$site_cdacatelh_thick
pfactor_red$subj_cdmdfrlh_thick <- pfactor_red$cdmdfrlh_thick - pfactor_red$site_cdmdfrlh_thick
pfactor_red$subj_cuneuslh_thick <- pfactor_red$cuneuslh_thick - pfactor_red$site_cuneuslh_thick
pfactor_red$subj_ehinallh_thick <- pfactor_red$ehinallh_thick - pfactor_red$site_ehinallh_thick
pfactor_red$subj_fusiformlh_thick <- pfactor_red$fusiformlh_thick - pfactor_red$site_fusiformlh_thick
pfactor_red$subj_ifpllh_thick <- pfactor_red$ifpllh_thick - pfactor_red$site_ifpllh_thick
pfactor_red$subj_iftmlh_thick <- pfactor_red$iftmlh_thick - pfactor_red$site_iftmlh_thick
pfactor_red$subj_ihcatelh_thick <- pfactor_red$ihcatelh_thick - pfactor_red$site_ihcatelh_thick
pfactor_red$subj_locclh_thick <- pfactor_red$locclh_thick - pfactor_red$site_locclh_thick
pfactor_red$subj_lobfrlh_thick <- pfactor_red$lobfrlh_thick - pfactor_red$site_lobfrlh_thick
pfactor_red$subj_linguallh_thick <- pfactor_red$linguallh_thick - pfactor_red$site_linguallh_thick
pfactor_red$subj_mobfrlh_thick <- pfactor_red$mobfrlh_thick - pfactor_red$site_mobfrlh_thick
pfactor_red$subj_mdtmlh_thick <- pfactor_red$mdtmlh_thick - pfactor_red$site_mdtmlh_thick
pfactor_red$subj_parahpallh_thick <- pfactor_red$parahpallh_thick - pfactor_red$site_parahpallh_thick
pfactor_red$subj_paracnlh_thick <- pfactor_red$paracnlh_thick - pfactor_red$site_paracnlh_thick
pfactor_red$subj_parsopclh_thick <- pfactor_red$parsopclh_thick - pfactor_red$site_parsopclh_thick
pfactor_red$subj_parsobislh_thick <- pfactor_red$parsobislh_thick - pfactor_red$site_parsobislh_thick
pfactor_red$subj_parstgrislh_thick <- pfactor_red$parstgrislh_thick - pfactor_red$site_parstgrislh_thick
pfactor_red$subj_pericclh_thick <- pfactor_red$pericclh_thick - pfactor_red$site_pericclh_thick
pfactor_red$subj_postcnlh_thick <- pfactor_red$postcnlh_thick - pfactor_red$site_postcnlh_thick
pfactor_red$subj_ptcatelh_thick <- pfactor_red$ptcatelh_thick - pfactor_red$site_ptcatelh_thick
pfactor_red$subj_precnlh_thick <- pfactor_red$precnlh_thick - pfactor_red$site_precnlh_thick
pfactor_red$subj_pclh_thick <- pfactor_red$pclh_thick - pfactor_red$site_pclh_thick
pfactor_red$subj_rracatelh_thick <- pfactor_red$rracatelh_thick - pfactor_red$site_rracatelh_thick
pfactor_red$subj_rrmdfrlh_thick <- pfactor_red$rrmdfrlh_thick - pfactor_red$site_rrmdfrlh_thick
pfactor_red$subj_sufrlh_thick <- pfactor_red$sufrlh_thick - pfactor_red$site_sufrlh_thick
pfactor_red$subj_supllh_thick <- pfactor_red$supllh_thick - pfactor_red$site_supllh_thick
pfactor_red$subj_sutmlh_thick <- pfactor_red$sutmlh_thick - pfactor_red$site_sutmlh_thick
pfactor_red$subj_smlh_thick <- pfactor_red$smlh_thick - pfactor_red$site_smlh_thick
pfactor_red$subj_frpolelh_thick <- pfactor_red$frpolelh_thick - pfactor_red$site_frpolelh_thick
pfactor_red$subj_tmpolelh_thick <- pfactor_red$tmpolelh_thick - pfactor_red$site_tmpolelh_thick
pfactor_red$subj_trvtmlh_thick <- pfactor_red$trvtmlh_thick - pfactor_red$site_trvtmlh_thick
pfactor_red$subj_insulalh_thick <- pfactor_red$insulalh_thick - pfactor_red$site_insulalh_thick
pfactor_red$subj_banksstsrh_thick <- pfactor_red$banksstsrh_thick - pfactor_red$site_banksstsrh_thick
pfactor_red$subj_cdacaterh_thick <- pfactor_red$cdacaterh_thick - pfactor_red$site_cdacaterh_thick
pfactor_red$subj_cdmdfrrh_thick <- pfactor_red$cdmdfrrh_thick - pfactor_red$site_cdmdfrrh_thick
pfactor_red$subj_cuneusrh_thick <- pfactor_red$cuneusrh_thick - pfactor_red$site_cuneusrh_thick
pfactor_red$subj_ehinalrh_thick <- pfactor_red$ehinalrh_thick - pfactor_red$site_ehinalrh_thick
pfactor_red$subj_fusiformrh_thick <- pfactor_red$fusiformrh_thick - pfactor_red$site_fusiformrh_thick
pfactor_red$subj_ifplrh_thick <- pfactor_red$ifplrh_thick - pfactor_red$site_ifplrh_thick
pfactor_red$subj_iftmrh_thick <- pfactor_red$iftmrh_thick - pfactor_red$site_iftmrh_thick
pfactor_red$subj_ihcaterh_thick <- pfactor_red$ihcaterh_thick - pfactor_red$site_ihcaterh_thick
pfactor_red$subj_loccrh_thick <- pfactor_red$loccrh_thick - pfactor_red$site_loccrh_thick
pfactor_red$subj_lobfrrh_thick <- pfactor_red$lobfrrh_thick - pfactor_red$site_lobfrrh_thick
pfactor_red$subj_lingualrh_thick <- pfactor_red$lingualrh_thick - pfactor_red$site_lingualrh_thick
pfactor_red$subj_mobfrrh_thick <- pfactor_red$mobfrrh_thick - pfactor_red$site_mobfrrh_thick
pfactor_red$subj_mdtmrh_thick <- pfactor_red$mdtmrh_thick - pfactor_red$site_mdtmrh_thick
pfactor_red$subj_parahpalrh_thick <- pfactor_red$parahpalrh_thick - pfactor_red$site_parahpalrh_thick
pfactor_red$subj_paracnrh_thick <- pfactor_red$paracnrh_thick - pfactor_red$site_paracnrh_thick
pfactor_red$subj_parsopcrh_thick <- pfactor_red$parsopcrh_thick - pfactor_red$site_parsopcrh_thick
pfactor_red$subj_parsobisrh_thick <- pfactor_red$parsobisrh_thick - pfactor_red$site_parsobisrh_thick
pfactor_red$subj_parstgrisrh_thick <- pfactor_red$parstgrisrh_thick - pfactor_red$site_parstgrisrh_thick
pfactor_red$subj_periccrh_thick <- pfactor_red$periccrh_thick - pfactor_red$site_periccrh_thick
pfactor_red$subj_postcnrh_thick <- pfactor_red$postcnrh_thick - pfactor_red$site_postcnrh_thick
pfactor_red$subj_ptcaterh_thick <- pfactor_red$ptcaterh_thick - pfactor_red$site_ptcaterh_thick
pfactor_red$subj_precnrh_thick <- pfactor_red$precnrh_thick - pfactor_red$site_precnrh_thick
pfactor_red$subj_pcrh_thick <- pfactor_red$pcrh_thick - pfactor_red$site_pcrh_thick
pfactor_red$subj_rracaterh_thick <- pfactor_red$rracaterh_thick - pfactor_red$site_rracaterh_thick
pfactor_red$subj_rrmdfrrh_thick <- pfactor_red$rrmdfrrh_thick - pfactor_red$site_rrmdfrrh_thick
pfactor_red$subj_sufrrh_thick <- pfactor_red$sufrrh_thick - pfactor_red$site_sufrrh_thick
pfactor_red$subj_suplrh_thick <- pfactor_red$suplrh_thick - pfactor_red$site_suplrh_thick
pfactor_red$subj_sutmrh_thick <- pfactor_red$sutmrh_thick - pfactor_red$site_sutmrh_thick
pfactor_red$subj_smrh_thick <- pfactor_red$smrh_thick - pfactor_red$site_smrh_thick
pfactor_red$subj_frpolerh_thick <- pfactor_red$smrh_thick - pfactor_red$site_smrh_thick
pfactor_red$subj_tmpolerh_thick <- pfactor_red$tmpolerh_thick - pfactor_red$site_tmpolerh_thick
pfactor_red$subj_trvtmrh_thick <- pfactor_red$trvtmrh_thick - pfactor_red$site_trvtmrh_thick
pfactor_red$subj_insularh_thick <- pfactor_red$insularh_thick - pfactor_red$site_insularh_thick

pfactor_red$subj_banksstslh_area <- pfactor_red$banksstslh_area - pfactor_red$site_banksstslh_area
pfactor_red$subj_cdacatelh_area <- pfactor_red$cdacatelh_area - pfactor_red$site_cdacatelh_area
pfactor_red$subj_cdmdfrlh_area <- pfactor_red$cdmdfrlh_area - pfactor_red$site_cdmdfrlh_area
pfactor_red$subj_cuneuslh_area <- pfactor_red$cuneuslh_area - pfactor_red$site_cuneuslh_area
pfactor_red$subj_ehinallh_area <- pfactor_red$ehinallh_area - pfactor_red$site_ehinallh_area
pfactor_red$subj_fusiformlh_area <- pfactor_red$fusiformlh_area - pfactor_red$site_fusiformlh_area
pfactor_red$subj_ifpllh_area <- pfactor_red$ifpllh_area - pfactor_red$site_ifpllh_area
pfactor_red$subj_iftmlh_area <- pfactor_red$iftmlh_area - pfactor_red$site_iftmlh_area
pfactor_red$subj_ihcatelh_area <- pfactor_red$ihcatelh_area - pfactor_red$site_ihcatelh_area
pfactor_red$subj_locclh_area <- pfactor_red$locclh_area - pfactor_red$site_locclh_area
pfactor_red$subj_lobfrlh_area <- pfactor_red$lobfrlh_area - pfactor_red$site_lobfrlh_area
pfactor_red$subj_linguallh_area <- pfactor_red$linguallh_area - pfactor_red$site_linguallh_area
pfactor_red$subj_mobfrlh_area <- pfactor_red$mobfrlh_area - pfactor_red$site_mobfrlh_area
pfactor_red$subj_mdtmlh_area <- pfactor_red$mdtmlh_area - pfactor_red$site_mdtmlh_area
pfactor_red$subj_parahpallh_area <- pfactor_red$parahpallh_area - pfactor_red$site_parahpallh_area
pfactor_red$subj_paracnlh_area <- pfactor_red$paracnlh_area - pfactor_red$site_paracnlh_area
pfactor_red$subj_parsopclh_area <- pfactor_red$parsopclh_area - pfactor_red$site_parsopclh_area
pfactor_red$subj_parsobislh_area <- pfactor_red$parsobislh_area - pfactor_red$site_parsobislh_area
pfactor_red$subj_parstgrislh_area <- pfactor_red$parstgrislh_area - pfactor_red$site_parstgrislh_area
pfactor_red$subj_pericclh_area <- pfactor_red$pericclh_area - pfactor_red$site_pericclh_area
pfactor_red$subj_postcnlh_area <- pfactor_red$postcnlh_area - pfactor_red$site_postcnlh_area
pfactor_red$subj_ptcatelh_area <- pfactor_red$ptcatelh_area - pfactor_red$site_ptcatelh_area
pfactor_red$subj_precnlh_area <- pfactor_red$precnlh_area - pfactor_red$site_precnlh_area
pfactor_red$subj_pclh_area <- pfactor_red$pclh_area - pfactor_red$site_pclh_area
pfactor_red$subj_rracatelh_area <- pfactor_red$rracatelh_area - pfactor_red$site_rracatelh_area
pfactor_red$subj_rrmdfrlh_area <- pfactor_red$rrmdfrlh_area - pfactor_red$site_rrmdfrlh_area
pfactor_red$subj_sufrlh_area <- pfactor_red$sufrlh_area - pfactor_red$site_sufrlh_area
pfactor_red$subj_supllh_area <- pfactor_red$supllh_area - pfactor_red$site_supllh_area
pfactor_red$subj_sutmlh_area <- pfactor_red$sutmlh_area - pfactor_red$site_sutmlh_area
pfactor_red$subj_smlh_area <- pfactor_red$smlh_area - pfactor_red$site_smlh_area
pfactor_red$subj_frpolelh_area <- pfactor_red$frpolelh_area - pfactor_red$site_frpolelh_area
pfactor_red$subj_tmpolelh_area <- pfactor_red$tmpolelh_area - pfactor_red$site_tmpolelh_area
pfactor_red$subj_trvtmlh_area <- pfactor_red$trvtmlh_area - pfactor_red$site_trvtmlh_area
pfactor_red$subj_insulalh_area <- pfactor_red$insulalh_area - pfactor_red$site_insulalh_area
pfactor_red$subj_banksstsrh_area <- pfactor_red$banksstsrh_area - pfactor_red$site_banksstsrh_area
pfactor_red$subj_cdacaterh_area <- pfactor_red$cdacaterh_area - pfactor_red$site_cdacaterh_area
pfactor_red$subj_cdmdfrrh_area <- pfactor_red$cdmdfrrh_area - pfactor_red$site_cdmdfrrh_area
pfactor_red$subj_cuneusrh_area <- pfactor_red$cuneusrh_area - pfactor_red$site_cuneusrh_area
pfactor_red$subj_ehinalrh_area <- pfactor_red$ehinalrh_area - pfactor_red$site_ehinalrh_area
pfactor_red$subj_fusiformrh_area <- pfactor_red$fusiformrh_area - pfactor_red$site_fusiformrh_area
pfactor_red$subj_ifplrh_area <- pfactor_red$ifplrh_area - pfactor_red$site_ifplrh_area
pfactor_red$subj_iftmrh_area <- pfactor_red$iftmrh_area - pfactor_red$site_iftmrh_area
pfactor_red$subj_ihcaterh_area <- pfactor_red$ihcaterh_area - pfactor_red$site_ihcaterh_area
pfactor_red$subj_loccrh_area <- pfactor_red$loccrh_area - pfactor_red$site_loccrh_area
pfactor_red$subj_lobfrrh_area <- pfactor_red$lobfrrh_area - pfactor_red$site_lobfrrh_area
pfactor_red$subj_lingualrh_area <- pfactor_red$lingualrh_area - pfactor_red$site_lingualrh_area
pfactor_red$subj_mobfrrh_area <- pfactor_red$mobfrrh_area - pfactor_red$site_mobfrrh_area
pfactor_red$subj_mdtmrh_area <- pfactor_red$mdtmrh_area - pfactor_red$site_mdtmrh_area
pfactor_red$subj_parahpalrh_area <- pfactor_red$parahpalrh_area - pfactor_red$site_parahpalrh_area
pfactor_red$subj_paracnrh_area <- pfactor_red$paracnrh_area - pfactor_red$site_paracnrh_area
pfactor_red$subj_parsopcrh_area <- pfactor_red$parsopcrh_area - pfactor_red$site_parsopcrh_area
pfactor_red$subj_parsobisrh_area <- pfactor_red$parsobisrh_area - pfactor_red$site_parsobisrh_area
pfactor_red$subj_parstgrisrh_area <- pfactor_red$parstgrisrh_area - pfactor_red$site_parstgrisrh_area
pfactor_red$subj_periccrh_area <- pfactor_red$periccrh_area - pfactor_red$site_periccrh_area
pfactor_red$subj_postcnrh_area <- pfactor_red$postcnrh_area - pfactor_red$site_postcnrh_area
pfactor_red$subj_ptcaterh_area <- pfactor_red$ptcaterh_area - pfactor_red$site_ptcaterh_area
pfactor_red$subj_precnrh_area <- pfactor_red$precnrh_area - pfactor_red$site_precnrh_area
pfactor_red$subj_pcrh_area <- pfactor_red$pcrh_area - pfactor_red$site_pcrh_area
pfactor_red$subj_rracaterh_area <- pfactor_red$rracaterh_area - pfactor_red$site_rracaterh_area
pfactor_red$subj_rrmdfrrh_area <- pfactor_red$rrmdfrrh_area - pfactor_red$site_rrmdfrrh_area
pfactor_red$subj_sufrrh_area <- pfactor_red$sufrrh_area - pfactor_red$site_sufrrh_area
pfactor_red$subj_suplrh_area <- pfactor_red$suplrh_area - pfactor_red$site_suplrh_area
pfactor_red$subj_sutmrh_area <- pfactor_red$sutmrh_area - pfactor_red$site_sutmrh_area
pfactor_red$subj_smrh_area <- pfactor_red$smrh_area - pfactor_red$site_smrh_area
pfactor_red$subj_frpolerh_area <- pfactor_red$smrh_area - pfactor_red$site_smrh_area
pfactor_red$subj_tmpolerh_area <- pfactor_red$tmpolerh_area - pfactor_red$site_tmpolerh_area
pfactor_red$subj_trvtmrh_area <- pfactor_red$trvtmrh_area - pfactor_red$site_trvtmrh_area
pfactor_red$subj_insularh_area <- pfactor_red$insularh_area - pfactor_red$site_insularh_area

#Create subject-centered subcortical volume variables
pfactor_red$subj_crbcortexlh_vol <- pfactor_red$crbcortexlh_vol - pfactor_red$site_crbcortexlh_vol
pfactor_red$subj_tplh_vol <- pfactor_red$tplh_vol - pfactor_red$site_tplh_vol
pfactor_red$subj_caudatelh_vol <- pfactor_red$caudatelh_vol - pfactor_red$site_caudatelh_vol
pfactor_red$subj_putamenlh_vol <- pfactor_red$putamenlh_vol - pfactor_red$site_putamenlh_vol
pfactor_red$subj_pallidumlh_vol <- pfactor_red$pallidumlh_vol - pfactor_red$site_pallidumlh_vol
pfactor_red$subj_bstem_vol <- pfactor_red$bstem_vol - pfactor_red$site_bstem_vol
pfactor_red$subj_hpuslh_vol <- pfactor_red$hpuslh_vol - pfactor_red$site_hpuslh_vol
pfactor_red$subj_amygdalalh_vol <- pfactor_red$amygdalalh_vol - pfactor_red$site_amygdalalh_vol
pfactor_red$subj_aal_vol <- pfactor_red$aal_vol - pfactor_red$site_aal_vol
pfactor_red$subj_vedclh_vol <- pfactor_red$vedclh_vol - pfactor_red$site_vedclh_vol
pfactor_red$subj_crbcortexrh_vol <- pfactor_red$crbcortexrh_vol - pfactor_red$site_crbcortexrh_vol
pfactor_red$subj_tprh_vol <- pfactor_red$tprh_vol - pfactor_red$site_tprh_vol
pfactor_red$subj_caudaterh_vol <- pfactor_red$caudaterh_vol - pfactor_red$site_caudaterh_vol
pfactor_red$subj_putamenrh_vol <- pfactor_red$putamenrh_vol - pfactor_red$site_putamenrh_vol
pfactor_red$subj_pallidumrh_vol <- pfactor_red$pallidumrh_vol - pfactor_red$site_pallidumrh_vol
pfactor_red$subj_hpusrh_vol <- pfactor_red$hpusrh_vol - pfactor_red$site_hpusrh_vol
pfactor_red$subj_amygdalarh_vol <- pfactor_red$amygdalarh_vol - pfactor_red$site_amygdalarh_vol
pfactor_red$subj_aar_vol <- pfactor_red$aar_vol - pfactor_red$site_aar_vol
pfactor_red$subj_vedcrh_vol <- pfactor_red$vedcrh_vol - pfactor_red$site_vedcrh_vol

write.csv(pfactor_red, "abcd_mlm_centered_vars_long_rel4_FINAL.csv") #write csv file


#Parcel-wise cortical volume analyses (68 cortical volume and area parcels, 19 subcortical volume parcels)
#conditional three-level growth model with parcel-wise structure measures predicting the psychopathology factor scores (looped) 
#covariates: sex, age, race/ethnicity dummies, scanner model dummies
#standardized betas, standard errors, and p-values for intercepts (i.e., intb, intse, intp) and slopes (i.e., slb, slse, slp) were saved in a csv file

vol1 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_banksstslh_vol_cent + scale(subj_banksstslh_vol) + site_banksstslh_vol_cent*wave + scale(subj_banksstslh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b1i <- as.data.frame(c(summary(vol1)$coefficients[14,1]))
vol_se1i <- as.data.frame(c(summary(vol1)$coefficients[14,2]))
vol_p1i <- as.data.frame(c(summary(vol1)$coefficients[14,5]))
vol_b1s <- as.data.frame(c(summary(vol1)$coefficients[16,1]))
vol_se1s <- as.data.frame(c(summary(vol1)$coefficients[16,2]))
vol_p1s <- as.data.frame(c(summary(vol1)$coefficients[16,5]))

vol2 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_cdacatelh_vol_cent + scale(subj_cdacatelh_vol) + site_cdacatelh_vol_cent*wave + scale(subj_cdacatelh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b2i <- as.data.frame(c(summary(vol2)$coefficients[14,1]))
vol_se2i <- as.data.frame(c(summary(vol2)$coefficients[14,2]))
vol_p2i <- as.data.frame(c(summary(vol2)$coefficients[14,5]))
vol_b2s <- as.data.frame(c(summary(vol2)$coefficients[16,1]))
vol_se2s <- as.data.frame(c(summary(vol2)$coefficients[16,2]))
vol_p2s <- as.data.frame(c(summary(vol2)$coefficients[16,5]))

vol3 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_cdmdfrlh_vol_cent + scale(subj_cdmdfrlh_vol) + site_cdmdfrlh_vol_cent*wave + scale(subj_cdmdfrlh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b3i <- as.data.frame(c(summary(vol3)$coefficients[14,1]))
vol_se3i <- as.data.frame(c(summary(vol3)$coefficients[14,2]))
vol_p3i <- as.data.frame(c(summary(vol3)$coefficients[14,5]))
vol_b3s <- as.data.frame(c(summary(vol3)$coefficients[16,1]))
vol_se3s <- as.data.frame(c(summary(vol3)$coefficients[16,2]))
vol_p3s <- as.data.frame(c(summary(vol3)$coefficients[16,5]))

vol4 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_cuneuslh_vol_cent + scale(subj_cuneuslh_vol) + site_cuneuslh_vol_cent*wave + scale(subj_cuneuslh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b4i <- as.data.frame(c(summary(vol4)$coefficients[14,1]))
vol_se4i <- as.data.frame(c(summary(vol4)$coefficients[14,2]))
vol_p4i <- as.data.frame(c(summary(vol4)$coefficients[14,5]))
vol_b4s <- as.data.frame(c(summary(vol4)$coefficients[16,1]))
vol_se4s <- as.data.frame(c(summary(vol4)$coefficients[16,2]))
vol_p4s <- as.data.frame(c(summary(vol4)$coefficients[16,5]))

vol5 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_ehinallh_vol_cent + scale(subj_ehinallh_vol) + site_ehinallh_vol_cent*wave + scale(subj_ehinallh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b5i <- as.data.frame(c(summary(vol5)$coefficients[14,1]))
vol_se5i <- as.data.frame(c(summary(vol5)$coefficients[14,2]))
vol_p5i <- as.data.frame(c(summary(vol5)$coefficients[14,5]))
vol_b5s <- as.data.frame(c(summary(vol5)$coefficients[16,1]))
vol_se5s <- as.data.frame(c(summary(vol5)$coefficients[16,2]))
vol_p5s <- as.data.frame(c(summary(vol5)$coefficients[16,5]))

vol6 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_fusiformlh_vol_cent + scale(subj_fusiformlh_vol) + site_fusiformlh_vol_cent*wave + scale(subj_fusiformlh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b6i <- as.data.frame(c(summary(vol6)$coefficients[14,1]))
vol_se6i <- as.data.frame(c(summary(vol6)$coefficients[14,2]))
vol_p6i <- as.data.frame(c(summary(vol6)$coefficients[14,5]))
vol_b6s <- as.data.frame(c(summary(vol6)$coefficients[16,1]))
vol_se6s <- as.data.frame(c(summary(vol6)$coefficients[16,2]))
vol_p6s <- as.data.frame(c(summary(vol6)$coefficients[16,5]))

vol7 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_ifpllh_vol_cent + scale(subj_ifpllh_vol) + site_ifpllh_vol_cent*wave + scale(subj_ifpllh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b7i <- as.data.frame(c(summary(vol7)$coefficients[14,1]))
vol_se7i <- as.data.frame(c(summary(vol7)$coefficients[14,2]))
vol_p7i <- as.data.frame(c(summary(vol7)$coefficients[14,5]))
vol_b7s <- as.data.frame(c(summary(vol7)$coefficients[16,1]))
vol_se7s <- as.data.frame(c(summary(vol7)$coefficients[16,2]))
vol_p7s <- as.data.frame(c(summary(vol7)$coefficients[16,5]))

vol8 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_iftmlh_vol_cent + scale(subj_iftmlh_vol) + site_iftmlh_vol_cent*wave + scale(subj_iftmlh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b8i <- as.data.frame(c(summary(vol8)$coefficients[14,1]))
vol_se8i <- as.data.frame(c(summary(vol8)$coefficients[14,2]))
vol_p8i <- as.data.frame(c(summary(vol8)$coefficients[14,5]))
vol_b8s <- as.data.frame(c(summary(vol8)$coefficients[16,1]))
vol_se8s <- as.data.frame(c(summary(vol8)$coefficients[16,2]))
vol_p8s <- as.data.frame(c(summary(vol8)$coefficients[16,5]))

vol9 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_ihcatelh_vol_cent + scale(subj_ihcatelh_vol) + site_ihcatelh_vol_cent*wave + scale(subj_ihcatelh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b9i <- as.data.frame(c(summary(vol9)$coefficients[14,1]))
vol_se9i <- as.data.frame(c(summary(vol9)$coefficients[14,2]))
vol_p9i <- as.data.frame(c(summary(vol9)$coefficients[14,5]))
vol_b9s <- as.data.frame(c(summary(vol9)$coefficients[16,1]))
vol_se9s <- as.data.frame(c(summary(vol9)$coefficients[16,2]))
vol_p9s <- as.data.frame(c(summary(vol9)$coefficients[16,5]))

vol10 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_locclh_vol_cent + scale(subj_locclh_vol) + site_locclh_vol_cent*wave + scale(subj_locclh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b10i <- as.data.frame(c(summary(vol10)$coefficients[14,1]))
vol_se10i <- as.data.frame(c(summary(vol10)$coefficients[14,2]))
vol_p10i <- as.data.frame(c(summary(vol10)$coefficients[14,5]))
vol_b10s <- as.data.frame(c(summary(vol10)$coefficients[16,1]))
vol_se10s <- as.data.frame(c(summary(vol10)$coefficients[16,2]))
vol_p10s <- as.data.frame(c(summary(vol10)$coefficients[16,5]))

vol11 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_lobfrlh_vol_cent + scale(subj_lobfrlh_vol) + site_lobfrlh_vol_cent*wave + scale(subj_lobfrlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b11i <- as.data.frame(c(summary(vol11)$coefficients[14,1]))
vol_se11i <- as.data.frame(c(summary(vol11)$coefficients[14,2]))
vol_p11i <- as.data.frame(c(summary(vol11)$coefficients[14,5]))
vol_b11s <- as.data.frame(c(summary(vol11)$coefficients[16,1]))
vol_se11s <- as.data.frame(c(summary(vol11)$coefficients[16,2]))
vol_p11s <- as.data.frame(c(summary(vol11)$coefficients[16,5]))

vol12 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_linguallh_vol_cent + scale(subj_linguallh_vol) + site_linguallh_vol_cent*wave + scale(subj_linguallh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b12i <- as.data.frame(c(summary(vol12)$coefficients[14,1]))
vol_se12i <- as.data.frame(c(summary(vol12)$coefficients[14,2]))
vol_p12i <- as.data.frame(c(summary(vol12)$coefficients[14,5]))
vol_b12s <- as.data.frame(c(summary(vol12)$coefficients[16,1]))
vol_se12s <- as.data.frame(c(summary(vol12)$coefficients[16,2]))
vol_p12s <- as.data.frame(c(summary(vol12)$coefficients[16,5]))

vol13 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_mobfrlh_vol_cent + scale(subj_mobfrlh_vol) + site_mobfrlh_vol_cent*wave + scale(subj_mobfrlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b13i <- as.data.frame(c(summary(vol13)$coefficients[14,1]))
vol_se13i <- as.data.frame(c(summary(vol13)$coefficients[14,2]))
vol_p13i <- as.data.frame(c(summary(vol13)$coefficients[14,5]))
vol_b13s <- as.data.frame(c(summary(vol13)$coefficients[16,1]))
vol_se13s <- as.data.frame(c(summary(vol13)$coefficients[16,2]))
vol_p13s <- as.data.frame(c(summary(vol13)$coefficients[16,5]))

vol14 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_mdtmlh_vol_cent + scale(subj_mdtmlh_vol) + site_mdtmlh_vol_cent*wave + scale(subj_mdtmlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b14i <- as.data.frame(c(summary(vol14)$coefficients[14,1]))
vol_se14i <- as.data.frame(c(summary(vol14)$coefficients[14,2]))
vol_p14i <- as.data.frame(c(summary(vol14)$coefficients[14,5]))
vol_b14s <- as.data.frame(c(summary(vol14)$coefficients[16,1]))
vol_se14s <- as.data.frame(c(summary(vol14)$coefficients[16,2]))
vol_p14s <- as.data.frame(c(summary(vol14)$coefficients[16,5]))

vol15 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_parahpallh_vol_cent + scale(subj_parahpallh_vol) + site_parahpallh_vol_cent*wave + scale(subj_parahpallh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b15i <- as.data.frame(c(summary(vol15)$coefficients[14,1]))
vol_se15i <- as.data.frame(c(summary(vol15)$coefficients[14,2]))
vol_p15i <- as.data.frame(c(summary(vol15)$coefficients[14,5]))
vol_b15s <- as.data.frame(c(summary(vol15)$coefficients[16,1]))
vol_se15s <- as.data.frame(c(summary(vol15)$coefficients[16,2]))
vol_p15s <- as.data.frame(c(summary(vol15)$coefficients[16,5]))

vol16 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_paracnlh_vol_cent + scale(subj_paracnlh_vol) + site_paracnlh_vol_cent*wave + scale(subj_paracnlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b16i <- as.data.frame(c(summary(vol16)$coefficients[14,1]))
vol_se16i <- as.data.frame(c(summary(vol16)$coefficients[14,2]))
vol_p16i <- as.data.frame(c(summary(vol16)$coefficients[14,5]))
vol_b16s <- as.data.frame(c(summary(vol16)$coefficients[16,1]))
vol_se16s <- as.data.frame(c(summary(vol16)$coefficients[16,2]))
vol_p16s <- as.data.frame(c(summary(vol16)$coefficients[16,5]))

vol17 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_parsopclh_vol_cent + scale(subj_parsopclh_vol) + site_parsopclh_vol_cent*wave + scale(subj_parsopclh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b17i <- as.data.frame(c(summary(vol17)$coefficients[14,1]))
vol_se17i <- as.data.frame(c(summary(vol17)$coefficients[14,2]))
vol_p17i <- as.data.frame(c(summary(vol17)$coefficients[14,5]))
vol_b17s <- as.data.frame(c(summary(vol17)$coefficients[16,1]))
vol_se17s <- as.data.frame(c(summary(vol17)$coefficients[16,2]))
vol_p17s <- as.data.frame(c(summary(vol17)$coefficients[16,5]))

vol18 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_parsobislh_vol_cent + scale(subj_parsobislh_vol) + site_parsobislh_vol_cent*wave + scale(subj_parsobislh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b18i <- as.data.frame(c(summary(vol18)$coefficients[14,1]))
vol_se18i <- as.data.frame(c(summary(vol18)$coefficients[14,2]))
vol_p18i <- as.data.frame(c(summary(vol18)$coefficients[14,5]))
vol_b18s <- as.data.frame(c(summary(vol18)$coefficients[16,1]))
vol_se18s <- as.data.frame(c(summary(vol18)$coefficients[16,2]))
vol_p18s <- as.data.frame(c(summary(vol18)$coefficients[16,5]))

vol19 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_parstgrislh_vol_cent + scale(subj_parstgrislh_vol) + site_parstgrislh_vol_cent*wave + scale(subj_parstgrislh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b19i <- as.data.frame(c(summary(vol19)$coefficients[14,1]))
vol_se19i <- as.data.frame(c(summary(vol19)$coefficients[14,2]))
vol_p19i <- as.data.frame(c(summary(vol19)$coefficients[14,5]))
vol_b19s <- as.data.frame(c(summary(vol19)$coefficients[16,1]))
vol_se19s <- as.data.frame(c(summary(vol19)$coefficients[16,2]))
vol_p19s <- as.data.frame(c(summary(vol19)$coefficients[16,5]))


vol20 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_pericclh_vol_cent + scale(subj_pericclh_vol) + site_pericclh_vol_cent*wave + scale(subj_pericclh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b20i <- as.data.frame(c(summary(vol20)$coefficients[14,1]))
vol_se20i <- as.data.frame(c(summary(vol20)$coefficients[14,2]))
vol_p20i <- as.data.frame(c(summary(vol20)$coefficients[14,5]))
vol_b20s <- as.data.frame(c(summary(vol20)$coefficients[16,1]))
vol_se20s <- as.data.frame(c(summary(vol20)$coefficients[16,2]))
vol_p20s <- as.data.frame(c(summary(vol20)$coefficients[16,5]))

vol21 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_postcnlh_vol_cent + scale(subj_postcnlh_vol) + site_postcnlh_vol_cent*wave + scale(subj_postcnlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b21i <- as.data.frame(c(summary(vol21)$coefficients[14,1]))
vol_se21i <- as.data.frame(c(summary(vol21)$coefficients[14,2]))
vol_p21i <- as.data.frame(c(summary(vol21)$coefficients[14,5]))
vol_b21s <- as.data.frame(c(summary(vol21)$coefficients[16,1]))
vol_se21s <- as.data.frame(c(summary(vol21)$coefficients[16,2]))
vol_p21s <- as.data.frame(c(summary(vol21)$coefficients[16,5]))

vol22 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_ptcatelh_vol_cent + scale(subj_ptcatelh_vol) + site_ptcatelh_vol_cent*wave + scale(subj_ptcatelh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b22i <- as.data.frame(c(summary(vol22)$coefficients[14,1]))
vol_se22i <- as.data.frame(c(summary(vol22)$coefficients[14,2]))
vol_p22i <- as.data.frame(c(summary(vol22)$coefficients[14,5]))
vol_b22s <- as.data.frame(c(summary(vol22)$coefficients[16,1]))
vol_se22s <- as.data.frame(c(summary(vol22)$coefficients[16,2]))
vol_p22s <- as.data.frame(c(summary(vol22)$coefficients[16,5]))

vol23 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_precnlh_vol_cent + scale(subj_precnlh_vol) + site_precnlh_vol_cent*wave + scale(subj_precnlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b23i <- as.data.frame(c(summary(vol23)$coefficients[14,1]))
vol_se23i <- as.data.frame(c(summary(vol23)$coefficients[14,2]))
vol_p23i <- as.data.frame(c(summary(vol23)$coefficients[14,5]))
vol_b23s <- as.data.frame(c(summary(vol23)$coefficients[16,1]))
vol_se23s <- as.data.frame(c(summary(vol23)$coefficients[16,2]))
vol_p23s <- as.data.frame(c(summary(vol23)$coefficients[16,5]))

vol24 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_pclh_vol_cent + scale(subj_pclh_vol) + site_pclh_vol_cent*wave + scale(subj_pclh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b24i <- as.data.frame(c(summary(vol24)$coefficients[14,1]))
vol_se24i <- as.data.frame(c(summary(vol24)$coefficients[14,2]))
vol_p24i <- as.data.frame(c(summary(vol24)$coefficients[14,5]))
vol_b24s <- as.data.frame(c(summary(vol24)$coefficients[16,1]))
vol_se24s <- as.data.frame(c(summary(vol24)$coefficients[16,2]))
vol_p24s <- as.data.frame(c(summary(vol24)$coefficients[16,5]))

vol25 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_rracatelh_vol_cent + scale(subj_rracatelh_vol) + site_rracatelh_vol_cent*wave + scale(subj_rracatelh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b25i <- as.data.frame(c(summary(vol25)$coefficients[14,1]))
vol_se25i <- as.data.frame(c(summary(vol25)$coefficients[14,2]))
vol_p25i <- as.data.frame(c(summary(vol25)$coefficients[14,5]))
vol_b25s <- as.data.frame(c(summary(vol25)$coefficients[16,1]))
vol_se25s <- as.data.frame(c(summary(vol25)$coefficients[16,2]))
vol_p25s <- as.data.frame(c(summary(vol25)$coefficients[16,5]))

vol26 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_rrmdfrlh_vol_cent + scale(subj_rrmdfrlh_vol) + site_rrmdfrlh_vol_cent*wave + scale(subj_rrmdfrlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b26i <- as.data.frame(c(summary(vol26)$coefficients[14,1]))
vol_se26i <- as.data.frame(c(summary(vol26)$coefficients[14,2]))
vol_p26i <- as.data.frame(c(summary(vol26)$coefficients[14,5]))
vol_b26s <- as.data.frame(c(summary(vol26)$coefficients[16,1]))
vol_se26s <- as.data.frame(c(summary(vol26)$coefficients[16,2]))
vol_p26s <- as.data.frame(c(summary(vol26)$coefficients[16,5]))

vol27 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_sufrlh_vol_cent + scale(subj_sufrlh_vol) + site_sufrlh_vol_cent*wave + scale(subj_sufrlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b27i <- as.data.frame(c(summary(vol27)$coefficients[14,1]))
vol_se27i <- as.data.frame(c(summary(vol27)$coefficients[14,2]))
vol_p27i <- as.data.frame(c(summary(vol27)$coefficients[14,5]))
vol_b27s <- as.data.frame(c(summary(vol27)$coefficients[16,1]))
vol_se27s <- as.data.frame(c(summary(vol27)$coefficients[16,2]))
vol_p27s <- as.data.frame(c(summary(vol27)$coefficients[16,5]))

vol28 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_supllh_vol_cent + scale(subj_supllh_vol) + site_supllh_vol_cent*wave + scale(subj_supllh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b28i <- as.data.frame(c(summary(vol28)$coefficients[14,1]))
vol_se28i <- as.data.frame(c(summary(vol28)$coefficients[14,2]))
vol_p28i <- as.data.frame(c(summary(vol28)$coefficients[14,5]))
vol_b28s <- as.data.frame(c(summary(vol28)$coefficients[16,1]))
vol_se28s <- as.data.frame(c(summary(vol28)$coefficients[16,2]))
vol_p28s <- as.data.frame(c(summary(vol28)$coefficients[16,5]))

vol29 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_sutmlh_vol_cent + scale(subj_sutmlh_vol) + site_sutmlh_vol_cent*wave + scale(subj_sutmlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b29i <- as.data.frame(c(summary(vol29)$coefficients[14,1]))
vol_se29i <- as.data.frame(c(summary(vol29)$coefficients[14,2]))
vol_p29i <- as.data.frame(c(summary(vol29)$coefficients[14,5]))
vol_b29s <- as.data.frame(c(summary(vol29)$coefficients[16,1]))
vol_se29s <- as.data.frame(c(summary(vol29)$coefficients[16,2]))
vol_p29s <- as.data.frame(c(summary(vol29)$coefficients[16,5]))

vol30 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_smlh_vol_cent + scale(subj_smlh_vol) + site_smlh_vol_cent*wave + scale(subj_smlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b30i <- as.data.frame(c(summary(vol30)$coefficients[14,1]))
vol_se30i <- as.data.frame(c(summary(vol30)$coefficients[14,2]))
vol_p30i <- as.data.frame(c(summary(vol30)$coefficients[14,5]))
vol_b30s <- as.data.frame(c(summary(vol30)$coefficients[16,1]))
vol_se30s <- as.data.frame(c(summary(vol30)$coefficients[16,2]))
vol_p30s <- as.data.frame(c(summary(vol30)$coefficients[16,5]))

vol31 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_frpolelh_vol_cent + scale(subj_frpolelh_vol) + site_frpolelh_vol_cent*wave + scale(subj_frpolelh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b31i <- as.data.frame(c(summary(vol31)$coefficients[14,1]))
vol_se31i <- as.data.frame(c(summary(vol31)$coefficients[14,2]))
vol_p31i <- as.data.frame(c(summary(vol31)$coefficients[14,5]))
vol_b31s <- as.data.frame(c(summary(vol31)$coefficients[16,1]))
vol_se31s <- as.data.frame(c(summary(vol31)$coefficients[16,2]))
vol_p31s <- as.data.frame(c(summary(vol31)$coefficients[16,5]))

vol32 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_tmpolelh_vol_cent + scale(subj_tmpolelh_vol) + site_tmpolelh_vol_cent*wave + scale(subj_tmpolelh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b32i <- as.data.frame(c(summary(vol32)$coefficients[14,1]))
vol_se32i <- as.data.frame(c(summary(vol32)$coefficients[14,2]))
vol_p32i <- as.data.frame(c(summary(vol32)$coefficients[14,5]))
vol_b32s <- as.data.frame(c(summary(vol32)$coefficients[16,1]))
vol_se32s <- as.data.frame(c(summary(vol32)$coefficients[16,2]))
vol_p32s <- as.data.frame(c(summary(vol32)$coefficients[16,5]))

vol33 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_trvtmlh_vol_cent + scale(subj_trvtmlh_vol) + site_trvtmlh_vol_cent*wave + scale(subj_trvtmlh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b33i <- as.data.frame(c(summary(vol33)$coefficients[14,1]))
vol_se33i <- as.data.frame(c(summary(vol33)$coefficients[14,2]))
vol_p33i <- as.data.frame(c(summary(vol33)$coefficients[14,5]))
vol_b33s <- as.data.frame(c(summary(vol33)$coefficients[16,1]))
vol_se33s <- as.data.frame(c(summary(vol33)$coefficients[16,2]))
vol_p33s <- as.data.frame(c(summary(vol33)$coefficients[16,5]))

vol34 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_insulalh_vol_cent + scale(subj_insulalh_vol) + site_insulalh_vol_cent*wave + scale(subj_insulalh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b34i <- as.data.frame(c(summary(vol34)$coefficients[14,1]))
vol_se34i <- as.data.frame(c(summary(vol34)$coefficients[14,2]))
vol_p34i <- as.data.frame(c(summary(vol34)$coefficients[14,5]))
vol_b34s <- as.data.frame(c(summary(vol34)$coefficients[16,1]))
vol_se34s <- as.data.frame(c(summary(vol34)$coefficients[16,2]))
vol_p34s <- as.data.frame(c(summary(vol34)$coefficients[16,5]))

vol35 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_banksstsrh_vol_cent + scale(subj_banksstsrh_vol) + site_banksstsrh_vol_cent*wave + scale(subj_banksstsrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b35i <- as.data.frame(c(summary(vol35)$coefficients[14,1]))
vol_se35i <- as.data.frame(c(summary(vol35)$coefficients[14,2]))
vol_p35i <- as.data.frame(c(summary(vol35)$coefficients[14,5]))
vol_b35s <- as.data.frame(c(summary(vol35)$coefficients[16,1]))
vol_se35s <- as.data.frame(c(summary(vol35)$coefficients[16,2]))
vol_p35s <- as.data.frame(c(summary(vol35)$coefficients[16,5]))

vol36 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_cdacaterh_vol_cent + scale(subj_cdacaterh_vol) + site_cdacaterh_vol_cent*wave + scale(subj_cdacaterh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b36i <- as.data.frame(c(summary(vol36)$coefficients[14,1]))
vol_se36i <- as.data.frame(c(summary(vol36)$coefficients[14,2]))
vol_p36i <- as.data.frame(c(summary(vol36)$coefficients[14,5]))
vol_b36s <- as.data.frame(c(summary(vol36)$coefficients[16,1]))
vol_se36s <- as.data.frame(c(summary(vol36)$coefficients[16,2]))
vol_p36s <- as.data.frame(c(summary(vol36)$coefficients[16,5]))

vol37 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_cdmdfrrh_vol_cent + scale(subj_cdmdfrrh_vol) + site_cdmdfrrh_vol_cent*wave + scale(subj_cdmdfrrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b37i <- as.data.frame(c(summary(vol37)$coefficients[14,1]))
vol_se37i <- as.data.frame(c(summary(vol37)$coefficients[14,2]))
vol_p37i <- as.data.frame(c(summary(vol37)$coefficients[14,5]))
vol_b37s <- as.data.frame(c(summary(vol37)$coefficients[16,1]))
vol_se37s <- as.data.frame(c(summary(vol37)$coefficients[16,2]))
vol_p37s <- as.data.frame(c(summary(vol37)$coefficients[16,5]))

vol38 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_cuneusrh_vol_cent + scale(subj_cuneusrh_vol) + site_cuneusrh_vol_cent*wave + scale(subj_cuneusrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b38i <- as.data.frame(c(summary(vol38)$coefficients[14,1]))
vol_se38i <- as.data.frame(c(summary(vol38)$coefficients[14,2]))
vol_p38i <- as.data.frame(c(summary(vol38)$coefficients[14,5]))
vol_b38s <- as.data.frame(c(summary(vol38)$coefficients[16,1]))
vol_se38s <- as.data.frame(c(summary(vol38)$coefficients[16,2]))
vol_p38s <- as.data.frame(c(summary(vol38)$coefficients[16,5]))

vol39 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_ehinalrh_vol_cent + scale(subj_ehinalrh_vol) + site_ehinalrh_vol_cent*wave + scale(subj_ehinalrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b39i <- as.data.frame(c(summary(vol39)$coefficients[14,1]))
vol_se39i <- as.data.frame(c(summary(vol39)$coefficients[14,2]))
vol_p39i <- as.data.frame(c(summary(vol39)$coefficients[14,5]))
vol_b39s <- as.data.frame(c(summary(vol39)$coefficients[16,1]))
vol_se39s <- as.data.frame(c(summary(vol39)$coefficients[16,2]))
vol_p39s <- as.data.frame(c(summary(vol39)$coefficients[16,5]))

vol40 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_fusiformrh_vol_cent + scale(subj_fusiformrh_vol) + site_fusiformrh_vol_cent*wave + scale(subj_fusiformrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b40i <- as.data.frame(c(summary(vol40)$coefficients[14,1]))
vol_se40i <- as.data.frame(c(summary(vol40)$coefficients[14,2]))
vol_p40i <- as.data.frame(c(summary(vol40)$coefficients[14,5]))
vol_b40s <- as.data.frame(c(summary(vol40)$coefficients[16,1]))
vol_se40s <- as.data.frame(c(summary(vol40)$coefficients[16,2]))
vol_p40s <- as.data.frame(c(summary(vol40)$coefficients[16,5]))

vol41 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_ifplrh_vol_cent + scale(subj_ifplrh_vol) + site_ifplrh_vol_cent*wave + scale(subj_ifplrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b41i <- as.data.frame(c(summary(vol41)$coefficients[14,1]))
vol_se41i <- as.data.frame(c(summary(vol41)$coefficients[14,2]))
vol_p41i <- as.data.frame(c(summary(vol41)$coefficients[14,5]))
vol_b41s <- as.data.frame(c(summary(vol41)$coefficients[16,1]))
vol_se41s <- as.data.frame(c(summary(vol41)$coefficients[16,2]))
vol_p41s <- as.data.frame(c(summary(vol41)$coefficients[16,5]))

vol42 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_iftmrh_vol_cent + scale(subj_iftmrh_vol) + site_iftmrh_vol_cent*wave + scale(subj_iftmrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b42i <- as.data.frame(c(summary(vol42)$coefficients[14,1]))
vol_se42i <- as.data.frame(c(summary(vol42)$coefficients[14,2]))
vol_p42i <- as.data.frame(c(summary(vol42)$coefficients[14,5]))
vol_b42s <- as.data.frame(c(summary(vol42)$coefficients[16,1]))
vol_se42s <- as.data.frame(c(summary(vol42)$coefficients[16,2]))
vol_p42s <- as.data.frame(c(summary(vol42)$coefficients[16,5]))

vol43 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_ihcaterh_vol_cent + scale(subj_ihcaterh_vol) + site_ihcaterh_vol_cent*wave + scale(subj_ihcaterh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b43i <- as.data.frame(c(summary(vol43)$coefficients[14,1]))
vol_se43i <- as.data.frame(c(summary(vol43)$coefficients[14,2]))
vol_p43i <- as.data.frame(c(summary(vol43)$coefficients[14,5]))
vol_b43s <- as.data.frame(c(summary(vol43)$coefficients[16,1]))
vol_se43s <- as.data.frame(c(summary(vol43)$coefficients[16,2]))
vol_p43s <- as.data.frame(c(summary(vol43)$coefficients[16,5]))

vol44 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_loccrh_vol_cent + scale(subj_loccrh_vol) + site_loccrh_vol_cent*wave + scale(subj_loccrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b44i <- as.data.frame(c(summary(vol44)$coefficients[14,1]))
vol_se44i <- as.data.frame(c(summary(vol44)$coefficients[14,2]))
vol_p44i <- as.data.frame(c(summary(vol44)$coefficients[14,5]))
vol_b44s <- as.data.frame(c(summary(vol44)$coefficients[16,1]))
vol_se44s <- as.data.frame(c(summary(vol44)$coefficients[16,2]))
vol_p44s <- as.data.frame(c(summary(vol44)$coefficients[16,5]))

vol45 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_lobfrrh_vol_cent + scale(subj_lobfrrh_vol) + site_lobfrrh_vol_cent*wave + scale(subj_lobfrrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b45i <- as.data.frame(c(summary(vol45)$coefficients[14,1]))
vol_se45i <- as.data.frame(c(summary(vol45)$coefficients[14,2]))
vol_p45i <- as.data.frame(c(summary(vol45)$coefficients[14,5]))
vol_b45s <- as.data.frame(c(summary(vol45)$coefficients[16,1]))
vol_se45s <- as.data.frame(c(summary(vol45)$coefficients[16,2]))
vol_p45s <- as.data.frame(c(summary(vol45)$coefficients[16,5]))

vol46 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_lingualrh_vol_cent + scale(subj_lingualrh_vol) + site_lingualrh_vol_cent*wave + scale(subj_lingualrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b46i <- as.data.frame(c(summary(vol46)$coefficients[14,1]))
vol_se46i <- as.data.frame(c(summary(vol46)$coefficients[14,2]))
vol_p46i <- as.data.frame(c(summary(vol46)$coefficients[14,5]))
vol_b46s <- as.data.frame(c(summary(vol46)$coefficients[16,1]))
vol_se46s <- as.data.frame(c(summary(vol46)$coefficients[16,2]))
vol_p46s <- as.data.frame(c(summary(vol46)$coefficients[16,5]))

vol47 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_mobfrrh_vol_cent + scale(subj_mobfrrh_vol) + site_mobfrrh_vol_cent*wave + scale(subj_mobfrrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b47i <- as.data.frame(c(summary(vol47)$coefficients[14,1]))
vol_se47i <- as.data.frame(c(summary(vol47)$coefficients[14,2]))
vol_p47i <- as.data.frame(c(summary(vol47)$coefficients[14,5]))
vol_b47s <- as.data.frame(c(summary(vol47)$coefficients[16,1]))
vol_se47s <- as.data.frame(c(summary(vol47)$coefficients[16,2]))
vol_p47s <- as.data.frame(c(summary(vol47)$coefficients[16,5]))

vol48 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_mdtmrh_vol_cent + scale(subj_mdtmrh_vol) + site_mdtmrh_vol_cent*wave + scale(subj_mdtmrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b48i <- as.data.frame(c(summary(vol48)$coefficients[14,1]))
vol_se48i <- as.data.frame(c(summary(vol48)$coefficients[14,2]))
vol_p48i <- as.data.frame(c(summary(vol48)$coefficients[14,5]))
vol_b48s <- as.data.frame(c(summary(vol48)$coefficients[16,1]))
vol_se48s <- as.data.frame(c(summary(vol48)$coefficients[16,2]))
vol_p48s <- as.data.frame(c(summary(vol48)$coefficients[16,5]))

vol49 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_parahpalrh_vol_cent + scale(subj_parahpalrh_vol) + site_parahpalrh_vol_cent*wave + scale(subj_parahpalrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b49i <- as.data.frame(c(summary(vol49)$coefficients[14,1]))
vol_se49i <- as.data.frame(c(summary(vol49)$coefficients[14,2]))
vol_p49i <- as.data.frame(c(summary(vol49)$coefficients[14,5]))
vol_b49s <- as.data.frame(c(summary(vol49)$coefficients[16,1]))
vol_se49s <- as.data.frame(c(summary(vol49)$coefficients[16,2]))
vol_p49s <- as.data.frame(c(summary(vol49)$coefficients[16,5]))

vol50 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_paracnrh_vol_cent + scale(subj_paracnrh_vol) + site_paracnrh_vol_cent*wave + scale(subj_paracnrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b50i <- as.data.frame(c(summary(vol50)$coefficients[14,1]))
vol_se50i <- as.data.frame(c(summary(vol50)$coefficients[14,2]))
vol_p50i <- as.data.frame(c(summary(vol50)$coefficients[14,5]))
vol_b50s <- as.data.frame(c(summary(vol50)$coefficients[16,1]))
vol_se50s <- as.data.frame(c(summary(vol50)$coefficients[16,2]))
vol_p50s <- as.data.frame(c(summary(vol50)$coefficients[16,5]))

vol51 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_parsopcrh_vol_cent + scale(subj_parsopcrh_vol) + site_parsopcrh_vol_cent*wave + scale(subj_parsopcrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b51i <- as.data.frame(c(summary(vol51)$coefficients[14,1]))
vol_se51i <- as.data.frame(c(summary(vol51)$coefficients[14,2]))
vol_p51i <- as.data.frame(c(summary(vol51)$coefficients[14,5]))
vol_b51s <- as.data.frame(c(summary(vol51)$coefficients[16,1]))
vol_se51s <- as.data.frame(c(summary(vol51)$coefficients[16,2]))
vol_p51s <- as.data.frame(c(summary(vol51)$coefficients[16,5]))

vol52 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_parsobisrh_vol_cent + scale(subj_parsobisrh_vol) + site_parsobisrh_vol_cent*wave + scale(subj_parsobisrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b52i <- as.data.frame(c(summary(vol52)$coefficients[14,1]))
vol_se52i <- as.data.frame(c(summary(vol52)$coefficients[14,2]))
vol_p52i <- as.data.frame(c(summary(vol52)$coefficients[14,5]))
vol_b52s <- as.data.frame(c(summary(vol52)$coefficients[16,1]))
vol_se52s <- as.data.frame(c(summary(vol52)$coefficients[16,2]))
vol_p52s <- as.data.frame(c(summary(vol52)$coefficients[16,5]))

vol53 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_parstgrisrh_vol_cent + scale(subj_parstgrisrh_vol) + site_parstgrisrh_vol_cent*wave + scale(subj_parstgrisrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b53i <- as.data.frame(c(summary(vol53)$coefficients[14,1]))
vol_se53i <- as.data.frame(c(summary(vol53)$coefficients[14,2]))
vol_p53i <- as.data.frame(c(summary(vol53)$coefficients[14,5]))
vol_b53s <- as.data.frame(c(summary(vol53)$coefficients[16,1]))
vol_se53s <- as.data.frame(c(summary(vol53)$coefficients[16,2]))
vol_p53s <- as.data.frame(c(summary(vol53)$coefficients[16,5]))

vol54 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_periccrh_vol_cent + scale(subj_periccrh_vol) + site_periccrh_vol_cent*wave + scale(subj_periccrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b54i <- as.data.frame(c(summary(vol54)$coefficients[14,1]))
vol_se54i <- as.data.frame(c(summary(vol54)$coefficients[14,2]))
vol_p54i <- as.data.frame(c(summary(vol54)$coefficients[14,5]))
vol_b54s <- as.data.frame(c(summary(vol54)$coefficients[16,1]))
vol_se54s <- as.data.frame(c(summary(vol54)$coefficients[16,2]))
vol_p54s <- as.data.frame(c(summary(vol54)$coefficients[16,5]))

vol55 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_postcnrh_vol_cent + scale(subj_postcnrh_vol) + site_postcnrh_vol_cent*wave + scale(subj_postcnrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b55i <- as.data.frame(c(summary(vol55)$coefficients[14,1]))
vol_se55i <- as.data.frame(c(summary(vol55)$coefficients[14,2]))
vol_p55i <- as.data.frame(c(summary(vol55)$coefficients[14,5]))
vol_b55s <- as.data.frame(c(summary(vol55)$coefficients[16,1]))
vol_se55s <- as.data.frame(c(summary(vol55)$coefficients[16,2]))
vol_p55s <- as.data.frame(c(summary(vol55)$coefficients[16,5]))

vol56 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_ptcaterh_vol_cent + scale(subj_ptcaterh_vol) + site_ptcaterh_vol_cent*wave + scale(subj_ptcaterh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b56i <- as.data.frame(c(summary(vol56)$coefficients[14,1]))
vol_se56i <- as.data.frame(c(summary(vol56)$coefficients[14,2]))
vol_p56i <- as.data.frame(c(summary(vol56)$coefficients[14,5]))
vol_b56s <- as.data.frame(c(summary(vol56)$coefficients[16,1]))
vol_se56s <- as.data.frame(c(summary(vol56)$coefficients[16,2]))
vol_p56s <- as.data.frame(c(summary(vol56)$coefficients[16,5]))

vol57 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_precnrh_vol_cent + scale(subj_precnrh_vol) + site_precnrh_vol_cent*wave + scale(subj_precnrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b57i <- as.data.frame(c(summary(vol57)$coefficients[14,1]))
vol_se57i <- as.data.frame(c(summary(vol57)$coefficients[14,2]))
vol_p57i <- as.data.frame(c(summary(vol57)$coefficients[14,5]))
vol_b57s <- as.data.frame(c(summary(vol57)$coefficients[16,1]))
vol_se57s <- as.data.frame(c(summary(vol57)$coefficients[16,2]))
vol_p57s <- as.data.frame(c(summary(vol57)$coefficients[16,5]))

vol58 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_pcrh_vol_cent + scale(subj_pcrh_vol) + site_pcrh_vol_cent*wave + scale(subj_pcrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b58i <- as.data.frame(c(summary(vol58)$coefficients[14,1]))
vol_se58i <- as.data.frame(c(summary(vol58)$coefficients[14,2]))
vol_p58i <- as.data.frame(c(summary(vol58)$coefficients[14,5]))
vol_b58s <- as.data.frame(c(summary(vol58)$coefficients[16,1]))
vol_se58s <- as.data.frame(c(summary(vol58)$coefficients[16,2]))
vol_p58s <- as.data.frame(c(summary(vol58)$coefficients[16,5]))

vol59 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_rracaterh_vol_cent + scale(subj_rracaterh_vol) + site_rracaterh_vol_cent*wave + scale(subj_rracaterh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b59i <- as.data.frame(c(summary(vol59)$coefficients[14,1]))
vol_se59i <- as.data.frame(c(summary(vol59)$coefficients[14,2]))
vol_p59i <- as.data.frame(c(summary(vol59)$coefficients[14,5]))
vol_b59s <- as.data.frame(c(summary(vol59)$coefficients[16,1]))
vol_se59s <- as.data.frame(c(summary(vol59)$coefficients[16,2]))
vol_p59s <- as.data.frame(c(summary(vol59)$coefficients[16,5]))

vol60 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_rrmdfrrh_vol_cent + scale(subj_rrmdfrrh_vol) + site_rrmdfrrh_vol_cent*wave + scale(subj_rrmdfrrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b60i <- as.data.frame(c(summary(vol60)$coefficients[14,1]))
vol_se60i <- as.data.frame(c(summary(vol60)$coefficients[14,2]))
vol_p60i <- as.data.frame(c(summary(vol60)$coefficients[14,5]))
vol_b60s <- as.data.frame(c(summary(vol60)$coefficients[16,1]))
vol_se60s <- as.data.frame(c(summary(vol60)$coefficients[16,2]))
vol_p60s <- as.data.frame(c(summary(vol60)$coefficients[16,5]))

vol61 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_sufrrh_vol_cent + scale(subj_sufrrh_vol) + site_sufrrh_vol_cent*wave + scale(subj_sufrrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b61i <- as.data.frame(c(summary(vol61)$coefficients[14,1]))
vol_se61i <- as.data.frame(c(summary(vol61)$coefficients[14,2]))
vol_p61i <- as.data.frame(c(summary(vol61)$coefficients[14,5]))
vol_b61s <- as.data.frame(c(summary(vol61)$coefficients[16,1]))
vol_se61s <- as.data.frame(c(summary(vol61)$coefficients[16,2]))
vol_p61s <- as.data.frame(c(summary(vol61)$coefficients[16,5]))

vol62 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_suplrh_vol_cent + scale(subj_suplrh_vol) + site_suplrh_vol_cent*wave + scale(subj_suplrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b62i <- as.data.frame(c(summary(vol62)$coefficients[14,1]))
vol_se62i <- as.data.frame(c(summary(vol62)$coefficients[14,2]))
vol_p62i <- as.data.frame(c(summary(vol62)$coefficients[14,5]))
vol_b62s <- as.data.frame(c(summary(vol62)$coefficients[16,1]))
vol_se62s <- as.data.frame(c(summary(vol62)$coefficients[16,2]))
vol_p62s <- as.data.frame(c(summary(vol62)$coefficients[16,5]))

vol63 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_sutmrh_vol_cent + scale(subj_sutmrh_vol) + site_sutmrh_vol_cent*wave + scale(subj_sutmrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b63i <- as.data.frame(c(summary(vol63)$coefficients[14,1]))
vol_se63i <- as.data.frame(c(summary(vol63)$coefficients[14,2]))
vol_p63i <- as.data.frame(c(summary(vol63)$coefficients[14,5]))
vol_b63s <- as.data.frame(c(summary(vol63)$coefficients[16,1]))
vol_se63s <- as.data.frame(c(summary(vol63)$coefficients[16,2]))
vol_p63s <- as.data.frame(c(summary(vol63)$coefficients[16,5]))

vol64 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_smrh_vol_cent + scale(subj_smrh_vol) + site_smrh_vol_cent*wave + scale(subj_smrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b64i <- as.data.frame(c(summary(vol64)$coefficients[14,1]))
vol_se64i <- as.data.frame(c(summary(vol64)$coefficients[14,2]))
vol_p64i <- as.data.frame(c(summary(vol64)$coefficients[14,5]))
vol_b64s <- as.data.frame(c(summary(vol64)$coefficients[16,1]))
vol_se64s <- as.data.frame(c(summary(vol64)$coefficients[16,2]))
vol_p64s <- as.data.frame(c(summary(vol64)$coefficients[16,5]))

vol65 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_frpolerh_vol_cent + scale(subj_frpolerh_vol) + site_frpolerh_vol_cent*wave + scale(subj_frpolerh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b65i <- as.data.frame(c(summary(vol65)$coefficients[14,1]))
vol_se65i <- as.data.frame(c(summary(vol65)$coefficients[14,2]))
vol_p65i <- as.data.frame(c(summary(vol65)$coefficients[14,5]))
vol_b65s <- as.data.frame(c(summary(vol65)$coefficients[16,1]))
vol_se65s <- as.data.frame(c(summary(vol65)$coefficients[16,2]))
vol_p65s <- as.data.frame(c(summary(vol65)$coefficients[16,5]))

vol66 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_tmpolerh_vol_cent + scale(subj_tmpolerh_vol) + site_tmpolerh_vol_cent*wave + scale(subj_tmpolerh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b66i <- as.data.frame(c(summary(vol66)$coefficients[14,1]))
vol_se66i <- as.data.frame(c(summary(vol66)$coefficients[14,2]))
vol_p66i <- as.data.frame(c(summary(vol66)$coefficients[14,5]))
vol_b66s <- as.data.frame(c(summary(vol66)$coefficients[16,1]))
vol_se66s <- as.data.frame(c(summary(vol66)$coefficients[16,2]))
vol_p66s <- as.data.frame(c(summary(vol66)$coefficients[16,5]))

vol67 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_trvtmrh_vol_cent + scale(subj_trvtmrh_vol) + site_trvtmrh_vol_cent*wave + scale(subj_trvtmrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b67i <- as.data.frame(c(summary(vol67)$coefficients[14,1]))
vol_se67i <- as.data.frame(c(summary(vol67)$coefficients[14,2]))
vol_p67i <- as.data.frame(c(summary(vol67)$coefficients[14,5]))
vol_b67s <- as.data.frame(c(summary(vol67)$coefficients[16,1]))
vol_se67s <- as.data.frame(c(summary(vol67)$coefficients[16,2]))
vol_p67s <- as.data.frame(c(summary(vol67)$coefficients[16,5]))

vol68 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_insularh_vol_cent + scale(subj_insularh_vol) + site_insularh_vol_cent*wave + scale(subj_insularh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
vol_b68i <- as.data.frame(c(summary(vol68)$coefficients[14,1]))
vol_se68i <- as.data.frame(c(summary(vol68)$coefficients[14,2]))
vol_p68i <- as.data.frame(c(summary(vol68)$coefficients[14,5]))
vol_b68s <- as.data.frame(c(summary(vol68)$coefficients[16,1]))
vol_se68s <- as.data.frame(c(summary(vol68)$coefficients[16,2]))
vol_p68s <- as.data.frame(c(summary(vol68)$coefficients[16,5]))

#create data frame with all parcel-wise cortical volume st. estimates, SEs, and p-values
vol_parcel <- data.frame(x=c("vol"))

newvolbi <- c(vol_b1i,vol_b2i,vol_b3i,vol_b4i,vol_b5i,vol_b6i,vol_b7i,vol_b8i,vol_b9i,vol_b10i,vol_b11i,vol_b12i,vol_b13i,vol_b14i,vol_b15i,
              vol_b16i,vol_b17i,vol_b18i,vol_b19i,vol_b20i,vol_b21i,vol_b22i,vol_b23i,vol_b24i,vol_b25i,vol_b26i,vol_b27i,vol_b28i,vol_b29i,
              vol_b30i,vol_b31i,vol_b32i,vol_b33i,vol_b34i,vol_b35i,vol_b36i,vol_b37i,vol_b38i,vol_b39i,vol_b40i,vol_b41i,vol_b42i,vol_b43i,
              vol_b44i,vol_b45i,vol_b46i,vol_b47i,vol_b48i,vol_b49i,vol_b50i,vol_b51i,vol_b52i,vol_b53i,vol_b54i,vol_b55i,vol_b56i,vol_b57i,
              vol_b58i,vol_b59i,vol_b60i,vol_b61i,vol_b62i,vol_b63i,vol_b64i,vol_b65i,vol_b66i,vol_b67i,vol_b68i)
vol_parcel <- cbind(vol_parcel,newvolbi)

newvolbs <- c(vol_b1s,vol_b2s,vol_b3s,vol_b4s,vol_b5s,vol_b6s,vol_b7s,vol_b8s,vol_b9s,vol_b10s,vol_b11s,vol_b12s,vol_b13s,vol_b14s,vol_b15s,
              vol_b16s,vol_b17s,vol_b18s,vol_b19s,vol_b20s,vol_b21s,vol_b22s,vol_b23s,vol_b24s,vol_b25s,vol_b26s,vol_b27s,vol_b28s,vol_b29s,
              vol_b30s,vol_b31s,vol_b32s,vol_b33s,vol_b34s,vol_b35s,vol_b36s,vol_b37s,vol_b38s,vol_b39s,vol_b40s,vol_b41s,vol_b42s,vol_b43s,
              vol_b44s,vol_b45s,vol_b46s,vol_b47s,vol_b48s,vol_b49s,vol_b50s,vol_b51s,vol_b52s,vol_b53s,vol_b54s,vol_b55s,vol_b56s,vol_b57s,
              vol_b58s,vol_b59s,vol_b60s,vol_b61s,vol_b62s,vol_b63s,vol_b64s,vol_b65s,vol_b66s,vol_b67s,vol_b68s)
vol_parcel <- cbind(vol_parcel,newvolbs)

newvolsei <- c(vol_se1i,vol_se2i,vol_se3i,vol_se4i,vol_se5i,vol_se6i,vol_se7i,vol_se8i,vol_se9i,vol_se10i,vol_se11i,vol_se12i,vol_se13i,
               vol_se14i,vol_se15i,vol_se16i,vol_se17i,vol_se18i,vol_se19i,vol_se20i,vol_se21i,vol_se22i,vol_se23i,vol_se24i,vol_se25i,
               vol_se26i,vol_se27i,vol_se28i,vol_se29i,vol_se30i,vol_se31i,vol_se32i,vol_se33i,vol_se34i,vol_se35i,vol_se36i,vol_se37i,
               vol_se38i,vol_se39i,vol_se40i,vol_se41i,vol_se42i,vol_se43i,vol_se44i,vol_se45i,vol_se46i,vol_se47i,vol_se48i,vol_se49i,
               vol_se50i,vol_se51i,vol_se52i,vol_se53i,vol_se54i,vol_se55i,vol_se56i,vol_se57i,vol_se58i,vol_se59i,vol_se60i,vol_se61i,
               vol_se62i,vol_se63i,vol_se64i,vol_se65i,vol_se66i,vol_se67i,vol_se68i)
vol_parcel <- cbind(vol_parcel,newvolsei)

newvolses <- c(vol_se1s,vol_se2s,vol_se3s,vol_se4s,vol_se5s,vol_se6s,vol_se7s,vol_se8s,vol_se9s,vol_se10s,vol_se11s,vol_se12s,vol_se13s,
               vol_se14s,vol_se15s,vol_se16s,vol_se17s,vol_se18s,vol_se19s,vol_se20s,vol_se21s,vol_se22s,vol_se23s,vol_se24s,vol_se25s,
               vol_se26s,vol_se27s,vol_se28s,vol_se29s,vol_se30s,vol_se31s,vol_se32s,vol_se33s,vol_se34s,vol_se35s,vol_se36s,vol_se37s,
               vol_se38s,vol_se39s,vol_se40s,vol_se41s,vol_se42s,vol_se43s,vol_se44s,vol_se45s,vol_se46s,vol_se47s,vol_se48s,vol_se49s,
               vol_se50s,vol_se51s,vol_se52s,vol_se53s,vol_se54s,vol_se55s,vol_se56s,vol_se57s,vol_se58s,vol_se59s,vol_se60s,vol_se61s,
               vol_se62s,vol_se63s,vol_se64s,vol_se65s,vol_se66s,vol_se67s,vol_se68s)
vol_parcel <- cbind(vol_parcel,newvolses)

newvolpi <- c(vol_p1i,vol_p2i,vol_p3i,vol_p4i,vol_p5i,vol_p6i,vol_p7i,vol_p8i,vol_p9i,vol_p10i,vol_p11i,vol_p12i,vol_p13i,vol_p14i,vol_p15i,
              vol_p16i,vol_p17i,vol_p18i,vol_p19i,vol_p20i,vol_p21i,vol_p22i,vol_p23i,vol_p24i,vol_p25i,vol_p26i,vol_p27i,vol_p28i,vol_p29i,
              vol_p30i,vol_p31i,vol_p32i,vol_p33i,vol_p34i,vol_p35i,vol_p36i,vol_p37i,vol_p38i,vol_p39i,vol_p40i,vol_p41i,vol_p42i,vol_p43i,
              vol_p44i,vol_p45i,vol_p46i,vol_p47i,vol_p48i,vol_p49i,vol_p50i,vol_p51i,vol_p52i,vol_p53i,vol_p54i,vol_p55i,vol_p56i,vol_p57i,
              vol_p58i,vol_p59i,vol_p60i,vol_p61i,vol_p62i,vol_p63i,vol_p64i,vol_p65i,vol_p66i,vol_p67i,vol_p68i)
vol_parcel <- cbind(vol_parcel,newvolpi)

newvolps <- c(vol_p1s,vol_p2s,vol_p3s,vol_p4s,vol_p5s,vol_p6s,vol_p7s,vol_p8s,vol_p9s,vol_p10s,vol_p11s,vol_p12s,vol_p13s,vol_p14s,vol_p15s,
              vol_p16s,vol_p17s,vol_p18s,vol_p19s,vol_p20s,vol_p21s,vol_p22s,vol_p23s,vol_p24s,vol_p25s,vol_p26s,vol_p27s,vol_p28s,vol_p29s,
              vol_p30s,vol_p31s,vol_p32s,vol_p33s,vol_p34s,vol_p35s,vol_p36s,vol_p37s,vol_p38s,vol_p39s,vol_p40s,vol_p41s,vol_p42s,vol_p43s,
              vol_p44s,vol_p45s,vol_p46s,vol_p47s,vol_p48s,vol_p49s,vol_p50s,vol_p51s,vol_p52s,vol_p53s,vol_p54s,vol_p55s,vol_p56s,vol_p57s,
              vol_p58s,vol_p59s,vol_p60s,vol_p61s,vol_p62s,vol_p63s,vol_p64s,vol_p65s,vol_p66s,vol_p67s,vol_p68s)
vol_parcel <- cbind(vol_parcel,newvolps)

names(vol_parcel) <- c('vol','banksstslh_vol_bi','cdacatelh_vol_bi','cdmdfrlh_vol_bi','cuneuslh_vol_bi','ehinallh_vol_bi','fusiformlh_vol_bi',
                       'ifpllh_vol_bi','iftmlh_vol_bi','ihcatelh_vol_bi','locclh_vol_bi','lobfrlh_vol_bi','linguallh_vol_bi',
                       'mobfrlh_vol_bi','mdtmlh_vol_bi','parahpallh_vol_bi','paracnlh_vol_bi','parsopclh_vol_bi','parsobislh_vol_bi',
                       'parstgrislh_vol_bi','pericclh_vol_bi','postcnlh_vol_bi','ptcatelh_vol_bi','precnlh_vol_bi','pclh_vol_bi',
                       'rracatelh_vol_bi','rrmdfrlh_vol_bi','sufrlh_vol_bi','supllh_vol_bi','sutmlh_vol_bi','smlh_vol_bi','frpolelh_vol_bi',
                       'tmpolelh_vol_bi','trvtmlh_vol_bi','insulalh_vol_bi','banksstsrh_vol_bi','cdacaterh_vol_bi','cdmdfrrh_vol_bi',
                       'cuneusrh_vol_bi','ehinalrh_vol_bi','fusiformrh_vol_bi','ifplrh_vol_bi','iftmrh_vol_bi','ihcaterh_vol_bi',
                       'loccrh_vol_bi','lobfrrh_vol_bi','lingualrh_vol_bi','mobfrrh_vol_bi','mdtmrh_vol_bi','parahpalrh_vol_bi',
                       'paracnrh_vol_bi','parsopcrh_vol_bi','parsobisrh_vol_bi','parstgrisrh_vol_bi','periccrh_vol_bi',
                       'postcnrh_vol_bi','ptcaterh_vol_bi','precnrh_vol_bi','pcrh_vol_bi','rracaterh_vol_bi','rrmdfrrh_vol_bi',
                       'sufrrh_vol_bi','suplrh_vol_bi','sutmrh_vol_bi','smrh_vol_bi','frpolerh_vol_bi',
                       'tmpolerh_vol_bi','trvtmrh_vol_bi','insularh_vol_bi',
                       'banksstslh_vol_bs','cdacatelh_vol_bs','cdmdfrlh_vol_bs','cuneuslh_vol_bs','ehinallh_vol_bs','fusiformlh_vol_bs',
                       'ifpllh_vol_bs','iftmlh_vol_bs','ihcatelh_vol_bs','locclh_vol_bs','lobfrlh_vol_bs','linguallh_vol_bs',
                       'mobfrlh_vol_bs','mdtmlh_vol_bs','parahpallh_vol_bs','paracnlh_vol_bs','parsopclh_vol_bs','parsobislh_vol_bs',
                       'parstgrislh_vol_bs','pericclh_vol_bs','postcnlh_vol_bs','ptcatelh_vol_bs','precnlh_vol_bs','pclh_vol_bs',
                       'rracatelh_vol_bs','rrmdfrlh_vol_bs','sufrlh_vol_bs','supllh_vol_bs','sutmlh_vol_bs','smlh_vol_bs','frpolelh_vol_bs',
                       'tmpolelh_vol_bs','trvtmlh_vol_bs','insulalh_vol_bs','banksstsrh_vol_bs','cdacaterh_vol_bs','cdmdfrrh_vol_bs',
                       'cuneusrh_vol_bs','ehinalrh_vol_bs','fusiformrh_vol_bs','ifplrh_vol_bs','iftmrh_vol_bs','ihcaterh_vol_bs',
                       'loccrh_vol_bs','lobfrrh_vol_bs','lingualrh_vol_bs','mobfrrh_vol_bs','mdtmrh_vol_bs','parahpalrh_vol_bs',
                       'paracnrh_vol_bs','parsopcrh_vol_bs','parsobisrh_vol_bs','parstgrisrh_vol_bs','periccrh_vol_bs',
                       'postcnrh_vol_bs','ptcaterh_vol_bs','precnrh_vol_bs','pcrh_vol_bs','rracaterh_vol_bs','rrmdfrrh_vol_bs',
                       'sufrrh_vol_bs','suplrh_vol_bs','sutmrh_vol_bs','smrh_vol_bs','frpolerh_vol_bs',
                       'tmpolerh_vol_bs','trvtmrh_vol_bs','insularh_vol_bs',
                       'banksstslh_vol_sei','cdacatelh_vol_sei','cdmdfrlh_vol_sei','cuneuslh_vol_sei','ehinallh_vol_sei','fusiformlh_vol_sei',
                       'ifpllh_vol_sei','iftmlh_vol_sei','ihcatelh_vol_sei','locclh_vol_sei','lobfrlh_vol_sei','linguallh_vol_sei',
                       'mobfrlh_vol_sei','mdtmlh_vol_sei','parahpallh_vol_sei','paracnlh_vol_sei','parsopclh_vol_sei',
                       'parsobislh_vol_sei','parstgrislh_vol_sei','pericclh_vol_sei','postcnlh_vol_sei','ptcatelh_vol_sei',
                       'precnlh_vol_sei','pclh_vol_sei','rracatelh_vol_sei','rrmdfrlh_vol_sei','sufrlh_vol_sei','supllh_vol_sei',
                       'sutmlh_vol_sei','smlh_vol_sei','frpolelh_vol_sei','tmpolelh_vol_sei','trvtmlh_vol_sei','insulalh_vol_sei',
                       'banksstsrh_vol_sei','cdacaterh_vol_sei','cdmdfrrh_vol_sei','cuneusrh_vol_sei','ehinalrh_vol_sei',
                       'fusiformrh_vol_sei','ifplrh_vol_sei','iftmrh_vol_sei','ihcaterh_vol_sei','loccrh_vol_sei','lobfrrh_vol_sei',
                       'lingualrh_vol_sei','mobfrrh_vol_sei','mdtmrh_vol_sei','parahpalrh_vol_sei','paracnrh_vol_sei',
                       'parsopcrh_vol_sei','parsobisrh_vol_sei','parstgrisrh_vol_sei','periccrh_vol_sei','postcnrh_vol_sei',
                       'ptcaterh_vol_sei','precnrh_vol_sei','pcrh_vol_sei','rracaterh_vol_sei','rrmdfrrh_vol_sei',
                       'sufrrh_vol_sei','suplrh_vol_sei','sutmrh_vol_sei','smrh_vol_sei','frpolerh_vol_sei',
                       'tmpolerh_vol_sei','trvtmrh_vol_sei','insularh_vol_sei',
                       'banksstslh_vol_ses','cdacatelh_vol_ses','cdmdfrlh_vol_ses','cuneuslh_vol_ses','ehinallh_vol_ses','fusiformlh_vol_ses',
                       'ifpllh_vol_ses','iftmlh_vol_ses','ihcatelh_vol_ses','locclh_vol_ses','lobfrlh_vol_ses','linguallh_vol_ses',
                       'mobfrlh_vol_ses','mdtmlh_vol_ses','parahpallh_vol_ses','paracnlh_vol_ses','parsopclh_vol_ses',
                       'parsobislh_vol_ses','parstgrislh_vol_ses','pericclh_vol_ses','postcnlh_vol_ses','ptcatelh_vol_ses',
                       'precnlh_vol_ses','pclh_vol_ses','rracatelh_vol_ses','rrmdfrlh_vol_ses','sufrlh_vol_ses','supllh_vol_ses',
                       'sutmlh_vol_ses','smlh_vol_ses','frpolelh_vol_ses','tmpolelh_vol_ses','trvtmlh_vol_ses','insulalh_vol_ses',
                       'banksstsrh_vol_ses','cdacaterh_vol_ses','cdmdfrrh_vol_ses','cuneusrh_vol_ses','ehinalrh_vol_ses',
                       'fusiformrh_vol_ses','ifplrh_vol_ses','iftmrh_vol_ses','ihcaterh_vol_ses','loccrh_vol_ses','lobfrrh_vol_ses',
                       'lingualrh_vol_ses','mobfrrh_vol_ses','mdtmrh_vol_ses','parahpalrh_vol_ses','paracnrh_vol_ses',
                       'parsopcrh_vol_ses','parsobisrh_vol_ses','parstgrisrh_vol_ses','periccrh_vol_ses','postcnrh_vol_ses',
                       'ptcaterh_vol_ses','precnrh_vol_ses','pcrh_vol_ses','rracaterh_vol_ses','rrmdfrrh_vol_ses',
                       'sufrrh_vol_ses','suplrh_vol_ses','sutmrh_vol_ses','smrh_vol_ses','frpolerh_vol_ses',
                       'tmpolerh_vol_ses','trvtmrh_vol_ses','insularh_vol_ses',
                       'banksstslh_vol_pi','cdacatelh_vol_pi','cdmdfrlh_vol_pi','cuneuslh_vol_pi','ehinallh_vol_pi','fusiformlh_vol_pi','ifpllh_vol_pi','iftmlh_vol_pi',
                       'ihcatelh_vol_pi','locclh_vol_pi','lobfrlh_vol_pi','linguallh_vol_pi','mobfrlh_vol_pi','mdtmlh_vol_pi',
                       'parahpallh_vol_pi','paracnlh_vol_pi','parsopclh_vol_pi','parsobislh_vol_pi','parstgrislh_vol_pi',
                       'pericclh_vol_pi','postcnlh_vol_pi','ptcatelh_vol_pi','precnlh_vol_pi','pclh_vol_pi','rracatelh_vol_pi',
                       'rrmdfrlh_vol_pi','sufrlh_vol_pi','supllh_vol_pi','sutmlh_vol_pi','smlh_vol_pi','frpolelh_vol_pi',
                       'tmpolelh_vol_pi','trvtmlh_vol_pi','insulalh_vol_pi','banksstsrh_vol_pi','cdacaterh_vol_pi','cdmdfrrh_vol_pi',
                       'cuneusrh_vol_pi','ehinalrh_vol_pi','fusiformrh_vol_pi','ifplrh_vol_pi','iftmrh_vol_pi','ihcaterh_vol_pi',
                       'loccrh_vol_pi','lobfrrh_vol_pi','lingualrh_vol_pi','mobfrrh_vol_pi','mdtmrh_vol_pi','parahpalrh_vol_pi',
                       'paracnrh_vol_pi','parsopcrh_vol_pi','parsobisrh_vol_pi','parstgrisrh_vol_pi','periccrh_vol_pi',
                       'postcnrh_vol_pi','ptcaterh_vol_pi','precnrh_vol_pi','pcrh_vol_pi','rracaterh_vol_pi','rrmdfrrh_vol_pi',
                       'sufrrh_vol_pi','suplrh_vol_pi','sutmrh_vol_pi','smrh_vol_pi','frpolerh_vol_pi',
                       'tmpolerh_vol_pi','trvtmrh_vol_pi','insularh_vol_pi',
                       'banksstslh_vol_ps','cdacatelh_vol_ps','cdmdfrlh_vol_ps','cuneuslh_vol_ps','ehinallh_vol_ps','fusiformlh_vol_ps','ifpllh_vol_ps','iftmlh_vol_ps',
                       'ihcatelh_vol_ps','locclh_vol_ps','lobfrlh_vol_ps','linguallh_vol_ps','mobfrlh_vol_ps','mdtmlh_vol_ps',
                       'parahpallh_vol_ps','paracnlh_vol_ps','parsopclh_vol_ps','parsobislh_vol_ps','parstgrislh_vol_ps',
                       'pericclh_vol_ps','postcnlh_vol_ps','ptcatelh_vol_ps','precnlh_vol_ps','pclh_vol_ps','rracatelh_vol_ps',
                       'rrmdfrlh_vol_ps','sufrlh_vol_ps','supllh_vol_ps','sutmlh_vol_ps','smlh_vol_ps','frpolelh_vol_ps',
                       'tmpolelh_vol_ps','trvtmlh_vol_ps','insulalh_vol_ps','banksstsrh_vol_ps','cdacaterh_vol_ps','cdmdfrrh_vol_ps',
                       'cuneusrh_vol_ps','ehinalrh_vol_ps','fusiformrh_vol_ps','ifplrh_vol_ps','iftmrh_vol_ps','ihcaterh_vol_ps',
                       'loccrh_vol_ps','lobfrrh_vol_ps','lingualrh_vol_ps','mobfrrh_vol_ps','mdtmrh_vol_ps','parahpalrh_vol_ps',
                       'paracnrh_vol_ps','parsopcrh_vol_ps','parsobisrh_vol_ps','parstgrisrh_vol_ps','periccrh_vol_ps',
                       'postcnrh_vol_ps','ptcaterh_vol_ps','precnrh_vol_ps','pcrh_vol_ps','rracaterh_vol_ps','rrmdfrrh_vol_ps',
                       'sufrrh_vol_ps','suplrh_vol_ps','sutmrh_vol_ps','smrh_vol_ps','frpolerh_vol_ps',
                       'tmpolerh_vol_ps','trvtmrh_vol_ps','insularh_vol_ps')

#calculate 95% CIs and create lower and upper bound variables 
vol_parcel$ci_lower_banksstslh_voli <- vol_parcel$banksstslh_vol_bi - 1.96*vol_parcel$banksstslh_vol_sei 
vol_parcel$ci_upper_banksstslh_voli <- vol_parcel$banksstslh_vol_bi + 1.96*vol_parcel$banksstslh_vol_sei
vol_parcel$ci_lower_cdacatelh_voli <- vol_parcel$cdacatelh_vol_bi - 1.96*vol_parcel$cdacatelh_vol_sei 
vol_parcel$ci_upper_cdacatelh_voli <- vol_parcel$cdacatelh_vol_bi + 1.96*vol_parcel$cdacatelh_vol_sei
vol_parcel$ci_lower_cdmdfrlh_voli <- vol_parcel$cdmdfrlh_vol_bi - 1.96*vol_parcel$cdmdfrlh_vol_sei 
vol_parcel$ci_upper_cdmdfrlh_voli <- vol_parcel$cdmdfrlh_vol_bi + 1.96*vol_parcel$cdmdfrlh_vol_sei
vol_parcel$ci_lower_cuneuslh_voli <- vol_parcel$cuneuslh_vol_bi - 1.96*vol_parcel$cuneuslh_vol_sei 
vol_parcel$ci_upper_cuneuslh_voli <- vol_parcel$cuneuslh_vol_bi + 1.96*vol_parcel$cuneuslh_vol_sei
vol_parcel$ci_lower_ehinallh_voli <- vol_parcel$ehinallh_vol_bi - 1.96*vol_parcel$ehinallh_vol_sei 
vol_parcel$ci_upper_ehinallh_volli <- vol_parcel$ehinallh_vol_bi + 1.96*vol_parcel$ehinallh_vol_sei
vol_parcel$ci_lower_fusiformlh_voli <- vol_parcel$fusiformlh_vol_bi - 1.96*vol_parcel$fusiformlh_vol_sei 
vol_parcel$ci_upper_fusiformlh_voli <- vol_parcel$fusiformlh_vol_bi + 1.96*vol_parcel$fusiformlh_vol_sei
vol_parcel$ci_lower_ifpllh_voli <- vol_parcel$ifpllh_vol_bi - 1.96*vol_parcel$ifpllh_vol_sei 
vol_parcel$ci_upper_ifpllh_voli <- vol_parcel$ifpllh_vol_bi + 1.96*vol_parcel$ifpllh_vol_sei
vol_parcel$ci_lower_iftmlh_voli <- vol_parcel$iftmlh_vol_bi - 1.96*vol_parcel$iftmlh_vol_sei
vol_parcel$ci_upper_iftmlh_voli <- vol_parcel$iftmlh_vol_bi + 1.96*vol_parcel$iftmlh_vol_sei
vol_parcel$ci_lower_ihcatelh_voli <- vol_parcel$ihcatelh_vol_bi - 1.96*vol_parcel$ihcatelh_vol_sei 
vol_parcel$ci_upper_ihcatelh_voli <- vol_parcel$ihcatelh_vol_bi + 1.96*vol_parcel$ihcatelh_vol_sei
vol_parcel$ci_lower_locclh_voli <- vol_parcel$locclh_vol_bi - 1.96*vol_parcel$locclh_vol_sei
vol_parcel$ci_upper_locclh_voli <- vol_parcel$locclh_vol_bi + 1.96*vol_parcel$locclh_vol_sei
vol_parcel$ci_lower_lobfrlh_voli <- vol_parcel$lobfrlh_vol_bi - 1.96*vol_parcel$lobfrlh_vol_sei 
vol_parcel$ci_upper_lobfrlh_voli <- vol_parcel$lobfrlh_vol_bi + 1.96*vol_parcel$lobfrlh_vol_sei
vol_parcel$ci_lower_linguallh_voli <- vol_parcel$linguallh_vol_bi - 1.96*vol_parcel$linguallh_vol_sei 
vol_parcel$ci_upper_linguallh_voli <- vol_parcel$linguallh_vol_bi + 1.96*vol_parcel$linguallh_vol_sei
vol_parcel$ci_lower_mobfrlh_voli <- vol_parcel$mobfrlh_vol_bi - 1.96*vol_parcel$mobfrlh_vol_sei 
vol_parcel$ci_upper_mobfrlh_voli <- vol_parcel$mobfrlh_vol_bi + 1.96*vol_parcel$mobfrlh_vol_sei
vol_parcel$ci_lower_mdtmlh_voli <- vol_parcel$mdtmlh_vol_bi - 1.96*vol_parcel$mdtmlh_vol_sei 
vol_parcel$ci_upper_mdtmlh_voli <- vol_parcel$mdtmlh_vol_bi + 1.96*vol_parcel$mdtmlh_vol_sei
vol_parcel$ci_lower_parahpallh_voli <- vol_parcel$parahpallh_vol_bi - 1.96*vol_parcel$parahpallh_vol_sei 
vol_parcel$ci_upper_parahpallh_voli <- vol_parcel$parahpallh_vol_bi + 1.96*vol_parcel$parahpallh_vol_sei
vol_parcel$ci_lower_paracnlh_voli <- vol_parcel$paracnlh_vol_bi - 1.96*vol_parcel$paracnlh_vol_sei 
vol_parcel$ci_upper_paracnlh_voli <- vol_parcel$paracnlh_vol_bi + 1.96*vol_parcel$paracnlh_vol_sei
vol_parcel$ci_lower_parsopclh_voli <- vol_parcel$parsopclh_vol_bi - 1.96*vol_parcel$parsopclh_vol_sei 
vol_parcel$ci_upper_parsopclh_voli <- vol_parcel$parsopclh_vol_bi + 1.96*vol_parcel$parsopclh_vol_sei
vol_parcel$ci_lower_parsobislh_voli <- vol_parcel$parsobislh_vol_bi - 1.96*vol_parcel$parsobislh_vol_sei 
vol_parcel$ci_upper_parsobislh_voli <- vol_parcel$parsobislh_vol_bi + 1.96*vol_parcel$parsobislh_vol_sei
vol_parcel$ci_lower_parstgrislh_voli <- vol_parcel$parstgrislh_vol_bi - 1.96*vol_parcel$parstgrislh_vol_sei 
vol_parcel$ci_upper_parstgrislh_voli <- vol_parcel$parstgrislh_vol_bi + 1.96*vol_parcel$parstgrislh_vol_sei
vol_parcel$ci_lower_pericclh_voli <- vol_parcel$pericclh_vol_bi - 1.96*vol_parcel$pericclh_vol_sei 
vol_parcel$ci_upper_pericclh_voli <- vol_parcel$pericclh_vol_bi + 1.96*vol_parcel$pericclh_vol_sei
vol_parcel$ci_lower_postcnlh_voli <- vol_parcel$postcnlh_vol_bi - 1.96*vol_parcel$postcnlh_vol_sei 
vol_parcel$ci_upper_postcnlh_voli <- vol_parcel$postcnlh_vol_bi + 1.96*vol_parcel$postcnlh_vol_sei
vol_parcel$ci_lower_ptcatelh_voli <- vol_parcel$ptcatelh_vol_bi - 1.96*vol_parcel$ptcatelh_vol_sei 
vol_parcel$ci_upper_ptcatelh_voli <- vol_parcel$ptcatelh_vol_bi + 1.96*vol_parcel$ptcatelh_vol_sei
vol_parcel$ci_lower_precnlh_voli <- vol_parcel$precnlh_vol_bi - 1.96*vol_parcel$precnlh_vol_sei 
vol_parcel$ci_upper_precnlh_voli <- vol_parcel$precnlh_vol_bi + 1.96*vol_parcel$precnlh_vol_sei
vol_parcel$ci_lower_pclh_voli <- vol_parcel$pclh_vol_bi - 1.96*vol_parcel$pclh_vol_sei 
vol_parcel$ci_upper_pclh_voli <- vol_parcel$pclh_vol_bi + 1.96*vol_parcel$pclh_vol_sei
vol_parcel$ci_lower_rracatelh_voli <- vol_parcel$rracatelh_vol_bi - 1.96*vol_parcel$rracatelh_vol_sei 
vol_parcel$ci_upper_rracatelh_voli <- vol_parcel$rracatelh_vol_bi + 1.96*vol_parcel$rracatelh_vol_sei
vol_parcel$ci_lower_rrmdfrlh_voli <- vol_parcel$rrmdfrlh_vol_bi - 1.96*vol_parcel$rrmdfrlh_vol_sei 
vol_parcel$ci_upper_rrmdfrlh_voli <- vol_parcel$rrmdfrlh_vol_bi + 1.96*vol_parcel$rrmdfrlh_vol_sei
vol_parcel$ci_lower_sufrlh_voli <- vol_parcel$sufrlh_vol_bi - 1.96*vol_parcel$sufrlh_vol_sei 
vol_parcel$ci_upper_sufrlh_voli <- vol_parcel$sufrlh_vol_bi + 1.96*vol_parcel$sufrlh_vol_sei
vol_parcel$ci_lower_supllh_voli <- vol_parcel$supllh_vol_bi - 1.96*vol_parcel$supllh_vol_sei 
vol_parcel$ci_upper_supllh_voli <- vol_parcel$supllh_vol_bi + 1.96*vol_parcel$supllh_vol_sei
vol_parcel$ci_lower_sutmlh_voli <- vol_parcel$sutmlh_vol_bi - 1.96*vol_parcel$sutmlh_vol_sei 
vol_parcel$ci_upper_sutmlh_voli <- vol_parcel$sutmlh_vol_bi + 1.96*vol_parcel$sutmlh_vol_sei
vol_parcel$ci_lower_smlh_voli <- vol_parcel$smlh_vol_bi - 1.96*vol_parcel$smlh_vol_sei 
vol_parcel$ci_upper_smlh_voli <- vol_parcel$smlh_vol_bi + 1.96*vol_parcel$smlh_vol_sei
vol_parcel$ci_lower_frpolelh_voli <- vol_parcel$frpolelh_vol_bi - 1.96*vol_parcel$frpolelh_vol_sei 
vol_parcel$ci_upper_frpolelh_voli <- vol_parcel$frpolelh_vol_bi + 1.96*vol_parcel$frpolelh_vol_sei
vol_parcel$ci_lower_tmpolelh_voli <- vol_parcel$tmpolelh_vol_bi - 1.96*vol_parcel$tmpolelh_vol_sei 
vol_parcel$ci_upper_tmpolelh_voli <- vol_parcel$tmpolelh_vol_bi + 1.96*vol_parcel$tmpolelh_vol_sei
vol_parcel$ci_lower_trvtmlh_voli <- vol_parcel$trvtmlh_vol_bi - 1.96*vol_parcel$trvtmlh_vol_sei 
vol_parcel$ci_upper_trvtmlh_voli <- vol_parcel$trvtmlh_vol_bi + 1.96*vol_parcel$trvtmlh_vol_sei
vol_parcel$ci_lower_insulalh_voli <- vol_parcel$insulalh_vol_bi - 1.96*vol_parcel$insulalh_vol_sei 
vol_parcel$ci_upper_insulalh_voli <- vol_parcel$insulalh_vol_bi + 1.96*vol_parcel$insulalh_vol_sei

vol_parcel$ci_lower_banksstsrh_voli <- vol_parcel$banksstsrh_vol_bi - 1.96*vol_parcel$banksstsrh_vol_sei 
vol_parcel$ci_upper_banksstsrh_voli <- vol_parcel$banksstsrh_vol_bi + 1.96*vol_parcel$banksstsrh_vol_sei
vol_parcel$ci_lower_cdacaterh_voli <- vol_parcel$cdacaterh_vol_bi - 1.96*vol_parcel$cdacaterh_vol_sei 
vol_parcel$ci_upper_cdacaterh_voli <- vol_parcel$cdacaterh_vol_bi + 1.96*vol_parcel$cdacaterh_vol_sei
vol_parcel$ci_lower_cdmdfrrh_voli <- vol_parcel$cdmdfrrh_vol_bi - 1.96*vol_parcel$cdmdfrrh_vol_sei 
vol_parcel$ci_upper_cdmdfrrh_voli <- vol_parcel$cdmdfrrh_vol_bi + 1.96*vol_parcel$cdmdfrrh_vol_sei
vol_parcel$ci_lower_cuneusrh_voli <- vol_parcel$cuneusrh_vol_bi - 1.96*vol_parcel$cuneusrh_vol_sei 
vol_parcel$ci_upper_cuneusrh_voli <- vol_parcel$cuneusrh_vol_bi + 1.96*vol_parcel$cuneusrh_vol_sei
vol_parcel$ci_lower_ehinalrh_voli <- vol_parcel$ehinalrh_vol_bi - 1.96*vol_parcel$ehinalrh_vol_sei 
vol_parcel$ci_upper_ehinalrh_volli <- vol_parcel$ehinalrh_vol_bi + 1.96*vol_parcel$ehinalrh_vol_sei
vol_parcel$ci_lower_fusiformrh_voli <- vol_parcel$fusiformrh_vol_bi - 1.96*vol_parcel$fusiformrh_vol_sei 
vol_parcel$ci_upper_fusiformrh_voli <- vol_parcel$fusiformrh_vol_bi + 1.96*vol_parcel$fusiformrh_vol_sei
vol_parcel$ci_lower_ifplrh_voli <- vol_parcel$ifplrh_vol_bi - 1.96*vol_parcel$ifplrh_vol_sei 
vol_parcel$ci_upper_ifplrh_voli <- vol_parcel$ifplrh_vol_bi + 1.96*vol_parcel$ifplrh_vol_sei
vol_parcel$ci_lower_iftmrh_voli <- vol_parcel$iftmrh_vol_bi - 1.96*vol_parcel$iftmrh_vol_sei 
vol_parcel$ci_upper_iftmrh_voli <- vol_parcel$iftmrh_vol_bi + 1.96*vol_parcel$iftmrh_vol_sei
vol_parcel$ci_lower_ihcaterh_voli <- vol_parcel$ihcaterh_vol_bi - 1.96*vol_parcel$ihcaterh_vol_sei 
vol_parcel$ci_upper_ihcaterh_voli <- vol_parcel$ihcaterh_vol_bi + 1.96*vol_parcel$ihcaterh_vol_sei
vol_parcel$ci_lower_loccrh_voli <- vol_parcel$loccrh_vol_bi - 1.96*vol_parcel$loccrh_vol_sei 
vol_parcel$ci_upper_loccrh_voli <- vol_parcel$loccrh_vol_bi + 1.96*vol_parcel$loccrh_vol_sei
vol_parcel$ci_lower_lobfrrh_voli <- vol_parcel$lobfrrh_vol_bi - 1.96*vol_parcel$lobfrrh_vol_sei 
vol_parcel$ci_upper_lobfrrh_voli <- vol_parcel$lobfrrh_vol_bi + 1.96*vol_parcel$lobfrrh_vol_sei
vol_parcel$ci_lower_lingualrh_voli <- vol_parcel$lingualrh_vol_bi - 1.96*vol_parcel$lingualrh_vol_sei 
vol_parcel$ci_upper_lingualrh_voli <- vol_parcel$lingualrh_vol_bi + 1.96*vol_parcel$lingualrh_vol_sei
vol_parcel$ci_lower_mobfrrh_voli <- vol_parcel$mobfrrh_vol_bi - 1.96*vol_parcel$mobfrrh_vol_sei 
vol_parcel$ci_upper_mobfrrh_voli <- vol_parcel$mobfrrh_vol_bi + 1.96*vol_parcel$mobfrrh_vol_sei
vol_parcel$ci_lower_mdtmrh_voli <- vol_parcel$mdtmrh_vol_bi - 1.96*vol_parcel$mdtmrh_vol_sei 
vol_parcel$ci_upper_mdtmrh_voli <- vol_parcel$mdtmrh_vol_bi + 1.96*vol_parcel$mdtmrh_vol_sei
vol_parcel$ci_lower_parahpalrh_voli <- vol_parcel$parahpalrh_vol_bi - 1.96*vol_parcel$parahpalrh_vol_sei 
vol_parcel$ci_upper_parahpalrh_voli <- vol_parcel$parahpalrh_vol_bi + 1.96*vol_parcel$parahpalrh_vol_sei
vol_parcel$ci_lower_paracnrh_voli <- vol_parcel$paracnrh_vol_bi - 1.96*vol_parcel$paracnrh_vol_sei 
vol_parcel$ci_upper_paracnrh_voli <- vol_parcel$paracnrh_vol_bi + 1.96*vol_parcel$paracnrh_vol_sei
vol_parcel$ci_lower_parsopcrh_voli <- vol_parcel$parsopcrh_vol_bi - 1.96*vol_parcel$parsopcrh_vol_sei 
vol_parcel$ci_upper_parsopcrh_voli <- vol_parcel$parsopcrh_vol_bi + 1.96*vol_parcel$parsopcrh_vol_sei
vol_parcel$ci_lower_parsobisrh_voli <- vol_parcel$parsobisrh_vol_bi - 1.96*vol_parcel$parsobisrh_vol_sei 
vol_parcel$ci_upper_parsobisrh_voli <- vol_parcel$parsobisrh_vol_bi + 1.96*vol_parcel$parsobisrh_vol_sei
vol_parcel$ci_lower_parstgrisrh_voli <- vol_parcel$parstgrisrh_vol_bi - 1.96*vol_parcel$parstgrisrh_vol_sei 
vol_parcel$ci_upper_parstgrisrh_voli <- vol_parcel$parstgrisrh_vol_bi + 1.96*vol_parcel$parstgrisrh_vol_sei
vol_parcel$ci_lower_periccrh_voli <- vol_parcel$periccrh_vol_bi - 1.96*vol_parcel$periccrh_vol_sei 
vol_parcel$ci_upper_periccrh_voli <- vol_parcel$periccrh_vol_bi + 1.96*vol_parcel$periccrh_vol_sei
vol_parcel$ci_lower_postcnrh_voli <- vol_parcel$postcnrh_vol_bi - 1.96*vol_parcel$postcnrh_vol_sei 
vol_parcel$ci_upper_postcnrh_voli <- vol_parcel$postcnrh_vol_bi + 1.96*vol_parcel$postcnrh_vol_sei
vol_parcel$ci_lower_ptcaterh_voli <- vol_parcel$ptcaterh_vol_bi - 1.96*vol_parcel$ptcaterh_vol_sei 
vol_parcel$ci_upper_ptcaterh_voli <- vol_parcel$ptcaterh_vol_bi + 1.96*vol_parcel$ptcaterh_vol_sei
vol_parcel$ci_lower_precnrh_voli <- vol_parcel$precnrh_vol_bi - 1.96*vol_parcel$precnrh_vol_sei 
vol_parcel$ci_upper_precnrh_voli <- vol_parcel$precnrh_vol_bi + 1.96*vol_parcel$precnrh_vol_sei
vol_parcel$ci_lower_pcrh_voli <- vol_parcel$pcrh_vol_bi - 1.96*vol_parcel$pcrh_vol_sei 
vol_parcel$ci_upper_pcrh_voli <- vol_parcel$pcrh_vol_bi + 1.96*vol_parcel$pcrh_vol_sei
vol_parcel$ci_lower_rracaterh_voli <- vol_parcel$rracaterh_vol_bi - 1.96*vol_parcel$rracaterh_vol_sei 
vol_parcel$ci_upper_rracaterh_voli <- vol_parcel$rracaterh_vol_bi + 1.96*vol_parcel$rracaterh_vol_sei
vol_parcel$ci_lower_rrmdfrrh_voli <- vol_parcel$rrmdfrrh_vol_bi - 1.96*vol_parcel$rrmdfrrh_vol_sei 
vol_parcel$ci_upper_rrmdfrrh_voli <- vol_parcel$rrmdfrrh_vol_bi + 1.96*vol_parcel$rrmdfrrh_vol_sei
vol_parcel$ci_lower_sufrrh_voli <- vol_parcel$sufrrh_vol_bi - 1.96*vol_parcel$sufrrh_vol_sei 
vol_parcel$ci_upper_sufrrh_voli <- vol_parcel$sufrrh_vol_bi + 1.96*vol_parcel$sufrrh_vol_sei
vol_parcel$ci_lower_suplrh_voli <- vol_parcel$suplrh_vol_bi - 1.96*vol_parcel$suplrh_vol_sei 
vol_parcel$ci_upper_suplrh_voli <- vol_parcel$suplrh_vol_bi + 1.96*vol_parcel$suplrh_vol_sei
vol_parcel$ci_lower_sutmrh_voli <- vol_parcel$sutmrh_vol_bi - 1.96*vol_parcel$sutmrh_vol_sei 
vol_parcel$ci_upper_sutmrh_voli <- vol_parcel$sutmrh_vol_bi + 1.96*vol_parcel$sutmrh_vol_sei
vol_parcel$ci_lower_smrh_voli <- vol_parcel$smrh_vol_bi - 1.96*vol_parcel$smrh_vol_sei 
vol_parcel$ci_upper_smrh_voli <- vol_parcel$smrh_vol_bi + 1.96*vol_parcel$smrh_vol_sei
vol_parcel$ci_lower_frpolerh_voli <- vol_parcel$frpolerh_vol_bi - 1.96*vol_parcel$frpolerh_vol_sei 
vol_parcel$ci_upper_frpolerh_voli <- vol_parcel$frpolerh_vol_bi + 1.96*vol_parcel$frpolerh_vol_sei
vol_parcel$ci_lower_tmpolerh_voli <- vol_parcel$tmpolerh_vol_bi - 1.96*vol_parcel$tmpolerh_vol_sei 
vol_parcel$ci_upper_tmpolerh_voli <- vol_parcel$tmpolerh_vol_bi + 1.96*vol_parcel$tmpolerh_vol_sei
vol_parcel$ci_lower_trvtmrh_voli <- vol_parcel$trvtmrh_vol_bi - 1.96*vol_parcel$trvtmrh_vol_sei 
vol_parcel$ci_upper_trvtmrh_voli <- vol_parcel$trvtmrh_vol_bi + 1.96*vol_parcel$trvtmrh_vol_sei
vol_parcel$ci_lower_insularh_voli <- vol_parcel$insularh_vol_bi - 1.96*vol_parcel$insularh_vol_sei 
vol_parcel$ci_upper_insularh_voli <- vol_parcel$insularh_vol_bi + 1.96*vol_parcel$insularh_vol_sei

write.csv(vol_parcel, "Parcel-Wise Cortical Volume MLM Analysis Output Standardized_FINAL.csv") #write csv file


#Parcel-wise cortical area analyses

area1 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_banksstslh_area_cent + scale(subj_banksstslh_area) + site_banksstslh_area_cent*wave + scale(subj_banksstslh_area)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b1i <- as.data.frame(c(summary(area1)$coefficients[14,1]))
area_se1i <- as.data.frame(c(summary(area1)$coefficients[14,2]))
area_p1i <- as.data.frame(c(summary(area1)$coefficients[14,5]))
area_b1s <- as.data.frame(c(summary(area1)$coefficients[16,1]))
area_se1s <- as.data.frame(c(summary(area1)$coefficients[16,2]))
area_p1s <- as.data.frame(c(summary(area1)$coefficients[16,5]))

area2 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_cdacatelh_area_cent + scale(subj_cdacatelh_area) + site_cdacatelh_area_cent*wave + scale(subj_cdacatelh_area)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b2i <- as.data.frame(c(summary(area2)$coefficients[14,1]))
area_se2i <- as.data.frame(c(summary(area2)$coefficients[14,2]))
area_p2i <- as.data.frame(c(summary(area2)$coefficients[14,5]))
area_b2s <- as.data.frame(c(summary(area2)$coefficients[16,1]))
area_se2s <- as.data.frame(c(summary(area2)$coefficients[16,2]))
area_p2s <- as.data.frame(c(summary(area2)$coefficients[16,5]))

area3 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_cdmdfrlh_area_cent + scale(subj_cdmdfrlh_area) + site_cdmdfrlh_area_cent*wave + scale(subj_cdmdfrlh_area)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b3i <- as.data.frame(c(summary(area3)$coefficients[14,1]))
area_se3i <- as.data.frame(c(summary(area3)$coefficients[14,2]))
area_p3i <- as.data.frame(c(summary(area3)$coefficients[14,5]))
area_b3s <- as.data.frame(c(summary(area3)$coefficients[16,1]))
area_se3s <- as.data.frame(c(summary(area3)$coefficients[16,2]))
area_p3s <- as.data.frame(c(summary(area3)$coefficients[16,5]))

area4 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_cuneuslh_area_cent + scale(subj_cuneuslh_area) + site_cuneuslh_area_cent*wave + scale(subj_cuneuslh_area)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b4i <- as.data.frame(c(summary(area4)$coefficients[14,1]))
area_se4i <- as.data.frame(c(summary(area4)$coefficients[14,2]))
area_p4i <- as.data.frame(c(summary(area4)$coefficients[14,5]))
area_b4s <- as.data.frame(c(summary(area4)$coefficients[16,1]))
area_se4s <- as.data.frame(c(summary(area4)$coefficients[16,2]))
area_p4s <- as.data.frame(c(summary(area4)$coefficients[16,5]))

area5 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_ehinallh_area_cent + scale(subj_ehinallh_area) + site_ehinallh_area_cent*wave + scale(subj_ehinallh_area)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b5i <- as.data.frame(c(summary(area5)$coefficients[14,1]))
area_se5i <- as.data.frame(c(summary(area5)$coefficients[14,2]))
area_p5i <- as.data.frame(c(summary(area5)$coefficients[14,5]))
area_b5s <- as.data.frame(c(summary(area5)$coefficients[16,1]))
area_se5s <- as.data.frame(c(summary(area5)$coefficients[16,2]))
area_p5s <- as.data.frame(c(summary(area5)$coefficients[16,5]))

area6 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_fusiformlh_area_cent + scale(subj_fusiformlh_area) + site_fusiformlh_area_cent*wave + scale(subj_fusiformlh_area)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b6i <- as.data.frame(c(summary(area6)$coefficients[14,1]))
area_se6i <- as.data.frame(c(summary(area6)$coefficients[14,2]))
area_p6i <- as.data.frame(c(summary(area6)$coefficients[14,5]))
area_b6s <- as.data.frame(c(summary(area6)$coefficients[16,1]))
area_se6s <- as.data.frame(c(summary(area6)$coefficients[16,2]))
area_p6s <- as.data.frame(c(summary(area6)$coefficients[16,5]))

area7 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_ifpllh_area_cent + scale(subj_ifpllh_area) + site_ifpllh_area_cent*wave + scale(subj_ifpllh_area)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b7i <- as.data.frame(c(summary(area7)$coefficients[14,1]))
area_se7i <- as.data.frame(c(summary(area7)$coefficients[14,2]))
area_p7i <- as.data.frame(c(summary(area7)$coefficients[14,5]))
area_b7s <- as.data.frame(c(summary(area7)$coefficients[16,1]))
area_se7s <- as.data.frame(c(summary(area7)$coefficients[16,2]))
area_p7s <- as.data.frame(c(summary(area7)$coefficients[16,5]))

area8 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_iftmlh_area_cent + scale(subj_iftmlh_area) + site_iftmlh_area_cent*wave + scale(subj_iftmlh_area)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b8i <- as.data.frame(c(summary(area8)$coefficients[14,1]))
area_se8i <- as.data.frame(c(summary(area8)$coefficients[14,2]))
area_p8i <- as.data.frame(c(summary(area8)$coefficients[14,5]))
area_b8s <- as.data.frame(c(summary(area8)$coefficients[16,1]))
area_se8s <- as.data.frame(c(summary(area8)$coefficients[16,2]))
area_p8s <- as.data.frame(c(summary(area8)$coefficients[16,5]))

area9 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_ihcatelh_area_cent + scale(subj_ihcatelh_area) + site_ihcatelh_area_cent*wave + scale(subj_ihcatelh_area)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b9i <- as.data.frame(c(summary(area9)$coefficients[14,1]))
area_se9i <- as.data.frame(c(summary(area9)$coefficients[14,2]))
area_p9i <- as.data.frame(c(summary(area9)$coefficients[14,5]))
area_b9s <- as.data.frame(c(summary(area9)$coefficients[16,1]))
area_se9s <- as.data.frame(c(summary(area9)$coefficients[16,2]))
area_p9s <- as.data.frame(c(summary(area9)$coefficients[16,5]))

area10 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_locclh_area_cent + scale(subj_locclh_area) + site_locclh_area_cent*wave + scale(subj_locclh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b10i <- as.data.frame(c(summary(area10)$coefficients[14,1]))
area_se10i <- as.data.frame(c(summary(area10)$coefficients[14,2]))
area_p10i <- as.data.frame(c(summary(area10)$coefficients[14,5]))
area_b10s <- as.data.frame(c(summary(area10)$coefficients[16,1]))
area_se10s <- as.data.frame(c(summary(area10)$coefficients[16,2]))
area_p10s <- as.data.frame(c(summary(area10)$coefficients[16,5]))

area11 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_lobfrlh_area_cent + scale(subj_lobfrlh_area) + site_lobfrlh_area_cent*wave + scale(subj_lobfrlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b11i <- as.data.frame(c(summary(area11)$coefficients[14,1]))
area_se11i <- as.data.frame(c(summary(area11)$coefficients[14,2]))
area_p11i <- as.data.frame(c(summary(area11)$coefficients[14,5]))
area_b11s <- as.data.frame(c(summary(area11)$coefficients[16,1]))
area_se11s <- as.data.frame(c(summary(area11)$coefficients[16,2]))
area_p11s <- as.data.frame(c(summary(area11)$coefficients[16,5]))

area12 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_linguallh_area_cent + scale(subj_linguallh_area) + site_linguallh_area_cent*wave + scale(subj_linguallh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b12i <- as.data.frame(c(summary(area12)$coefficients[14,1]))
area_se12i <- as.data.frame(c(summary(area12)$coefficients[14,2]))
area_p12i <- as.data.frame(c(summary(area12)$coefficients[14,5]))
area_b12s <- as.data.frame(c(summary(area12)$coefficients[16,1]))
area_se12s <- as.data.frame(c(summary(area12)$coefficients[16,2]))
area_p12s <- as.data.frame(c(summary(area12)$coefficients[16,5]))

area13 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_mobfrlh_area_cent + scale(subj_mobfrlh_area) + site_mobfrlh_area_cent*wave + scale(subj_mobfrlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b13i <- as.data.frame(c(summary(area13)$coefficients[14,1]))
area_se13i <- as.data.frame(c(summary(area13)$coefficients[14,2]))
area_p13i <- as.data.frame(c(summary(area13)$coefficients[14,5]))
area_b13s <- as.data.frame(c(summary(area13)$coefficients[16,1]))
area_se13s <- as.data.frame(c(summary(area13)$coefficients[16,2]))
area_p13s <- as.data.frame(c(summary(area13)$coefficients[16,5]))

area14 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_mdtmlh_area_cent + scale(subj_mdtmlh_area) + site_mdtmlh_area_cent*wave + scale(subj_mdtmlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b14i <- as.data.frame(c(summary(area14)$coefficients[14,1]))
area_se14i <- as.data.frame(c(summary(area14)$coefficients[14,2]))
area_p14i <- as.data.frame(c(summary(area14)$coefficients[14,5]))
area_b14s <- as.data.frame(c(summary(area14)$coefficients[16,1]))
area_se14s <- as.data.frame(c(summary(area14)$coefficients[16,2]))
area_p14s <- as.data.frame(c(summary(area14)$coefficients[16,5]))

area15 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_parahpallh_area_cent + scale(subj_parahpallh_area) + site_parahpallh_area_cent*wave + scale(subj_parahpallh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b15i <- as.data.frame(c(summary(area15)$coefficients[14,1]))
area_se15i <- as.data.frame(c(summary(area15)$coefficients[14,2]))
area_p15i <- as.data.frame(c(summary(area15)$coefficients[14,5]))
area_b15s <- as.data.frame(c(summary(area15)$coefficients[16,1]))
area_se15s <- as.data.frame(c(summary(area15)$coefficients[16,2]))
area_p15s <- as.data.frame(c(summary(area15)$coefficients[16,5]))

area16 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_paracnlh_area_cent + scale(subj_paracnlh_area) + site_paracnlh_area_cent*wave + scale(subj_paracnlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b16i <- as.data.frame(c(summary(area16)$coefficients[14,1]))
area_se16i <- as.data.frame(c(summary(area16)$coefficients[14,2]))
area_p16i <- as.data.frame(c(summary(area16)$coefficients[14,5]))
area_b16s <- as.data.frame(c(summary(area16)$coefficients[16,1]))
area_se16s <- as.data.frame(c(summary(area16)$coefficients[16,2]))
area_p16s <- as.data.frame(c(summary(area16)$coefficients[16,5]))

area17 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_parsopclh_area_cent + scale(subj_parsopclh_area) + site_parsopclh_area_cent*wave + scale(subj_parsopclh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b17i <- as.data.frame(c(summary(area17)$coefficients[14,1]))
area_se17i <- as.data.frame(c(summary(area17)$coefficients[14,2]))
area_p17i <- as.data.frame(c(summary(area17)$coefficients[14,5]))
area_b17s <- as.data.frame(c(summary(area17)$coefficients[16,1]))
area_se17s <- as.data.frame(c(summary(area17)$coefficients[16,2]))
area_p17s <- as.data.frame(c(summary(area17)$coefficients[16,5]))

area18 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_parsobislh_area_cent + scale(subj_parsobislh_area) + site_parsobislh_area_cent*wave + scale(subj_parsobislh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b18i <- as.data.frame(c(summary(area18)$coefficients[14,1]))
area_se18i <- as.data.frame(c(summary(area18)$coefficients[14,2]))
area_p18i <- as.data.frame(c(summary(area18)$coefficients[14,5]))
area_b18s <- as.data.frame(c(summary(area18)$coefficients[16,1]))
area_se18s <- as.data.frame(c(summary(area18)$coefficients[16,2]))
area_p18s <- as.data.frame(c(summary(area18)$coefficients[16,5]))

area19 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_parstgrislh_area_cent + scale(subj_parstgrislh_area) + site_parstgrislh_area_cent*wave + scale(subj_parstgrislh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b19i <- as.data.frame(c(summary(area19)$coefficients[14,1]))
area_se19i <- as.data.frame(c(summary(area19)$coefficients[14,2]))
area_p19i <- as.data.frame(c(summary(area19)$coefficients[14,5]))
area_b19s <- as.data.frame(c(summary(area19)$coefficients[16,1]))
area_se19s <- as.data.frame(c(summary(area19)$coefficients[16,2]))
area_p19s <- as.data.frame(c(summary(area19)$coefficients[16,5]))


area20 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_pericclh_area_cent + scale(subj_pericclh_area) + site_pericclh_area_cent*wave + scale(subj_pericclh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b20i <- as.data.frame(c(summary(area20)$coefficients[14,1]))
area_se20i <- as.data.frame(c(summary(area20)$coefficients[14,2]))
area_p20i <- as.data.frame(c(summary(area20)$coefficients[14,5]))
area_b20s <- as.data.frame(c(summary(area20)$coefficients[16,1]))
area_se20s <- as.data.frame(c(summary(area20)$coefficients[16,2]))
area_p20s <- as.data.frame(c(summary(area20)$coefficients[16,5]))

area21 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_postcnlh_area_cent + scale(subj_postcnlh_area) + site_postcnlh_area_cent*wave + scale(subj_postcnlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b21i <- as.data.frame(c(summary(area21)$coefficients[14,1]))
area_se21i <- as.data.frame(c(summary(area21)$coefficients[14,2]))
area_p21i <- as.data.frame(c(summary(area21)$coefficients[14,5]))
area_b21s <- as.data.frame(c(summary(area21)$coefficients[16,1]))
area_se21s <- as.data.frame(c(summary(area21)$coefficients[16,2]))
area_p21s <- as.data.frame(c(summary(area21)$coefficients[16,5]))

area22 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_ptcatelh_area_cent + scale(subj_ptcatelh_area) + site_ptcatelh_area_cent*wave + scale(subj_ptcatelh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b22i <- as.data.frame(c(summary(area22)$coefficients[14,1]))
area_se22i <- as.data.frame(c(summary(area22)$coefficients[14,2]))
area_p22i <- as.data.frame(c(summary(area22)$coefficients[14,5]))
area_b22s <- as.data.frame(c(summary(area22)$coefficients[16,1]))
area_se22s <- as.data.frame(c(summary(area22)$coefficients[16,2]))
area_p22s <- as.data.frame(c(summary(area22)$coefficients[16,5]))

area23 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_precnlh_area_cent + scale(subj_precnlh_area) + site_precnlh_area_cent*wave + scale(subj_precnlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b23i <- as.data.frame(c(summary(area23)$coefficients[14,1]))
area_se23i <- as.data.frame(c(summary(area23)$coefficients[14,2]))
area_p23i <- as.data.frame(c(summary(area23)$coefficients[14,5]))
area_b23s <- as.data.frame(c(summary(area23)$coefficients[16,1]))
area_se23s <- as.data.frame(c(summary(area23)$coefficients[16,2]))
area_p23s <- as.data.frame(c(summary(area23)$coefficients[16,5]))

area24 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_pclh_area_cent + scale(subj_pclh_area) + site_pclh_area_cent*wave + scale(subj_pclh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b24i <- as.data.frame(c(summary(area24)$coefficients[14,1]))
area_se24i <- as.data.frame(c(summary(area24)$coefficients[14,2]))
area_p24i <- as.data.frame(c(summary(area24)$coefficients[14,5]))
area_b24s <- as.data.frame(c(summary(area24)$coefficients[16,1]))
area_se24s <- as.data.frame(c(summary(area24)$coefficients[16,2]))
area_p24s <- as.data.frame(c(summary(area24)$coefficients[16,5]))

area25 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_rracatelh_area_cent + scale(subj_rracatelh_area) + site_rracatelh_area_cent*wave + scale(subj_rracatelh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b25i <- as.data.frame(c(summary(area25)$coefficients[14,1]))
area_se25i <- as.data.frame(c(summary(area25)$coefficients[14,2]))
area_p25i <- as.data.frame(c(summary(area25)$coefficients[14,5]))
area_b25s <- as.data.frame(c(summary(area25)$coefficients[16,1]))
area_se25s <- as.data.frame(c(summary(area25)$coefficients[16,2]))
area_p25s <- as.data.frame(c(summary(area25)$coefficients[16,5]))

area26 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_rrmdfrlh_area_cent + scale(subj_rrmdfrlh_area) + site_rrmdfrlh_area_cent*wave + scale(subj_rrmdfrlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b26i <- as.data.frame(c(summary(area26)$coefficients[14,1]))
area_se26i <- as.data.frame(c(summary(area26)$coefficients[14,2]))
area_p26i <- as.data.frame(c(summary(area26)$coefficients[14,5]))
area_b26s <- as.data.frame(c(summary(area26)$coefficients[16,1]))
area_se26s <- as.data.frame(c(summary(area26)$coefficients[16,2]))
area_p26s <- as.data.frame(c(summary(area26)$coefficients[16,5]))

area27 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_sufrlh_area_cent + scale(subj_sufrlh_area) + site_sufrlh_area_cent*wave + scale(subj_sufrlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b27i <- as.data.frame(c(summary(area27)$coefficients[14,1]))
area_se27i <- as.data.frame(c(summary(area27)$coefficients[14,2]))
area_p27i <- as.data.frame(c(summary(area27)$coefficients[14,5]))
area_b27s <- as.data.frame(c(summary(area27)$coefficients[16,1]))
area_se27s <- as.data.frame(c(summary(area27)$coefficients[16,2]))
area_p27s <- as.data.frame(c(summary(area27)$coefficients[16,5]))

area28 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_supllh_area_cent + scale(subj_supllh_area) + site_supllh_area_cent*wave + scale(subj_supllh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b28i <- as.data.frame(c(summary(area28)$coefficients[14,1]))
area_se28i <- as.data.frame(c(summary(area28)$coefficients[14,2]))
area_p28i <- as.data.frame(c(summary(area28)$coefficients[14,5]))
area_b28s <- as.data.frame(c(summary(area28)$coefficients[16,1]))
area_se28s <- as.data.frame(c(summary(area28)$coefficients[16,2]))
area_p28s <- as.data.frame(c(summary(area28)$coefficients[16,5]))

area29 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_sutmlh_area_cent + scale(subj_sutmlh_area) + site_sutmlh_area_cent*wave + scale(subj_sutmlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b29i <- as.data.frame(c(summary(area29)$coefficients[14,1]))
area_se29i <- as.data.frame(c(summary(area29)$coefficients[14,2]))
area_p29i <- as.data.frame(c(summary(area29)$coefficients[14,5]))
area_b29s <- as.data.frame(c(summary(area29)$coefficients[16,1]))
area_se29s <- as.data.frame(c(summary(area29)$coefficients[16,2]))
area_p29s <- as.data.frame(c(summary(area29)$coefficients[16,5]))

area30 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_smlh_area_cent + scale(subj_smlh_area) + site_smlh_area_cent*wave + scale(subj_smlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b30i <- as.data.frame(c(summary(area30)$coefficients[14,1]))
area_se30i <- as.data.frame(c(summary(area30)$coefficients[14,2]))
area_p30i <- as.data.frame(c(summary(area30)$coefficients[14,5]))
area_b30s <- as.data.frame(c(summary(area30)$coefficients[16,1]))
area_se30s <- as.data.frame(c(summary(area30)$coefficients[16,2]))
area_p30s <- as.data.frame(c(summary(area30)$coefficients[16,5]))

area31 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_frpolelh_area_cent + scale(subj_frpolelh_area) + site_frpolelh_area_cent*wave + scale(subj_frpolelh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b31i <- as.data.frame(c(summary(area31)$coefficients[14,1]))
area_se31i <- as.data.frame(c(summary(area31)$coefficients[14,2]))
area_p31i <- as.data.frame(c(summary(area31)$coefficients[14,5]))
area_b31s <- as.data.frame(c(summary(area31)$coefficients[16,1]))
area_se31s <- as.data.frame(c(summary(area31)$coefficients[16,2]))
area_p31s <- as.data.frame(c(summary(area31)$coefficients[16,5]))

area32 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_tmpolelh_area_cent + scale(subj_tmpolelh_area) + site_tmpolelh_area_cent*wave + scale(subj_tmpolelh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b32i <- as.data.frame(c(summary(area32)$coefficients[14,1]))
area_se32i <- as.data.frame(c(summary(area32)$coefficients[14,2]))
area_p32i <- as.data.frame(c(summary(area32)$coefficients[14,5]))
area_b32s <- as.data.frame(c(summary(area32)$coefficients[16,1]))
area_se32s <- as.data.frame(c(summary(area32)$coefficients[16,2]))
area_p32s <- as.data.frame(c(summary(area32)$coefficients[16,5]))

area33 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_trvtmlh_area_cent + scale(subj_trvtmlh_area) + site_trvtmlh_area_cent*wave + scale(subj_trvtmlh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b33i <- as.data.frame(c(summary(area33)$coefficients[14,1]))
area_se33i <- as.data.frame(c(summary(area33)$coefficients[14,2]))
area_p33i <- as.data.frame(c(summary(area33)$coefficients[14,5]))
area_b33s <- as.data.frame(c(summary(area33)$coefficients[16,1]))
area_se33s <- as.data.frame(c(summary(area33)$coefficients[16,2]))
area_p33s <- as.data.frame(c(summary(area33)$coefficients[16,5]))

area34 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_insulalh_area_cent + scale(subj_insulalh_area) + site_insulalh_area_cent*wave + scale(subj_insulalh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b34i <- as.data.frame(c(summary(area34)$coefficients[14,1]))
area_se34i <- as.data.frame(c(summary(area34)$coefficients[14,2]))
area_p34i <- as.data.frame(c(summary(area34)$coefficients[14,5]))
area_b34s <- as.data.frame(c(summary(area34)$coefficients[16,1]))
area_se34s <- as.data.frame(c(summary(area34)$coefficients[16,2]))
area_p34s <- as.data.frame(c(summary(area34)$coefficients[16,5]))

area35 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_banksstsrh_area_cent + scale(subj_banksstsrh_area) + site_banksstsrh_area_cent*wave + scale(subj_banksstsrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b35i <- as.data.frame(c(summary(area35)$coefficients[14,1]))
area_se35i <- as.data.frame(c(summary(area35)$coefficients[14,2]))
area_p35i <- as.data.frame(c(summary(area35)$coefficients[14,5]))
area_b35s <- as.data.frame(c(summary(area35)$coefficients[16,1]))
area_se35s <- as.data.frame(c(summary(area35)$coefficients[16,2]))
area_p35s <- as.data.frame(c(summary(area35)$coefficients[16,5]))

area36 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_cdacaterh_area_cent + scale(subj_cdacaterh_area) + site_cdacaterh_area_cent*wave + scale(subj_cdacaterh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b36i <- as.data.frame(c(summary(area36)$coefficients[14,1]))
area_se36i <- as.data.frame(c(summary(area36)$coefficients[14,2]))
area_p36i <- as.data.frame(c(summary(area36)$coefficients[14,5]))
area_b36s <- as.data.frame(c(summary(area36)$coefficients[16,1]))
area_se36s <- as.data.frame(c(summary(area36)$coefficients[16,2]))
area_p36s <- as.data.frame(c(summary(area36)$coefficients[16,5]))

area37 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_cdmdfrrh_area_cent + scale(subj_cdmdfrrh_area) + site_cdmdfrrh_area_cent*wave + scale(subj_cdmdfrrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b37i <- as.data.frame(c(summary(area37)$coefficients[14,1]))
area_se37i <- as.data.frame(c(summary(area37)$coefficients[14,2]))
area_p37i <- as.data.frame(c(summary(area37)$coefficients[14,5]))
area_b37s <- as.data.frame(c(summary(area37)$coefficients[16,1]))
area_se37s <- as.data.frame(c(summary(area37)$coefficients[16,2]))
area_p37s <- as.data.frame(c(summary(area37)$coefficients[16,5]))

area38 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_cuneusrh_area_cent + scale(subj_cuneusrh_area) + site_cuneusrh_area_cent*wave + scale(subj_cuneusrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b38i <- as.data.frame(c(summary(area38)$coefficients[14,1]))
area_se38i <- as.data.frame(c(summary(area38)$coefficients[14,2]))
area_p38i <- as.data.frame(c(summary(area38)$coefficients[14,5]))
area_b38s <- as.data.frame(c(summary(area38)$coefficients[16,1]))
area_se38s <- as.data.frame(c(summary(area38)$coefficients[16,2]))
area_p38s <- as.data.frame(c(summary(area38)$coefficients[16,5]))

area39 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_ehinalrh_area_cent + scale(subj_ehinalrh_area) + site_ehinalrh_area_cent*wave + scale(subj_ehinalrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b39i <- as.data.frame(c(summary(area39)$coefficients[14,1]))
area_se39i <- as.data.frame(c(summary(area39)$coefficients[14,2]))
area_p39i <- as.data.frame(c(summary(area39)$coefficients[14,5]))
area_b39s <- as.data.frame(c(summary(area39)$coefficients[16,1]))
area_se39s <- as.data.frame(c(summary(area39)$coefficients[16,2]))
area_p39s <- as.data.frame(c(summary(area39)$coefficients[16,5]))

area40 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_fusiformrh_area_cent + scale(subj_fusiformrh_area) + site_fusiformrh_area_cent*wave + scale(subj_fusiformrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b40i <- as.data.frame(c(summary(area40)$coefficients[14,1]))
area_se40i <- as.data.frame(c(summary(area40)$coefficients[14,2]))
area_p40i <- as.data.frame(c(summary(area40)$coefficients[14,5]))
area_b40s <- as.data.frame(c(summary(area40)$coefficients[16,1]))
area_se40s <- as.data.frame(c(summary(area40)$coefficients[16,2]))
area_p40s <- as.data.frame(c(summary(area40)$coefficients[16,5]))

area41 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_ifplrh_area_cent + scale(subj_ifplrh_area) + site_ifplrh_area_cent*wave + scale(subj_ifplrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b41i <- as.data.frame(c(summary(area41)$coefficients[14,1]))
area_se41i <- as.data.frame(c(summary(area41)$coefficients[14,2]))
area_p41i <- as.data.frame(c(summary(area41)$coefficients[14,5]))
area_b41s <- as.data.frame(c(summary(area41)$coefficients[16,1]))
area_se41s <- as.data.frame(c(summary(area41)$coefficients[16,2]))
area_p41s <- as.data.frame(c(summary(area41)$coefficients[16,5]))

area42 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_iftmrh_area_cent + scale(subj_iftmrh_area) + site_iftmrh_area_cent*wave + scale(subj_iftmrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b42i <- as.data.frame(c(summary(area42)$coefficients[14,1]))
area_se42i <- as.data.frame(c(summary(area42)$coefficients[14,2]))
area_p42i <- as.data.frame(c(summary(area42)$coefficients[14,5]))
area_b42s <- as.data.frame(c(summary(area42)$coefficients[16,1]))
area_se42s <- as.data.frame(c(summary(area42)$coefficients[16,2]))
area_p42s <- as.data.frame(c(summary(area42)$coefficients[16,5]))

area43 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_ihcaterh_area_cent + scale(subj_ihcaterh_area) + site_ihcaterh_area_cent*wave + scale(subj_ihcaterh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b43i <- as.data.frame(c(summary(area43)$coefficients[14,1]))
area_se43i <- as.data.frame(c(summary(area43)$coefficients[14,2]))
area_p43i <- as.data.frame(c(summary(area43)$coefficients[14,5]))
area_b43s <- as.data.frame(c(summary(area43)$coefficients[16,1]))
area_se43s <- as.data.frame(c(summary(area43)$coefficients[16,2]))
area_p43s <- as.data.frame(c(summary(area43)$coefficients[16,5]))

area44 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_loccrh_area_cent + scale(subj_loccrh_area) + site_loccrh_area_cent*wave + scale(subj_loccrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b44i <- as.data.frame(c(summary(area44)$coefficients[14,1]))
area_se44i <- as.data.frame(c(summary(area44)$coefficients[14,2]))
area_p44i <- as.data.frame(c(summary(area44)$coefficients[14,5]))
area_b44s <- as.data.frame(c(summary(area44)$coefficients[16,1]))
area_se44s <- as.data.frame(c(summary(area44)$coefficients[16,2]))
area_p44s <- as.data.frame(c(summary(area44)$coefficients[16,5]))

area45 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_lobfrrh_area_cent + scale(subj_lobfrrh_area) + site_lobfrrh_area_cent*wave + scale(subj_lobfrrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b45i <- as.data.frame(c(summary(area45)$coefficients[14,1]))
area_se45i <- as.data.frame(c(summary(area45)$coefficients[14,2]))
area_p45i <- as.data.frame(c(summary(area45)$coefficients[14,5]))
area_b45s <- as.data.frame(c(summary(area45)$coefficients[16,1]))
area_se45s <- as.data.frame(c(summary(area45)$coefficients[16,2]))
area_p45s <- as.data.frame(c(summary(area45)$coefficients[16,5]))

area46 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_lingualrh_area_cent + scale(subj_lingualrh_area) + site_lingualrh_area_cent*wave + scale(subj_lingualrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b46i <- as.data.frame(c(summary(area46)$coefficients[14,1]))
area_se46i <- as.data.frame(c(summary(area46)$coefficients[14,2]))
area_p46i <- as.data.frame(c(summary(area46)$coefficients[14,5]))
area_b46s <- as.data.frame(c(summary(area46)$coefficients[16,1]))
area_se46s <- as.data.frame(c(summary(area46)$coefficients[16,2]))
area_p46s <- as.data.frame(c(summary(area46)$coefficients[16,5]))

area47 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_mobfrrh_area_cent + scale(subj_mobfrrh_area) + site_mobfrrh_area_cent*wave + scale(subj_mobfrrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b47i <- as.data.frame(c(summary(area47)$coefficients[14,1]))
area_se47i <- as.data.frame(c(summary(area47)$coefficients[14,2]))
area_p47i <- as.data.frame(c(summary(area47)$coefficients[14,5]))
area_b47s <- as.data.frame(c(summary(area47)$coefficients[16,1]))
area_se47s <- as.data.frame(c(summary(area47)$coefficients[16,2]))
area_p47s <- as.data.frame(c(summary(area47)$coefficients[16,5]))

area48 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_mdtmrh_area_cent + scale(subj_mdtmrh_area) + site_mdtmrh_area_cent*wave + scale(subj_mdtmrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b48i <- as.data.frame(c(summary(area48)$coefficients[14,1]))
area_se48i <- as.data.frame(c(summary(area48)$coefficients[14,2]))
area_p48i <- as.data.frame(c(summary(area48)$coefficients[14,5]))
area_b48s <- as.data.frame(c(summary(area48)$coefficients[16,1]))
area_se48s <- as.data.frame(c(summary(area48)$coefficients[16,2]))
area_p48s <- as.data.frame(c(summary(area48)$coefficients[16,5]))

area49 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_parahpalrh_area_cent + scale(subj_parahpalrh_area) + site_parahpalrh_area_cent*wave + scale(subj_parahpalrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b49i <- as.data.frame(c(summary(area49)$coefficients[14,1]))
area_se49i <- as.data.frame(c(summary(area49)$coefficients[14,2]))
area_p49i <- as.data.frame(c(summary(area49)$coefficients[14,5]))
area_b49s <- as.data.frame(c(summary(area49)$coefficients[16,1]))
area_se49s <- as.data.frame(c(summary(area49)$coefficients[16,2]))
area_p49s <- as.data.frame(c(summary(area49)$coefficients[16,5]))

area50 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_paracnrh_area_cent + scale(subj_paracnrh_area) + site_paracnrh_area_cent*wave + scale(subj_paracnrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b50i <- as.data.frame(c(summary(area50)$coefficients[14,1]))
area_se50i <- as.data.frame(c(summary(area50)$coefficients[14,2]))
area_p50i <- as.data.frame(c(summary(area50)$coefficients[14,5]))
area_b50s <- as.data.frame(c(summary(area50)$coefficients[16,1]))
area_se50s <- as.data.frame(c(summary(area50)$coefficients[16,2]))
area_p50s <- as.data.frame(c(summary(area50)$coefficients[16,5]))

area51 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_parsopcrh_area_cent + scale(subj_parsopcrh_area) + site_parsopcrh_area_cent*wave + scale(subj_parsopcrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b51i <- as.data.frame(c(summary(area51)$coefficients[14,1]))
area_se51i <- as.data.frame(c(summary(area51)$coefficients[14,2]))
area_p51i <- as.data.frame(c(summary(area51)$coefficients[14,5]))
area_b51s <- as.data.frame(c(summary(area51)$coefficients[16,1]))
area_se51s <- as.data.frame(c(summary(area51)$coefficients[16,2]))
area_p51s <- as.data.frame(c(summary(area51)$coefficients[16,5]))

area52 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_parsobisrh_area_cent + scale(subj_parsobisrh_area) + site_parsobisrh_area_cent*wave + scale(subj_parsobisrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b52i <- as.data.frame(c(summary(area52)$coefficients[14,1]))
area_se52i <- as.data.frame(c(summary(area52)$coefficients[14,2]))
area_p52i <- as.data.frame(c(summary(area52)$coefficients[14,5]))
area_b52s <- as.data.frame(c(summary(area52)$coefficients[16,1]))
area_se52s <- as.data.frame(c(summary(area52)$coefficients[16,2]))
area_p52s <- as.data.frame(c(summary(area52)$coefficients[16,5]))

area53 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_parstgrisrh_area_cent + scale(subj_parstgrisrh_area) + site_parstgrisrh_area_cent*wave + scale(subj_parstgrisrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b53i <- as.data.frame(c(summary(area53)$coefficients[14,1]))
area_se53i <- as.data.frame(c(summary(area53)$coefficients[14,2]))
area_p53i <- as.data.frame(c(summary(area53)$coefficients[14,5]))
area_b53s <- as.data.frame(c(summary(area53)$coefficients[16,1]))
area_se53s <- as.data.frame(c(summary(area53)$coefficients[16,2]))
area_p53s <- as.data.frame(c(summary(area53)$coefficients[16,5]))

area54 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_periccrh_area_cent + scale(subj_periccrh_area) + site_periccrh_area_cent*wave + scale(subj_periccrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b54i <- as.data.frame(c(summary(area54)$coefficients[14,1]))
area_se54i <- as.data.frame(c(summary(area54)$coefficients[14,2]))
area_p54i <- as.data.frame(c(summary(area54)$coefficients[14,5]))
area_b54s <- as.data.frame(c(summary(area54)$coefficients[16,1]))
area_se54s <- as.data.frame(c(summary(area54)$coefficients[16,2]))
area_p54s <- as.data.frame(c(summary(area54)$coefficients[16,5]))

area55 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_postcnrh_area_cent + scale(subj_postcnrh_area) + site_postcnrh_area_cent*wave + scale(subj_postcnrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b55i <- as.data.frame(c(summary(area55)$coefficients[14,1]))
area_se55i <- as.data.frame(c(summary(area55)$coefficients[14,2]))
area_p55i <- as.data.frame(c(summary(area55)$coefficients[14,5]))
area_b55s <- as.data.frame(c(summary(area55)$coefficients[16,1]))
area_se55s <- as.data.frame(c(summary(area55)$coefficients[16,2]))
area_p55s <- as.data.frame(c(summary(area55)$coefficients[16,5]))

area56 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_ptcaterh_area_cent + scale(subj_ptcaterh_area) + site_ptcaterh_area_cent*wave + scale(subj_ptcaterh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b56i <- as.data.frame(c(summary(area56)$coefficients[14,1]))
area_se56i <- as.data.frame(c(summary(area56)$coefficients[14,2]))
area_p56i <- as.data.frame(c(summary(area56)$coefficients[14,5]))
area_b56s <- as.data.frame(c(summary(area56)$coefficients[16,1]))
area_se56s <- as.data.frame(c(summary(area56)$coefficients[16,2]))
area_p56s <- as.data.frame(c(summary(area56)$coefficients[16,5]))

area57 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_precnrh_area_cent + scale(subj_precnrh_area) + site_precnrh_area_cent*wave + scale(subj_precnrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b57i <- as.data.frame(c(summary(area57)$coefficients[14,1]))
area_se57i <- as.data.frame(c(summary(area57)$coefficients[14,2]))
area_p57i <- as.data.frame(c(summary(area57)$coefficients[14,5]))
area_b57s <- as.data.frame(c(summary(area57)$coefficients[16,1]))
area_se57s <- as.data.frame(c(summary(area57)$coefficients[16,2]))
area_p57s <- as.data.frame(c(summary(area57)$coefficients[16,5]))

area58 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_pcrh_area_cent + scale(subj_pcrh_area) + site_pcrh_area_cent*wave + scale(subj_pcrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b58i <- as.data.frame(c(summary(area58)$coefficients[14,1]))
area_se58i <- as.data.frame(c(summary(area58)$coefficients[14,2]))
area_p58i <- as.data.frame(c(summary(area58)$coefficients[14,5]))
area_b58s <- as.data.frame(c(summary(area58)$coefficients[16,1]))
area_se58s <- as.data.frame(c(summary(area58)$coefficients[16,2]))
area_p58s <- as.data.frame(c(summary(area58)$coefficients[16,5]))

area59 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_rracaterh_area_cent + scale(subj_rracaterh_area) + site_rracaterh_area_cent*wave + scale(subj_rracaterh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b59i <- as.data.frame(c(summary(area59)$coefficients[14,1]))
area_se59i <- as.data.frame(c(summary(area59)$coefficients[14,2]))
area_p59i <- as.data.frame(c(summary(area59)$coefficients[14,5]))
area_b59s <- as.data.frame(c(summary(area59)$coefficients[16,1]))
area_se59s <- as.data.frame(c(summary(area59)$coefficients[16,2]))
area_p59s <- as.data.frame(c(summary(area59)$coefficients[16,5]))

area60 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_rrmdfrrh_area_cent + scale(subj_rrmdfrrh_area) + site_rrmdfrrh_area_cent*wave + scale(subj_rrmdfrrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b60i <- as.data.frame(c(summary(area60)$coefficients[14,1]))
area_se60i <- as.data.frame(c(summary(area60)$coefficients[14,2]))
area_p60i <- as.data.frame(c(summary(area60)$coefficients[14,5]))
area_b60s <- as.data.frame(c(summary(area60)$coefficients[16,1]))
area_se60s <- as.data.frame(c(summary(area60)$coefficients[16,2]))
area_p60s <- as.data.frame(c(summary(area60)$coefficients[16,5]))

area61 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_sufrrh_area_cent + scale(subj_sufrrh_area) + site_sufrrh_area_cent*wave + scale(subj_sufrrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b61i <- as.data.frame(c(summary(area61)$coefficients[14,1]))
area_se61i <- as.data.frame(c(summary(area61)$coefficients[14,2]))
area_p61i <- as.data.frame(c(summary(area61)$coefficients[14,5]))
area_b61s <- as.data.frame(c(summary(area61)$coefficients[16,1]))
area_se61s <- as.data.frame(c(summary(area61)$coefficients[16,2]))
area_p61s <- as.data.frame(c(summary(area61)$coefficients[16,5]))

area62 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_suplrh_area_cent + scale(subj_suplrh_area) + site_suplrh_area_cent*wave + scale(subj_suplrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b62i <- as.data.frame(c(summary(area62)$coefficients[14,1]))
area_se62i <- as.data.frame(c(summary(area62)$coefficients[14,2]))
area_p62i <- as.data.frame(c(summary(area62)$coefficients[14,5]))
area_b62s <- as.data.frame(c(summary(area62)$coefficients[16,1]))
area_se62s <- as.data.frame(c(summary(area62)$coefficients[16,2]))
area_p62s <- as.data.frame(c(summary(area62)$coefficients[16,5]))

area63 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_sutmrh_area_cent + scale(subj_sutmrh_area) + site_sutmrh_area_cent*wave + scale(subj_sutmrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b63i <- as.data.frame(c(summary(area63)$coefficients[14,1]))
area_se63i <- as.data.frame(c(summary(area63)$coefficients[14,2]))
area_p63i <- as.data.frame(c(summary(area63)$coefficients[14,5]))
area_b63s <- as.data.frame(c(summary(area63)$coefficients[16,1]))
area_se63s <- as.data.frame(c(summary(area63)$coefficients[16,2]))
area_p63s <- as.data.frame(c(summary(area63)$coefficients[16,5]))

area64 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_smrh_area_cent + scale(subj_smrh_area) + site_smrh_area_cent*wave + scale(subj_smrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b64i <- as.data.frame(c(summary(area64)$coefficients[14,1]))
area_se64i <- as.data.frame(c(summary(area64)$coefficients[14,2]))
area_p64i <- as.data.frame(c(summary(area64)$coefficients[14,5]))
area_b64s <- as.data.frame(c(summary(area64)$coefficients[16,1]))
area_se64s <- as.data.frame(c(summary(area64)$coefficients[16,2]))
area_p64s <- as.data.frame(c(summary(area64)$coefficients[16,5]))

area65 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_frpolerh_area_cent + scale(subj_frpolerh_area) + site_frpolerh_area_cent*wave + scale(subj_frpolerh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b65i <- as.data.frame(c(summary(area65)$coefficients[14,1]))
area_se65i <- as.data.frame(c(summary(area65)$coefficients[14,2]))
area_p65i <- as.data.frame(c(summary(area65)$coefficients[14,5]))
area_b65s <- as.data.frame(c(summary(area65)$coefficients[16,1]))
area_se65s <- as.data.frame(c(summary(area65)$coefficients[16,2]))
area_p65s <- as.data.frame(c(summary(area65)$coefficients[16,5]))

area66 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_tmpolerh_area_cent + scale(subj_tmpolerh_area) + site_tmpolerh_area_cent*wave + scale(subj_tmpolerh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b66i <- as.data.frame(c(summary(area66)$coefficients[14,1]))
area_se66i <- as.data.frame(c(summary(area66)$coefficients[14,2]))
area_p66i <- as.data.frame(c(summary(area66)$coefficients[14,5]))
area_b66s <- as.data.frame(c(summary(area66)$coefficients[16,1]))
area_se66s <- as.data.frame(c(summary(area66)$coefficients[16,2]))
area_p66s <- as.data.frame(c(summary(area66)$coefficients[16,5]))

area67 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_trvtmrh_area_cent + scale(subj_trvtmrh_area) + site_trvtmrh_area_cent*wave + scale(subj_trvtmrh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b67i <- as.data.frame(c(summary(area67)$coefficients[14,1]))
area_se67i <- as.data.frame(c(summary(area67)$coefficients[14,2]))
area_p67i <- as.data.frame(c(summary(area67)$coefficients[14,5]))
area_b67s <- as.data.frame(c(summary(area67)$coefficients[16,1]))
area_se67s <- as.data.frame(c(summary(area67)$coefficients[16,2]))
area_p67s <- as.data.frame(c(summary(area67)$coefficients[16,5]))

area68 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_insularh_area_cent + scale(subj_insularh_area) + site_insularh_area_cent*wave + scale(subj_insularh_area)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
area_b68i <- as.data.frame(c(summary(area68)$coefficients[14,1]))
area_se68i <- as.data.frame(c(summary(area68)$coefficients[14,2]))
area_p68i <- as.data.frame(c(summary(area68)$coefficients[14,5]))
area_b68s <- as.data.frame(c(summary(area68)$coefficients[16,1]))
area_se68s <- as.data.frame(c(summary(area68)$coefficients[16,2]))
area_p68s <- as.data.frame(c(summary(area68)$coefficients[16,5]))

#create data frame with all parcel-wise cortical areaume st. estimates, SEs, and p-values
area_parcel <- data.frame(x=c("area"))

newareabi <- c(area_b1i,area_b2i,area_b3i,area_b4i,area_b5i,area_b6i,area_b7i,area_b8i,area_b9i,area_b10i,area_b11i,area_b12i,area_b13i,area_b14i,area_b15i,
               area_b16i,area_b17i,area_b18i,area_b19i,area_b20i,area_b21i,area_b22i,area_b23i,area_b24i,area_b25i,area_b26i,area_b27i,area_b28i,area_b29i,
               area_b30i,area_b31i,area_b32i,area_b33i,area_b34i,area_b35i,area_b36i,area_b37i,area_b38i,area_b39i,area_b40i,area_b41i,area_b42i,area_b43i,
               area_b44i,area_b45i,area_b46i,area_b47i,area_b48i,area_b49i,area_b50i,area_b51i,area_b52i,area_b53i,area_b54i,area_b55i,area_b56i,area_b57i,
               area_b58i,area_b59i,area_b60i,area_b61i,area_b62i,area_b63i,area_b64i,area_b65i,area_b66i,area_b67i,area_b68i)
area_parcel <- cbind(area_parcel,newareabi)

newareabs <- c(area_b1s,area_b2s,area_b3s,area_b4s,area_b5s,area_b6s,area_b7s,area_b8s,area_b9s,area_b10s,area_b11s,area_b12s,area_b13s,area_b14s,area_b15s,
               area_b16s,area_b17s,area_b18s,area_b19s,area_b20s,area_b21s,area_b22s,area_b23s,area_b24s,area_b25s,area_b26s,area_b27s,area_b28s,area_b29s,
               area_b30s,area_b31s,area_b32s,area_b33s,area_b34s,area_b35s,area_b36s,area_b37s,area_b38s,area_b39s,area_b40s,area_b41s,area_b42s,area_b43s,
               area_b44s,area_b45s,area_b46s,area_b47s,area_b48s,area_b49s,area_b50s,area_b51s,area_b52s,area_b53s,area_b54s,area_b55s,area_b56s,area_b57s,
               area_b58s,area_b59s,area_b60s,area_b61s,area_b62s,area_b63s,area_b64s,area_b65s,area_b66s,area_b67s,area_b68s)
area_parcel <- cbind(area_parcel,newareabs)

newareasei <- c(area_se1i,area_se2i,area_se3i,area_se4i,area_se5i,area_se6i,area_se7i,area_se8i,area_se9i,area_se10i,area_se11i,area_se12i,area_se13i,
                area_se14i,area_se15i,area_se16i,area_se17i,area_se18i,area_se19i,area_se20i,area_se21i,area_se22i,area_se23i,area_se24i,area_se25i,
                area_se26i,area_se27i,area_se28i,area_se29i,area_se30i,area_se31i,area_se32i,area_se33i,area_se34i,area_se35i,area_se36i,area_se37i,
                area_se38i,area_se39i,area_se40i,area_se41i,area_se42i,area_se43i,area_se44i,area_se45i,area_se46i,area_se47i,area_se48i,area_se49i,
                area_se50i,area_se51i,area_se52i,area_se53i,area_se54i,area_se55i,area_se56i,area_se57i,area_se58i,area_se59i,area_se60i,area_se61i,
                area_se62i,area_se63i,area_se64i,area_se65i,area_se66i,area_se67i,area_se68i)
area_parcel <- cbind(area_parcel,newareasei)

newareases <- c(area_se1s,area_se2s,area_se3s,area_se4s,area_se5s,area_se6s,area_se7s,area_se8s,area_se9s,area_se10s,area_se11s,area_se12s,area_se13s,
                area_se14s,area_se15s,area_se16s,area_se17s,area_se18s,area_se19s,area_se20s,area_se21s,area_se22s,area_se23s,area_se24s,area_se25s,
                area_se26s,area_se27s,area_se28s,area_se29s,area_se30s,area_se31s,area_se32s,area_se33s,area_se34s,area_se35s,area_se36s,area_se37s,
                area_se38s,area_se39s,area_se40s,area_se41s,area_se42s,area_se43s,area_se44s,area_se45s,area_se46s,area_se47s,area_se48s,area_se49s,
                area_se50s,area_se51s,area_se52s,area_se53s,area_se54s,area_se55s,area_se56s,area_se57s,area_se58s,area_se59s,area_se60s,area_se61s,
                area_se62s,area_se63s,area_se64s,area_se65s,area_se66s,area_se67s,area_se68s)
area_parcel <- cbind(area_parcel,newareases)

newareapi <- c(area_p1i,area_p2i,area_p3i,area_p4i,area_p5i,area_p6i,area_p7i,area_p8i,area_p9i,area_p10i,area_p11i,area_p12i,area_p13i,area_p14i,area_p15i,
               area_p16i,area_p17i,area_p18i,area_p19i,area_p20i,area_p21i,area_p22i,area_p23i,area_p24i,area_p25i,area_p26i,area_p27i,area_p28i,area_p29i,
               area_p30i,area_p31i,area_p32i,area_p33i,area_p34i,area_p35i,area_p36i,area_p37i,area_p38i,area_p39i,area_p40i,area_p41i,area_p42i,area_p43i,
               area_p44i,area_p45i,area_p46i,area_p47i,area_p48i,area_p49i,area_p50i,area_p51i,area_p52i,area_p53i,area_p54i,area_p55i,area_p56i,area_p57i,
               area_p58i,area_p59i,area_p60i,area_p61i,area_p62i,area_p63i,area_p64i,area_p65i,area_p66i,area_p67i,area_p68i)
area_parcel <- cbind(area_parcel,newareapi)

newareaps <- c(area_p1s,area_p2s,area_p3s,area_p4s,area_p5s,area_p6s,area_p7s,area_p8s,area_p9s,area_p10s,area_p11s,area_p12s,area_p13s,area_p14s,area_p15s,
               area_p16s,area_p17s,area_p18s,area_p19s,area_p20s,area_p21s,area_p22s,area_p23s,area_p24s,area_p25s,area_p26s,area_p27s,area_p28s,area_p29s,
               area_p30s,area_p31s,area_p32s,area_p33s,area_p34s,area_p35s,area_p36s,area_p37s,area_p38s,area_p39s,area_p40s,area_p41s,area_p42s,area_p43s,
               area_p44s,area_p45s,area_p46s,area_p47s,area_p48s,area_p49s,area_p50s,area_p51s,area_p52s,area_p53s,area_p54s,area_p55s,area_p56s,area_p57s,
               area_p58s,area_p59s,area_p60s,area_p61s,area_p62s,area_p63s,area_p64s,area_p65s,area_p66s,area_p67s,area_p68s)
area_parcel <- cbind(area_parcel,newareaps)

names(area_parcel) <- c('area','banksstslh_area_bi','cdacatelh_area_bi','cdmdfrlh_area_bi','cuneuslh_area_bi','ehinallh_area_bi','fusiformlh_area_bi',
                        'ifpllh_area_bi','iftmlh_area_bi','ihcatelh_area_bi','locclh_area_bi','lobfrlh_area_bi','linguallh_area_bi',
                        'mobfrlh_area_bi','mdtmlh_area_bi','parahpallh_area_bi','paracnlh_area_bi','parsopclh_area_bi','parsobislh_area_bi',
                        'parstgrislh_area_bi','pericclh_area_bi','postcnlh_area_bi','ptcatelh_area_bi','precnlh_area_bi','pclh_area_bi',
                        'rracatelh_area_bi','rrmdfrlh_area_bi','sufrlh_area_bi','supllh_area_bi','sutmlh_area_bi','smlh_area_bi','frpolelh_area_bi',
                        'tmpolelh_area_bi','trvtmlh_area_bi','insulalh_area_bi','banksstsrh_area_bi','cdacaterh_area_bi','cdmdfrrh_area_bi',
                        'cuneusrh_area_bi','ehinalrh_area_bi','fusiformrh_area_bi','ifplrh_area_bi','iftmrh_area_bi','ihcaterh_area_bi',
                        'loccrh_area_bi','lobfrrh_area_bi','lingualrh_area_bi','mobfrrh_area_bi','mdtmrh_area_bi','parahpalrh_area_bi',
                        'paracnrh_area_bi','parsopcrh_area_bi','parsobisrh_area_bi','parstgrisrh_area_bi','periccrh_area_bi',
                        'postcnrh_area_bi','ptcaterh_area_bi','precnrh_area_bi','pcrh_area_bi','rracaterh_area_bi','rrmdfrrh_area_bi',
                        'sufrrh_area_bi','suplrh_area_bi','sutmrh_area_bi','smrh_area_bi','frpolerh_area_bi',
                        'tmpolerh_area_bi','trvtmrh_area_bi','insularh_area_bi',
                        'banksstslh_area_bs','cdacatelh_area_bs','cdmdfrlh_area_bs','cuneuslh_area_bs','ehinallh_area_bs','fusiformlh_area_bs',
                        'ifpllh_area_bs','iftmlh_area_bs','ihcatelh_area_bs','locclh_area_bs','lobfrlh_area_bs','linguallh_area_bs',
                        'mobfrlh_area_bs','mdtmlh_area_bs','parahpallh_area_bs','paracnlh_area_bs','parsopclh_area_bs','parsobislh_area_bs',
                        'parstgrislh_area_bs','pericclh_area_bs','postcnlh_area_bs','ptcatelh_area_bs','precnlh_area_bs','pclh_area_bs',
                        'rracatelh_area_bs','rrmdfrlh_area_bs','sufrlh_area_bs','supllh_area_bs','sutmlh_area_bs','smlh_area_bs','frpolelh_area_bs',
                        'tmpolelh_area_bs','trvtmlh_area_bs','insulalh_area_bs','banksstsrh_area_bs','cdacaterh_area_bs','cdmdfrrh_area_bs',
                        'cuneusrh_area_bs','ehinalrh_area_bs','fusiformrh_area_bs','ifplrh_area_bs','iftmrh_area_bs','ihcaterh_area_bs',
                        'loccrh_area_bs','lobfrrh_area_bs','lingualrh_area_bs','mobfrrh_area_bs','mdtmrh_area_bs','parahpalrh_area_bs',
                        'paracnrh_area_bs','parsopcrh_area_bs','parsobisrh_area_bs','parstgrisrh_area_bs','periccrh_area_bs',
                        'postcnrh_area_bs','ptcaterh_area_bs','precnrh_area_bs','pcrh_area_bs','rracaterh_area_bs','rrmdfrrh_area_bs',
                        'sufrrh_area_bs','suplrh_area_bs','sutmrh_area_bs','smrh_area_bs','frpolerh_area_bs',
                        'tmpolerh_area_bs','trvtmrh_area_bs','insularh_area_bs',
                        'banksstslh_area_sei','cdacatelh_area_sei','cdmdfrlh_area_sei','cuneuslh_area_sei','ehinallh_area_sei','fusiformlh_area_sei',
                        'ifpllh_area_sei','iftmlh_area_sei','ihcatelh_area_sei','locclh_area_sei','lobfrlh_area_sei','linguallh_area_sei',
                        'mobfrlh_area_sei','mdtmlh_area_sei','parahpallh_area_sei','paracnlh_area_sei','parsopclh_area_sei',
                        'parsobislh_area_sei','parstgrislh_area_sei','pericclh_area_sei','postcnlh_area_sei','ptcatelh_area_sei',
                        'precnlh_area_sei','pclh_area_sei','rracatelh_area_sei','rrmdfrlh_area_sei','sufrlh_area_sei','supllh_area_sei',
                        'sutmlh_area_sei','smlh_area_sei','frpolelh_area_sei','tmpolelh_area_sei','trvtmlh_area_sei','insulalh_area_sei',
                        'banksstsrh_area_sei','cdacaterh_area_sei','cdmdfrrh_area_sei','cuneusrh_area_sei','ehinalrh_area_sei',
                        'fusiformrh_area_sei','ifplrh_area_sei','iftmrh_area_sei','ihcaterh_area_sei','loccrh_area_sei','lobfrrh_area_sei',
                        'lingualrh_area_sei','mobfrrh_area_sei','mdtmrh_area_sei','parahpalrh_area_sei','paracnrh_area_sei',
                        'parsopcrh_area_sei','parsobisrh_area_sei','parstgrisrh_area_sei','periccrh_area_sei','postcnrh_area_sei',
                        'ptcaterh_area_sei','precnrh_area_sei','pcrh_area_sei','rracaterh_area_sei','rrmdfrrh_area_sei',
                        'sufrrh_area_sei','suplrh_area_sei','sutmrh_area_sei','smrh_area_sei','frpolerh_area_sei',
                        'tmpolerh_area_sei','trvtmrh_area_sei','insularh_area_sei',
                        'banksstslh_area_ses','cdacatelh_area_ses','cdmdfrlh_area_ses','cuneuslh_area_ses','ehinallh_area_ses','fusiformlh_area_ses',
                        'ifpllh_area_ses','iftmlh_area_ses','ihcatelh_area_ses','locclh_area_ses','lobfrlh_area_ses','linguallh_area_ses',
                        'mobfrlh_area_ses','mdtmlh_area_ses','parahpallh_area_ses','paracnlh_area_ses','parsopclh_area_ses',
                        'parsobislh_area_ses','parstgrislh_area_ses','pericclh_area_ses','postcnlh_area_ses','ptcatelh_area_ses',
                        'precnlh_area_ses','pclh_area_ses','rracatelh_area_ses','rrmdfrlh_area_ses','sufrlh_area_ses','supllh_area_ses',
                        'sutmlh_area_ses','smlh_area_ses','frpolelh_area_ses','tmpolelh_area_ses','trvtmlh_area_ses','insulalh_area_ses',
                        'banksstsrh_area_ses','cdacaterh_area_ses','cdmdfrrh_area_ses','cuneusrh_area_ses','ehinalrh_area_ses',
                        'fusiformrh_area_ses','ifplrh_area_ses','iftmrh_area_ses','ihcaterh_area_ses','loccrh_area_ses','lobfrrh_area_ses',
                        'lingualrh_area_ses','mobfrrh_area_ses','mdtmrh_area_ses','parahpalrh_area_ses','paracnrh_area_ses',
                        'parsopcrh_area_ses','parsobisrh_area_ses','parstgrisrh_area_ses','periccrh_area_ses','postcnrh_area_ses',
                        'ptcaterh_area_ses','precnrh_area_ses','pcrh_area_ses','rracaterh_area_ses','rrmdfrrh_area_ses',
                        'sufrrh_area_ses','suplrh_area_ses','sutmrh_area_ses','smrh_area_ses','frpolerh_area_ses',
                        'tmpolerh_area_ses','trvtmrh_area_ses','insularh_area_ses',
                        'banksstslh_area_pi','cdacatelh_area_pi','cdmdfrlh_area_pi','cuneuslh_area_pi','ehinallh_area_pi','fusiformlh_area_pi','ifpllh_area_pi','iftmlh_area_pi',
                        'ihcatelh_area_pi','locclh_area_pi','lobfrlh_area_pi','linguallh_area_pi','mobfrlh_area_pi','mdtmlh_area_pi',
                        'parahpallh_area_pi','paracnlh_area_pi','parsopclh_area_pi','parsobislh_area_pi','parstgrislh_area_pi',
                        'pericclh_area_pi','postcnlh_area_pi','ptcatelh_area_pi','precnlh_area_pi','pclh_area_pi','rracatelh_area_pi',
                        'rrmdfrlh_area_pi','sufrlh_area_pi','supllh_area_pi','sutmlh_area_pi','smlh_area_pi','frpolelh_area_pi',
                        'tmpolelh_area_pi','trvtmlh_area_pi','insulalh_area_pi','banksstsrh_area_pi','cdacaterh_area_pi','cdmdfrrh_area_pi',
                        'cuneusrh_area_pi','ehinalrh_area_pi','fusiformrh_area_pi','ifplrh_area_pi','iftmrh_area_pi','ihcaterh_area_pi',
                        'loccrh_area_pi','lobfrrh_area_pi','lingualrh_area_pi','mobfrrh_area_pi','mdtmrh_area_pi','parahpalrh_area_pi',
                        'paracnrh_area_pi','parsopcrh_area_pi','parsobisrh_area_pi','parstgrisrh_area_pi','periccrh_area_pi',
                        'postcnrh_area_pi','ptcaterh_area_pi','precnrh_area_pi','pcrh_area_pi','rracaterh_area_pi','rrmdfrrh_area_pi',
                        'sufrrh_area_pi','suplrh_area_pi','sutmrh_area_pi','smrh_area_pi','frpolerh_area_pi',
                        'tmpolerh_area_pi','trvtmrh_area_pi','insularh_area_pi',
                        'banksstslh_area_ps','cdacatelh_area_ps','cdmdfrlh_area_ps','cuneuslh_area_ps','ehinallh_area_ps','fusiformlh_area_ps','ifpllh_area_ps','iftmlh_area_ps',
                        'ihcatelh_area_ps','locclh_area_ps','lobfrlh_area_ps','linguallh_area_ps','mobfrlh_area_ps','mdtmlh_area_ps',
                        'parahpallh_area_ps','paracnlh_area_ps','parsopclh_area_ps','parsobislh_area_ps','parstgrislh_area_ps',
                        'pericclh_area_ps','postcnlh_area_ps','ptcatelh_area_ps','precnlh_area_ps','pclh_area_ps','rracatelh_area_ps',
                        'rrmdfrlh_area_ps','sufrlh_area_ps','supllh_area_ps','sutmlh_area_ps','smlh_area_ps','frpolelh_area_ps',
                        'tmpolelh_area_ps','trvtmlh_area_ps','insulalh_area_ps','banksstsrh_area_ps','cdacaterh_area_ps','cdmdfrrh_area_ps',
                        'cuneusrh_area_ps','ehinalrh_area_ps','fusiformrh_area_ps','ifplrh_area_ps','iftmrh_area_ps','ihcaterh_area_ps',
                        'loccrh_area_ps','lobfrrh_area_ps','lingualrh_area_ps','mobfrrh_area_ps','mdtmrh_area_ps','parahpalrh_area_ps',
                        'paracnrh_area_ps','parsopcrh_area_ps','parsobisrh_area_ps','parstgrisrh_area_ps','periccrh_area_ps',
                        'postcnrh_area_ps','ptcaterh_area_ps','precnrh_area_ps','pcrh_area_ps','rracaterh_area_ps','rrmdfrrh_area_ps',
                        'sufrrh_area_ps','suplrh_area_ps','sutmrh_area_ps','smrh_area_ps','frpolerh_area_ps',
                        'tmpolerh_area_ps','trvtmrh_area_ps','insularh_area_ps')

#calculate 95% CIs and create lower and upper bound variables 
area_parcel$ci_lower_banksstslh_areai <- area_parcel$banksstslh_area_bi - 1.96*area_parcel$banksstslh_area_sei 
area_parcel$ci_upper_banksstslh_areai <- area_parcel$banksstslh_area_bi + 1.96*area_parcel$banksstslh_area_sei
area_parcel$ci_lower_cdacatelh_areai <- area_parcel$cdacatelh_area_bi - 1.96*area_parcel$cdacatelh_area_sei 
area_parcel$ci_upper_cdacatelh_areai <- area_parcel$cdacatelh_area_bi + 1.96*area_parcel$cdacatelh_area_sei
area_parcel$ci_lower_cdmdfrlh_areai <- area_parcel$cdmdfrlh_area_bi - 1.96*area_parcel$cdmdfrlh_area_sei 
area_parcel$ci_upper_cdmdfrlh_areai <- area_parcel$cdmdfrlh_area_bi + 1.96*area_parcel$cdmdfrlh_area_sei
area_parcel$ci_lower_cuneuslh_areai <- area_parcel$cuneuslh_area_bi - 1.96*area_parcel$cuneuslh_area_sei 
area_parcel$ci_upper_cuneuslh_areai <- area_parcel$cuneuslh_area_bi + 1.96*area_parcel$cuneuslh_area_sei
area_parcel$ci_lower_ehinallh_areai <- area_parcel$ehinallh_area_bi - 1.96*area_parcel$ehinallh_area_sei 
area_parcel$ci_upper_ehinallh_areali <- area_parcel$ehinallh_area_bi + 1.96*area_parcel$ehinallh_area_sei
area_parcel$ci_lower_fusiformlh_areai <- area_parcel$fusiformlh_area_bi - 1.96*area_parcel$fusiformlh_area_sei 
area_parcel$ci_upper_fusiformlh_areai <- area_parcel$fusiformlh_area_bi + 1.96*area_parcel$fusiformlh_area_sei
area_parcel$ci_lower_ifpllh_areai <- area_parcel$ifpllh_area_bi - 1.96*area_parcel$ifpllh_area_sei 
area_parcel$ci_upper_ifpllh_areai <- area_parcel$ifpllh_area_bi + 1.96*area_parcel$ifpllh_area_sei
area_parcel$ci_lower_iftmlh_areai <- area_parcel$iftmlh_area_bi - 1.96*area_parcel$iftmlh_area_sei
area_parcel$ci_upper_iftmlh_areai <- area_parcel$iftmlh_area_bi + 1.96*area_parcel$iftmlh_area_sei
area_parcel$ci_lower_ihcatelh_areai <- area_parcel$ihcatelh_area_bi - 1.96*area_parcel$ihcatelh_area_sei 
area_parcel$ci_upper_ihcatelh_areai <- area_parcel$ihcatelh_area_bi + 1.96*area_parcel$ihcatelh_area_sei
area_parcel$ci_lower_locclh_areai <- area_parcel$locclh_area_bi - 1.96*area_parcel$locclh_area_sei
area_parcel$ci_upper_locclh_areai <- area_parcel$locclh_area_bi + 1.96*area_parcel$locclh_area_sei
area_parcel$ci_lower_lobfrlh_areai <- area_parcel$lobfrlh_area_bi - 1.96*area_parcel$lobfrlh_area_sei 
area_parcel$ci_upper_lobfrlh_areai <- area_parcel$lobfrlh_area_bi + 1.96*area_parcel$lobfrlh_area_sei
area_parcel$ci_lower_linguallh_areai <- area_parcel$linguallh_area_bi - 1.96*area_parcel$linguallh_area_sei 
area_parcel$ci_upper_linguallh_areai <- area_parcel$linguallh_area_bi + 1.96*area_parcel$linguallh_area_sei
area_parcel$ci_lower_mobfrlh_areai <- area_parcel$mobfrlh_area_bi - 1.96*area_parcel$mobfrlh_area_sei 
area_parcel$ci_upper_mobfrlh_areai <- area_parcel$mobfrlh_area_bi + 1.96*area_parcel$mobfrlh_area_sei
area_parcel$ci_lower_mdtmlh_areai <- area_parcel$mdtmlh_area_bi - 1.96*area_parcel$mdtmlh_area_sei 
area_parcel$ci_upper_mdtmlh_areai <- area_parcel$mdtmlh_area_bi + 1.96*area_parcel$mdtmlh_area_sei
area_parcel$ci_lower_parahpallh_areai <- area_parcel$parahpallh_area_bi - 1.96*area_parcel$parahpallh_area_sei 
area_parcel$ci_upper_parahpallh_areai <- area_parcel$parahpallh_area_bi + 1.96*area_parcel$parahpallh_area_sei
area_parcel$ci_lower_paracnlh_areai <- area_parcel$paracnlh_area_bi - 1.96*area_parcel$paracnlh_area_sei 
area_parcel$ci_upper_paracnlh_areai <- area_parcel$paracnlh_area_bi + 1.96*area_parcel$paracnlh_area_sei
area_parcel$ci_lower_parsopclh_areai <- area_parcel$parsopclh_area_bi - 1.96*area_parcel$parsopclh_area_sei 
area_parcel$ci_upper_parsopclh_areai <- area_parcel$parsopclh_area_bi + 1.96*area_parcel$parsopclh_area_sei
area_parcel$ci_lower_parsobislh_areai <- area_parcel$parsobislh_area_bi - 1.96*area_parcel$parsobislh_area_sei 
area_parcel$ci_upper_parsobislh_areai <- area_parcel$parsobislh_area_bi + 1.96*area_parcel$parsobislh_area_sei
area_parcel$ci_lower_parstgrislh_areai <- area_parcel$parstgrislh_area_bi - 1.96*area_parcel$parstgrislh_area_sei 
area_parcel$ci_upper_parstgrislh_areai <- area_parcel$parstgrislh_area_bi + 1.96*area_parcel$parstgrislh_area_sei
area_parcel$ci_lower_pericclh_areai <- area_parcel$pericclh_area_bi - 1.96*area_parcel$pericclh_area_sei 
area_parcel$ci_upper_pericclh_areai <- area_parcel$pericclh_area_bi + 1.96*area_parcel$pericclh_area_sei
area_parcel$ci_lower_postcnlh_areai <- area_parcel$postcnlh_area_bi - 1.96*area_parcel$postcnlh_area_sei 
area_parcel$ci_upper_postcnlh_areai <- area_parcel$postcnlh_area_bi + 1.96*area_parcel$postcnlh_area_sei
area_parcel$ci_lower_ptcatelh_areai <- area_parcel$ptcatelh_area_bi - 1.96*area_parcel$ptcatelh_area_sei 
area_parcel$ci_upper_ptcatelh_areai <- area_parcel$ptcatelh_area_bi + 1.96*area_parcel$ptcatelh_area_sei
area_parcel$ci_lower_precnlh_areai <- area_parcel$precnlh_area_bi - 1.96*area_parcel$precnlh_area_sei 
area_parcel$ci_upper_precnlh_areai <- area_parcel$precnlh_area_bi + 1.96*area_parcel$precnlh_area_sei
area_parcel$ci_lower_pclh_areai <- area_parcel$pclh_area_bi - 1.96*area_parcel$pclh_area_sei 
area_parcel$ci_upper_pclh_areai <- area_parcel$pclh_area_bi + 1.96*area_parcel$pclh_area_sei
area_parcel$ci_lower_rracatelh_areai <- area_parcel$rracatelh_area_bi - 1.96*area_parcel$rracatelh_area_sei 
area_parcel$ci_upper_rracatelh_areai <- area_parcel$rracatelh_area_bi + 1.96*area_parcel$rracatelh_area_sei
area_parcel$ci_lower_rrmdfrlh_areai <- area_parcel$rrmdfrlh_area_bi - 1.96*area_parcel$rrmdfrlh_area_sei 
area_parcel$ci_upper_rrmdfrlh_areai <- area_parcel$rrmdfrlh_area_bi + 1.96*area_parcel$rrmdfrlh_area_sei
area_parcel$ci_lower_sufrlh_areai <- area_parcel$sufrlh_area_bi - 1.96*area_parcel$sufrlh_area_sei 
area_parcel$ci_upper_sufrlh_areai <- area_parcel$sufrlh_area_bi + 1.96*area_parcel$sufrlh_area_sei
area_parcel$ci_lower_supllh_areai <- area_parcel$supllh_area_bi - 1.96*area_parcel$supllh_area_sei 
area_parcel$ci_upper_supllh_areai <- area_parcel$supllh_area_bi + 1.96*area_parcel$supllh_area_sei
area_parcel$ci_lower_sutmlh_areai <- area_parcel$sutmlh_area_bi - 1.96*area_parcel$sutmlh_area_sei 
area_parcel$ci_upper_sutmlh_areai <- area_parcel$sutmlh_area_bi + 1.96*area_parcel$sutmlh_area_sei
area_parcel$ci_lower_smlh_areai <- area_parcel$smlh_area_bi - 1.96*area_parcel$smlh_area_sei 
area_parcel$ci_upper_smlh_areai <- area_parcel$smlh_area_bi + 1.96*area_parcel$smlh_area_sei
area_parcel$ci_lower_frpolelh_areai <- area_parcel$frpolelh_area_bi - 1.96*area_parcel$frpolelh_area_sei 
area_parcel$ci_upper_frpolelh_areai <- area_parcel$frpolelh_area_bi + 1.96*area_parcel$frpolelh_area_sei
area_parcel$ci_lower_tmpolelh_areai <- area_parcel$tmpolelh_area_bi - 1.96*area_parcel$tmpolelh_area_sei 
area_parcel$ci_upper_tmpolelh_areai <- area_parcel$tmpolelh_area_bi + 1.96*area_parcel$tmpolelh_area_sei
area_parcel$ci_lower_trvtmlh_areai <- area_parcel$trvtmlh_area_bi - 1.96*area_parcel$trvtmlh_area_sei 
area_parcel$ci_upper_trvtmlh_areai <- area_parcel$trvtmlh_area_bi + 1.96*area_parcel$trvtmlh_area_sei
area_parcel$ci_lower_insulalh_areai <- area_parcel$insulalh_area_bi - 1.96*area_parcel$insulalh_area_sei 
area_parcel$ci_upper_insulalh_areai <- area_parcel$insulalh_area_bi + 1.96*area_parcel$insulalh_area_sei

area_parcel$ci_lower_banksstsrh_areai <- area_parcel$banksstsrh_area_bi - 1.96*area_parcel$banksstsrh_area_sei 
area_parcel$ci_upper_banksstsrh_areai <- area_parcel$banksstsrh_area_bi + 1.96*area_parcel$banksstsrh_area_sei
area_parcel$ci_lower_cdacaterh_areai <- area_parcel$cdacaterh_area_bi - 1.96*area_parcel$cdacaterh_area_sei 
area_parcel$ci_upper_cdacaterh_areai <- area_parcel$cdacaterh_area_bi + 1.96*area_parcel$cdacaterh_area_sei
area_parcel$ci_lower_cdmdfrrh_areai <- area_parcel$cdmdfrrh_area_bi - 1.96*area_parcel$cdmdfrrh_area_sei 
area_parcel$ci_upper_cdmdfrrh_areai <- area_parcel$cdmdfrrh_area_bi + 1.96*area_parcel$cdmdfrrh_area_sei
area_parcel$ci_lower_cuneusrh_areai <- area_parcel$cuneusrh_area_bi - 1.96*area_parcel$cuneusrh_area_sei 
area_parcel$ci_upper_cuneusrh_areai <- area_parcel$cuneusrh_area_bi + 1.96*area_parcel$cuneusrh_area_sei
area_parcel$ci_lower_ehinalrh_areai <- area_parcel$ehinalrh_area_bi - 1.96*area_parcel$ehinalrh_area_sei 
area_parcel$ci_upper_ehinalrh_areali <- area_parcel$ehinalrh_area_bi + 1.96*area_parcel$ehinalrh_area_sei
area_parcel$ci_lower_fusiformrh_areai <- area_parcel$fusiformrh_area_bi - 1.96*area_parcel$fusiformrh_area_sei 
area_parcel$ci_upper_fusiformrh_areai <- area_parcel$fusiformrh_area_bi + 1.96*area_parcel$fusiformrh_area_sei
area_parcel$ci_lower_ifplrh_areai <- area_parcel$ifplrh_area_bi - 1.96*area_parcel$ifplrh_area_sei 
area_parcel$ci_upper_ifplrh_areai <- area_parcel$ifplrh_area_bi + 1.96*area_parcel$ifplrh_area_sei
area_parcel$ci_lower_iftmrh_areai <- area_parcel$iftmrh_area_bi - 1.96*area_parcel$iftmrh_area_sei 
area_parcel$ci_upper_iftmrh_areai <- area_parcel$iftmrh_area_bi + 1.96*area_parcel$iftmrh_area_sei
area_parcel$ci_lower_ihcaterh_areai <- area_parcel$ihcaterh_area_bi - 1.96*area_parcel$ihcaterh_area_sei 
area_parcel$ci_upper_ihcaterh_areai <- area_parcel$ihcaterh_area_bi + 1.96*area_parcel$ihcaterh_area_sei
area_parcel$ci_lower_loccrh_areai <- area_parcel$loccrh_area_bi - 1.96*area_parcel$loccrh_area_sei 
area_parcel$ci_upper_loccrh_areai <- area_parcel$loccrh_area_bi + 1.96*area_parcel$loccrh_area_sei
area_parcel$ci_lower_lobfrrh_areai <- area_parcel$lobfrrh_area_bi - 1.96*area_parcel$lobfrrh_area_sei 
area_parcel$ci_upper_lobfrrh_areai <- area_parcel$lobfrrh_area_bi + 1.96*area_parcel$lobfrrh_area_sei
area_parcel$ci_lower_lingualrh_areai <- area_parcel$lingualrh_area_bi - 1.96*area_parcel$lingualrh_area_sei 
area_parcel$ci_upper_lingualrh_areai <- area_parcel$lingualrh_area_bi + 1.96*area_parcel$lingualrh_area_sei
area_parcel$ci_lower_mobfrrh_areai <- area_parcel$mobfrrh_area_bi - 1.96*area_parcel$mobfrrh_area_sei 
area_parcel$ci_upper_mobfrrh_areai <- area_parcel$mobfrrh_area_bi + 1.96*area_parcel$mobfrrh_area_sei
area_parcel$ci_lower_mdtmrh_areai <- area_parcel$mdtmrh_area_bi - 1.96*area_parcel$mdtmrh_area_sei 
area_parcel$ci_upper_mdtmrh_areai <- area_parcel$mdtmrh_area_bi + 1.96*area_parcel$mdtmrh_area_sei
area_parcel$ci_lower_parahpalrh_areai <- area_parcel$parahpalrh_area_bi - 1.96*area_parcel$parahpalrh_area_sei 
area_parcel$ci_upper_parahpalrh_areai <- area_parcel$parahpalrh_area_bi + 1.96*area_parcel$parahpalrh_area_sei
area_parcel$ci_lower_paracnrh_areai <- area_parcel$paracnrh_area_bi - 1.96*area_parcel$paracnrh_area_sei 
area_parcel$ci_upper_paracnrh_areai <- area_parcel$paracnrh_area_bi + 1.96*area_parcel$paracnrh_area_sei
area_parcel$ci_lower_parsopcrh_areai <- area_parcel$parsopcrh_area_bi - 1.96*area_parcel$parsopcrh_area_sei 
area_parcel$ci_upper_parsopcrh_areai <- area_parcel$parsopcrh_area_bi + 1.96*area_parcel$parsopcrh_area_sei
area_parcel$ci_lower_parsobisrh_areai <- area_parcel$parsobisrh_area_bi - 1.96*area_parcel$parsobisrh_area_sei 
area_parcel$ci_upper_parsobisrh_areai <- area_parcel$parsobisrh_area_bi + 1.96*area_parcel$parsobisrh_area_sei
area_parcel$ci_lower_parstgrisrh_areai <- area_parcel$parstgrisrh_area_bi - 1.96*area_parcel$parstgrisrh_area_sei 
area_parcel$ci_upper_parstgrisrh_areai <- area_parcel$parstgrisrh_area_bi + 1.96*area_parcel$parstgrisrh_area_sei
area_parcel$ci_lower_periccrh_areai <- area_parcel$periccrh_area_bi - 1.96*area_parcel$periccrh_area_sei 
area_parcel$ci_upper_periccrh_areai <- area_parcel$periccrh_area_bi + 1.96*area_parcel$periccrh_area_sei
area_parcel$ci_lower_postcnrh_areai <- area_parcel$postcnrh_area_bi - 1.96*area_parcel$postcnrh_area_sei 
area_parcel$ci_upper_postcnrh_areai <- area_parcel$postcnrh_area_bi + 1.96*area_parcel$postcnrh_area_sei
area_parcel$ci_lower_ptcaterh_areai <- area_parcel$ptcaterh_area_bi - 1.96*area_parcel$ptcaterh_area_sei 
area_parcel$ci_upper_ptcaterh_areai <- area_parcel$ptcaterh_area_bi + 1.96*area_parcel$ptcaterh_area_sei
area_parcel$ci_lower_precnrh_areai <- area_parcel$precnrh_area_bi - 1.96*area_parcel$precnrh_area_sei 
area_parcel$ci_upper_precnrh_areai <- area_parcel$precnrh_area_bi + 1.96*area_parcel$precnrh_area_sei
area_parcel$ci_lower_pcrh_areai <- area_parcel$pcrh_area_bi - 1.96*area_parcel$pcrh_area_sei 
area_parcel$ci_upper_pcrh_areai <- area_parcel$pcrh_area_bi + 1.96*area_parcel$pcrh_area_sei
area_parcel$ci_lower_rracaterh_areai <- area_parcel$rracaterh_area_bi - 1.96*area_parcel$rracaterh_area_sei 
area_parcel$ci_upper_rracaterh_areai <- area_parcel$rracaterh_area_bi + 1.96*area_parcel$rracaterh_area_sei
area_parcel$ci_lower_rrmdfrrh_areai <- area_parcel$rrmdfrrh_area_bi - 1.96*area_parcel$rrmdfrrh_area_sei 
area_parcel$ci_upper_rrmdfrrh_areai <- area_parcel$rrmdfrrh_area_bi + 1.96*area_parcel$rrmdfrrh_area_sei
area_parcel$ci_lower_sufrrh_areai <- area_parcel$sufrrh_area_bi - 1.96*area_parcel$sufrrh_area_sei 
area_parcel$ci_upper_sufrrh_areai <- area_parcel$sufrrh_area_bi + 1.96*area_parcel$sufrrh_area_sei
area_parcel$ci_lower_suplrh_areai <- area_parcel$suplrh_area_bi - 1.96*area_parcel$suplrh_area_sei 
area_parcel$ci_upper_suplrh_areai <- area_parcel$suplrh_area_bi + 1.96*area_parcel$suplrh_area_sei
area_parcel$ci_lower_sutmrh_areai <- area_parcel$sutmrh_area_bi - 1.96*area_parcel$sutmrh_area_sei 
area_parcel$ci_upper_sutmrh_areai <- area_parcel$sutmrh_area_bi + 1.96*area_parcel$sutmrh_area_sei
area_parcel$ci_lower_smrh_areai <- area_parcel$smrh_area_bi - 1.96*area_parcel$smrh_area_sei 
area_parcel$ci_upper_smrh_areai <- area_parcel$smrh_area_bi + 1.96*area_parcel$smrh_area_sei
area_parcel$ci_lower_frpolerh_areai <- area_parcel$frpolerh_area_bi - 1.96*area_parcel$frpolerh_area_sei 
area_parcel$ci_upper_frpolerh_areai <- area_parcel$frpolerh_area_bi + 1.96*area_parcel$frpolerh_area_sei
area_parcel$ci_lower_tmpolerh_areai <- area_parcel$tmpolerh_area_bi - 1.96*area_parcel$tmpolerh_area_sei 
area_parcel$ci_upper_tmpolerh_areai <- area_parcel$tmpolerh_area_bi + 1.96*area_parcel$tmpolerh_area_sei
area_parcel$ci_lower_trvtmrh_areai <- area_parcel$trvtmrh_area_bi - 1.96*area_parcel$trvtmrh_area_sei 
area_parcel$ci_upper_trvtmrh_areai <- area_parcel$trvtmrh_area_bi + 1.96*area_parcel$trvtmrh_area_sei
area_parcel$ci_lower_insularh_areai <- area_parcel$insularh_area_bi - 1.96*area_parcel$insularh_area_sei 
area_parcel$ci_upper_insularh_areai <- area_parcel$insularh_area_bi + 1.96*area_parcel$insularh_area_sei

write.csv(area_parcel, "Parcel-Wise Cortical Area MLM Analysis Output Standardized_FINAL.csv") #write csv file

#Parcel-wise cortical thickness analyses

thick1 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_banksstslh_thick_cent + scale(subj_banksstslh_thick) + site_banksstslh_thick_cent*wave + scale(subj_banksstslh_thick)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b1i <- as.data.frame(c(summary(thick1)$coefficients[14,1]))
thick_se1i <- as.data.frame(c(summary(thick1)$coefficients[14,2]))
thick_p1i <- as.data.frame(c(summary(thick1)$coefficients[14,5]))
thick_b1s <- as.data.frame(c(summary(thick1)$coefficients[16,1]))
thick_se1s <- as.data.frame(c(summary(thick1)$coefficients[16,2]))
thick_p1s <- as.data.frame(c(summary(thick1)$coefficients[16,5]))

thick2 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_cdacatelh_thick_cent + scale(subj_cdacatelh_thick) + site_cdacatelh_thick_cent*wave + scale(subj_cdacatelh_thick)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b2i <- as.data.frame(c(summary(thick2)$coefficients[14,1]))
thick_se2i <- as.data.frame(c(summary(thick2)$coefficients[14,2]))
thick_p2i <- as.data.frame(c(summary(thick2)$coefficients[14,5]))
thick_b2s <- as.data.frame(c(summary(thick2)$coefficients[16,1]))
thick_se2s <- as.data.frame(c(summary(thick2)$coefficients[16,2]))
thick_p2s <- as.data.frame(c(summary(thick2)$coefficients[16,5]))

thick3 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_cdmdfrlh_thick_cent + scale(subj_cdmdfrlh_thick) + site_cdmdfrlh_thick_cent*wave + scale(subj_cdmdfrlh_thick)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b3i <- as.data.frame(c(summary(thick3)$coefficients[14,1]))
thick_se3i <- as.data.frame(c(summary(thick3)$coefficients[14,2]))
thick_p3i <- as.data.frame(c(summary(thick3)$coefficients[14,5]))
thick_b3s <- as.data.frame(c(summary(thick3)$coefficients[16,1]))
thick_se3s <- as.data.frame(c(summary(thick3)$coefficients[16,2]))
thick_p3s <- as.data.frame(c(summary(thick3)$coefficients[16,5]))

thick4 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_cuneuslh_thick_cent + scale(subj_cuneuslh_thick) + site_cuneuslh_thick_cent*wave + scale(subj_cuneuslh_thick)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b4i <- as.data.frame(c(summary(thick4)$coefficients[14,1]))
thick_se4i <- as.data.frame(c(summary(thick4)$coefficients[14,2]))
thick_p4i <- as.data.frame(c(summary(thick4)$coefficients[14,5]))
thick_b4s <- as.data.frame(c(summary(thick4)$coefficients[16,1]))
thick_se4s <- as.data.frame(c(summary(thick4)$coefficients[16,2]))
thick_p4s <- as.data.frame(c(summary(thick4)$coefficients[16,5]))

thick5 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_ehinallh_thick_cent + scale(subj_ehinallh_thick) + site_ehinallh_thick_cent*wave + scale(subj_ehinallh_thick)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b5i <- as.data.frame(c(summary(thick5)$coefficients[14,1]))
thick_se5i <- as.data.frame(c(summary(thick5)$coefficients[14,2]))
thick_p5i <- as.data.frame(c(summary(thick5)$coefficients[14,5]))
thick_b5s <- as.data.frame(c(summary(thick5)$coefficients[16,1]))
thick_se5s <- as.data.frame(c(summary(thick5)$coefficients[16,2]))
thick_p5s <- as.data.frame(c(summary(thick5)$coefficients[16,5]))

thick6 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_fusiformlh_thick_cent + scale(subj_fusiformlh_thick) + site_fusiformlh_thick_cent*wave + scale(subj_fusiformlh_thick)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b6i <- as.data.frame(c(summary(thick6)$coefficients[14,1]))
thick_se6i <- as.data.frame(c(summary(thick6)$coefficients[14,2]))
thick_p6i <- as.data.frame(c(summary(thick6)$coefficients[14,5]))
thick_b6s <- as.data.frame(c(summary(thick6)$coefficients[16,1]))
thick_se6s <- as.data.frame(c(summary(thick6)$coefficients[16,2]))
thick_p6s <- as.data.frame(c(summary(thick6)$coefficients[16,5]))

thick7 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_ifpllh_thick_cent + scale(subj_ifpllh_thick) + site_ifpllh_thick_cent*wave + scale(subj_ifpllh_thick)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b7i <- as.data.frame(c(summary(thick7)$coefficients[14,1]))
thick_se7i <- as.data.frame(c(summary(thick7)$coefficients[14,2]))
thick_p7i <- as.data.frame(c(summary(thick7)$coefficients[14,5]))
thick_b7s <- as.data.frame(c(summary(thick7)$coefficients[16,1]))
thick_se7s <- as.data.frame(c(summary(thick7)$coefficients[16,2]))
thick_p7s <- as.data.frame(c(summary(thick7)$coefficients[16,5]))

thick8 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_iftmlh_thick_cent + scale(subj_iftmlh_thick) + site_iftmlh_thick_cent*wave + scale(subj_iftmlh_thick)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b8i <- as.data.frame(c(summary(thick8)$coefficients[14,1]))
thick_se8i <- as.data.frame(c(summary(thick8)$coefficients[14,2]))
thick_p8i <- as.data.frame(c(summary(thick8)$coefficients[14,5]))
thick_b8s <- as.data.frame(c(summary(thick8)$coefficients[16,1]))
thick_se8s <- as.data.frame(c(summary(thick8)$coefficients[16,2]))
thick_p8s <- as.data.frame(c(summary(thick8)$coefficients[16,5]))

thick9 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                 site_ihcatelh_thick_cent + scale(subj_ihcatelh_thick) + site_ihcatelh_thick_cent*wave + scale(subj_ihcatelh_thick)*wave +
                 (1 + wave|site_id/id), 
               data = pfactor_red, 
               control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b9i <- as.data.frame(c(summary(thick9)$coefficients[14,1]))
thick_se9i <- as.data.frame(c(summary(thick9)$coefficients[14,2]))
thick_p9i <- as.data.frame(c(summary(thick9)$coefficients[14,5]))
thick_b9s <- as.data.frame(c(summary(thick9)$coefficients[16,1]))
thick_se9s <- as.data.frame(c(summary(thick9)$coefficients[16,2]))
thick_p9s <- as.data.frame(c(summary(thick9)$coefficients[16,5]))

thick10 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_locclh_thick_cent + scale(subj_locclh_thick) + site_locclh_thick_cent*wave + scale(subj_locclh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b10i <- as.data.frame(c(summary(thick10)$coefficients[14,1]))
thick_se10i <- as.data.frame(c(summary(thick10)$coefficients[14,2]))
thick_p10i <- as.data.frame(c(summary(thick10)$coefficients[14,5]))
thick_b10s <- as.data.frame(c(summary(thick10)$coefficients[16,1]))
thick_se10s <- as.data.frame(c(summary(thick10)$coefficients[16,2]))
thick_p10s <- as.data.frame(c(summary(thick10)$coefficients[16,5]))

thick11 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_lobfrlh_thick_cent + scale(subj_lobfrlh_thick) + site_lobfrlh_thick_cent*wave + scale(subj_lobfrlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b11i <- as.data.frame(c(summary(thick11)$coefficients[14,1]))
thick_se11i <- as.data.frame(c(summary(thick11)$coefficients[14,2]))
thick_p11i <- as.data.frame(c(summary(thick11)$coefficients[14,5]))
thick_b11s <- as.data.frame(c(summary(thick11)$coefficients[16,1]))
thick_se11s <- as.data.frame(c(summary(thick11)$coefficients[16,2]))
thick_p11s <- as.data.frame(c(summary(thick11)$coefficients[16,5]))

thick12 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_linguallh_thick_cent + scale(subj_linguallh_thick) + site_linguallh_thick_cent*wave + scale(subj_linguallh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b12i <- as.data.frame(c(summary(thick12)$coefficients[14,1]))
thick_se12i <- as.data.frame(c(summary(thick12)$coefficients[14,2]))
thick_p12i <- as.data.frame(c(summary(thick12)$coefficients[14,5]))
thick_b12s <- as.data.frame(c(summary(thick12)$coefficients[16,1]))
thick_se12s <- as.data.frame(c(summary(thick12)$coefficients[16,2]))
thick_p12s <- as.data.frame(c(summary(thick12)$coefficients[16,5]))

thick13 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_mobfrlh_thick_cent + scale(subj_mobfrlh_thick) + site_mobfrlh_thick_cent*wave + scale(subj_mobfrlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b13i <- as.data.frame(c(summary(thick13)$coefficients[14,1]))
thick_se13i <- as.data.frame(c(summary(thick13)$coefficients[14,2]))
thick_p13i <- as.data.frame(c(summary(thick13)$coefficients[14,5]))
thick_b13s <- as.data.frame(c(summary(thick13)$coefficients[16,1]))
thick_se13s <- as.data.frame(c(summary(thick13)$coefficients[16,2]))
thick_p13s <- as.data.frame(c(summary(thick13)$coefficients[16,5]))

thick14 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_mdtmlh_thick_cent + scale(subj_mdtmlh_thick) + site_mdtmlh_thick_cent*wave + scale(subj_mdtmlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b14i <- as.data.frame(c(summary(thick14)$coefficients[14,1]))
thick_se14i <- as.data.frame(c(summary(thick14)$coefficients[14,2]))
thick_p14i <- as.data.frame(c(summary(thick14)$coefficients[14,5]))
thick_b14s <- as.data.frame(c(summary(thick14)$coefficients[16,1]))
thick_se14s <- as.data.frame(c(summary(thick14)$coefficients[16,2]))
thick_p14s <- as.data.frame(c(summary(thick14)$coefficients[16,5]))

thick15 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_parahpallh_thick_cent + scale(subj_parahpallh_thick) + site_parahpallh_thick_cent*wave + scale(subj_parahpallh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b15i <- as.data.frame(c(summary(thick15)$coefficients[14,1]))
thick_se15i <- as.data.frame(c(summary(thick15)$coefficients[14,2]))
thick_p15i <- as.data.frame(c(summary(thick15)$coefficients[14,5]))
thick_b15s <- as.data.frame(c(summary(thick15)$coefficients[16,1]))
thick_se15s <- as.data.frame(c(summary(thick15)$coefficients[16,2]))
thick_p15s <- as.data.frame(c(summary(thick15)$coefficients[16,5]))

thick16 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_paracnlh_thick_cent + scale(subj_paracnlh_thick) + site_paracnlh_thick_cent*wave + scale(subj_paracnlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b16i <- as.data.frame(c(summary(thick16)$coefficients[14,1]))
thick_se16i <- as.data.frame(c(summary(thick16)$coefficients[14,2]))
thick_p16i <- as.data.frame(c(summary(thick16)$coefficients[14,5]))
thick_b16s <- as.data.frame(c(summary(thick16)$coefficients[16,1]))
thick_se16s <- as.data.frame(c(summary(thick16)$coefficients[16,2]))
thick_p16s <- as.data.frame(c(summary(thick16)$coefficients[16,5]))

thick17 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_parsopclh_thick_cent + scale(subj_parsopclh_thick) + site_parsopclh_thick_cent*wave + scale(subj_parsopclh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b17i <- as.data.frame(c(summary(thick17)$coefficients[14,1]))
thick_se17i <- as.data.frame(c(summary(thick17)$coefficients[14,2]))
thick_p17i <- as.data.frame(c(summary(thick17)$coefficients[14,5]))
thick_b17s <- as.data.frame(c(summary(thick17)$coefficients[16,1]))
thick_se17s <- as.data.frame(c(summary(thick17)$coefficients[16,2]))
thick_p17s <- as.data.frame(c(summary(thick17)$coefficients[16,5]))

thick18 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_parsobislh_thick_cent + scale(subj_parsobislh_thick) + site_parsobislh_thick_cent*wave + scale(subj_parsobislh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b18i <- as.data.frame(c(summary(thick18)$coefficients[14,1]))
thick_se18i <- as.data.frame(c(summary(thick18)$coefficients[14,2]))
thick_p18i <- as.data.frame(c(summary(thick18)$coefficients[14,5]))
thick_b18s <- as.data.frame(c(summary(thick18)$coefficients[16,1]))
thick_se18s <- as.data.frame(c(summary(thick18)$coefficients[16,2]))
thick_p18s <- as.data.frame(c(summary(thick18)$coefficients[16,5]))

thick19 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_parstgrislh_thick_cent + scale(subj_parstgrislh_thick) + site_parstgrislh_thick_cent*wave + scale(subj_parstgrislh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b19i <- as.data.frame(c(summary(thick19)$coefficients[14,1]))
thick_se19i <- as.data.frame(c(summary(thick19)$coefficients[14,2]))
thick_p19i <- as.data.frame(c(summary(thick19)$coefficients[14,5]))
thick_b19s <- as.data.frame(c(summary(thick19)$coefficients[16,1]))
thick_se19s <- as.data.frame(c(summary(thick19)$coefficients[16,2]))
thick_p19s <- as.data.frame(c(summary(thick19)$coefficients[16,5]))


thick20 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_pericclh_thick_cent + scale(subj_pericclh_thick) + site_pericclh_thick_cent*wave + scale(subj_pericclh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b20i <- as.data.frame(c(summary(thick20)$coefficients[14,1]))
thick_se20i <- as.data.frame(c(summary(thick20)$coefficients[14,2]))
thick_p20i <- as.data.frame(c(summary(thick20)$coefficients[14,5]))
thick_b20s <- as.data.frame(c(summary(thick20)$coefficients[16,1]))
thick_se20s <- as.data.frame(c(summary(thick20)$coefficients[16,2]))
thick_p20s <- as.data.frame(c(summary(thick20)$coefficients[16,5]))

thick21 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_postcnlh_thick_cent + scale(subj_postcnlh_thick) + site_postcnlh_thick_cent*wave + scale(subj_postcnlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b21i <- as.data.frame(c(summary(thick21)$coefficients[14,1]))
thick_se21i <- as.data.frame(c(summary(thick21)$coefficients[14,2]))
thick_p21i <- as.data.frame(c(summary(thick21)$coefficients[14,5]))
thick_b21s <- as.data.frame(c(summary(thick21)$coefficients[16,1]))
thick_se21s <- as.data.frame(c(summary(thick21)$coefficients[16,2]))
thick_p21s <- as.data.frame(c(summary(thick21)$coefficients[16,5]))

thick22 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_ptcatelh_thick_cent + scale(subj_ptcatelh_thick) + site_ptcatelh_thick_cent*wave + scale(subj_ptcatelh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b22i <- as.data.frame(c(summary(thick22)$coefficients[14,1]))
thick_se22i <- as.data.frame(c(summary(thick22)$coefficients[14,2]))
thick_p22i <- as.data.frame(c(summary(thick22)$coefficients[14,5]))
thick_b22s <- as.data.frame(c(summary(thick22)$coefficients[16,1]))
thick_se22s <- as.data.frame(c(summary(thick22)$coefficients[16,2]))
thick_p22s <- as.data.frame(c(summary(thick22)$coefficients[16,5]))

thick23 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_precnlh_thick_cent + scale(subj_precnlh_thick) + site_precnlh_thick_cent*wave + scale(subj_precnlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b23i <- as.data.frame(c(summary(thick23)$coefficients[14,1]))
thick_se23i <- as.data.frame(c(summary(thick23)$coefficients[14,2]))
thick_p23i <- as.data.frame(c(summary(thick23)$coefficients[14,5]))
thick_b23s <- as.data.frame(c(summary(thick23)$coefficients[16,1]))
thick_se23s <- as.data.frame(c(summary(thick23)$coefficients[16,2]))
thick_p23s <- as.data.frame(c(summary(thick23)$coefficients[16,5]))

thick24 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_pclh_thick_cent + scale(subj_pclh_thick) + site_pclh_thick_cent*wave + scale(subj_pclh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b24i <- as.data.frame(c(summary(thick24)$coefficients[14,1]))
thick_se24i <- as.data.frame(c(summary(thick24)$coefficients[14,2]))
thick_p24i <- as.data.frame(c(summary(thick24)$coefficients[14,5]))
thick_b24s <- as.data.frame(c(summary(thick24)$coefficients[16,1]))
thick_se24s <- as.data.frame(c(summary(thick24)$coefficients[16,2]))
thick_p24s <- as.data.frame(c(summary(thick24)$coefficients[16,5]))

thick25 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_rracatelh_thick_cent + scale(subj_rracatelh_thick) + site_rracatelh_thick_cent*wave + scale(subj_rracatelh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b25i <- as.data.frame(c(summary(thick25)$coefficients[14,1]))
thick_se25i <- as.data.frame(c(summary(thick25)$coefficients[14,2]))
thick_p25i <- as.data.frame(c(summary(thick25)$coefficients[14,5]))
thick_b25s <- as.data.frame(c(summary(thick25)$coefficients[16,1]))
thick_se25s <- as.data.frame(c(summary(thick25)$coefficients[16,2]))
thick_p25s <- as.data.frame(c(summary(thick25)$coefficients[16,5]))

thick26 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_rrmdfrlh_thick_cent + scale(subj_rrmdfrlh_thick) + site_rrmdfrlh_thick_cent*wave + scale(subj_rrmdfrlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b26i <- as.data.frame(c(summary(thick26)$coefficients[14,1]))
thick_se26i <- as.data.frame(c(summary(thick26)$coefficients[14,2]))
thick_p26i <- as.data.frame(c(summary(thick26)$coefficients[14,5]))
thick_b26s <- as.data.frame(c(summary(thick26)$coefficients[16,1]))
thick_se26s <- as.data.frame(c(summary(thick26)$coefficients[16,2]))
thick_p26s <- as.data.frame(c(summary(thick26)$coefficients[16,5]))

thick27 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_sufrlh_thick_cent + scale(subj_sufrlh_thick) + site_sufrlh_thick_cent*wave + scale(subj_sufrlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b27i <- as.data.frame(c(summary(thick27)$coefficients[14,1]))
thick_se27i <- as.data.frame(c(summary(thick27)$coefficients[14,2]))
thick_p27i <- as.data.frame(c(summary(thick27)$coefficients[14,5]))
thick_b27s <- as.data.frame(c(summary(thick27)$coefficients[16,1]))
thick_se27s <- as.data.frame(c(summary(thick27)$coefficients[16,2]))
thick_p27s <- as.data.frame(c(summary(thick27)$coefficients[16,5]))

thick28 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_supllh_thick_cent + scale(subj_supllh_thick) + site_supllh_thick_cent*wave + scale(subj_supllh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b28i <- as.data.frame(c(summary(thick28)$coefficients[14,1]))
thick_se28i <- as.data.frame(c(summary(thick28)$coefficients[14,2]))
thick_p28i <- as.data.frame(c(summary(thick28)$coefficients[14,5]))
thick_b28s <- as.data.frame(c(summary(thick28)$coefficients[16,1]))
thick_se28s <- as.data.frame(c(summary(thick28)$coefficients[16,2]))
thick_p28s <- as.data.frame(c(summary(thick28)$coefficients[16,5]))

thick29 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_sutmlh_thick_cent + scale(subj_sutmlh_thick) + site_sutmlh_thick_cent*wave + scale(subj_sutmlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b29i <- as.data.frame(c(summary(thick29)$coefficients[14,1]))
thick_se29i <- as.data.frame(c(summary(thick29)$coefficients[14,2]))
thick_p29i <- as.data.frame(c(summary(thick29)$coefficients[14,5]))
thick_b29s <- as.data.frame(c(summary(thick29)$coefficients[16,1]))
thick_se29s <- as.data.frame(c(summary(thick29)$coefficients[16,2]))
thick_p29s <- as.data.frame(c(summary(thick29)$coefficients[16,5]))

thick30 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_smlh_thick_cent + scale(subj_smlh_thick) + site_smlh_thick_cent*wave + scale(subj_smlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b30i <- as.data.frame(c(summary(thick30)$coefficients[14,1]))
thick_se30i <- as.data.frame(c(summary(thick30)$coefficients[14,2]))
thick_p30i <- as.data.frame(c(summary(thick30)$coefficients[14,5]))
thick_b30s <- as.data.frame(c(summary(thick30)$coefficients[16,1]))
thick_se30s <- as.data.frame(c(summary(thick30)$coefficients[16,2]))
thick_p30s <- as.data.frame(c(summary(thick30)$coefficients[16,5]))

thick31 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_frpolelh_thick_cent + scale(subj_frpolelh_thick) + site_frpolelh_thick_cent*wave + scale(subj_frpolelh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b31i <- as.data.frame(c(summary(thick31)$coefficients[14,1]))
thick_se31i <- as.data.frame(c(summary(thick31)$coefficients[14,2]))
thick_p31i <- as.data.frame(c(summary(thick31)$coefficients[14,5]))
thick_b31s <- as.data.frame(c(summary(thick31)$coefficients[16,1]))
thick_se31s <- as.data.frame(c(summary(thick31)$coefficients[16,2]))
thick_p31s <- as.data.frame(c(summary(thick31)$coefficients[16,5]))

thick32 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_tmpolelh_thick_cent + scale(subj_tmpolelh_thick) + site_tmpolelh_thick_cent*wave + scale(subj_tmpolelh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b32i <- as.data.frame(c(summary(thick32)$coefficients[14,1]))
thick_se32i <- as.data.frame(c(summary(thick32)$coefficients[14,2]))
thick_p32i <- as.data.frame(c(summary(thick32)$coefficients[14,5]))
thick_b32s <- as.data.frame(c(summary(thick32)$coefficients[16,1]))
thick_se32s <- as.data.frame(c(summary(thick32)$coefficients[16,2]))
thick_p32s <- as.data.frame(c(summary(thick32)$coefficients[16,5]))

thick33 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_trvtmlh_thick_cent + scale(subj_trvtmlh_thick) + site_trvtmlh_thick_cent*wave + scale(subj_trvtmlh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b33i <- as.data.frame(c(summary(thick33)$coefficients[14,1]))
thick_se33i <- as.data.frame(c(summary(thick33)$coefficients[14,2]))
thick_p33i <- as.data.frame(c(summary(thick33)$coefficients[14,5]))
thick_b33s <- as.data.frame(c(summary(thick33)$coefficients[16,1]))
thick_se33s <- as.data.frame(c(summary(thick33)$coefficients[16,2]))
thick_p33s <- as.data.frame(c(summary(thick33)$coefficients[16,5]))

thick34 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_insulalh_thick_cent + scale(subj_insulalh_thick) + site_insulalh_thick_cent*wave + scale(subj_insulalh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b34i <- as.data.frame(c(summary(thick34)$coefficients[14,1]))
thick_se34i <- as.data.frame(c(summary(thick34)$coefficients[14,2]))
thick_p34i <- as.data.frame(c(summary(thick34)$coefficients[14,5]))
thick_b34s <- as.data.frame(c(summary(thick34)$coefficients[16,1]))
thick_se34s <- as.data.frame(c(summary(thick34)$coefficients[16,2]))
thick_p34s <- as.data.frame(c(summary(thick34)$coefficients[16,5]))

thick35 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_banksstsrh_thick_cent + scale(subj_banksstsrh_thick) + site_banksstsrh_thick_cent*wave + scale(subj_banksstsrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b35i <- as.data.frame(c(summary(thick35)$coefficients[14,1]))
thick_se35i <- as.data.frame(c(summary(thick35)$coefficients[14,2]))
thick_p35i <- as.data.frame(c(summary(thick35)$coefficients[14,5]))
thick_b35s <- as.data.frame(c(summary(thick35)$coefficients[16,1]))
thick_se35s <- as.data.frame(c(summary(thick35)$coefficients[16,2]))
thick_p35s <- as.data.frame(c(summary(thick35)$coefficients[16,5]))

thick36 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_cdacaterh_thick_cent + scale(subj_cdacaterh_thick) + site_cdacaterh_thick_cent*wave + scale(subj_cdacaterh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b36i <- as.data.frame(c(summary(thick36)$coefficients[14,1]))
thick_se36i <- as.data.frame(c(summary(thick36)$coefficients[14,2]))
thick_p36i <- as.data.frame(c(summary(thick36)$coefficients[14,5]))
thick_b36s <- as.data.frame(c(summary(thick36)$coefficients[16,1]))
thick_se36s <- as.data.frame(c(summary(thick36)$coefficients[16,2]))
thick_p36s <- as.data.frame(c(summary(thick36)$coefficients[16,5]))

thick37 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_cdmdfrrh_thick_cent + scale(subj_cdmdfrrh_thick) + site_cdmdfrrh_thick_cent*wave + scale(subj_cdmdfrrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b37i <- as.data.frame(c(summary(thick37)$coefficients[14,1]))
thick_se37i <- as.data.frame(c(summary(thick37)$coefficients[14,2]))
thick_p37i <- as.data.frame(c(summary(thick37)$coefficients[14,5]))
thick_b37s <- as.data.frame(c(summary(thick37)$coefficients[16,1]))
thick_se37s <- as.data.frame(c(summary(thick37)$coefficients[16,2]))
thick_p37s <- as.data.frame(c(summary(thick37)$coefficients[16,5]))

thick38 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_cuneusrh_thick_cent + scale(subj_cuneusrh_thick) + site_cuneusrh_thick_cent*wave + scale(subj_cuneusrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b38i <- as.data.frame(c(summary(thick38)$coefficients[14,1]))
thick_se38i <- as.data.frame(c(summary(thick38)$coefficients[14,2]))
thick_p38i <- as.data.frame(c(summary(thick38)$coefficients[14,5]))
thick_b38s <- as.data.frame(c(summary(thick38)$coefficients[16,1]))
thick_se38s <- as.data.frame(c(summary(thick38)$coefficients[16,2]))
thick_p38s <- as.data.frame(c(summary(thick38)$coefficients[16,5]))

thick39 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_ehinalrh_thick_cent + scale(subj_ehinalrh_thick) + site_ehinalrh_thick_cent*wave + scale(subj_ehinalrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b39i <- as.data.frame(c(summary(thick39)$coefficients[14,1]))
thick_se39i <- as.data.frame(c(summary(thick39)$coefficients[14,2]))
thick_p39i <- as.data.frame(c(summary(thick39)$coefficients[14,5]))
thick_b39s <- as.data.frame(c(summary(thick39)$coefficients[16,1]))
thick_se39s <- as.data.frame(c(summary(thick39)$coefficients[16,2]))
thick_p39s <- as.data.frame(c(summary(thick39)$coefficients[16,5]))

thick40 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_fusiformrh_thick_cent + scale(subj_fusiformrh_thick) + site_fusiformrh_thick_cent*wave + scale(subj_fusiformrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b40i <- as.data.frame(c(summary(thick40)$coefficients[14,1]))
thick_se40i <- as.data.frame(c(summary(thick40)$coefficients[14,2]))
thick_p40i <- as.data.frame(c(summary(thick40)$coefficients[14,5]))
thick_b40s <- as.data.frame(c(summary(thick40)$coefficients[16,1]))
thick_se40s <- as.data.frame(c(summary(thick40)$coefficients[16,2]))
thick_p40s <- as.data.frame(c(summary(thick40)$coefficients[16,5]))

thick41 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_ifplrh_thick_cent + scale(subj_ifplrh_thick) + site_ifplrh_thick_cent*wave + scale(subj_ifplrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b41i <- as.data.frame(c(summary(thick41)$coefficients[14,1]))
thick_se41i <- as.data.frame(c(summary(thick41)$coefficients[14,2]))
thick_p41i <- as.data.frame(c(summary(thick41)$coefficients[14,5]))
thick_b41s <- as.data.frame(c(summary(thick41)$coefficients[16,1]))
thick_se41s <- as.data.frame(c(summary(thick41)$coefficients[16,2]))
thick_p41s <- as.data.frame(c(summary(thick41)$coefficients[16,5]))

thick42 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_iftmrh_thick_cent + scale(subj_iftmrh_thick) + site_iftmrh_thick_cent*wave + scale(subj_iftmrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b42i <- as.data.frame(c(summary(thick42)$coefficients[14,1]))
thick_se42i <- as.data.frame(c(summary(thick42)$coefficients[14,2]))
thick_p42i <- as.data.frame(c(summary(thick42)$coefficients[14,5]))
thick_b42s <- as.data.frame(c(summary(thick42)$coefficients[16,1]))
thick_se42s <- as.data.frame(c(summary(thick42)$coefficients[16,2]))
thick_p42s <- as.data.frame(c(summary(thick42)$coefficients[16,5]))

thick43 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_ihcaterh_thick_cent + scale(subj_ihcaterh_thick) + site_ihcaterh_thick_cent*wave + scale(subj_ihcaterh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b43i <- as.data.frame(c(summary(thick43)$coefficients[14,1]))
thick_se43i <- as.data.frame(c(summary(thick43)$coefficients[14,2]))
thick_p43i <- as.data.frame(c(summary(thick43)$coefficients[14,5]))
thick_b43s <- as.data.frame(c(summary(thick43)$coefficients[16,1]))
thick_se43s <- as.data.frame(c(summary(thick43)$coefficients[16,2]))
thick_p43s <- as.data.frame(c(summary(thick43)$coefficients[16,5]))

thick44 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_loccrh_thick_cent + scale(subj_loccrh_thick) + site_loccrh_thick_cent*wave + scale(subj_loccrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b44i <- as.data.frame(c(summary(thick44)$coefficients[14,1]))
thick_se44i <- as.data.frame(c(summary(thick44)$coefficients[14,2]))
thick_p44i <- as.data.frame(c(summary(thick44)$coefficients[14,5]))
thick_b44s <- as.data.frame(c(summary(thick44)$coefficients[16,1]))
thick_se44s <- as.data.frame(c(summary(thick44)$coefficients[16,2]))
thick_p44s <- as.data.frame(c(summary(thick44)$coefficients[16,5]))

thick45 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_lobfrrh_thick_cent + scale(subj_lobfrrh_thick) + site_lobfrrh_thick_cent*wave + scale(subj_lobfrrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b45i <- as.data.frame(c(summary(thick45)$coefficients[14,1]))
thick_se45i <- as.data.frame(c(summary(thick45)$coefficients[14,2]))
thick_p45i <- as.data.frame(c(summary(thick45)$coefficients[14,5]))
thick_b45s <- as.data.frame(c(summary(thick45)$coefficients[16,1]))
thick_se45s <- as.data.frame(c(summary(thick45)$coefficients[16,2]))
thick_p45s <- as.data.frame(c(summary(thick45)$coefficients[16,5]))

thick46 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_lingualrh_thick_cent + scale(subj_lingualrh_thick) + site_lingualrh_thick_cent*wave + scale(subj_lingualrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b46i <- as.data.frame(c(summary(thick46)$coefficients[14,1]))
thick_se46i <- as.data.frame(c(summary(thick46)$coefficients[14,2]))
thick_p46i <- as.data.frame(c(summary(thick46)$coefficients[14,5]))
thick_b46s <- as.data.frame(c(summary(thick46)$coefficients[16,1]))
thick_se46s <- as.data.frame(c(summary(thick46)$coefficients[16,2]))
thick_p46s <- as.data.frame(c(summary(thick46)$coefficients[16,5]))

thick47 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_mobfrrh_thick_cent + scale(subj_mobfrrh_thick) + site_mobfrrh_thick_cent*wave + scale(subj_mobfrrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b47i <- as.data.frame(c(summary(thick47)$coefficients[14,1]))
thick_se47i <- as.data.frame(c(summary(thick47)$coefficients[14,2]))
thick_p47i <- as.data.frame(c(summary(thick47)$coefficients[14,5]))
thick_b47s <- as.data.frame(c(summary(thick47)$coefficients[16,1]))
thick_se47s <- as.data.frame(c(summary(thick47)$coefficients[16,2]))
thick_p47s <- as.data.frame(c(summary(thick47)$coefficients[16,5]))

thick48 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_mdtmrh_thick_cent + scale(subj_mdtmrh_thick) + site_mdtmrh_thick_cent*wave + scale(subj_mdtmrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b48i <- as.data.frame(c(summary(thick48)$coefficients[14,1]))
thick_se48i <- as.data.frame(c(summary(thick48)$coefficients[14,2]))
thick_p48i <- as.data.frame(c(summary(thick48)$coefficients[14,5]))
thick_b48s <- as.data.frame(c(summary(thick48)$coefficients[16,1]))
thick_se48s <- as.data.frame(c(summary(thick48)$coefficients[16,2]))
thick_p48s <- as.data.frame(c(summary(thick48)$coefficients[16,5]))

thick49 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_parahpalrh_thick_cent + scale(subj_parahpalrh_thick) + site_parahpalrh_thick_cent*wave + scale(subj_parahpalrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b49i <- as.data.frame(c(summary(thick49)$coefficients[14,1]))
thick_se49i <- as.data.frame(c(summary(thick49)$coefficients[14,2]))
thick_p49i <- as.data.frame(c(summary(thick49)$coefficients[14,5]))
thick_b49s <- as.data.frame(c(summary(thick49)$coefficients[16,1]))
thick_se49s <- as.data.frame(c(summary(thick49)$coefficients[16,2]))
thick_p49s <- as.data.frame(c(summary(thick49)$coefficients[16,5]))

thick50 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_paracnrh_thick_cent + scale(subj_paracnrh_thick) + site_paracnrh_thick_cent*wave + scale(subj_paracnrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b50i <- as.data.frame(c(summary(thick50)$coefficients[14,1]))
thick_se50i <- as.data.frame(c(summary(thick50)$coefficients[14,2]))
thick_p50i <- as.data.frame(c(summary(thick50)$coefficients[14,5]))
thick_b50s <- as.data.frame(c(summary(thick50)$coefficients[16,1]))
thick_se50s <- as.data.frame(c(summary(thick50)$coefficients[16,2]))
thick_p50s <- as.data.frame(c(summary(thick50)$coefficients[16,5]))

thick51 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_parsopcrh_thick_cent + scale(subj_parsopcrh_thick) + site_parsopcrh_thick_cent*wave + scale(subj_parsopcrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b51i <- as.data.frame(c(summary(thick51)$coefficients[14,1]))
thick_se51i <- as.data.frame(c(summary(thick51)$coefficients[14,2]))
thick_p51i <- as.data.frame(c(summary(thick51)$coefficients[14,5]))
thick_b51s <- as.data.frame(c(summary(thick51)$coefficients[16,1]))
thick_se51s <- as.data.frame(c(summary(thick51)$coefficients[16,2]))
thick_p51s <- as.data.frame(c(summary(thick51)$coefficients[16,5]))

thick52 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_parsobisrh_thick_cent + scale(subj_parsobisrh_thick) + site_parsobisrh_thick_cent*wave + scale(subj_parsobisrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b52i <- as.data.frame(c(summary(thick52)$coefficients[14,1]))
thick_se52i <- as.data.frame(c(summary(thick52)$coefficients[14,2]))
thick_p52i <- as.data.frame(c(summary(thick52)$coefficients[14,5]))
thick_b52s <- as.data.frame(c(summary(thick52)$coefficients[16,1]))
thick_se52s <- as.data.frame(c(summary(thick52)$coefficients[16,2]))
thick_p52s <- as.data.frame(c(summary(thick52)$coefficients[16,5]))

thick53 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_parstgrisrh_thick_cent + scale(subj_parstgrisrh_thick) + site_parstgrisrh_thick_cent*wave + scale(subj_parstgrisrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b53i <- as.data.frame(c(summary(thick53)$coefficients[14,1]))
thick_se53i <- as.data.frame(c(summary(thick53)$coefficients[14,2]))
thick_p53i <- as.data.frame(c(summary(thick53)$coefficients[14,5]))
thick_b53s <- as.data.frame(c(summary(thick53)$coefficients[16,1]))
thick_se53s <- as.data.frame(c(summary(thick53)$coefficients[16,2]))
thick_p53s <- as.data.frame(c(summary(thick53)$coefficients[16,5]))

thick54 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_periccrh_thick_cent + scale(subj_periccrh_thick) + site_periccrh_thick_cent*wave + scale(subj_periccrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b54i <- as.data.frame(c(summary(thick54)$coefficients[14,1]))
thick_se54i <- as.data.frame(c(summary(thick54)$coefficients[14,2]))
thick_p54i <- as.data.frame(c(summary(thick54)$coefficients[14,5]))
thick_b54s <- as.data.frame(c(summary(thick54)$coefficients[16,1]))
thick_se54s <- as.data.frame(c(summary(thick54)$coefficients[16,2]))
thick_p54s <- as.data.frame(c(summary(thick54)$coefficients[16,5]))

thick55 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_postcnrh_thick_cent + scale(subj_postcnrh_thick) + site_postcnrh_thick_cent*wave + scale(subj_postcnrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b55i <- as.data.frame(c(summary(thick55)$coefficients[14,1]))
thick_se55i <- as.data.frame(c(summary(thick55)$coefficients[14,2]))
thick_p55i <- as.data.frame(c(summary(thick55)$coefficients[14,5]))
thick_b55s <- as.data.frame(c(summary(thick55)$coefficients[16,1]))
thick_se55s <- as.data.frame(c(summary(thick55)$coefficients[16,2]))
thick_p55s <- as.data.frame(c(summary(thick55)$coefficients[16,5]))

thick56 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_ptcaterh_thick_cent + scale(subj_ptcaterh_thick) + site_ptcaterh_thick_cent*wave + scale(subj_ptcaterh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b56i <- as.data.frame(c(summary(thick56)$coefficients[14,1]))
thick_se56i <- as.data.frame(c(summary(thick56)$coefficients[14,2]))
thick_p56i <- as.data.frame(c(summary(thick56)$coefficients[14,5]))
thick_b56s <- as.data.frame(c(summary(thick56)$coefficients[16,1]))
thick_se56s <- as.data.frame(c(summary(thick56)$coefficients[16,2]))
thick_p56s <- as.data.frame(c(summary(thick56)$coefficients[16,5]))

thick57 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_precnrh_thick_cent + scale(subj_precnrh_thick) + site_precnrh_thick_cent*wave + scale(subj_precnrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b57i <- as.data.frame(c(summary(thick57)$coefficients[14,1]))
thick_se57i <- as.data.frame(c(summary(thick57)$coefficients[14,2]))
thick_p57i <- as.data.frame(c(summary(thick57)$coefficients[14,5]))
thick_b57s <- as.data.frame(c(summary(thick57)$coefficients[16,1]))
thick_se57s <- as.data.frame(c(summary(thick57)$coefficients[16,2]))
thick_p57s <- as.data.frame(c(summary(thick57)$coefficients[16,5]))

thick58 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_pcrh_thick_cent + scale(subj_pcrh_thick) + site_pcrh_thick_cent*wave + scale(subj_pcrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b58i <- as.data.frame(c(summary(thick58)$coefficients[14,1]))
thick_se58i <- as.data.frame(c(summary(thick58)$coefficients[14,2]))
thick_p58i <- as.data.frame(c(summary(thick58)$coefficients[14,5]))
thick_b58s <- as.data.frame(c(summary(thick58)$coefficients[16,1]))
thick_se58s <- as.data.frame(c(summary(thick58)$coefficients[16,2]))
thick_p58s <- as.data.frame(c(summary(thick58)$coefficients[16,5]))

thick59 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_rracaterh_thick_cent + scale(subj_rracaterh_thick) + site_rracaterh_thick_cent*wave + scale(subj_rracaterh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b59i <- as.data.frame(c(summary(thick59)$coefficients[14,1]))
thick_se59i <- as.data.frame(c(summary(thick59)$coefficients[14,2]))
thick_p59i <- as.data.frame(c(summary(thick59)$coefficients[14,5]))
thick_b59s <- as.data.frame(c(summary(thick59)$coefficients[16,1]))
thick_se59s <- as.data.frame(c(summary(thick59)$coefficients[16,2]))
thick_p59s <- as.data.frame(c(summary(thick59)$coefficients[16,5]))

thick60 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_rrmdfrrh_thick_cent + scale(subj_rrmdfrrh_thick) + site_rrmdfrrh_thick_cent*wave + scale(subj_rrmdfrrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b60i <- as.data.frame(c(summary(thick60)$coefficients[14,1]))
thick_se60i <- as.data.frame(c(summary(thick60)$coefficients[14,2]))
thick_p60i <- as.data.frame(c(summary(thick60)$coefficients[14,5]))
thick_b60s <- as.data.frame(c(summary(thick60)$coefficients[16,1]))
thick_se60s <- as.data.frame(c(summary(thick60)$coefficients[16,2]))
thick_p60s <- as.data.frame(c(summary(thick60)$coefficients[16,5]))

thick61 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_sufrrh_thick_cent + scale(subj_sufrrh_thick) + site_sufrrh_thick_cent*wave + scale(subj_sufrrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b61i <- as.data.frame(c(summary(thick61)$coefficients[14,1]))
thick_se61i <- as.data.frame(c(summary(thick61)$coefficients[14,2]))
thick_p61i <- as.data.frame(c(summary(thick61)$coefficients[14,5]))
thick_b61s <- as.data.frame(c(summary(thick61)$coefficients[16,1]))
thick_se61s <- as.data.frame(c(summary(thick61)$coefficients[16,2]))
thick_p61s <- as.data.frame(c(summary(thick61)$coefficients[16,5]))

thick62 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_suplrh_thick_cent + scale(subj_suplrh_thick) + site_suplrh_thick_cent*wave + scale(subj_suplrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b62i <- as.data.frame(c(summary(thick62)$coefficients[14,1]))
thick_se62i <- as.data.frame(c(summary(thick62)$coefficients[14,2]))
thick_p62i <- as.data.frame(c(summary(thick62)$coefficients[14,5]))
thick_b62s <- as.data.frame(c(summary(thick62)$coefficients[16,1]))
thick_se62s <- as.data.frame(c(summary(thick62)$coefficients[16,2]))
thick_p62s <- as.data.frame(c(summary(thick62)$coefficients[16,5]))

thick63 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_sutmrh_thick_cent + scale(subj_sutmrh_thick) + site_sutmrh_thick_cent*wave + scale(subj_sutmrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b63i <- as.data.frame(c(summary(thick63)$coefficients[14,1]))
thick_se63i <- as.data.frame(c(summary(thick63)$coefficients[14,2]))
thick_p63i <- as.data.frame(c(summary(thick63)$coefficients[14,5]))
thick_b63s <- as.data.frame(c(summary(thick63)$coefficients[16,1]))
thick_se63s <- as.data.frame(c(summary(thick63)$coefficients[16,2]))
thick_p63s <- as.data.frame(c(summary(thick63)$coefficients[16,5]))

thick64 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_smrh_thick_cent + scale(subj_smrh_thick) + site_smrh_thick_cent*wave + scale(subj_smrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b64i <- as.data.frame(c(summary(thick64)$coefficients[14,1]))
thick_se64i <- as.data.frame(c(summary(thick64)$coefficients[14,2]))
thick_p64i <- as.data.frame(c(summary(thick64)$coefficients[14,5]))
thick_b64s <- as.data.frame(c(summary(thick64)$coefficients[16,1]))
thick_se64s <- as.data.frame(c(summary(thick64)$coefficients[16,2]))
thick_p64s <- as.data.frame(c(summary(thick64)$coefficients[16,5]))

thick65 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_frpolerh_thick_cent + scale(subj_frpolerh_thick) + site_frpolerh_thick_cent*wave + scale(subj_frpolerh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b65i <- as.data.frame(c(summary(thick65)$coefficients[14,1]))
thick_se65i <- as.data.frame(c(summary(thick65)$coefficients[14,2]))
thick_p65i <- as.data.frame(c(summary(thick65)$coefficients[14,5]))
thick_b65s <- as.data.frame(c(summary(thick65)$coefficients[16,1]))
thick_se65s <- as.data.frame(c(summary(thick65)$coefficients[16,2]))
thick_p65s <- as.data.frame(c(summary(thick65)$coefficients[16,5]))

thick66 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_tmpolerh_thick_cent + scale(subj_tmpolerh_thick) + site_tmpolerh_thick_cent*wave + scale(subj_tmpolerh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b66i <- as.data.frame(c(summary(thick66)$coefficients[14,1]))
thick_se66i <- as.data.frame(c(summary(thick66)$coefficients[14,2]))
thick_p66i <- as.data.frame(c(summary(thick66)$coefficients[14,5]))
thick_b66s <- as.data.frame(c(summary(thick66)$coefficients[16,1]))
thick_se66s <- as.data.frame(c(summary(thick66)$coefficients[16,2]))
thick_p66s <- as.data.frame(c(summary(thick66)$coefficients[16,5]))

thick67 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_trvtmrh_thick_cent + scale(subj_trvtmrh_thick) + site_trvtmrh_thick_cent*wave + scale(subj_trvtmrh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b67i <- as.data.frame(c(summary(thick67)$coefficients[14,1]))
thick_se67i <- as.data.frame(c(summary(thick67)$coefficients[14,2]))
thick_p67i <- as.data.frame(c(summary(thick67)$coefficients[14,5]))
thick_b67s <- as.data.frame(c(summary(thick67)$coefficients[16,1]))
thick_se67s <- as.data.frame(c(summary(thick67)$coefficients[16,2]))
thick_p67s <- as.data.frame(c(summary(thick67)$coefficients[16,5]))

thick68 <- lmer(scaleint ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                  site_insularh_thick_cent + scale(subj_insularh_thick) + site_insularh_thick_cent*wave + scale(subj_insularh_thick)*wave +
                  (1 + wave|site_id/id), 
                data = pfactor_red, 
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
thick_b68i <- as.data.frame(c(summary(thick68)$coefficients[14,1]))
thick_se68i <- as.data.frame(c(summary(thick68)$coefficients[14,2]))
thick_p68i <- as.data.frame(c(summary(thick68)$coefficients[14,5]))
thick_b68s <- as.data.frame(c(summary(thick68)$coefficients[16,1]))
thick_se68s <- as.data.frame(c(summary(thick68)$coefficients[16,2]))
thick_p68s <- as.data.frame(c(summary(thick68)$coefficients[16,5]))

#create data frame with all parcel-wise cortical thickume st. estimates, SEs, and p-values
thick_parcel <- data.frame(x=c("thick"))

newthickbi <- c(thick_b1i,thick_b2i,thick_b3i,thick_b4i,thick_b5i,thick_b6i,thick_b7i,thick_b8i,thick_b9i,thick_b10i,thick_b11i,thick_b12i,thick_b13i,thick_b14i,thick_b15i,
                thick_b16i,thick_b17i,thick_b18i,thick_b19i,thick_b20i,thick_b21i,thick_b22i,thick_b23i,thick_b24i,thick_b25i,thick_b26i,thick_b27i,thick_b28i,thick_b29i,
                thick_b30i,thick_b31i,thick_b32i,thick_b33i,thick_b34i,thick_b35i,thick_b36i,thick_b37i,thick_b38i,thick_b39i,thick_b40i,thick_b41i,thick_b42i,thick_b43i,
                thick_b44i,thick_b45i,thick_b46i,thick_b47i,thick_b48i,thick_b49i,thick_b50i,thick_b51i,thick_b52i,thick_b53i,thick_b54i,thick_b55i,thick_b56i,thick_b57i,
                thick_b58i,thick_b59i,thick_b60i,thick_b61i,thick_b62i,thick_b63i,thick_b64i,thick_b65i,thick_b66i,thick_b67i,thick_b68i)
thick_parcel <- cbind(thick_parcel,newthickbi)

newthickbs <- c(thick_b1s,thick_b2s,thick_b3s,thick_b4s,thick_b5s,thick_b6s,thick_b7s,thick_b8s,thick_b9s,thick_b10s,thick_b11s,thick_b12s,thick_b13s,thick_b14s,thick_b15s,
                thick_b16s,thick_b17s,thick_b18s,thick_b19s,thick_b20s,thick_b21s,thick_b22s,thick_b23s,thick_b24s,thick_b25s,thick_b26s,thick_b27s,thick_b28s,thick_b29s,
                thick_b30s,thick_b31s,thick_b32s,thick_b33s,thick_b34s,thick_b35s,thick_b36s,thick_b37s,thick_b38s,thick_b39s,thick_b40s,thick_b41s,thick_b42s,thick_b43s,
                thick_b44s,thick_b45s,thick_b46s,thick_b47s,thick_b48s,thick_b49s,thick_b50s,thick_b51s,thick_b52s,thick_b53s,thick_b54s,thick_b55s,thick_b56s,thick_b57s,
                thick_b58s,thick_b59s,thick_b60s,thick_b61s,thick_b62s,thick_b63s,thick_b64s,thick_b65s,thick_b66s,thick_b67s,thick_b68s)
thick_parcel <- cbind(thick_parcel,newthickbs)

newthicksei <- c(thick_se1i,thick_se2i,thick_se3i,thick_se4i,thick_se5i,thick_se6i,thick_se7i,thick_se8i,thick_se9i,thick_se10i,thick_se11i,thick_se12i,thick_se13i,
                 thick_se14i,thick_se15i,thick_se16i,thick_se17i,thick_se18i,thick_se19i,thick_se20i,thick_se21i,thick_se22i,thick_se23i,thick_se24i,thick_se25i,
                 thick_se26i,thick_se27i,thick_se28i,thick_se29i,thick_se30i,thick_se31i,thick_se32i,thick_se33i,thick_se34i,thick_se35i,thick_se36i,thick_se37i,
                 thick_se38i,thick_se39i,thick_se40i,thick_se41i,thick_se42i,thick_se43i,thick_se44i,thick_se45i,thick_se46i,thick_se47i,thick_se48i,thick_se49i,
                 thick_se50i,thick_se51i,thick_se52i,thick_se53i,thick_se54i,thick_se55i,thick_se56i,thick_se57i,thick_se58i,thick_se59i,thick_se60i,thick_se61i,
                 thick_se62i,thick_se63i,thick_se64i,thick_se65i,thick_se66i,thick_se67i,thick_se68i)
thick_parcel <- cbind(thick_parcel,newthicksei)

newthickses <- c(thick_se1s,thick_se2s,thick_se3s,thick_se4s,thick_se5s,thick_se6s,thick_se7s,thick_se8s,thick_se9s,thick_se10s,thick_se11s,thick_se12s,thick_se13s,
                 thick_se14s,thick_se15s,thick_se16s,thick_se17s,thick_se18s,thick_se19s,thick_se20s,thick_se21s,thick_se22s,thick_se23s,thick_se24s,thick_se25s,
                 thick_se26s,thick_se27s,thick_se28s,thick_se29s,thick_se30s,thick_se31s,thick_se32s,thick_se33s,thick_se34s,thick_se35s,thick_se36s,thick_se37s,
                 thick_se38s,thick_se39s,thick_se40s,thick_se41s,thick_se42s,thick_se43s,thick_se44s,thick_se45s,thick_se46s,thick_se47s,thick_se48s,thick_se49s,
                 thick_se50s,thick_se51s,thick_se52s,thick_se53s,thick_se54s,thick_se55s,thick_se56s,thick_se57s,thick_se58s,thick_se59s,thick_se60s,thick_se61s,
                 thick_se62s,thick_se63s,thick_se64s,thick_se65s,thick_se66s,thick_se67s,thick_se68s)
thick_parcel <- cbind(thick_parcel,newthickses)

newthickpi <- c(thick_p1i,thick_p2i,thick_p3i,thick_p4i,thick_p5i,thick_p6i,thick_p7i,thick_p8i,thick_p9i,thick_p10i,thick_p11i,thick_p12i,thick_p13i,thick_p14i,thick_p15i,
                thick_p16i,thick_p17i,thick_p18i,thick_p19i,thick_p20i,thick_p21i,thick_p22i,thick_p23i,thick_p24i,thick_p25i,thick_p26i,thick_p27i,thick_p28i,thick_p29i,
                thick_p30i,thick_p31i,thick_p32i,thick_p33i,thick_p34i,thick_p35i,thick_p36i,thick_p37i,thick_p38i,thick_p39i,thick_p40i,thick_p41i,thick_p42i,thick_p43i,
                thick_p44i,thick_p45i,thick_p46i,thick_p47i,thick_p48i,thick_p49i,thick_p50i,thick_p51i,thick_p52i,thick_p53i,thick_p54i,thick_p55i,thick_p56i,thick_p57i,
                thick_p58i,thick_p59i,thick_p60i,thick_p61i,thick_p62i,thick_p63i,thick_p64i,thick_p65i,thick_p66i,thick_p67i,thick_p68i)
thick_parcel <- cbind(thick_parcel,newthickpi)

newthickps <- c(thick_p1s,thick_p2s,thick_p3s,thick_p4s,thick_p5s,thick_p6s,thick_p7s,thick_p8s,thick_p9s,thick_p10s,thick_p11s,thick_p12s,thick_p13s,thick_p14s,thick_p15s,
                thick_p16s,thick_p17s,thick_p18s,thick_p19s,thick_p20s,thick_p21s,thick_p22s,thick_p23s,thick_p24s,thick_p25s,thick_p26s,thick_p27s,thick_p28s,thick_p29s,
                thick_p30s,thick_p31s,thick_p32s,thick_p33s,thick_p34s,thick_p35s,thick_p36s,thick_p37s,thick_p38s,thick_p39s,thick_p40s,thick_p41s,thick_p42s,thick_p43s,
                thick_p44s,thick_p45s,thick_p46s,thick_p47s,thick_p48s,thick_p49s,thick_p50s,thick_p51s,thick_p52s,thick_p53s,thick_p54s,thick_p55s,thick_p56s,thick_p57s,
                thick_p58s,thick_p59s,thick_p60s,thick_p61s,thick_p62s,thick_p63s,thick_p64s,thick_p65s,thick_p66s,thick_p67s,thick_p68s)
thick_parcel <- cbind(thick_parcel,newthickps)

names(thick_parcel) <- c('thick','banksstslh_thick_bi','cdacatelh_thick_bi','cdmdfrlh_thick_bi','cuneuslh_thick_bi','ehinallh_thick_bi','fusiformlh_thick_bi',
                         'ifpllh_thick_bi','iftmlh_thick_bi','ihcatelh_thick_bi','locclh_thick_bi','lobfrlh_thick_bi','linguallh_thick_bi',
                         'mobfrlh_thick_bi','mdtmlh_thick_bi','parahpallh_thick_bi','paracnlh_thick_bi','parsopclh_thick_bi','parsobislh_thick_bi',
                         'parstgrislh_thick_bi','pericclh_thick_bi','postcnlh_thick_bi','ptcatelh_thick_bi','precnlh_thick_bi','pclh_thick_bi',
                         'rracatelh_thick_bi','rrmdfrlh_thick_bi','sufrlh_thick_bi','supllh_thick_bi','sutmlh_thick_bi','smlh_thick_bi','frpolelh_thick_bi',
                         'tmpolelh_thick_bi','trvtmlh_thick_bi','insulalh_thick_bi','banksstsrh_thick_bi','cdacaterh_thick_bi','cdmdfrrh_thick_bi',
                         'cuneusrh_thick_bi','ehinalrh_thick_bi','fusiformrh_thick_bi','ifplrh_thick_bi','iftmrh_thick_bi','ihcaterh_thick_bi',
                         'loccrh_thick_bi','lobfrrh_thick_bi','lingualrh_thick_bi','mobfrrh_thick_bi','mdtmrh_thick_bi','parahpalrh_thick_bi',
                         'paracnrh_thick_bi','parsopcrh_thick_bi','parsobisrh_thick_bi','parstgrisrh_thick_bi','periccrh_thick_bi',
                         'postcnrh_thick_bi','ptcaterh_thick_bi','precnrh_thick_bi','pcrh_thick_bi','rracaterh_thick_bi','rrmdfrrh_thick_bi',
                         'sufrrh_thick_bi','suplrh_thick_bi','sutmrh_thick_bi','smrh_thick_bi','frpolerh_thick_bi',
                         'tmpolerh_thick_bi','trvtmrh_thick_bi','insularh_thick_bi',
                         'banksstslh_thick_bs','cdacatelh_thick_bs','cdmdfrlh_thick_bs','cuneuslh_thick_bs','ehinallh_thick_bs','fusiformlh_thick_bs',
                         'ifpllh_thick_bs','iftmlh_thick_bs','ihcatelh_thick_bs','locclh_thick_bs','lobfrlh_thick_bs','linguallh_thick_bs',
                         'mobfrlh_thick_bs','mdtmlh_thick_bs','parahpallh_thick_bs','paracnlh_thick_bs','parsopclh_thick_bs','parsobislh_thick_bs',
                         'parstgrislh_thick_bs','pericclh_thick_bs','postcnlh_thick_bs','ptcatelh_thick_bs','precnlh_thick_bs','pclh_thick_bs',
                         'rracatelh_thick_bs','rrmdfrlh_thick_bs','sufrlh_thick_bs','supllh_thick_bs','sutmlh_thick_bs','smlh_thick_bs','frpolelh_thick_bs',
                         'tmpolelh_thick_bs','trvtmlh_thick_bs','insulalh_thick_bs','banksstsrh_thick_bs','cdacaterh_thick_bs','cdmdfrrh_thick_bs',
                         'cuneusrh_thick_bs','ehinalrh_thick_bs','fusiformrh_thick_bs','ifplrh_thick_bs','iftmrh_thick_bs','ihcaterh_thick_bs',
                         'loccrh_thick_bs','lobfrrh_thick_bs','lingualrh_thick_bs','mobfrrh_thick_bs','mdtmrh_thick_bs','parahpalrh_thick_bs',
                         'paracnrh_thick_bs','parsopcrh_thick_bs','parsobisrh_thick_bs','parstgrisrh_thick_bs','periccrh_thick_bs',
                         'postcnrh_thick_bs','ptcaterh_thick_bs','precnrh_thick_bs','pcrh_thick_bs','rracaterh_thick_bs','rrmdfrrh_thick_bs',
                         'sufrrh_thick_bs','suplrh_thick_bs','sutmrh_thick_bs','smrh_thick_bs','frpolerh_thick_bs',
                         'tmpolerh_thick_bs','trvtmrh_thick_bs','insularh_thick_bs',
                         'banksstslh_thick_sei','cdacatelh_thick_sei','cdmdfrlh_thick_sei','cuneuslh_thick_sei','ehinallh_thick_sei','fusiformlh_thick_sei',
                         'ifpllh_thick_sei','iftmlh_thick_sei','ihcatelh_thick_sei','locclh_thick_sei','lobfrlh_thick_sei','linguallh_thick_sei',
                         'mobfrlh_thick_sei','mdtmlh_thick_sei','parahpallh_thick_sei','paracnlh_thick_sei','parsopclh_thick_sei',
                         'parsobislh_thick_sei','parstgrislh_thick_sei','pericclh_thick_sei','postcnlh_thick_sei','ptcatelh_thick_sei',
                         'precnlh_thick_sei','pclh_thick_sei','rracatelh_thick_sei','rrmdfrlh_thick_sei','sufrlh_thick_sei','supllh_thick_sei',
                         'sutmlh_thick_sei','smlh_thick_sei','frpolelh_thick_sei','tmpolelh_thick_sei','trvtmlh_thick_sei','insulalh_thick_sei',
                         'banksstsrh_thick_sei','cdacaterh_thick_sei','cdmdfrrh_thick_sei','cuneusrh_thick_sei','ehinalrh_thick_sei',
                         'fusiformrh_thick_sei','ifplrh_thick_sei','iftmrh_thick_sei','ihcaterh_thick_sei','loccrh_thick_sei','lobfrrh_thick_sei',
                         'lingualrh_thick_sei','mobfrrh_thick_sei','mdtmrh_thick_sei','parahpalrh_thick_sei','paracnrh_thick_sei',
                         'parsopcrh_thick_sei','parsobisrh_thick_sei','parstgrisrh_thick_sei','periccrh_thick_sei','postcnrh_thick_sei',
                         'ptcaterh_thick_sei','precnrh_thick_sei','pcrh_thick_sei','rracaterh_thick_sei','rrmdfrrh_thick_sei',
                         'sufrrh_thick_sei','suplrh_thick_sei','sutmrh_thick_sei','smrh_thick_sei','frpolerh_thick_sei',
                         'tmpolerh_thick_sei','trvtmrh_thick_sei','insularh_thick_sei',
                         'banksstslh_thick_ses','cdacatelh_thick_ses','cdmdfrlh_thick_ses','cuneuslh_thick_ses','ehinallh_thick_ses','fusiformlh_thick_ses',
                         'ifpllh_thick_ses','iftmlh_thick_ses','ihcatelh_thick_ses','locclh_thick_ses','lobfrlh_thick_ses','linguallh_thick_ses',
                         'mobfrlh_thick_ses','mdtmlh_thick_ses','parahpallh_thick_ses','paracnlh_thick_ses','parsopclh_thick_ses',
                         'parsobislh_thick_ses','parstgrislh_thick_ses','pericclh_thick_ses','postcnlh_thick_ses','ptcatelh_thick_ses',
                         'precnlh_thick_ses','pclh_thick_ses','rracatelh_thick_ses','rrmdfrlh_thick_ses','sufrlh_thick_ses','supllh_thick_ses',
                         'sutmlh_thick_ses','smlh_thick_ses','frpolelh_thick_ses','tmpolelh_thick_ses','trvtmlh_thick_ses','insulalh_thick_ses',
                         'banksstsrh_thick_ses','cdacaterh_thick_ses','cdmdfrrh_thick_ses','cuneusrh_thick_ses','ehinalrh_thick_ses',
                         'fusiformrh_thick_ses','ifplrh_thick_ses','iftmrh_thick_ses','ihcaterh_thick_ses','loccrh_thick_ses','lobfrrh_thick_ses',
                         'lingualrh_thick_ses','mobfrrh_thick_ses','mdtmrh_thick_ses','parahpalrh_thick_ses','paracnrh_thick_ses',
                         'parsopcrh_thick_ses','parsobisrh_thick_ses','parstgrisrh_thick_ses','periccrh_thick_ses','postcnrh_thick_ses',
                         'ptcaterh_thick_ses','precnrh_thick_ses','pcrh_thick_ses','rracaterh_thick_ses','rrmdfrrh_thick_ses',
                         'sufrrh_thick_ses','suplrh_thick_ses','sutmrh_thick_ses','smrh_thick_ses','frpolerh_thick_ses',
                         'tmpolerh_thick_ses','trvtmrh_thick_ses','insularh_thick_ses',
                         'banksstslh_thick_pi','cdacatelh_thick_pi','cdmdfrlh_thick_pi','cuneuslh_thick_pi','ehinallh_thick_pi','fusiformlh_thick_pi','ifpllh_thick_pi','iftmlh_thick_pi',
                         'ihcatelh_thick_pi','locclh_thick_pi','lobfrlh_thick_pi','linguallh_thick_pi','mobfrlh_thick_pi','mdtmlh_thick_pi',
                         'parahpallh_thick_pi','paracnlh_thick_pi','parsopclh_thick_pi','parsobislh_thick_pi','parstgrislh_thick_pi',
                         'pericclh_thick_pi','postcnlh_thick_pi','ptcatelh_thick_pi','precnlh_thick_pi','pclh_thick_pi','rracatelh_thick_pi',
                         'rrmdfrlh_thick_pi','sufrlh_thick_pi','supllh_thick_pi','sutmlh_thick_pi','smlh_thick_pi','frpolelh_thick_pi',
                         'tmpolelh_thick_pi','trvtmlh_thick_pi','insulalh_thick_pi','banksstsrh_thick_pi','cdacaterh_thick_pi','cdmdfrrh_thick_pi',
                         'cuneusrh_thick_pi','ehinalrh_thick_pi','fusiformrh_thick_pi','ifplrh_thick_pi','iftmrh_thick_pi','ihcaterh_thick_pi',
                         'loccrh_thick_pi','lobfrrh_thick_pi','lingualrh_thick_pi','mobfrrh_thick_pi','mdtmrh_thick_pi','parahpalrh_thick_pi',
                         'paracnrh_thick_pi','parsopcrh_thick_pi','parsobisrh_thick_pi','parstgrisrh_thick_pi','periccrh_thick_pi',
                         'postcnrh_thick_pi','ptcaterh_thick_pi','precnrh_thick_pi','pcrh_thick_pi','rracaterh_thick_pi','rrmdfrrh_thick_pi',
                         'sufrrh_thick_pi','suplrh_thick_pi','sutmrh_thick_pi','smrh_thick_pi','frpolerh_thick_pi',
                         'tmpolerh_thick_pi','trvtmrh_thick_pi','insularh_thick_pi',
                         'banksstslh_thick_ps','cdacatelh_thick_ps','cdmdfrlh_thick_ps','cuneuslh_thick_ps','ehinallh_thick_ps','fusiformlh_thick_ps','ifpllh_thick_ps','iftmlh_thick_ps',
                         'ihcatelh_thick_ps','locclh_thick_ps','lobfrlh_thick_ps','linguallh_thick_ps','mobfrlh_thick_ps','mdtmlh_thick_ps',
                         'parahpallh_thick_ps','paracnlh_thick_ps','parsopclh_thick_ps','parsobislh_thick_ps','parstgrislh_thick_ps',
                         'pericclh_thick_ps','postcnlh_thick_ps','ptcatelh_thick_ps','precnlh_thick_ps','pclh_thick_ps','rracatelh_thick_ps',
                         'rrmdfrlh_thick_ps','sufrlh_thick_ps','supllh_thick_ps','sutmlh_thick_ps','smlh_thick_ps','frpolelh_thick_ps',
                         'tmpolelh_thick_ps','trvtmlh_thick_ps','insulalh_thick_ps','banksstsrh_thick_ps','cdacaterh_thick_ps','cdmdfrrh_thick_ps',
                         'cuneusrh_thick_ps','ehinalrh_thick_ps','fusiformrh_thick_ps','ifplrh_thick_ps','iftmrh_thick_ps','ihcaterh_thick_ps',
                         'loccrh_thick_ps','lobfrrh_thick_ps','lingualrh_thick_ps','mobfrrh_thick_ps','mdtmrh_thick_ps','parahpalrh_thick_ps',
                         'paracnrh_thick_ps','parsopcrh_thick_ps','parsobisrh_thick_ps','parstgrisrh_thick_ps','periccrh_thick_ps',
                         'postcnrh_thick_ps','ptcaterh_thick_ps','precnrh_thick_ps','pcrh_thick_ps','rracaterh_thick_ps','rrmdfrrh_thick_ps',
                         'sufrrh_thick_ps','suplrh_thick_ps','sutmrh_thick_ps','smrh_thick_ps','frpolerh_thick_ps',
                         'tmpolerh_thick_ps','trvtmrh_thick_ps','insularh_thick_ps')

#calculate 95% CIs and create lower and upper bound variables 
thick_parcel$ci_lower_banksstslh_thicks <- thick_parcel$banksstslh_thick_bs - 1.96*thick_parcel$banksstslh_thick_ses 
thick_parcel$ci_upper_banksstslh_thicks <- thick_parcel$banksstslh_thick_bs + 1.96*thick_parcel$banksstslh_thick_ses
thick_parcel$ci_lower_cdacatelh_thicks <- thick_parcel$cdacatelh_thick_bs - 1.96*thick_parcel$cdacatelh_thick_ses 
thick_parcel$ci_upper_cdacatelh_thicks <- thick_parcel$cdacatelh_thick_bs + 1.96*thick_parcel$cdacatelh_thick_ses
thick_parcel$ci_lower_cdmdfrlh_thicks <- thick_parcel$cdmdfrlh_thick_bs - 1.96*thick_parcel$cdmdfrlh_thick_ses 
thick_parcel$ci_upper_cdmdfrlh_thicks <- thick_parcel$cdmdfrlh_thick_bs + 1.96*thick_parcel$cdmdfrlh_thick_ses
thick_parcel$ci_lower_cuneuslh_thicks <- thick_parcel$cuneuslh_thick_bs - 1.96*thick_parcel$cuneuslh_thick_ses 
thick_parcel$ci_upper_cuneuslh_thicks <- thick_parcel$cuneuslh_thick_bs + 1.96*thick_parcel$cuneuslh_thick_ses
thick_parcel$ci_lower_ehinallh_thicks <- thick_parcel$ehinallh_thick_bs - 1.96*thick_parcel$ehinallh_thick_ses 
thick_parcel$ci_upper_ehinallh_thicks <- thick_parcel$ehinallh_thick_bs + 1.96*thick_parcel$ehinallh_thick_ses
thick_parcel$ci_lower_fusiformlh_thicks <- thick_parcel$fusiformlh_thick_bs - 1.96*thick_parcel$fusiformlh_thick_ses 
thick_parcel$ci_upper_fusiformlh_thicks <- thick_parcel$fusiformlh_thick_bs + 1.96*thick_parcel$fusiformlh_thick_ses
thick_parcel$ci_lower_ifpllh_thicks <- thick_parcel$ifpllh_thick_bs - 1.96*thick_parcel$ifpllh_thick_ses 
thick_parcel$ci_upper_ifpllh_thicks <- thick_parcel$ifpllh_thick_bs + 1.96*thick_parcel$ifpllh_thick_ses
thick_parcel$ci_lower_iftmlh_thicks <- thick_parcel$iftmlh_thick_bs - 1.96*thick_parcel$iftmlh_thick_ses
thick_parcel$ci_upper_iftmlh_thicks <- thick_parcel$iftmlh_thick_bs + 1.96*thick_parcel$iftmlh_thick_ses
thick_parcel$ci_lower_ihcatelh_thicks <- thick_parcel$ihcatelh_thick_bs - 1.96*thick_parcel$ihcatelh_thick_ses 
thick_parcel$ci_upper_ihcatelh_thicks <- thick_parcel$ihcatelh_thick_bs + 1.96*thick_parcel$ihcatelh_thick_ses
thick_parcel$ci_lower_locclh_thicks <- thick_parcel$locclh_thick_bs - 1.96*thick_parcel$locclh_thick_ses
thick_parcel$ci_upper_locclh_thicks <- thick_parcel$locclh_thick_bs + 1.96*thick_parcel$locclh_thick_ses
thick_parcel$ci_lower_lobfrlh_thicks <- thick_parcel$lobfrlh_thick_bs - 1.96*thick_parcel$lobfrlh_thick_ses 
thick_parcel$ci_upper_lobfrlh_thicks <- thick_parcel$lobfrlh_thick_bs + 1.96*thick_parcel$lobfrlh_thick_ses
thick_parcel$ci_lower_linguallh_thicks <- thick_parcel$linguallh_thick_bs - 1.96*thick_parcel$linguallh_thick_ses 
thick_parcel$ci_upper_linguallh_thicks <- thick_parcel$linguallh_thick_bs + 1.96*thick_parcel$linguallh_thick_ses
thick_parcel$ci_lower_mobfrlh_thicks <- thick_parcel$mobfrlh_thick_bs - 1.96*thick_parcel$mobfrlh_thick_ses 
thick_parcel$ci_upper_mobfrlh_thicks <- thick_parcel$mobfrlh_thick_bs + 1.96*thick_parcel$mobfrlh_thick_ses
thick_parcel$ci_lower_mdtmlh_thicks <- thick_parcel$mdtmlh_thick_bs - 1.96*thick_parcel$mdtmlh_thick_ses 
thick_parcel$ci_upper_mdtmlh_thicks <- thick_parcel$mdtmlh_thick_bs + 1.96*thick_parcel$mdtmlh_thick_ses
thick_parcel$ci_lower_parahpallh_thicks <- thick_parcel$parahpallh_thick_bs - 1.96*thick_parcel$parahpallh_thick_ses 
thick_parcel$ci_upper_parahpallh_thicks <- thick_parcel$parahpallh_thick_bs + 1.96*thick_parcel$parahpallh_thick_ses
thick_parcel$ci_lower_paracnlh_thicks <- thick_parcel$paracnlh_thick_bs - 1.96*thick_parcel$paracnlh_thick_ses 
thick_parcel$ci_upper_paracnlh_thicks <- thick_parcel$paracnlh_thick_bs + 1.96*thick_parcel$paracnlh_thick_ses
thick_parcel$ci_lower_parsopclh_thicks <- thick_parcel$parsopclh_thick_bs - 1.96*thick_parcel$parsopclh_thick_ses 
thick_parcel$ci_upper_parsopclh_thicks <- thick_parcel$parsopclh_thick_bs + 1.96*thick_parcel$parsopclh_thick_ses
thick_parcel$ci_lower_parsobislh_thicks <- thick_parcel$parsobislh_thick_bs - 1.96*thick_parcel$parsobislh_thick_ses 
thick_parcel$ci_upper_parsobislh_thicks <- thick_parcel$parsobislh_thick_bs + 1.96*thick_parcel$parsobislh_thick_ses
thick_parcel$ci_lower_parstgrislh_thicks <- thick_parcel$parstgrislh_thick_bs - 1.96*thick_parcel$parstgrislh_thick_ses 
thick_parcel$ci_upper_parstgrislh_thicks <- thick_parcel$parstgrislh_thick_bs + 1.96*thick_parcel$parstgrislh_thick_ses
thick_parcel$ci_lower_pericclh_thicks <- thick_parcel$pericclh_thick_bs - 1.96*thick_parcel$pericclh_thick_ses 
thick_parcel$ci_upper_pericclh_thicks <- thick_parcel$pericclh_thick_bs + 1.96*thick_parcel$pericclh_thick_ses
thick_parcel$ci_lower_postcnlh_thicks <- thick_parcel$postcnlh_thick_bs - 1.96*thick_parcel$postcnlh_thick_ses 
thick_parcel$ci_upper_postcnlh_thicks <- thick_parcel$postcnlh_thick_bs + 1.96*thick_parcel$postcnlh_thick_ses
thick_parcel$ci_lower_ptcatelh_thicks <- thick_parcel$ptcatelh_thick_bs - 1.96*thick_parcel$ptcatelh_thick_ses 
thick_parcel$ci_upper_ptcatelh_thicks <- thick_parcel$ptcatelh_thick_bs + 1.96*thick_parcel$ptcatelh_thick_ses
thick_parcel$ci_lower_precnlh_thicks <- thick_parcel$precnlh_thick_bs - 1.96*thick_parcel$precnlh_thick_ses 
thick_parcel$ci_upper_precnlh_thicks <- thick_parcel$precnlh_thick_bs + 1.96*thick_parcel$precnlh_thick_ses
thick_parcel$ci_lower_pclh_thicks <- thick_parcel$pclh_thick_bs - 1.96*thick_parcel$pclh_thick_ses 
thick_parcel$ci_upper_pclh_thicks <- thick_parcel$pclh_thick_bs + 1.96*thick_parcel$pclh_thick_ses
thick_parcel$ci_lower_rracatelh_thicks <- thick_parcel$rracatelh_thick_bs - 1.96*thick_parcel$rracatelh_thick_ses 
thick_parcel$ci_upper_rracatelh_thicks <- thick_parcel$rracatelh_thick_bs + 1.96*thick_parcel$rracatelh_thick_ses
thick_parcel$ci_lower_rrmdfrlh_thicks <- thick_parcel$rrmdfrlh_thick_bs - 1.96*thick_parcel$rrmdfrlh_thick_ses 
thick_parcel$ci_upper_rrmdfrlh_thicks <- thick_parcel$rrmdfrlh_thick_bs + 1.96*thick_parcel$rrmdfrlh_thick_ses
thick_parcel$ci_lower_sufrlh_thicks <- thick_parcel$sufrlh_thick_bs - 1.96*thick_parcel$sufrlh_thick_ses 
thick_parcel$ci_upper_sufrlh_thicks <- thick_parcel$sufrlh_thick_bs + 1.96*thick_parcel$sufrlh_thick_ses
thick_parcel$ci_lower_supllh_thicks <- thick_parcel$supllh_thick_bs - 1.96*thick_parcel$supllh_thick_ses 
thick_parcel$ci_upper_supllh_thicks <- thick_parcel$supllh_thick_bs + 1.96*thick_parcel$supllh_thick_ses
thick_parcel$ci_lower_sutmlh_thicks <- thick_parcel$sutmlh_thick_bs - 1.96*thick_parcel$sutmlh_thick_ses 
thick_parcel$ci_upper_sutmlh_thicks <- thick_parcel$sutmlh_thick_bs + 1.96*thick_parcel$sutmlh_thick_ses
thick_parcel$ci_lower_smlh_thicks <- thick_parcel$smlh_thick_bs - 1.96*thick_parcel$smlh_thick_ses 
thick_parcel$ci_upper_smlh_thicks <- thick_parcel$smlh_thick_bs + 1.96*thick_parcel$smlh_thick_ses
thick_parcel$ci_lower_frpolelh_thicks <- thick_parcel$frpolelh_thick_bs - 1.96*thick_parcel$frpolelh_thick_ses 
thick_parcel$ci_upper_frpolelh_thicks <- thick_parcel$frpolelh_thick_bs + 1.96*thick_parcel$frpolelh_thick_ses
thick_parcel$ci_lower_tmpolelh_thicks <- thick_parcel$tmpolelh_thick_bs - 1.96*thick_parcel$tmpolelh_thick_ses 
thick_parcel$ci_upper_tmpolelh_thicks <- thick_parcel$tmpolelh_thick_bs + 1.96*thick_parcel$tmpolelh_thick_ses
thick_parcel$ci_lower_trvtmlh_thicks <- thick_parcel$trvtmlh_thick_bs - 1.96*thick_parcel$trvtmlh_thick_ses 
thick_parcel$ci_upper_trvtmlh_thicks <- thick_parcel$trvtmlh_thick_bs + 1.96*thick_parcel$trvtmlh_thick_ses
thick_parcel$ci_lower_insulalh_thicks <- thick_parcel$insulalh_thick_bs - 1.96*thick_parcel$insulalh_thick_ses 
thick_parcel$ci_upper_insulalh_thicks <- thick_parcel$insulalh_thick_bs + 1.96*thick_parcel$insulalh_thick_ses

thick_parcel$ci_lower_banksstsrh_thicks <- thick_parcel$banksstsrh_thick_bs - 1.96*thick_parcel$banksstsrh_thick_ses 
thick_parcel$ci_upper_banksstsrh_thicks <- thick_parcel$banksstsrh_thick_bs + 1.96*thick_parcel$banksstsrh_thick_ses
thick_parcel$ci_lower_cdacaterh_thicks <- thick_parcel$cdacaterh_thick_bs - 1.96*thick_parcel$cdacaterh_thick_ses 
thick_parcel$ci_upper_cdacaterh_thicks <- thick_parcel$cdacaterh_thick_bs + 1.96*thick_parcel$cdacaterh_thick_ses
thick_parcel$ci_lower_cdmdfrrh_thicks <- thick_parcel$cdmdfrrh_thick_bs - 1.96*thick_parcel$cdmdfrrh_thick_ses 
thick_parcel$ci_upper_cdmdfrrh_thicks <- thick_parcel$cdmdfrrh_thick_bs + 1.96*thick_parcel$cdmdfrrh_thick_ses
thick_parcel$ci_lower_cuneusrh_thicks <- thick_parcel$cuneusrh_thick_bs - 1.96*thick_parcel$cuneusrh_thick_ses 
thick_parcel$ci_upper_cuneusrh_thicks <- thick_parcel$cuneusrh_thick_bs + 1.96*thick_parcel$cuneusrh_thick_ses
thick_parcel$ci_lower_ehinalrh_thicks <- thick_parcel$ehinalrh_thick_bs - 1.96*thick_parcel$ehinalrh_thick_ses 
thick_parcel$ci_upper_ehinalrh_thicks <- thick_parcel$ehinalrh_thick_bs + 1.96*thick_parcel$ehinalrh_thick_ses
thick_parcel$ci_lower_fusiformrh_thicks <- thick_parcel$fusiformrh_thick_bs - 1.96*thick_parcel$fusiformrh_thick_ses 
thick_parcel$ci_upper_fusiformrh_thicks <- thick_parcel$fusiformrh_thick_bs + 1.96*thick_parcel$fusiformrh_thick_ses
thick_parcel$ci_lower_ifplrh_thicks <- thick_parcel$ifplrh_thick_bs - 1.96*thick_parcel$ifplrh_thick_ses 
thick_parcel$ci_upper_ifplrh_thicks <- thick_parcel$ifplrh_thick_bs + 1.96*thick_parcel$ifplrh_thick_ses
thick_parcel$ci_lower_iftmrh_thicks <- thick_parcel$iftmrh_thick_bs - 1.96*thick_parcel$iftmrh_thick_ses 
thick_parcel$ci_upper_iftmrh_thicks <- thick_parcel$iftmrh_thick_bs + 1.96*thick_parcel$iftmrh_thick_ses
thick_parcel$ci_lower_ihcaterh_thicks <- thick_parcel$ihcaterh_thick_bs - 1.96*thick_parcel$ihcaterh_thick_ses 
thick_parcel$ci_upper_ihcaterh_thicks <- thick_parcel$ihcaterh_thick_bs + 1.96*thick_parcel$ihcaterh_thick_ses
thick_parcel$ci_lower_loccrh_thicks <- thick_parcel$loccrh_thick_bs - 1.96*thick_parcel$loccrh_thick_ses 
thick_parcel$ci_upper_loccrh_thicks <- thick_parcel$loccrh_thick_bs + 1.96*thick_parcel$loccrh_thick_ses
thick_parcel$ci_lower_lobfrrh_thicks <- thick_parcel$lobfrrh_thick_bs - 1.96*thick_parcel$lobfrrh_thick_ses 
thick_parcel$ci_upper_lobfrrh_thicks <- thick_parcel$lobfrrh_thick_bs + 1.96*thick_parcel$lobfrrh_thick_ses
thick_parcel$ci_lower_lingualrh_thicks <- thick_parcel$lingualrh_thick_bs - 1.96*thick_parcel$lingualrh_thick_ses 
thick_parcel$ci_upper_lingualrh_thicks <- thick_parcel$lingualrh_thick_bs + 1.96*thick_parcel$lingualrh_thick_ses
thick_parcel$ci_lower_mobfrrh_thicks <- thick_parcel$mobfrrh_thick_bs - 1.96*thick_parcel$mobfrrh_thick_ses 
thick_parcel$ci_upper_mobfrrh_thicks <- thick_parcel$mobfrrh_thick_bs + 1.96*thick_parcel$mobfrrh_thick_ses
thick_parcel$ci_lower_mdtmrh_thicks <- thick_parcel$mdtmrh_thick_bs - 1.96*thick_parcel$mdtmrh_thick_ses 
thick_parcel$ci_upper_mdtmrh_thicks <- thick_parcel$mdtmrh_thick_bs + 1.96*thick_parcel$mdtmrh_thick_ses
thick_parcel$ci_lower_parahpalrh_thicks <- thick_parcel$parahpalrh_thick_bs - 1.96*thick_parcel$parahpalrh_thick_ses 
thick_parcel$ci_upper_parahpalrh_thicks <- thick_parcel$parahpalrh_thick_bs + 1.96*thick_parcel$parahpalrh_thick_ses
thick_parcel$ci_lower_paracnrh_thicks <- thick_parcel$paracnrh_thick_bs - 1.96*thick_parcel$paracnrh_thick_ses 
thick_parcel$ci_upper_paracnrh_thicks <- thick_parcel$paracnrh_thick_bs + 1.96*thick_parcel$paracnrh_thick_ses
thick_parcel$ci_lower_parsopcrh_thicks <- thick_parcel$parsopcrh_thick_bs - 1.96*thick_parcel$parsopcrh_thick_ses 
thick_parcel$ci_upper_parsopcrh_thicks <- thick_parcel$parsopcrh_thick_bs + 1.96*thick_parcel$parsopcrh_thick_ses
thick_parcel$ci_lower_parsobisrh_thicks <- thick_parcel$parsobisrh_thick_bs - 1.96*thick_parcel$parsobisrh_thick_ses 
thick_parcel$ci_upper_parsobisrh_thicks <- thick_parcel$parsobisrh_thick_bs + 1.96*thick_parcel$parsobisrh_thick_ses
thick_parcel$ci_lower_parstgrisrh_thicks <- thick_parcel$parstgrisrh_thick_bs - 1.96*thick_parcel$parstgrisrh_thick_ses 
thick_parcel$ci_upper_parstgrisrh_thicks <- thick_parcel$parstgrisrh_thick_bs + 1.96*thick_parcel$parstgrisrh_thick_ses
thick_parcel$ci_lower_periccrh_thicks <- thick_parcel$periccrh_thick_bs - 1.96*thick_parcel$periccrh_thick_ses 
thick_parcel$ci_upper_periccrh_thicks <- thick_parcel$periccrh_thick_bs + 1.96*thick_parcel$periccrh_thick_ses
thick_parcel$ci_lower_postcnrh_thicks <- thick_parcel$postcnrh_thick_bs - 1.96*thick_parcel$postcnrh_thick_ses 
thick_parcel$ci_upper_postcnrh_thicks <- thick_parcel$postcnrh_thick_bs + 1.96*thick_parcel$postcnrh_thick_ses
thick_parcel$ci_lower_ptcaterh_thicks <- thick_parcel$ptcaterh_thick_bs - 1.96*thick_parcel$ptcaterh_thick_ses 
thick_parcel$ci_upper_ptcaterh_thicks <- thick_parcel$ptcaterh_thick_bs + 1.96*thick_parcel$ptcaterh_thick_ses
thick_parcel$ci_lower_precnrh_thicks <- thick_parcel$precnrh_thick_bs - 1.96*thick_parcel$precnrh_thick_ses 
thick_parcel$ci_upper_precnrh_thicks <- thick_parcel$precnrh_thick_bs + 1.96*thick_parcel$precnrh_thick_ses
thick_parcel$ci_lower_pcrh_thicks <- thick_parcel$pcrh_thick_bs - 1.96*thick_parcel$pcrh_thick_ses 
thick_parcel$ci_upper_pcrh_thicks <- thick_parcel$pcrh_thick_bs + 1.96*thick_parcel$pcrh_thick_ses
thick_parcel$ci_lower_rracaterh_thicks <- thick_parcel$rracaterh_thick_bs - 1.96*thick_parcel$rracaterh_thick_ses 
thick_parcel$ci_upper_rracaterh_thicks <- thick_parcel$rracaterh_thick_bs + 1.96*thick_parcel$rracaterh_thick_ses
thick_parcel$ci_lower_rrmdfrrh_thicks <- thick_parcel$rrmdfrrh_thick_bs - 1.96*thick_parcel$rrmdfrrh_thick_ses 
thick_parcel$ci_upper_rrmdfrrh_thicks <- thick_parcel$rrmdfrrh_thick_bs + 1.96*thick_parcel$rrmdfrrh_thick_ses
thick_parcel$ci_lower_sufrrh_thicks <- thick_parcel$sufrrh_thick_bs - 1.96*thick_parcel$sufrrh_thick_ses 
thick_parcel$ci_upper_sufrrh_thicks <- thick_parcel$sufrrh_thick_bs + 1.96*thick_parcel$sufrrh_thick_ses
thick_parcel$ci_lower_suplrh_thicks <- thick_parcel$suplrh_thick_bs - 1.96*thick_parcel$suplrh_thick_ses 
thick_parcel$ci_upper_suplrh_thicks <- thick_parcel$suplrh_thick_bs + 1.96*thick_parcel$suplrh_thick_ses
thick_parcel$ci_lower_sutmrh_thicks <- thick_parcel$sutmrh_thick_bs - 1.96*thick_parcel$sutmrh_thick_ses 
thick_parcel$ci_upper_sutmrh_thicks <- thick_parcel$sutmrh_thick_bs + 1.96*thick_parcel$sutmrh_thick_ses
thick_parcel$ci_lower_smrh_thicks <- thick_parcel$smrh_thick_bs - 1.96*thick_parcel$smrh_thick_ses 
thick_parcel$ci_upper_smrh_thicks <- thick_parcel$smrh_thick_bs + 1.96*thick_parcel$smrh_thick_ses
thick_parcel$ci_lower_frpolerh_thicks <- thick_parcel$frpolerh_thick_bs - 1.96*thick_parcel$frpolerh_thick_ses 
thick_parcel$ci_upper_frpolerh_thicks <- thick_parcel$frpolerh_thick_bs + 1.96*thick_parcel$frpolerh_thick_ses
thick_parcel$ci_lower_tmpolerh_thicks <- thick_parcel$tmpolerh_thick_bs - 1.96*thick_parcel$tmpolerh_thick_ses 
thick_parcel$ci_upper_tmpolerh_thicks <- thick_parcel$tmpolerh_thick_bs + 1.96*thick_parcel$tmpolerh_thick_ses
thick_parcel$ci_lower_trvtmrh_thicks <- thick_parcel$trvtmrh_thick_bs - 1.96*thick_parcel$trvtmrh_thick_ses 
thick_parcel$ci_upper_trvtmrh_thicks <- thick_parcel$trvtmrh_thick_bs + 1.96*thick_parcel$trvtmrh_thick_ses
thick_parcel$ci_lower_insularh_thicks <- thick_parcel$insularh_thick_bs - 1.96*thick_parcel$insularh_thick_ses 
thick_parcel$ci_upper_insularh_thicks <- thick_parcel$insularh_thick_bs + 1.96*thick_parcel$insularh_thick_ses

write.csv(thick_parcel, "Parcel-Wise Cortical Thickness Standardized_FINAL.csv") #write csv file

interact_plot(thick55, pred = wave, modx = subj_postcnrh_thick, 
              modx.values=c(-0.14,0,0.14),
              x.label = "Wave", y.label = "INT Factor Scores")

#Parcel-wise subcortical volume analyses

sub1 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_crbcortexlh_vol_cent + scale(subj_crbcortexlh_vol) + site_crbcortexlh_vol_cent*wave + scale(subj_crbcortexlh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b1i <- as.data.frame(c(summary(sub1)$coefficients[14,1]))
sub_se1i <- as.data.frame(c(summary(sub1)$coefficients[14,2]))
sub_p1i <- as.data.frame(c(summary(sub1)$coefficients[14,5]))
sub_b1s <- as.data.frame(c(summary(sub1)$coefficients[16,1]))
sub_se1s <- as.data.frame(c(summary(sub1)$coefficients[16,2]))
sub_p1s <- as.data.frame(c(summary(sub1)$coefficients[16,5]))


sub2 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_tplh_vol_cent + scale(subj_tplh_vol) + site_tplh_vol_cent*wave + scale(subj_tplh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b2i <- as.data.frame(c(summary(sub2)$coefficients[14,1]))
sub_se2i <- as.data.frame(c(summary(sub2)$coefficients[14,2]))
sub_p2i <- as.data.frame(c(summary(sub2)$coefficients[14,5]))
sub_b2s <- as.data.frame(c(summary(sub2)$coefficients[16,1]))
sub_se2s <- as.data.frame(c(summary(sub2)$coefficients[16,2]))
sub_p2s <- as.data.frame(c(summary(sub2)$coefficients[16,5]))

sub3 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_caudatelh_vol_cent + scale(subj_caudatelh_vol) + site_caudatelh_vol_cent*wave + scale(subj_caudatelh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b3i <- as.data.frame(c(summary(sub3)$coefficients[14,1]))
sub_se3i <- as.data.frame(c(summary(sub3)$coefficients[14,2]))
sub_p3i <- as.data.frame(c(summary(sub3)$coefficients[14,5]))
sub_b3s <- as.data.frame(c(summary(sub3)$coefficients[16,1]))
sub_se3s <- as.data.frame(c(summary(sub3)$coefficients[16,2]))
sub_p3s <- as.data.frame(c(summary(sub3)$coefficients[16,5]))

sub4 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_putamenlh_vol_cent + scale(subj_putamenlh_vol) + site_putamenlh_vol_cent*wave + scale(subj_putamenlh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b4i <- as.data.frame(c(summary(sub4)$coefficients[14,1]))
sub_se4i <- as.data.frame(c(summary(sub4)$coefficients[14,2]))
sub_p4i <- as.data.frame(c(summary(sub4)$coefficients[14,5]))
sub_b4s <- as.data.frame(c(summary(sub4)$coefficients[16,1]))
sub_se4s <- as.data.frame(c(summary(sub4)$coefficients[16,2]))
sub_p4s <- as.data.frame(c(summary(sub4)$coefficients[16,5]))

sub5 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_pallidumlh_vol_cent + scale(subj_pallidumlh_vol) + site_pallidumlh_vol_cent*wave + scale(subj_pallidumlh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b5i <- as.data.frame(c(summary(sub5)$coefficients[14,1]))
sub_se5i <- as.data.frame(c(summary(sub5)$coefficients[14,2]))
sub_p5i <- as.data.frame(c(summary(sub5)$coefficients[14,5]))
sub_b5s <- as.data.frame(c(summary(sub5)$coefficients[16,1]))
sub_se5s <- as.data.frame(c(summary(sub5)$coefficients[16,2]))
sub_p5s <- as.data.frame(c(summary(sub5)$coefficients[16,5]))

sub6 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_bstem_vol_cent + scale(subj_bstem_vol) + site_bstem_vol_cent*wave + scale(subj_bstem_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b6i <- as.data.frame(c(summary(sub6)$coefficients[14,1]))
sub_se6i <- as.data.frame(c(summary(sub6)$coefficients[14,2]))
sub_p6i <- as.data.frame(c(summary(sub6)$coefficients[14,5]))
sub_b6s <- as.data.frame(c(summary(sub6)$coefficients[16,1]))
sub_se6s <- as.data.frame(c(summary(sub6)$coefficients[16,2]))
sub_p6s <- as.data.frame(c(summary(sub6)$coefficients[16,5]))

sub7 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_hpuslh_vol_cent + scale(subj_hpuslh_vol) + site_hpuslh_vol_cent*wave + scale(subj_hpuslh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b7i <- as.data.frame(c(summary(sub7)$coefficients[14,1]))
sub_se7i <- as.data.frame(c(summary(sub7)$coefficients[14,2]))
sub_p7i <- as.data.frame(c(summary(sub7)$coefficients[14,5]))
sub_b7s <- as.data.frame(c(summary(sub7)$coefficients[16,1]))
sub_se7s <- as.data.frame(c(summary(sub7)$coefficients[16,2]))
sub_p7s <- as.data.frame(c(summary(sub7)$coefficients[16,5]))

sub8 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_amygdalalh_vol_cent + scale(subj_amygdalalh_vol) + site_amygdalalh_vol_cent*wave + scale(subj_amygdalalh_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b8i <- as.data.frame(c(summary(sub8)$coefficients[14,1]))
sub_se8i <- as.data.frame(c(summary(sub8)$coefficients[14,2]))
sub_p8i <- as.data.frame(c(summary(sub8)$coefficients[14,5]))
sub_b8s <- as.data.frame(c(summary(sub8)$coefficients[16,1]))
sub_se8s <- as.data.frame(c(summary(sub8)$coefficients[16,2]))
sub_p8s <- as.data.frame(c(summary(sub8)$coefficients[16,5]))

sub9 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
               site_aal_vol_cent + scale(subj_aal_vol) + site_aal_vol_cent*wave + scale(subj_aal_vol)*wave +
               (1 + wave|site_id/id), 
             data = pfactor_red, 
             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b9i <- as.data.frame(c(summary(sub9)$coefficients[14,1]))
sub_se9i <- as.data.frame(c(summary(sub9)$coefficients[14,2]))
sub_p9i <- as.data.frame(c(summary(sub9)$coefficients[14,5]))
sub_b9s <- as.data.frame(c(summary(sub9)$coefficients[16,1]))
sub_se9s <- as.data.frame(c(summary(sub9)$coefficients[16,2]))
sub_p9s <- as.data.frame(c(summary(sub9)$coefficients[16,5]))

sub10 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_vedclh_vol_cent + scale(subj_vedclh_vol) + site_vedclh_vol_cent*wave + scale(subj_vedclh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b10i <- as.data.frame(c(summary(sub10)$coefficients[14,1]))
sub_se10i <- as.data.frame(c(summary(sub10)$coefficients[14,2]))
sub_p10i <- as.data.frame(c(summary(sub10)$coefficients[14,5]))
sub_b10s <- as.data.frame(c(summary(sub10)$coefficients[16,1]))
sub_se10s <- as.data.frame(c(summary(sub10)$coefficients[16,2]))
sub_p10s <- as.data.frame(c(summary(sub10)$coefficients[16,5]))

sub11 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_crbcortexrh_vol_cent + scale(subj_crbcortexrh_vol) + site_crbcortexrh_vol_cent*wave + scale(subj_crbcortexrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b11i <- as.data.frame(c(summary(sub11)$coefficients[14,1]))
sub_se11i <- as.data.frame(c(summary(sub11)$coefficients[14,2]))
sub_p11i <- as.data.frame(c(summary(sub11)$coefficients[14,5]))
sub_b11s <- as.data.frame(c(summary(sub11)$coefficients[16,1]))
sub_se11s <- as.data.frame(c(summary(sub11)$coefficients[16,2]))
sub_p11s <- as.data.frame(c(summary(sub11)$coefficients[16,5]))

sub12 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_tprh_vol_cent + scale(subj_tprh_vol) + site_tprh_vol_cent*wave + scale(subj_tprh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b12i <- as.data.frame(c(summary(sub12)$coefficients[14,1]))
sub_se12i <- as.data.frame(c(summary(sub12)$coefficients[14,2]))
sub_p12i <- as.data.frame(c(summary(sub12)$coefficients[14,5]))
sub_b12s <- as.data.frame(c(summary(sub12)$coefficients[16,1]))
sub_se12s <- as.data.frame(c(summary(sub12)$coefficients[16,2]))
sub_p12s <- as.data.frame(c(summary(sub12)$coefficients[16,5]))

sub13 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_caudaterh_vol_cent + scale(subj_caudaterh_vol) + site_caudaterh_vol_cent*wave + scale(subj_caudaterh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b13i <- as.data.frame(c(summary(sub13)$coefficients[14,1]))
sub_se13i <- as.data.frame(c(summary(sub13)$coefficients[14,2]))
sub_p13i <- as.data.frame(c(summary(sub13)$coefficients[14,5]))
sub_b13s <- as.data.frame(c(summary(sub13)$coefficients[16,1]))
sub_se13s <- as.data.frame(c(summary(sub13)$coefficients[16,2]))
sub_p13s <- as.data.frame(c(summary(sub13)$coefficients[16,5]))

sub14 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_putamenrh_vol_cent + scale(subj_putamenrh_vol) + site_putamenrh_vol_cent*wave + scale(subj_putamenrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b14i <- as.data.frame(c(summary(sub14)$coefficients[14,1]))
sub_se14i <- as.data.frame(c(summary(sub14)$coefficients[14,2]))
sub_p14i <- as.data.frame(c(summary(sub14)$coefficients[14,5]))
sub_b14s <- as.data.frame(c(summary(sub14)$coefficients[16,1]))
sub_se14s <- as.data.frame(c(summary(sub14)$coefficients[16,2]))
sub_p14s <- as.data.frame(c(summary(sub14)$coefficients[16,5]))

sub15 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_pallidumrh_vol_cent + scale(subj_pallidumrh_vol) + site_pallidumrh_vol_cent*wave + scale(subj_pallidumrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b15i <- as.data.frame(c(summary(sub15)$coefficients[14,1]))
sub_se15i <- as.data.frame(c(summary(sub15)$coefficients[14,2]))
sub_p15i <- as.data.frame(c(summary(sub15)$coefficients[14,5]))
sub_b15s <- as.data.frame(c(summary(sub15)$coefficients[16,1]))
sub_se15s <- as.data.frame(c(summary(sub15)$coefficients[16,2]))
sub_p15s <- as.data.frame(c(summary(sub15)$coefficients[16,5]))

sub16 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_hpusrh_vol_cent + scale(subj_hpusrh_vol) + site_hpusrh_vol_cent*wave + scale(subj_hpusrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b16i <- as.data.frame(c(summary(sub16)$coefficients[14,1]))
sub_se16i <- as.data.frame(c(summary(sub16)$coefficients[14,2]))
sub_p16i <- as.data.frame(c(summary(sub16)$coefficients[14,5]))
sub_b16s <- as.data.frame(c(summary(sub16)$coefficients[16,1]))
sub_se16s <- as.data.frame(c(summary(sub16)$coefficients[16,2]))
sub_p16s <- as.data.frame(c(summary(sub16)$coefficients[16,5]))

sub17 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_amygdalarh_vol_cent + scale(subj_amygdalarh_vol) + site_amygdalarh_vol_cent*wave + scale(subj_amygdalarh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b17i <- as.data.frame(c(summary(sub17)$coefficients[14,1]))
sub_se17i <- as.data.frame(c(summary(sub17)$coefficients[14,2]))
sub_p17i <- as.data.frame(c(summary(sub17)$coefficients[14,5]))
sub_b17s <- as.data.frame(c(summary(sub17)$coefficients[16,1]))
sub_se17s <- as.data.frame(c(summary(sub17)$coefficients[16,2]))
sub_p17s <- as.data.frame(c(summary(sub17)$coefficients[16,5]))

sub18 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_aar_vol_cent + scale(subj_aar_vol) + site_aar_vol_cent*wave + scale(subj_aar_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b18i <- as.data.frame(c(summary(sub18)$coefficients[14,1]))
sub_se18i <- as.data.frame(c(summary(sub18)$coefficients[14,2]))
sub_p18i <- as.data.frame(c(summary(sub18)$coefficients[14,5]))
sub_b18s <- as.data.frame(c(summary(sub18)$coefficients[16,1]))
sub_se18s <- as.data.frame(c(summary(sub18)$coefficients[16,2]))
sub_p18s <- as.data.frame(c(summary(sub18)$coefficients[16,5]))

sub19 <- lmer(scalep ~ wave + sex + age + black + asian + hisp + other + achieva + discovery + ingenia + prisma +  
                site_vedcrh_vol_cent + scale(subj_vedcrh_vol) + site_vedcrh_vol_cent*wave + scale(subj_vedcrh_vol)*wave +
                (1 + wave|site_id/id), 
              data = pfactor_red, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
sub_b19i <- as.data.frame(c(summary(sub19)$coefficients[14,1]))
sub_se19i <- as.data.frame(c(summary(sub19)$coefficients[14,2]))
sub_p19i <- as.data.frame(c(summary(sub19)$coefficients[14,5]))
sub_b19s <- as.data.frame(c(summary(sub19)$coefficients[16,1]))
sub_se19s <- as.data.frame(c(summary(sub19)$coefficients[16,2]))
sub_p19s <- as.data.frame(c(summary(sub19)$coefficients[16,5]))

#create data frame with all parcel-wise cortical volume st. estimates, SEs, and p-values
sub_parcel <- data.frame(x=c("sub"))

newsubbi <- c(sub_b1i,sub_b2i,sub_b3i,sub_b4i,sub_b5i,sub_b6i,sub_b7i,sub_b8i,sub_b9i,sub_b10i,sub_b11i,sub_b12i,sub_b13i,sub_b14i,sub_b15i,
              sub_b16i,sub_b17i,sub_b18i,sub_b19i)
sub_parcel <- cbind(sub_parcel,newsubbi)

newsubbs <- c(sub_b1s,sub_b2s,sub_b3s,sub_b4s,sub_b5s,sub_b6s,sub_b7s,sub_b8s,sub_b9s,sub_b10s,sub_b11s,sub_b12s,sub_b13s,sub_b14s,sub_b15s,
              sub_b16s,sub_b17s,sub_b18s,sub_b19s)
sub_parcel <- cbind(sub_parcel,newsubbs)

newsubsei <- c(sub_se1i,sub_se2i,sub_se3i,sub_se4i,sub_se5i,sub_se6i,sub_se7i,sub_se8i,sub_se9i,sub_se10i,sub_se11i,sub_se12i,sub_se13i,
               sub_se14i,sub_se15i,sub_se16i,sub_se17i,sub_se18i,sub_se19i)
sub_parcel <- cbind(sub_parcel,newsubsei)

newsubses <- c(sub_se1s,sub_se2s,sub_se3s,sub_se4s,sub_se5s,sub_se6s,sub_se7s,sub_se8s,sub_se9s,sub_se10s,sub_se11s,sub_se12s,sub_se13s,
               sub_se14s,sub_se15s,sub_se16s,sub_se17s,sub_se18s,sub_se19s)
sub_parcel <- cbind(sub_parcel,newsubses)

newsubpi <- c(sub_p1i,sub_p2i,sub_p3i,sub_p4i,sub_p5i,sub_p6i,sub_p7i,sub_p8i,sub_p9i,sub_p10i,sub_p11i,sub_p12i,sub_p13i,sub_p14i,sub_p15i,
              sub_p16i,sub_p17i,sub_p18i,sub_p19i)
sub_parcel <- cbind(sub_parcel,newsubpi)

newsubps <- c(sub_p1s,sub_p2s,sub_p3s,sub_p4s,sub_p5s,sub_p6s,sub_p7s,sub_p8s,sub_p9s,sub_p10s,sub_p11s,sub_p12s,sub_p13s,sub_p14s,sub_p15s,
              sub_p16s,sub_p17s,sub_p18s,sub_p19s)
sub_parcel <- cbind(sub_parcel,newsubps)

names(sub_parcel) <- c('sub','crbcortexlh_vol_bi','tplh_vol_bi','caudatelh_vol_bi','putamenlh_vol_bi','pallidumlh_vol_bi',
                       'bstem_vol_bi','hpuslh_vol_bi','amygdalalh_vol_bi','aal_vol_bi','vedclh_vol_bi','crbcortexrh_vol_bi',
                       'tprh_vol_bi','caudaterh_vol_bi','putamenrh_vol_bi','pallidumrh_vol_bi','hpusrh_vol_bi',
                       'amygdalarh_vol_bi','aar_vol_bi','vedcrh_vol_bi',
                       'crbcortexlh_vol_bs','tplh_vol_bs','caudatelh_vol_bs','putamenlh_vol_bs','pallidumlh_vol_bs',
                       'bstem_vol_bs','hpuslh_vol_bs','amygdalalh_vol_bs','aal_vol_bs','vedclh_vol_bs','crbcortexrh_vol_bs',
                       'tprh_vol_bs','caudaterh_vol_bs','putamenrh_vol_bs','pallidumrh_vol_bs','hpusrh_vol_bs',
                       'amygdalarh_vol_bs','aar_vol_bs','vedcrh_vol_bs',
                       'crbcortexlh_vol_sei','tplh_vol_sei','caudatelh_vol_sei','putamenlh_vol_sei','pallidumlh_vol_sei',
                       'bstem_vol_sei','hpuslh_vol_sei','amygdalalh_vol_sei','aal_vol_sei','vedclh_vol_sei','crbcortexrh_vol_sei',
                       'tprh_vol_sei','caudaterh_vol_sei','putamenrh_vol_sei','pallidumrh_vol_sei','hpusrh_vol_sei',
                       'amygdalarh_vol_sei','aar_vol_sei','vedcrh_vol_sei',
                       'crbcortexlh_vol_ses','tplh_vol_ses','caudatelh_vol_ses','putamenlh_vol_ses','pallidumlh_vol_ses',
                       'bstem_vol_ses','hpuslh_vol_ses','amygdalalh_vol_ses','aal_vol_ses','vedclh_vol_ses','crbcortexrh_vol_ses',
                       'tprh_vol_ses','caudaterh_vol_ses','putamenrh_vol_ses','pallidumrh_vol_ses','hpusrh_vol_ses',
                       'amygdalarh_vol_ses','aar_vol_ses','vedcrh_vol_ses',
                       'crbcortexlh_vol_pi','tplh_vol_pi','caudatelh_vol_pi','putamenlh_vol_pi','pallidumlh_vol_pi',
                       'bstem_vol_pi','hpuslh_vol_pi','amygdalalh_vol_pi','aal_vol_pi','vedclh_vol_pi','crbcortexrh_vol_pi',
                       'tprh_vol_pi','caudaterh_vol_pi','putamenrh_vol_pi','pallidumrh_vol_pi','hpusrh_vol_pi',
                       'amygdalarh_vol_pi','aar_vol_pi','vedcrh_vol_pi',
                       'crbcortexlh_vol_ps','tplh_vol_ps','caudatelh_vol_ps','putamenlh_vol_ps','pallidumlh_vol_ps',
                       'bstem_vol_ps','hpuslh_vol_ps','amygdalalh_vol_ps','aal_vol_ps','vedclh_vol_ps','crbcortexrh_vol_ps',
                       'tprh_vol_ps','caudaterh_vol_ps','putamenrh_vol_ps','pallidumrh_vol_ps','hpusrh_vol_ps',
                       'amygdalarh_vol_ps','aar_vol_ps','vedcrh_vol_ps')

#calculate 95% CIs and create lower and upper bound variables 
sub_parcel$ci_lower_crbcortexlh_voli <- sub_parcel$crbcortexlh_vol_bi - 1.96*sub_parcel$crbcortexlh_vol_sei
sub_parcel$ci_upper_crbcortexlh_voli <- sub_parcel$crbcortexlh_vol_bi + 1.96*sub_parcel$crbcortexlh_vol_sei
sub_parcel$ci_lower_tplh_voli <- sub_parcel$tplh_vol_bi - 1.96*sub_parcel$tplh_vol_sei
sub_parcel$ci_upper_tplh_voli <- sub_parcel$tplh_vol_bi + 1.96*sub_parcel$tplh_vol_sei
sub_parcel$ci_lower_caudatelh_voli <- sub_parcel$caudatelh_vol_bi - 1.96*sub_parcel$caudatelh_vol_sei 
sub_parcel$ci_upper_caudatelh_voli <- sub_parcel$caudatelh_vol_bi + 1.96*sub_parcel$caudatelh_vol_sei
sub_parcel$ci_lower_putamenlh_voli <- sub_parcel$putamenlh_vol_bi - 1.96*sub_parcel$putamenlh_vol_sei 
sub_parcel$ci_upper_putamenlh_voli <- sub_parcel$putamenlh_vol_bi + 1.96*sub_parcel$putamenlh_vol_sei
sub_parcel$ci_lower_pallidumlh_voli <- sub_parcel$pallidumlh_vol_bi - 1.96*sub_parcel$pallidumlh_vol_sei 
sub_parcel$ci_upper_pallidumlh_voli <- sub_parcel$pallidumlh_vol_bi + 1.96*sub_parcel$pallidumlh_vol_sei
sub_parcel$ci_lower_bstem_voli <- sub_parcel$bstem_vol_bi - 1.96*sub_parcel$bstem_vol_sei 
sub_parcel$ci_upper_bstem_voli <- sub_parcel$bstem_vol_bi + 1.96*sub_parcel$bstem_vol_sei
sub_parcel$ci_lower_hpuslh_voli <- sub_parcel$hpuslh_vol_bi - 1.96*sub_parcel$hpuslh_vol_sei 
sub_parcel$ci_upper_hpuslh_voli <- sub_parcel$hpuslh_vol_bi + 1.96*sub_parcel$hpuslh_vol_sei
sub_parcel$ci_lower_amygdalalh_voli <- sub_parcel$amygdalalh_vol_bi - 1.96*sub_parcel$amygdalalh_vol_sei 
sub_parcel$ci_upper_amygdalalh_voli <- sub_parcel$amygdalalh_vol_bi + 1.96*sub_parcel$amygdalalh_vol_sei
sub_parcel$ci_lower_aal_voli <- sub_parcel$aal_vol_bi - 1.96*sub_parcel$aal_vol_sei 
sub_parcel$ci_upper_aal_voli <- sub_parcel$aal_vol_bi + 1.96*sub_parcel$aal_vol_sei
sub_parcel$ci_lower_vedclh_voli <- sub_parcel$vedclh_vol_bi - 1.96*sub_parcel$vedclh_vol_sei 
sub_parcel$ci_upper_vedclh_voli <- sub_parcel$vedclh_vol_bi + 1.96*sub_parcel$vedclh_vol_sei

sub_parcel$ci_lower_crbcortexrh_voli <- sub_parcel$crbcortexrh_vol_bi - 1.96*sub_parcel$crbcortexrh_vol_sei 
sub_parcel$ci_upper_crbcortexrh_voli <- sub_parcel$crbcortexrh_vol_bi + 1.96*sub_parcel$crbcortexrh_vol_sei
sub_parcel$ci_lower_tprh_voli <- sub_parcel$tprh_vol_bi - 1.96*sub_parcel$tprh_vol_sei 
sub_parcel$ci_upper_tprh_voli <- sub_parcel$tprh_vol_bi + 1.96*sub_parcel$tprh_vol_sei
sub_parcel$ci_lower_caudaterh_voli <- sub_parcel$caudaterh_vol_bi - 1.96*sub_parcel$caudaterh_vol_sei
sub_parcel$ci_upper_caudaterh_voli <- sub_parcel$caudaterh_vol_bi + 1.96*sub_parcel$caudaterh_vol_sei
sub_parcel$ci_lower_putamenrh_voli <- sub_parcel$putamenrh_vol_bi - 1.96*sub_parcel$putamenrh_vol_sei 
sub_parcel$ci_upper_putamenrh_voli <- sub_parcel$putamenrh_vol_bi + 1.96*sub_parcel$putamenrh_vol_sei
sub_parcel$ci_lower_pallidumrh_voli <- sub_parcel$pallidumrh_vol_bi - 1.96*sub_parcel$pallidumrh_vol_sei 
sub_parcel$ci_upper_pallidumrh_voli <- sub_parcel$pallidumrh_vol_bi + 1.96*sub_parcel$pallidumrh_vol_sei
sub_parcel$ci_lower_hpusrh_voli <- sub_parcel$hpusrh_vol_bi - 1.96*sub_parcel$hpusrh_vol_sei 
sub_parcel$ci_upper_hpusrh_voli <- sub_parcel$hpusrh_vol_bi + 1.96*sub_parcel$hpusrh_vol_sei
sub_parcel$ci_lower_amygdalarh_voli <- sub_parcel$amygdalarh_vol_bi - 1.96*sub_parcel$amygdalarh_vol_sei 
sub_parcel$ci_upper_amygdalarh_voli <- sub_parcel$amygdalarh_vol_bi + 1.96*sub_parcel$amygdalarh_vol_sei
sub_parcel$ci_lower_aar_voli <- sub_parcel$aar_vol_bi - 1.96*sub_parcel$aar_vol_sei 
sub_parcel$ci_upper_aar_voli <- sub_parcel$aar_vol_bi + 1.96*sub_parcel$aar_vol_sei
sub_parcel$ci_lower_vedcrh_voli <- sub_parcel$vedcrh_vol_bi - 1.96*sub_parcel$vedcrh_vol_sei 
sub_parcel$ci_upper_vedcrh_voli <- sub_parcel$vedcrh_vol_bi + 1.96*sub_parcel$vedcrh_vol_sei

write.csv(sub_parcel, "Parcel-Wise Subcortical Volume MLM Analysis Output Standardized_FINAL.csv") #write csv file


#ggseg figure of parcel-wise cortical volume analyses with p showing std. betas - for parcels that survived FDR correction (Figure 2A)

vol_results_fdr= data.frame(cbind(region=c("bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",
                                           "fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital",
                                           "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","parahippocampal",
                                           "paracentral","pars opercularis","pars orbitalis","pars triangularis","pericalcarine",
                                           "postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
                                           "rostral middle frontal","superior frontal","superior parietal","superior temporal",
                                           "supramarginal","temporal pole","transverse temporal","insula",
                                           "bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",
                                           "fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital",
                                           "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","parahippocampal",
                                           "paracentral","pars opercularis","pars orbitalis","pars triangularis","pericalcarine",
                                           "postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
                                           "rostral middle frontal","superior frontal","superior parietal","superior temporal",
                                           "supramarginal","frontal pole","temporal pole","transverse temporal","insula"),
                                  stdb=c(-0.042368321,-0.035717618,-0.036385733,-0.025642445,-0.03143749,-0.064598625,-0.069576686,
                                         -0.054125201,-0.053063053,-0.042422373,-0.063394951,-0.051823327,-0.035497437,-0.074020718,
                                         -0.046714585,-0.047766861,-0.022548963,-0.042677898,-0.031600652,-0.02184436,-0.07081975,
                                         -0.056360759,-0.081063624,-0.057722969,-0.053262367,-0.058546438,-0.066479253,-0.060252816,
                                         -0.060426257,-0.077828238,-0.041944366,-0.02965012,-0.062519618,-0.036433259,-0.049468731,
                                         -0.038449068,-0.028045668,-0.042602437,-0.0730509,-0.07345152,-0.073528626,-0.033691576,
                                         -0.040178517,-0.061678671,-0.03090881,-0.043071623,-0.081208478,-0.038049519,-0.058252336,
                                         -0.029326333,-0.049818597,-0.02290476,-0.024347901,-0.077672815,-0.059497585,-0.076860362,
                                         -0.054229699,-0.045887035,-0.058750679,-0.060038727,-0.045835286,-0.056136484,-0.068509275,
                                         -0.068412512,-0.032373698,-0.038495449,-0.067274413),
                                  hemi=c("left","left","left","left","left","left","left","left","left","left","left","left","left","left",
                                         "left","left","left","left","left","left","left","left","left","left","left","left","left","left",
                                         "left","left","left","left","left","right","right","right","right","right","right","right",
                                         "right","right","right","right","right","right","right","right","right","right","right","right","right",
                                         "right","right","right","right","right","right","right","right","right","right","right","right","right",
                                         "right")),
                            stringsAsFactors=F)

vol_results_fdr %>% 
  ggseg(mapping=aes(fill=as.numeric(stdb)),
        colour="black",size=.6) +
  scale_fill_gradient(low = "yellow",
                      high = "red") +
  ggtitle("Parcel-Wise Cortical Volume")

#ggseg figure of parcel-wise surface area analyses with p showing std. betas - for parcels that survived FDR correction (Figure 2B)

area_results_fdr = data.frame(cbind(region=c("bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",
                                             "fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital",
                                             "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","parahippocampal",
                                             "paracentral","pars opercularis","pars orbitalis","pars triangularis","pericalcarine",
                                             "postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
                                             "rostral middle frontal","superior frontal","superior parietal","superior temporal",
                                             "supramarginal","frontal pole","temporal pole","transverse temporal","insula",
                                             "bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",
                                             "fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital",
                                             "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","parahippocampal",
                                             "paracentral","pars opercularis","pars orbitalis","pars triangularis","pericalcarine",
                                             "postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
                                             "rostral middle frontal","superior frontal","superior parietal","superior temporal",
                                             "supramarginal","frontal pole","temporal pole","transverse temporal","insula"),
                                    stdb=c(-0.038278791,-0.044845218,-0.042424554,-0.036515479,-0.030829005,-0.067368816,-0.066182309,
                                           -0.052949451,-0.058380138,-0.043662973,-0.065922201,-0.048812821,-0.044014255,-0.06448414,
                                           -0.046843378,-0.04401945,-0.033647862,-0.056584523,-0.046453091,-0.029500029,-0.079751855,
                                           -0.068358722,-0.073826869,-0.063881976,-0.061525938,-0.066855365,-0.072314495,-0.056852908,
                                           -0.071361287,-0.077119706,-0.058310678,-0.053063398,-0.031987585,-0.067686192,-0.043065692,
                                           -0.061373984,-0.048279067,-0.038015022,-0.040751586,-0.080862887,-0.076608991,-0.052949451,
                                           -0.042308527,-0.04172251,-0.065465506,-0.032785757,-0.068806165,-0.084900452,-0.045006604,
                                           -0.055122552,-0.036044644,-0.069159196,-0.035125562,-0.029212916,-0.080899515,-0.064136644,
                                           -0.070306278,-0.061630446,-0.060028754,-0.068942733,-0.070410928,-0.047015608,-0.065977443,
                                           -0.070331422,-0.070333964,-0.047115717,-0.042272271,-0.058706186),
                                    hemi=c("left","left","left","left","left","left","left","left","left","left","left","left","left","left",
                                           "left","left","left","left","left","left","left","left","left","left","left","left","left","left",
                                           "left","left","left","left","left","left","right","right","right","right","right","right","right",
                                           "right","right","right","right","right","right","right","right","right","right","right","right","right",
                                           "right","right","right","right","right","right","right","right","right","right","right","right","right",
                                           "right")),
                              stringsAsFactors=F)

area_results_fdr %>% 
  ggseg(mapping=aes(fill=as.numeric(stdb)), 
        colour="black",size=.6) +
  scale_fill_gradient(low = "yellow",
                      high = "red") +
  ggtitle("Parcel-Wise Cortical Surface Area")

ggseg(atlas=aseg,mapping=aes(fill=region))


#ggseg figure of parcel-wise cortical thickness analyses with INT showing std. betas - for parcels that survived FDR correction (Figure S2)

thick_results_fdr = data.frame(cbind(region=c("inferior temporal","lingual","middle temporal","parahippocampal","paracentral","postcentral","precentral",
                                              "superior parietal","cuneus","lingual","postcentral","superior parietal"),
                                     stdb=c(0.01397659,0.013431228,0.013937159,0.013453852,0.014844924,0.018167898,0.015874702,
                                            0.013916114,0.013694289,0.016436235,0.016635508,0.013832886),
                                     hemi=c("left","left","left","left","left","left","left","left","right","right","right","right")),
                               stringsAsFactors=F)

thick_results_fdr %>% 
  ggseg(mapping=aes(fill=as.numeric(stdb)), position="stacked",
        colour="black",size=.6) +
  scale_fill_gradient(high = "yellow",
                      low = "red") +
  ggtitle("Parcel-Wise Cortical Thickness")


#forest plot of parcel-wise subcortical volume analysis with p (Figure 2C)

sub_forest <- data.frame(variable=c("Left Cerebellar Cortex","Left Thalamus Proper","Left Caudate","Left Putamen","Left Pallidum",
                                    "Brainstem","Left Hippocampus","Left Amygdala","Left Accumbens","Left Ventral DC",
                                    "Right Cerebellar Cortex","Right Thalamus Proper","Right Caudate","Right Putamen","Right Pallidum",
                                    "Right Hippocampus","Right Amygdala","Right Accumbens","Right Ventral DC"),
                         stdb=c(-0.062526529,-0.068174585,-0.050167024,-0.04306801,-0.043305729,-0.069430328,-0.068697041,-0.053180425,
                                -0.04888736,-0.064164359,-0.055479243,-0.055929446,-0.046133413,-0.050503415,-0.030368098,-0.061715518,
                                -0.055990013,-0.05952229,-0.058611321),
                         ub=c(-0.040013875,-0.046446858,-0.02947583,-0.021562199,-0.022104379,-0.047509596,-0.046879632,-0.031469191,
                              -0.027605907,-0.04211372,-0.032795663,-0.034034863,-0.025415504,-0.028962497,-0.009133533,-0.040122157,
                              -0.034160325,-0.038574215,-0.03650152),
                         lb=c(-0.085039183,-0.089902312,-0.070858217,-0.064573821,-0.064507078,-0.09135106,-0.09051445,-0.07489166,
                              -0.070168813,-0.086214997,-0.078162823,-0.07782403,-0.066851322,-0.072044334,-0.051602664,-0.083308879,
                              -0.077819702,-0.080470366,-0.080721123))
tabletext <- cbind(
  c(expr("Left Cerebellar Cortex"),expr("Left Thalamus Proper"),expr("Left Caudate"),expr("Left Putamen"),expr("Left Pallidum"),
    expr("Brainstem"),expr("Left Hippocampus"),expr("Left Amygdala"),expr("Left Accumbens"),expr("Left Ventral Diencephalon"),
    expr("Right Cerebellar Cortex"),expr("Right Thalamus Proper"),expr("Right Caudate"),expr("Right Putamen"),expr("Right Pallidum"),
    expr("Right Hippocampus"),expr("Right Amygdala"),expr("Right Accumbens"),expr("Right Ventral Diencephalon")) 
)
png(paste("Parcel_wise_subcortical_volume_forest_stdb.png",sep=""),height=6,width=5,res=300,units="in")
print( forestplot(tabletext,
                  mean =sub_forest$stdb,
                  lower=sub_forest$lb,
                  upper=sub_forest$ub,
                  boxsize=.15,
                  txt_gp = fpTxtGp( ticks=gpar(fontfamily="", cex=1), xlab=gpar(fontfamily="", cex=1.25) ),
                  xlab="Standardized Beta" ) )
dev.off()

