############################################################################################
# Analysis of qPCR Data
#############################################################################################

# Date : 21 November 2023
# Author : Ari Fina Bintarti

# Install packages
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
install.packages("ggpubr")
install.packages("car")
install.packages("agricolae")
install.packages("multcompView")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("ggtext")
install.packages("sjmisc") 
install.packages("sjPlot")
install.packages("MASS")
install.packages("FSA")
install.packages('mvtnorm', dep = TRUE)
install.packages("rcompanion")
install.packages("onewaytests")
install.packages("PerformanceAnalytics")
install.packages("gvlma")
install.packages("ggpmisc")
install.packages("fitdistrplus")
install.packages('BiocManager')
install.packages("dplyr")
install.packages("lme4")
install.packages("nlme")
install.packages("lmerTest")
install.packages("multcomp")
install.packages("ape")
install.packages("devtools", dependencies = TRUE)
BiocManager::install("phyloseq")
install.packages("datarium")
install.packages("rstatix")
install.packages("export")
library(devtools)
library(multcomp)
library(car)
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggtext)
library(reshape)
library(ggpubr)
library(car)
library(agricolae)
library(multcompView)
library(grid)
library(gridExtra)
library(sjmisc)
library(sjPlot)
library(MASS)
library(FSA)
library(rcompanion)
library(onewaytests)
library(ggsignif)
library(PerformanceAnalytics)
library(gvlma)
library(ggpmisc)
library(tibble)
library(fitdistrplus)
library(lme4)
library(nlme)
library(ape)
library(lmerTest)
#library(phyloseq)
library(datarium)
library(rstatix)
library(export)
#install.packages("afex") 
#install.packages("qqplotr")
library(afex)
library(performance)
library(qqplotr)
library(emmeans)
library(reshape2)

###########################################################################
# 1. Response variable: BACTERIAL AND ARCHAEAL COPY NUMBER
###########################################################################

#### 1a. Analyses of Bulk Soil - Copies per Gram Dry weight of Soil ####
setwd('/Users/arifinabintarti/Documents/France/microservices/')
#setwd('D:/Fina/INRAE_Project/microservices/')
qPCR <- read.csv("qPCR_results_LP_stat.csv")
qPCR.BS <- qPCR[1:120,]
str(qPCR.BS)
qPCR.BS$nr <- as.factor(qPCR.BS$nr)
qPCR.BS$plot <- as.factor(qPCR.BS$plot)
qPCR.BS$block <- as.factor(qPCR.BS$block)
qPCR.BS$irrigation <- as.factor(qPCR.BS$irrigation)
qPCR.BS$fertilization <- as.factor(qPCR.BS$fertilization)
qPCR.BS$type <- as.factor(qPCR.BS$type)
qPCR.BS$sampling.date <- as.factor(qPCR.BS$sampling.date)
qPCR.BS$var3 <- as.factor(qPCR.BS$var3)
qPCR.BS$x <- as.factor(qPCR.BS$x)
qPCR.BS$rep <- as.factor(qPCR.BS$rep)

qPCR.BS$sampling.date <- factor(qPCR.BS$sampling.date, levels = c("28/04/2022", "1/6/22", "5/7/22", "20/07/2022", "13/09/2022"),
                                labels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"))
str(qPCR.BS)

# perform log transformation
qPCR.BS$AOA_logDWS <- log10(qPCR.BS$AOA_nbc_per_g_DW_soil)
qPCR.BS$AOB_logDWS <- log10(qPCR.BS$AOB_nbc_per_g_DW_soil)
qPCR.BS$ComA_logDWS <- log10(qPCR.BS$ComA_nbc_per_g_DW_soil)
qPCR.BS$ComB_logDWS <- log10(qPCR.BS$ComB_nbc_per_g_DW_soil)
qPCR.BS$Tot_logDWS <- log10(qPCR.BS$Tot_nbc_per_g_DW_soil)
# perform arcsin root square
qPCR.BS$AOA_16.arc.ratio <- asin(sqrt(qPCR.BS$AOA_16S_ratio_percent / 100))
qPCR.BS$AOB_16.arc.ratio <- asin(sqrt(qPCR.BS$AOB_16S_ratio_percent / 100))
qPCR.BS$ComA_16.arc.ratio <- asin(sqrt(qPCR.BS$ComA_16S_ratio_percent / 100))
qPCR.BS$ComB_16.arc.ratio <- asin(sqrt(qPCR.BS$ComB_16S_ratio_percent / 100))
qPCR.BS$AOA_AOB.arc.ratio <- asin(sqrt(qPCR.BS$AOA_AOB_ratio / 100))
qPCR.BS$ComA_ComB.arc.ratio <- asin(sqrt(qPCR.BS$ComA_ComB_ratio / 100))

setwd('/Users/arifinabintarti/Documents/France/microservices/')
write.csv(qPCR.BS, file = "qPCR.BS.csv")
#qPCR.BS$AOA_logngDNA <- log10(qPCR.BS$AOA_nbc_per_ngDNA)
#qPCR.BS$AOB_logngDNA <- log10(qPCR.BS$AOB_nbc_per_ngDNA)

########################################################################################

#### 1. AOA Abundance per gram DWS ####

# anova test for non-transformed AOA 
# Linear mixed model
library(lmerTest)
# variation in the intercept among levels of f:g (the interaction between f and g)
# Intercept varying among block within sampling date.
# so we need a random effect enumerated by all of block by sampling date combinations
aoa.dws <- lmerTest::lmer(AOA_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date+
                       (1|sampling.date:block), data=qPCR.BS, na.action=na.omit)
anova(aoa.dws)
# test assumptions:
library(DHARMa)
shapiro.test(resid(aoa.dws)) # not normal
plot(simulateResiduals(aoa.dws)) # okay
#*** Need to transform the data
# anova test for transformed AOA 
aoa.log.dws <- lmerTest::lmer(AOA_logDWS ~ irrigation*fertilization*sampling.date+
                       (1|sampling.date:block), data=qPCR.BS, na.action=na.omit)
anova(aoa.log.dws)
aoa.log.dws2 <- lmerTest::lmer(AOA_logDWS ~ irrigation*fertilization*sampling.date+
                       (1|plot), data=qPCR.BS, na.action=na.omit)
anova(aoa.log.dws2)

aoa.log.dws.emm2 <- emmeans(aoa.log.dws2, ~ irrigation | fertilization*sampling.date)
aoa.log.dws.pair2 <- pairs(aoa.log.dws.emm2)
aoa.log.dws.pair2
# test assumptions:
shapiro.test(resid(aoa.log.dws)) # normal
plot(simulateResiduals(aoa.log.dws)) # okay
# Pairwise comparison:
aoa.log.dws.emm <- emmeans(aoa.log.dws2, ~ irrigation | fertilization*sampling.date)
aoa.log.dws.pair <- pairs(aoa.log.dws.emm)
aoa.log.dws.pair
aoa.log.dws.pair.DF <- as.data.frame(summary(aoa.log.dws.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
#write.csv(aoa.log.dws.pair.DF, file = "AOA_log_dws_pair.csv")


#### 2. AOB Abundance per gram DWS ####

# anova test for transformed AOB DWS
aob.dws <- lmerTest::lmer(AOB_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date+
                       (1|sampling.date:block), data=qPCR.BS, na.action=na.omit)
anova(aob.dws)
# test assumptions:
shapiro.test(resid(aob.dws)) # normal
plot(simulateResiduals(aob.dws)) # not that okay
#*** Need to transform the data
# anova test for transformed AOA 
aob.log.dws <- lmerTest::lmer(AOB_logDWS ~ irrigation*fertilization*sampling.date+
                       (1|sampling.date:block), data=qPCR.BS, na.action=na.omit)
anova(aob.log.dws)

aob.log.dws2 <- lmerTest::lmer(AOB_logDWS ~ irrigation*fertilization*sampling.date+
                       (1|plot), data=qPCR.BS, na.action=na.omit)
anova(aob.log.dws2)
# test assumptions:
shapiro.test(resid(aob.log.dws)) # normal
plot(simulateResiduals(aob.log.dws)) # better
# Pairwise comparison:
aob.log.dws.emm <- emmeans(aob.log.dws, ~ irrigation | fertilization*sampling.date)
aob.log.dws.pair <- pairs(aob.log.dws.emm)
aob.log.dws.pair
aob.log.dws.pair.DF <- as.data.frame(summary(aob.log.dws.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
#write.csv(aob.log.dws.pair.DF, file = "AOB_log_dws_pair.csv")


#### 3. COMAMMOX A Abundance per gram DWS ####

# anova test for non-transformed Comammox A DWS
comA.dws <- lmerTest::lmer(ComA_nbc_per_g_DW_soil ~ irrigation * fertilization * sampling.date + 
                          (1|sampling.date:block), data = qPCR.BS, na.action = na.omit)
anova(comA.dws)
# test assumptions:
shapiro.test(resid(comA.dws)) # very not normal
plot(simulateResiduals(comA.dws)) # not good
#*** Need to transform the data
# anova test for transformed Comammox A 
comA.log.dws <- lmerTest::lmer(ComA_logDWS ~ irrigation*fertilization*sampling.date+
                       (1|sampling.date:block), data=qPCR.BS, na.action=na.omit)
anova(comA.log.dws)

comA.log.dws2 <- lmerTest::lmer(ComA_logDWS ~ irrigation*fertilization*sampling.date+
                       (1|plot), data=qPCR.BS, na.action=na.omit)
anova(comA.log.dws2)

# test assumptions:
shapiro.test(resid(comA.log.dws)) # normal
plot(simulateResiduals(comA.log.dws)) # better
# Pairwise comparison:
comA.log.dws.emm <- emmeans(comA.log.dws, ~ irrigation | fertilization*sampling.date)
comA.log.dws.pair <- pairs(comA.log.dws.emm)
comA.log.dws.pair
comA.log.dws.pair.DF <- as.data.frame(summary(comA.log.dws.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
#write.csv(comA.log.dws.pair.DF, file = "ComA_log_dws_pair.csv")


#### 4. COMAMMOX B Abundance per gram DWS ####

# anova test for non-transformed Comammox B DWS
comB.dws <- lmerTest::lmer(ComB_nbc_per_g_DW_soil ~ irrigation * fertilization * sampling.date + 
                          (1|sampling.date:block), data = qPCR.BS, na.action = na.omit)
anova(comB.dws)
# test assumptions:
shapiro.test(resid(comB.dws)) # very not normal
plot(simulateResiduals(comB.dws)) # okay
#*** Need to transform the data
# anova test for transformed Comammox B 
comB.log.dws <- lmerTest::lmer(ComB_logDWS ~ irrigation*fertilization*sampling.date+
                       (1|sampling.date:block), data=qPCR.BS, na.action=na.omit)
anova(comB.log.dws)
comB.log.dws2 <- lmerTest::lmer(ComB_logDWS ~ irrigation*fertilization*sampling.date+
                       (1|plot), data=qPCR.BS, na.action=na.omit)
anova(comB.log.dws2)
0# test assumptions:
shapiro.test(resid(comB.log.dws)) # normal
plot(simulateResiduals(comB.log.dws)) # good
# Pairwise comparison:
comB.log.dws.emm <- emmeans(comB.log.dws, ~ irrigation | fertilization*sampling.date,lmer.df = "satterthwaite")
comB.log.dws.pair <- pairs(comB.log.dws.emm)
comB.log.dws.pair.DF <- as.data.frame(summary(comB.log.dws.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
#write.csv(comB.log.dws.pair.DF, file = "ComB_log_dws_pair.csv")

#### 5. 16S ####

# Linear mixed model
library(lmerTest)
# Intercept varying among block within sampling date.
# so we need a random effect enumerated by all of block by sampling date combinations
t.16S <- lmerTest::lmer(Tot_logDWS ~ irrigation*fertilization*sampling.date+
                       (1|sampling.date:block), data=qPCR.BS, na.action=na.omit)
anova(t.16S)

t.16S2 <- lmerTest::lmer(Tot_logDWS ~ irrigation*fertilization*sampling.date+(1|plot),
                         data=qPCR.BS, na.action=na.omit)
anova(t.16S2)
# lmer with unbalanced
t.16S3 = lmer(Tot_logDWS ~ irrigation*fertilization*sampling.date+(1|plot),
   contrasts = list(irrigation="contr.sum",fertilization="contr.sum",sampling.date="contr.sum"),data=qPCR.BS)
anova(t.16S3)
Anova(t.16S3, test="F", type="III") # similar with anova(t.16S3)
# pairwise
tot.log.dws.emm2 <- emmeans(t.16S2, ~ irrigation | fertilization*sampling.date)
tot.log.dws.pair2 <- pairs(tot.log.dws.emm2)
tot.log.dws.pair2

# checkin with aov
summary(aov(Tot_logDWS~irrigation*fertilization*sampling.date+
             Error(block/(irrigation*fertilization*sampling.date)), data=qPCR.BS))
summary(aov(Tot_logDWS~irrigation*fertilization*sampling.date+
             Error(plot), data=qPCR.BS))
summary(aov(Tot_logDWS~irrigation*fertilization*sampling.date+
             Error(block/sampling.date), data=qPCR.BS))



# test assumptions:
library(DHARMa)
shapiro.test(resid(t.16S)) # normal
plot(simulateResiduals(t.16S)) # okay
# Homogeneity of variance is not required at the nested level in linear mixed models. 
# Mixed models don’t assume homoskedasticity in the classical way. 
# Pairwise comparison:
library(emmeans)
tot.log.dws.emm <- emmeans(t.16S, ~ irrigation | fertilization*sampling.date)
tot.log.dws.pair <- pairs(tot.log.dws.emm)
tot.log.dws.pair.DF <- as.data.frame(summary(tot.log.dws.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil/')
#write.csv(tot.log.dws.pair.DF, file = "16S_log_dws_pair.csv") ### CHANGED


#### 6. AOA/16S RATIO ####

### anova test for AOA/16S Ratio Percentage
# Linear mixed model
t.aoa.16.rat.percent <- lmerTest::lmer(AOA_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                   (1|block:sampling.date), data=qPCR.BS, na.action=na.omit)
anova(t.aoa.16.rat.percent)
#specifies a separate intercept for each subject. 
t.aoa.16.rat.percent2 <- lmerTest::lmer(AOA_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                   (1|plot), data=qPCR.BS, na.action=na.omit) 
anova(t.aoa.16.rat.percent2)

#In order to allow for the change over time to differ across participants 
#(i.e. to explicitly model individual differences in change over time), 
#you also need to allow for the effect of Day to be random
t.aoa.16.rat.percent3 <- lmerTest::lmer(AOA_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                   (sampling.date|plot), data=qPCR.BS, na.action=na.omit,control=lmerControl(check.nobs.vs.nRE="ignore")) 
anova(t.aoa.16.rat.percent3)
# with anova (results )
summary(aov(AOA_16S_ratio_percent~irrigation*fertilization*sampling.date+
             Error(plot/(irrigation*fertilization*sampling.date)), data=qPCR.BS))



# test assumption
shapiro.test(resid(t.aoa.16.rat.percent)) # not normal
plot(simulateResiduals(t.aoa.16.rat.percent)) # not good
# ** Transformation needed
### anova test for AOA/16S Ratio Arcsin
t.aoa <- lmerTest::lmer(AOA_16.arc.ratio ~ irrigation*fertilization*sampling.date+
                   (1|block:sampling.date), data=qPCR.BS, na.action=na.omit)
anova(t.aoa)
t.aoa2 <- lmerTest::lmer(AOA_16.arc.ratio ~ irrigation*fertilization*sampling.date+
                   (1|plot), data=qPCR.BS, na.action=na.omit)
anova(t.aoa2)
# with anova (results are similar with lmer using random:  (1|plot) )
summary(aov(AOA_16.arc.ratio~irrigation*fertilization*sampling.date+
             Error(plot/(irrigation*fertilization*sampling.date)), data=qPCR.BS))
# test assumption
shapiro.test(resid(t.aoa)) # slightly not normal
plot(simulateResiduals(t.aoa)) # okay
# Pairwise comparison:
AOA_16.arc.ratio.emm <- emmeans(t.aoa, ~ irrigation | fertilization*sampling.date)
AOA_16.arc.ratio.pair <- pairs(AOA_16.arc.ratio.emm)
AOA_16.arc.ratio.pair.DF <- as.data.frame(summary(AOA_16.arc.ratio.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
#write.csv(AOA_16.arc.ratio.pair.DF, file = "AOA_16_arc_ratio_pair.csv")

#### 7. AOB/16S RATIO ####

### anova test for AOB/16S Ratio Percentage
# Linear mixed model
t.aob.16.rat.percent <- lmerTest::lmer(AOB_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                   (1|block:sampling.date), data=qPCR.BS, na.action=na.omit, REML=F)
anova(t.aob.16.rat.percent)

t.aob2 <- lmerTest::lmer(AOB_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                   (1|plot), contrasts = list(irrigation="contr.sum",fertilization="contr.sum",sampling.date="contr.sum"),
                   data=qPCR.BS)
anova(t.aob2)
anova(t.aob2,t.aob.16.rat.percent) #boundary (singular) fit: see help('isSingular') --> happens when the model too ccomplex
# t.aob.16.rat.percent is better based on AIC value

# test assumptions:
shapiro.test(resid(t.aob.16.rat.percent)) # normal
plot(simulateResiduals(t.aob.16.rat.percent)) # okay
# ** Transformation is not needed
# Pairwise comparison:
AOB_16.percent.ratio.emm <- emmeans(t.aob.16.rat.percent, ~ irrigation | fertilization*sampling.date)
AOB_16.percent.ratio.pair <- pairs(AOB_16.percent.ratio.emm)
AOB_16.percent.ratio.pair.DF <- as.data.frame(summary(AOB_16.percent.ratio.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil')
#write.csv(AOB_16.percent.ratio.pair.DF, file = "AOB_16_percent_ratio_pair.csv")
#**just checking
### anova test for AOB/16S Ratio Arcsin
# Linear mixed model
t.aob.16.rat.arcsin <- lmerTest::lmer(AOB_16.arc.ratio ~ irrigation*fertilization*sampling.date+
                   (1|block:sampling.date), data=qPCR.BS, na.action=na.omit)
anova(t.aob.16.rat.arcsin)
# test assumptions:
shapiro.test(resid(t.aob.16.rat.arcsin)) # normal
plot(simulateResiduals(t.aob.16.rat.arcsin)) # actually it is better, but i decided to go without transformation
# no change of the anova test and the emmeans test

#### 8. COMAMMOX A/16S RATIO ####

### anova test for Comammox A/16S Ratio Percentage
# Linear mixed model
t.comA.16.rat.percent <- lmerTest::lmer(ComA_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                   (1|block:sampling.date), data=qPCR.BS, na.action=na.omit)
anova(t.comA.16.rat.percent)
# test assumptions:
shapiro.test(resid(t.comA.16.rat.percent)) # normal
plot(simulateResiduals(t.comA.16.rat.percent)) # okay
# Pairwise comparison:
ComA_16.percent.ratio.emm <- emmeans(t.comA.16.rat.percent, ~ irrigation | fertilization*sampling.date)
ComA_16.percent.ratio.pair <- pairs(ComA_16.percent.ratio.emm)
ComA_16.percent.ratio.pair.DF <- as.data.frame(ComA_16.percent.ratio.pair)
### anova test for ComA/16S Ratio Arcsin
t.comA <- lmerTest::lmer(ComA_16.arc.ratio ~ irrigation*fertilization*sampling.date+
                       (1|block:sampling.date), data=qPCR.BS, na.action=na.omit)
anova(t.comA)
t.comA2 <- lmerTest::lmer(ComA_16.arc.ratio ~ irrigation*fertilization*sampling.date+
                       (1|plot), data=qPCR.BS, contrasts = list(irrigation="contr.sum",
                                                                fertilization="contr.sum",sampling.date="contr.sum"))
anova(t.comA2)
ComA_16.arc.ratio.emm2 <- emmeans(t.comA2, ~ irrigation | fertilization*sampling.date)
ComA_16.arc.ratio.pair2 <- pairs(ComA_16.arc.ratio.emm2)
# test assumptions:
shapiro.test(resid(t.comA)) # normal (much better)
plot(simulateResiduals(t.comA)) # okay
# Pairwise comparison:
ComA_16.arc.ratio.emm <- emmeans(t.comA, ~ irrigation | fertilization*sampling.date)
ComA_16.arc.ratio.pair <- pairs(ComA_16.arc.ratio.emm)
ComA_16.arc.ratio.pair.DF <- as.data.frame(summary(ComA_16.arc.ratio.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
#write.csv(ComA_16.arc.ratio.pair.DF, file = "ComA_16_arc_ratio_pair.csv")

#### 9. COMAMMOX B /16S RATIO ####

### anova test for Comammox B /16S Ratio Percentage
# Linear mixed model
t.comB.16.rat.percent <- lmerTest::lmer(ComB_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                   (1|block:sampling.date), data=qPCR.BS, na.action=na.omit)
anova(t.comB.16.rat.percent)
# test assumptions:
shapiro.test(resid(t.comB.16.rat.percent)) # very not normal
plot(simulateResiduals(t.comB.16.rat.percent)) # not okay
#** I need data transformation
### anova test for ComB/16S Ratio Arcsin
t.comB=lmerTest::lmer(ComB_16.arc.ratio ~ irrigation*fertilization*sampling.date+
                        (1|block:sampling.date), data=qPCR.BS,
                        contrasts = list(irrigation="contr.sum",fertilization="contr.sum",sampling.date="contr.sum"))
anova(t.comB)
Anova(t.comB, test="F", type="III")


t.comB2=lmerTest::lmer(ComB_16.arc.ratio ~ irrigation*fertilization*sampling.date+(1|block)+
                        (1|block:sampling.date), data=qPCR.BS,contrasts = list(irrigation="contr.sum",
                                                                fertilization="contr.sum",sampling.date="contr.sum"))
anova(t.comB, t.comB2)
anova(t.comB2)
car::Anova(t.comB2,test="F", type="III")

# test with anova (the data is unbalanced due to the 16S data has one missing value)

# test assumptions:
shapiro.test(resid(t.comB)) # still not normal
plot(simulateResiduals(t.comB)) # much better
# Pairwise comparison:
ComB_16.arc.ratio.emm <- emmeans(t.comB2, ~ irrigation | fertilization*sampling.date)
ComB_16.arc.ratio.pair <- pairs(ComB_16.arc.ratio.emm)
ComB_16.arc.ratio.pair.DF <- as.data.frame(summary(ComB_16.arc.ratio.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
#write.csv(ComB_16.arc.ratio.pair.DF, file = "ComB_16_arc_ratio_pair.csv")


#### 10. AOA/AOB RATIO ####

### anova test for AOA/AOB Ratio without transformation
AOA_AOB_rat <- lmerTest::lmer(AOA_AOB_ratio ~ irrigation*fertilization*sampling.date+
                   (1|block:sampling.date), data=qPCR.BS, na.action=na.omit)
anova(AOA_AOB_rat)
# test assumptions:
shapiro.test(resid(AOA_AOB_rat)) # very not normal
plot(simulateResiduals(AOA_AOB_rat)) # very not okay
#** I need data transformation
### anova test for AOA/AOB Ratio arcsin sqrt transformed
t.AOA_AOB <- lmerTest::lmer(AOA_AOB.arc.ratio ~ irrigation*fertilization*sampling.date+
                          (1|block:sampling.date), data=qPCR.BS, na.action=na.omit)
anova(t.AOA_AOB)
# test assumptions:
shapiro.test(resid(t.AOA_AOB)) # normal
plot(simulateResiduals(t.AOA_AOB)) # slightly not okay
# Pairwise comparison:
AOA_AOB.arc.ratio.emm <- emmeans(t.AOA_AOB, ~ irrigation | fertilization*sampling.date)
AOA_AOB.arc.ratio.pair <- pairs(AOA_AOB.arc.ratio.emm)
AOA_AOB.arc.ratio.pair.DF <- as.data.frame(summary(AOA_AOB.arc.ratio.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
#write.csv(AOA_AOB.arc.ratio.pair.DF, file = "AOA_AOB_arc_ratio_pair.csv")


#### 11. ComA/ComB RATIO ####

### anova test for ComA/ComB Ratio without transformation
t.ComA_ComB <- lmerTest::lmer(ComA_ComB_ratio ~ irrigation*fertilization*sampling.date+
                       (1|block:sampling.date), data=qPCR.BS, na.action=na.omit)
anova(t.ComA_ComB)
# test assumptions:
shapiro.test(resid(t.ComA_ComB)) # normal
plot(simulateResiduals(t.ComA_ComB)) # good
# Pairwise comparison:
ComA_ComB.ratio.emm <- emmeans(t.ComA_ComB, ~ irrigation | fertilization*sampling.date)
ComA_ComB.ratio.pair <- pairs(ComA_ComB.ratio.emm)
ComA_ComB.ratio.pair.DF <- as.data.frame(summary(ComA_ComB.ratio.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil')
#write.csv(ComA_ComB.ratio.pair.DF, file = "ComA_ComB_ratio_pair.csv")


########################################################################################################################################
# Plotting

qPCR.BS$x
qPCR.BS.ed <- qPCR.BS %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")

#### 1. Bulk Soil - AOA Abundance per gram DWS ####

aoa.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,7e+08)+
  ylab(bquote(bold('AOA'~(copies~g^-1~dry~soil))))+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
 facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.position =  "none",
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
aoa.cop.pwc.plot
#setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("AOA_gDWS.BS.tiff",
       aoa.cop.pwc.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


# 2. Bulk Soil - AOB

aob.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOB_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,7e+08)+
  ylab(bquote(bold('AOB'~(copies~g^-1~dry~soil))))+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('AOB'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.position =  "none",
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")#+ggtitle("A. Bulk Soil")
aob.cop.pwc.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
aob.dws.emm.rstat <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(AOB_nbc_per_g_DW_soil  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aob.log.dws)
aob.dws.emm.rstat
# add x y position
aob.log.dws.xy <- aob.dws.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
aob.cop.pwc.plot2 <- aob.cop.pwc.plot + 
  stat_pvalue_manual(aob.log.dws.xy,
                     #step.increase = 1,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)
  #scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
#aob.cop.pwc.plot3 <- aob.cop.pwc.plot2 +ylim(0,3.2e+10)
aob.cop.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("AOB_gDWS.BS.tiff",
       aob.cop.pwc.plot2, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)



# 3. Bulk Soil - Comammox A

comA.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=ComA_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,7e+08)+
  ylab(bquote(bold('Comammox A'~(copies~g^-1~dry~soil))))+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox A'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.position =  "none",
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")#+ggtitle("A. Bulk Soil")
comA.cop.pwc.plot
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("COM_A_gDWS.BS.tiff",
       comA.cop.pwc.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


# 4. Bulk Soil - Comammox B

comB.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=ComB_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,7e+07)+
  ylab(bquote(bold('Comammox B'~(copies~g^-1~dry~soil))))+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.position =  "none",
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")#+ggtitle("A. Bulk Soil")
comB.cop.pwc.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
comB.dws.emm.rstat <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(ComB_nbc_per_g_DW_soil  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = comB.log.dws)
comB.dws.emm.rstat
# add x y position
comB.log.dws.xy <- comB.dws.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
comB.cop.pwc.plot2 <- comB.cop.pwc.plot + 
  stat_pvalue_manual(comB.log.dws.xy,
                     #step.increase = 1,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)
comB.cop.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("COM_B_gDWS.BS.tiff",
       comB.cop.pwc.plot2, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)



# 4. Bulk Soil - 16S

tot.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=Tot_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  ylab(bquote(bold('16S'~(copies~g^-1~dry~soil))))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab('16S')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_alpha_manual(values = c(1, 0.5),
                     #labels=c('Control', 'Drought'),
                     #guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
tot.cop.pwc.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
tot.log.dws.emm.rstat <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(Tot_nbc_per_g_DW_soil  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.16S)
tot.log.dws.emm.rstat 
# add x y position
tot.log.dws.xy <- tot.log.dws.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
tot.cop.pwc.plot2 <- tot.cop.pwc.plot + 
  stat_pvalue_manual(tot.log.dws.xy,
                     #step.increase = 1,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
tot.cop.pwc.plot3 <- tot.cop.pwc.plot2 +ylim(0,3.2e+10)

setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("16Scopies_BS.tiff",
       tot.cop.pwc.plot3, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)



library(patchwork)
copiesdws.All <-  aoa.cop.pwc.plot / aob.cop.pwc.plot / comA.cop.pwc.plot / comB.cop.pwc.plot
copiesdws.All

# 5. Bulk Soil - AOA/16S RATIO

AOA_16S.rat.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_16S_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('AOA/16S (%)')+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
AOA_16S.rat.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
AOA_16S.rat.emm.rstat <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(AOA_16S_ratio_percent  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.aoa)
AOA_16S.rat.emm.rstat 
# add x y position
AOA_16S.rat.emm.xy <- AOA_16S.rat.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
AOA_16S.rat.emm.pwc.plot2 <- AOA_16S.rat.plot + 
  stat_pvalue_manual(AOA_16S.rat.emm.xy,
                     #step.increase = 1,
                     #label = "p.adj.signif",size=3.5,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
AOA_16S.rat.plot3 <- AOA_16S.rat.emm.pwc.plot2 +   ylim(0,10)
AOA_16S.rat.plot3 
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("AOA_16S_ratio_BS.tiff",
       AOA_16S.rat.plot3, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)



#6. Bulk Soil - AOB/16S RATIO

AOB_16S.rat.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOB_16S_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  #ylim(0,8.9)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  ylab('AOB/16S (%)')+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
AOB_16S.rat.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
aob_16S_percent_rat_emm <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(AOB_16S_ratio_percent ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.aob)
aob_16S_percent.xy <- aob_16S_percent_rat_emm %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
AOB_16S.rat.plot2 <- AOB_16S.rat.plot + 
  stat_pvalue_manual(aob_16S_percent.xy,label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
 scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
AOB_16S.rat.plot2
AOB_16S.rat.plot3 <- AOB_16S.rat.plot2 +   ylim(0,5)
AOB_16S.rat.plot3 
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("AOB_16S_ratio_BS.tiff",
       AOB_16S.rat.plot3, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


# 7. Bulk Soil - ComammoxA/16S RATIO

comA_16S.rat.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=ComA_16S_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('Comammox A/16S (%)')+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
comA_16S.rat.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
comA_16S.rat.emm.rstat <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(ComA_16S_ratio_percent  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.comA)
comA_16S.rat.emm.rstat 
# add x y position
comA_16S.rat.emm.xy <- comA_16S.rat.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
comA_16S.rat.emm.pwc.plot2 <- comA_16S.rat.plot + 
  stat_pvalue_manual(comA_16S.rat.emm.xy,
                     #step.increase = 1,
                     #label = "p.adj.signif",size=3.5,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
comA_16S.rat.plot3 <- comA_16S.rat.emm.pwc.plot2 +   ylim(0,1.2)
comA_16S.rat.plot3 
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("comA_16S_ratio_BS.tiff",
       comA_16S.rat.plot3, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


# 8. Bulk Soil - ComammoxB/16S RATIO

comB_16S.rat.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=ComB_16S_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('Comammox B/16S (%)')+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
comB_16S.rat.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
comB_16S.rat.emm.rstat <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(ComB_16S_ratio_percent  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.comB)
comB_16S.rat.emm.rstat 
# add x y position
comB_16S.rat.emm.xy <- comB_16S.rat.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
comB_16S.rat.emm.pwc.plot2 <- comB_16S.rat.plot + 
  stat_pvalue_manual(comB_16S.rat.emm.xy,
                     #step.increase = 1,
                     #label = "p.adj.signif",size=3.5,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
comB_16S.rat.plot3 <- comB_16S.rat.emm.pwc.plot2 +   ylim(0,0.25)
comB_16S.rat.plot3 
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("comB_16S_ratio_BS.tiff",
       comB_16S.rat.plot3, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)

# 9. Bulk Soil - AOA/AOB RATIO

AOA_AOB.rat.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_AOB_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('AOA/AOB (%)')+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
AOA_AOB.rat.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
AOA_AOB.rat.emm.rstat <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(AOA_AOB_ratio_percent  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.AOA_AOB)
AOA_AOB.rat.emm.rstat 
# add x y position
AOA_AOB.rat.emm.xy <- AOA_AOB.rat.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
AOA_AOB.rat.emm.pwc.plot2 <- AOA_AOB.rat.plot + 
  stat_pvalue_manual(AOA_AOB.rat.emm.xy,
                     #step.increase = 1,
                     #label = "p.adj.signif",size=3.5,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
AOA_AOB.rat.plot3 <- AOA_AOB.rat.emm.pwc.plot2 +   ylim(0,500)
AOA_AOB.rat.plot3 
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("AOA_AOB_ratio_BS.tiff",
       AOA_AOB.rat.plot3, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)

# 9. Bulk Soil - comA/comB RATIO

comA_comB.rat.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=ComA_ComB_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('Coma A/Coma B (%)')+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
comA_comB.rat.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
comA_comB.rat.emm.rstat <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(ComA_ComB_ratio_percent  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.ComA_ComB)
comA_comB.rat.emm.rstat 
# add x y position
comA_comB.rat.emm.xy <- comA_comB.rat.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
comA_comB.rat.emm.pwc.plot2 <- comA_comB.rat.plot + 
  stat_pvalue_manual(comA_comB.rat.emm.xy,
                     #step.increase = 1,
                     #label = "p.adj.signif",size=3.5,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
comA_comB.rat.plot3 <- comA_comB.rat.emm.pwc.plot2 +   ylim(0,1650)
comA_comB.rat.plot3 
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("comA_comB_ratio_BS.tiff",
       comA_comB.rat.plot3, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)
#_________________________________________________________________________________
# line chart AOA Copies per gram DWS Bulk Soil 

qPCR.BS.ed

# tidy up the data frame
qPCR.BS.sum <- qPCR.BS.ed %>%
  group_by(irrigation, fertilization,x, sampling.date,var3) %>%
  summarize(AOA_nbc_per_g_DW_soil=mean(AOA_nbc_per_g_DW_soil),
            #sd.aoa.DWS=sd(AOA_nbc_per_g_DW_soil, na.rm = TRUE),
            AOA_nbc_per_ngDNA=mean(AOA_nbc_per_ngDNA),
            #sd.aoa.ngDNA=sd(AOA_nbc_per_ngDNA, na.rm = TRUE),
            AOB_nbc_per_g_DW_soil=mean(AOB_nbc_per_g_DW_soil),
            #sd.aob.DWS=sd(AOB_nbc_per_g_DW_soil, na.rm = TRUE),
            AOB_nbc_per_ngDNA=mean(AOB_nbc_per_ngDNA))
            #sd.aob.ngDNA=sd(AOB_nbc_per_ngDNA, na.rm = TRUE))
qPCR.BS.sum
str(qPCR.BS.sum)
sd.aoa.DWS <- sd(qPCR.BS.ed$AOA_nbc_per_g_DW_soil)



qPCR.BS.ed$lm_pred_val <- predict(aoa.BS.copies.mod,newdata = qPCR.BS.ed,
                              interval =  "confidence"
) %>% as.data.frame()


ggplot(qPCR.BS.ed, aes(x = sampling.date, y = AOA_nbc_per_g_DW_soil, linetype=irrigation)) +
  geom_jitter(position = position_jitter(0.2),
               color = "darkgray") + 
  geom_line(linewidth=1.15,aes(group = x,col=fertilization), data = qPCR.BS.sum) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #geom_errorbar(aes(ymin = AOA_nbc_per_g_DW_soil-sd.aoa.DWS, ymax = AOA_nbc_per_g_DW_soil+sd.aoa.DWS),
                #position = position_dodge(.5),
                #data = qPCR.BS.sum, width = 0.1) +
  geom_point(data = qPCR.BS.sum, size = 2)+
  geom_ribbon(aes(ymin = AOA_nbc_per_g_DW_soil - sd.aoa.DWS, ymax = AOA_nbc_per_g_DW_soil + sd.aoa.DWS,
                  group = x, fill=fertilization),data = qPCR.BS.sum,
              linetype=0, alpha=0.07) +
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  ylab(bquote(~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("A. Bulk Soil")

#_______________________________________________________________________________
# AOA abundance per ng DNA

set.seed(13)
aoa.BS.copies.ngDNA.mod <- lm(qPCR.BS$AOA_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date, data=qPCR.BS)
anova(aoa.BS.copies.ngDNA.mod)

aoa.ngDNA.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = fertilization, alpha=irrigation))+
  theme_bw() +
  labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote(~italic(amoA)~'gene'~(copies~ng^-1~DNA)))+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")
aoa.ngDNA.pwc.plot

#install.packages("patchwork")
#library(patchwork)
aoa.copies.bulk.plot <- aoa.cop.pwc.plot+aoa.ngDNA.pwc.plot + plot_layout(nrow = 2)
aoa.copies.bulk.plot

########################################################################################################################################

#### 1a. Analyses of Rhizosphere - Copies per ng of DNA ####

#setwd('D:/Fina/INRAE_Project/microservices/')
setwd('/Users/arifinabintarti/Documents/France/microservices/')
qPCR <- read.csv("qPCR_results_LP_stat.csv")
qPCR.RS <- qPCR[121:192,]
str(qPCR.RS)
qPCR.RS$nr <- as.factor(qPCR.RS$nr)
qPCR.RS$plot <- as.factor(qPCR.RS$plot)
qPCR.RS$block <- as.factor(qPCR.RS$block)
qPCR.RS$irrigation <- as.factor(qPCR.RS$irrigation)
qPCR.RS$fertilization <- as.factor(qPCR.RS$fertilization)
qPCR.RS$type <- as.factor(qPCR.RS$type)
qPCR.RS$sampling.date <- as.factor(qPCR.RS$sampling.date)
qPCR.RS$var3 <- as.factor(qPCR.RS$var3)
qPCR.RS$x <- as.factor(qPCR.RS$x)
qPCR.RS$rep <- as.factor(qPCR.RS$rep)
qPCR.RS$sampling.date <- factor(qPCR.RS$sampling.date, levels = c("28/04/2022", "1/6/22", "5/7/22"),
                                labels = c("Apr 28th", "Jun 1st", "Jul 5th"))

##subset
aoa.rh.M <- qPCR.RS[which(qPCR.RS$fertilization == "M"),]
aoa.rh.M
aoa.rh.D <- qPCR.RS[which(qPCR.RS$fertilization == "D"),]
aoa.rh.D
aoa.rh.K <- qPCR.RS[which(qPCR.RS$fertilization == "K"),]
aoa.rh.K
aoa.rh.D1 <- qPCR.RS[which(qPCR.RS$fertilization == "D"
                            & qPCR.RS$sampling.date == "04-28-22" ),]
aoa.rh.D2 <- qPCR.RS[which(qPCR.RS$fertilization == "D"
                            & qPCR.RS$sampling.date == "06-01-22" ),]
aoa.rh.D3 <- qPCR.RS[which(qPCR.RS$fertilization == "D"
                            & qPCR.RS$sampling.date == "07-05-22" ),]
aoa.rh.K1 <- qPCR.RS[which(qPCR.RS$fertilization == "K"
                           & qPCR.RS$sampling.date == "04-28-22" ),]
aoa.rh.K2 <- qPCR.RS[which(qPCR.RS$fertilization == "K"
                           & qPCR.RS$sampling.date == "06-01-22" ),]
aoa.rh.K3 <- qPCR.RS[which(qPCR.RS$fertilization == "K"
                           & qPCR.RS$sampling.date == "07-05-22" ),]
aoa.rh.M1 <- qPCR.RS[which(qPCR.RS$fertilization == "M"
                           & qPCR.RS$sampling.date == "04-28-22" ),]
aoa.rh.M2 <- qPCR.RS[which(qPCR.RS$fertilization == "M"
                           & qPCR.RS$sampling.date == "06-01-22" ),]
aoa.rh.M3 <- qPCR.RS[which(qPCR.RS$fertilization == "M"
                           & qPCR.RS$sampling.date == "07-05-22" ),]


# perform log transformation
qPCR.RS$AOA_logDNA <- log10(qPCR.RS$AOA_nbc_per_ngDNA)
qPCR.RS$AOB_logDNA <- log10(qPCR.RS$AOB_nbc_per_ngDNA)
qPCR.RS$ComA_logDNA <- log10(qPCR.RS$ComA_nbc_per_ngDNA)
qPCR.RS$ComB_logDNA <- log10(qPCR.RS$ComB_nbc_per_ngDNA)
qPCR.RS$Tot_logDNA <- log10(qPCR.RS$Tot_nbc_per_ngDNA)
# perform arcsin root square
qPCR.RS$AOA_16.arc.ratio.rh <- asin(sqrt(qPCR.RS$AOA_16S_ratio_percent / 100))
qPCR.RS$AOB_16.arc.ratio.rh <- asin(sqrt(qPCR.RS$AOB_16S_ratio_percent / 100))
qPCR.RS$ComA_16.arc.ratio.rh <- asin(sqrt(qPCR.RS$ComA_16S_ratio_percent / 100))
qPCR.RS$ComB_16.arc.ratio.rh <- asin(sqrt(qPCR.RS$ComB_16S_ratio_percent / 100))
qPCR.RS$AOA_AOB.arc.ratio.rh <- asin(sqrt(qPCR.RS$AOA_AOB_ratio / 100))
qPCR.RS$ComA_ComB.arc.ratio.rh <- asin(sqrt(qPCR.RS$ComA_ComB_ratio / 100))


##### 1. AOA in Rhizosphere

# Anova test for transformed AOA nbc per ngDNA Log transformed Data

aoa.log.dna.aov <- aov_ez("plot", "AOA_logDNA", qPCR.RS, 
                          within = "sampling.date",
                          between = c("fertilization","irrigation"),
                          type = 2,
                          return = afex_options("return_aov"),
                          anova_table = list(correction="none"))
aoa.log.dna.aov

# Three-Way Mixed (Split-Plot) ANOVA 
aoa.log.dna.aov2 <- anova_test(
  data = qPCR.RS, dv = AOA_logDNA, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aoa.log.dna.aov2)
aoa.log.dna.aov2.df <- get_anova_table(aoa.log.dna.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOA_log_DNA")
table2csv(x=aoa.log.dna.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(aoa.log.dna.aov) # good
check_sphericity(aoa.log.dna.aov) # good
aoa.log.dna.is_norm <- check_normality(aoa.log.dna.aov)
shapiro.test(qPCR.RS$AOA_logDNA) # NORMAL
shapiro.test(qPCR.RS$AOA_nbc_per_ngDNA) # NORMAL
plot(aoa.log.dna.is_norm, type = "qq")
plot(aoa.log.dna.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOA_logDNA", ggtheme = theme_bw()) # better
ggqqplot(qPCR.RS, "AOA_nbc_per_ngDNA", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(aoa.log.dna.aov))
# Pairwise comparison:
aoa.log.dna.emm <- emmeans(aoa.log.dna.aov, ~ irrigation | fertilization*sampling.date)
aoa.log.dna.pair <- pairs(aoa.log.dna.emm)
aoa.log.dna.pair.DF <- as.data.frame(summary(aoa.log.dna.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(aoa.log.dna.pair.DF, file = "AOA_log_DNA_pair.csv")



# Anova test for non transformed AOA nbc per ngDNA 

aoa.ngDNA.aov <- aov_ez("plot", "AOA_nbc_per_ngDNA", qPCR.RS, 
                          within = "sampling.date",
                          between = c("fertilization","irrigation"),
                          type = 2,
                          return = afex_options("return_aov"),
                          anova_table = list(correction="none"))
aoa.ngDNA.aov

# Three-Way Mixed (Split-Plot) ANOVA 
aoa.ngDNA.aov2 <- anova_test(
  data = qPCR.RS, dv = AOA_nbc_per_ngDNA, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aoa.ngDNA.aov2)
aoa.ngDNA.aov2.df <- get_anova_table(aoa.ngDNA.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOA_ngDNA")
table2csv(x=aoa.ngDNA.aov2.df,file=filen, digits = 1, digitspvals = 3)

#another model
t.aoa.rh=lmerTest::lmer(AOA_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date+
                       (1|block:sampling.date), 
                     data=qPCR.RS, na.action=na.omit)

anova(t.aoa.rh)


# test assumptions:
check_homogeneity(aoa.ngDNA.aov) # good
check_sphericity(aoa.ngDNA.aov) # good
aoa.ngDNA.is_norm <- check_normality(aoa.ngDNA.aov)
shapiro.test(qPCR.RS$AOA_nbc_per_ngDNA) # NORMAL
plot(aoa.ngDNA.is_norm, type = "qq")
plot(aoa.ngDNA.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOA_nbc_per_ngDNA", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(aoa.ngDNA.aov))
# Pairwise comparison:
aoa.ngDNA.emm <- emmeans(aoa.ngDNA.aov, ~ irrigation | fertilization*sampling.date)
aoa.ngDNA.pair <- pairs(aoa.ngDNA.emm)
aoa.ngDNA.pair.DF <- as.data.frame(summary(aoa.ngDNA.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(aoa.ngDNA.pair.DF, file = "AOA_ngDNA_pair.csv")


##### 2. AOB in Rhizosphere

# Anova test for transformed AOA nbc per ngDNA Log transformed Data

aob.log.dna.aov <- aov_ez("plot", "AOB_logDNA", qPCR.RS, 
                          within = "sampling.date",
                          between = c("irrigation","fertilization"),
                          type = 2,
                          return = afex_options("return_aov"),
                          anova_table = list(correction= "none"))
aob.log.dna.aov

# Three-Way Mixed (Split-Plot) ANOVA 
aob.log.dna.aov2 <- anova_test(
  data = qPCR.RS, dv = AOB_logDNA, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aob.log.dna.aov2)
aob.log.dna.aov2.df <- get_anova_table(aob.log.dna.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOB_log_DNA")
table2csv(x=aob.log.dna.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(aob.log.dna.aov) # good
check_sphericity(aob.log.dna.aov) # good
aob.log.dna.is_norm <- check_normality(aob.log.dna.aov)
shapiro.test(qPCR.RS$AOB_logDNA) # NORMAL
shapiro.test(qPCR.RS$AOB_nbc_per_ngDNA) # NORMAL
plot(aob.log.dna.is_norm, type = "qq")
plot(aob.log.dna.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOB_logDNA", ggtheme = theme_bw()) # better
ggqqplot(qPCR.RS, "AOB_nbc_per_ngDNA", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(aob.log.dna.aov))
# Pairwise comparison:
aob.log.dna.emm <- emmeans(aob.log.dna.aov, ~ irrigation | fertilization*sampling.date)
aob.log.dna.pair <- pairs(aob.log.dna.emm)
aob.log.dna.pair.DF <- as.data.frame(summary(aob.log.dna.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(aob.log.dna.pair.DF, file = "AOB_log_DNA_pair.csv")



# Anova test for non transformed AOA nbc per ngDNA 

aob.ngDNA.aov <- aov_ez("plot", "AOB_nbc_per_ngDNA", qPCR.RS, 
                        within = "sampling.date",
                        between = c("irrigation","fertilization"),
                        type = 2,
                        return = afex_options("return_aov"),
                        anova_table = list(correction="none"))
aob.ngDNA.aov

# Three-Way Mixed (Split-Plot) ANOVA 
aob.ngDNA.aov2 <- anova_test(
  data = qPCR.RS, dv = AOB_nbc_per_ngDNA, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aob.ngDNA.aov2)
aob.ngDNA.aov2.df <- get_anova_table(aob.ngDNA.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOB_ngDNA")
table2csv(x=aob.ngDNA.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(aob.ngDNA.aov) # good
check_sphericity(aob.ngDNA.aov) # good
aob.ngDNA.is_norm <- check_normality(aob.ngDNA.aov)
shapiro.test(qPCR.RS$AOB_nbc_per_ngDNA) # NORMAL
plot(aob.ngDNA.is_norm, type = "qq")
plot(aob.ngDNA.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOB_nbc_per_ngDNA", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(aob.ngDNA.aov))
# Pairwise comparison:
aob.ngDNA.emm <- emmeans(aob.ngDNA.aov, ~ irrigation | fertilization*sampling.date)
aob.ngDNA.pair <- pairs(aob.ngDNA.emm)
aob.ngDNA.pair.DF <- as.data.frame(summary(aob.ngDNA.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(aob.ngDNA.pair.DF, file = "AOB_ngDNA_pair.csv")


##### 3. Comammox A in Rhizosphere

# Anova test for non transformed Comammox A nbc per ngDNA 

ComA.ngDNA.aov <- aov_ez("plot", "ComA_nbc_per_ngDNA", qPCR.RS, 
                        within = "sampling.date",
                        between = c("irrigation","fertilization"),
                        correction_aov = "GG",
                        type = 3,
                        return = afex_options("return_aov"),
                        anova_table = list())
ComA.ngDNA.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComA.ngDNA.aov2 <- anova_test(
  data = qPCR.RS, dv = ComA_nbc_per_ngDNA, wid = plot, type = 3,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(ComA.ngDNA.aov2)
ComA.ngDNA.aov2.df <- get_anova_table(ComA.ngDNA.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/ComA_ngDNA")
table2csv(x=ComA.ngDNA.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComA.ngDNA.aov) # SLIGHTLY NOT HOMOGEN
check_sphericity(ComA.ngDNA.aov) # good
ComA.ngDNA.is_norm <- check_normality(ComA.ngDNA.aov)
shapiro.test(qPCR.RS$ComA_nbc_per_ngDNA) # slightly NOT NORMAL
plot(ComA.ngDNA.is_norm, type = "density")
plot(ComA.ngDNA.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "ComA_nbc_per_ngDNA", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(ComA.ngDNA.aov))
# Pairwise comparison:
ComA.ngDNA.emm <- emmeans(ComA.ngDNA.aov, ~ irrigation | fertilization*sampling.date)
ComA.ngDNA.pair <- pairs(ComA.ngDNA.emm)
ComA.ngDNA.pair.DF <- as.data.frame(summary(ComA.ngDNA.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(ComA.ngDNA.pair.DF, file = "ComA_ngDNA_pair.csv")


# Anova test for log transformed Comammox A nbc per ngDNA 

ComA.log.dna.aov <- aov_ez("plot", "ComA_logDNA", qPCR.RS, 
                         within = "sampling.date",
                         between = c("irrigation","fertilization"),
                         type = 3,
                         return = afex_options("return_aov"),
                         anova_table = list())
ComA.log.dna.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComA.log.dna.aov2 <- anova_test(
  data = qPCR.RS, dv = ComA_logDNA, wid = plot, type = 3,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(ComA.log.dna.aov2)
ComA.log.dna.aov2.df <- get_anova_table(ComA.log.dna.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/ComA_log_DNA")
table2csv(x=ComA.log.dna.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComA.log.dna.aov) # VERY NOT HOMOGEN
check_sphericity(ComA.log.dna.aov) # good
ComA.log.dna.is_norm <- check_normality(ComA.log.dna.aov)
shapiro.test(qPCR.RS$ComA_logDNA) # NORMAL
plot(ComA.log.dna.is_norm, type = "density")
plot(ComA.log.dna.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "ComA_logDNA", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(ComA.log.dna.aov))
# Pairwise comparison:
ComA.log.dna.emm <- emmeans(ComA.log.dna.aov, ~ irrigation | fertilization*sampling.date)
ComA.log.dna.pair <- pairs(ComA.log.dna.emm)
ComA.log.dna.pair.DF <- as.data.frame(summary(ComA.log.dna.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(ComA.log.dna.pair.DF, file = "ComA_log_DNA_pair.csv")

##### 4. Comammox B in Rhizosphere

# Anova test for non transformed Comammox B nbc per ngDNA 

ComB.ngDNA.aov <- aov_ez("plot", "ComB_nbc_per_ngDNA", qPCR.RS, 
                         within = "sampling.date",
                         between = c("irrigation","fertilization"),
                         type = 2,
                         return = afex_options("return_aov"),
                         anova_table = list(correction="none"))
ComB.ngDNA.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComB.ngDNA.aov2 <- anova_test(
  data = qPCR.RS, dv = ComB_nbc_per_ngDNA, wid = plot, type = 3,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(ComB.ngDNA.aov2)
ComB.ngDNA.aov2.df <- get_anova_table(ComB.ngDNA.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/ComB_ngDNA")
table2csv(x=ComB.ngDNA.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComB.ngDNA.aov) # good
check_sphericity(ComB.ngDNA.aov) # good
ComB.ngDNA.is_norm <- check_normality(ComB.ngDNA.aov)
shapiro.test(qPCR.RS$ComB_nbc_per_ngDNA) # NORMAL
plot(ComB.ngDNA.is_norm, type = "qq")
plot(ComB.ngDNA.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "ComB_nbc_per_ngDNA", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(ComB.ngDNA.aov))
# Pairwise comparison:
ComB.ngDNA.emm <- emmeans(ComB.ngDNA.aov, ~ irrigation | fertilization*sampling.date)
ComB.ngDNA.pair <- pairs(ComB.ngDNA.emm)
ComB.ngDNA.pair.DF <- as.data.frame(summary(ComB.ngDNA.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(ComB.ngDNA.pair.DF, file = "ComB_ngDNA_pair.csv")


#### 5. 16S in Rhizosphere ####

# linear mixed model for non transformed 16S nbc per ngDNA 
#t.tot.rh1=lmerTest::lmer(Tot_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date+
                          #(1|block/plot),data=qPCR.RS, na.action=na.omit, REML=F)
t.tot.rh2 <- lmerTest::lmer(Tot_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date+
                          (sampling.date|block), data=qPCR.RS, na.action=na.omit)
anova(t.tot.rh2)
# test assumption
shapiro.test(resid(t.tot.rh2)) # normal
plot(simulateResiduals(t.tot.rh2)) # very good
# ** Transformation is not needed
# Pairwise comparison:
Tot.ngDNA.emm <- emmeans(t.tot.rh2, ~ irrigation | fertilization*sampling.date)
Tot.ngDNA.pair <- pairs(Tot.ngDNA.emm)
Tot.ngDNA.pair.DF <- as.data.frame(summary(Tot.ngDNA.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
#write.csv(Tot.ngDNA.pair.DF, file = "16S_ngDNA_pair.csv")


#### 6. AOA/16S Ratio in Rhizosphere ####


# linear mixed model test for AOA/16S Percent Ratio
t.aoa_16S_percent.rh <- lmerTest::lmer(AOA_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                        (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)

anova(t.aoa_16S_percent.rh)
# test assumption
shapiro.test(resid(t.aoa_16S_percent.rh)) # normal
plot(simulateResiduals(t.aoa_16S_percent.rh)) # okay
#*** No need to do transformation
# Pairwise comparison:
AOA_16S.percent.rat.rh.emm <- emmeans(t.aoa_16S_percent.rh, ~ irrigation | fertilization*sampling.date)
AOA_16S.percent.rat.rh.pair <- pairs(AOA_16S.percent.rat.rh.emm)
AOA_16S.percent.rat.rh.DF <- as.data.frame(summary(AOA_16S.percent.rat.rh.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
#write.csv(AOA_16S.percent.rat.rh.DF, file = "AOA_16S.percent.rat.pair.rh.csv")
#Just checking
# linear mixed model test for AOA/16S Arcsin 
t.aoa_16S.rh <- lmerTest::lmer(AOA_16.arc.ratio.rh ~ irrigation*fertilization*sampling.date+
                (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)

anova(t.aoa_16S.rh)
# test assumption
shapiro.test(resid(t.aoa_16S.rh)) # normal
plot(simulateResiduals(t.aoa_16S.rh)) # okay


##### 7. AOB/16S Ratio in Rhizosphere


# linear mixed model test for AOB/16S Percent Ratio
t.aob_16S_percent.rh <- lmerTest::lmer(AOB_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                        (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)
anova(t.aob_16S_percent.rh)
# test assumption
shapiro.test(resid(t.aob_16S_percent.rh)) # normal
plot(simulateResiduals(t.aob_16S_percent.rh)) # good
#*** No need to do transformation
# Pairwise comparison:
AOB_16S.percent.rat.rh.emm <- emmeans(t.aob_16S_percent.rh, ~ irrigation | fertilization*sampling.date)
AOB_16S.percent.rat.rh.pair <- pairs(AOB_16S.percent.rat.rh.emm)
AOB_16S.percent.rat.rh.DF <- as.data.frame(summary(AOB_16S.percent.rat.rh.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
#write.csv(AOB_16S.percent.rat.rh.DF, file = "AOB_16S.percent.rat.pair.rh.csv")
# Just Checking
# linear mixed model test for AOB/16S Arcsin 
t.aob_16S_arcs.rh <- lmerTest::lmer(AOB_16.arc.ratio.rh ~ irrigation*fertilization*sampling.date+
                     (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)
anova(t.aob_16S_arcs.rh)
# test assumption
shapiro.test(resid(t.aob_16S_arcs.rh)) # normal
plot(simulateResiduals(t.aob_16S_arcs.rh)) # good


#### 7. ComA/16S Ratio in Rhizosphere ####


# Linear mixed model for ComA/16S Ratio in percent
t.comA_16S_percent.rh <- lmerTest::lmer(ComA_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                         (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)
anova(t.comA_16S_percent.rh)
0# test assumption
shapiro.test(resid(t.comA_16S_percent.rh)) # normal
plot(simulateResiduals(t.comA_16S_percent.rh)) # good
#*** No need to do transformation
# Pairwise comparison:
ComA_16S.ratio.rh.emm <- emmeans(t.comA_16S_percent.rh, ~ irrigation | fertilization*sampling.date)
ComA_16S.ratio.rh.pair <- pairs(ComA_16S.ratio.rh.emm)
ComA_16S.ratio.rh.DF <- as.data.frame(summary(ComA_16S.ratio.rh.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
#write.csv(ComA_16S.ratio.rh.DF, file = "ComA_16S_percent_ratio.pair.rh.csv")
# Just Checking
# linear mixed model for arcsin square root transformed ComA/16S Ratio
t.comA_16S_arcs.rh <- lmerTest::lmer(ComA_16.arc.ratio.rh ~ irrigation*fertilization*sampling.date+
                      (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)

anova(t.comA_16S_arcs.rh)
# test assumption
shapiro.test(resid(t.comA_16S_arcs.rh)) # normal
plot(simulateResiduals(t.comA_16S_arcs.rh)) # good


#### 8. ComB/16S Ratio in Rhizosphere ####

# Linear mixed model test for ComB/16S Ratio in percent
t.comB_16S_percent.rh <- lmerTest::lmer(ComB_16S_ratio_percent ~ irrigation*fertilization*sampling.date+
                         (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)
anova(t.comB_16S_percent.rh)
# test assumption
shapiro.test(resid(t.comB_16S_percent.rh)) # normal
plot(simulateResiduals(t.comB_16S_percent.rh)) # good
#*** No need to do data transformation
# Pairwise comparison:
ComB_16S.ratio.rh.emm <- emmeans(t.comB_16S_percent.rh, ~ irrigation | fertilization*sampling.date)
ComB_16S.ratio.rh.pair <- pairs(ComB_16S.ratio.rh.emm)
ComB_16S.ratio.rh.DF <- as.data.frame(summary(ComB_16S.ratio.rh.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
#write.csv(ComB_16S.ratio.rh.DF, file = "ComB_16S_percent_ratio.pair.rh.csv")
# Just Checking
# Linear mixed model test for ComB/16S Arcsin Ratio 
t.comB_16S_arc.rh <- lmerTest::lmer(ComB_16.arc.ratio.rh ~ irrigation*fertilization*sampling.date+
                         (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)

anova(t.comB_16S_arc.rh)
# test assumption
shapiro.test(resid(t.comB_16S_arc.rh)) # normal
plot(simulateResiduals(t.comB_16S_arc.rh)) # good


#### 9. AOA/AOB Ratio in Rhizosphere ####

# Linear mixed model test for AOA/AOB Ratio on non transformed data
t.AOA_AOB_percent.rh <- lmerTest::lmer(AOA_AOB_ratio ~ irrigation*fertilization*sampling.date+
                        (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)

anova(t.AOA_AOB_percent.rh)
# test assumption
shapiro.test(resid(t.AOA_AOB_percent.rh)) # normal
plot(simulateResiduals(t.AOA_AOB_percent.rh)) # okay
#***No need data transformation
# Pairwise comparison:
AOA_AOB.ratio.rh.emm <- emmeans(t.AOA_AOB_percent.rh, ~ irrigation | fertilization*sampling.date)
AOA_AOB.ratio.rh.pair <- pairs(AOA_AOB.ratio.rh.emm)
AOA_AOB.ratio.rh.DF <- as.data.frame(summary(AOA_AOB.ratio.rh.pair))
AOA_AOB.ratio.rh.DF
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
#write.csv(AOA_AOB.ratio.rh.DF, file = "AOA_AOB_ratio.pair.rh.csv")
# Just Checking
# Linear mixed model test for AOA/AOB Ratio on arcsin square root transformed data
t.AOA_AOB_arcsin.rh <- lmerTest::lmer(AOA_AOB.arc.ratio.rh ~ irrigation*fertilization*sampling.date+
                       (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)
anova(t.AOA_AOB_arcsin.rh)
# test assumption
shapiro.test(resid(t.AOA_AOB_arcsin.rh)) # normal
plot(simulateResiduals(t.AOA_AOB_arcsin.rh)) # okay


#### 9. ComA/ComB Ratio in Rhizosphere ####


# Linear mixed model test for ComA/ComB Ratio on non transformed data
t.ComA_ComB_percent.rh <- lmerTest::lmer(ComA_ComB_ratio ~ irrigation*fertilization*sampling.date+
                          (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)
anova(t.ComA_ComB_percent.rh)
# test assumption
shapiro.test(resid(t.ComA_ComB_percent.rh)) # not normal
plot(simulateResiduals(t.ComA_ComB_percent.rh)) # looks good actually
#*** No need to do data transsformation
# Pairwise comparison:
ComA_ComB.ratio.rh.emm <- emmeans(t.ComA_ComB_percent.rh, ~ irrigation | fertilization*sampling.date)
ComA_ComB.ratio.rh.pair <- pairs(ComA_ComB.ratio.rh.emm)
ComA_ComB.ratio.rh.DF <- as.data.frame(summary(ComA_ComB.ratio.rh.pair))
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
#write.csv(ComA_ComB.ratio.rh.DF, file = "ComA_ComB_ratio.pair.rh.csv")
# Just checking
# linear mixed model test for ComA_ComB Ratio on arcsin square root transformed data
t.ComA_ComB_arcsin.rh <- lmerTest::lmer(ComA_ComB.arc.ratio.rh ~ irrigation*fertilization*sampling.date+
                         (1|block:sampling.date), data=qPCR.RS, na.action=na.omit)
anova(t.ComA_ComB_arcsin.rh)
# test assumption
shapiro.test(resid(t.ComA_ComB_arcsin.rh)) # still not normal
plot(simulateResiduals(t.ComA_ComB_arcsin.rh)) # looks good 





#########################################################################################################################################
# check stats
aoa.RS.copies <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  get_summary_stats(AOA_nbc_per_ngDNA, type = "mean_sd")

aoa.RS.copies.plot <- ggboxplot(
  qPCR.RS, x = "irrigation", y = "AOA_nbc_per_ngDNA",
  color = "fertilization", palette = "jco",
  facet.by =  "sampling.date")
aoa.RS.copies.plot

# check assumption (outliers)
aoa.RS.copies.out <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  identify_outliers(AOA_nbc_per_ngDNA) # no extreme outliers
# Saphiro-Wilk for normality
aoa.RS.copies.SW <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(AOA_nbc_per_ngDNA)
ggqqplot(qPCR.RS, "AOA_nbc_per_ngDNA", ggtheme = theme_bw()) +
  facet_grid(sampling.date ~ fertilization, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.RS.copies.Lave <- qPCR.RS %>%
  group_by(sampling.date) %>%
  levene_test(AOA_nbc_per_ngDNA ~ irrigation*fertilization)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aoa.RS.copies.aov <- anova_test(
  data = qPCR.RS, dv = AOA_nbc_per_ngDNA, wid = nr ,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aoa.RS.copies.aov)

# TWO-WAY
res.aov2.rh <- aov(AOA_nbc_per_ngDNA ~ irrigation + fertilization, data = qPCR.RS)
summary(res.aov2.rh)
# THREE WAY
res.aov3.rh <- qPCR.RS %>% anova_test(AOA_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date)
res.aov3.rh
# Model Fit
set.seed(13)
aoa.RS.copies.lm.mod <- lm(qPCR.RS$AOA_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date, data=qPCR.RS)
anova(aoa.RS.copies.lm.mod) # similar to type II anova in three-way 
# Pairwise comparisons
aoa.rh.K %>%
  pairwise_t_test(AOA_nbc_per_ngDNA ~ irrigation, p.adjust.method = "BH")
# Pairwise comparisons
library(emmeans)
pwc.aoa.ngDNA.rhizo <- qPCR.RS %>%
  group_by(fertilization, sampling.date) %>%
  emmeans_test(AOA_nbc_per_ngDNA ~ irrigation, p.adjust.method = "BH") %>%
  select(-df, -statistic, -p) # Remove details
pwc.aoa.ngDNA.rhizo # significant if only i grouped by fertilization not significant per sampling date


qPCR.RS$x
qPCR.RS.ed <- qPCR.RS %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")
# plot
aoa.ngDNA.rhizo.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=AOA_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = fertilization, alpha=irrigation))+
  theme_bw() +
  labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote('amoA gene'~(copies~ng^-1~DNA)))+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ ggtitle("B. Rhizosphere")
aoa.ngDNA.rhizo.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_copies.Rhizo.tiff",
       aoa.ngDNA.rhizo.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


aoa.copies.all.plot <- aoa.cop.pwc.plot+aoa.ngDNA.rhizo.plot + plot_layout(nrow = 2)
aoa.copies.all.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_copies.ALL.tiff",
       aoa.copies.all.plot, device = "tiff",
       width = 12, height =11, 
       units= "in", dpi = 600)


############################################################################################################################
# 2. Response variable: AOB COPY NUMBER
############################################################################################################################
# 1. Bulk Soil

aob.BS.copies <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  get_summary_stats(AOB_nbc_per_g_DW_soil, type = "mean_sd")


aob.BS.copies.plot <- ggboxplot(
  qPCR.BS, x = "irrigation", y = "AOB_nbc_per_g_DW_soil",
  color = "fertilization", palette = "jco",
  facet.by =  "sampling.date")
aob.BS.copies.plot

# check assumption (outliers)
aob.BS.copies.out <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  identify_outliers(AOB_nbc_per_g_DW_soil) # no extreme outliers
# Saphiro-Wilk for normality
aob.BS.copies.SW <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(AOB_nbc_per_g_DW_soil)
ggqqplot(qPCR.BS, "AOB_nbc_per_g_DW_soil", ggtheme = theme_bw()) +
  facet_grid(sampling.date ~ fertilization, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.BS.copies.Lave <- qPCR.BS %>%
  group_by(sampling.date) %>%
  levene_test(AOB_nbc_per_g_DW_soil ~ irrigation*fertilization)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aob.BS.copies.aov <- anova_test(
  data = qPCR.BS, type=3, dv = AOB_nbc_per_g_DW_soil, wid = plot,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aob.BS.copies.aov)

aob.bulk.DW <- aov(AOB_nbc_per_g_DW_soil~irrigation, qPCR.BS)
anova(aob.bulk.DW)




# Model Fit
set.seed(13)
aob.BS.copies.mod <- lm(qPCR.BS$AOB_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date, data=qPCR.BS)
anova(aob.BS.copies.mod)

set.seed(13)
aob.BS.copies.mod.lmer <- lmerTest::lmer(qPCR.BS$AOB_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date+(1|plot), data=qPCR.BS)
anova(aob.BS.copies.mod.lmer)

# plot
aob.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOB_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = fertilization, alpha=irrigation))+
  theme_bw() +
  labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote('amoA gene'~(copies~g^-1~dry~soil)))+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
aob.cop.pwc.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_copies.Bulk.tiff",
       aob.cop.pwc.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)

#______________________________________________________________________________
# 2. Rhizosphere

# check assumption (outliers)
aob.RS.copies.out <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  identify_outliers(AOB_nbc_per_ngDNA) # no extreme outliers
# Saphiro-Wilk for normality
aob.RS.copies.SW <- qPCR.RS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(AOB_nbc_per_ngDNA)
ggqqplot(qPCR.RS, "AOB_nbc_per_ngDNA", ggtheme = theme_bw()) +
  facet_grid(sampling.date ~ fertilization, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.RS.copies.Lave <- qPCR.RS %>%
  group_by(sampling.date) %>%
  levene_test(AOB_nbc_per_ngDNA ~ irrigation*fertilization)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aob.RS.copies.aov <- anova_test(
  data = qPCR.RS, dv = AOB_nbc_per_ngDNA, wid = plot,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aob.RS.copies.aov)
# Model Fit


set.seed(13)
aob.RS.copies.lm.mod <- lm(qPCR.RS$AOB_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date, data=qPCR.RS)
anova(aob.RS.copies.lm.mod)


aob.RS.copies.mod <- lmerTest::lmer(qPCR.RS$AOB_nbc_per_ngDNA ~ irrigation*fertilization + (1|plot), data=qPCR.RS)
anova(aob.RS.copies.mod)


aob.test.pw.rhizo <- qPCR.RS %>%
  group_by(fertilization, sampling.date) %>%
  emmeans_test(AOB_nbc_per_ngDNA ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95)
aob.test.pw.rhizo

# plot
aob.ngDNA.rhizo.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=AOB_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = fertilization, alpha=irrigation))+
  theme_bw() +
  labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote('amoA gene'~(copies~ng^-1~DNA)))+
  scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ ggtitle("B. Rhizosphere")
aob.ngDNA.rhizo.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_copies.Rhizo.tiff",
       aob.ngDNA.rhizo.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


aob.copies.all.plot <- aob.cop.pwc.plot+aob.ngDNA.rhizo.plot + plot_layout(nrow = 2)
aob.copies.all.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOB_copies.ALL.tiff",
       aob.copies.all.plot, device = "tiff",
       width = 12, height =11, 
       units= "in", dpi = 600)