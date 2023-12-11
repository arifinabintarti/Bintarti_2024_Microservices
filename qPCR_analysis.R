#############################################################################################
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
library(phyloseq)
library(datarium)
library(rstatix)
library(export)

###########################################################################
# 1. Response variable: BACTERIAL AND ARCHAEAL COPY NUMBER
###########################################################################

#### 1a. Analyses of Bulk Soil - Copies per Gram Dry weight of Soil ####
#setwd('/Users/arifinabintarti/Documents/France/microservices/')
setwd('D:/Fina/INRAE_Project/microservices/')
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


# Check ANOVA assumptions: 

# 1. AOA 

# Boxplot
qPCR.BS.AOA_nbc_dw.plot <- ggboxplot(
  qPCR.BS, x = "sampling.date", y = "AOA_nbc_per_g_DW_soil",
  color = "irrigation", palette = "jco",
  facet.by =  "fertilization")
qPCR.BS.AOA_nbc_dw.plot
# check assumption (outliers)
qPCR.BS.AOA_nbc_dws.out <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  identify_outliers(AOA_logDWS) # no extreme outliers
# Saphiro-Wilk for normality
qPCR.BS.AOA_nbc_dws.SW <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(AOA_nbc_per_g_DW_soil)

ggqqplot(qPCR.BS, "AOA_nbc_per_g_DW_soil", ggtheme = theme_bw())+
  facet_grid(sampling.date ~ irrigation*fertilization, labeller = "label_both")#All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
ggqqplot(qPCR.BS, "AOA_nbc_per_g_DW_soil", ggtheme = theme_bw())
# Lavene test
qPCR.BS.AOA_nbc_dws.Lave <- qPCR.BS %>%
  group_by(fertilization,sampling.date) %>%
  levene_test(AOA_nbc_per_g_DW_soil ~ irrigation)


# Three-Way Mixed (Split-Plot) ANOVA 
AOA_nbc_dws.aov <- anova_test(
  data = qPCR.BS, dv = AOA_nbc_per_g_DW_soil, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(AOA_nbc_dws.aov)

#install.packages("afex") 
#install.packages("qqplotr")
library(afex)
library(performance)
library(qqplotr)
library(emmeans)
library(reshape2)

data(obk.long, package = "afex")

o1 <- aov_ez("plot", "AOA_nbc_per_g_DW_soil", qPCR.BS, 
             within = "sampling.date",
             between = c("fertilization","irrigation"),
             type = 2,
             return = afex_options("return_aov"),
             anova_table = list())
o1
check_homogeneity(o1)
check_sphericity(o1)
is_norm <- check_normality(o1)
shapiro.test(qPCR.BS$AOA_nbc_per_g_DW_soil) # NOT NORMAL
plot(is_norm)
plot(is_norm, type = "qq")
plot(is_norm, type = "qq", detrend = TRUE)
knitr::kable(nice(o1))
m1 <- emmeans(o1, ~ irrigation | fertilization*sampling.date)
pairs(m1)
#pwc <- test(pairs(m1), by=NULL, adjust="bh")
#pairs(m1, simple="each",adjust = "bh")

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

#qPCR.BS$AOA_logngDNA <- log10(qPCR.BS$AOA_nbc_per_ngDNA)
#qPCR.BS$AOB_logngDNA <- log10(qPCR.BS$AOB_nbc_per_ngDNA)

# anova test for transformed AOA DWS
aoa.log.dws.aov <- aov_ez("plot", "AOA_logDWS", qPCR.BS, 
             within = "sampling.date",
             between = c("fertilization","irrigation"),
             type = 2,
             return = afex_options("return_aov"),
             anova_table = list(correction = "none"))
aoa.log.dws.aov

# Three-Way Mixed (Split-Plot) ANOVA 
aoa.log.dws.aov2 <- anova_test(
  data = qPCR.BS, dv = AOA_logDWS, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aoa.log.dws.aov2)
aoa.log.dws.aov2.df <- get_anova_table(aoa.log.dws.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil/AOA_log_dws")
table2csv(x=aoa.log.dws.aov2.df,file=filen, digits = 1, digitspvals = 3)


# try another model
set.seed(13)
aoa.log.dws.lmer <- lmerTest::lmer(AOA_logDWS ~ irrigation * fertilization * sampling.date + 
                                    (1|plot), data = qPCR.BS, na.action = na.omit)
anova(aoa.log.dws.lmer)
shapiro.test(resid(aoa.log.dws.lmer)) #normal
library(DHARMa) 
plot(simulateResiduals(aoa.log.dws.lmer))

# test assumptions:
check_homogeneity(aoa.log.dws.aov) # slightly not homogen, but anova is robust in balanced data anyway
check_sphericity(aoa.log.dws.aov) # OK
aoa.log.dws.is_norm <- check_normality(aoa.log.dws.aov)
shapiro.test(qPCR.BS$AOA_logDWS) # NORMAL
plot(aoa.log.dws.is_norm, type = "qq")
plot(aoa.log.dws.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "AOA_logDWS", ggtheme = theme_bw())

# Pairwise comparison:
aoa.log.dws.emm <- emmeans(aoa.log.dws.lmer, ~ irrigation | fertilization*sampling.date)
aoa.log.dws.pair <- pairs(aoa.log.dws.emm)
aoa.log.dws.pair.DF <- as.data.frame(summary(aoa.log.dws.pair)) # nothing signif
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil')
write.csv(aoa.log.dws.pair.DF, file = "AOA_log_dws_pair.csv")

# 2. AOB

# anova test for transformed AOB DWS
aob.log.dws.aov <- aov_ez("plot", "AOB_logDWS", qPCR.BS, 
                          within = "sampling.date",
                          between = c("fertilization","irrigation"),
                          type = 2,
                          return = afex_options("return_aov"),
                          anova_table = list())
aob.log.dws.aov

# Three-Way Mixed (Split-Plot) ANOVA 
aob.log.dws.aov2 <- anova_test(data = qPCR.BS, 
                               dv = AOB_logDWS, 
                               wid = plot, 
                               type = 2,
                               within = sampling.date, 
                               between = c(irrigation, fertilization))
get_anova_table(aob.log.dws.aov2) # sphericity is automatically corrected in anova_test
aob.log.dws.aov2.df <- get_anova_table(aob.log.dws.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/AOB_log_dws")
table2csv(x=aob.log.dws.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(aob.log.dws.aov) # good
check_sphericity(aob.log.dws.aov) #Sphericity violated
aob.log.dws.is_norm <- check_normality(aob.log.dws.aov)
plot(aob.log.dws.is_norm, type = "qq")
plot(aob.log.dws.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "AOB_nbc_per_g_DW_soil", ggtheme = theme_bw())
ggqqplot(qPCR.BS, "AOB_logDWS", ggtheme = theme_bw())
hist(qPCR.BS$AOB_nbc_per_g_DW_soil)
hist(qPCR.BS$AOB_logDWS)
shapiro.test(qPCR.BS$AOB_nbc_per_g_DW_soil) # not NORMAL
shapiro.test(qPCR.BS$AOB_logDWS) # not NORMAL

# try another model
set.seed(13)
aob.log.dws.lmer <- lmerTest::lmer(AOB_logDWS ~ irrigation * fertilization * sampling.date + 
                                    (1|plot), data = qPCR.BS, na.action = na.omit)
anova(aob.log.dws.lmer)
shapiro.test(resid(aob.log.dws.lmer)) #normal
library(DHARMa) 
plot(simulateResiduals(aob.log.dws.lmer))

# Pairwise comparison:
aob.log.dws.emm <- emmeans(aob.log.dws.lmer, ~ irrigation | fertilization*sampling.date)
aob.log.dws.pair <- pairs(aob.log.dws.emm)
aob.log.dws.pair
aob.log.dws.pair.DF <- as.data.frame(summary(aob.log.dws.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(aob.log.dws.pair.DF, file = "AOB_log_dws_pair.csv")



# anova test for non transformed AOB DWS

aob.dws.aov <- aov_ez("plot", "AOB_nbc_per_g_DW_soil", qPCR.BS, 
                          within = "sampling.date",
                          between = c("fertilization","irrigation"),
                          type = 2,
                          return = afex_options("return_aov"),
                          anova_table = list(correction="none"))
aob.dws.aov

# Three-Way Mixed (Split-Plot) ANOVA 
aob.dws.aov2 <- anova_test(data = qPCR.BS, 
                               dv = AOB_nbc_per_g_DW_soil, 
                               wid = plot, 
                               type = 2,
                               within = sampling.date, 
                               between = c(irrigation, fertilization))
get_anova_table(aob.dws.aov2) # sphericity is automatically corrected in anova_test
aob.dws.aov2.df <- get_anova_table(aob.dws.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil/AOB_dws")
table2csv(x=aob.dws.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(aob.dws.aov) # good
check_sphericity(aob.dws.aov) # good
aob.dws.is_norm <- check_normality(aob.dws.aov)
plot(aob.dws.is_norm, type = "qq")
plot(aob.dws.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "AOB_nbc_per_g_DW_soil", ggtheme = theme_bw())
hist(qPCR.BS$AOB_nbc_per_g_DW_soil)

# try another model
set.seed(13)
aob.dws.lmer <- lmerTest::lmer(AOB_nbc_per_g_DW_soil ~ irrigation * fertilization * sampling.date + 
                                    (1|plot), data = qPCR.BS, na.action = na.omit)
anova(aob.dws.lmer)
shapiro.test(resid(aob.dws.lmer)) #normal
plot(simulateResiduals(aob.dws.lmer))

# tidy anova table
knitr::kable(nice(aob.dws.aov))
# Pairwise comparison:
aob.dws.emm <- emmeans(aob.dws.lmer, ~ irrigation | fertilization*sampling.date)
aob.dws.pair <- pairs(aob.dws.emm)
aob.dws.pair
aob.dws.pair.DF <- as.data.frame(summary(aob.dws.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil')
write.csv(aob.dws.pair.DF, file = "AOB_dws_pair.csv")

# 3. COMAMMOX A

# anova test for transformed Comammox A DWS
comA.log.dws.aov <- aov_ez("plot", "ComA_logDWS", qPCR.BS, 
                          within = "sampling.date",
                          between = c("fertilization","irrigation"),
                          type = 2,
                          return = afex_options("return_aov"),
                          anova_table = list())
comA.log.dws.aov

# Three-Way Mixed (Split-Plot) ANOVA 
comA.log.dws.aov2 <- anova_test(data = qPCR.BS, 
                               dv = ComA_logDWS, 
                               wid = plot, 
                               type = 2,
                               within = sampling.date, 
                               between = c(irrigation, fertilization))
get_anova_table(comA.log.dws.aov2)# sphericity is automatically corrected in anova_test
comA.log.dws.aov2.df <- get_anova_table(comA.log.dws.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/ComA_log_dws")
table2csv(x=comA.log.dws.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(comA.log.dws.aov) # good
check_sphericity(comA.log.dws.aov) #Sphericity violated
comA.log.dws.is_norm <- check_normality(comA.log.dws.aov)
plot(comA.log.dws.is_norm, type = "qq")
plot(comA.log.dws.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "ComA_nbc_per_g_DW_soil", ggtheme = theme_bw()) # violated
ggqqplot(qPCR.BS, "ComA_logDWS", ggtheme = theme_bw()) #better

# try another model
set.seed(13)
comA.log.dws.lmer <- lmerTest::lmer(ComA_logDWS ~ irrigation * fertilization * sampling.date + 
                                    (1|plot), data = qPCR.BS, na.action = na.omit)
anova(comA.log.dws.lmer)
shapiro.test(resid(comA.log.dws.lmer)) #normal
library(DHARMa) 
plot(simulateResiduals(comA.log.dws.lmer))

# tidy anova table
knitr::kable(nice(comA.log.dws.aov))
# Pairwise comparison:
comA.log.dws.emm <- emmeans(comA.log.dws.lmer, ~ irrigation | fertilization*sampling.date, lmer.df = "satterthwaite")
comA.log.dws.pair <- pairs(comA.log.dws.emm)
comA.log.dws.pair
comA.log.dws.pair.DF <- as.data.frame(summary(comA.log.dws.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(comA.log.dws.pair.DF, file = "ComA_log_dws_pair.csv")

# 4. COMAMMOX B

# anova test for transformed Comammox B DWS
comB.log.dws.aov <- aov_ez("plot", "ComB_logDWS", qPCR.BS, 
                           within = "sampling.date",
                           between = c("fertilization","irrigation"),
                           type = 2,
                           correction_aov = "GG",
                           return = afex_options("return_aov"),
                           anova_table = list())
comB.log.dws.aov

# Three-Way Mixed (Split-Plot) ANOVA 
comB.log.dws.aov2 <- anova_test(data = qPCR.BS, 
                                dv = ComB_logDWS, 
                                wid = plot, 
                                type = 2,
                                within = sampling.date, 
                                between = c(irrigation, fertilization))
get_anova_table(comB.log.dws.aov2)
comB.log.dws.aov2.df <- get_anova_table(comB.log.dws.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/ComB_log_dws")
table2csv(x=comB.log.dws.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(comB.log.dws.aov) # good
check_sphericity(comB.log.dws.aov) #Sphericity violated
comB.log.dws.is_norm <- check_normality(comB.log.dws.aov)
plot(comB.log.dws.is_norm, type = "qq")
plot(comB.log.dws.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "ComB_nbc_per_g_DW_soil", ggtheme = theme_bw())
ggqqplot(qPCR.BS, "ComB_logDWS", ggtheme = theme_bw()) 
hist(qPCR.BS$ComB_nbc_per_g_DW_soil)
hist(qPCR.BS$ComB_logDWS)
# Saphiro-Wilk for normality
qPCR.BS.ComB_logDWS.SW <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(ComB_logDWS)


# try another model
set.seed(13)
comB.log.dws.lmer <- lmerTest::lmer(ComB_logDWS ~ irrigation * fertilization * sampling.date + 
                                    (1|plot), data = qPCR.BS, na.action = na.omit)
anova(comB.log.dws.lmer)
shapiro.test(resid(comB.log.dws.lmer)) #normal
plot(simulateResiduals(comB.log.dws.lmer))

# tidy anova table
knitr::kable(nice(comB.log.dws.aov))
# Pairwise comparison:
comB.log.dws.emm <- emmeans(comB.log.dws.lmer, ~ irrigation | fertilization*sampling.date,lmer.df = "satterthwaite")
comB.log.dws.pair <- pairs(comB.log.dws.emm)
comB.log.dws.pair.DF <- as.data.frame(summary(comB.log.dws.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(comB.log.dws.pair.DF, file = "ComB_log_dws_pair.csv")

# 5. 16S

# anova test for transformed The Total microbial per DWS
tot.log.dws.aov <- aov_ez("plot", "Tot_logDWS", qPCR.BS, 
                           within = "sampling.date",
                           between = c("fertilization","irrigation"),
                           type = 2,
                           return = afex_options("return_aov"),
                           anova_table = list())
tot.log.dws.aov

# Three-Way Mixed (Split-Plot) ANOVA 
tot.log.dws.aov2 <- anova_test(data = qPCR.BS, 
                                dv = Tot_logDWS, 
                                wid = plot, 
                                type = 2,
                                within = c(sampling.date, fertilization),
                                between = irrigation)
get_anova_table(tot.log.dws.aov2)

pwc <- qPCR.BS %>%
  group_by(sampling.date, fertilization) %>%
  pairwise_t_test(
    Tot_logDWS ~ irrigation, paired = TRUE, 
    p.adjust.method = "fdr")
pwc
pwc <- pwc %>% add_xy_position(x = "sampling.date")

tot.log.dws.aov2.df <- get_anova_table(tot.log.dws.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/16S_log_dws")
table2csv(x=tot.log.dws.aov2.df,file=filen, digits = 1, digitspvals = 3)

# try another model
#tot.log.dws.lme <- lme(Tot_logDWS ~ irrigation*fertilization*sampling.date, random=~1|rep/fertilization/irrigation, data=qPCR.BS,na.action = na.omit, method = "REML")
#anova(tot.log.dws.lme)
#tot.log.dws.lmer <- lmerTest::lmer(Tot_logDWS ~ irrigation*fertilization*sampling.date+(1|plot), data=qPCR.BS,na.action = na.omit, REML=F)
#anova(tot.log.dws.lmer)
#tot.log.dws.lmer2 <- lmerTest::lmer(Tot_logDWS ~ irrigation*fertilization*sampling.date+(1|block), data=qPCR.BS,na.action = na.omit, REML=F)
#anova(tot.log.dws.lmer2)
tot.log.dws.lmer <- lmerTest::lmer(Tot_logDWS ~ irrigation*fertilization*sampling.date+
                                    (sampling.date|rep), data=qPCR.BS, na.action=na.omit)

tot.log.dws.lmer <- lmerTest::lmer(Tot_logDWS ~ irrigation*fertilization*sampling.date+
                                    (1|block/plot), data=qPCR.BS, na.action=na.omit)

#t=afex::mixed(Tot_logDWS ~ irrigation*fertilization*sampling.date+
                                    #(1|plot), data=qPCR.BS, na.action=na.omit)

anova(tot.log.dws.lmer)
shapiro.test(resid(tot.log.dws.lmer)) #normal
plot(simulateResiduals(tot.log.dws.lmer))

# test assumptions:
check_homogeneity(tot.log.dws.aov) # good
check_sphericity(tot.log.dws.aov) #Sphericity violated
tot.log.dws.is_norm <- check_normality(tot.log.dws.aov)
plot(tot.log.dws.is_norm, type = "qq")
plot(tot.log.dws.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "Tot_nbc_per_g_DW_soil", ggtheme = theme_bw())
ggqqplot(qPCR.BS, "Tot_logDWS", ggtheme = theme_bw()) 
hist(qPCR.BS$Tot_nbc_per_g_DW_soil)
hist(qPCR.BS$Tot_logDWS)
# Saphiro-Wilk for normality
qPCR.BS.Tot_logDWS.SW <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(Tot_logDWS)

# tidy anova table
knitr::kable(nice(tot.log.dws.aov))
# Pairwise comparison:
tot.log.dws.emm <- emmeans(t, ~ irrigation | fertilization*sampling.date,mode = "satterthwaite")
tot.log.dws.pair <- pairs(tot.log.dws.emm)
tot.log.dws.pair.DF <- as.data.frame(summary(tot.log.dws.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(tot.log.dws.pair.DF, file = "16S_log_dws_pair.csv")

# 6. AOA/16S RATIO

### anova test for AOA/16S Ratio Percentage

AOA_16.percent.ratio.aov <- aov_ez("plot", "AOA_16S_ratio_percent", qPCR.BS, 
                          within = "sampling.date",
                          between = c("fertilization","irrigation"),
                          type = 2,
                          correction_aov = "GG",
                          return = afex_options("return_aov"),
                          anova_table = list())
AOA_16.percent.ratio.aov 
# test assumptions:
check_homogeneity(AOA_16.percent.ratio.aov) # good
check_sphericity(AOA_16.percent.ratio.aov) #Sphericity violated
AOA_16.percent.ratio.is_norm <- check_normality(AOA_16.percent.ratio.aov)
plot(AOA_16.percent.ratio.is_norm, type = "qq")
plot(AOA_16.percent.ratio.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "AOA_16S_ratio_percent", ggtheme = theme_bw())
hist(qPCR.BS$AOA_16S_ratio_percent)
# Saphiro-Wilk for normality
#qPCR.BS.AOA_16S_ratio_arc.SW <- qPCR.BS %>%
  #group_by(irrigation, fertilization, sampling.date) %>%
shapiro_test(qPCR.BS$AOA_16S_ratio_percent)

# try another model
AOA_16.percent.ratio.lme <- lme(AOA_16S_ratio_percent~ irrigation*fertilization*sampling.date, random=~1|rep/sampling.date/fertilization/irrigation, data=qPCR.BS,method="REML",na.action = na.omit)
anova(AOA_16.percent.ratio.lme)

AOA_16.percent.ratio.lmer <- lmerTest::lmer(AOA_16S_ratio_percent ~ irrigation*fertilization*sampling.date+(1|plot), data=qPCR.BS,na.action = na.omit)
anova(AOA_16.percent.ratio.lmer)



# tidy anova table
knitr::kable(nice(AOA_16.percent.ratio.aov))
# Pairwise comparison:
AOA_16.percent.ratio.emm <- emmeans(AOA_16.percent.ratio.lmer, ~ irrigation | fertilization*sampling.date,model="multivariate")
AOA_16.percent.ratio.pair <- pairs(AOA_16.percent.ratio.emm)
AOA_16.percent.ratio.pair.DF <- as.data.frame(AOA_16.percent.ratio.pair)

### anova test for AOA/16S Ratio Arcsin

AOA_16.arc.ratio.aov <- aov_ez("plot", "AOA_16.arc.ratio", qPCR.BS, 
                                   within = "sampling.date",
                                   between = c("fertilization","irrigation"),
                                   type = 2,
                                   return = afex_options("return_aov"),
                                   anova_table = list(correction = "GG"))
AOA_16.arc.ratio.aov

# Three-Way Mixed (Split-Plot) ANOVA 
AOA_16.arc.ratio.aov2 <- anova_test(data = qPCR.BS, 
                               dv = AOA_16.arc.ratio, 
                               wid = plot, 
                               type = 2,
                               within = sampling.date, 
                               between = c(irrigation, fertilization))
get_anova_table(AOA_16.arc.ratio.aov2)
AOA_16.arc.ratio.aov2.df <- get_anova_table(AOA_16.arc.ratio.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/AOA_16_arc_ratio")
table2csv(x=AOA_16.arc.ratio.aov2.df,file=filen, digits = 1, digitspvals = 3)


# test assumptions:
check_homogeneity(AOA_16.arc.ratio.aov) # good
check_sphericity(AOA_16.arc.ratio.aov) #Sphericity violated
AOA_16.arc.ratio.is_norm <- check_normality(AOA_16.arc.ratio.aov)
ggqqplot(qPCR.BS, "AOA_16.arc.ratio", ggtheme = theme_bw())
plot(AOA_16.arc.ratio.is_norm, type = "qq")
plot(AOA_16.arc.ratio.is_norm, type = "qq", detrend = TRUE)
# NORMAL
# tidy anova table
knitr::kable(nice(AOA_16.arc.ratio.aov))
# Pairwise comparison:
AOA_16.arc.ratio.emm <- emmeans(AOA_16.arc.ratio.aov, ~ irrigation | fertilization*sampling.date)
AOA_16.arc.ratio.pair <- pairs(AOA_16.arc.ratio.emm)
AOA_16.arc.ratio.pair.DF <- as.data.frame(summary(AOA_16.arc.ratio.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(AOA_16.arc.ratio.pair.DF, file = "AOA_16_arc_ratio_pair.csv")

# 7. AOB/16S RATIO

### anova test for AOB/16S Ratio Percentage

AOB_16.percent.ratio.aov <- aov_ez("plot", "AOB_16S_ratio_percent", qPCR.BS, 
                                   within = "sampling.date",
                                   between = c("fertilization","irrigation"),
                                   type = 2,
                                   return = afex_options("return_aov"),
                                   anova_table = list(correction = "GG"))
AOB_16.percent.ratio.aov

#aov_car(AOB_16S_ratio_percent ~ irrigation * fertilization + Error(plot/(sampling.date)),
        #data = qPCR.BS)

# Three-Way Mixed (Split-Plot) ANOVA 
AOB_16.percent.ratio.aov2 <- anova_test(data = qPCR.BS.noNA, 
                                    dv = AOB_16S_ratio_percent, 
                                    wid = plot, 
                                    type = 2,
                                    within = sampling.date, 
                                    between = c(irrigation, fertilization))
get_anova_table(AOB_16.percent.ratio.aov2)
AOB_16.percent.ratio.aov2.df <- get_anova_table(AOB_16.percent.ratio.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil/AOB_16_percent_ratio")
table2csv(x=AOB_16.percent.ratio.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(AOB_16.percent.ratio.aov) # good
check_sphericity(AOB_16.percent.ratio.aov) #Sphericity violated
AOB_16.percent.ratio.is_norm <- check_normality(AOB_16.percent.ratio.aov)
plot(AOB_16.percent.ratio.is_norm, type = "qq")
plot(AOB_16.percent.ratio.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "AOB_16S_ratio_percent", ggtheme = theme_bw())
ggqqplot(qPCR.BS, "AOB_16.arc.ratio", ggtheme = theme_bw())
hist(qPCR.BS$AOB_16S_ratio_percent)
hist(qPCR.BS$AOB_16.arc.ratio)
shapiro_test(qPCR.BS$AOB_16.arc.ratio)

# NORMAL
# tidy anova table
knitr::kable(nice(AOB_16.percent.ratio.aov))
# Pairwise comparison:
AOB_16.percent.ratio.emm <- emmeans(AOB_16.percent.ratio.aov, ~ irrigation | fertilization*sampling.date)
AOB_16.percent.ratio.pair <- pairs(AOB_16.percent.ratio.emm)
AOB_16.percent.ratio.pair.DF <- as.data.frame(summary(AOB_16.percent.ratio.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil')
write.csv(AOB_16.percent.ratio.pair.DF, file = "AOB_16_percent_ratio_pair.csv")


### anova test for AOB/16S Ratio Arcsin

AOB_16.arc.ratio.aov <- aov_ez("plot", "AOB_16.arc.ratio", qPCR.BS, 
                               within = "sampling.date",
                               between = c("fertilization","irrigation"),
                               type = 2,
                               return = afex_options("return_aov"),
                               na.rm=T,
                               anova_table = list())
AOB_16.arc.ratio.aov



qPCR.BS$ferblo <- factor(qPCR.BS$fertilization:qPCR.BS$block)
qPCR.BS$irriblo <- factor(qPCR.BS$irrigation:qPCR.BS$block)
#Three-Way Mixed (Split-Plot) ANOVA 
AOB_16.arc.ratio.aov2 <- anova_test(data = qPCR.BS, 
                                    dv = AOB_16.arc.ratio, 
                                    wid = plot, 
                                    type = 3,
                                    within = sampling.date, 
                                    between = c(irrigation, fertilization))
get_anova_table(AOB_16.arc.ratio.aov2)
AOB_16.arc.ratio.aov2.df <- get_anova_table(AOB_16.arc.ratio.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/AOB_16_arc_ratio")
table2csv(x=AOB_16.arc.ratio.aov2.df,file=filen, digits = 1, digitspvals = 3)

# try another model
me=lme(AOB_16S_ratio_percent ~ irrigation*fertilization*sampling.date, random=~1|rep/sampling.date/fertilization/irrigation, data=qPCR.BS,na.action = na.omit)
anova(me)


#_________________________________________________________________________________________________________________________________________________________________
AOB_16.arc.ratio.aov3 <- lmerTest::lmer(qPCR.BS$AOB_16.arc.ratio ~ irrigation*fertilization*sampling.date+(1|plot), data=qPCR.BS,
                                        na.action = na.omit)
anova(AOB_16.arc.ratio.aov3)
shapiro.test(resid(AOB_16.arc.ratio.aov3)) #ok
# Lavene test
AOB_16.arc.ratio.Lave <- qPCR.BS %>%
  group_by(sampling.date,fertilization) %>%
  levene_test(AOB_16.arc.ratio ~ irrigation)


t1=lmerTest::lmer(AOB_16S_ratio_percent ~ irrigation*fertilization*sampling.date+(1|rep:sampling.date), 
                  data = qPCR.BS,na.action = na.omit)
anova(t1)
AIC(t1)
library(lmerTest)
t <- lmerTest::lmer(AOB_16S_ratio_percent ~ irrigation*fertilization*sampling.date+(1|plot), data = qPCR.BS, 
                    na.action = na.omit, REML=F)
anova(t)
AIC(t)

t2=lmerTest::lmer(AOB_16S_ratio_percent ~ irrigation*fertilization*sampling.date+(1|block:sampling.date), data = qPCR.BS, 
                    na.action = na.omit)
#______________________________________________________________________________________________________________________







# test assumptions:
check_homogeneity(AOB_16.arc.ratio.aov) # good
check_sphericity(AOB_16.arc.ratio.aov) #Sphericity violated
AOB_16.arc.ratio.is_norm <- check_normality(AOB_16.arc.ratio.aov)
plot(AOB_16.arc.ratio.is_norm, type = "qq")
plot(AOB_16.arc.ratio.is_norm, type = "qq", detrend = TRUE)
# NORMAL
# tidy anova table
knitr::kable(nice(AOB_16.arc.ratio.aov))
# Pairwise comparison:
AOB_16.arc.ratio.emm <- emmeans(me, ~ irrigation | fertilization*sampling.date)
AOB_16.arc.ratio.pair <- pairs(AOB_16.arc.ratio.emm)
AOB_16.arc.ratio.pair.DF <- as.data.frame(summary(AOB_16.arc.ratio.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(AOB_16.arc.ratio.pair.DF, file = "AOB_16_arc_ratio_pair.csv")

# 8. COMAMMOX A/16S RATIO

### anova test for Comammox A/16S Ratio Percentage

ComA_16.percent.ratio.aov <- aov_ez("plot", "ComA_16S_ratio_percent", qPCR.BS, 
                                   within = "sampling.date",
                                   between = c("fertilization","irrigation"),
                                   type = 2,
                                   #correction_aov = "GG",
                                   return = afex_options("return_aov"),
                                   anova_table = list())
ComA_16.percent.ratio.aov

# test assumptions:
check_homogeneity(ComA_16.percent.ratio.aov) # slightly not homogen
check_sphericity(ComA_16.percent.ratio.aov) # good
ComA_16.percent.ratio.is_norm <- check_normality(ComA_16.percent.ratio.aov)
plot(ComA_16.percent.ratio.is_norm, type = "qq")
plot(ComA_16.percent.ratio.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "ComA_16S_ratio_percent", ggtheme = theme_bw())
ggqqplot(qPCR.BS, "ComA_16.arc.ratio", ggtheme = theme_bw())

# NORMAL
# tidy anova table
knitr::kable(nice(ComA_16.percent.ratio.aov))
# Pairwise comparison:
ComA_16.percent.ratio.emm <- emmeans(ComA_16.percent.ratio.aov, ~ irrigation | fertilization*sampling.date)
ComA_16.percent.ratio.pair <- pairs(ComA_16.percent.ratio.emm)
ComA_16.percent.ratio.pair.DF <- as.data.frame(ComA_16.percent.ratio.pair)

### anova test for ComA/16S Ratio Arcsin

ComA_16.arc.ratio.aov <- aov_ez("plot", "ComA_16.arc.ratio", qPCR.BS, 
                               within = "sampling.date",
                               between = c("fertilization","irrigation"),
                               type = 2,
                               #correction_aov = "GG",
                               return = afex_options("return_aov"),
                               anova_table = list())
ComA_16.arc.ratio.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComA_16.arc.ratio.aov2 <- anova_test(data = qPCR.BS, 
                                    dv = ComA_16.arc.ratio, 
                                    wid = plot, 
                                    type = 2,
                                    within = sampling.date, 
                                    between = c(irrigation, fertilization))
get_anova_table(ComA_16.arc.ratio.aov2)
ComA_16.arc.ratio.aov2.df <- get_anova_table(ComA_16.arc.ratio.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/ComA_16_arc_ratio")
table2csv(x=ComA_16.arc.ratio.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComA_16.arc.ratio.aov) # good
check_sphericity(ComA_16.arc.ratio.aov) # good
ComA_16.arc.ratio.is_norm <- check_normality(ComA_16.arc.ratio.aov)
plot(ComA_16.arc.ratio.is_norm, type = "qq")
plot(ComA_16.arc.ratio.is_norm, type = "qq", detrend = TRUE)
# NORMAL
# tidy anova table
knitr::kable(nice(ComA_16.arc.ratio.aov))
# Pairwise comparison:
ComA_16.arc.ratio.emm <- emmeans(ComA_16.arc.ratio.aov, ~ irrigation | fertilization*sampling.date)
ComA_16.arc.ratio.pair <- pairs(ComA_16.arc.ratio.emm)
ComA_16.arc.ratio.pair.DF <- as.data.frame(summary(ComA_16.arc.ratio.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(ComA_16.arc.ratio.pair.DF, file = "ComA_16_arc_ratio_pair.csv")

# 9. COMAMMOX B /16S RATIO

### anova test for Comammox B /16S Ratio Percentage

ComB_16.percent.ratio.aov <- aov_ez("plot", "ComB_16S_ratio_percent", qPCR.BS, 
                                    within = "sampling.date",
                                    between = c("fertilization","irrigation"),
                                    type = 2,
                                    correction_aov = "GG",
                                    return = afex_options("return_aov"),
                                    anova_table = list())
ComB_16.percent.ratio.aov

# test assumptions:
check_homogeneity(ComB_16.percent.ratio.aov) # homogen
check_sphericity(ComB_16.percent.ratio.aov) # violated
ComB_16.percent.ratio.is_norm <- check_normality(ComB_16.percent.ratio.aov)
plot(ComB_16.percent.ratio.is_norm, type = "qq")
plot(ComB_16.percent.ratio.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "ComB_16S_ratio_percent", ggtheme = theme_bw())
ggqqplot(qPCR.BS, "ComB_16.arc.ratio", ggtheme = theme_bw())
hist(qPCR.BS$ComB_16S_ratio_percent)

# NORMAL
# tidy anova table
knitr::kable(nice(ComB_16.percent.ratio.aov))
# Pairwise comparison:
ComB_16.percent.ratio.emm <- emmeans(ComB_16.percent.ratio.aov, ~ irrigation | fertilization*sampling.date)
ComB_16.percent.ratio.pair <- pairs(ComB_16.percent.ratio.emm)
ComB_16.percent.ratio.pair.DF <- as.data.frame(ComB_16.percent.ratio.pair)

### anova test for ComB/16S Ratio Arcsin

ComB_16.arc.ratio.aov <- aov_ez("plot", "ComB_16.arc.ratio", qPCR.BS, 
                                within = "sampling.date",
                                between = c("fertilization","irrigation"),
                                type = 2,
                                correction_aov = "GG",
                                return = afex_options("return_aov"),
                                anova_table = list())
ComB_16.arc.ratio.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComB_16.arc.ratio.aov2 <- anova_test(data = qPCR.BS, 
                                     dv = ComB_16.arc.ratio, 
                                     wid = plot, 
                                     type = 2,
                                     within = sampling.date, 
                                     between = c(irrigation, fertilization))
get_anova_table(ComB_16.arc.ratio.aov2)
ComB_16.arc.ratio.aov2.df <- get_anova_table(ComB_16.arc.ratio.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/ComB_16_arc_ratio")
table2csv(x=ComB_16.arc.ratio.aov2.df,file=filen, digits = 1, digitspvals = 3)


# test assumptions:
check_homogeneity(ComB_16.arc.ratio.aov) # good
check_sphericity(ComB_16.arc.ratio.aov) # violated
ComB_16.arc.ratio.is_norm <- check_normality(ComB_16.arc.ratio.aov)
plot(ComB_16.arc.ratio.is_norm, type = "qq")
plot(ComB_16.arc.ratio.is_norm, type = "qq", detrend = TRUE)
# NORMAL
# tidy anova table
knitr::kable(nice(ComB_16.arc.ratio.aov))
# Pairwise comparison:
ComB_16.arc.ratio.emm <- emmeans(ComB_16.arc.ratio.aov, ~ irrigation | fertilization*sampling.date)
ComB_16.arc.ratio.pair <- pairs(ComB_16.arc.ratio.emm)
ComB_16.arc.ratio.pair.DF <- as.data.frame(summary(ComB_16.arc.ratio.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(ComB_16.arc.ratio.pair.DF, file = "ComB_16_arc_ratio_pair.csv")


# 10. AOA/AOB RATIO

### anova test for AOA/AOB Ratio without transformation

AOA_AOB.ratio.aov <- aov_ez("plot", "AOA_AOB_ratio", qPCR.BS, 
                                   within = "sampling.date",
                                   between = c("fertilization","irrigation"),
                                   type = 2,
                                   correction_aov = "GG",
                                   return = afex_options("return_aov"),
                                   anova_table = list())
AOA_AOB.ratio.aov

# test assumptions:
check_homogeneity(AOA_AOB.ratio.aov) # good
check_sphericity(AOA_AOB.ratio.aov) #Sphericity violated
AOA_AOB.ratio.is_norm <- check_normality(AOA_AOB.ratio.aov)
plot(AOA_AOB.ratio.is_norm, type = "qq")
plot(AOA_AOB.ratio.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "AOA_AOB_ratio", ggtheme = theme_bw())
shapiro_test(qPCR.BS$AOA_AOB_ratio)
hist(qPCR.BS$AOA_AOB_ratio)

#facet_grid(sampling.date ~ irrigation*fertilization, labeller = "label_both")
# NORMAL
# tidy anova table
knitr::kable(nice(AOA_AOB.ratio.aov))
# Pairwise comparison:
AOA_AOB.ratio.emm <- emmeans(AOA_AOB.ratio.aov, ~ irrigation | fertilization*sampling.date)
AOA_AOB.ratio.pair <- pairs(AOA_AOB.ratio.emm)
AOA_AOB.ratio.pair.DF <- as.data.frame(AOA_AOB.ratio.pair)

### anova test for AOA/AOB Ratio arcsin sqrt transformed

AOA_AOB.arc.ratio.aov <- aov_ez("plot", "AOA_AOB.arc.ratio", qPCR.BS, 
                            within = "sampling.date",
                            between = c("fertilization","irrigation"),
                            type = 2,
                            correction_aov = "GG",
                            return = afex_options("return_aov"),
                            anova_table = list())
AOA_AOB.arc.ratio.aov

# Three-Way Mixed (Split-Plot) ANOVA 
AOA_AOB.arc.ratio.aov2 <- anova_test(data = qPCR.BS, 
                                    dv = AOA_AOB.arc.ratio, 
                                    wid = plot, 
                                    type = 2,
                                    within = sampling.date, 
                                    between = c(irrigation, fertilization))
get_anova_table(AOA_AOB.arc.ratio.aov2)
AOA_AOB.arc.ratio.aov2.df <- get_anova_table(AOA_AOB.arc.ratio.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/AOA_AOB_arc_ratio")
table2csv(x=AOA_AOB.arc.ratio.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(AOA_AOB.arc.ratio.aov) # good
check_sphericity(AOA_AOB.arc.ratio.aov) #Sphericity violated
AOA_AOB.arc.ratio.is_norm <- check_normality(AOA_AOB.arc.ratio.aov)
plot(AOA_AOB.arc.ratio.is_norm, type = "qq")
plot(AOA_AOB.arc.ratio.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.BS, "AOA_AOB.arc.ratio", ggtheme = theme_bw())
#facet_grid(sampling.date ~ irrigation*fertilization, labeller = "label_both")

# tidy anova table
knitr::kable(nice(AOA_AOB.arc.ratio.aov))
# Pairwise comparison:
AOA_AOB.arc.ratio.emm <- emmeans(AOA_AOB.arc.ratio.aov, ~ irrigation | fertilization*sampling.date)
AOA_AOB.arc.ratio.pair <- pairs(AOA_AOB.arc.ratio.emm)
AOA_AOB.arc.ratio.pair.DF <- as.data.frame(summary(AOA_AOB.arc.ratio.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(AOA_AOB.arc.ratio.pair.DF, file = "AOA_AOB_arc_ratio_pair.csv")


# 11. ComA/ComB RATIO

### anova test for ComA/ComB Ratio without transformation

ComA_ComB.ratio.aov <- aov_ez("plot", "ComA_ComB_ratio", qPCR.BS, 
                            within = "sampling.date",
                            between = c("fertilization","irrigation"),
                            type = 3, # there is interaction effect
                            return = afex_options("return_aov"),
                            anova_table = list(correction = "none"))
ComA_ComB.ratio.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComA_ComB.ratio.aov2 <- anova_test(data = qPCR.BS, 
                                       dv = ComA_ComB_ratio, 
                                       wid = plot, 
                                       type = 3, # there is interaction
                                       within = sampling.date, 
                                       between = c(irrigation, fertilization))
get_anova_table(ComA_ComB.ratio.aov2)
ComA_ComB.ratio.aov2.df <- get_anova_table(ComA_ComB.ratio.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil/ComA_ComB_ratio")
table2csv(x=ComA_ComB.ratio.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComA_ComB.ratio.aov) # good
check_sphericity(ComA_ComB.ratio.aov) #good
ComA_ComB.ratio.is_norm <- check_normality(ComA_ComB.ratio.aov)
ggqqplot(qPCR.BS, "ComA_ComB_ratio", ggtheme = theme_bw())
  #facet_grid(sampling.date ~ irrigation*fertilization, labeller = "label_both")
plot(ComA_ComB.ratio.is_norm, type = "qq")
plot(ComA_ComB.ratio.is_norm, type = "qq", detrend = TRUE)
shapiro_test(qPCR.BS$ComA_ComB_ratio) #good

# tidy anova table
knitr::kable(nice(ComA_ComB.ratio.aov))
t=nice(ComA_ComB.ratio.aov)
# Pairwise comparison:
ComA_ComB.ratio.emm <- emmeans(ComA_ComB.ratio.aov, ~ irrigation | fertilization*sampling.date)
ComA_ComB.ratio.pair <- pairs(ComA_ComB.ratio.emm)
ComA_ComB.ratio.pair.DF <- as.data.frame(summary(ComA_ComB.ratio.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Bulk Soil')
write.csv(ComA_ComB.ratio.pair.DF, file = "ComA_ComB_ratio_pair.csv")


### anova test for ComA/ComB Ratio arcsin transformed

ComA_ComB.arc.ratio.aov <- aov_ez("plot", "ComA_ComB.arc.ratio", qPCR.BS, 
                              within = "sampling.date",
                              between = c("fertilization","irrigation"),
                              type = 3,
                              #correction_aov = "GG",
                              return = afex_options("return_aov"),
                              anova_table = list())
ComA_ComB.arc.ratio.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComA_ComB.arc.ratio.aov2 <- anova_test(data = qPCR.BS, 
                                     dv = ComA_ComB.arc.ratio, 
                                     wid = plot, 
                                     type = 3, # there is interaction
                                     within = sampling.date, 
                                     between = c(irrigation, fertilization))
get_anova_table(ComA_ComB.arc.ratio.aov2)
ComA_ComB.arc.ratio.aov2.df <- get_anova_table(ComA_ComB.arc.ratio.aov2)
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/ComA_ComB_arc_ratio")
table2csv(x=ComA_ComB.arc.ratio.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComA_ComB.arc.ratio.aov) # good
check_sphericity(ComA_ComB.arc.ratio.aov) #Sphericity violated
ComA_ComB.arc.ratio.is_norm <- check_normality(ComA_ComB.arc.ratio.aov)
ggqqplot(qPCR.BS, "ComA_ComB.arc.ratio", ggtheme = theme_bw())
#facet_grid(sampling.date ~ irrigation*fertilization, labeller = "label_both")
plot(ComA_ComB.arc.ratio.is_norm, type = "qq")
plot(ComA_ComB.arc.ratio.is_norm, type = "qq", detrend = TRUE)
# NORMAL
# tidy anova table
knitr::kable(nice(ComA_ComB.arc.ratio.aov))
# Pairwise comparison:
ComA_ComB.arc.ratio.emm <- emmeans(ComA_ComB.arc.ratio.aov, ~ irrigation | fertilization*sampling.date)
ComA_ComB.arc.ratio.pair <- pairs(ComA_ComB.arc.ratio.emm)
ComA_ComB.arc.ratio.pair.DF <- as.data.frame(summary(ComA_ComB.arc.ratio.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat')
write.csv(ComA_ComB.arc.ratio.pair.DF, file = "ComA_ComB_arc_ratio_pair.csv")



##############################################################################################################################################
aoa.BS.copies <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  get_summary_stats(AOA_nbc_per_g_DW_soil, type = "mean_sd")


aoa.BS.copies.plot <- ggboxplot(
  qPCR.BS, x = "irrigation", y = "AOA_nbc_per_g_DW_soil",
  color = "fertilization", palette = "jco",
  facet.by =  "sampling.date")
aoa.BS.copies.plot

# check assumption (outliers)
aoa.BS.copies.out <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  identify_outliers(AOA_nbc_per_g_DW_soil) # no extreme outliers
# Saphiro-Wilk for normality
aoa.BS.copies.SW <- qPCR.BS %>%
  group_by(irrigation, fertilization, sampling.date) %>%
  shapiro_test(AOA_nbc_per_g_DW_soil)
ggqqplot(qPCR.BS, "AOA_nbc_per_g_DW_soil", ggtheme = theme_bw()) +
  facet_grid(sampling.date ~ fertilization, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.BS.copies.Lave <- qPCR.BS %>%
  group_by(sampling.date) %>%
  levene_test(AOA_nbc_per_g_DW_soil ~ irrigation*fertilization)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aoa.BS.copies.aov <- anova_test(
  data = qPCR.BS, dv = AOA_nbc_per_g_DW_soil, wid = plot,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aoa.BS.copies.aov)
# TWO-WAY
res.aov2 <- aov(AOA_nbc_per_g_DW_soil ~ irrigation + fertilization, data = qPCR.BS)
summary(res.aov2)
# THREE WAY
res.aov3 <- qPCR.BS %>% anova_test(AOA_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date)
res.aov3

res.aov3.ngDNA <- qPCR.BS %>% anova_test(AOA_nbc_per_ngDNA ~ irrigation*fertilization*sampling.date)
res.aov3.ngDNA

aoa.BS.ngDNA.aov <- anova_test(
  data = qPCR.BS, dv = AOA_nbc_per_g_DW_soil, wid = plot,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(aoa.BS.ngDNA.aov)
# model
aoa.BS.copies.mod <- lm(qPCR.BS$AOA_nbc_per_g_DW_soil ~ irrigation*fertilization*sampling.date, data=qPCR.BS)
anova(aoa.BS.copies.mod)
# pairwise comparisons
aoa.test.pw <- qPCR.BS %>%
  group_by(fertilization, sampling.date) %>%
  emmeans_test(AOA_nbc_per_g_DW_soil ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.BS.copies.mod)
aoa.test.pw
# Pairwise comparisons
library(emmeans)
pwc.aoa.DWS <- qPCR.BS %>%
  group_by(fertilization, sampling.date) %>%
  emmeans_test(AOA_nbc_per_g_DW_soil ~ irrigation, p.adjust.method = "BH") %>%
  select(-df, -statistic, -p) # Remove details
pwc.aoa.DWS # either grouped by fertilization or both by fertilization and per sampling date, drought is not significant

##subset
aoa.cop.M <- qPCR.BS[which(qPCR.BS$fertilization == "M"),]
aoa.cop.M
aoa.cop.D <- qPCR.BS[which(qPCR.BS$fertilization == "D"),]
aoa.cop.D
aoa.cop.K <- qPCR.BS[which(qPCR.BS$fertilization == "K"),]
aoa.cop.K
aoa.cop.M.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M)
summary(aoa.cop.M.aov)
TukeyHSD(aoa.cop.M.aov, conf.level=.95)

aoa.cop.D.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.D)
summary(aoa.cop.D.aov)
TukeyHSD(aoa.cop.D.aov, conf.level=.95)

aoa.cop.K.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.K)
summary(aoa.cop.K.aov)
TukeyHSD(aoa.cop.K.aov, conf.level=.95)

aoa.cop.M1 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "04-28-22" ),]
aoa.cop.M1.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M1)
summary(aoa.cop.M1.aov)
TukeyHSD(aoa.cop.M1.aov, conf.level=.95)

aoa.cop.M2 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "06-01-22" ),]
aoa.cop.M2.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M2)
summary(aoa.cop.M2.aov)
TukeyHSD(aoa.cop.M2.aov, conf.level=.95)

aoa.cop.M3 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "07-05-22" ),]
aoa.cop.M3.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M3)
summary(aoa.cop.M3.aov)
TukeyHSD(aoa.cop.M3.aov, conf.level=.95)

aoa.cop.M4 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "07-20-22" ),]
aoa.cop.M4.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M4)
summary(aoa.cop.M4.aov)
TukeyHSD(aoa.cop.M4.aov, conf.level=.95)

aoa.cop.M5 <- qPCR.BS[which(qPCR.BS$fertilization == "D"
                            & qPCR.BS$sampling.date == "09-13-22" ),]
aoa.cop.M5.aov <- aov(AOA_nbc_per_g_DW_soil ~ irrigation, data = aoa.cop.M5)
summary(aoa.cop.M5.aov)
TukeyHSD(aoa.cop.M5.aov, conf.level=.95)




########################################################################################################################################
# Plotting

qPCR.BS$x
qPCR.BS.ed <- qPCR.BS %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")

# 1. Bulk Soil - AOA

aoa.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  labs(y="AOA copy number per g dry soil", fill='Treatment')+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
aoa.cop.pwc.plot

aoa.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,7e+08)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('AOA'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('AOA')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("A. Bulk Soil")
aoa.cop.pwc.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_copies.Bulk.tiff",
       aoa.cop.pwc.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


# 2. Bulk Soil - AOB

aob.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOB_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,7e+08)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('AOB'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('AOB')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        #strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")#+ggtitle("A. Bulk Soil")
aob.cop.pwc.plot


# 3. Bulk Soil - Comammox A

comA.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=ComA_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,7e+08)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox A'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('Comammox A')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        #strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")#+ggtitle("A. Bulk Soil")
comA.cop.pwc.plot

# 4. Bulk Soil - Comammox B

comB.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=ComB_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,7e+08)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('Comammox B')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        #strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")#+ggtitle("A. Bulk Soil")
comB.cop.pwc.plot

# 4. Bulk Soil - 16S

tot.cop.pwc.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=Tot_nbc_per_g_DW_soil)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,3.2e+10)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  ylab(bquote('16S rRNA gene'~(copies~g^-1~dry~soil)))+
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

setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("16Scopies_BS.tiff",
       tot.cop.pwc.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)



library(patchwork)
copiesdws.All <-  aoa.cop.pwc.plot / aob.cop.pwc.plot / comA.cop.pwc.plot / comB.cop.pwc.plot
copiesdws.All

# 5. Bulk Soil - AOA/16S RATIO

AOA_16S.rat.plot <- ggplot(qPCR.BS.ed, aes(x=sampling.date, y=AOA_16S_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  ylim(0,8.9)+
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

setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("AOA_16S_ratio_BS.tiff",
       AOA_16S.rat.plot, device = "tiff",
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
        panel.grid.minor = element_blank())
AOB_16S.rat.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
aob_16S_percent_rat_emm <- qPCR.BS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(AOB_16S_ratio_percent ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t)
aob_16S_percent_rat_emm
aob_16S_percent.xy <- aob_16S_percent_rat_emm %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
aob_16S_percent.xy
# plotting the pairwise comparisons among treatments (emmeans results)
AOB_16S.rat.plot2 <- AOB_16S.rat.plot + 
  stat_pvalue_manual(aob_16S_percent.xy,label = "p.adj.signif", size=4, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 0, color = "blue",
                     tip.length = 0.01, hide.ns = TRUE)+
 scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
AOB_16S.rat.plot2










AOB_16.arc.ratio.pair.DF




























#__________________________________________________________________________________
# line chart AOA Copies per gram DWS Bulk Soil 

qPCR.BS.ed

# tidy up the data frame
qPCR.BS.sum <- qPCR.BS.ed %>%
  group_by(irrigation, fertilization,x, sampling.date,var3) %>%
  summarize(AOA_nbc_per_g_DW_soil=mean(AOA_nbc_per_g_DW_soil),
            sd.aoa.DWS=sd(AOA_nbc_per_g_DW_soil, na.rm = TRUE),
            AOA_nbc_per_ngDNA=mean(AOA_nbc_per_ngDNA),
            sd.aoa.ngDNA=sd(AOA_nbc_per_ngDNA, na.rm = TRUE),
            AOB_nbc_per_g_DW_soil=mean(AOB_nbc_per_g_DW_soil),
            sd.aob.DWS=sd(AOB_nbc_per_g_DW_soil, na.rm = TRUE),
            AOB_nbc_per_ngDNA=mean(AOB_nbc_per_ngDNA),
            sd.aob.ngDNA=sd(AOB_nbc_per_ngDNA, na.rm = TRUE))
qPCR.BS.sum
str(qPCR.BS.sum)

qPCR.BS.ed$lm_pred_val <- predict(aoa.BS.copies.mod,newdata = qPCR.BS.ed,
                              interval =  "confidence"
) %>% as.data.frame()


ggplot(qPCR.BS.ed, aes(x = sampling.date, y = AOA_nbc_per_g_DW_soil, linetype=irrigation)) +
  geom_jitter(position = position_jitter(0.2),
               color = "darkgray") + 
  geom_line(linewidth=1.15,aes(group = x,col=fertilization), data = qPCR.BS.sum) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"))+
  geom_errorbar(aes(ymin = AOA_nbc_per_g_DW_soil-sd.aoa.DWS, ymax = AOA_nbc_per_g_DW_soil+sd.aoa.DWS),
                data = qPCR.BS.sum, width = 0.2) +
  geom_point(data = qPCR.BS.sum, size = 2)+
  #geom_ribbon(aes(ymin = AOA_nbc_per_g_DW_soil+lm_pred_val$lwr, ymax = AOA_nbc_per_g_DW_soil-lm_pred_val$upr), alpha = 0.1)+
  ylab(bquote(~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  theme_bw() +
  theme(legend.title = element_blank(),
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

setwd('D:/Fina/INRAE_Project/microservices/')
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
qPCR.RS$sampling.date <- factor(qPCR.RS$sampling.date, levels = c("28/04/2022", "01/06/2022", "05/07/2022"),
                                labels = c("04-28-22", "06-01-22", "07-05-22"))

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


##### 5. 16S in Rhizosphere

# Anova test for non transformed 16S nbc per ngDNA 

Tot.ngDNA.aov <- aov_ez("plot", "Tot_nbc_per_ngDNA", qPCR.RS, 
                         within = "sampling.date",
                         between = c("irrigation","fertilization"),
                         #correction_aov = "GG",
                         type = 2,
                         return = afex_options("return_aov"),
                         anova_table = list(correction="none"))
Tot.ngDNA.aov

# Three-Way Mixed (Split-Plot) ANOVA 
Tot.ngDNA.aov2 <- anova_test(
  data = qPCR.RS, dv = Tot_nbc_per_ngDNA, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(Tot.ngDNA.aov2)
Tot.ngDNA.aov2.df <- get_anova_table(Tot.ngDNA.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/16S_ngDNA")
table2csv(x=Tot.ngDNA.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(Tot.ngDNA.aov) # good
check_sphericity(Tot.ngDNA.aov) # good
Tot.ngDNA.is_norm <- check_normality(Tot.ngDNA.aov)
shapiro.test(qPCR.RS$Tot_nbc_per_ngDNA) # NORMAL
plot(Tot.ngDNA.is_norm, type = "qq")
plot(Tot.ngDNA.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "Tot_nbc_per_ngDNA", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(Tot.ngDNA.aov))
# Pairwise comparison:
Tot.ngDNA.emm <- emmeans(Tot.ngDNA.aov, ~ irrigation | fertilization*sampling.date)
Tot.ngDNA.pair <- pairs(Tot.ngDNA.emm)
Tot.ngDNA.pair.DF <- as.data.frame(summary(Tot.ngDNA.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(Tot.ngDNA.pair.DF, file = "16S_ngDNA_pair.csv")


##### 6. AOA/16S Ratio in Rhizosphere

# Anova test for arcsin square root transformed AOA/16S Ratio

AOA_16S.arc.rat.rh.aov <- aov_ez("plot", "AOA_16.arc.ratio.rh", qPCR.RS, 
                        within = "sampling.date",
                        between = c("irrigation","fertilization"),
                        #correction_aov = "GG",
                        type = 2,
                        return = afex_options("return_aov"),
                        anova_table = list())
AOA_16S.arc.rat.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
AOA_16S.arc.rat.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = AOA_16.arc.ratio.rh, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(AOA_16S.arc.rat.rh.aov2)
AOA_16S.arc.rat.rh.aov2.df <- get_anova_table(AOA_16S.arc.rat.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOA_16S_arc_ratio")
table2csv(x=AOA_16S.arc.rat.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(AOA_16S.arc.rat.rh.aov) # good
check_sphericity(AOA_16S.arc.rat.rh.aov) # good
AOA_16S.arc.rat.rh.is_norm <- check_normality(AOA_16S.arc.rat.rh.aov)
shapiro.test(qPCR.RS$AOA_16.arc.ratio.rh) # NORMAL
plot(AOA_16S.arc.rat.rh.is_norm, type = "qq")
plot(AOA_16S.arc.rat.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOA_16.arc.ratio.rh", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(AOA_16S.arc.rat.rh.aov))
# Pairwise comparison:
AOA_16S.arc.rat.rh.emm <- emmeans(AOA_16S.arc.rat.rh.aov, ~ irrigation | fertilization*sampling.date)
AOA_16S.arc.rat.rh.pair <- pairs(AOA_16S.arc.rat.rh.emm)
AOA_16S.arc.rat.rh.DF <- as.data.frame(summary(AOA_16S.arc.rat.rh.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(AOA_16S.arc.rat.rh.DF, file = "AOA_16S.arc.rat.pair.rh.csv")

# Anova test for AOA/16S Ratio in percent

AOA_16S.ratio.rh.aov <- aov_ez("plot", "AOA_16S_ratio_percent", qPCR.RS, 
                                 within = "sampling.date",
                                 between = c("irrigation","fertilization"),
                                 #correction_aov = "GG",
                                 type = 2,
                                 return = afex_options("return_aov"),
                                 anova_table = list(correction="none"))
AOA_16S.ratio.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
AOA_16S.ratio.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = AOA_16S_ratio_percent, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(AOA_16S.ratio.rh.aov2)
AOA_16S.ratio.rh.aov2.df <- get_anova_table(AOA_16S.ratio.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOA_16S_percent_ratio")
table2csv(x=AOA_16S.ratio.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(AOA_16S.ratio.rh.aov) # good
check_sphericity(AOA_16S.ratio.rh.aov) # good
AOA_16S.ratio.rh.is_norm <- check_normality(AOA_16S.ratio.rh.aov)
shapiro.test(qPCR.RS$AOA_16S_ratio_percent) # NORMAL
plot(AOA_16S.ratio.rh.is_norm, type = "qq")
plot(AOA_16S.ratio.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOA_16S_ratio_percent", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(AOA_16S.ratio.rh.aov))
# Pairwise comparison:
AOA_16S.ratio.rh.emm <- emmeans(AOA_16S.ratio.rh.aov, ~ irrigation | fertilization*sampling.date)
AOA_16S.ratio.rh.pair <- pairs(AOA_16S.ratio.rh.emm)
AOA_16S.ratio.rh.DF <- as.data.frame(summary(AOA_16S.ratio.rh.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(AOA_16S.ratio.rh.DF, file = "AOA_16S_percent_ratio.pair.rh.csv")


##### 7. AOB/16S Ratio in Rhizosphere

# Anova test for arcsin square root transformed AOB/16S Ratio

AOB_16S.arc.rat.rh.aov <- aov_ez("plot", "AOB_16.arc.ratio.rh", qPCR.RS, 
                                 within = "sampling.date",
                                 between = c("irrigation","fertilization"),
                                 type = 3,
                                 return = afex_options("return_aov"),
                                 anova_table = list(correction="none"))
AOB_16S.arc.rat.rh.aov

str# Three-Way Mixed (Split-Plot) ANOVA 
AOB_16S.arc.rat.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = AOB_16.arc.ratio.rh, wid = plot, type = 3,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(AOB_16S.arc.rat.rh.aov2)
AOB_16S.arc.rat.rh.aov2.df <- get_anova_table(AOB_16S.arc.rat.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOB_16S_arc_ratio")
table2csv(x=AOB_16S.arc.rat.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(AOB_16S.arc.rat.rh.aov) # good
check_sphericity(AOB_16S.arc.rat.rh.aov) # good
AOB_16S.arc.rat.rh.is_norm <- check_normality(AOB_16S.arc.rat.rh.aov)
shapiro.test(qPCR.RS$AOB_16.arc.ratio.rh) # NORMAL
plot(AOB_16S.arc.rat.rh.is_norm, type = "qq")
plot(AOB_16S.arc.rat.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOB_16.arc.ratio.rh", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(AOB_16S.arc.rat.rh.aov))
# Pairwise comparison:
AOB_16S.arc.rat.rh.emm <- emmeans(AOB_16S.arc.rat.rh.aov, ~ irrigation | fertilization*sampling.date)
AOB_16S.arc.rat.rh.pair <- pairs(AOB_16S.arc.rat.rh.emm)
AOB_16S.arc.rat.rh.DF <- as.data.frame(summary(AOB_16S.arc.rat.rh.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(AOB_16S.arc.rat.rh.DF, file = "AOB_16S.arc.rat.pair.rh.csv")


# Anova test for AOB/16S Ratio in percent

AOB_16S.ratio.rh.aov <- aov_ez("plot", "AOB_16S_ratio_percent", qPCR.RS, 
                               within = "sampling.date",
                               between = c("irrigation","fertilization"),
                               type = 3,
                               return = afex_options("return_aov"),
                               anova_table = list())
AOB_16S.ratio.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
AOB_16S.ratio.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = AOB_16S_ratio_percent, wid = plot, type = 3,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(AOB_16S.ratio.rh.aov2)
AOB_16S.ratio.rh.aov2.df <- get_anova_table(AOB_16S.ratio.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOB_16S_percent_ratio")
table2csv(x=AOB_16S.ratio.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(AOB_16S.ratio.rh.aov) # good
check_sphericity(AOB_16S.ratio.rh.aov) # not good
AOB_16S.ratio.rh.is_norm <- check_normality(AOB_16S.ratio.rh.aov)
shapiro.test(qPCR.RS$AOB_16S_ratio_percent) # NORMAL
plot(AOB_16S.ratio.rh.is_norm, type = "qq")
plot(AOB_16S.ratio.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOB_16S_ratio_percent", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(AOB_16S.ratio.rh.aov))
# Pairwise comparison:
AOB_16S.ratio.rh.emm <- emmeans(AOB_16S.ratio.rh.aov, ~ irrigation | fertilization*sampling.date)
AOB_16S.ratio.rh.pair <- pairs(AOB_16S.ratio.rh.emm)
AOB_16S.ratio.rh.DF <- as.data.frame(summary(AOB_16S.ratio.rh.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(AOB_16S.ratio.rh.DF, file = "AOB_16S_percent_ratio.pair.rh.csv")


##### 7. ComA/16S Ratio in Rhizosphere

# Anova test for ComA/16S Ratio in percent

ComA_16S.ratio.rh.aov <- aov_ez("plot", "ComA_16S_ratio_percent", qPCR.RS, 
                               within = "sampling.date",
                               between = c("irrigation","fertilization"),
                               type = 2,
                               return = afex_options("return_aov"),
                               anova_table = list(correction="none"))
ComA_16S.ratio.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComA_16S.ratio.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = ComA_16S_ratio_percent, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(ComA_16S.ratio.rh.aov2)
ComA_16S.ratio.rh.aov2.df <- get_anova_table(ComA_16S.ratio.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/ComA_16S_percent_ratio")
table2csv(x=ComA_16S.ratio.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComA_16S.ratio.rh.aov) # SLIGHTLY NOT HOMOGEN
check_sphericity(ComA_16S.ratio.rh.aov) # good
ComA_16S.ratio.rh.is_norm <- check_normality(ComA_16S.ratio.rh.aov)
shapiro.test(qPCR.RS$ComA_16S_ratio_percent) # NOT NORMAL
plot(ComA_16S.ratio.rh.is_norm, type = "density")
plot(ComA_16S.ratio.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "ComA_16S_ratio_percent", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(ComA_16S.ratio.rh.aov))
# Pairwise comparison:
ComA_16S.ratio.rh.emm <- emmeans(ComA_16S.ratio.rh.aov, ~ irrigation | fertilization*sampling.date)
ComA_16S.ratio.rh.pair <- pairs(ComA_16S.ratio.rh.emm)
ComA_16S.ratio.rh.DF <- as.data.frame(summary(ComA_16S.ratio.rh.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(ComA_16S.ratio.rh.DF, file = "ComA_16S_percent_ratio.pair.rh.csv")


# Anova test for arcsin square root transformed ComA16S Ratio

ComA_16S.arc.rat.rh.aov <- aov_ez("plot", "ComA_16.arc.ratio.rh", qPCR.RS, 
                                 within = "sampling.date",
                                 between = c("irrigation","fertilization"),
                                 type = 2,
                                 return = afex_options("return_aov"),
                                 anova_table = list(correction="none"))
ComA_16S.arc.rat.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComA_16S.arc.rat.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = ComA_16.arc.ratio.rh, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(ComA_16S.arc.rat.rh.aov2)
ComA_16S.arc.rat.rh.aov2.df <- get_anova_table(ComA_16S.arc.rat.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/ComA_16S_arc_ratio")
table2csv(x=ComA_16S.arc.rat.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComA_16S.arc.rat.rh.aov) # good
check_sphericity(ComA_16S.arc.rat.rh.aov) # good
ComA_16S.arc.rat.rh.is_norm <- check_normality(ComA_16S.arc.rat.rh.aov)
shapiro.test(qPCR.RS$ComA_16.arc.ratio.rh) # NORMAL
plot(ComA_16S.arc.rat.rh.is_norm, type = "density")
plot(ComA_16S.arc.rat.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "ComA_16.arc.ratio.rh", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(ComA_16S.arc.rat.rh.aov))
# Pairwise comparison:
ComA_16S.arc.rat.rh.emm <- emmeans(ComA_16S.arc.rat.rh.aov, ~ irrigation | fertilization*sampling.date)
ComA_16S.arc.rat.rh.pair <- pairs(ComA_16S.arc.rat.rh.emm)
ComA_16S.arc.rat.rh.DF <- as.data.frame(summary(ComA_16S.arc.rat.rh.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(ComA_16S.arc.rat.rh.DF, file = "ComA_16S.arc.rat.pair.rh.csv")

##### 8. ComB/16S Ratio in Rhizosphere

# Anova test for ComB/16S Ratio in percent

ComB_16S.ratio.rh.aov <- aov_ez("plot", "ComB_16S_ratio_percent", qPCR.RS, 
                                within = "sampling.date",
                                between = c("irrigation","fertilization"),
                                type = 2,
                                return = afex_options("return_aov"),
                                anova_table = list(correction="none"))
ComB_16S.ratio.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComB_16S.ratio.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = ComB_16S_ratio_percent, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(ComB_16S.ratio.rh.aov2)
ComB_16S.ratio.rh.aov2.df <- get_anova_table(ComB_16S.ratio.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/ComB_16S_percent_ratio")
table2csv(x=ComB_16S.ratio.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComB_16S.ratio.rh.aov) # good
check_sphericity(ComB_16S.ratio.rh.aov) # good
ComB_16S.ratio.rh.is_norm <- check_normality(ComB_16S.ratio.rh.aov)
shapiro.test(qPCR.RS$ComB_16S_ratio_percent) # NORMAL
plot(ComB_16S.ratio.rh.is_norm, type = "density")
plot(ComB_16S.ratio.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "ComB_16S_ratio_percent", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(ComB_16S.ratio.rh.aov))
# Pairwise comparison:
ComB_16S.ratio.rh.emm <- emmeans(ComB_16S.ratio.rh.aov, ~ irrigation | fertilization*sampling.date)
ComB_16S.ratio.rh.pair <- pairs(ComB_16S.ratio.rh.emm)
ComB_16S.ratio.rh.DF <- as.data.frame(summary(ComB_16S.ratio.rh.pair))
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(ComB_16S.ratio.rh.DF, file = "ComB_16S_percent_ratio.pair.rh.csv")

##### 9. AOA/AOB Ratio in Rhizosphere

# Anova test for AOA/AOB Ratio on non transformed data

AOA_AOB.ratio.rh.aov <- aov_ez("plot", "AOA_AOB_ratio", qPCR.RS, 
                                within = "sampling.date",
                                between = c("irrigation","fertilization"),
                                type = 2,
                                return = afex_options("return_aov"),
                                anova_table = list(correction="none"))
AOA_AOB.ratio.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
AOA_AOB.ratio.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = AOA_AOB_ratio_percent, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(AOA_AOB.ratio.rh.aov2)
AOA_AOB.ratio.rh.aov2.df <- get_anova_table(AOA_AOB.ratio.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOA_AOB_ratio")
table2csv(x=AOA_AOB.ratio.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(AOA_AOB.ratio.rh.aov) # good
check_sphericity(AOA_AOB.ratio.rh.aov) # good
AOA_AOB.ratio.rh.is_norm <- check_normality(AOA_AOB.ratio.rh.aov)
shapiro.test(qPCR.RS$AOA_AOB_ratio_percent) # NOT NORMAL
plot(AOA_AOB.ratio.rh.is_norm, type = "density")
plot(AOA_AOB.ratio.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOA_AOB_ratio_percent", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(AOA_AOB.ratio.rh.aov))
# Pairwise comparison:
AOA_AOB.ratio.rh.emm <- emmeans(AOA_AOB.ratio.rh.aov, ~ irrigation | fertilization*sampling.date)
AOA_AOB.ratio.rh.pair <- pairs(AOA_AOB.ratio.rh.emm)
AOA_AOB.ratio.rh.DF <- as.data.frame(summary(AOA_AOB.ratio.rh.pair))
AOA_AOB.ratio.rh.DF
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(AOA_AOB.ratio.rh.DF, file = "AOA_AOB_ratio.pair.rh.csv")


# Anova test for AOA/AOB Ratio on arcsin square root transformed data

AOA_AOB.arcs.ratio.rh.aov <- aov_ez("plot", "AOA_AOB.arc.ratio.rh", qPCR.RS, 
                               within = "sampling.date",
                               between = c("irrigation","fertilization"),
                               type = 2,
                               return = afex_options("return_aov"),
                               anova_table = list(correction="none"))
AOA_AOB.arcs.ratio.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
AOA_AOB.arcs.ratio.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = AOA_AOB.arc.ratio.rh, wid = plot, type = 2,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(AOA_AOB.arcs.ratio.rh.aov2)
AOA_AOB.arcs.ratio.rh.aov2.df <- get_anova_table(AOA_AOB.arcs.ratio.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/AOA_AOB_arc_ratio")
table2csv(x=AOA_AOB.arcs.ratio.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(AOA_AOB.arcs.ratio.rh.aov) # good
check_sphericity(AOA_AOB.arcs.ratio.rh.aov) # good
AOA_AOB.arcs.ratio.rh.is_norm <- check_normality(AOA_AOB.arcs.ratio.rh.aov)
shapiro.test(qPCR.RS$AOA_AOB.arc.ratio.rh) # NORMAL
plot(AOA_AOB.arcs.ratio.rh.is_norm, type = "density")
plot(AOA_AOB.arcs.ratio.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "AOA_AOB.arc.ratio.rh", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(AOA_AOB.arcs.ratio.rh.aov))
# Pairwise comparison:
AOA_AOB.arcs.ratio.rh.emm <- emmeans(AOA_AOB.arcs.ratio.rh.aov, ~ irrigation | fertilization*sampling.date)
AOA_AOB.arcs.ratio.rh.pair <- pairs(AOA_AOB.arcs.ratio.rh.emm)
AOA_AOB.arcs.ratio.rh.DF <- as.data.frame(summary(AOA_AOB.arcs.ratio.rh.pair))
AOA_AOB.arcs.ratio.rh.DF
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(AOA_AOB.arcs.ratio.rh.DF, file = "AOA_AOB.arcs.ratio.pair.rh.csv")


##### 9. ComA/ComB Ratio in Rhizosphere

# Anova test for ComA/ComB Ratio on non transformed data

ComA_ComB.ratio.rh.aov <- aov_ez("plot", "ComA_ComB_ratio", qPCR.RS, 
                               within = "sampling.date",
                               between = c("irrigation","fertilization"),
                               type = 3,
                               return = afex_options("return_aov"),
                               anova_table = list(correction="none"))
ComA_ComB.ratio.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComA_ComB.ratio.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = ComA_ComB_ratio_percent, wid = plot, type = 3,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(ComA_ComB.ratio.rh.aov2)
ComA_ComB.ratio.rh.aov2.df <- get_anova_table(ComA_ComB.ratio.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/ComA_ComB_ratio")
table2csv(x=ComA_ComB.ratio.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComA_ComB.ratio.rh.aov) # good
check_sphericity(ComA_ComB.ratio.rh.aov) # good
ComA_ComB.ratio.rh.is_norm <- check_normality(ComA_ComB.ratio.rh.aov)
shapiro.test(qPCR.RS$ComA_ComB_ratio_percent) # NORMAL
plot(ComA_ComB.ratio.rh.is_norm, type = "density")
plot(ComA_ComB.ratio.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "ComA_ComB_ratio_percent", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(ComA_ComB.ratio.rh.aov))
# Pairwise comparison:
ComA_ComB.ratio.rh.emm <- emmeans(ComA_ComB.ratio.rh.aov, ~ irrigation | fertilization*sampling.date)
ComA_ComB.ratio.rh.pair <- pairs(ComA_ComB.ratio.rh.emm)
ComA_ComB.ratio.rh.DF <- as.data.frame(summary(ComA_ComB.ratio.rh.pair))
ComA_ComB.ratio.rh.DF
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(ComA_ComB.ratio.rh.DF, file = "ComA_ComB_ratio.pair.rh.csv")


# Anova test for ComA_ComB Ratio on arcsin square root transformed data

ComA_ComB.arcs.ratio.rh.aov <- aov_ez("plot", "ComA_ComB.arc.ratio.rh", qPCR.RS, 
                                    within = "sampling.date",
                                    between = c("irrigation","fertilization"),
                                    type = 3,
                                    return = afex_options("return_aov"),
                                    anova_table = list(correction="none"))
ComA_ComB.arcs.ratio.rh.aov

# Three-Way Mixed (Split-Plot) ANOVA 
ComA_ComB.arcs.ratio.rh.aov2 <- anova_test(
  data = qPCR.RS, dv = ComA_ComB.arc.ratio.rh, wid = plot, type = 3,
  within = sampling.date, between = c(irrigation, fertilization))
get_anova_table(ComA_ComB.arcs.ratio.rh.aov2)
ComA_ComB.arcs.ratio.rh.aov2.df <- get_anova_table(ComA_ComB.arcs.ratio.rh.aov2)
#setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/')
filen <- paste("D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere/ComA_ComB_arc_ratio")
table2csv(x=ComA_ComB.arcs.ratio.rh.aov2.df,file=filen, digits = 1, digitspvals = 3)

# test assumptions:
check_homogeneity(ComA_ComB.arcs.ratio.rh.aov) # good
check_sphericity(ComA_ComB.arcs.ratio.rh.aov) # good
ComA_ComB.arcs.ratio.rh.is_norm <- check_normality(ComA_ComB.arcs.ratio.rh.aov)
shapiro.test(qPCR.RS$ComA_ComB.arc.ratio.rh) # NORMAL
plot(ComA_ComB.arcs.ratio.rh.is_norm, type = "density")
plot(ComA_ComB.arcs.ratio.rh.is_norm, type = "qq", detrend = TRUE)
ggqqplot(qPCR.RS, "ComA_ComB.arc.ratio.rh", ggtheme = theme_bw())

# tidy anova table
knitr::kable(nice(ComA_ComB.arcs.ratio.rh.aov))
# Pairwise comparison:
ComA_ComB.arcs.ratio.rh.emm <- emmeans(ComA_ComB.arcs.ratio.rh.aov, ~ irrigation | fertilization*sampling.date)
ComA_ComB.arcs.ratio.rh.pair <- pairs(ComA_ComB.arcs.ratio.rh.emm)
ComA_ComB.arcs.ratio.rh.DF <- as.data.frame(summary(ComA_ComB.arcs.ratio.rh.pair))
ComA_ComB.arcs.ratio.rh.DF
setwd('D:/Fina/INRAE_Project/microservices/qPCR_stat/Rhizosphere')
write.csv(ComA_ComB.arcs.ratio.rh.DF, file = "ComA_ComB.arcs.ratio.pair.rh.csv")


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

