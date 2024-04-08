###########################################################################
#################### ALPHA DIVERSITY ANALYSIS: COM ########################
###########################################################################

# Author: Ari Fina Bintarti
# Date: 02/08/2023

#install.packages("datarium")
#install.packages("rstatix")
library(datarium)
library(rstatix)
###########################################################################
# 1. Response variable: Richness
###########################################################################
# 1a. Analyses of Bulk Soil
com.meta.bulk <- com.meta.df.ed[1:118,]
str(com.meta.bulk)
com.meta.bulk$SampleID<-factor(com.meta.bulk$SampleID)
com.meta.bulk$PlotID<-factor(com.meta.bulk$PlotID)
com.meta.bulk$Irrigation<-factor(com.meta.bulk$Irrigation)
com.meta.bulk$Block<-factor(com.meta.bulk$Block)
com.meta.bulk$x<-factor(com.meta.bulk$x)
com.meta.bulk$rep<-factor(com.meta.bulk$rep)
com.meta.bulk$rep2<-factor(com.meta.bulk$rep2)
com.meta.bulk$var2<-factor(com.meta.bulk$var2)
com.meta.bulk$var3<-factor(com.meta.bulk$var3)
com.meta.bulk$sqrtRich <- sqrt(com.meta.bulk$Richness)

com.meta.bulk.sum.rich <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
com.bulk.sum.rich.plot <- ggboxplot(
  com.meta.bulk, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.bulk.sum.rich.plot
# check assumption (outliers)
com.bulk.rich.out <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
com.bulk.rich.SW <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(com.meta.bulk, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.bulk.rich.Lave <- com.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Repeated-measured ANOVA 
set.seed(13)
com.bulk.rich.aov2 <- anova_test(
  data = com.meta.bulk, type=3, dv = Richness, wid = rep2,
  #between=Date,
  within = c(Irrigation,Treatment,Date))
get_anova_table(com.bulk.rich.aov2)

# Two-way ANOVA at each fertilization level
com.bulk.rich.aov2way <-com.meta.bulk %>%
  group_by(Date) %>%
  anova_test(dv = Richness, wid = rep2, within = c(Irrigation, Treatment), type=3)
get_anova_table(com.bulk.rich.aov2way)
# Effect of drought at each fert X date
com.drought.effect <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  anova_test(dv = Richness, wid = rep2, within = Irrigation)
com.drought.effect
# Pairwise comparisons
com.BS.pwc <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  pairwise_t_test(Richness ~ Irrigation, paired = F, p.adjust.method = "BH") 
com.BS.pwc

com.BS.rich.aov2 <- aov(Richness ~ Irrigation*Treatment*Date+
                         Error(rep2/(Irrigation*Treatment*Date)), data=com.meta.bulk)
summary(com.BS.rich.aov2)

# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.pwc.trt <- com.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95)
View(com.pwc.trt)
# 2. between irrigation:
com.pwc.irr <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95)
View(com.pwc.irr)
# lmer model
com.rich.bulk.lmer <- lmerTest::lmer(com.meta.bulk$Richness ~ Irrigation*Treatment*Date+(1|Block)+(1|Block:Date),
                                     data=com.meta.bulk,contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"))
#boundary (singular) fit: see help('isSingular')
car::Anova(com.rich.bulk.lmer, test="F", type="III") 
# test assumptions:
library(DHARMa)
shapiro.test(resid(com.rich.bulk.lmer)) # normal
plot(simulateResiduals(com.rich.bulk.lmer)) # okay

# Test only 4 time points:
com.meta.bulk.4 <- com.meta.bulk[1:94,]
str(com.meta.bulk.4)
com.meta.bulk.4$SampleID<-factor(com.meta.bulk.4$SampleID)
com.meta.bulk.4$PlotID<-factor(com.meta.bulk.4$PlotID)
com.meta.bulk.4$Irrigation<-factor(com.meta.bulk.4$Irrigation)
com.meta.bulk.4$x<-factor(com.meta.bulk.4$x)
com.meta.bulk.4$rep<-factor(com.meta.bulk.4$rep)
com.meta.bulk.4$rep2<-factor(com.meta.bulk.4$rep2)
com.meta.bulk.4$var2<-factor(com.meta.bulk.4$var2)
com.meta.bulk.4$var3<-factor(com.meta.bulk.4$var3)
com.meta.bulk.4$Date<-factor(com.meta.bulk.4$Date)

# check assumption (outliers)
com.bulk.rich.out.4 <- com.meta.bulk.4 %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(sqrtRich) # no extreme outliers
View(com.bulk.rich.out.4)
# Saphiro-Wilk for normality
com.bulk.rich.SW.4 <- com.meta.bulk.4 %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(sqrtRich)
View(com.bulk.rich.SW.4)
ggqqplot(com.meta.bulk.4, "sqrtRich", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.bulk.rich.Lave.4 <- com.meta.bulk.4 %>%
  group_by(Date) %>%
  levene_test(sqrtRich ~ Irrigation*Treatment)
View(com.bulk.rich.Lave.4)

set.seed(13)
# Three-way repeated measured anova
com.bulk.rich.aov.4 <- anova_test(
  data = com.meta.bulk.4, type=3, dv = Richness, wid = rep2,
  #between=Date,
  within = c(Irrigation, Treatment, Date))
get_anova_table(com.bulk.rich.aov.4)
# lmer model
com.rich.bulk.lmer.4 <- lmerTest::lmer(com.meta.bulk.4$Richness ~ Irrigation*Treatment*Date+(1|Block)+(1|Block:Date),
                                     data=com.meta.bulk.4,contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"))
#boundary (singular) fit: see help('isSingular')
car::Anova(com.rich.bulk.lmer.4, test="F", type="III") 
# test assumptions:
library(DHARMa)
shapiro.test(resid(com.rich.bulk.lmer.4)) # normal
plot(simulateResiduals(com.rich.bulk.lmer.4)) # okay
# Two-way ANOVA 
two.way <- com.meta.bulk %>%
  group_by(Date) %>%
  anova_test(dv = Richness, wid = rep, within = c(Irrigation, Treatment))
get_anova_table(two.way)

pwc <- com.meta.bulk %>%
  group_by(Treatment, Date) %>%
  emmeans_test(Richness ~ Irrigation, p.adjust.method = "BH") %>%
  select(-df, -statistic, -p) # Remove details


# Test # Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.pwc.trt.4 <- com.meta.bulk.4 %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95)
View(com.pwc.trt.4)
# 2. between irrigation:
com.pwc.irr.4 <- com.meta.bulk.4 %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95)
View(com.pwc.irr.4)

#Method 3 
#model3 <- lme(Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, method="REML", data=com.meta.bulk)
#anova(model3)
#model1 <- aov(com.meta.bulk$Richness ~ Irrigation*Treatment*Date + Error(PlotID/Irrigation*Treatment*Date) , data=com.meta.bulk)
summary(model1)

############################################################################################################
# Model Fit
set.seed(13)
#install.packages("lme4", type = "source")
library(lme4)
library(lmerTest)
com.rich.bulk.lmer <- lmerTest::lmer(Richness ~ Irrigation*Treatment*Date +(1|Block)+(1|Block:Date),
                                    contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"),data=com.meta.bulk)
#boundary (singular) fit: see help('isSingular')
Anova(com.rich.bulk.lmer, test="F", type="III") 

com.rich.bulk.mod2 <- lmerTest::lmer(com.meta.bulk$Richness ~ Irrigation*Treatment*Date+(1|rep), contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"),
                                     data=com.meta.bulk)
anova(com.rich.bulk.mod2)
Anova(com.rich.bulk.mod2, test="F", type="III") 
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.rich.bulk <- com.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.rich.bulk.mod)
# 2. between irrigation:
com.emm.rich.irri.bulk <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.rich.bulk.mod2)
com.emm.rich.irri.bulk
# Comammox Richness Summary
com.rich.sum <- com.meta.bulk %>%
  group_by(Irrigation, Treatment) %>%
  get_summary_stats(Richness, type = "mean_sd")
View(com.rich.sum)
com.rich.sum.date <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
View(com.rich.sum.date)
com.rich.sum.plot <- ggboxplot(
  com.meta.bulk, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.rich.sum.plot

######################################################################
# 1b. Analyses of rhizosphere Soil

com.meta.rh <- com.meta.df.ed[119:190,]
view(com.meta.rh)
com.meta.rh$SampleID<-factor(com.meta.rh$SampleID)
com.meta.rh$PlotID<-factor(com.meta.rh$PlotID)
com.meta.rh$Irrigation<-factor(com.meta.rh$Irrigation)
com.meta.rh$Block<-factor(com.meta.rh$Block)
com.meta.rh$x<-factor(com.meta.rh$x)
com.meta.rh$rep<-factor(com.meta.rh$rep)
com.meta.rh$rep2<-factor(com.meta.rh$rep2)
com.meta.rh$var2<-factor(com.meta.rh$var2)
com.meta.rh$var3<-factor(com.meta.rh$var3)
com.meta.rh$Date <- factor(com.meta.rh$Date)

com.meta.rh.sum.rich <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
com.rh.sum.rich.plot <- ggboxplot(
  com.meta.rh, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.rh.sum.rich.plot
# check assumption (outliers)
com.rh.rich.out <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
com.rh.rich.SW <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(com.meta.rh, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.rh.rich.Lave <- com.meta.rh %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
com.rh.rich.aov <- anova_test(
  data = com.meta.rh, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.rh.rich.aov)

# Three-way repeated measured anova
com.rh.rich.aov <- anova_test(
  data = com.meta.rh, type=3, dv = Richness, wid = rep2,
  #between=Date,
  within = c(Irrigation, Treatment, Date))
get_anova_table(com.rh.rich.aov)

COM.RS.Rich.aov2 <- aov(Richness ~ Irrigation*Treatment*Date + Error(rep2/(Irrigation*Treatment*Date)), data=com.meta.rh)
summary(COM.RS.Rich.aov2)

# Model Fit
set.seed(13)
# Test Method 3 
com.RS.rich.lmer <- lmerTest::lmer(com.meta.rh$Richness ~ Irrigation*Treatment*Date +(1|Block)+(1|Block:Date), 
                           data=com.meta.rh,contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"))
car::Anova(com.RS.rich.lmer, test="F", type="III") 

# pairwise comparisons
# 1. between fertilization treatment:
com.emm.rich.rh <- com.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.rich.rhizo.mod)
# 2. between irrigation:
com.emm.rich.irri.rh <- com.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.RS.rich.lmer)
com.emm.rich.irri.rh 

###########################################################################
# 2. Response variable: Shannon
###########################################################################
# 2a. Analyses of Bulk Soil
com.meta.bulk.sum.sha <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
com.bulk.sum.sha.plot <- ggboxplot(
  com.meta.bulk, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.bulk.sum.sha.plot
# check assumption (outliers)
com.bulk.sha.out <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
View(com.bulk.sha.out)
# Saphiro-Wilk for normality
com.bulk.sha.SW <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
View(com.bulk.sha.SW)
ggqqplot(com.meta.bulk, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.bulk.sha.Lave <- com.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
View(com.bulk.sha.Lave)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
com.bulk.sha.aov <- anova_test(
  data = com.meta.bulk, type=2, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.bulk.sha.aov)
# Three-way repeated measured anova
com.bulk.sha.aov <- anova_test(
  data = com.meta.bulk, type=3, dv = Shannon, wid = rep2,
  #between=Date,
  within = c(Irrigation, Treatment, Date))
get_anova_table(com.bulk.sha.aov)
# Two-way ANOVA at each fertilization level
com.bulk.sha.aov2way <-com.meta.bulk %>%
  group_by(Date) %>%
  anova_test(dv = Shannon, wid = rep2, within = c(Irrigation, Treatment), type=3)
get_anova_table(com.bulk.sha.aov2way)
# Effect of drought at each fert X date
com.sha.drought.effect <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  anova_test(dv = Shannon, wid = rep2, within = Irrigation)
com.sha.drought.effect
# Pairwise comparisons
com.sha.BS.pwc <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  pairwise_t_test(Shannon ~ Irrigation, paired = F, p.adjust.method = "BH") 
com.sha.BS.pwc



# Analysis of 4 time points
# check assumption (outliers)
com.bulk.sha.out.4 <- com.meta.bulk.4 %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
View(com.bulk.sha.out.4)
# Saphiro-Wilk for normality
com.bulk.sha.SW.4 <- com.meta.bulk.4 %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
View(com.bulk.sha.SW.4)
ggqqplot(com.meta.bulk.4, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.bulk.sha.Lave.4 <- com.meta.bulk.4 %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
View(com.bulk.sha.Lave.4)

set.seed(13)
# Three-way repeated measured anova for 4 time points
com.bulk.sha.aov.4 <- anova_test(
  data = com.meta.bulk.4, type=3, dv = Shannon, wid = rep2,
  #between=Date,
  within = c(Irrigation, Treatment, Date))
get_anova_table(com.bulk.sha.aov.4)



############################################################################################################
# Model Fit
set.seed(13)
# lmer model
com.sha.bulk.lmer <- lmerTest::lmer(com.meta.bulk$Shannon ~ Irrigation*Treatment*Date+(1|Block)+(1|Block:Date),
                                     data=com.meta.bulk,contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"))
#boundary (singular) fit: see help('isSingular')
car::Anova(com.sha.bulk.lmer, test="F", type="III") 
# lmer model 4 time pointss
com.sha.bulk.lmer.4 <- lmerTest::lmer(com.meta.bulk.4$Shannon ~ Irrigation*Treatment*Date+(1|Block)+(1|Block:Date),
                                     data=com.meta.bulk.4,contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"))
#boundary (singular) fit: see help('isSingular')
car::Anova(com.sha.bulk.lmer.4, test="F", type="III") 
# test assumptions:
library(DHARMa)
hist(com.meta.bulk$Shannon) # looks fine
shapiro.test(resid(com.sha.bulk.lmer)) # slightly not normal
plot(simulateResiduals(com.sha.bulk.lmer)) # okay
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.sha.bulk <- com.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Shannon ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.sha.bulk.mod)
# 2. between irrigation:
com.emm.sha.irri.bulk <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Shannon ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95,model=com.sha.bulk.lmer)
com.emm.sha.irri.bulk

# Comammox Shannon Summary
com.sha.sum <- com.meta.bulk %>%
  group_by(Irrigation, Treatment) %>%
  get_summary_stats(Shannon, type = "mean_sd")
View(com.sha.sum)
com.sha.sum.date <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
View(com.sha.sum.date)
com.sha.sum.plot <- ggboxplot(
  com.meta.bulk, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.sha.sum.plot
#############################################################################################################

##################################################################
# 2b. Analyses of Rhizosphere
com.meta.rh.sum.sha <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
com.rh.sum.sha.plot <- ggboxplot(
  com.meta.rh, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.rh.sum.sha.plot
# check assumption (outliers)
com.rh.sha.out <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
# Saphiro-Wilk for normality
com.rh.sha.SW <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
ggqqplot(com.meta.rh, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.rh.sha.Lave <- com.min.meta.rh %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(com.meta.rh)
com.rh.sha.aov <- anova_test(
  data = com.meta.rh, type=3, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.rh.sha.aov)


# Three-way repeated measured anova
com.rh.sha.aov <- anova_test(
  data = com.meta.rh, type=3, dv = Shannon, wid = rep2,
  #between=Date,
  within = c(Irrigation, Treatment, Date))
get_anova_table(com.rh.sha.aov)
############################################################################################################
# Model Fit
set.seed(13)
com.sha.rh.mod <- lmerTest::lmer(com.meta.rh$Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.rh)
anova(com.sha.rh.mod, type = 2)

com.sha.rh.mod2 <- lmerTest::lmer(com.meta.rh$Shannon ~ Irrigation*Treatment*Date+(1|Block:Date), data=com.meta.rh, na.action=na.omit)
anova(com.sha.rh.mod2)

# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.sha.rh <- com.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Shannon ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.sha.rh.mod2)
# 2. between irrigation:
com.emm.sha.irri.rh <- com.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Shannon ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.sha.rh.mod2)
com.emm.sha.irri.rh

#############################################################################################################

###########################################################################
# 3. Response variable: Inverse Simpson
###########################################################################
# 4a. Analyses of Bulk Soil
com.meta.bulk.invsum.simp <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
com.bulk.sum.invsimp.plot <- ggboxplot(
  com.meta.bulk, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.bulk.sum.invsimp.plot
# check assumption (outliers)
com.bulk.invsimp.out <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
com.bulk.invsimp.SW <- com.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(com.meta.bulk, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.bulk.invsimp.Lave <- com.meta.bulk %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
com.bulk.invsimp.aov <- anova_test(
  data = com.meta.bulk, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.bulk.invsimp.aov)

# Three-Way Repeated-Measures ANOVA
com.bulk.simp.aov <- anova_test(
  data = com.meta.bulk, type=3,dv = InvSimpson, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(com.bulk.simp.aov)


# Three-Way Repeated-Measures ANOVA 4 tIME pOINTS
com.bulk.simp.aov.4 <- anova_test(
  data = com.meta.bulk.4, type=3,dv = InvSimpson, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(com.bulk.simp.aov.4)


############################################################################################################
# Model Fit
set.seed(13)
com.invsimp.bulk.mod <- lmerTest::lmer(com.meta.bulk$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.bulk)
anova(com.invsimp.bulk.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.invsimp.bulk <- com.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(InvSimpson ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.invsimp.bulk.mod)
# 2. between irrigation:
com.emm.invsimp.irri.bulk <- com.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(InvSimpson ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.invsimp.bulk.mod)
#############################################################################################################

##################################################################
# 3b. Analyses of Rhizosphere
com.meta.rh.sum.invsimp <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
com.rh.sum.invsimp.plot <- ggboxplot(
  com.meta.rh, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
com.rh.sum.invsimp.plot
# check assumption (outliers)
com.rh.invsimp.out <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
com.rh.invsimp.SW <- com.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(com.meta.rh, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
com.rh.invsimp.Lave <- com.meta.rh %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
com.rh.invsimp.aov <- anova_test(
  data = com.min.meta.rh, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(com.rh.invsimp.aov)


# Three-way repeated measured anova
com.rh.simp.aov <- anova_test(
  data = com.meta.rh, type=3, dv = InvSimpson, wid = rep2,
  #between=Date,
  within = c(Irrigation, Treatment, Date))
get_anova_table(com.rh.simp.aov)
############################################################################################################
# Model Fit
set.seed(13)
com.invsimp.rh.mod <- lmerTest::lmer(com.meta.rh$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=com.meta.rh)
anova(com.invsimp.rh.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
com.emm.invsimp.rh <- com.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(InvSimpson ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.invsimp.rh.mod)
# 2. between irrigation:
com.emm.invsimp.irri.rh <- com.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(InvSimpson ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = com.invsimp.rh.mod)
#############################################################################################################













