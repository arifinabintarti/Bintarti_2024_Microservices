###########################################################################
#################### ALPHA DIVERSITY ANALYSIS: AOA ########################
###########################################################################

# Author: Ari Fina Bintarti
# Date: 28/07/2023

#install.packages("datarium")
#install.packages("rstatix")
library(datarium)
library(rstatix)
library(DHARMa)
###########################################################################
# 1. Response variable: Richness
###########################################################################
# 1a. Analyses of Bulk Soil
aoa.meta.bulk <- aoa.meta.df[1:120,]
str(aoa.meta.bulk)
aoa.meta.bulk$SampleID<-factor(aoa.meta.bulk$SampleID)
aoa.meta.bulk$PlotID<-factor(aoa.meta.bulk$PlotID)
aoa.meta.bulk$Irrigation<-factor(aoa.meta.bulk$Irrigation)
aoa.meta.bulk$Block<-factor(aoa.meta.bulk$Block)
aoa.meta.bulk$x<-factor(aoa.meta.bulk$x)
aoa.meta.bulk$rep<-factor(aoa.meta.bulk$rep)
aoa.meta.bulk$rep2<-factor(aoa.meta.bulk$rep2)
aoa.meta.bulk$var2<-factor(aoa.meta.bulk$var2)
aoa.meta.bulk$var3<-factor(aoa.meta.bulk$var3)
aoa.meta.bulk$AOA_sqrtRich <- sqrt(aoa.meta.bulk$Richness)


aoa.meta.bulk.sum.rich <- aoa.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
aoa.bulk.sum.rich.plot <- ggboxplot(
  aoa.meta.bulk, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aoa.bulk.sum.rich.plot
# check assumption (outliers)
aoa.bulk.rich.out <- aoa.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
aoa.bulk.rich.SW <- aoa.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(aoa.meta.bulk, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.bulk.rich.Lave <- aoa.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aoa.bulk.rich.aov <- anova_test(
  data = aoa.meta.bulk, type=2, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aoa.bulk.rich.aov)

# Three-Way Repeated-Measures ANOVA
aoa.bulk.rich.aov <- anova_test(
  data = aoa.meta.bulk, type=2,dv = Richness, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.bulk.rich.aov)

# Test only 4 time points:
aoa.meta.bulk.4 <- aoa.meta.bulk[1:96,]
str(aoa.meta.bulk.4)
aoa.meta.bulk.4$SampleID<-factor(aoa.meta.bulk.4$SampleID)
aoa.meta.bulk.4$PlotID<-factor(aoa.meta.bulk.4$PlotID)
aoa.meta.bulk.4$Irrigation<-factor(aoa.meta.bulk.4$Irrigation)
aoa.meta.bulk.4$Block<-factor(aoa.meta.bulk.4$Block)
aoa.meta.bulk.4$Date<-factor(aoa.meta.bulk.4$Date)
aoa.meta.bulk.4$x<-factor(aoa.meta.bulk.4$x)
aoa.meta.bulk.4$rep<-factor(aoa.meta.bulk.4$rep)
aoa.meta.bulk.4$rep2<-factor(aoa.meta.bulk.4$rep2)
aoa.meta.bulk.4$var2<-factor(aoa.meta.bulk.4$var2)
aoa.meta.bulk.4$var3<-factor(aoa.meta.bulk.4$var3)

# Three-Way Repeated-Measures ANOVA
aoa.bulk.rich.aov.4 <- anova_test(
  data = aoa.meta.bulk.4, type=2,dv = Richness, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.bulk.rich.aov.4)



# Test Method 3 
#model3 <- lme(Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, method="REML", data=aoa.meta.bulk)
#anova(model3)
#model1 <- aov(aob.meta.bulk$Richness ~ Irrigation*Treatment*Date + Error(PlotID/Irrigation*Treatment) , data=aoa.meta.bulk)
#summary(model1)

############################################################################################################
# Model Fit
library(lmerTest)
set.seed(13)
aoa.rich.bulk.mod <- lmerTest::lmer(aoa.meta.bulk$Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=aoa.meta.bulk)
anova(aoa.rich.bulk.mod, type = 3)

aoa.rich.bulk.mod2 <- lmerTest::lmer(aoa.meta.bulk$Richness ~ Irrigation*Treatment*Date+(1|Block)+(1|Block:Date),
                                     data=aoa.meta.bulk,contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"))
car::Anova(aoa.rich.bulk.mod2, test="F", type="III") 
# test assumptions:
library(DHARMa)
shapiro.test(resid(aoa.rich.bulk.mod2)) # not normal
plot(simulateResiduals(aoa.rich.bulk.mod2)) # okay
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
aoa.emm.rich.bulk <- aoa.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.rich.bulk.mod2)

# 2. between irrigation:
aoa.emm.rich.irri.bulk <- aoa.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.rich.bulk.mod2)
aoa.emm.rich.irri.bulk
# Pairwise comparison:
library(emmeans)
aoa.rich.emm <- emmeans(aoa.rich.bulk.mod2, ~ Treatment)
aoa.rich.emm.pair <- pairs(aoa.rich.emm)
aoa.rich.emm.pair
aoa.rich.emm.pair.DF <- as.data.frame(summary(aoa.rich.emm.pair))
######################################################################
# 1b. Analyses of rhizosphere Soil

aoa.meta.rh <- aoa.meta.df[121:192,]
str(aoa.meta.rh)
View(aoa.meta.rh)
aoa.meta.rh$SampleID<-factor(aoa.meta.rh$SampleID)
aoa.meta.rh$PlotID<-factor(aoa.meta.rh$PlotID)
aoa.meta.rh$Irrigation<-factor(aoa.meta.rh$Irrigation)
aoa.meta.rh$Block<-factor(aoa.meta.rh$Block)
aoa.meta.rh$x<-factor(aoa.meta.rh$x)
aoa.meta.rh$rep<-factor(aoa.meta.rh$rep)
aoa.meta.rh$rep2<-factor(aoa.meta.rh$rep2)
aoa.meta.rh$var2<-factor(aoa.meta.rh$var2)
aoa.meta.rh$var3<-factor(aoa.meta.rh$var3)
aoa.meta.rh$Date <- factor(aoa.meta.rh$Date)

aoa.meta.rh.sum.rich <- aoa.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
aoa.rh.sum.rich.plot <- ggboxplot(
  aoa.meta.rh, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aoa.rh.sum.rich.plot
# check assumption (outliers)
aoa.rh.rich.out <- aoa.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
aoa.rh.rich.SW <- aoa.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(aoa.meta.rh, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.rh.rich.Lave <- aoa.meta.rh %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aoa.rh.rich.aov <- anova_test(
  data = aoa.meta.rh, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aoa.rh.rich.aov)


# Three-Way Repeated-Measures ANOVA
aoa.rich.RS.aov1 <- anova_test(
  data = aoa.meta.rh, type=2,dv = Richness, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.rich.RS.aov1)
# other model
aoa.rich.RS.aov2 <- aov(aoa.meta.rh$Richness ~ Irrigation*Treatment*Date + 
               Error(rep2/(Irrigation*Treatment*Date)), data=aoa.meta.rh)
summary(aoa.rich.RS.aov2)
hist(aoa.meta.rh$Richness) # fine

# Model Fit
set.seed(13)
aoa.rich.rhizo.mod <- lmerTest::lmer(aoa.meta.rh$Richness ~ Irrigation*Treatment*Date +(1|Block)+(1|Block:Date), 
                                     data=aoa.meta.rh,contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"))
car::Anova(aoa.rich.rhizo.mod, test="F", type="III") 
shapiro.test(resid(aoa.rich.rhizo.mod)) # not normal
plot(simulateResiduals(aoa.rich.rhizo.mod)) # okay

# pairwise comparisons
# 1. between fertilization treatment:
aoa.emm.rich.rh <- aoa.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.rich.rhizo.mod)
# 2. between irrigation:
aoa.emm.rich.irri.rh <- aoa.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.rich.RS.aov2)
aoa.emm.rich.irri.rh
# pairwise
aoa.rich.emm.rh <- emmeans(aoa.rich.rhizo.mod2, ~ Treatment)
aoa.rich.emm.rh.pair <- pairs(aoa.rich.emm.rh)
aoa.rich.emm.rh.pair
###########################################################################
# 2. Response variable: Shannon
###########################################################################
# 2a. Analyses of Bulk Soil
aoa.meta.bulk.sum.sha <- aoa.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
aoa.bulk.sum.sha.plot <- ggboxplot(
  aoa.meta.bulk, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aoa.bulk.sum.sha.plot
# check assumption (outliers)
aoa.bulk.sha.out <- aoa.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
view(aoa.bulk.sha.out)
# Saphiro-Wilk for normality
aoa.bulk.sha.SW <- aoa.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
view(aoa.bulk.sha.SW)
ggqqplot(aoa.meta.bulk, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.bulk.sha.Lave <- aoa.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aoa.bulk.sha.aov <- anova_test(
  data = aoa.meta.bulk, type=2, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aoa.bulk.sha.aov)

# Three-Way Repeated-Measures ANOVA
aoa.bulk.sha.aov <- anova_test(
  data = aoa.meta.bulk, type=3,dv = Shannon, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.bulk.sha.aov)

# 4 Time Points:
# Three-Way Repeated-Measures ANOVA
aoa.bulk.sha.aov.4 <- anova_test(
  data = aoa.meta.bulk.4, type=3,dv = Shannon, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.bulk.sha.aov.4)
# Test Method 3 
aoa.bulk.sha.lmer4 <- lme4::lmer(aoa.meta.bulk.4$Shannon ~ Irrigation*Treatment*Date +(1|Block)+(1|Block:Date), 
                           data=aoa.meta.bulk.4,contrasts = list(Irrigation="contr.sum",Treatment="contr.sum",Date="contr.sum"))
car::Anova(aoa.bulk.sha.lmer4, test="F", type="III") 
shapiro.test(resid(aoa.bulk.sha.lmer4)) # 
plot(simulateResiduals(aoa.bulk.sha.lmer4)) # 

############################################################################################################
# Model Fit
set.seed(13)
#anova(aoa.sha.bulk.mod)

aoa.sha.bulk.mod2 <- lmerTest::lmer(aoa.meta.bulk$Shannon ~ Irrigation*Treatment*Date +(1|Date:Block), data=aoa.meta.bulk)
anova(aoa.sha.bulk.mod2)

shapiro.test(resid(aoa.sha.bulk.mod2)) # not normal
plot(simulateResiduals(aoa.sha.bulk.mod2)) # okay
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
aoa.emm.sha.bulk <- aoa.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Shannon ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.sha.bulk.mod)
# 2. between irrigation:
aoa.emm.sha.irri.bulk <- aoa.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Shannon ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.sha.bulk.mod2)
aoa.emm.sha.irri.bulk
# pairwise
aoa.sha.emm <- emmeans(aoa.sha.bulk.mod2, ~ Treatment)
aoa.sha.emm.pair <- pairs(aoa.sha.emm)
aoa.sha.emm.pair
#############################################################################################################

##################################################################
# 2b. Analyses of Rhizosphere
aoa.meta.rh.sum.sha <- aoa.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
aoa.rh.sum.sha.plot <- ggboxplot(
  aoa.meta.rh, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aoa.rh.sum.sha.plot
# check assumption (outliers)
aoa.rh.sha.out <- aoa.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
# Saphiro-Wilk for normality
aoa.rh.sha.SW <- aoa.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
ggqqplot(aoa.meta.rh, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.rh.sha.Lave <- aoa.min.meta.rh %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(aoa.meta.rh)
aoa.rh.sha.aov <- anova_test(
  data = aoa.meta.rh, type=3, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aoa.rh.sha.aov)

# Three-Way Repeated-Measures ANOVA
aoa.rh.sha.aov <- anova_test(
  data = aoa.meta.rh, type=3,dv = Shannon, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.rh.sha.aov)
############################################################################################################
# Model Fit
set.seed(13)
aoa.sha.rh.mod <- lmerTest::lmer(aoa.meta.rh$Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=aoa.meta.rh)
anova(aoa.sha.rh.mod, type = 2)

aoa.sha.rhizo.mod2 <- lmerTest::lmer(aoa.meta.rh$Shannon ~ Irrigation*Treatment*Date +(1|Block:Date), data=aoa.meta.rh)
anova(aoa.sha.rhizo.mod2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
aoa.emm.sha.rh <- aoa.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Shannon ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.sha.rh.mod)
# 2. between irrigation:
aoa.emm.sha.irri.rh <- aoa.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Shannon ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.sha.rh.mod)
#############################################################################################################

###########################################################################
# 3. Response variable: Inverse Simpson
###########################################################################
# 4a. Analyses of Bulk Soil
aoa.meta.bulk.invsum.simp <- aoa.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
aoa.bulk.sum.invsimp.plot <- ggboxplot(
  aoa.meta.bulk, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aoa.bulk.sum.invsimp.plot
# check assumption (outliers)
aoa.bulk.invsimp.out <- aoa.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
aoa.bulk.invsimp.SW <- aoa.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(aoa.meta.bulk, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.bulk.invsimp.Lave <- aoa.meta.bulk %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aoa.bulk.invsimp.aov <- anova_test(
  data = aoa.meta.bulk, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aoa.bulk.invsimp.aov)

# Three-Way Repeated-Measures ANOVA
aoa.bulk.simp.aov <- anova_test(
  data = aoa.meta.bulk, type=3,dv = InvSimpson, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.bulk.simp.aov)


# Three-Way Repeated-Measures ANOVA 4 tIME pOINTS
aoa.bulk.simp.aov.4 <- anova_test(
  data = aoa.meta.bulk.4, type=3,dv = InvSimpson, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.bulk.simp.aov.4)
############################################################################################################
# Model Fit
set.seed(13)
aoa.invsimp.bulk.mod <- lmerTest::lmer(aoa.meta.bulk$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=aoa.meta.bulk)
anova(aoa.invsimp.bulk.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
aoa.emm.invsimp.bulk <- aoa.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(InvSimpson ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.invsimp.bulk.mod)
# 2. between irrigation:
aoa.emm.invsimp.irri.bulk <- aoa.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(InvSimpson ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.invsimp.bulk.mod)
#############################################################################################################

##################################################################
# 3b. Analyses of Rhizosphere
aoa.meta.rh.sum.invsimp <- aoa.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
aoa.rh.sum.invsimp.plot <- ggboxplot(
  aoa.meta.rh, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aoa.rh.sum.invsimp.plot
# check assumption (outliers)
aoa.rh.invsimp.out <- aoa.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
aoa.rh.invsimp.SW <- aoa.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(aoa.meta.rh, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aoa.rh.invsimp.Lave <- aoa.meta.rh %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aoa.rh.invsimp.aov <- anova_test(
  data = aoa.min.meta.rh, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aoa.rh.invsimp.aov)

# Three-Way Repeated-Measures ANOVA
aoa.rh.simp.aov <- anova_test(
  data = aoa.meta.rh, type=2,dv = InvSimpson, wid = rep2,
  within = c(Irrigation, Treatment, Date))
get_anova_table(aoa.rh.simp.aov)
############################################################################################################
# Model Fit
set.seed(13)
aoa.invsimp.rh.mod <- lmerTest::lmer(aoa.meta.rh$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=aoa.meta.rh)
anova(aoa.invsimp.rh.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
aoa.emm.invsimp.rh <- aoa.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(InvSimpson ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.invsimp.rh.mod)
# 2. between irrigation:
aoa.emm.invsimp.irri.rh <- aoa.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(InvSimpson ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = aoa.invsimp.rh.mod)
#############################################################################################################













