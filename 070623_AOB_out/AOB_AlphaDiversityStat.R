###########################################################################
#################### ALPHA DIVERSITY ANALYSIS: AOB ########################
###########################################################################

# Author: Ari Fina Bintarti
# Date: 03/07/2023

install.packages("datarium")
install.packages("rstatix")
library(datarium)
library(rstatix)
###########################################################################
# 1. Response variable: Richness
###########################################################################
# 1a. Analyses of Bulk Soil
aob.meta.bulk <- aob.meta.df.sub[1:119,]
str(aob.meta.bulk)
aob.meta.bulk.sum.rich <- aob.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
aob.bulk.sum.rich.plot <- ggboxplot(
  aob.meta.bulk, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.bulk.sum.rich.plot
# check assumption (outliers)
aob.bulk.rich.out <- aob.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
aob.bulk.rich.SW <- aob.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(aob.meta.bulk, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.bulk.rich.Lave <- aob.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
aob.bulk.rich.aov <- anova_test(
  data = aob.meta.bulk, type=2, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.bulk.rich.aov)
# Test Method 3 
model3 <- lme(Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, data=aob.meta.bulk, na.action = na.omit)
anova(model3) ### lme() IS SIMILAR to lme4 lmer() ###
#model1 <- aov(aob.meta.bulk$Richness ~ Irrigation*Treatment*Date + Error(PlotID/Irrigation*Treatment) , data=aob.meta.bulk)
#summary(model1)

############################################################################################################
# Model Fit
set.seed(13)
rich.bulk.mod <- lmerTest::lmer(aob.meta.bulk$Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.bulk, na.action = na.omit)
anova(rich.bulk.mod,  type="II")
library(lme4)
library(afex)
rich.bulk.mod2 <- lme4::lmer(aob.meta.bulk$Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.bulk, na.action = na.omit)
anova(rich.bulk.mod2)
anova(model1)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
emm.rich.bulk <- aob.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = rich.bulk.mod)
# 2. between irrigation:
emm.rich.irri.bulk <- aob.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = rich.bulk.mod)
#############################################################################################################

# Check other methods:
# pairwise comparisons
# library(emmeans)
emm.bulk <- emmeans(rich.bulk.mod, ~ Irrigation*Treatment*Date)
con <- contrast(emm.bulk, "pairwise", simple = "Irrigation")
# summary(con, adjust="fdr")
meta04 <- aob.meta.bulk[which(aob.meta.bulk$Date == "04-28-22"),]
rich.bulk.mod1 <- lmer(meta04$Richness ~ Irrigation*Treatment + (1|var2), data=meta04)
sum.mod = summary(rich.bulk.mod1)$coefficients
rownames(sum.mod)


# 1. between irrigation treatments
emm.irri.bulk <- emmeans(rich.bulk.mod1, pairwise~Irrigation|Treatment)
sum=summary(emm.irri.bulk)
sum.con=sum[["contrasts"]]
tmp = unlist(strsplit(as.character(sum.con$contrast)," - "))
tmp
sum.con[,"a"] <- tmp[seq(1,length(tmp),by=2)]
sum.con[,"b"] <- tmp[seq(2,length(tmp),by=2)]
v = sum.con[sum.con$a == "1a" | sum.con$b == "1a" ,]

con.irri.bulk <- contrast(emm.irri.bulk, "pairwise", adjust="fdr")
#df.con.irri.bulk <-  as.data.frame(con.irri.bulk)

# 2. among fertilization treatments
#emm.trt.bulk <- emmeans(rich.bulk.mod,~Treatment|Irrigation*Date, lmer.df = "satterthwaite")
#con.trt.bulk <- contrast(emm.trt.bulk, "pairwise", adjust="fdr")
#df.con.trt.bulk <-  as.data.frame(con.trt.bulk)

# adding significance codes
#df.con.trt.bulk$p.adj.signif <- ifelse(df.con.trt.bulk$p.value>0.05,"ns", 
                                       #ifelse(df.con.trt.bulk$p.value>=0.01&df.con.trt.bulk$p.value<=0.05,"*",
                                       #ifelse(df.con.trt.bulk$p.value>=0.001&df.con.trt.bulk$p.value<=0.01,"**",
                                       #ifelse(df.con.trt.bulk$p.value>=0&df.con.trt.bulk$p.value<=0.001,"***", "ns"))))
#df.con.trt.bulk.tid <- separate(data = df.con.trt.bulk, col = contrast, into = c("group1", "group2"), sep = "\\-")
#df.con.trt.bulk.tid$Type <- "Bulk Soil"

#df.con.irri.bulk$p.adj.signif <- ifelse(df.con.irri.bulk$p.value>0.05,"ns", 
                                        #ifelse(df.con.irri.bulk$p.value>=0.01&df.con.irri.bulk$p.value<=0.05,"*",
                                        #ifelse(df.con.irri.bulk$p.value>=0.001&df.con.irri.bulk$p.value<=0.01,"**",
                                        #ifelse(df.con.irri.bulk$p.value>=0&df.con.irri.bulk$p.value<=0.001,"***", "ns"))))
#df.con.irri.bulk.tid <- separate(data = df.con.irri.bulk, col = contrast, into = c("group1", "group2"), sep = "\\-")
#df.con.irri.bulk.tid$Type <- "Bulk Soil"

######################################################################
# 1b. Analyses of rhizosphere Soil
aob.meta.rh <- aob.meta.df.sub[120:191,]
str(aob.meta.rh)
aob.meta.rh.sum.rich <- aob.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Richness, type = "mean_sd")
aob.rh.sum.rich.plot <- ggboxplot(
  aob.meta.rh, x = "Irrigation", y = "Richness",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.rh.sum.rich.plot
# check assumption (outliers)
aob.rh.rich.out <- aob.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Richness) # no extreme outliers
# Saphiro-Wilk for normality
aob.rh.rich.SW <- aob.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Richness)
ggqqplot(aob.meta.rh, "Richness", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.rh.rich.Lave <- aob.meta.rh %>%
  group_by(Date) %>%
  levene_test(Richness ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aob.rh.rich.aov <- anova_test(
  data = aob.meta.rh, dv = Richness, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.rh.rich.aov)
# other model
# model1 <- aov(aob.meta.rh$Richness ~ Irrigation*Treatment*Date + Error(PlotID/Irrigation*Treatment), data=aob.meta.rh)
# summary(model1)
# Test Method 3
# model3 <- lme(Richness ~ Irrigation*Treatment*Date, random=~1 | PlotID, method="ML", data=aob.meta.rh)
# anova(model3)

# Model Fit
set.seed(13)
rich.rhizo.mod <- lmerTest::lmer(aob.meta.rh$Richness ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.rh)
#library(pbkrtest)
#anova(rich.rhizo.mod, type = 2, ddf= "Kenward-Roger") # similar with below
anova(rich.rhizo.mod, type = 2)
# pairwise comparisons
# 1. between fertilization treatment:
emm.rich.rh <- aob.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Richness ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = rich.rhizo.mod)
# 2. between irrigation:
emm.rich.irri.rh <- aob.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Richness ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = rich.rhizo.mod)
###########################################################################
# 2. Response variable: Shannon
###########################################################################
# 2a. Analyses of Bulk Soil
aob.meta.bulk.sum.sha <- aob.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
aob.bulk.sum.sha.plot <- ggboxplot(
  aob.meta.bulk, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.bulk.sum.sha.plot
# check assumption (outliers)
aob.bulk.sha.out <- aob.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
# Saphiro-Wilk for normality
aob.bulk.sha.SW <- aob.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
ggqqplot(aob.meta.bulk, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.bulk.sha.Lave <- aob.meta.bulk %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aob.bulk.sha.aov <- anova_test(
  data = aob.meta.bulk, type=2, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.bulk.sha.aov)
############################################################################################################
# Model Fit
set.seed(13)
sha.bulk.mod <- lmerTest::lmer(aob.meta.bulk$Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.bulk)
anova(sha.bulk.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
emm.sha.bulk <- aob.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Shannon ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = sha.bulk.mod)
# 2. between irrigation:
emm.sha.irri.bulk <- aob.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Shannon ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = sha.bulk.mod)
#############################################################################################################

##################################################################
# 2b. Analyses of Rhizosphere
aob.meta.rh.sum.sha <- aob.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(Shannon, type = "mean_sd")
aob.rh.sum.sha.plot <- ggboxplot(
  aob.meta.rh, x = "Irrigation", y = "Shannon",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.rh.sum.sha.plot
# check assumption (outliers)
aob.rh.sha.out <- aob.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(Shannon) # no extreme outliers
# Saphiro-Wilk for normality
aob.rh.sha.SW <- aob.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(Shannon)
ggqqplot(aob.meta.rh, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.rh.sha.Lave <- aob.min.meta.rh %>%
  group_by(Date) %>%
  levene_test(Shannon ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
str(aob.meta.rh)
aob.rh.sha.aov <- anova_test(
  data = aob.meta.rh, type=3, dv = Shannon, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.rh.sha.aov)
############################################################################################################
# Model Fit
set.seed(13)
sha.rh.mod <- lmerTest::lmer(aob.meta.rh$Shannon ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.rh)
anova(sha.rh.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
emm.sha.rh <- aob.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(Shannon ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = sha.rh.mod)
# 2. between irrigation:
emm.sha.irri.rh <- aob.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(Shannon ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = sha.rh.mod)
#############################################################################################################

###########################################################################
# 3. Response variable: Inverse Simpson
###########################################################################
# 4a. Analyses of Bulk Soil
aob.meta.bulk.invsum.simp <- aob.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
aob.bulk.sum.invsimp.plot <- ggboxplot(
  aob.meta.bulk, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.bulk.sum.invsimp.plot
# check assumption (outliers)
aob.bulk.invsimp.out <- aob.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
aob.bulk.invsimp.SW <- aob.meta.bulk %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(aob.meta.bulk, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.bulk.invsimp.Lave <- aob.meta.bulk %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aob.bulk.invsimp.aov <- anova_test(
  data = aob.meta.bulk, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.bulk.invsimp.aov)
############################################################################################################
# Model Fit
set.seed(13)
invsimp.bulk.mod <- lmerTest::lmer(aob.meta.bulk$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.bulk)
anova(invsimp.bulk.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
emm.invsimp.bulk <- aob.meta.bulk %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(InvSimpson ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = invsimp.bulk.mod)
# 2. between irrigation:
emm.invsimp.irri.bulk <- aob.meta.bulk %>%
  group_by(Date, Treatment) %>%
  emmeans_test(InvSimpson ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = invsimp.bulk.mod)
#############################################################################################################

##################################################################
# 3b. Analyses of Rhizosphere
aob.meta.rh.sum.invsimp <- aob.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  get_summary_stats(InvSimpson, type = "mean_sd")
aob.rh.sum.invsimp.plot <- ggboxplot(
  aob.meta.rh, x = "Irrigation", y = "InvSimpson",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
aob.rh.sum.invsimp.plot
# check assumption (outliers)
aob.rh.invsimp.out <- aob.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  identify_outliers(InvSimpson) # no extreme outliers
# Saphiro-Wilk for normality
aob.rh.invsimp.SW <- aob.meta.rh %>%
  group_by(Irrigation, Treatment, Date) %>%
  shapiro_test(InvSimpson)
ggqqplot(aob.meta.rh, "InvSimpson", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
aob.rh.invsimp.Lave <- aob.meta.rh %>%
  group_by(Date) %>%
  levene_test(InvSimpson ~ Irrigation*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
aob.rh.invsimp.aov <- anova_test(
  data = aob.min.meta.rh, dv = InvSimpson, wid = PlotID,
  within = Date, between = c(Irrigation, Treatment))
get_anova_table(aob.rh.invsimp.aov)
############################################################################################################
# Model Fit
set.seed(13)
invsimp.rh.mod <- lmerTest::lmer(aob.meta.rh$InvSimpson ~ Irrigation*Treatment*Date +(1|PlotID), data=aob.meta.rh)
anova(invsimp.rh.mod, type = 2)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
emm.invsimp.rh <- aob.meta.rh %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(InvSimpson ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = invsimp.rh.mod)
# 2. between irrigation:
emm.invsimp.irri.rh <- aob.meta.rh %>%
  group_by(Date, Treatment) %>%
  emmeans_test(InvSimpson ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = invsimp.rh.mod)
#############################################################################################################













