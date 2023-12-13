################################################################################
# Microservices Project: Soil Properties Statistical Analysis
################################################################################
# Date : 23 August 2023
# Author : Ari Fina Bintarti

library(datarium)
library(rstatix)

# set working directory
setwd('/Users/arifinabintarti/Documents/France/microservices/')
setwd('D:/Fina/INRAE_Project/microservices/')
wd <- print(getwd())
# load the metadata
soilprop <- read.csv('Soil_SamplingPlotData_DOK2022_edited.csv')
soilprop$Date <- factor(soilprop$Date, levels = c("28/04/2022", "01/06/2022", "05/07/2022", "20/07/2022", "13/09/2022"),
                          labels = c("04-28-22", "06-01-22", "07-05-22", "07-20-22", "09-13-22"))
soilprop$Treatment <- factor(soilprop$Treatment, levels = c("D", "K", "M"),
                                    labels = c("Biodynamic (D)", "Conventional (K)", "Mineral fertilized (M)"))
###########################################################################

# 1. Response variable: GWC (Gravimetric water content in water content (g) / dry soil (g))
str(soilprop)
gwc.sum <- soilprop %>%
  group_by(Type, Treatment, Date) %>%
  get_summary_stats(GWC, type = "mean_sd")
gwc.sum.plot <- ggboxplot(
  soilprop, x = "Type", y = "GWC",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
gwc.sum.plot
# check assumption (outliers)
gwc.out <- soilprop %>%
  group_by(Type, Treatment, Date) %>%
  identify_outliers(GWC) # no extreme outliers
# Saphiro-Wilk for normality
gwc.SW <- soilprop %>%
  group_by(Type, Treatment, Date) %>%
  shapiro_test(GWC)
ggqqplot(soilprop, "GWC", ggtheme = theme_bw()) +
  facet_grid(Date ~ Treatment, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
gwc.Lave <- soilprop %>%
  group_by(Date) %>%
  levene_test(GWC ~ Type*Treatment)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
gwc.3aov <- anova_test(
  data = soilprop, type=3, dv = GWC, wid = PlotID,
  within = Date, between = c(Type, Treatment))
get_anova_table(gwc.3aov)
# Other models 
gwc.lme <- lme(GWC ~ Type*Treatment*Date, random=~1 | PlotID, method="ML", data=soilprop)
anova(gwc.lme, type = "sequential") # similar results as the three way mixed ANOVA above!!!
gwc.aov <- aov(soilprop$GWC ~ Type*Treatment*Date + Error(PlotID/Type*Treatment) , data=soilprop)
summary(gwc.aov) # # similar results as the three way mixed ANOVA above!!!

# Model Fit
set.seed(13)
gwc.lmer <- lmerTest::lmer(soilprop$GWC ~ Type*Treatment*Date +(1|PlotID), data=soilprop)
anova(gwc.lmer, type = 3)

# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
gwc.pwc.trt <- soilprop %>%
  group_by(Date, Type) %>%
  emmeans_test(GWC ~ Treatment, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = gwc.lmer)
# 2. between irrigation:
gwc.pwc.irr <- soilprop %>%
  group_by(Date, Treatment) %>%
  emmeans_test(GWC ~ Type, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = gwc.lmer)

# GWC: plotting the significance between irrigation within treatment and date

aob.sha.pwc.plot <- ggplot(aob.meta.df.sub, aes(x=Irrigation, y=Shannon)) +
  geom_boxplot(aes(fill = Treatment))+
  theme_bw() +
  labs(y="AOB Shannon")+
  labs(pattern="Irrigation")+
  scale_fill_viridis(discrete=T)+
  facet_grid(Type~ Date,scales="free_x")+
  theme(legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
aob.sha.pwc.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
xy.sha.bulk <- emm.sha.bulk %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
xy.sha.rh <- emm.sha.rh %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
# #combine two data frames and adding 'Type'
df.xy.sha.bulk <- as.data.frame(xy.sha.bulk)
df.xy.sha.rh <- as.data.frame(xy.sha.rh)
df.xy.sha.all <- rbind(df.xy.sha.bulk, df.xy.sha.rh) 
df.xy.sha.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
# plotting the pairwise comparisons among treatments (emmeans results)
aob.sha.pwc.plot2 <- aob.sha.pwc.plot + 
  stat_pvalue_manual(df.xy.sha.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.sha.pwc.plot2
# save figure
setwd('D:/Fina/INRAE_Project/microservices_fig/AOB')
ggsave("AOB_sha_boxplot.tiff",
       aob.sha.pwc.plot2, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)







#############################################################################################################






























