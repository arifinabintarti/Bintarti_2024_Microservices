################################################################################
# Microservices Project: Soil Properties Statistical Analysis
################################################################################
# Date : 23 August 2023
# Author : Ari Fina Bintarti

library(datarium)
library(rstatix)

# set working directory
setwd('/Users/arifinabintarti/Documents/France/microservices/')
#setwd('D:/Fina/INRAE_Project/microservices/')
wd <- print(getwd())
# load the metadata
soilprop <- read.csv('Soil_SamplingPlotData_DOK2022_edited.csv')
soilprop$Date <- factor(soilprop$Date, levels = c("28/04/2022", "1/6/22", "5/7/22", "20/07/2022", "13/09/2022"),
                          labels = c("Apr", "Jun", "Jul5", "Jul20", "Sep"))
soilprop$Fertilization <- factor(soilprop$Fertilization, levels = c("D", "K", "M"),
                                    labels = c("BIODYN", "CONFYM", "CONMIN"))
str(soilprop)
soilprop$SampleID<-factor(soilprop$SampleID)
soilprop$PlotID<-factor(soilprop$PlotID)
soilprop$Irrigation<-factor(soilprop$Irrigation)
soilprop$Block<-factor(soilprop$Block)
soilprop$x<-factor(soilprop$x)
soilprop$rep<-factor(soilprop$rep)
soilprop$rep2<-factor(soilprop$rep2)
soilprop$var3<-factor(soilprop$var3)
soilprop$Depth<-factor(soilprop$Depth)
soilprop$Mg_mgkg<-as.numeric(soilprop$Mg_mgkg)

###########################################################################

# 1. Response variable: GWC (Gravimetric water content in water content (g) / dry soil (g))
str(soilprop)
gwc.sum.ave <- soilprop %>%
  group_by(Irrigation, Fertilization) %>%
  get_summary_stats(GWC, type = "mean_sd")
View(gwc.sum.ave)
gwc.sum <- soilprop %>%
  group_by(Date, Irrigation, Fertilization) %>%
  get_summary_stats(GWC, type = "mean_sd")
View(gwc.sum)
gwc.sum.plot <- ggboxplot(
  soilprop, x = "Type", y = "GWC",
  color = "Treatment", palette = "jco",
  facet.by =  "Date")
gwc.sum.plot
# check assumption (outliers)
gwc.out <- soilprop %>%
  group_by(Irrigation, Fertilization, Date) %>%
  identify_outliers(GWC) # no extreme outliers
View(gwc.out)
# Saphiro-Wilk for normality
gwc.SW <- soilprop %>%
  group_by(Irrigation, Fertilization, Date) %>%
  shapiro_test(GWC)
View(gwc.SW)
ggqqplot(soilprop, "GWC", ggtheme = theme_bw()) +
  facet_grid(Date ~ Fertilization, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
gwc.Lave <- soilprop %>%
  group_by(Date) %>%
  levene_test(GWC ~ Irrigation*Fertilization)
View(gwc.Lave)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Mixed (Split-Plot) ANOVA 
set.seed(13)
gwc.3aov <- anova_test(
  data = soilprop, type=3, dv = GWC, wid = PlotID,
  within = Date, between = c(Type, Treatment))
get_anova_table(gwc.3aov)
# Three-Way Repeated-Measures ANOVA
gwc.3rmaov <- anova_test(
  data = soilprop, type=3,dv = GWC, wid = rep2,
  within = c(Irrigation, Fertilization, Date))
get_anova_table(gwc.3rmaov)
soilprop$sqrtGWC <- sqrt(soilprop$GWC)
#similarly
gwc.3rmaov2 <- aov(GWC ~ Irrigation*Fertilization*Date + Error(rep2/(Irrigation*Fertilization*Date)), data=soilprop)
summary(gwc.3rmaov2)
# Other models 
gwc.lme <- lme(GWC ~ Irrigation*Fertilization*Date, random=~1 | PlotID, data=soilprop)
anova(gwc.lme, type = "sequential") # similar results as the three way mixed ANOVA above!!!
gwc.aov <- aov(soilprop$GWC ~ Irrigation*Fertilization*Date + Error(Block/Irrigation*Fertilization*Date) , data=soilprop)
summary(gwc.aov) # # similar results as the three way mixed ANOVA above!!!
# Model Fit
set.seed(13)
gwc.lmer <- lmerTest::lmer(soilprop$GWC ~ Irrigation*Fertilization*Date +(1|Block)+(1|Block:Date), 
                           data=soilprop,contrasts = list(Irrigation="contr.sum",Fertilization="contr.sum",Date="contr.sum"))
car::Anova(gwc.lmer, test="F", type="III") 
# test assumptions:
hist(soilprop$GWC) #not bad
shapiro.test(resid(gwc.lmer)) # not normal
plot(simulateResiduals(gwc.lmer)) # okay
plot(gwc.lmer, which = 3)
# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between fertilization treatment:
gwc.pwc.trt <- soilprop %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(GWC ~ Fertilization, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = gwc.lmer)
View(gwc.pwc.trt)
# 2. between irrigation:
gwc.pwc.irr <- soilprop %>%
  group_by(Date, Fertilization) %>%
  emmeans_test(GWC ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = gwc.3rmaov2)
View(gwc.pwc.irr)
# Test only 4 time points:
soilprop.4 <- soilprop[1:96,]
str(soilprop.4)
soilprop.4$SampleID<-factor(soilprop.4$SampleID)
soilprop.4$PlotID<-factor(soilprop.4$PlotID)
soilprop.4$Irrigation<-factor(soilprop.4$Irrigation)
soilprop.4$Block<-factor(soilprop.4$Block)
soilprop.4$x<-factor(soilprop.4$x)
soilprop.4$rep<-factor(soilprop.4$rep)
soilprop.4$rep2<-factor(soilprop.4$rep2)
soilprop.4$var3<-factor(soilprop.4$var3)
soilprop.4$Date<-factor(soilprop.4$Date)
# Three-Way Repeated-Measures ANOVA
gwc.3rmaov.4 <- anova_test(
  data = soilprop.4, type=3,dv = GWC, wid = rep2,
  within = c(Irrigation, Fertilization, Date))
get_anova_table(gwc.3rmaov.4)
# Model Fit
set.seed(13)
gwc.lmer.4 <- lmerTest::lmer(soilprop.4$GWC ~ Irrigation*Fertilization*Date +(1|Block)+(1|Block:Date), 
                           data=soilprop.4,contrasts = list(Irrigation="contr.sum",Fertilization="contr.sum",Date="contr.sum"))
anova(gwc.lmer.4, type = 3)
car::Anova(gwc.lmer.4, test="F", type="III") 

# GWC: plotting the significance between irrigation within treatment and date
gwc.pwc.plot <- ggplot(soilprop, aes(x=Date, y=GWC)) +
  geom_boxplot(aes(fill = Fertilization))+
  theme_bw() +
  labs(y="AOB Shannon")+
  labs(pattern="Irrigation")+
  scale_fill_viridis(discrete=T)+
  #facet_grid(Type~ Date,scales="free_x")+
  theme(legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
gwc.pwc.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
xy.gwc <- emm.sha.bulk %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
xy.sha.rh <- emm.sha.rh %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
# #combine two data frames and adding 'Type'
df.xy.sha.bulk <- as.data.frame(xy.sha.bulk)
df.xy.sha.rh <- as.data.frame(xy.sha.rh)
df.xy.sha.all <- rbind(df.xy.sha.bulk, df.xy.sha.rh) 
df.xy.sha.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
#############################################################################################################

# 1. Response variable: NH4 (Ammonium content)
# check assumption (outliers)
nh4.out <- soilprop %>%
  group_by(Irrigation, Fertilization, Date) %>%
  identify_outliers(log10NH4) # two extreme outliers
View(nh4.out)
# Saphiro-Wilk for normality
nh4.SW <- soilprop %>%
  group_by(Irrigation, Fertilization, Date) %>%
  shapiro_test(log10NH4) #  good
View(nh4.SW)
ggqqplot(soilprop, "log10NH4", ggtheme = theme_bw()) +
  facet_grid(Date ~ Fertilization, labeller = "label_both") # good
nh4.Lave <- soilprop %>%
  group_by(Date) %>%
  levene_test(log10NH4 ~ Irrigation*Fertilization) # good
View(nh4.Lave)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Repeated-Measures ANOVA
nh4.3rmaov <- anova_test(
  data = soilprop, type=3,dv = log10NH4, wid = rep2,
  within = c(Irrigation, Fertilization, Date))
get_anova_table(nh4.3rmaov)
#similarly
nh4.3rmaov2 <- aov(log10NH4 ~ Irrigation*Fertilization*Date + Error(rep2/(Irrigation*Fertilization*Date)), data=soilprop)
summary(nh4.3rmaov2)
# Model Fit
set.seed(13)
nh4.lmer <- lmerTest::lmer(soilprop$log10NH4 ~ Irrigation*Fertilization*Date +(1|Block)+(1|Block:Date), 
                           data=soilprop,contrasts = list(Irrigation="contr.sum",Fertilization="contr.sum",Date="contr.sum"))
anova(nh4.lmer, type = 3)
car::Anova(nh4.lmer, test="F", type="III") 
# test assumptions:
hist(soilprop$NH4) # highly right-skewed
shapiro.test(resid(nh4.lmer)) # normal
plot(simulateResiduals(nh4.lmer)) # okay
## MUST transform the data
## adding constant val of 1 and log
## log10
soilprop$log10NH4 <- log10(soilprop$NH4+1)
hist(soilprop$log10NH4) # better
## log
soilprop$logNH4 <- log(soilprop$NH4+1)
hist(soilprop$logNH4) # better
# Model Fit
set.seed(13)
nh4.log.lmer.4 <- lmerTest::lmer(soilprop.4$log10NH4 ~ Irrigation*Fertilization*Date +(1|Block)+(1|Block:Date), 
                           data=soilprop.4,contrasts = list(Irrigation="contr.sum",Fertilization="contr.sum",Date="contr.sum"))
anova(nh4.log.lmer.4, type = 3)
car::Anova(nh4.log.lmer.4, test="F", type="III") 
shapiro.test(resid(nh4.log.lmer)) # normal
plot(simulateResiduals(nh4.log.lmer)) # okay
# Three-Way Repeated-Measures ANOVA: 4 Time Points
soilprop.4$log10NH4 <- log10(soilprop.4$NH4+1) # add x+a =1
nh4.log.3rmaov.4 <- anova_test(
  data = soilprop.4, type=3,dv = log10NH4, wid = rep2,
  within = c(Irrigation, Fertilization, Date))
get_anova_table(nh4.log.3rmaov.4)
#similarly
nh4.3rmaov.4 <- aov(log10NH4 ~ Irrigation*Fertilization*Date + Error(rep2/(Irrigation*Fertilization*Date)), data=soilprop.4)
summary(nh4.3rmaov.4)
# Two-way ANOVA at each fertilization level
nh4.log.2way <- soilprop.4 %>%
  group_by(Fertilization) %>%
  anova_test(dv = log10NH4, wid = rep2, within = c(Irrigation, Date))
get_anova_table(nh4.log.2way)
# Effect of drought at each fert X date
drought.effect <- soilprop.4 %>%
  group_by(Date, Fertilization) %>%
  anova_test(dv = log10NH4, wid = rep2, within = Irrigation)
drought.effect
# pairwise comparisons
# 1. between fertilization treatment:
nh4.crop.emm <- soilprop.4 %>%
  group_by(Date, Irrigation) %>%
  emmeans_test(log10NH4 ~ Fertilization, 
               p.adjust.method = "BH", 
               conf.level = 0.95)
View(nh4.crop.emm)
# 2. between irrigation:
nh4.irri.emm <- soilprop %>%
  group_by(Date, Fertilization) %>%
  emmeans_test(log10NH4 ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model=nh4.3rmaov2)
View(nh4.irri.emm) 
# NH4 Summary
nh4.sum <- soilprop %>%
  group_by(Irrigation, Fertilization) %>%
  get_summary_stats(NH4, type = "mean_sd")
View(nh4.sum)
nh4.sum.date <- soilprop %>%
  group_by(Irrigation, Fertilization, Date) %>%
  get_summary_stats(NH4, type = "mean_sd")
View(nh4.sum.date)
nh4.sum.plot <- ggboxplot(
  soilprop, x = "Irrigation", y = "NH4",
  color = "Fertilization", palette = "jco",
  facet.by =  "Date")
nh4.sum.plot


#############################################################################################################
# 1. Response variable: NO3 (Ammonium content)
# check assumption (outliers)
no3.out <- soilprop %>%
  group_by(Irrigation, Fertilization, Date) %>%
  identify_outliers(NO3) # no extreme outliers
View(no3.out)
# Saphiro-Wilk for normality
no3.SW <- soilprop %>%
  group_by(Irrigation, Fertilization, Date) %>%
  shapiro_test(NO3) # okay
View(no3.SW)
ggqqplot(soilprop, "NO3", ggtheme = theme_bw()) +
  facet_grid(Date ~ Fertilization, labeller = "label_both") #not good
no3.Lave <- soilprop %>%
  group_by(Date) %>%
  levene_test(logNO3 ~ Irrigation*Fertilization) #violated
View(no3.Lave)
#If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# Three-Way Repeated-Measures ANOVA
no3.3rmaov <- anova_test(
  data = soilprop, type=3,dv = logNO3, wid = rep2,
  within = c(Irrigation, Fertilization, Date))
get_anova_table(no3.3rmaov)
#similarly
no3.3rmaov2 <- aov(logNO3 ~ Irrigation*Fertilization*Date + Error(rep2/(Irrigation*Fertilization*Date)), data=soilprop)
summary(no3.3rmaov2)
# Model Fit
set.seed(13)
no3.lmer <- lmerTest::lmer(soilprop$sqrtNO3 ~ Irrigation*Fertilization*Date +(1|Block)+(1|Block:Date), 
                           data=soilprop,contrasts = list(Irrigation="contr.sum",Fertilization="contr.sum",Date="contr.sum"))
anova(no3.lmer, type = 3)
car::Anova(no3.lmer, test="F", type="III") 
# test assumptions:
hist(soilprop$NO3) # highly right-skewed
shapiro.test(resid(no3.lmer)) # not normal
plot(simulateResiduals(no3.lmer)) # not okay
## MUST transform the data
soilprop$logNO3 <- log(soilprop$NO3)
hist(soilprop$logNO3) # better
hist(soilprop$pH)
# Three-Way Repeated-Measures ANOVA: 4 Time Points
soilprop.4$logNO3 <- log(soilprop.4$NO3)
hist(soilprop.4$logNO3)
no3.3rmaov.4 <- anova_test(
  data = soilprop.4, type=3,dv = logNO3, wid = rep2,
  within = c(Irrigation, Fertilization, Date))
get_anova_table(no3.3rmaov.4)
# model with lmer
no3.lmer.4 <- lmerTest::lmer(soilprop.4$logNO3 ~ Irrigation*Fertilization*Date +(1|Block)+(1|Block:Date), 
                           data=soilprop.4,contrasts = list(Irrigation="contr.sum",Fertilization="contr.sum",Date="contr.sum"))
car::Anova(no3.lmer.4, test="F", type="III") 
# Two-way ANOVA at each fertilization level
no3.log.2way <- soilprop.4 %>%
  group_by(Date) %>%
  anova_test(dv = logNO3, wid = rep2, within = c(Irrigation, Fertilization))
get_anova_table(no3.log.2way)
# Effect of drought at each fert X date
no3.drought.effect <- soilprop.4 %>%
  group_by(Date, Fertilization) %>%
  anova_test(dv = logNO3, wid = rep2, within = Irrigation)
no3.drought.effect
# pairwise comparisons
# 1. between fertilization treatment:
no3.crop.emm <- soilprop %>%
  group_by(Irrigation) %>%
  emmeans_test(logNO3 ~ Fertilization, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model=no3.3rmaov2)
View(no3.crop.emm)
# 2. between irrigation:
no3.irri.emm <- soilprop %>%
  group_by(Date, Fertilization) %>%
  emmeans_test(logNO3 ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model=no3.3rmaov2)
View(no3.irri.emm) 
# NO3 Summary
no3.sum.ave <- soilprop %>%
  group_by(Irrigation, Fertilization) %>%
  get_summary_stats(NO3, type = "mean_sd")
View(no3.sum.ave)
# NO3Summary per date
no3.sum<- soilprop %>%
  group_by(Irrigation, Fertilization, Date) %>%
  get_summary_stats(NO3, type = "mean_sd")
View(no3.sum)

no3.sum.plot <- ggboxplot(
  soilprop, x = "Irrigation", y = "NO3",
  color = "Fertilization", palette = "jco",
  facet.by =  "Date")
no3.sum.plot

###############################################################################################################################

# PLOTS

# GWC
# adding ANOVA ressults of 4 samples:

soilprop$x <- factor(soilprop$x, levels = c("cont.D", "rain.D", "cont.K","rain.K","cont.M","rain.M"))
                             
stat_text.gwc <- data.frame(Date = 0.5, GWC = 0.05,Fertilization="BIODYN", label="D ***\nC *\nT ***\nD x T ***")

gwc.plot <- ggplot(soilprop , aes(x=Date, y=GWC)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control','Biodyn-drought', 'Confym-control','Confym-drought',
                             'Conmin-control', 'Conmin-drought'))+
  #ylab(bquote(bold('Gravimetric Water Content'~(g~g^-1~dry~soil))))+
  ylab(bquote('Gravimetric Water Content'~(g~g^-1~dry~soil)))+
  #labs(subtitle = "")+
  facet_wrap(~ Fertilization)+
  #scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        #strip.text = element_blank(),
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_text(size=18),
        plot.subtitle = element_text(size=20, face = "bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 geom_vline(xintercept = 3.4, linetype="dashed", colour="darkgrey") +
 annotate(geom = "text", x = 3.6, y = 0, hjust = 0, size = 4, label = "Rewetting", color = "gray25")+
 geom_label(data = stat_text.gwc,label=stat_text.gwc$label,hjust=0, colour="black", size=4, fontface="bold")
 gwc.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
gwc.pwc.irr <- soilprop %>%
  group_by(Date, Fertilization) %>%
  emmeans_test(GWC ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = gwc.3rmaov2)
View(gwc.pwc.irr)
# add x y position
gwc.pwc.xy <- gwc.pwc.irr  %>% 
  add_xy_position(x = "Date", dodge = 0.8) 
# plotting the pairwise comparisons among treatments (emmeans results)
 gwc.plot2 <-  gwc.plot + 
  stat_pvalue_manual(gwc.pwc.xy,x = "Date", y.position = 0.4,
                     #step.increase = 1,
                     #label = "p = {scales::pvalue(p.adj)}",size=3, 
                     label = "p.adj.signif",size=5,
                     #bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     #bracket.shorten = 1, color = "black",
                     tip.length = 0.01, hide.ns = F)
 gwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Fig.3dpi300.tiff",
       gwc.plot2, device = "tiff",bg="white",
       width = 11.5, height = 6.5, 
       units= "in", dpi = 300, compression="lzw")
ggsave("Fig.3.2dpi300.tiff",
       gwc.plot2, device = "tiff",bg="white",
       width = 11.5, height = 8, 
       #width = 11.5, height = 6.5,
       units= "in", dpi = 300, compression="lzw")



# NH4
# adding ANOVA results of 4 samples:

stat_text.NH4 <- data.frame(Date = 0.5, NH4 = 20,Fertilization="BIODYN", label="D **\nC ***\nT *\nD x C ***\nD x T ***\nC x T *\nD x C x T **")

NH4.plot <- ggplot(soilprop , aes(x=Date, y=NH4)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control','Biodyn-drought', 'Confym-control','Confym-drought',
                             'Conmin-control', 'Conmin-drought'))+
  ylab('NH<sub>4<sup>+ </sub>content<br>(mg Kg<sup>-1</sup>soil dry matter)')+
  #ylab(bquote(NH4^{"+"}~(mg~Kg^-1~soil~dry~matter)))+
  labs(title = "A")+
  facet_wrap(~ Fertilization)+
  #scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        #strip.text = element_blank(),
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_markdown(size=18),
        plot.title = element_text(size=20, face = "bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 geom_vline(xintercept = 3.4, linetype="dashed", colour="darkgrey") +
 #annotate(geom = "text", x = 3.6, y = -1, hjust = 0, size = 4, label = "Rewetting", color = "gray25")+
 geom_label(data = stat_text.NH4,label=stat_text.NH4$label,hjust=0, colour="black", size=4, fontface="bold")
NH4.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
nh4.irri.emm <- soilprop %>%
  group_by(Date, Fertilization) %>%
  emmeans_test(log10NH4 ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model=nh4.3rmaov2)
# add x y position
nh4.pwc.xy <- nh4.irri.emm  %>% 
  add_xy_position(x = "Date", dodge = 0.8) 
# plotting the pairwise comparisons among treatments (emmeans results)
NH4.plot2 <-  NH4.plot + 
  stat_pvalue_manual(nh4.pwc.xy,x = "Date", y.position = 33,
                     #step.increase = 1,
                     #label = "p = {scales::pvalue(p.adj)}",size=3, 
                     label = "p.adj.signif",size=5,
                     #bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     #bracket.shorten = 1, color = "black",
                     tip.length = 0.01, hide.ns = F)
NH4.plot2

# NO3
# adding ANOVA results of 4 samples:

stat_text.NO3 <- data.frame(Date = 0.5, NO3 = 40,Fertilization="BIODYN", label="D ***\nT ***\nD x C ***\nD x T ***\nC x T *")

NO3.plot <- ggplot(soilprop , aes(x=Date, y=NO3)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control','Biodyn-drought', 'Confym-control','Confym-drought',
                             'Conmin-control', 'Conmin-drought'))+
  ylab('NO<sub>3<sup>- </sub>content<br>(mg Kg<sup>-1</sup>soil dry matter)')+
  #ylab(bquote(NO3^{"-"}~(mg~Kg^-1~soil~dry~matter)))+
  labs(title = "B")+
  facet_wrap(~ Fertilization)+
  #scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.text = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_markdown(size=18),
        plot.title = element_text(size=20, face = "bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 geom_vline(xintercept = 3.4, linetype="dashed", colour="darkgrey") +
 #annotate(geom = "text", x = 3.6, y = -1, hjust = 0, size = 4, label = "Rewetting", color = "gray25")+
 geom_label(data = stat_text.NO3,label=stat_text.NO3$label,hjust=0, colour="black", size=4, fontface="bold")
NO3.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
no3.irri.emm <- soilprop %>%
  group_by(Date, Fertilization) %>%
  emmeans_test(logNO3 ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model=no3.3rmaov2)
# add x y position
no3.pwc.xy <- no3.irri.emm  %>% 
  add_xy_position(x = "Date", dodge = 0.8) 
# plotting the pairwise comparisons among treatments (emmeans results)
NO3.plot2 <-  NO3.plot + 
  stat_pvalue_manual(no3.pwc.xy,x = "Date", y.position = 65,
                     #step.increase = 1,
                     #label = "p = {scales::pvalue(p.adj)}",size=3, 
                     label = "p.adj.signif",size=5,
                     #bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     #bracket.shorten = 1, color = "black",
                     tip.length = 0.01, hide.ns = F)
NO3.plot2
 

N.pools <- (NH4.plot2 / NO3.plot2 / n2o.plot2)+
 plot_layout(guides = "collect") & 
 theme(legend.position = 'bottom',legend.title = element_blank())
N.pools

setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Fig.4.dpi300.tiff",
       N.pools, device = "tiff",bg="white",
       width = 11.5, height = 9.5, 
       units= "in", dpi = 300, compression="lzw")
###################################################################################################################################################################
####################################################################################################################################################################
# set working directory
setwd('/Users/arifinabintarti/Documents/France/microservices/')
# setwd('D:/Fina/INRAE_Project/microservices/')
wd <- print(getwd())
# load the whole N2O data
N2O.dat <- as.data.frame(read.csv('N20_data_ed.csv'))
view(N2O.dat)
# subset only data with flux value (1 value per sample)
N2O.dat.sub <- N2O.dat[N2O.dat$TIME_POINT == '1',]
view(N2O.dat.sub)
# make a unique time
N2O.dat.sub$Date <- factor(N2O.dat.sub$DSAMP, levels = c("4/5/22", "4/12/22", "4/19/22", "4/26/22",
                                                  "5/3/22","5/10/22", "5/17/22","5/24/22","5/31/22",
                                                  "6/10/22","6/14/22","6/21/22","6/28/22",
                                                  "7/5/22","7/15/22", "7/20/22"),
                          labels = c("Apr","Apr","Apr","Apr",
                                     "Jun","Jun","Jun","Jun","Jun",
                                     "Jul5","Jul5","Jul5","Jul5",
                                     "Jul5","Jul20","Jul20"))
# calculate mean
N2O.dat.sub.mean <- N2O.dat.sub %>%
  group_by(PLOT, fertilization, irrigation, Date) %>%
  summarise(Mean.N2O = mean(flux_nmol_m2_per_s, na.rm = TRUE)) %>%
  arrange(Date, PLOT)
N2O.dat.sub.mean <- N2O.dat.sub.mean %>% 
  rename(PlotID = PLOT)
view(N2O.dat.sub.mean)
N2O.dat.sub.mean <- as.data.frame(N2O.dat.sub.mean)
str(N2O.dat.sub.mean)
N2O.dat.sub.mean$PlotID<- factor(N2O.dat.sub.mean$PlotID)
N2O.dat.sub.mean$fertilization <- factor(N2O.dat.sub.mean$fertilization)
N2O.dat.sub.mean$irrigation <- factor(N2O.dat.sub.mean$irrigation)
# test data distribution
plotNormalHistogram(N2O.dat.sub.mean$Mean.N2O) # right skewed
plotNormalHistogram(log(N2O.dat.sub.mean$Mean.N2O)) #NaNs produced, right-skewed
plotNormalHistogram(log(N2O.dat.sub.mean$Mean.N2O+2)) # right skewed
# same with the log
# Cube root transformation (better)
T_cub = sign(N2O.dat.sub.mean$Mean.N2O) * abs(N2O.dat.sub.mean$Mean.N2O)^(1/3)
plotNormalHistogram(T_cub) # looks okay
N2O.dat.sub.mean$N2O.cubrt <- T_cub
skewness(soilprop.n2o$N2O.cubrt)
# tukey'ss ladder of power principle (not suitable)
T_tuk = transformTukey(N2O.dat.sub.mean$Mean,
                     plotit=FALSE)  #lambda of –7.8
plotNormalHistogram(T_tuk)
qqnorm(T_cub,
       ylab="Sample Quantiles for N20")
qqline(T_cub,
       col="red")
# box-cox (not suitable)
library(MASS)
Box = boxcox(Mean ~ irrigation*fertilization*Date,
             data = N2O.dat.sub.mean,
             lambda = seq(-6,6,0.1))
Cox = data.frame(Box$x, Box$y)
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
Cox2[1,]
lambda = Cox2[1, "Box.x"]
N20.mean.box = (N2O.dat.sub.mean$Mean ^ lambda - 1)/lambda 
plotNormalHistogram(N20.mean.box)
qqnorm(N20.mean.box,
       ylab="Sample Quantiles for N20")
# subset soilproperties data (remove the month of "Sep")
soilprop.nosep <- soilprop[1:96,]
view(soilprop.nosep)
# adding N2O data into soil properties data
soilprop.n2o <- soilprop.nosep
soilprop.n2o$mean.N2Oflux <- N2O.dat.sub.mean$Mean.N2O
soilprop.n2o$N2Oflux.cbrt <- N2O.dat.sub.mean$N2O.cubrt
str(soilprop.n2o)
soilprop.n2o$SampleID <- factor(soilprop.n2o$SampleID)
soilprop.n2o$PlotID <- factor(soilprop.n2o$PlotID)
soilprop.n2o$Block <- factor(soilprop.n2o$Block)
soilprop.n2o$Irrigation <- factor(soilprop.n2o$Irrigation)
soilprop.n2o$Fertilization <- factor(soilprop.n2o$Fertilization)
soilprop.n2o$Date <- factor(soilprop.n2o$Date)
soilprop.n2o$rep2 <- factor(soilprop.n2o$rep2)
soilprop.n2o$var3 <- factor(soilprop.n2o$var3)
soilprop.n2o$x <- factor(soilprop.n2o$x)
setwd('/Users/arifinabintarti/Documents/France/microservices/')
#write.csv(soilprop.n2o, file = "soilprop.n2o.csv")

# combine 2 dataframes
view(soilprop)
view(soilprop.n2o)
soilprop <- rownames_to_column(soilprop, var = "SampleID")
soilprop.n2o <- rownames_to_column(soilprop.n2o, var = "SampleID")
just.n2o <- soilprop.n2o[,c(1,31)]
head(just.n2o)
# left_join
soilprop.comp <- left_join(soilprop,just.n2o, by="SampleID")
view(soilprop.comp)



###########################################################
# statistical test
# Repeated measures ANOVA is fairly robust insofar as P-values are scientifically similar for raw data (with no outlier removal) and cube roots.
# Three-Way Repeated-Measures ANOVA with transformed data
n2o.3rmaov.cbrt <- anova_test(
  data = soilprop.n2o, type=3,dv = N2Oflux.cbrt, wid = rep2,
  within = c(Irrigation, Fertilization, Date))
get_anova_table(n2o.3rmaov.cbrt)
#similarly
n2o.3rmaov2.cbrt2 <- aov(N2Oflux.cbrt ~ Irrigation*Fertilization*Date + Error(rep2/(Irrigation*Fertilization*Date)), data=soilprop.n2o)
summary(n2o.3rmaov2.cbrt2)
# Model Fit
set.seed(13)
N2O.lmer.cbrt <- lmerTest::lmer(N2Oflux.cbrt ~ Irrigation*Fertilization*Date + (1|Block)+(1|Block:Date), 
                           data=soilprop.n2o, contrasts = list(Irrigation="contr.sum",Fertilization="contr.sum",Date="contr.sum"))
car::Anova(N2O.lmer.cbrt, test="F", type="III") # singularity warning

#N2O.dat.ed$x <- factor(N2O.dat.ed$x, levels = c("cont.D", "rain.D", "cont.K","rain.K","cont.M","rain.M"))
# check assumption (outliers)
N2O.out <- soilprop.n2o %>%
  group_by(Irrigation, Fertilization, Date) %>%
  identify_outliers(N2Oflux.cbrt) # no extreme outliers
View(N2O.out)
# Saphiro-Wilk for normality
N2O.SW <- soilprop.n2o %>%
  group_by(Irrigation, Fertilization, Date) %>%
  shapiro_test(N2Oflux.cbrt)
View(N2O.SW)
ggqqplot(soilprop.n2o, "N2Oflux.cbrt", ggtheme = theme_bw()) 
  #facet_grid(~ Fertilization*Irrigation, labeller = "label_both") #All the points fall approximately along the reference line, for each cell. So we can assume normality of the data
# Lavene test
N2O.Lave <- soilprop.n2o %>%
  group_by(Date) %>%
  levene_test(mean.N2Oflux.cbrt ~ Irrigation*Fertilization)
View(N2O.Lave)
# glmmtmb
#N2O.lmer <- glmer(mean.N2Oflux ~ Irrigation*Fertilization*Date + (1|Block)+(1|Block:Date), 
                  #data= soilprop.n2o, family = "gaussian",
                  #control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)),
                  #contrasts = list(Irrigation="contr.sum",Fertilization="contr.sum",Date="contr.sum"))
#car::Anova(N2O.lmer)
# test assumptions:
hist(soilprop.n2o$mean.N2Oflux)
hist(log(soilprop.n2o$mean.N2Oflux+1))
hist(residuals(N2O.lmer.cbrt)) # looks good
plot(simulateResiduals(N2O.lmer))
plot(simulateResiduals(N2O.lmer.cbrt)) # much better
shapiro.test(resid(N2O.lmer.cbrt)) # not good

# Fit pairwise comparisons
# Performs pairwise comparisons between groups using the estimated marginal means. Pipe-friendly wrapper around the functions emmeans() + contrast() from the emmeans package,
# 1. between irrigation:
# using aov mod
n2o.pwc.irr <- soilprop.n2o %>%
  group_by(Date, Fertilization) %>%
  emmeans_test(mean.N2Oflux ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = n2o.3rmaov2.cbrt2)
View(n2o.pwc.irr)
# using lmer mod
n2o.pwc.irr2 <- soilprop.n2o %>%
  group_by(Date, Fertilization) %>%
  emmeans_test(N2Oflux.cbrt ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = N2O.lmer.cbrt)
View(n2o.pwc.irr2)

# plot N2O
stat_text.n2o <- data.frame(Date = 0.5, mean.N2Oflux = 2,Fertilization="BIODYN", label="D *\nT ***\nD x T **\nC x T ***\nD x C x T **")


n2o.plot <- ggplot(soilprop.comp , aes(x=Date, y=mean.N2Oflux)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control','Biodyn-drought', 'Confym-control','Confym-drought',
                             'Conmin-control', 'Conmin-drought'))+
  ylab('Average N<sub>2</sub>O flux<br>(nmol m<sup>-2</sup>s<sup>-1)')+
  #ylab('Comammmox B abundance<br>(copies g<sup>-1</sup>dry soil)')+
  labs(subtitle = "C")+
  facet_wrap(~ Fertilization)+
  #scale_y_continuous(limits = c(0, NA))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        strip.text = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, hjust = 1, angle = 45),
        axis.title.y = element_markdown(size=18),
        plot.subtitle = element_text(size=20, face = "bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 geom_vline(xintercept = 3.4, linetype="dashed", colour="darkgrey") +
 annotate(geom = "text", x = 3.5, y = -0.23, hjust = 0, size = 4, label = "Rewetting", color = "gray25")+
 geom_label(data = stat_text.n2o,label=stat_text.n2o$label,hjust=0, colour="black", size=4, fontface="bold")
n2o.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
n2o.pwc.irr <- soilprop.n2o %>%
  group_by(Date, Fertilization) %>%
  emmeans_test(mean.N2Oflux ~ Irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = n2o.3rmaov2.cbrt2) 
# add x y position
n2o.pwc.xy <- n2o.pwc.irr  %>% 
  add_xy_position(x = "Date", dodge = 0.8) 
# plotting the pairwise comparisons among treatments (emmeans results)
n2o.plot2 <-  n2o.plot + 
  stat_pvalue_manual(n2o.pwc.xy,x = "Date", y.position = 3.2,
                     label = "p.adj.signif",size=5,
                     tip.length = 0.01, hide.ns = F)
n2o.plot2 # this one is using lmer
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Fig.7.dpi300.tiff",
       n2o.plot2, device = "tiff",bg="white",
       width = 15, height = 9.5, 
       units= "in", dpi = 300, compression="lzw")

N.pools <- (NH4.plot2 / NO3.plot2 / n2o.plot2)+
 plot_layout(guides = "collect") & 
 theme(legend.position = 'bottom',legend.title = element_blank())
N.pools

setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Fig.4.2.dpi300.tiff",
       N.pools, device = "tiff",bg="white",
       width = 11.5, height = 12, 
       units= "in", dpi = 300, compression="lzw")




