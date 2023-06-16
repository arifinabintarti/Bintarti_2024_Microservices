#################################### Analysis of amoA gene of AOB Illumina MiSeq Data #####################################
##
# Date : 30 May 2023
# Author : Ari Fina BINTARTI

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
install.packages("car")
install.packages("multcomp")
install.packages("ape")
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

# SET THE WORKING DIRECTORY
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB.ASV-analysis')
wd <- print(getwd())
# load the asv table
aob.asv <- read.table('annotated.AOB.ASVs.counts.tsv', sep='\t', header=T, row.names = 1, check.names = FALSE)
dim(aob.asv) # 1338  192
sort(colSums(aob.asv, na.rm = FALSE, dims = 1), decreasing = F) # there are no asv that does not exist in at least one sample.
# load the taxonomy table
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/')
aob.tax <- read.csv("besthit.diamond.output.curateddb.AOB.ASVs.csv")
dim(aob.tax) # 1338
# load the metadata
setwd('/Users/arifinabintarti/Documents/France/microservices/')
meta_micro <- read.csv("meta_microservices.csv")
# load phylogenetic tree (nwk file)
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB.Phylogenetic-analysis/')
aob.tre <- ape::read.tree("tree.AOB.nwk")

############################################################################
# rarefaction curve
set.seed(13)
rarecurve(t(aob.asv), step=50, cex=0.5, lwd=2, ylab="ASV", label=F)
#BiocManager::install("phyloseq")
library(phyloseq)

## make a phyloseq object of the asv table, taxonomy table, metadata

# re-order the rownames of the asv table to match the colnames of the metadata.
re_order <- match(rownames(meta_micro), colnames(aob.asv))
aob.asv.ord  <- aob.asv[ ,re_order]
aob.asv.physeq = otu_table(aob.asv.ord, taxa_are_rows = TRUE) # asv table
sample_names(aob.asv.physeq)
# adding "S" for sample names to avoid possible problem later on
sample_names(aob.asv.physeq) <- paste0("S", sample_names(aob.asv.physeq))

# phyloseq object of the taxonomy table
aob.tax <- column_to_rownames(aob.tax, var = "ASVid")
aob.tax.physeq = tax_table(as.matrix(aob.tax)) # taxonomy table
 
# phyloseq object of the metadata
rownames(meta_micro) <- sample_names(aob.asv.physeq)
aob.meta.physeq <- sample_data(meta_micro)# meta data
sample_names(aob.meta.physeq)

# make phyloseq object
aob.physeq <- merge_phyloseq(aob.asv.physeq,aob.tax.physeq,aob.meta.physeq)
aob.physeq

# run the ggrare function attached in the file "generating_rarecurve.r"
set.seed(13)
aob.rare <- ggrare(aob.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
Palette <- c("#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(aob.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aob.rare <- aob.rare + 
 theme_bw()+
 scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
 scale_size_manual(values = 60)+
 #scale_fill_manual()+
 labs(title = "(a) AOB", )+
 theme( strip.text.x = element_text(size=14, face='bold'),
        axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x = element_text(size=15,face="bold"),
        legend.position = "right",
        legend.title = element_text(size=15),
        legend.text = element_text(size = 13),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
        ylab("Number of ASVs")+xlab("Reads")
 
plot.aob.rare

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')

ggsave("AOB_rarecurve.jpg",
       plot.aob.rare, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)

sort(sample_sums(aob.physeq), decreasing = F)


# rarefy to minimum sequencing depth
set.seed(13)
aob.rare.min.physeq <- rarefy_even_depth(aob.physeq, sample.size = min(sample_sums(aob.physeq)),
  rngseed = 13, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
sort(sample_sums(aob.rare.min), decreasing = F) # 178 OTUs were removed because they are no longer present in any sample after random subsampling
                                                # no sample removed

# run the ggrare function attached in the file "generating_rarecurve.r"
aob.rare.min <- ggrare(aob.rare.min.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
Palette <- c("#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(aob.rare.min.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aob.rare.min <- aob.rare.min + 
 theme_bw()+
 scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
 scale_size_manual(values = 60)+
 #scale_fill_manual()+
 labs(title = "(a) AOB", )+
 theme( strip.text.x = element_text(size=14, face='bold'),
        axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x = element_text(size=15,face="bold"),
        legend.position = "right",
        legend.title = element_text(size=15),
        legend.text = element_text(size = 13),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
        ylab("Number of ASVs")+xlab("Reads")
 
plot.aob.rare.min

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')

ggsave("AOB_rarecurve_min.jpg",
       plot.aob.rare.min, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)

# rarefy to 2k of sequencing depth
set.seed(13)
aob.rare.2k.physeq <- rarefy_even_depth(aob.physeq, sample.size = 2000,
  rngseed = 13, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
aob.rare.2k.physeq
sort(sample_sums(aob.rare.2k.physeq), decreasing = F) # 102 OTUs were removed because they are no longer present in any sample after random subsampling
                                                # 5 samples removed
# subset samples with sample sum less than 2000 reads
aob.physeq2k <- prune_samples(sample_sums(aob.physeq) < 2000, aob.physeq)
sort(sample_sums(aob.physeq2k), decreasing = F)



############################################################################

# Calculate the alpha diversity (Richness and Pielou's evenness, we also calculates Shannon index) 

#### AOB alpha diversity using rarefied data to the minimum sequencing depth (reads=858) ###

colSums(otu_table(aob.rare.min.physeq))
aob.asv.min <- as.data.frame(otu_table(aob.rare.min.physeq))
dim(aob.asv.min) # 1160 ASVs
aob.asv.min_pa <- 1*(aob.asv.min>0)
aob.min.s <- specnumber(aob.asv.min, MARGIN = 2) # richness
aob.min.richness <- as.data.frame(aob.min.s) 
aob.min.h <- diversity(t(aob.asv.min), index = 'shannon') # Shannon index
aob.min.shannon <- as.data.frame(aob.min.h)
aob.min.pielou <- aob.min.h/log(aob.min.s) # Pielou's evenness
aob.min.evenness <- as.data.frame(aob.min.pielou)
aob.min.d <- diversity(t(aob.asv.min), index = 'simpson') # Simpson index
aob.min.simpson <- as.data.frame(aob.min.d)
aob.min.meta <- data.frame(meta_micro) # make data frame of the map data
aob.min.meta$Richness <- aob.min.s
aob.min.meta$Shannon <- aob.min.h
aob.min.meta$Pielou <- aob.min.pielou
aob.min.meta$Simpson <- aob.min.d
# Statistical Analyses: Alpha Diversity
aob.min.meta$Irrigation<-as.factor(aob.min.meta$Irrigation)
aob.min.meta$Treatment<-as.factor(aob.min.meta$Treatment)
aob.min.meta$Date<-as.factor(aob.min.meta$Date)
aob.min.meta$Type<-as.factor(aob.min.meta$Type)

# 1. Richness

# Anova mixed model (no interaction)
aob.min.mod <- lmer(aob.min.meta$Richness ~ Irrigation + Treatment + Type + Date + (1|PlotID), data=aob.min.meta)
car::Anova(aob.min.mod, type = "II")
# Anova mixed model (with interaction)
aob.min.mod2 <- lmer(aob.min.meta$Richness ~ Irrigation + Treatment + Irrigation*Treatment + Date + Type + (1|PlotID), data=aob.min.meta)
car::Anova(aob.min.mod2, type = "III")
# model evaluation
# normality
plot(aob.min.mod)
qqnorm(resid(aob.min.mod))
qqline(resid(aob.min.mod)) #There is some deviation from from the expected normal line towards the tails, but overall the line looks straight and therefore pretty normal and suggests that the assumption is not violated
# Generate residual and predicted values
aob.min.mod_resids <- residuals(aob.min.mod)
aob.min.mod_preds <- predict(aob.min.mod)
plot(aob.min.mod_resids ~ aob.min.mod_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(aob.min.mod_resids) #data errors are normally distributed
kurtosis(resid(aob.min.mod),method = 'sample') #3.6 (in the range of -7 to 7)
# homoscedasticity
leveneTest(residuals(aob.min.mod) ~ factor(aob.min.meta$PlotID)) #Since the p value is greater than 0.05, we can say that the variance of the residuals is equal and therefore the assumption of homoscedasticity is met
hist(aob.min.mod_resids)
# posthoc
library(multcomp)
post.rich.trt <- glht(aob.min.mod, linfct = mcp(Treatment = "Tukey"), test = adjusted("BH"))
mcs.rich.trt = summary(post.rich.trt,
              test=adjusted("BH"))
# get the significant letter
label.rich.trt <- cld(mcs.rich.trt,
    level=0.05,
    decreasing=TRUE)
label.rich.trt.df <- as.data.frame(label.rich.trt$mcletters$Letters)
names(label.rich.trt.df)[names(label.rich.trt.df) == "label.rich.trt$mcletters$Letters"] <- "Letter"
label.rich.trt.df <- rownames_to_column(label.rich.trt.df, "Treatment")
label.rich.trt.df
# calculate the max value of the ricchness and put them in the same data frame of the group and the significant letter
aob.sum.rich.trt <- aob.min.meta %>%
  group_by(Treatment) %>% 
  summarize(max.rich=max(Richness))
aob.sum.rich.trt2 <- left_join(label.rich.trt.df,aob.sum.rich.trt, by='Treatment')
aob.sum.rich.trt2

# line chart of AOB richness 
aob.min.meta.df <- data.frame(meta_micro)
aob.min.meta.df$Richness <- aob.min.s
aob.min.meta.df$Shannon <- aob.min.h
aob.min.meta.df$Pielou <- aob.min.pielou
aob.min.meta.df$Simpson <- aob.min.d
aob.min.meta.df$SampleID<-as.factor(aob.min.meta.df$SampleID)
aob.min.meta.df$PlotID<-as.factor(aob.min.meta$PlotID)
aob.min.meta.df$Date  <- as.Date(aob.min.meta.df$Date , "%m/%d/%Y")
str(aob.min.meta.df)
# tidy up the data frame
aob.min.meta.df.tidy <- aob.min.meta.df %>%
                             group_by(Irrigation, Treatment, Type, Date) %>%
                             summarize(Mean.Rich=mean(Richness),
                                       Mean.Sha=mean(Shannon),
                                       Mean.Simp=mean(Simpson))
#setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/')
#write.csv(aob.min.meta.df.tidy, file = "aob.min.meta.df.tidy.csv")
aob.min.meta.df.tidy.ed <- read.csv("aob.min.meta.df.tidy.csv")
aob.min.meta.df.tidy.ed.bulk <- aob.min.meta.df.tidy.ed[]
#install.packages("rcartocolor")
library(rcartocolor)
carto_pal(n = NULL, 'Safe')
display_carto_pal(7, "Vivid")
carto_pal(n = NULL, 'Vivid')
color.trt <- c(D="#E58606", K="#5D69B1", M="#52BCA3")
           #CDRS="#99C945", CCMRS="#CC61B0", CKRS="#24796C")
#palette.colors(n = NULL, "Polychrome 36")
install.packages("ggnewscale")
library(ggnewscale)
aob.min.meta.df.tidy.ed$Date  <- as.Date(aob.min.meta.df.tidy.ed$Date , "%m/%d/%Y")
aob.min.rich.plot <- ggplot(aob.min.meta.df.tidy.ed, aes(x = Date, y = Mean.Rich, linetype=Irrigation))+
                             geom_line(size=1.15, aes(group = var2, col=Treatment))+
                             ylim(0, 110)+
                             theme_bw() +
                             scale_color_manual(values = color.trt, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             labs(x="Time Point", y="AOB Richness")+
                             theme(legend.position = "right",
                             axis.text.y = element_text(size = 16),
                             axis.text.x = element_text(size = 16),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.title.y = element_text(size=18,face="bold"),
                             axis.title.x = element_text(size=18,face="bold"),
                             legend.title = element_text(size=15, face='bold'),
                             legend.text=element_text(size=14))+
                             new_scale_color() +
                             geom_point(alpha=0.5,aes(col=Type), size=4)+
                             #scale_shape_manual(values=c(16,17), labels = c("Bulk Soil", "Rhizosphere"))
                             scale_colour_manual(values=c('black', 'red'), labels = c("Bulk Soil", "Rhizosphere"))
                             
aob.min.rich.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_rich.eps",
       aob.min.rich.plot, device = "eps",
       width = 10, height = 7, 
       units= "in", dpi = 600)

# Box plot of AOB Richness

aob.min.meta.df

#install.packages("ggpattern")
#options(timeout=600)
#install.packages("sf")
library(sf)
library(ggpattern)

soiltype <- c("Bulk Soil", "Rhizosphere")
names(soiltype) <- c("BS", "RS")
date <- c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13")
names(date) <- c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13")
label <- c("Biodynamic", "Conventional", "Mineral fertilized")
set.seed(13)
aob.min.rich.boxplot <- ggplot(aob.min.meta.df, aes(x=Treatment, y=Richness, pattern = Irrigation, fill = Treatment)) +
               geom_boxplot(aes(fill = Treatment))+
               #geom_jitter(position = position_jitter(seed=13), alpha=0.3)+
               theme_bw() +
               labs(x="Treatment", y="AOB Richness")+
               expand_limits(y = 0)+
               #labs(pattern="Irrigation")+
               scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
               #ylim(0,2.4)+
               scale_pattern_manual(values = c(Rainout = "stripe", Control = "none"))+
               facet_grid(Type ~ Date, labeller=labeller(Date=date,Type=soiltype),scales = 'free_x')+
               geom_boxplot_pattern(position = position_dodge(preserve = "single"), 
                     color = "black", pattern_fill = "black", 
                     pattern_angle = 45, pattern_density = 0.1, 
                     pattern_spacing = 0.018, 
                     pattern_key_scale_factor = 0.6) +
                     theme(legend.position="right",
                     legend.key.size = unit(1, 'cm'),
                     legend.title = element_text(size=15, face='bold'),
                     legend.text = element_text(size=15),
                     axis.text.y = element_text(size = 18),
                     axis.text.x = element_blank(),
                     strip.text.x = element_blank(),
                     strip.text.y = element_text(size=18),
                     plot.title = element_text(size = 20,face = 'bold'),
                     axis.title.x =element_blank(),
                     axis.title.y = element_markdown(size=18,face="bold"),
                     plot.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
                     guides(pattern = guide_legend(override.aes = list(fill = "white")),fill = "none")
                     #fill = guide_legend(override.aes = list(pattern = "none")))

aob.min.rich.boxplot
library(grid)
library(gridExtra)
aob.min.rich.boxplot.grob <- ggplotGrob(aob.min.rich.boxplot)
names(aob.min.rich.boxplot.grob$grobs)
aob.min.rich.boxplot.grob$layout



# 2. Shannon

# Anova mixed model (no interaction)
aob.min.sha.mod <- lmer(aob.min.meta$Shannon ~ Irrigation + Treatment + Type + Date + (1|PlotID), data=aob.min.meta)
car::Anova(aob.min.sha.mod, type = "II")
# Anova mixed model (with interaction)
aob.min.sha.mod2 <- lmer(aob.min.meta$Shannon ~ Irrigation + Treatment + Irrigation*Treatment + Date + Type + (1|PlotID), data=aob.min.meta)
car::Anova(aob.min.sha.mod2, type = "III")
# model evaluation
# normality
plot(aob.min.sha.mod)
qqnorm(resid(aob.min.sha.mod))
qqline(resid(aob.min.sha.mod)) #There is some deviation from from the expected normal line towards the tails, but overall the line looks straight and therefore pretty normal and suggests that the assumption is not violated
# Generate residual and predicted values
aob.min.sha.mod_resids <- residuals(aob.min.sha.mod)
aob.min.sha.mod_preds <- predict(aob.min.sha.mod)
plot(aob.min.sha.mod_resids ~ aob.min.sha.mod_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)                            
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(aob.min.sha.mod_resids) #data errors are not normally distributed
kurtosis(resid(aob.min.sha.mod),method = 'sample') #3.6 (in the range of -7 to 7)
# homoscedasticity
leveneTest(residuals(aob.min.sha.mod) ~ factor(aob.min.meta$PlotID)) #Since the p value is greater than 0.05, we can say that the variance of the residuals is equal and therefore the assumption of homoscedasticity is met
hist(aob.min.sha.mod_resids)








aob.min.meta.df %>% 
group_by(gear,carb) %>% 
mutate(mean=mean(disp)) %>% 
ggplot(aes(gear,mean,group=carb,color=carb))
+geom_line()+geom_point()




str(aob.min.meta)
































# 1. Bulk Soil: 
# compare alpha diversity between control and shelter, and among managements systems within treatment
aob.min.meta.bulk <- aob.min.meta[1:120,] # select only bulk soil sample
aob.min.meta$Irrigation<-as.factor(aob.min.meta$Irrigation)
aob.min.meta$Treatment<-as.factor(aob.min.meta$Treatment)
aob.min.meta$Date<-as.factor(aob.min.meta$Date)
# Richness: fit nested ANOVA
set.seed(13)
aob.min.nest.rich <- aov(aob.min.meta.bulk$Richness ~ aob.min.meta.bulk$Treatment / aob.min.meta.bulk$Irrigation)
summary(aob.min.nest.rich) # (p-value Irri & Treat = 0.839 & 1.03e-10 ***, F-val Irri & Treat = 0.42 & 0.047), no significant effect of irrigation and management systems to AOB richness
# double check the effect of irrigation on AOB richness
set.seed(13)
aob.rich.trt <- lm(aob.meta.bulk$Richness ~ Treatment, data=aob.meta.bulk, na.action=na.exclude)
aob.rich.trt
drop1(aob.rich.trt,~.,test="F") # not significant, similar p-val and F-val
# Shannon: fit nested ANOVA
set.seed(13)
aob.nest.sha <- aov(aob.meta.bulk$Shannon ~ aob.meta.bulk$Treatment / aob.meta.bulk$Irrigation)
summary(aob.nest.sha) # (p-value Irri & Treat = 0.35 & 0.95, F-val Irri & Treat = 1.1 & 0.048), no significant effect of irrigation and management systems to AOB shannon diversity index
# Evenness: fit nested ANOVA
set.seed(13)
aob.nest.pie <- aov(aob.meta.bulk$Pielou ~ aob.meta.bulk$Treatment / aob.meta.bulk$Irrigation)
summary(aob.nest.pie) # (p-value Irri & Treat = 0.78 & 0.71, F-val Irri & Treat = 0.35 & 0.33), no significant effect of irrigation and management systems to AOB evenness.

aob.mod <- lmer(aob.min.meta.bulk$Richness ~ Irrigation * Date + Treatment + (1|PlotID), data=aob.min.meta.bulk)
car::Anova(aob.mod, type = "III")









# 2. Rhizosphere Soil: 
# compare alpha diversity between control and shelter, and among managements systems within treatment
aob.meta.rhizo <- aob.meta[121:192,] # select only rhizosphere soil sample
aob.meta.rhizo$Irrigation<-as.factor(aob.meta.rhizo$Irrigation)
aob.meta.rhizo$Treatment<-as.factor(aob.meta.rhizo$Treatment)
# Richness: fit nested ANOVA
set.seed(13)
aob.nest.rich.rz <- aov(aob.meta.rhizo$Richness ~ aob.meta.rhizo$Treatment / aob.meta.rhizo$Irrigation)
summary(aob.nest.rich.rz) # (p-value Irri & Treat = 0.21 & 0.66, F-val Irri & Treat = 1.58 & 0.59), no significant effect of irrigation and management systems to AOB richness
# Shannon: fit nested ANOVA
set.seed(13)
aob.nest.sha.rz <- aov(aob.meta.rhizo$Shannon ~ aob.meta.rhizo$Treatment / aob.meta.rhizo$Irrigation)
summary(aob.nest.sha.rz) # (p-value Irri & Treat = 0.27 & 0.12, F-val Irri & Treat = 1.2 & 1.85), no significant effect of irrigation and management systems to AOB shannon diversity index
# Evenness: fit nested ANOVA
set.seed(13)
aob.nest.pie.rz <- aov(aob.meta.rhizo$Pielou ~ aob.meta.rhizo$Irrigation / aob.meta.rhizo$Treatment)
summary(aob.nest.pie.rz) # (p-value Irri & Treat = 0.38 & 0.29, F-val Irri & Treat = 0.75 & 1.25), no significant effect of irrigation and management systems to AOB evenness.

### Conclusion: Irrigation and management systems have no effect on the AOB Alpha diversity ###






