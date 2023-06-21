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

# Build a phylogenetic tree using the ASV Table
BiocManager::install("DECIPHER")
library(DECIPHER) 
library(phangorn); packageVersion("phangorn")






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
sort(sample_sums(aob.rare.min.physeq), decreasing = F) # 178 OTUs were removed because they are no longer present in any sample after random subsampling
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
aob.min.inv.d <- diversity(t(aob.asv.min), index = 'invsimpson')
aob.min.meta <- data.frame(meta_micro) # make data frame of the map data
aob.min.meta$Richness <- aob.min.s
aob.min.meta$Shannon <- aob.min.h
aob.min.meta$Pielou <- aob.min.pielou
aob.min.meta$Simpson <- aob.min.d
aob.min.meta$InvSimpson <- aob.min.inv.d
# Statistical Analyses: Alpha Diversity
aob.min.meta$Irrigation<-as.factor(aob.min.meta$Irrigation)
aob.min.meta$Treatment<-as.factor(aob.min.meta$Treatment)
aob.min.meta$Date<-as.factor(aob.min.meta$Date)
aob.min.meta$Type<-as.factor(aob.min.meta$Type)

# 1. Richness

# Anova mixed model (no interaction) = all samples
set.seed(13)
str(aob.min.meta.df)
aob.min.meta.df$Irrigation <- as.factor(aob.min.meta.df$Irrigation)
aob.min.meta.df$Date <- as.factor(aob.min.meta.df$Date)
aob.min.mod <- lmer(aob.min.meta.df$Richness ~ Irrigation + Treatment + Type + Date + (1|PlotID), data=aob.min.meta.df)
car::Anova(aob.min.mod, type = "II")
# Anova mixed model (with interaction) = all samples
set.seed(13)
str(aob.min.meta.df)
aob.min.mod2 <- lmer(aob.min.meta.df$Richness ~ Irrigation*Treatment*Date*Type + (1|PlotID), data=aob.min.meta.df)
car::Anova(aob.min.mod2, type = "III")

# Line plot of AOB richness 
aob.min.meta.df <- data.frame(meta_micro)
aob.min.meta.df$Richness <- aob.min.s
aob.min.meta.df$Shannon <- aob.min.h
aob.min.meta.df$Pielou <- aob.min.pielou
aob.min.meta.df$Simpson <- aob.min.d
aob.min.meta.df$InvSimpson <- aob.min.inv.d
aob.min.meta.df$SampleID<-as.factor(aob.min.meta.df$SampleID)
aob.min.meta.df$PlotID<-as.factor(aob.min.meta$PlotID)
aob.min.meta.df$Date  <- as.Date(aob.min.meta.df$Date , "%m/%d/%Y")
str(aob.min.meta.df)
# tidy up the data frame
aob.min.meta.df.tidy <- aob.min.meta.df %>%
                             group_by(Irrigation, Treatment, Type, Date) %>%
                             summarize(Mean.Rich=mean(Richness),
                                       Mean.Sha=mean(Shannon),
                                       Mean.Simp=mean(Simpson),
                                       Mean.invsimp=mean(InvSimpson))
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/')
#write.csv(aob.min.meta.df.tidy, file = "aob.min.meta.df.tidy2.csv")
aob.min.meta.df.tidy.ed <- read.csv("aob.min.meta.df.tidy.csv")
#install.packages("rcartocolor")
library(rcartocolor)
carto_pal(n = NULL, 'Safe')
display_carto_pal(7, "Vivid")
carto_pal(n = NULL, 'Vivid')
color.trt <- c(D="#E58606", K="#5D69B1", M="#52BCA3")
#install.packages("ggnewscale")
library(ggnewscale)
aob.min.meta.df.tidy.ed$Type <- factor(aob.min.meta.df.tidy.ed$Type, levels = c("BS", "RS"),
                  labels = c("Bulk Soil", "Rhizosphere")
                  )
aob.min.meta.df.tidy.ed$Date  <- as.Date(aob.min.meta.df.tidy.ed$Date , "%m/%d/%Y")
aob.min.rich.plot <- ggplot(aob.min.meta.df.tidy.ed, aes(x = Date, y = Mean.Rich, linetype=Irrigation))+
                             geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
                             facet_wrap(~ Type, strip.position="top", nrow = 1)+
                             #expand_limits(y = 0)+
                             theme_bw() +
                             scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             #scale_color_manual(values = color.trt, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             labs(x="Time Point", y="AOB Richness")+
                             theme(legend.position = "bottom",
                             strip.text = element_text(size=18),
                             axis.text.y = element_text(size = 16),
                             axis.text.x = element_text(size = 16),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.title.y = element_text(size=18,face="bold"),
                             axis.title.x = element_blank(),
                             legend.title = element_text(size=15, face='bold'),
                             legend.text=element_text(size=14))+
                             #new_scale_color() +
                             geom_point(alpha=0.5,aes(col=Treatment), size=4)+
                             theme(legend.direction = "horizontal", legend.box = "vertical")
                             #scale_shape_manual(values=c(16,17), labels = c("Bulk Soil", "Rhizosphere"))
                             #scale_colour_manual(values=c('black', 'red'))
                             
aob.min.rich.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_rich3.eps",
       aob.min.rich.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

dev.off()

#Line plot of AOB Shannon
#install.packages("scales")
library(scales)
str(aob.min.meta.df.tidy.ed)
aob.min.sha.plot <- ggplot(aob.min.meta.df.tidy.ed, aes(x = Date, y = Mean.Sha, linetype=Irrigation))+
                             geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
                             facet_wrap(~ Type, strip.position="top", nrow = 1)+
                             theme_bw() +
                             scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             #scale_color_manual(values = color.trt, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             labs(x="Time Point", y="AOB Shannon")+
                             theme(legend.position = "bottom",
                             strip.text = element_text(size=18),
                             axis.text.y = element_text(size = 16),
                             axis.text.x = element_text(size = 16),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.title.y = element_text(size=18,face="bold"),
                             axis.title.x = element_blank(),
                             legend.title = element_text(size=15, face='bold'),
                             legend.text=element_text(size=14))+
                             #new_scale_color() +
                             geom_point(alpha=0.5,aes(col=Treatment), size=4)+
                             theme(legend.direction = "horizontal", legend.box = "vertical")
                             #scale_shape_manual(values=c(16,17), labels = c("Bulk Soil", "Rhizosphere"))
                             #scale_colour_manual(values=c('black', 'red'))+
                             #scale_y_continuous(labels = label_number(accuracy = 0.05))
                             
aob.min.sha.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_sha2.eps",
       aob.min.sha.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

#Line plot of AOB Simpson
str(aob.min.meta.df.tidy.ed)
aob.min.simp.plot <- ggplot(aob.min.meta.df.tidy.ed, aes(x = Date, y = Mean.Simp, linetype=Irrigation))+
                             geom_line(linewidth=1, aes(group = var2, col=Treatment))+
                             #ylim(0.9,0.97)+
                             facet_wrap(~ Type, strip.position="top", nrow = 1)+
                             #expand_limits(y = 0.9)+
                             theme_bw() +
                             scale_color_manual(values = color.trt, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             labs(x="Time Point", y="AOB Simpson")+
                             theme(legend.position = "bottom",
                             strip.text = element_text(size=18),
                             axis.text.y = element_text(size = 16),
                             axis.text.x = element_text(size = 16),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.title.y = element_text(size=18,face="bold"),
                             axis.title.x = element_blank(),
                             legend.title = element_text(size=15, face='bold'),
                             legend.text=element_text(size=14))+
                             #new_scale_color() +
                             geom_point(alpha=0.5,aes(col=Treatment), size=4)+
                             theme(legend.direction = "horizontal", legend.box = "vertical")
                             #scale_color_manual(values = color.trt)
                             #scale_shape_manual(values=c(16,17), labels = c("Bulk Soil", "Rhizosphere"))
                             #scale_colour_manual(values=c('black', 'red'))
                             #scale_y_continuous(labels = label_number(accuracy = 0.01))
                             
aob.min.simp.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_simp.eps",
       aob.min.simp.plot, device = cairo_ps,
       width = 8, height = 7, 
       units= "in", dpi = 600)

#Line plot of AOB Inverse Simpson
str(aob.min.meta.df.tidy.ed)
aob.min.invsimp.plot <- ggplot(aob.min.meta.df.tidy.ed, aes(x = Date, y = Mean.invsimp, linetype=Irrigation))+
                             geom_line(linewidth=1, aes(group = var2, col=Treatment))+
                             #ylim(0.9,0.97)+
                             facet_wrap(~ Type, strip.position="top", nrow = 1)+
                             #expand_limits(y = 0.9)+
                             theme_bw() +
                             scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             #scale_color_manual(values = color.trt, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
                             labs(x="Time Point", y="AOB Inverse Simpson")+
                             theme(legend.position = "bottom",
                             strip.text = element_text(size=18),
                             axis.text.y = element_text(size = 16),
                             axis.text.x = element_text(size = 16),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.title.y = element_text(size=18,face="bold"),
                             axis.title.x = element_blank(),
                             legend.title = element_text(size=15, face='bold'),
                             legend.text=element_text(size=14))+
                             #new_scale_color() +
                             geom_point(alpha=0.5,aes(col=Treatment), size=4)+
                             theme(legend.direction = "horizontal", legend.box = "vertical")
                             #scale_color_manual(values = color.trt)
                             #scale_shape_manual(values=c(16,17), labels = c("Bulk Soil", "Rhizosphere"))
                             #scale_colour_manual(values=c('black', 'red'))
                             #scale_y_continuous(labels = label_number(accuracy = 0.01))
                             
aob.min.invsimp.plot                            
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_invsimp.eps",
       aob.min.invsimp.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

# Box plot of AOB Richness

#install.packages("ggpattern")
#install.packages("sf")
library(sf)
library(ggpattern)
soiltype <- c("Bulk Soil", "Rhizosphere")
names(soiltype) <- c("BS", "RS")
aob.min.meta.df$Type <- factor(aob.min.meta.df$Type, levels = c("BS", "RS"),
                  labels = c("Bulk Soil", "Rhizosphere"))
aob.min.meta.df$Treatment <- factor(aob.min.meta.df$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral fertilized"))
set.seed(13)
#aob.min.rich.boxplot <- ggplot(aob.min.meta.df, aes(x=Treatment, y=Richness, pattern = Irrigation, fill = Treatment)) +
               #geom_boxplot(aes(fill = Treatment))+
               #geom_jitter(position = position_jitter(seed=13), alpha=0.3)+
               #theme_bw() +
               #labs(x="Treatment", y="AOB Richness")+
               #expand_limits(y = 0)+
               #labs(pattern="Irrigation")+
               #scale_colour_viridis(discrete=T)+
               #scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
               #scale_pattern_manual(values = c(Rainout = "stripe", Control = "none"))+
               #facet_grid(Type~Date)+
               #geom_boxplot_pattern(position = position_dodge(preserve = "single"), 
                     #color = "black", pattern_fill = "black", 
                     #pattern_angle = 45, pattern_density = 0.05, 
                     #pattern_spacing = 0.018, 
                     #pattern_key_scale_factor = 0.6)+
                     #theme(legend.title = element_text(size=15, face='bold'),
                     #legend.text = element_text(size=15),
                     #strip.text = element_text(size=18),
                     #axis.text.y = element_text(size = 18),
                     #axis.title.y = element_text(size=18,face="bold"),
                     #axis.text.x = element_blank(),
                     #axis.title.x =element_blank(),
                     #axis.ticks.x = element_blank(),
                     #plot.background = element_blank(),
                     #panel.grid.major = element_blank(),
                     #panel.grid.minor = element_blank())+
                     #guides(pattern = guide_legend(override.aes = list(fill = "white")),
                     #fill = guide_legend(override.aes = list(pattern = "none")))
#aob.min.rich.boxplot
   
#setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
#ggsave("AOB_min_rich_box.eps",
       #aob.min.rich.boxplot, device = "eps",
       #width = 12.5, height = 7, 
       #units= "in", dpi = 600)

# Richness: plotting the significance across treatment
aob.min.rich.pwc.plot <- ggplot(aob.min.meta, aes(x=Irrigation, y=Richness)) +
               geom_boxplot(aes(fill=Treatment))+
               #geom_jitter(position = position_jitter(seed=13), alpha=0.3)+
               theme_bw() +
               labs(y="AOB Richness")+
               expand_limits(y = 0)+
               labs(pattern="Irrigation")+
               scale_fill_viridis(discrete=T)+
               #scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
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
aob.min.rich.pwc.plot
# adding asterix from the stats analyses
aob.rich.pwc.all <- aob.min.meta %>%
  group_by(Type, Date, Irrigation) %>%
  pairwise_t_test(Richness ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif)
# plotting the significance
aob.rich.pwc.all <- aob.rich.pwc.all %>% add_xy_position(x = "Irrigation", dodge = 0.8)
aob.min.rich.pwc.plot2 <- aob.min.rich.pwc.plot + 
 stat_pvalue_manual(aob.rich.pwc.all,label = "p.adj.signif", tip.length = 0, hide.ns = TRUE)+
 scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.min.rich.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_rich_all_test.eps",
       aob.min.rich.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# Shannon: plotting the significance across treatment
aob.min.sha.pwc.plot <- ggplot(aob.min.meta, aes(x=Irrigation, y=Shannon)) +
               geom_boxplot(aes(fill = Treatment))+
               #geom_jitter(position = position_jitter(seed=13), alpha=0.3)+
               theme_bw() +
               labs(y="AOB Shannon")+
               #expand_limits(y = 0)+
               labs(pattern="Irrigation")+
               scale_fill_viridis(discrete=T)+
               #scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
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
aob.min.sha.pwc.plot
# adding asterix from the stats analyses
aob.sha.pwc.all <- aob.min.meta %>%
  group_by(Type, Date, Irrigation) %>%
  pairwise_t_test(Shannon ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif)
# plotting the significance
aob.sha.pwc.all <- aob.sha.pwc.all %>% add_xy_position(x = "Irrigation", dodge = 0.8)
aob.min.sha.pwc.plot2 <- aob.min.sha.pwc.plot + 
 stat_pvalue_manual(aob.sha.pwc.all,label = "p.adj.signif", tip.length = 0, hide.ns = TRUE)+
 scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.min.sha.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_sha_all.eps",
       aob.min.sha.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# Simpson: plotting the significance across treatment
aob.min.simp.pwc.plot <- ggplot(aob.min.meta, aes(x=Irrigation, y=Simpson)) +
               geom_boxplot(aes(fill = Treatment))+
               #geom_jitter(position = position_jitter(seed=13), alpha=0.3)+
               theme_bw() +
               labs(y="AOB Simpson")+
               #expand_limits(y = 0)+
               labs(pattern="Irrigation")+
               scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
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
aob.min.simp.pwc.plot
# adding asterix from the stats analyses
aob.simp.pwc.all <- aob.min.meta %>%
  group_by(Type, Date, Irrigation) %>%
  pairwise_t_test(Simpson ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif)
# plotting the significance
aob.simp.pwc.all <- aob.simp.pwc.all %>% add_xy_position(x = "Irrigation", dodge = 0.8)
aob.min.simp.pwc.plot2 <- aob.min.simp.pwc.plot + 
 stat_pvalue_manual(aob.simp.pwc.all,label = "p.adj.signif", tip.length = 0, hide.ns = TRUE)+
 scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.min.simp.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_simp_all.eps",
       aob.min.simp.pwc.plot2, device = "eps",
       width = 16, height =7, 
       units= "in", dpi = 600)


# Inverse Simpson: plotting the significance across treatment
aob.min.invsimp.pwc.plot <- ggplot(aob.min.meta, aes(x=Irrigation, y=InvSimpson)) +
               geom_boxplot(aes(fill = Treatment))+
               #geom_jitter(position = position_jitter(seed=13), alpha=0.3)+
               theme_bw() +
               labs(y="AOB Inverse Simpson")+
               #expand_limits(y = 0)+
               labs(pattern="Irrigation")+
               scale_fill_viridis(discrete=T)+
               #scale_fill_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
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
aob.min.invsimp.pwc.plot
# adding asterix from the stats analyses
aob.invsimp.pwc.all <- aob.min.meta %>%
  group_by(Type, Date, Irrigation) %>%
  pairwise_t_test(InvSimpson ~ Treatment, p.adjust.method = "BH") %>%
  select(-p, -p.signif)
# plotting the significance
aob.invsimp.pwc.all <- aob.invsimp.pwc.all %>% add_xy_position(x = "Irrigation", dodge = 0.8)
aob.min.invsimp.pwc.plot2 <- aob.min.invsimp.pwc.plot + 
 stat_pvalue_manual(aob.invsimp.pwc.all,label = "p.adj.signif", tip.length = 0, hide.ns = TRUE)+
 scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aob.min.invsimp.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_min_invsimp_all.eps",
       aob.min.invsimp.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

###################################################################################
# Beta Diversity Analyses on Rarefied Data: AOB
###################################################################################
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
#aob.asv.min_PA <- 1*(aob.asv.min>0)
#aob.asv.min_PA
aob.asv.min_dist <- vegdist(t(aob.asv.min), method = "bray")
# jaccard
aob.asv.min_dist_jac <- vegdist(t(aob.asv.min), binary = TRUE, method = "jaccard") 
# Weighted UniFrac
aob.wUF_dist <- UniFrac(aob.rare.min.physeq, weighted = TRUE)

# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
aob.asv.min_pcoa <- cmdscale(aob.asv.min_dist, eig=T)
# jaccard
aob.asv.min_pcoa.jac <- cmdscale(aob.asv.min_dist_jac, eig=T)

# scores of PC1 and PC2
ax1.scores <- aob.asv.min_pcoa$points[,1]
ax2.scores <- aob.asv.min_pcoa$points[,2] 
# jaccard
ax1.scores.j <- aob.asv.min_pcoa.jac$points[,1]
ax2.scores.j <- aob.asv.min_pcoa.jac$points[,2]

#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

#calculate percent variance explained, then add to plot
ax1 <- aob.asv.min_pcoa$eig[1]/sum(aob.asv.min_pcoa$eig)
ax2 <- aob.asv.min_pcoa$eig[2]/sum(aob.asv.min_pcoa$eig)
aob.map.pcoa <- cbind(aob.min.meta,ax1.scores,ax2.scores)
# jaccard
ax1.j <- aob.asv.min_pcoa.jac$eig[1]/sum(aob.asv.min_pcoa.jac$eig)
ax2.j <- aob.asv.min_pcoa.jac$eig[2]/sum(aob.asv.min_pcoa.jac$eig)
aob.map.pcoa.j <- cbind(aob.min.meta,ax1.scores.j,ax2.scores.j)

# simple plot
aob.pcoa_plot <- plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
# PCoA Plot 
#require("ggrepel")
library(ggrepel)
library(viridis)
set.seed(13)
aob.pcoa_plot <- ggplot(data = aob.map.pcoa, aes(x=ax1.scores, y=ax2.scores))+
            theme_bw()+
            geom_point(data = aob.map.pcoa, aes(x = ax1.scores, y = ax2.scores, col=Treatment, shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB PCoA Plot (Bray-Curtis)")+
            theme(legend.position="right",
            legend.title = element_text(size=15, face='bold'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            #plot.subtitle = element_text(size = 20, face = 'bold'),
            axis.text=element_text(size=16), 
            axis.title=element_text(size=17,face="bold"),
            legend.text=element_text(size=15),
            legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_PCoA_rare.tiff",
       aob.pcoa_plot, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)
# jaccard
set.seed(13)
aob.pcoa_plot.jac <- ggplot(data = aob.map.pcoa.j, aes(x=ax1.scores.j, y=ax2.scores.j))+
            theme_bw()+
            geom_point(data = aob.map.pcoa.j, aes(x = ax1.scores.j, y = ax2.scores.j, col=Treatment, shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.j,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2.j,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB PCoA Plot (Bray-Curtis)")+
            theme(legend.position="right",
            legend.title = element_text(size=15, face='bold'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            #plot.subtitle = element_text(size = 20, face = 'bold'),
            axis.text=element_text(size=16), 
            axis.title=element_text(size=17,face="bold"),
            legend.text=element_text(size=15),
            legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot.jac
#aob.pcoa_plot+geom_text_repel(aes(label = SampleID),size = 3, max.overlaps = Inf)
# Calculated the statistical analysis of beta diversity using PERMANOVA
set.seed(13)
aob.asv.min_dist
aob.adonis <- adonis2(aob.asv.min_dist ~ Irrigation*Treatment*Date*Type, data=aob.min.meta, 
                 permutation=999,
                 method="bray", 
                 strata = NULL)
aob.adonis
set.seed(13)
aob.adonis.jac <- adonis2(aob.asv.min_dist_jac ~ Irrigation*Treatment*Date*Type, data=aob.min.meta, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
aob.adonis.jac
############################################################################
# beta diversity on unrarefied data
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
aob.asv_PA <- 1*(aob.asv>0)
aob.asv_PA
aob.asv_dist <- vegdist(t(aob.asv_PA), binary = TRUE, method = "bray")
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
aob.asv_pcoa <- cmdscale(aob.asv_dist, eig=T)
# scores of PC1 and PC2
ax1.scores.x <- aob.asv_pcoa$points[,1]
ax2.scores.x <- aob.asv_pcoa$points[,2] 
#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
#calculate percent variance explained, then add to plot
ax1.x <- aob.asv_pcoa$eig[1]/sum(aob.asv_pcoa$eig)
ax2.x <- aob.asv_pcoa$eig[2]/sum(aob.asv_pcoa$eig)
aob.map.pcoa.x <- cbind(aob.min.meta,ax1.scores.x,ax2.scores.x)
# simple plot
aob.pcoa_plot.x <- plot(ax1.scores.x, ax2.scores.x, xlab=paste("PCoA1: ",round(ax1.x,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2.x,3)*100,"% var. explained", sep=""))
# PCoA Plot 
#require("ggrepel")
library(ggrepel)
library(viridis)
set.seed(13)
aob.pcoa_plot.x<- ggplot(data = aob.map.pcoa.x, aes(x=ax1.scores.x, y=ax2.scores.x))+
            theme_bw()+
            geom_point(data = aob.map.pcoa.x, aes(x = ax1.scores.x, y = ax2.scores.x, col=Treatment, shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.x,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2.x,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB PCoA Plot")+
            theme(legend.position="right",
            legend.title = element_text(size=15, face='bold'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            #plot.subtitle = element_text(size = 20, face = 'bold'),
            axis.text=element_text(size=16), 
            axis.title=element_text(size=17,face="bold"),
            legend.text=element_text(size=15),
            legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot.x
set.seed(13)
aob.adonis.x <- adonis2(aob.asv_dist ~ Irrigation*Treatment*Date, data=aob.min.meta, 
                 permutation=999,
                 method="bray", 
                 strata = NULL)
aob.adonis.x

























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






