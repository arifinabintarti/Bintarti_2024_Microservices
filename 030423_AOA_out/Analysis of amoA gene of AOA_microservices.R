#############################################################################################
# Analysis of amoA gene of AOA Illumina MiSeq Data 
#############################################################################################
# Date : 13 July 2023
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
install.packages("devtools", dependencies = TRUE)
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
library(devtools)
library(phyloseq)

# SET THE WORKING DIRECTORY
setwd('/Users/arifinabintarti/Documents/France/microservices/030423_AOA_out/AOA.ASV-analysis')
setwd('D:/Fina/INRAE_Project/microservices/030423_AOA_out/AOA.ASV-analysis')
wd <- print(getwd())
# load the asv table
aoa.asv <- read.table('annotated.AOA.ASVs.counts.tsv', sep='\t', header=T, row.names = 1, check.names = FALSE)
dim(aoa.asv) # 646  192
sort(rowSums(aoa.asv, na.rm = FALSE, dims = 1), decreasing = F) # there are no asv that does not exist in at least one sample.
# load the taxonomy table
setwd('/Users/arifinabintarti/Documents/France/microservices/030423_AOA_out/')
setwd('D:/Fina/INRAE_Project/microservices/030423_AOA_out/')
aoa.tax <- read.table("besthit.diamond.output.curateddb.AOA.ASVs.edited.csv", sep = ';', header=T)
dim(aoa.tax) # 646
# load the metadata
setwd('/Users/arifinabintarti/Documents/France/microservices/')
setwd('D:/Fina/INRAE_Project/microservices/')
meta_micro <- read.csv("meta_microservices.csv")
# load phylogenetic tree (nwk file)
setwd('D:/Fina/INRAE_Project/microservices/030423_AOA_out/AOA-rooted-tree/')
AOA_rooted_tree <- ape::read.tree("tree.nwk")
AOA_rooted_tree

############################################################################
# rarefaction curve
set.seed(13)
rarecurve(t(aoa.asv), step=50, cex=0.5, lwd=2, ylab="ASV", label=F)
#BiocManager::install("phyloseq")

## make a phyloseq object of the asv table, taxonomy table, metadata

# re-order the rownames of the asv table to match the colnames of the metadata.
re_order <- match(rownames(meta_micro), colnames(aoa.asv))
aoa.asv.ord  <- aoa.asv[ ,re_order]
aoa.asv.physeq = otu_table(aoa.asv.ord, taxa_are_rows = TRUE) # asv table
sample_names(aoa.asv.physeq)
# adding "S" for sample names to avoid possible problem later on
sample_names(aoa.asv.physeq) <- paste0("S", sample_names(aoa.asv.physeq))

# phyloseq object of the taxonomy table
aoa.tax <- column_to_rownames(aoa.tax, var = "ASVid")
#row.names(aoa.tax) <- aoa.tax$ASVid
aoa.tax.physeq = tax_table(as.matrix(aoa.tax)) # taxonomy table
 
# phyloseq object of the metadata
meta_micro$Date <- factor(meta_micro$Date, levels = c("4/28/22", "06/01/2022", "07/05/2022", "7/20/22", "9/13/22"),
                          labels = c("04-28-22", "06-01-22", "07-05-22", "07-20-22", "09-13-22"))
rownames(meta_micro) <- sample_names(aoa.asv.physeq)
aoa.meta.physeq <- sample_data(meta_micro)# meta data
sample_names(aoa.meta.physeq)

# read the rooted tree
setwd('D:/Fina/INRAE_Project/microservices/030423_AOA_out/AOA-rooted-tree/')
AOA_rooted_tree <- ape::read.tree("tree.nwk")

# make phyloseq object
aoa.physeq <- merge_phyloseq(aoa.asv.physeq,aoa.tax.physeq,aoa.meta.physeq,AOA_rooted_tree)
aoa.physeq
sample_data(aoa.physeq)$SampleID <- paste0("S", sample_data(aoa.physeq)$SampleID)
sample_data(aoa.physeq)

# run the ggrare function attached in the file "generating_rarecurve.r"
aoa.rare <- ggrare(aoa.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
Palette <- c("#FF7F00", "#662F00")
names(Palette) <- levels(sample_data(aoa.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aoa.rare <- aoa.rare + 
 theme_bw()+
 scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
 scale_size_manual(values = 60)+
 labs(title = "(b) AOA", )+
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
 
plot.aoa.rare

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_rarecurve.jpg",
       plot.aoa.rare, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)

setwd('D:/Fina/INRAE_Project/microservices_fig/AOA/')
ggsave("AOA_rarecurve.tiff",
       plot.aoa.rare, device = "tiff",
       width = 10, height = 7, 
       units= "in", dpi = 600)

# rarefy to minimum sequencing depth (minimum reads = 3832 reads)
set.seed(13)
aoa.rare.min.physeq <- rarefy_even_depth(aoa.physeq, sample.size = min(sample_sums(aoa.physeq)),
  rngseed = 13, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
aoa.rare.min.physeq
sort(sample_sums(aoa.rare.min.physeq), decreasing = F) # 54 OTUs were removed because they are no longer present in any sample after random subsampling
                                                # no sample removed
sort(rowSums(otu_table(aoa.rare.min.physeq), na.rm = FALSE, dims = 1), decreasing = F)

# run the ggrare function attached in the file "generating_rarecurve.r"
aoa.rare.min <- ggrare(aoa.rare.min.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)
#set up your own color palette
Palette <- c("#FF7F00", "#662F00")
names(Palette) <- levels(sample_data(aoa.rare.min.physeq)$Type)
Palette
legend_title <- "Sample Type"
# plot after rarefaction
library(ggtext)
plot.aoa.rare.min <- aoa.rare.min + 
                     theme_bw()+
                     scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
                     scale_size_manual(values = 60)+
                     labs(title = "(b) AOA", )+
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
                            ylab("Number of ASVs") + xlab("Reads")
 plot.aoa.rare.min

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_rarecurve_min.jpg",
       plot.aoa.rare.min, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA/')
ggsave("AOA_rarecurve_min.tiff",
       plot.aoa.rare.min, device = "tiff",
       width = 10, height = 7, 
       units= "in", dpi = 600)

########################################################################################################
# Calculate the alpha diversity (Richness, Shannon, and Inverse Simpson) 
########################################################################################################

colSums(otu_table(aoa.rare.min.physeq))
aoarare.asv.df <- as.data.frame(otu_table(aoa.rare.min.physeq))
dim(aoarare.asv.df) # 592 ASVs, 192 samples
aoarare.asv.df_pa <- 1*(aoarare.asv.df>0)
aoa.s <- specnumber(aoarare.asv.df, MARGIN = 2) # richness
aoa.richness <- as.data.frame(aoa.s) 
aoa.h <- diversity(t(aoarare.asv.df), index = 'shannon') # Shannon index
aoa.shannon <- as.data.frame(aoa.h)
aoa.d <- diversity(t(aoarare.asv.df), index = 'simpson') # Simpson index
aoa.simpson <- as.data.frame(aoa.d)
aoa.inv.d <- diversity(t(aoarare.asv.df), index = 'invsimpson')

# 1. Richness

# Line plot of AOA richness 
aoa.meta.df <- data.frame(meta_micro)
aoa.meta.df$Richness <- aoa.s
aoa.meta.df$Shannon <- aoa.h
aoa.meta.df$Simpson <- aoa.d
aoa.meta.df$InvSimpson <- aoa.inv.d
#aob.min.meta.df$Date  <- as.Date(aob.min.meta.df$Date , "%m/%d/%Y")
str(aoa.meta.df)
aoa.meta.df$Type <- factor(aoa.meta.df$Type, levels = c("BS", "RS"),
                               labels = c("Bulk Soil", "Rhizosphere"))
aoa.meta.df$Treatment <- factor(aoa.meta.df$Treatment, levels = c("D", "K", "M"),
                                    labels = c("Biodynamic", "Conventional", "Mineral fertilized"))
aoa.meta.df$SampleID<-as.factor(aoa.meta.df$SampleID)
aoa.meta.df$PlotID<-as.factor(aoa.meta.df$PlotID)
aoa.meta.df$Irrigation<-as.factor(aoa.meta.df$Irrigation)
# tidy up the data frame
aoa.meta.df.tidy <- aoa.meta.df %>%
  group_by(Irrigation, Treatment, Date, Type, var2,var3) %>%
  summarize(Mean.Rich=mean(Richness),
            Mean.Sha=mean(Shannon),
            Mean.Simp=mean(Simpson),
            Mean.invsimp=mean(InvSimpson))
str(aoa.meta.df.tidy)

# plotting line chart
#remotes::install_github("Nowosad/rcartocolor")
library(rcartocolor)
carto_pal(n = NULL, 'Safe')
carto_pal(n = NULL, 'Vivid')
display_carto_pal(7, "Vivid")
display_carto_pal(12, "Safe")
color.trt <- c(D="#E58606", K="#5D69B1", M="#52BCA3")
#install.packages("ggnewscale")
library(ggnewscale)
#install.packages("viridis")
library(viridis)
aoa.meta.df.tidy$Date <- factor(aoa.meta.df.tidy$Date, levels = unique(aoa.meta.df.tidy$Date))
aoa.rich.plot <- ggplot(aoa.meta.df.tidy, aes(x = Date, y = Mean.Rich, linetype=Irrigation))+
  geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
  facet_wrap(~ Type, strip.position="top", nrow = 1)+
  theme_bw() +
  scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
  labs(x="Time Point", y="AOA Richness")+
  theme(legend.position = "bottom",
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 16),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        legend.title = element_text(size=15, face='bold'),
        legend.text=element_text(size=14))+
  geom_point(alpha=0.5,aes(col=Treatment), size=4)+
  theme(legend.direction = "horizontal", legend.box = "vertical")
aoa.rich.plot                            

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_min_rich3.eps",
       aoa.min.rich.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_rich.tiff",
       aoa.rich.plot, device = tiff,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

# 2. Shannon

#Line plot of AOA Shannon
library(scales)
aoa.sha.plot <- ggplot(aoa.meta.df.tidy, aes(x = Date, y = Mean.Sha, linetype=Irrigation))+
  geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
  facet_wrap(~ Type, strip.position="top", nrow = 1)+
  theme_bw() +
  scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
  labs(x="Time Point", y="AOA Shannon")+
  theme(legend.position = "bottom",
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 16),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        legend.title = element_text(size=15, face='bold'),
        legend.text=element_text(size=14))+
  geom_point(alpha=0.5,aes(col=Treatment), size=4)+
  theme(legend.direction = "horizontal", legend.box = "vertical")
aoa.sha.plot  

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_min_sha2.eps",
       aoa.min.sha.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_sha.tiff",
       aoa.sha.plot, device = tiff,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

#. 3. Inverse Simpson

#Line plot of AOA Inverse Simpson
aoa.invsimp.plot <- ggplot(aoa.meta.df.tidy, aes(x = Date, y = Mean.invsimp, linetype=Irrigation))+
  geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
  facet_wrap(~ Type, strip.position="top", nrow = 1)+
  theme_bw() +
  scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
  labs(x="Time Point", y="AOA Inverse Simpson")+
  theme(legend.position = "bottom",
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 16),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        legend.title = element_text(size=15, face='bold'),
        legend.text=element_text(size=14))+
  geom_point(alpha=0.5,aes(col=Treatment), size=4)+
  theme(legend.direction = "horizontal", legend.box = "vertical")
aoa.invsimp.plot 

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_min_invsimp.eps",
       aoa.min.invsimp.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_invsimp.tiff",
       aoa.invsimp.plot, device = tiff,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

#############################################################################

# 1. Richness

# Box plot of AOB Richness

#install.packages('remotes')
#remotes::install_github("coolbutuseless/ggpattern")
#install.packages("rstatix")
#install.packages("sf")
library(remotes)
library(rstatix)
library(sf)
library(ggpattern)

# Richness: plotting the pairwise comparisons across treatment 
aoa.rich.pwc.plot <- ggplot(aoa.meta.df, aes(x=Irrigation, y=Richness)) +
  geom_boxplot(aes(fill=Treatment))+
  theme_bw() +
  labs(y="AOA Richness")+
  expand_limits(y = 0)+
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
aoa.rich.pwc.plot

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_min_rich_mean_boxplot.eps",
       aoa.min.rich.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_rich_mean_boxplot.tiff",
       aoa.rich.pwc.plot, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# richness between irrigations
aoa.rich.pwc.irri.plot <- ggplot(aoa.meta.df, aes(x=Date, y=Richness)) +
  geom_boxplot(aes(group = var3, fill = Irrigation))+
  theme_bw() +
  labs(y="AOA Richness")+
  scale_fill_manual(values = c("#996035","#F2DACD"))+
  facet_grid(Type~ Treatment,scales="free_x")+
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
aoa.rich.pwc.irri.plot

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_rich_irri_boxplot.eps",
       aoa.rich.pwc.irri.plot, device = "eps",
       width = 10, height =5.5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_rich_irri_boxplot.tiff",
       aoa.rich.pwc.irri.plot, device = "tiff",
       width = 10, height =5.5, 
       units= "in", dpi = 600)

# 2. Shannon

# Shannon: plotting the significance across treatment

aoa.sha.pwc.plot <- ggplot(aoa.meta.df, aes(x=Irrigation, y=Shannon)) +
  geom_boxplot(aes(fill = Treatment))+
  theme_bw() +
  labs(y="AOA Shannon")+
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
aoa.sha.pwc.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
aoa.xy.sha.bulk <- aoa.emm.sha.bulk %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
aoa.xy.sha.rh <- aoa.emm.sha.rh %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
# #combine two data frames and adding 'Type'
aoa.df.xy.sha.bulk <- as.data.frame(aoa.xy.sha.bulk)
aoa.df.xy.sha.rh <- as.data.frame(aoa.xy.sha.rh)
aoa.df.xy.sha.all <- rbind(aoa.df.xy.sha.bulk, aoa.df.xy.sha.rh) 
aoa.df.xy.sha.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
# plotting the pairwise comparisons among treatments (emmeans results)
aoa.sha.pwc.plot2 <- aoa.sha.pwc.plot + 
  stat_pvalue_manual(aoa.df.xy.sha.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aoa.sha.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_min_sha_boxplot.eps",
       aoa.min.sha.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_sha_boxplot.tiff",
       aoa.sha.pwc.plot2, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# shannon between irrigations

aoa.sha.pwc.irri.plot <- ggplot(aoa.meta.df, aes(x=Date, y=Shannon)) +
  geom_boxplot(aes(group = var3, fill = Irrigation))+
  theme_bw() +
  labs(y="AOA Shannon")+
  scale_fill_manual(values = c("#996035","#F2DACD"))+
  facet_grid(Type~ Treatment,scales="free_x")+
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
aoa.sha.pwc.irri.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_sha_irri_boxplot.eps",
       aoa.sha.pwc.irri.plot, device = "eps",
       width = 10, height =5.5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_sha_irri_boxplot.tiff",
       aoa.sha.pwc.irri.plot, device = "tiff",
       width = 10, height =5.5, 
       units= "in", dpi = 600)

# Inverse Simpson

# Inverse Simpson: plotting the significance across treatment
aoa.invsimp.pwc.plot <- ggplot(aoa.meta.df, aes(x=Irrigation, y=InvSimpson)) +
  geom_boxplot(aes(fill = Treatment))+
  theme_bw() +
  labs(y="AOA Inverse Simpson")+
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
aoa.invsimp.pwc.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
aoa.xy.invsimp.bulk <- aoa.emm.invsimp.bulk %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
aoa.xy.invsimp.rh <- aoa.emm.invsimp.rh %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
# #combine two data frames and adding 'Type'
aoa.df.xy.invsimp.bulk <- as.data.frame(aoa.xy.invsimp.bulk)
aoa.df.xy.invsimp.rh <- as.data.frame(aoa.xy.invsimp.rh)
aoa.df.xy.invsimp.all <- rbind(aoa.df.xy.invsimp.bulk, aoa.df.xy.invsimp.rh) 
aoa.df.xy.invsimp.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
# plotting the pairwise comparisons among treatments (emmeans results)
aoa.invsimp.pwc.plot2 <- aoa.invsimp.pwc.plot + 
  stat_pvalue_manual(aoa.df.xy.invsimp.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
aoa.invsimp.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_min_invsimp_all.eps",
       aoa.min.invsimp.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_invsimp_all.tiff",
       aoa.invsimp.pwc.plot2, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# inverse simpson between irrigations
aoa.invsimp.pwc.irri.plot <- ggplot(aoa.meta.df, aes(x=Date, y=InvSimpson)) +
  geom_boxplot(aes(group = var3, fill = Irrigation))+
  theme_bw() +
  labs(y="AOA Inverse Simpson")+
  scale_fill_manual(values = c("#996035","#F2DACD"))+
  facet_grid(Type~ Treatment,scales="free_x")+
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
aoa.invsimp.pwc.irri.plot

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_invsimp_irri_boxplot.eps",
       aoa.invsimp.pwc.irri.plot, device = "eps",
       width = 10, height =5.5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_invsimp_irri_boxplot.tiff",
       aoa.invsimp.pwc.irri.plot, device = "tiff",
       width = 10, height =5.5, 
       units= "in", dpi = 600)

###################################################################################
# Beta Diversity Analyses on Rarefied Data: AOA
###################################################################################

# FOR ALL SAMPLES
# 1. Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
aoa.asv.rare <- as.data.frame(otu_table(aoa.rare.min.physeq))
# Bray-Curtis using rarefied data:
aob_dist_bc <- vegdist(t(aob.asv.rare1k), method = "bray")
# Jaccard using rarefied data:
aob.dist_jac <- vegdist(t(aob.asv.rare1k), binary = TRUE, method = "jaccard")
# Weighted UniFrac using rarefied data:
aob.wUF_dist <- UniFrac(aob.rare.1282.seq, weighted=TRUE, normalized = TRUE)
aob.wUF_dist
# Unweighted UniFrac using rarefied data:
aob.uwUF_dist <- UniFrac(aob.rare.1282.seq, weighted=FALSE, normalized = TRUE)
aob.uwUF_dist

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis using rarefied data:
aob_pcoa.bc <- cmdscale(aob_dist_bc, eig=T)
# jaccard using rarefied data:
aob_pcoa.jac <- cmdscale(aob.dist_jac, eig=T)
# Weighted UniFrac using rarefied data:
aob_pcoa.wUF <- cmdscale(aob.wUF_dist, eig=T)
# Unweighted UniFrac using rarefied data:
aob_pcoa.uwUF <- cmdscale(aob.uwUF_dist, eig=T)

# 3. scores of PC1 and PC2

# bray-curtis:
ax1.scores <- aob_pcoa.bc$points[,1]
ax2.scores <- aob_pcoa.bc$points[,2] 
# jaccard:
ax1.scores.j <- aob_pcoa.jac$points[,1]
ax2.scores.j <- aob_pcoa.jac$points[,2]
# Weighted UniFrac using rarefied data:
ax1.scores.wUF <- aob_pcoa.wUF$points[,1]
ax2.scores.wUF <- aob_pcoa.wUF$points[,2]
# Unweighted UniFrac using rarefied data:
ax1.scores.uwUF <- aob_pcoa.uwUF$points[,1]
ax2.scores.uwUF <- aob_pcoa.uwUF$points[,2]

#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# 4. calculate percent variance explained, then add to plot

#Bray-curtis:
ax1 <- aob_pcoa.bc$eig[1]/sum(aob_pcoa.bc$eig)
ax2 <- aob_pcoa.bc$eig[2]/sum(aob_pcoa.bc$eig)
aob.map.pcoa <- cbind(aob.meta.df.sub,ax1.scores,ax2.scores)
# jaccard
ax1.j <- aob_pcoa.jac$eig[1]/sum(aob_pcoa.jac$eig)
ax2.j <- aob_pcoa.jac$eig[2]/sum(aob_pcoa.jac$eig)
aob.map.pcoa.j <- cbind(aob.meta.df.sub,ax1.scores.j,ax2.scores.j)
# Weighted UniFrac using rarefied:
ax1.wUF <- aob_pcoa.wUF$eig[1]/sum(aob_pcoa.wUF$eig)
ax2.wUF <- aob_pcoa.wUF$eig[2]/sum(aob_pcoa.wUF$eig)
aob.map.pcoa.wUF <- cbind(aob.meta.df.sub,ax1.scores.wUF,ax2.scores.wUF)
# unweighted UniFrac
ax1.uwUF <- aob_pcoa.uwUF$eig[1]/sum(aob_pcoa.uwUF$eig)
ax2.uwUF <- aob_pcoa.uwUF$eig[2]/sum(aob_pcoa.uwUF$eig)
aob.map.pcoa.uwUF <- cbind(aob.meta.df.sub,ax1.scores.uwUF,ax2.scores.uwUF)

# 5. PCoA Plot 

#require("ggrepel")
library(ggrepel)
library(viridis)

# a. Bray-Curtis:

set.seed(13)
aob.pcoa_plot <- ggplot(data = aob.map.pcoa, aes(x=ax1.scores, y=ax2.scores))+
  theme_bw()+
  geom_point(data = aob.map.pcoa, aes(x = ax1.scores, y = ax2.scores, col=Treatment, shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "AOB PCoA Plot (Bray-Curtis)")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
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

# b. Jaccard:

set.seed(13)
aob.pcoa_plot.jac <- ggplot(data = aob.map.pcoa.j, aes(x=ax1.scores.j, y=ax2.scores.j, color=Treatment))+
  theme_bw()+
  geom_point(data = aob.map.pcoa.j, aes(x = ax1.scores.j, y = ax2.scores.j, shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.j,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.j,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "AOB PCoA Plot (Bray-Curtis)")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot.jac
jaccard

# c. Weighted UniFrac using rarefied data:

aob.pcoa_plot.wUF <- ggplot(data = aob.map.pcoa.wUF, aes(x=ax1.scores.wUF, y=ax2.scores.wUF, color=Treatment))+
  theme_bw()+
  geom_point(data = aob.map.pcoa.wUF, aes(x = ax1.scores.wUF, y = ax2.scores.wUF , shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PC1:\n",round(ax1.wUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PC2:\n",round(ax2.wUF,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "AOB Weighted UniFrac (rarefied)")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot.wUF

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_WeightedUniFrac_rarefied.tiff",
       aob.pcoa_plot.wUF.rare, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# d. UnWeighted UniFrac using rarefied data:

aob.pcoa_plot.uwUF <- ggplot(data = aob.map.pcoa.uwUF, aes(x=ax1.scores.uwUF, y=ax2.scores.uwUF, color=Treatment))+
  theme_bw()+
  geom_point(data = aob.map.pcoa.uwUF, aes(x = ax1.scores.uwUF, y = ax2.scores.uwUF, shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PC1:\n",round(ax1.uwUF,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PC2:\n",round(ax2.uwUF,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "AOB Unweighted UniFrac (rarefied)")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))
aob.pcoa_plot.uwUF
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_UnweightedUniFrac_rarefied.tiff",
       aob.pcoa_plot.uwUF.rare, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

######################################################################################
# SEPARATE BETWEEN BULK SOIL AND RHIZOSPHERE
######################################################################################

# BULK SOIL

# 1. Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)

# Bray-Curtis - Bulk Soil :
aoarare.asv.df
aoa.asv.bulk <- aoarare.asv.df[,1:120]
aoa.asv.bulk1 <- aoa.asv.bulk[rowSums(aoa.asv.bulk)>0,]
sort(rowSums(aoa.asv.bulk1, na.rm = FALSE, dims = 1), decreasing = FALSE)
aoa.bulk_dist_bc <- vegdist(t(aoa.asv.bulk1), method = "bray")
# jaccard - Bulk Soil :
aoa.bulk_dist_jac <- vegdist(t(aoa.asv.bulk1), binary = TRUE, method = "jaccard")
# Weighted UniFrac (rarefied) - Bulk Soil:
aoa.physeq_bulk <- subset_samples(aoa.rare.min.physeq, Type=="BS")
aoa.physeq_bulk1 <- prune_taxa(taxa_sums(aoa.physeq_bulk)>0, aoa.physeq_bulk)
aoa.physeq_bulk1
sort(taxa_sums(aoa.physeq_bulk1), decreasing =F) #checking
aoa.bulk_dist_wUF <- UniFrac(aoa.physeq_bulk1, weighted=TRUE, normalized = TRUE)
aoa.bulk_dist_wUF
# Unweighted UniFrac (rarefied) -  Bulk Soil:
aoa.bulk_dist_uwUF <- UniFrac(aoa.physeq_bulk1, weighted=FALSE, normalized = TRUE)
aoa.bulk_dist_uwUF

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis - Bulk Soil:
aoa.bulk_pcoa_bc <- cmdscale(aoa.bulk_dist_bc, eig=T)
# Jaccard - Bulk Soil:
aoa.bulk_pcoa_jac <- cmdscale(aoa.bulk_dist_jac, eig=T)
# Weighted UniFrac - Bulk Soil:
aoa.bulk_pcoa_wUF <- cmdscale(aoa.bulk_dist_wUF, eig=T)
# Unweighted UniFrac - Bulk Soil:
aoa.bulk_pcoa.uwUF <- cmdscale(aoa.bulk_dist_uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis - Bulk Soil:
ax1.scores.bulk <- aoa.bulk_pcoa_bc$points[,1]
ax2.scores.bulk <- aoa.bulk_pcoa_bc$points[,2] 
# Jaccard - Bulk Soil:
ax1.scores.j.bulk <- aoa.bulk_pcoa_jac$points[,1]
ax2.scores.j.bulk <- aoa.bulk_pcoa_jac$points[,2]
# Weighted UniFrac - Bulk Soil:
ax1.scores.wUF.bulk <- aoa.bulk_pcoa_wUF$points[,1]
ax2.scores.wUF.bulk <- aoa.bulk_pcoa_wUF$points[,2]
# Unweighted UniFrac - Bulk Soil:
ax1.scores.uwUF.bulk <- aoa.bulk_pcoa.uwUF$points[,1]
ax2.scores.uwUF.bulk <- aoa.bulk_pcoa.uwUF$points[,2]

#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# 4. Envfit
env.aoa.bulk <- aoa.meta.df[1:120,c(13:19, 22,26:28)]
str(env.aoa.bulk)
env.aoa.bulk <- env.aoa.bulk %>% mutate_at(colnames(env.aoa.bulk), as.numeric)
# bray-curtis
set.seed(13)
env_fit.aoa.bc.bulk <- envfit(aoa.bulk_pcoa_bc, env.aoa.bulk, na.rm=TRUE)
env_fit.aoa.bc.bulk
# Jaccard 
set.seed(13)
env_fit.aoa.jac <- envfit(aoa.bulk_pcoa_jac, env.aoa.bulk, na.rm=TRUE)
# Weighted UniFrac
set.seed(13)
env_fit.aoa.wuF <- envfit(aoa.bulk_pcoa_wUF, env.aoa.bulk, na.rm=TRUE)
# UnWeighted UniFrac
set.seed(13)
env_fit.aoa.uwuF <- envfit(aoa.bulk_pcoa.uwUF, env.aoa.bulk, na.rm=TRUE)


# 5. calculate percent variance explained, then add to plot
aoa.meta.bulk <- aoa.meta.df[1:120,]
# Bray-curtis - Bulk Soil:
ax1.bulk <- aoa.bulk_pcoa_bc$eig[1]/sum(aoa.bulk_pcoa_bc$eig)
ax2.bulk <- aoa.bulk_pcoa_bc$eig[2]/sum(aoa.bulk_pcoa_bc$eig)
aoa.map.pcoa.bulk <- cbind(aoa.meta.bulk,ax1.scores.bulk,ax2.scores.bulk)
# Jaccard - Bulk Soil:
ax1.j.bulk <- aoa.bulk_pcoa_jac$eig[1]/sum(aoa.bulk_pcoa_jac$eig)
ax2.j.bulk <- aoa.bulk_pcoa_jac$eig[2]/sum(aoa.bulk_pcoa_jac$eig)
aoa.map.pcoa.j.bulk <- cbind(aoa.meta.bulk,ax1.scores.j.bulk,ax2.scores.j.bulk)
# Weighted UniFrac - Bulk Soil:
ax1.wUF.bulk <- aoa.bulk_pcoa_wUF$eig[1]/sum(aoa.bulk_pcoa_wUF$eig)
ax2.wUF.bulk <- aoa.bulk_pcoa_wUF$eig[2]/sum(aoa.bulk_pcoa_wUF$eig)
aoa.map.pcoa.wUF.bulk <- cbind(aoa.meta.bulk,ax1.scores.wUF.bulk,ax2.scores.wUF.bulk)
# Unweighted UniFrac - Bulk Soil:
ax1.uwUF.bulk <- aoa.bulk_pcoa.uwUF$eig[1]/sum(aoa.bulk_pcoa.uwUF$eig)
ax2.uwUF.bulk <- aoa.bulk_pcoa.uwUF$eig[2]/sum(aoa.bulk_pcoa.uwUF$eig)
aoa.map.pcoa.uwUF.bulk <- cbind(aoa.meta.bulk,ax1.scores.uwUF.bulk,ax2.scores.uwUF.bulk)

#################################################################################################

# RHIZOSPHERE

# Bray-Curtis - Rhizosphere :
aoa.asv.rh <- aoarare.asv.df[,121:192]
aoa.asv.rh1 <- aoa.asv.rh[rowSums(aoa.asv.rh)>0,]
sort(rowSums(aoa.asv.rh1, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(aoa.asv.rh1) 
aoa.rh_dist_bc <- vegdist(t(aoa.asv.rh1), method = "bray")
# jaccard - Rhizosphere :
aoa.rh_dist_jac <- vegdist(t(aoa.asv.rh1), binary = TRUE, method = "jaccard") 
# Weighted UniFrac (rarefied) - Rhizosphere :
aoa.physeq_rh <- subset_samples(aoa.rare.min.physeq, Type=="RS")
aoa.physeq_rh1 <- prune_taxa(taxa_sums(aoa.physeq_rh)>0, aoa.physeq_rh)
aoa.physeq_rh1
sort(taxa_sums(aoa.physeq_rh1), decreasing =F) #checking
aoa.rh_dist_wUF <- UniFrac(aoa.physeq_rh1, weighted=TRUE, normalized = TRUE)
aoa.rh_dist_wUF
# Unweighted UniFrac (rarefied) -  Rhizosphere :
aoa.rh_dist_uwUF <- UniFrac(aoa.physeq_rh1, weighted=FALSE, normalized = TRUE)
aoa.rh_dist_uwUF

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis - Rhizosphere :
aoa.rh_pcoa_bc <- cmdscale(aoa.rh_dist_bc, eig=T)
# Jaccard - Rhizosphere :
aoa.rh_pcoa_jac <- cmdscale(aoa.rh_dist_jac, eig=T)
# Weighted UniFrac - Rhizosphere :
aoa.rh_pcoa_wUF <- cmdscale(aoa.rh_dist_wUF, eig=T)
# Unweighted UniFrac - Rhizosphere :
aoa.rh_pcoa.uwUF <- cmdscale(aoa.rh_dist_uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis - Rhizosphere :
ax1.scores.rh <- aoa.rh_pcoa_bc$points[,1]
ax2.scores.rh <- aoa.rh_pcoa_bc$points[,2] 
# Jaccard - Rhizosphere :
ax1.scores.j.rh <- aoa.rh_pcoa_jac$points[,1]
ax2.scores.j.rh <- aoa.rh_pcoa_jac$points[,2]
# Weighted UniFrac - Rhizosphere :
ax1.scores.wUF.rh <- aoa.rh_pcoa_wUF$points[,1]
ax2.scores.wUF.rh <- aoa.rh_pcoa_wUF$points[,2]
# Unweighted UniFrac - Rhizosphere :
ax1.scores.uwUF.rh <- aoa.rh_pcoa.uwUF$points[,1]
ax2.scores.uwUF.rh <- aoa.rh_pcoa.uwUF$points[,2]

#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# 5. calculate percent variance explained, then add to plot
aoa.meta.rh <- aoa.meta.df[121:192,]
# Bray-curtis - Rhizosphere :
ax1.rh <- aoa.rh_pcoa_bc$eig[1]/sum(aoa.rh_pcoa_bc$eig)
ax2.rh <- aoa.rh_pcoa_bc$eig[2]/sum(aoa.rh_pcoa_bc$eig)
aoa.map.pcoa.rh <- cbind(aoa.meta.rh,ax1.scores.rh,ax2.scores.rh)
# Jaccard - Rhizosphere :
ax1.j.rh <- aoa.rh_pcoa_jac$eig[1]/sum(aoa.rh_pcoa_jac$eig)
ax2.j.rh <- aoa.rh_pcoa_jac$eig[2]/sum(aoa.rh_pcoa_jac$eig)
aoa.map.pcoa.j.rh <- cbind(aoa.meta.rh,ax1.scores.j.rh,ax2.scores.j.rh)
# Weighted UniFrac - Rhizosphere :
ax1.wUF.rh <- aoa.rh_pcoa_wUF$eig[1]/sum(aoa.rh_pcoa_wUF$eig)
ax2.wUF.rh <- aoa.rh_pcoa_wUF$eig[2]/sum(aoa.rh_pcoa_wUF$eig)
aoa.map.pcoa.wUF.rh <- cbind(aoa.meta.rh,ax1.scores.wUF.rh,ax2.scores.wUF.rh)
# Unweighted UniFrac - Rhizosphere :
ax1.uwUF.rh <- aoa.rh_pcoa.uwUF$eig[1]/sum(aoa.rh_pcoa.uwUF$eig)
ax2.uwUF.rh <- aoa.rh_pcoa.uwUF$eig[2]/sum(aoa.rh_pcoa.uwUF$eig)
aoa.map.pcoa.uwUF.rh <- cbind(aoa.meta.rh,ax1.scores.uwUF.rh,ax2.scores.uwUF.rh)

# 6. PCoA Plot 

#require("ggrepel")
library(ggrepel)
install.packages("viridis")
library(viridis)

# A. Bray-Curtis - Bulk Soil :
set.seed(33)
A.bc <- as.list(env_fit.aoa.bc.bulk$vectors) #shortcutting ef$vectors
pvals.bc<-as.data.frame(A.bc$pvals) #creating the dataframe
#environment scores (vectors scaled by R2 values)
env.scores1.bc <- as.data.frame(scores(env_fit.aoa.bc.bulk, display="vectors"))
env.scores2.bc <- cbind(env.scores1.bc, pvals)
env.scores3.bc <- cbind(env.scores2.bc,Variable=rownames(env.scores2.bc))
env.scores4.bc <- subset(env.scores3.bc,pvals<0.05)
set.seed(33)
mult <-.53

aoa.pcoa_bulk.plot <- ggplot(data = aoa.map.pcoa.bulk, aes(x=ax1.scores.bulk, y=ax2.scores.bulk, colour=Treatment))+
  theme_bw()+
  geom_point(data = aoa.map.pcoa.bulk, aes(x = ax1.scores.bulk, y = ax2.scores.bulk, shape=Irrigation),size=5, alpha= 0.6)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.bulk,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.bulk,3)*100,"% var. explained", sep=""))+
  geom_segment(data=env.scores4.bc,
               aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), 
               arrow = arrow(length = unit(0.3, "cm")),
               colour = "grey",inherit.aes = FALSE)+
  geom_text_repel(data = env.scores4.bc,
                  aes(x = mult*Dim1, y = mult*Dim2, label = Variable),
                  size = 5,fontface="bold",
                  position=position_jitter(width=0.03,height=0.001), inherit.aes = FALSE)+
  labs(colour = "Treatment",  title = "A. Bulk Soil")+
  theme(legend.position="none",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))+
  stat_ellipse()
aoa.pcoa_bulk.plot

# B. Bray-Curtis - Rhizosphere :
aoa.pcoa_rh.plot <- ggplot(data = aoa.map.pcoa.rh, aes(x=ax1.scores.rh, y=ax2.scores.rh, colour=Treatment))+
  theme_bw()+
  geom_point(data = aoa.map.pcoa.rh, aes(x = ax1.scores.rh, y = ax2.scores.rh, shape=Irrigation),size=5, alpha= 0.6)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.rh,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.rh,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "B. Rhizosphere")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))+
  stat_ellipse()
aoa.pcoa_rh.plot

#install.packages("patchwork")
library(patchwork)

aoa.bray.plot.envfit <- aoa.pcoa_bulk.plot |  aoa.pcoa_rh.plot
aoa.bray.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("aoa.bray.tiff",
       aoa.bray.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("aoa.bray.tiff",
       aoa.bray.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
ggsave("aoa.bray.envfit.tiff",
       aoa.bray.plot.envfit, device = "tiff",
       width = 16, height = 6, 
       units= "in", dpi = 600)

# A. Jaccard - Bulk Soil :
aoa.pcoa_bulk.jac <- ggplot(data = aoa.map.pcoa.j.bulk, aes(x=ax1.scores.j.bulk, y=ax2.scores.j.bulk, colour=Treatment))+
  theme_bw()+
  geom_point(data = aoa.map.pcoa.j.bulk, aes(x = ax1.scores.j.bulk, y = ax2.scores.j.bulk, shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.j.bulk,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.j.bulk,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "A. Bulk Soil")+
  theme(legend.position="none",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))+
  stat_ellipse()
aoa.pcoa_bulk.jac

# B. Jaccard - Rhizosphere :
aoa.pcoa_rh.jac <- ggplot(data = aoa.map.pcoa.j.rh, aes(x=ax1.scores.j.rh, y=ax2.scores.j.rh, colour=Treatment))+
  theme_bw()+
  geom_point(data = aoa.map.pcoa.j.rh, aes(x = ax1.scores.j.rh, y = ax2.scores.j.rh, shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.j.rh,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.j.rh,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "B. Rhizosphere")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))+
  stat_ellipse()
aoa.pcoa_rh.jac

aoa.jac.plot <- aoa.pcoa_bulk.jac |  aoa.pcoa_rh.jac
aoa.jac.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("aoa.jac.tiff",
       aoa.jac.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("aoa.jac.tiff",
       aoa.jac.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)

# A. Weighted UniFrac - Bulk Soil :
set.seed(33)
A.wu <- as.list(env_fit.aoa.wuF$vectors) #shortcutting ef$vectors
pvals.wu<-as.data.frame(A.wu$pvals) #creating the dataframe
#environment scores (vectors scaled by R2 values)
env.scores1.wu <- as.data.frame(scores(env_fit.aoa.wuF, display="vectors"))
env.scores2.wu <- cbind(env.scores1.wu, pvals)
env.scores3.wu <- cbind(env.scores2.wu,Variable=rownames(env.scores2.wu))
env.scores4.wu <- subset(env.scores3.wu,pvals<0.05)
set.seed(33)
mult <-.3

aoa.pcoa_bulk.wUF <- ggplot(data = aoa.map.pcoa.wUF.bulk, aes(x=ax1.scores.wUF.bulk, y=ax2.scores.wUF.bulk, colour=Treatment))+
  theme_bw()+
  geom_point(data = aoa.map.pcoa.wUF.bulk, aes(x = ax1.scores.wUF.bulk, y = ax2.scores.wUF.bulk, shape=Irrigation),size=5, alpha= 0.6)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.wUF.bulk,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.wUF.bulk,3)*100,"% var. explained", sep=""))+
  geom_segment(data=env.scores4.wu,
               aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), 
               arrow = arrow(length = unit(0.3, "cm")),
               colour = "grey",inherit.aes = FALSE)+
  geom_text_repel(data = env.scores4.wu,
                  aes(x = mult*Dim1, y = mult*Dim2, label = Variable),
                  size = 5,fontface="bold",
                  position=position_jitter(width=0.03,height=0.001), inherit.aes = FALSE)+
  labs(colour = "Treatment",  title = "A. Bulk Soil")+
  theme(legend.position="none",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))+
  stat_ellipse()
aoa.pcoa_bulk.wUF

# B. Weighted UniFrac - Rhizosphere :
aoa.pcoa_rh.wUF <- ggplot(data = aoa.map.pcoa.wUF.rh, aes(x=ax1.scores.wUF.rh, y=ax2.scores.wUF.rh, colour=Treatment))+
  theme_bw()+
  geom_point(data = aoa.map.pcoa.wUF.rh, aes(x = ax1.scores.wUF.rh, y = ax2.scores.wUF.rh, shape=Irrigation),size=5, alpha= 0.6)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.wUF.rh,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.wUF.rh,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "B. Rhizosphere")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))+
  stat_ellipse()
aoa.pcoa_rh.wUF

aoa.wUF.plot.envfit <- aoa.pcoa_bulk.wUF |  aoa.pcoa_rh.wUF
aoa.wUF.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("aoa.wUF.tiff",
       aoa.wUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("aoa.wUF.tiff",
       aoa.wUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
ggsave("aoa.wUF.envfit.tiff",
       aoa.wUF.plot.envfit, device = "tiff",
       width = 16, height = 6, 
       units= "in", dpi = 600)

# A. Unweighted UniFrac - Bulk Soil :
aoa.pcoa_bulk.uwUF <- ggplot(data = aoa.map.pcoa.uwUF.bulk, aes(x=ax1.scores.uwUF.bulk, y=ax2.scores.uwUF.bulk, colour=Treatment))+
  theme_bw()+
  geom_point(data = aoa.map.pcoa.uwUF.bulk, aes(x = ax1.scores.uwUF.bulk, y = ax2.scores.uwUF.bulk, shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.uwUF.bulk,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.uwUF.bulk,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "A. Bulk Soil")+
  theme(legend.position="none",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))+
  stat_ellipse()
aoa.pcoa_bulk.uwUF
#aob.pcoa_bulk.uwUF.id <- aob.pcoa_bulk.uwUF+geom_text_repel(aes(label = SampleID),size = 3, max.overlaps = Inf)
#aob.pcoa_bulk.uwUF.id
setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_UnweightedUniFrac_bulk_id.tiff",
       aoa.pcoa_bulk.uwUF.id, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# B. Unweighted UniFrac - Rhizosphere :
aoa.pcoa_rh.uwUF <- ggplot(data = aoa.map.pcoa.uwUF.rh, aes(x=ax1.scores.uwUF.rh, y=ax2.scores.uwUF.rh, colour=Treatment))+
  theme_bw()+
  geom_point(data = aoa.map.pcoa.uwUF.rh, aes(x = ax1.scores.uwUF.rh, y = ax2.scores.uwUF.rh, shape=Irrigation),size=5, alpha= 0.8)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.uwUF.rh,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.uwUF.rh,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "B. Rhizosphere")+
  theme(legend.position="right",
        legend.title = element_text(size=15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face="bold"),
        axis.text=element_text(size=16), 
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=15),
        legend.spacing.x = unit(0.05, 'cm'))+
  stat_ellipse()
aoa.pcoa_rh.uwUF

aoa.uwUF.plot <- aoa.pcoa_bulk.uwUF |  aoa.pcoa_rh.uwUF
aoa.uwUF.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("aoa.uwUF.tiff",
       aoa.uwUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("aoa.uwUF.tiff",
       aoa.uwUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)

############################################################################################
# PERMANOVA FOR BULK SOIL AND RHIZOSPHERE
############################################################################################

# A. Bray-Curtis - Bulk Soil : 
set.seed(13)
aoa.adonis.bulk <- adonis2(aoa.bulk_dist_bc ~ Irrigation*Treatment*Date, data=aoa.meta.bulk, 
                           permutation=999,
                           method="bray", 
                           strata = NULL) # only treatment is significant
aoa.adonis.bulk

set.seed(13)
aoa.adonis.bulk.irri <- adonis2(aoa.bulk_dist_bc ~ Irrigation, data=aoa.meta.bulk, 
                                permutation=999,
                                method="bray", 
                                strata = NULL) # not significant
aoa.adonis.bulk.irri

set.seed(13)
aoa.adonis.bulk.irri2 <- adonis2(aoa.bulk_dist_bc ~ Irrigation*Treatment, data=aoa.meta.bulk, 
                                 permutation=perm,method="bray") # not significant
aoa.adonis.bulk.irri2
perm = how(nperm = 999, within = Within(type="free"), 
           plots = with(aoa.meta.bulk, Plots(strata=Date, type="free")))
           #blocks = aoa.meta.bulk$Date)

set.seed(13)
aoa.adonis.bulk.trt <- adonis2(aoa.bulk_dist_bc ~ Treatment, data=aoa.meta.bulk, 
                               permutation=999,
                               method="bray", 
                               strata = NULL) # significant (p val = 0.001***)
aoa.adonis.bulk.trt

set.seed(13)
aoa.adonis.bulk.date <- adonis2(aoa.bulk_dist_bc ~ Date, data=aoa.meta.bulk, 
                                permutation=999,
                                method="bray", 
                                strata = NULL) # not significant
aoa.adonis.bulk.date

# B. Bray-Curtis - Rhizosphere : 
set.seed(13)
aoa.adonis.rh <- adonis2(aoa.rh_dist_bc ~ Irrigation*Treatment*Date, data=aoa.meta.rh, 
                         permutation=999,
                         method="bray", 
                         strata = NULL) # only treatment is significant
aoa.adonis.rh

set.seed(13)
aoa.adonis.rh.irri <- adonis2(aoa.rh_dist_bc ~ Irrigation, data=aoa.meta.rh, 
                              permutation=999,
                              method="bray", 
                              strata = NULL) # not significant
aoa.adonis.rh.irri

set.seed(13)
aoa.adonis.rh.irri2 <- adonis2(aoa.rh_dist_bc ~ Irrigation, data=aoa.meta.rh, 
                               permutation=999,
                               method="bray", 
                               strata = aoa.meta.rh$Treatment) # not significant
aoa.adonis.rh.irri2

set.seed(13)
aoa.adonis.rh.trt <- adonis2(aoa.rh_dist_bc ~ Treatment, data=aoa.meta.rh, 
                             permutation=999,
                             method="bray", 
                             strata = NULL) # treatment is significant ( p val = 0.001***)
aoa.adonis.rh.trt

set.seed(13)
aoa.adonis.rh.date <- adonis2(aoa.rh_dist_bc ~ Date, data=aoa.meta.rh, 
                              permutation=999,
                              method="bray", 
                              strata = NULL) # not significant
aoa.adonis.rh.date

# A. Jaccard - Bulk Soil : 
set.seed(13)
aoa.adonis.jac.bulk <- adonis2(aoa.bulk_dist_jac ~ Irrigation*Treatment*Date, data=aoa.meta.bulk, 
                               permutation=999,
                               method="jaccard", 
                               strata = NULL)
aoa.adonis.jac.bulk
# B. Jaccard - Rhizosphere : 
set.seed(13)
aoa.adonis.jac.rh <- adonis2(aoa.rh_dist_jac ~ Irrigation*Treatment*Date, data=aoa.meta.rh, 
                             permutation=999,
                             method="jaccard", 
                             strata = NULL)
aoa.adonis.jac.rh

# A. Weighted UniFrac - Bulk Soil : 
set.seed(13)
aoa.adonis.wuF.bulk <- adonis2(aoa.bulk_dist_wUF ~ Irrigation*Treatment*Date, data=aoa.meta.bulk, 
                               permutation=999, 
                               strata = NULL)
aoa.adonis.wuF.bulk

# B. Weighted UniFrac - Rhizosphere : 
set.seed(13)
aoa.adonis.wuF.rh <- adonis2(aoa.rh_dist_wUF ~ Irrigation*Treatment*Date, data=aoa.meta.rh, 
                             permutation=999, 
                             strata = NULL)
aoa.adonis.wuF.rh

# A. Unweighted UniFrac - Bulk Soil : 
set.seed(13)
aoa.adonis.uwuF.bulk <- adonis2(aoa.bulk_dist_uwUF ~ Irrigation*Treatment*Date, data=aoa.meta.bulk, 
                                permutation=999, 
                                strata = NULL)
aoa.adonis.uwuF.bulk
# B. Unweighted UniFrac - Rhizosphere : 
set.seed(13)
aoa.adonis.uwuF.rh <- adonis2(aoa.rh_dist_uwUF ~ Irrigation*Treatment*Date, data=aoa.meta.rh, 
                              permutation=999, 
                              strata = NULL)
aoa.adonis.uwuF.rh

########################################################################################
# Pairwise comparison analyses across treatments and between irrigation within date
########################################################################################
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# A. Pairwise Adonis Among Treatment (overall)

# 1. Bray-Curtis - Bulk Soil:
set.seed(13)
aoa.pw.bulk.trt_bc <- pairwiseAdonis::pairwise.adonis(aoa.bulk_dist_bc, 
                                                  aoa.meta.bulk$Treatment,
                                                  p.adjust.m = "BH")
aoa.pw.bulk.trt_bc # all pairwise comparisons are significant (p val =0.001**)

# 2. weighted UniFrac - Bulk Soil:
set.seed(13)
aoa.pw.bulk.trt_wUF <- pairwiseAdonis::pairwise.adonis(aoa.bulk_dist_wUF, 
                                                   aoa.meta.bulk$Treatment,
                                                   p.adjust.m = "BH")
aoa.pw.bulk.trt_wUF # all pairwise comparisons are significant (p val =0.001**)

# 3. Unweighted UniFrac - Bulk Soil:
set.seed(13)
aoa.pw.bulk.trt_uwUF <- pairwiseAdonis::pairwise.adonis(aoa.bulk_dist_uwUF, 
                                                    aoa.meta.bulk$Treatment,
                                                    p.adjust.m = "BH")
aoa.pw.bulk.trt_uwUF # all pairwise comparisons are significant (p val =0.001**)

# B. Pairwise Adonis Among Date

# 1. Bray-Curtis - Bulk Soil:
set.seed(13)
aoa.pw.bulk.dat_bc <- pairwiseAdonis::pairwise.adonis(aoa.bulk_dist_bc, 
                                                  aoa.meta.bulk$Date,
                                                  p.adjust.m = "BH")
aoa.pw.bulk.dat_bc # none are significant 

# 2. weighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.dat_wUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_wUF, 
                                                   aob.meta.bulk$Date,
                                                   p.adjust.m = "BH")
pw.bulk.dat_wUF # none are significant 

# 3. Unweighted UniFrac - Bulk Soil:
set.seed(13)
pw.bulk.dat_uwUF <- pairwiseAdonis::pairwise.adonis(aob.bulk_dist_uwUF, 
                                                    aob.meta.bulk$Date,
                                                    p.adjust.m = "BH")
pw.bulk.dat_uwUF # none are significant 

# C. Pairwise Adonis Between Irrigation within Treatment and Date

# separate the aob.asv table by date and treatment:

# 1.04-28-2022
# mineral:
m04.asv <- aob.asv.bulk1 %>% select(S1, S2, S9, S10, S15, S16, S23, S24)
m04.asv1 <- m04.asv[rowSums(m04.asv)>0,]
# biodynamic:
d04.asv <- aob.asv.bulk1 %>% select(S3, S4, S12, S13, S14, S21, S22)
d04.asv1 <- d04.asv[rowSums(d04.asv)>0,]
# conventional:
k04.asv <- aob.asv.bulk1 %>% select(S5, S6, S7, S8, S17, S18, S19, S20)
k04.asv1 <- k04.asv[rowSums(k04.asv)>0,]

# 2. 06-01-2022
# mineral:
m06.asv <- aob.asv.bulk1 %>% select(S25, S26, S33, S34, S39, S40, S47, S48)
m06.asv1 <- m06.asv[rowSums(m06.asv)>0,]
# biodynamic:
d06.asv <- aob.asv.bulk1 %>% select(S27, S28, S35, S36, S37, S38, S45, S46)
d06.asv1 <- d06.asv[rowSums(d06.asv)>0,]
# conventional:
k06.asv <- aob.asv.bulk1 %>% select(S29, S30, S31, S32, S41, S42, S43, S44)
k06.asv1 <- k06.asv[rowSums(k06.asv)>0,]

# 3. 07-05-2022
# mineral:
m0705.asv <- aob.asv.bulk1 %>% select(S49, S50, S57, S58, S63, S64, S71, S72)
m0705.asv1 <- m0705.asv[rowSums(m0705.asv)>0,]
# biodynamic:
d0705.asv <- aob.asv.bulk1 %>% select(S51, S52, S59, S60, S61, S62, S69, S70)
d0705.asv1 <- d0705.asv[rowSums(d0705.asv)>0,]
# conventional:
k0705.asv <- aob.asv.bulk1 %>% select(S53, S54, S55, S56, S65, S66, S67, S68)
k0705.asv1 <- k0705.asv[rowSums(k0705.asv)>0,]

# 4. 07-20-2022
# mineral:
m0720.asv <- aob.asv.bulk1 %>% select(S73, S74, S81, S82, S87, S88, S95, S96)
m0720.asv1 <- m0720.asv[rowSums(m0720.asv)>0,]
# biodynamic:
d0720.asv <- aob.asv.bulk1 %>% select(S75, S76, S83, S84, S85, S86, S93, S94)
d0720.asv1 <- d0720.asv[rowSums(d0720.asv)>0,]
# conventional:
k0720.asv <- aob.asv.bulk1 %>% select(S77, S78, S79, S80, S89, S90, S91, S92)
k0720.asv1 <- k0720.asv[rowSums(k0720.asv)>0,]

# 5. 09-13-2022
# mineral:
m09.asv <- aob.asv.bulk1 %>% select(S97, S98, S105, S106, S111, S112, S119, S120)
m09.asv1 <- m09.asv[rowSums(m09.asv)>0,]
# biodynamic:
d09.asv <- aob.asv.bulk1 %>% select(S99, S100, S107, S108, S109, S110, S117, S118)
d09.asv1 <- d09.asv[rowSums(d09.asv)>0,]
# conventional:
k09.asv <- aob.asv.bulk1 %>% select(S101, S102, S103, S104, S113, S114, S115, S116)
k09.asv1 <- k09.asv[rowSums(k09.asv)>0,]

# separate the metadata by date and treatment:
# 1. 
m04.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "04-28-22" 
                                & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d04.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "04-28-22" 
                                & aob.meta.bulk$Treatment == "Biodynamic"),]
k04.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "04-28-22" 
                                & aob.meta.bulk$Treatment == "Conventional"),]
# 2.  
m06.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "06-01-22" 
                                & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d06.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "06-01-22" 
                                & aob.meta.bulk$Treatment == "Biodynamic"),]
k06.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "06-01-22" 
                                & aob.meta.bulk$Treatment == "Conventional"),]
# 3. 
m0705.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-05-22" 
                                  & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d0705.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-05-22" 
                                  & aob.meta.bulk$Treatment == "Biodynamic"),]
k0705.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-05-22" 
                                  & aob.meta.bulk$Treatment == "Conventional"),]
# 4. 
m0720.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-20-22" 
                                  & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d0720.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-20-22" 
                                  & aob.meta.bulk$Treatment == "Biodynamic"),]
k0720.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "07-20-22" 
                                  & aob.meta.bulk$Treatment == "Conventional"),]
# 5. 
m09.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "09-13-22" 
                                & aob.meta.bulk$Treatment == "Mineral fertilized"),]
d09.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "09-13-22" 
                                & aob.meta.bulk$Treatment == "Biodynamic"),]
k09.meta <- aob.meta.bulk[which(aob.meta.bulk$Date == "09-13-22" 
                                & aob.meta.bulk$Treatment == "Conventional"),]

# calculate distance matrix - Bray Curtis
# 1.
m04dist_bc <- vegdist(t(m04.asv1), method = "bray")
d04dist_bc <- vegdist(t(d04.asv1), method = "bray")
k04dist_bc <- vegdist(t(k04.asv1), method = "bray")
# 2.
m06dist_bc <- vegdist(t(m06.asv1), method = "bray")
d06dist_bc <- vegdist(t(d06.asv1), method = "bray")
k06dist_bc <- vegdist(t(k06.asv1), method = "bray")
# 3.
m0705dist_bc <- vegdist(t(m0705.asv1), method = "bray")
d0705dist_bc <- vegdist(t(d0705.asv1), method = "bray")
k0705dist_bc <- vegdist(t(k0705.asv1), method = "bray")
# 4.
m0720dist_bc <- vegdist(t(m0720.asv1), method = "bray")
d0720dist_bc <- vegdist(t(d0720.asv1), method = "bray")
k0720dist_bc <- vegdist(t(k0720.asv1), method = "bray")
# 5.
m09dist_bc <- vegdist(t(m09.asv1), method = "bray")
d09dist_bc <- vegdist(t(d09.asv1), method = "bray")
k09dist_bc <- vegdist(t(k09.asv1), method = "bray")

# pairwise adonis between irrigation within treatment and date
# 1.
set.seed(13)
m04pwc.bc <- pairwiseAdonis::pairwise.adonis(m04dist_bc, m04.meta$Irrigation, p.adjust.m = "BH")
m04pwc.bc # not significant
set.seed(13)
d04pwc.bc <- pairwiseAdonis::pairwise.adonis(d04dist_bc, d04.meta$Irrigation, p.adjust.m = "BH")
d04pwc.bc # not significant
set.seed(13)
k04pwc.bc <- pairwiseAdonis::pairwise.adonis(k04dist_bc, k04.meta$Irrigation, p.adjust.m = "BH")
k04pwc.bc # not significant

# 2.
set.seed(13)
m06pwc.bc <- pairwiseAdonis::pairwise.adonis(m06dist_bc, m06.meta$Irrigation, p.adjust.m = "BH")
m06pwc.bc # not significant
set.seed(13)
d06pwc.bc <- pairwiseAdonis::pairwise.adonis(d06dist_bc, d06.meta$Irrigation, p.adjust.m = "BH")
d06pwc.bc # not significant
set.seed(13)
k06pwc.bc <- pairwiseAdonis::pairwise.adonis(k06dist_bc, k06.meta$Irrigation, p.adjust.m = "BH")
k06pwc.bc # not significant

# 3.
set.seed(13)
m0705pwc.bc <- pairwiseAdonis::pairwise.adonis(m0705dist_bc, m0705.meta$Irrigation, p.adjust.m = "BH")
m0705pwc.bc # not significant
set.seed(13)
d0705pwc.bc <- pairwiseAdonis::pairwise.adonis(d0705dist_bc, d0705.meta$Irrigation, p.adjust.m = "BH")
d0705pwc.bc # not significant
set.seed(13)
k0705pwc.bc <- pairwiseAdonis::pairwise.adonis(k0705dist_bc, k0705.meta$Irrigation, p.adjust.m = "BH")
k0705pwc.bc # not significant

# 4.
set.seed(13)
m0720pwc.bc <- pairwiseAdonis::pairwise.adonis(m0720dist_bc, m0720.meta$Irrigation, p.adjust.m = "BH")
m0720pwc.bc # not significant
set.seed(13)
d0720pwc.bc <- pairwiseAdonis::pairwise.adonis(d0720dist_bc, d0720.meta$Irrigation, p.adjust.m = "BH")
d0720pwc.bc # not significant
set.seed(13)
k0720pwc.bc <- pairwiseAdonis::pairwise.adonis(k0720dist_bc, k0720.meta$Irrigation, p.adjust.m = "BH")
k0720pwc.bc # not significant

# 5.
set.seed(13)
m09pwc.bc <- pairwiseAdonis::pairwise.adonis(m09dist_bc, m09.meta$Irrigation, p.adjust.m = "BH")
m09pwc.bc # not significant
set.seed(13)
d09pwc.bc <- pairwiseAdonis::pairwise.adonis(d09dist_bc, d09.meta$Irrigation, p.adjust.m = "BH")
d09pwc.bc # not significant
set.seed(13)
k09pwc.bc <- pairwiseAdonis::pairwise.adonis(k09dist_bc, k09.meta$Irrigation, p.adjust.m = "BH")
k09pwc.bc # not significant

# calculate distance matrix - Weighted UniFrac

# 1.
m04dist_bc <- vegdist(t(m04.asv1), method = "bray")
d04dist_bc <- vegdist(t(d04.asv1), method = "bray")
k04dist_bc <- vegdist(t(k04.asv1), method = "bray")
# 2.
m06dist_bc <- vegdist(t(m06.asv1), method = "bray")
d06dist_bc <- vegdist(t(d06.asv1), method = "bray")
k06dist_bc <- vegdist(t(k06.asv1), method = "bray")
# 3.
m0705dist_bc <- vegdist(t(m0705.asv1), method = "bray")
d0705dist_bc <- vegdist(t(d0705.asv1), method = "bray")
k0705dist_bc <- vegdist(t(k0705.asv1), method = "bray")
# 4.
m0720dist_bc <- vegdist(t(m0720.asv1), method = "bray")
d0720dist_bc <- vegdist(t(d0720.asv1), method = "bray")
k0720dist_bc <- vegdist(t(k0720.asv1), method = "bray")
# 5.
m09dist_bc <- vegdist(t(m09.asv1), method = "bray")
d09dist_bc <- vegdist(t(d09.asv1), method = "bray")
k09dist_bc <- vegdist(t(k09.asv1), method = "bray")

#########################################################################################
# AOA Community Composition
########################################################################################
# Phyloseq object of rarefied data and unrarefied data:
# 1. rarefied data
aoa.rare.min.physeq
tax_table(aoa.rare.min.physeq)
# merge taxa by species
aoa.sp <- tax_glom(aoa.rare.min.physeq, taxrank = "Sub_Clade2", NArm = F)
aoa.sp.ra <- transform_sample_counts(aoa.sp, function(x) x/sum(x))
sample_data(aoa.sp.ra)

aoa.sp.df <- psmelt(aoa.sp.ra) %>%
  group_by(var2, Type, Date, Treatment, Irrigation, Order, Clade, Sub_Clade2) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

aoa.order <- tax_glom(aoa.rare.min.physeq, taxrank = "Order", NArm = F)
aoa.order.ra <- transform_sample_counts(aoa.order, function(x) x/sum(x))
aoa.abund.trt <- psmelt(aoa.order.ra) %>%
  group_by(Type, Treatment, Order) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
aoa.cla <- tax_glom(aoa.rare.min.physeq, taxrank = "Clade", NArm = F)
aoa.cla.ra <- transform_sample_counts(aoa.cla, function(x) x/sum(x))
aoa.abund.trt.cla <- psmelt(aoa.cla.ra) %>%
  group_by(Type, Treatment, Clade) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)



#install.packages("Polychrome")
library(Polychrome)
# build-in color palette
#install.packages("colorBlindness")
library(colorBlindness)
displayAvailablePalette(color="white")
SteppedSequential5Steps
str(aoa.sp.df)
#install.packages("ggh4x")
library(ggh4x)
library(microshades)
microshades_palettes
microshades_cvd_palettes

grad.col.aoa <- c("#990F0F", "#B22C2C", "#CC5151", "#E57E7E", "#FFB2B2",
                  "#99540F", "#B26F2C", "#CC8E51", "#E5B17E", "#FFD8B2",
                  "#6B990F", "#85B22C", "#A3CC51", "#C3E57E", "#E5FFB2",
                  "#0F6B99", "#2C85B2", "#51A3CC", "#7EC3E5", "#B2E5FF",
                  "#260F99", "#422CB2", "#6551CC", "#8F7EE5", "#BFB2FF",
                  "#148F77", "#009E73", "#43BA8F", "#48C9B0", "#A3E4D7",
                  "#7D3560", "#A1527F", "#CC79A7", "#E794C1", "#EFB6D6")
                  

aoa.sp.df$Type <- factor(aoa.sp.df$Type, levels = c("BS", "RS"),
                         labels = c("Bulk Soil", "Rhizosphere"))
aoa.sp.df$Treatment <- factor(aoa.sp.df$Treatment, levels = c("D", "K", "M"),
                              labels = c("Biodynamic", "Conventional", "Mineral"))
aoa.sp.df$Sub_Clade2[aoa.sp.df$Mean<0.001] <- "Other (less than 0.1%)"
aoa.sp.df$Sub_Clade2[is.na(aoa.sp.df$Sub_Clade2)] <- "Other (less than 0.1%)"


legend <- "AOA Taxa"
library(scales)
library(scales)
library(forcats)
library(dplyr)
aoa.df <- aoa.sp.df
aoa.df$Sub_Clade <- paste(aoa.df$Order,aoa.df$Sub_Clade, sep = "-")

aoa.sp.plot <- ggplot(aoa.sp.df, aes(x=interaction(Date, Irrigation), y=Mean, fill=Sub_Clade2)) + 
  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(legend, values=grad.col.aoa)+
  facet_nested(~Type+Treatment, nest_line = element_line(linetype = 1), scales="free")+
  theme(legend.direction = "vertical",legend.position="right") + 
  guides(fill=guide_legend(ncol=1))+
  labs(y= "Mean Relative Abundance")+
  theme(plot.title = element_text(size = rel(1.5), face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=13, face = "bold"),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title =element_text(size=15,face="bold"),
        legend.text=element_text(size = 12),
        legend.title = element_text(size=13, face = "bold"),
        panel.grid = element_blank(), 
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA,linewidth= 0.2))+
  scale_y_continuous(expand = c(0,0))+ guides(x="axis_nested")
aoa.sp.plot

setwd('/Users/arifinabintarti/Documents/France/Figures/AOA/')
ggsave("AOA_meanRA_barplot2.eps",
       aoa.sp.plot, device = "eps",
       width = 15, height =6, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/AOA')
ggsave("AOA_meanRA_barplot.tiff",
       aoa.sp.plot, device = "tiff",
       width = 15, height =6, 
       units= "in", dpi = 600)



