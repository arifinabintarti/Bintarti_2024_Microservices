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

# read the rooted tree
setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/AOB-rooted-tree/')
AOB_rooted_tree <- ape::read.tree("tree.nwk")

# make phyloseq object
aob.physeq <- merge_phyloseq(aob.asv.physeq,aob.tax.physeq,aob.meta.physeq,AOB_rooted_tree)
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
  rngseed =13, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
sort(sample_sums(aob.rare.min.physeq), decreasing = F) # 178 OTUs were removed because they are no longer present in any sample after random subsampling
                                                # no sample removed

# run the ggrare function attached in the file "generating_rarecurve.r"
#aob.rare.min <- ggrare(aob.rare.min.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)
test <- ggrare(aob.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

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

######################################################################################################

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
aob.min.meta$SampleID <- as.factor(aob.min.meta$SampleID)
aob.min.meta$PlotID <- as.factor(aob.min.meta$PlotID)
aob.min.meta$Irrigation<-as.factor(aob.min.meta$Irrigation)
aob.min.meta$Treatment<-as.factor(aob.min.meta$Treatment)
aob.min.meta$Date<-as.factor(aob.min.meta$Date)
aob.min.meta$Type<-as.factor(aob.min.meta$Type)

# 1. Richness

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
aob.min.meta.df$Type <- factor(aob.min.meta.df$Type, levels = c("BS", "RS"),
                  labels = c("Bulk Soil", "Rhizosphere"))
aob.min.meta.df$Treatment <- factor(aob.min.meta.df$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral fertilized"))
# tidy up the data frame
aob.min.meta.df.tidy <- aob.min.meta.df %>%
                             group_by(Irrigation, Treatment, Date,  Type, var2) %>%
                             summarize(Mean.Rich=mean(Richness),
                                       Mean.Sha=mean(Shannon),
                                       Mean.Simp=mean(Simpson),
                                       Mean.invsimp=mean(InvSimpson))
str(aob.min.meta.df.tidy)
#setwd('/Users/arifinabintarti/Documents/France/microservices/070623_AOB_out/')
#write.csv(aob.min.meta.df.tidy, file = "aob.min.meta.df.tidy2.csv")
#aob.min.meta.df.tidy.ed <- read.csv("aob.min.meta.df.tidy.csv")
#install.packages("rcartocolor")
library(rcartocolor)
carto_pal(n = NULL, 'Safe')
display_carto_pal(7, "Vivid")
carto_pal(n = NULL, 'Vivid')
color.trt <- c(D="#E58606", K="#5D69B1", M="#52BCA3")
#install.packages("ggnewscale")
library(ggnewscale)
#aob.min.meta.df$Date  <- as.Date(aob.min.meta.df$Date , "%m/%d/%Y")
aob.min.rich.plot <- ggplot(aob.min.meta.df.tidy, aes(x = Date, y = Mean.Rich, linetype=Irrigation))+
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

# 2. Shannon

#Line plot of AOB Shannon
#install.packages("scales")
library(scales)
aob.min.sha.plot <- ggplot(aob.min.meta.df.tidy, aes(x = Date, y = Mean.Sha, linetype=Irrigation))+
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

# 3. Simpson

#Line plot of AOB Simpson
aob.min.simp.plot <- ggplot(aob.min.meta.df.tidy, aes(x = Date, y = Mean.Simp, linetype=Irrigation))+
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

#. 4. Inverse Simpson

#Line plot of AOB Inverse Simpson
aob.min.invsimp.plot <- ggplot(aob.min.meta.df.tidy, aes(x = Date, y = Mean.invsimp, linetype=Irrigation))+
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

#############################################################################

# 1. Richness

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


# 2. Shannon

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


# Simpson

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


# Inverse Simpson

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

# 1. Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)

#aob.asv.min_PA <- 1*(aob.asv.min>0)
#aob.asv.min_PA
# Bray-Curtis using rarefied data:
aob.asv.min_dist <- vegdist(t(aob.asv.min), method = "bray")
# jaccard using rarefied data:
aob.asv.min_dist_jac <- vegdist(t(aob.asv.min), binary = TRUE, method = "jaccard") 
# Weighted UniFrac on unrarefied data:
aob.wUF_dist <- UniFrac(aob.physeq, weighted=TRUE, normalized = TRUE)
aob.wUF_dist
# Weighted UniFrac using rarefied data:
aob.rare_tree = phy_tree(aob.rare.min.physeq)
sprintf("Is tree binary: %s", is.binary(aob.rare_tree)) ## Check if binary (dichotomy) & multifurcating (polytomy) trees (please check: https://github.com/joey711/phyloseq/issues/1643)
phy_tree(aob.rare.min.physeq) = multi2di(aob.rare_tree)#If FASLE, randomly resolve polytomies and replace tree in "aob.rare.min.physeq"
sprintf("Is tree binary: %s", is.binary(phy_tree(aob.rare.min.physeq)))
aob.wUF.rare_dist <- UniFrac(aob.rare.min.physeq, weighted=TRUE, normalized = TRUE)
aob.wUF.rare_dist
# Unweighted UniFrac on unrarefied data:
aob.uwUF_dist <- UniFrac(aob.physeq, weighted=FALSE, normalized = TRUE)
aob.uwUF_dist
# Unweighted UniFrac using rarefied data:
aob.uwUF.rare_dist <- UniFrac(aob.rare.min.physeq, weighted=FALSE, normalized = TRUE)
aob.uwUF.rare_dist

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis using rarefied data:
aob.asv.min_pcoa <- cmdscale(aob.asv.min_dist, eig=T)
# jaccard using rarefied data:
aob.asv.min_pcoa.jac <- cmdscale(aob.asv.min_dist_jac, eig=T)
# Weighted UniFrac using unrarefied data:
aob.asv.min_pcoa.wUF <- cmdscale(aob.wUF_dist, eig=T)
# Weighted UniFrac using rarefied data:
aob.rare_pcoa.wUF <- cmdscale(aob.wUF.rare_dist, eig=T)
# Unweighted UniFrac using unrarefied data:
aob.asv.min_pcoa.uwUF <- cmdscale(aob.uwUF_dist, eig=T)
# Unweighted UniFrac using rarefied data:
aob.rare_pcoa.uwUF <- cmdscale(aob.uwUF.rare_dist, eig=T)

# 3. scores of PC1 and PC2

# bray-curtis:
ax1.scores <- aob.asv.min_pcoa$points[,1]
ax2.scores <- aob.asv.min_pcoa$points[,2] 
# jaccard:
ax1.scores.j <- aob.asv.min_pcoa.jac$points[,1]
ax2.scores.j <- aob.asv.min_pcoa.jac$points[,2]
# Weighted UniFrac:
ax1.scores.wUF <- aob.asv.min_pcoa.wUF$points[,1]
ax2.scores.wUF <- aob.asv.min_pcoa.wUF$points[,2]
# Weighted UniFrac using rarefied data:
ax1.scores.wUF.rare <- aob.rare_pcoa.wUF$points[,1]
ax2.scores.wUF.rare <- aob.rare_pcoa.wUF$points[,2]
# unweighted UniFrac:
ax1.scores.uwUF <- aob.asv.min_pcoa.uwUF$points[,1]
ax2.scores.uwUF <- aob.asv.min_pcoa.uwUF$points[,2]
# Unweighted UniFrac using rarefied data:
ax1.scores.uwUF.rare <- aob.rare_pcoa.uwUF$points[,1]
ax2.scores.uwUF.rare <- aob.rare_pcoa.uwUF$points[,2]

#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# 4. calculate percent variance explained, then add to plot

#Bray-curtis:
ax1 <- aob.asv.min_pcoa$eig[1]/sum(aob.asv.min_pcoa$eig)
ax2 <- aob.asv.min_pcoa$eig[2]/sum(aob.asv.min_pcoa$eig)
aob.map.pcoa <- cbind(aob.min.meta,ax1.scores,ax2.scores)
# jaccard
ax1.j <- aob.asv.min_pcoa.jac$eig[1]/sum(aob.asv.min_pcoa.jac$eig)
ax2.j <- aob.asv.min_pcoa.jac$eig[2]/sum(aob.asv.min_pcoa.jac$eig)
aob.map.pcoa.j <- cbind(aob.min.meta,ax1.scores.j,ax2.scores.j)
# Weighted UniFrac
ax1.wUF.rare <- aob.rare_pcoa.wUF$eig[1]/sum(aob.rare_pcoa.wUF$eig)
ax2.wUF.rare <- aob.rare_pcoa.wUF$eig[2]/sum(aob.rare_pcoa.wUF$eig)
aob.map.pcoa.wUF.rare <- cbind(aob.min.meta,ax1.scores.wUF.rare,ax2.scores.wUF.rare)
# Weighted UniFrac using rarefied data:
ax1.wUF <- aob.asv.min_pcoa.wUF$eig[1]/sum(aob.asv.min_pcoa.wUF$eig)
ax2.wUF <- aob.asv.min_pcoa.wUF$eig[2]/sum(aob.asv.min_pcoa.wUF$eig)
aob.map.pcoa.wUF <- cbind(aob.min.meta,ax1.scores.wUF,ax2.scores.wUF)
# unweighted UniFrac
ax1.uwUF <- aob.asv.min_pcoa.uwUF$eig[1]/sum(aob.asv.min_pcoa.uwUF$eig)
ax2.uwUF <- aob.asv.min_pcoa.uwUF$eig[2]/sum(aob.asv.min_pcoa.uwUF$eig)
aob.map.pcoa.uwUF <- cbind(aob.min.meta,ax1.scores.uwUF,ax2.scores.uwUF)
# Unweighted UniFrac using rarefied data:
ax1.uwUF.rare <- aob.rare_pcoa.uwUF$eig[1]/sum(aob.rare_pcoa.uwUF$eig)
ax2.uwUF.rare <- aob.rare_pcoa.uwUF$eig[2]/sum(aob.rare_pcoa.uwUF$eig)
aob.map.pcoa.uwUF.rare <- cbind(aob.min.meta,ax1.scores.uwUF.rare,ax2.scores.uwUF.rare)

# simple plot
aob.pcoa_plot <- plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))

# 5. PCoA Plot 

#require("ggrepel")
library(ggrepel)
library(viridis)

# a. Bray-Curtis:

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

# b. Jaccard:

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
 jaccard
set.seed(13)

# c. Weighted UniFrac:

aob.pcoa_plot.wUF <- ggplot(data = aob.map.pcoa.wUF, aes(x=ax1.scores.wUF, y=ax2.scores.wUF, color=Treatment))+
            theme_bw()+
            geom_point(data = aob.map.pcoa.wUF, aes(x = ax1.scores.wUF, y = ax2.scores.wUF , shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PC1:\n",round(ax1.wUF,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PC2:\n",round(ax2.wUF,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB Weighted UniFrac")+
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
            legend.spacing.x = unit(0.05, 'cm'))+
            stat_ellipse()
aob.pcoa_plot.wUF
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_WeightedUniFrac.tiff",
       aob.pcoa_plot.wUF, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# d. Weighted UniFrac ussing rarefied data:

aob.pcoa_plot.wUF.rare <- ggplot(data = aob.map.pcoa.wUF.rare, aes(x=ax1.scores.wUF.rare, y=ax2.scores.wUF.rare, color=Treatment))+
            theme_bw()+
            geom_point(data = aob.map.pcoa.wUF.rare, aes(x = ax1.scores.wUF.rare, y = ax2.scores.wUF.rare , shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PC1:\n",round(ax1.wUF.rare,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PC2:\n",round(ax2.wUF.rare,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB Weighted UniFrac (rarefied)")+
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
            legend.spacing.x = unit(0.05, 'cm'))+
            stat_ellipse()
aob.pcoa_plot.wUF.rare
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_WeightedUniFrac_rarefied.tiff",
       aob.pcoa_plot.wUF.rare, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# e. UnWeighted UniFrac:

aob.pcoa_plot.uwUF <- ggplot(data = aob.map.pcoa.uwUF, aes(x=ax1.scores.uwUF, y=ax2.scores.uwUF, color=Treatment))+
            theme_bw()+
            geom_point(data = aob.map.pcoa.uwUF, aes(x = ax1.scores.uwUF, y = ax2.scores.uwUF , shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PC1:\n",round(ax1.uwUF,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PC2:\n",round(ax2.uwUF,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB Unweighted UniFrac")+
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
            legend.spacing.x = unit(0.05, 'cm'))+
            stat_ellipse()
aob.pcoa_plot.uwUF
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_UnweightedUniFrac.tiff",
       aob.pcoa_plot.uwUF, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# f. UnWeighted UniFrac using rarefied data:

aob.pcoa_plot.uwUF.rare <- ggplot(data = aob.map.pcoa.uwUF.rare, aes(x=ax1.scores.uwUF.rare, y=ax2.scores.uwUF.rare, color=Treatment))+
            theme_bw()+
            geom_point(data = aob.map.pcoa.uwUF.rare, aes(x = ax1.scores.uwUF.rare, y = ax2.scores.uwUF.rare, shape=Irrigation),size=5, alpha= 0.8)+
            #scale_colour_manual(values = c("#FC8D62","#8DA0CB","#66C2A5"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PC1:\n",round(ax1.uwUF.rare,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PC2:\n",round(ax2.uwUF.rare,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Treatment",  title = "AOB Unweighted UniFrac (rarefied)")+
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
            legend.spacing.x = unit(0.05, 'cm'))+
            stat_ellipse()
aob.pcoa_plot.uwUF.rare
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_UnweightedUniFrac_rarefied.tiff",
       aob.pcoa_plot.uwUF.rare, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)
###################################################################################
# PERMANOVA #

# Calculated the statistical analysis of beta diversity using PERMANOVA
set.seed(13)
# Bray-Curtis:
aob.asv.min_dist
aob.adonis <- adonis2(aob.asv.min_dist ~ Irrigation*Treatment*Date*Type, data=aob.min.meta, 
                 permutation=999,
                 method="bray", 
                 strata = NULL)
aob.adonis
aob.adonis.irri <- adonis2(aob.asv.min_dist ~ Irrigation, data=aob.min.meta, 
                 permutation=999,
                 method="bray", 
                 strata = NULL)
aob.adonis.irri #not signif
aob.adonis.irrixtrt <- adonis2(aob.asv.min_dist ~ Irrigation*Treatment, data=aob.min.meta, 
                 permutation=999,
                 method="bray", 
                 strata = NULL)
aob.adonis.irrixtrt #Irrigation             1    0.252 0.00667  1.6570  0.071 .  
                    #Treatment              2    8.766 0.23215 28.8371  0.001 ***
                    #Irrigation:Treatment   2    0.472 0.01251  1.5537  0.050 *  

# Jaccard:
set.seed(13)
aob.adonis.jac <- adonis2(aob.asv.min_dist_jac ~ Irrigation*Treatment*Date*Type, data=aob.min.meta, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
aob.adonis.jac

# Weighted UniFrac:
set.seed(13)
aob.adonis.wuF <- adonis2(aob.wUF_dist ~ Irrigation*Treatment*Date*Type, data=aob.min.meta, 
                 permutation=999, 
                 strata = NULL)
aob.adonis.wuF

# Weighted UniFrac using rarefied data:
set.seed(13)
aob.adonis.wuF.rare <- adonis2(aob.wUF.rare_dist ~ Irrigation*Treatment*Date*Type, data=aob.min.meta, 
                 permutation=999, 
                 strata = NULL)
aob.adonis.wuF.rare
aob.adonis.wuF.rare.irrixtrt <- adonis2(aob.wUF.rare_dist ~ Irrigation*Treatment, data=aob.min.meta, 
                 permutation=999,
                 strata = NULL)
aob.adonis.wuF.rare.irrixtrt #                      Df SumOfSqs      R2       F Pr(>F)    
                             #Irrigation             1  0.00616 0.00679  2.1888  0.092 .  
                             #Treatment              2  0.36905 0.40656 65.5750  0.001 ***
                             #Irrigation:Treatment   2  0.00914 0.01007  1.6248  0.136    

#Unweighted UniFrac:
set.seed(13)
aob.adonis.uwuF <- adonis2(aob.uwUF_dist ~ Irrigation*Treatment*Date*Type, data=aob.min.meta, 
                 permutation=999, 
                 strata = NULL)
aob.adonis.uwuF
#Unweighted UniFrac using rarefied data:
set.seed(13)
aob.adonis.uwuF.rare <- adonis2(aob.uwUF.rare_dist ~ Irrigation*Treatment*Date*Type, data=aob.min.meta, 
                 permutation=999, 
                 strata = NULL)
aob.adonis.uwuF.rare

############################################################################
# testing beta diversity analyses on unrarefied data
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

###################################################################################
# AOB Community Composition

# Phyloseq object of rarefied data and unrarefied data:
# 1. rarefied data
aob.rare.min.physeq
tax_table(aob.rare.min.physeq)
# merge taxa by species
aob.sp <- tax_glom(aob.rare.min.physeq, taxrank = "Species", NArm = F)
aob.sp.ra <- transform_sample_counts(aob.sp, function(x) x/sum(x))
sample_data(aob.sp.ra)

aob.sp.df <- psmelt(aob.sp.ra) %>%
  group_by(var2, Type, Date, Treatment, Irrigation, Genus, Species) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

colours <- ColourPalleteMulti(aob.sp.df, "Genus", "Species")
colours
install.packages("Polychrome")
library(Polychrome)
# build-in color palette
mycol = glasbey.colors(21)
mycol=c("#FE8F42","#FFFFFF","#0000FF","#FF0000","#009FFF","#00FF00","#000033",
        
        "#FF00B6","#005300","#9A4D42","#00FFBE","#783FC1","#F1085C","#DD00FF",
        
        "#201A01","#1F9698","#FFD300","#FFACFD","#B1CC71", "#720055", "#766C95")

sp_col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777", "#44AAAA", "#77CCCC", "#737373", "#117744", "#88CCAA",
          "#777711", "#fc8d59", "#fed976", 
          "#771122", "#AA4455", "#DD7788","#774411", "#AA7744", "#DDAA77")


str(aob.sp.df)

install.packages("ggh4x")
library(ggh4x)

aob.sp.df$Date <- as.Date(aob.sp.df$Date , "%m/%d/%Y")
aob.sp.df$Type <- factor(aob.sp.df$Type, levels = c("BS", "RS"),
                  labels = c("Bulk Soil", "Rhizosphere")
                  )
aob.sp.df$Treatment <- factor(aob.sp.df$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral"))
legend <- "AOB Taxa"
library(scales)
#x_cols <- rep(hue_pal()(length(unique(interaction(aob.sp.df$Date, aob.sp.df$Irrigation)))))
aob.sp.plot <- ggplot(aob.sp.df, aes(x=interaction(Date, Irrigation), y=Mean, fill=Species)) + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(legend, values=sp_col)+
                     #scale_fill_manual(values=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0"))+
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
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
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           #guides(fill=guide_legend(nrow=2,byrow=TRUE))
                           scale_y_continuous(expand = c(0,0))+
                           guides(x="axis_nested")
                           
aob.sp.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_meanRA_barplot2.eps",
       aob.sp.plot, device = "eps",
       width = 15, height =6, 
       units= "in", dpi = 600)

# 2. unrarefied data
aob.physeq
# merge taxa by species
aob.sp.unrare <- tax_glom(aob.physeq, taxrank = "Species", NArm = F)
aob.sp.unrare.ra <- transform_sample_counts(aob.sp.unrare, function(x) x/sum(x))
aob.sp.unrare.ra

aob.sp.unrare.df <- psmelt(aob.sp.unrare.ra) %>%
  group_by(var2, Type, Date, Treatment, Irrigation, Genus, Species) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

str(aob.sp.unrare.df)
sp_col_x=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777", "#44AAAA", "#77CCCC", "#737373", "#117744", "#88CCAA",
          "#777711", "#fc8d59", "#fed976", 
          "#771122", "#AA4455", "#DD7788","#774411", "#AA7744", "#DDAA77","black")
install.packages("ggh4x")
library(ggh4x)

aob.sp.unrare.df$Date <- as.Date(aob.sp.unrare.df$Date , "%m/%d/%Y")
aob.sp.unrare.df$Type <- factor(aob.sp.unrare.df$Type, levels = c("BS", "RS"),
                  labels = c("Bulk Soil", "Rhizosphere")
                  )
aob.sp.unrare.df$Treatment <- factor(aob.sp.unrare.df$Treatment, levels = c("D", "K", "M"),
                  labels = c("Biodynamic", "Conventional", "Mineral"))
legend <- "AOB Taxa"
library(scales)
#x_cols <- rep(hue_pal()(length(unique(interaction(aob.sp.df$Date, aob.sp.df$Irrigation)))))
aob.sp.unrare.plot <- ggplot(aob.sp.unrare.df, aes(x=interaction(Date, Irrigation), y=Mean, fill=Species)) + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(legend, values=sp_col_x)+
                     #scale_fill_manual(values=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0"))+
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
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
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           #guides(fill=guide_legend(nrow=2,byrow=TRUE))
                           scale_y_continuous(expand = c(0,0))+
                           guides(x="axis_nested")
                           
aob.sp.unrare.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')
ggsave("AOB_meanRA_unrare_barplot.eps",
       aob.sp.unrare.plot, device = "eps",
       width = 15, height =6, 
       units= "in", dpi = 600)

###################################################################################################
# Compare the analyses with rarefaction to 1282 reads (remove sample with reads < 1000)
###################################################################################################
# ASV Table
sort(colSums(aob.asv, na.rm = FALSE, dims = 1), decreasing = F)
set.seed(13)
aob.rare.1282.seq <- rarefy_even_depth(aob.physeq, sample.size = 1282,
  rngseed = 13, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
aob.rare.1282.seq # 1 samples removed (S11), 116 ASVs were removed
# run the ggrare function attached in the file "generating_rarecurve.r"
set.seed(13)
aob.rare.1282 <- ggrare(aob.rare.1282.seq, step = 1, color = "Type", label = "SampleID", se = FALSE)
#set up your own color palette
Palette <- c("#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(aob.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.aob.rare.1282 <- aob.rare.1282 + 
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
 
plot.aob.rare.1282

setwd('/Users/arifinabintarti/Documents/France/Figures/AOB/')

ggsave("AOB_rarecurve_1282.jpg",
       plot.aob.rare.1282, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)

sort(sample_sums(aob.physeq), decreasing = F)


# rarefy to minimum sequencing depth
set.seed(13)
aob.rare.min.physeq <- rarefy_even_depth(aob.physeq, sample.size = min(sample_sums(aob.physeq)),
  rngseed = 13, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
sort(sample_sums(aob.rare.min.physeq), decreasing = F) # 172 OTUs were removed because they are no longer present in any sample after random subsampling
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


















