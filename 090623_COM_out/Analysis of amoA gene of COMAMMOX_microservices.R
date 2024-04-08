#############################################################################################
# Analysis of amoA gene of Comammox Illumina MiSeq Data 
#############################################################################################
# Date : 1 August 2023
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
setwd('/Users/arifinabintarti/Documents/France/microservices/090623_COM_out/COM.ASV-analysis')
#setwd('D:/Fina/INRAE_Project/microservices/090623_COM_out/COM.ASV-analysis')
wd <- print(getwd())
# load the asv table
com.asv <- read.table('annotated.COM.ASVs.counts.tsv', sep='\t', header=T, row.names = 1, check.names = FALSE)
dim(com.asv) # 686 192
# remove the bad sample (sample # 26) from the OTU table
com.asv.sub <- com.asv[, -which(names(com.asv) == "26" )]
sort(rowSums(com.asv.sub, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(com.asv.sub)
# load the taxonomy table
setwd('/Users/arifinabintarti/Documents/France/microservices/090623_COM_out/')
#setwd('D:/Fina/INRAE_Project/microservices/090623_COM_out/')
com.tax <- read.table("besthit.diamond.output.curateddb.COM.ASVs.edited.csv", sep = ';', header=T)
dim(com.tax) # 680 
# load the metadata
setwd('/Users/arifinabintarti/Documents/France/microservices/')
#setwd('D:/Fina/INRAE_Project/microservices/')
meta_micro <- read.csv("meta_microservices.csv")
# remove the bad sample (sample # 26) from the metadata
meta_micro_sub <- meta_micro[-26,]
# load phylogenetic tree (nwk file)
setwd('/Users/arifinabintarti/Documents/France/microservices/090623_COM_out/COM-rooted-tree/')
#setwd('D:/Fina/INRAE_Project/microservices/090623_COM_out/COM-rooted-tree/')
COM_rooted_tree <- ape::read.tree("tree.nwk")

############################################################################
# rarefaction curve
set.seed(13)
rarecurve(t(com.asv), step=50, cex=0.5, lwd=2, ylab="ASV", label=F)
#BiocManager::install("phyloseq")

## make a phyloseq object of the asv table, taxonomy table, metadata

# re-order the rownames of the asv table to match the colnames of the metadata.
re_order <- match(rownames(meta_micro_sub), colnames(com.asv.sub))
com.asv.ord  <- com.asv.sub[ ,re_order]
com.asv.physeq = otu_table(com.asv.ord, taxa_are_rows = TRUE) # asv table
sample_names(com.asv.physeq)
# adding "S" for sample names to avoid possible problem later on
sample_names(com.asv.physeq) <- paste0("S", sample_names(com.asv.physeq))

# phyloseq object of the taxonomy table
com.tax <- column_to_rownames(com.tax, var = "ASVid")
#row.names(com.tax) <- com.tax$ASVid
com.tax.physeq = phyloseq::tax_table(as.matrix(com.tax)) # taxonomy table

# phyloseq object of the metadata
meta_micro_sub$Date <- factor(meta_micro_sub$Date, levels = c("4/28/22", "6/1/22", "7/5/22", "7/20/22", "9/13/22"),
                          labels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"))
rownames(meta_micro_sub) <- sample_names(com.asv.physeq)
com.meta.physeq <- sample_data(meta_micro_sub)# meta data
sample_names(com.meta.physeq)

# read the rooted tree
setwd('/Users/arifinabintarti/Documents/France/microservices/090623_COM_out/COM-rooted-tree/')
#setwd('D:/Fina/INRAE_Project/microservices/090623_COM_out/COM-rooted-tree/')
COM_rooted_tree <- ape::read.tree("tree.nwk")

# make phyloseq object
com.physeq <- merge_phyloseq(com.asv.physeq,com.tax.physeq,com.meta.physeq,COM_rooted_tree)
com.physeq # 653 taxa
sample_data(com.physeq)$SampleID <- paste0("S", sample_data(com.physeq)$SampleID)
sample_data(com.physeq)

# run the ggrare function attached in the file "generating_rarecurve.r"
com.rare <- ggrare(com.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)

#set up your own color palette
Palette <- c("#FF7F00", "#662F00")
names(Palette) <- levels(sample_data(com.physeq)$Type)
Palette
legend_title <- "Sample Type"

library(ggtext)
plot.com.rare <- com.rare + 
  theme_bw()+
  scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
  scale_size_manual(values = 60)+
  labs(title = "(c) COMAMMOX", )+
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
plot.com.rare

setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_rarecurve.jpg",
       plot.com.rare, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)

setwd('D:/Fina/INRAE_Project/microservices_fig/COM/')
ggsave("COM_rarecurve.tiff",
       plot.com.rare, device = "tiff",
       width = 10, height = 7, 
       units= "in", dpi = 600)

# rarefy to minimum sequencing depth (minimum reads =  reads)
set.seed(333)
com.rare.min.physeq <- rarefy_even_depth(com.physeq, sample.size = 5242,
                                         rngseed = 333, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
com.rare.min.physeq #632 taxa
sort(sample_sums(com.rare.min.physeq), decreasing = F) # 21 OTUs were removed because they are no longer present in any sample after random subsampling
# 1 sample removed (S52)
sort(rowSums(otu_table(com.rare.min.physeq), na.rm = FALSE, dims = 1), decreasing = F)

# run the ggrare function attached in the file "generating_rarecurve.r"
com.rare.min <- ggrare(com.rare.min.physeq, step = 1, color = "Type", label = "SampleID", se = FALSE)
#set up your own color palette
Palette <- c("#FF7F00", "#662F00")
names(Palette) <- levels(sample_data(com.rare.min.physeq)$Type)
Palette
legend_title <- "Sample Type"
# plot after rarefaction
library(ggtext)
plot.com.rare.min <- com.rare.min + 
  theme_bw()+
  scale_color_manual(legend_title,values = Palette, labels = c("Bulk Soil", "Rhizosphere"))+
  scale_size_manual(values = 60)+
  labs(title = "(c) COMAMMOX", )+
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
plot.com.rare.min

setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_rarecurve_min.jpg",
       plot.com.rare.min, device = "jpg",
       width = 10, height = 7, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM/')
ggsave("COM_rarecurve_min.tiff",
       plot.com.rare.min, device = "tiff",
       width = 10, height = 7, 
       units= "in", dpi = 600)

########################################################################################################
# Calculate the alpha diversity (Richness, Shannon, and Inverse Simpson) 
########################################################################################################

colSums(otu_table(com.rare.min.physeq))
com.rare.asv.df <- as.data.frame(otu_table(com.rare.min.physeq))
dim(com.rare.asv.df) # 632 ASVs, 190 samples
com.rare.asv.df_pa <- 1*(com.rare.asv.df>0)
com.s <- specnumber(com.rare.asv.df, MARGIN = 2) # richness
com.richness <- as.data.frame(com.s) 
com.h <- diversity(t(com.rare.asv.df), index = 'shannon') # Shannon index
com.shannon <- as.data.frame(com.h)
com.d <- diversity(t(com.rare.asv.df), index = 'simpson') # Simpson index
com.simpson <- as.data.frame(com.d)
com.inv.d <- diversity(t(com.rare.asv.df), index = 'invsimpson')

# 1. Richness

# Line plot of COM richness 
meta_micro_sub2 <- meta_micro_sub %>% filter(SampleID != "52")
com.meta.df <- data.frame(meta_micro_sub2)
com.meta.df$Richness <- com.s
com.meta.df$Shannon <- com.h
com.meta.df$Simpson <- com.d
com.meta.df$InvSimpson <- com.inv.d
#com.min.meta.df$Date  <- as.Date(com.meta.df$Date , "%m/%d/%Y")
str(com.meta.df)
com.meta.df$Type <- factor(com.meta.df$Type, levels = c("BS", "RS"),
                           labels = c("Bulk Soil", "Rhizosphere"))
com.meta.df$Treatment <- factor(com.meta.df$Treatment, levels = c("D", "K", "M"),
                                labels = c("BIODYN", "CONFYM", "CONMIN"))
com.meta.df$SampleID<-factor(com.meta.df$SampleID)
com.meta.df$PlotID<-factor(com.meta.df$PlotID)
com.meta.df$Irrigation<-factor(com.meta.df$Irrigation)
com.meta.df$Block<-factor(com.meta.df$Block)
com.meta.df$x<-factor(com.meta.df$x)
com.meta.df$DxD<-factor(com.meta.df$DxD)
com.meta.df$period<-factor(com.meta.df$period)
com.meta.df$rep2<-factor(com.meta.df$rep2)
com.meta.df$var2<-factor(com.meta.df$var2)
com.meta.df$var3<-factor(com.meta.df$var3)
com.meta.df[sapply(com.meta.df, is.character)] <- 
 lapply(com.meta.df[sapply(com.meta.df, is.character)], as.numeric)
com.meta.df[sapply(com.meta.df, is.integer)] <- 
 lapply(com.meta.df[sapply(com.meta.df, is.integer)], as.numeric)


# tidy up the data frame
com.meta.df.tidy <- com.meta.df %>%
  group_by(Irrigation, Treatment, Date, Type, var2,var3) %>%
  summarize(Mean.Rich=mean(Richness),
            Mean.Sha=mean(Shannon),
            Mean.Simp=mean(Simpson),
            Mean.invsimp=mean(InvSimpson))
str(com.meta.df.tidy)

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
com.meta.df.tidy$Date <- factor(com.meta.df.tidy$Date, levels = unique(com.meta.df.tidy$Date))
com.rich.plot <- ggplot(com.meta.df.tidy, aes(x = Date, y = Mean.Rich, linetype=Irrigation))+
  geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
  facet_wrap(~ Type, strip.position="top", nrow = 1)+
  theme_bw() +
  scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
  labs(x="Time Point", y="Comammox Richness")+
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
com.rich.plot                            

setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_min_rich3.eps",
       com.min.rich.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("COM_rich.tiff",
       com.rich.plot, device = tiff,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

# 2. Shannon

#Line plot of COM Shannon
library(scales)
com.sha.plot <- ggplot(com.meta.df.tidy, aes(x = Date, y = Mean.Sha, linetype=Irrigation))+
  geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
  facet_wrap(~ Type, strip.position="top", nrow = 1)+
  theme_bw() +
  scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
  labs(x="Time Point", y="Comammox Shannon")+
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
com.sha.plot  

setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_min_sha2.eps",
       com.min.sha.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("COM_sha.tiff",
       com.sha.plot, device = tiff,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

#. 3. Inverse Simpson

#Line plot of COM Inverse Simpson
com.invsimp.plot <- ggplot(com.meta.df.tidy, aes(x = Date, y = Mean.invsimp, linetype=Irrigation))+
  geom_line(linewidth=1.15, aes(group = var2, col=Treatment))+
  facet_wrap(~ Type, strip.position="top", nrow = 1)+
  theme_bw() +
  scale_colour_viridis(discrete=T, labels = c("Biodynamic", "Conventional", "Mineral fertilized"))+
  labs(x="Time Point", y="Comammox Inverse Simpson")+
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
com.invsimp.plot 

setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_min_invsimp.eps",
       com.min.invsimp.plot, device = cairo_ps,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("COM_invsimp.tiff",
       com.invsimp.plot, device = tiff,
       width = 9, height = 7.5, 
       units= "in", dpi = 600)

#############################################################################

# 1. Richness

# Box plot of COM Richness

#install.packages('remotes')
#remotes::install_github("coolbutuseless/ggpattern")
#install.packages("rstatix")
#install.packages("sf")
library(remotes)
library(rstatix)
library(sf)
library(ggpattern)

com.meta.df$x
com.meta.df.ed <- com.meta.df %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")

# Richness: plotting the pairwise comparisons across treatment 
com.rich.pwc.plot <- ggplot(com.meta.df, aes(x=Irrigation, y=Richness)) +
  geom_boxplot(aes(fill=Treatment))+
  theme_bw() +
  labs(y="Comammox Richness")+
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
com.rich.pwc.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
com.xy.rich.bulk <- com.emm.rich.bulk %>% 
add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
com.xy.rich.rh <- com.emm.rich.rh %>% 
add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
#combine two data frames and adding 'Type'
com.df.xy.rich.bulk <- as.data.frame(com.xy.rich.bulk)
com.df.xy.rich.rh <- as.data.frame(com.xy.rich.rh)
com.df.xy.rich.all <- rbind(com.df.xy.rich.bulk, com.df.xy.rich.rh) 
com.df.xy.rich.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
# plotting the pairwise comparisons among treatments (emmeans results)
com.rich.pwc.plot2 <- com.rich.pwc.plot + 
stat_pvalue_manual(com.df.xy.rich.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
com.rich.pwc.plot2

setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_min_rich_mean_boxplot.eps",
       com.min.rich.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("COM_rich_mean_boxplot.tiff",
       com.rich.pwc.plot2, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# richness between irrigations
com.rich.pwc.irri.plot <- ggplot(com.meta.df.ed, aes(x=Date, y=Richness)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  labs(y="Comammox Richness")+
  facet_grid(Type~ Treatment,scales="free_x")+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
com.rich.pwc.irri.plot

setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_rich_irri_boxplot.tiff",
       com.rich.pwc.irri.plot, device = "tiff",
       width = 10, height =5.5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("COM_rich_irri_boxplot.tiff",
       com.rich.pwc.irri.plot, device = "tiff",
       width = 10, height =5.5, 
       units= "in", dpi = 600)

# 2. Shannon

# Shannon: plotting the significance across treatment

com.sha.pwc.plot <- ggplot(com.meta.df, aes(x=Irrigation, y=Shannon)) +
  geom_boxplot(aes(fill = Treatment))+
  theme_bw() +
  labs(y="Comammox Shannon")+
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
com.sha.pwc.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
com.xy.sha.bulk <- com.emm.sha.bulk %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
com.xy.sha.rh <- com.emm.sha.rh %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
# combine two data frames and adding 'Type'
com.df.xy.sha.bulk <- as.data.frame(com.xy.sha.bulk)
com.df.xy.sha.rh <- as.data.frame(com.xy.sha.rh)
com.df.xy.sha.all <- rbind(com.df.xy.sha.bulk, com.df.xy.sha.rh) 
com.df.xy.sha.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
# plotting the pairwise comparisons among treatments (emmeans results)
com.sha.pwc.plot2 <- com.sha.pwc.plot + 
  stat_pvalue_manual(com.df.xy.sha.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
com.sha.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_min_sha_boxplot.eps",
       com.min.sha.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("COM_sha_boxplot.tiff",
       com.sha.pwc.plot2, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# shannon between irrigations

com.sha.pwc.irri.plot <- ggplot(com.meta.df.ed, aes(x=Date, y=Shannon)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  labs(y="Comammox Shannon")+
  facet_grid(Type~ Treatment,scales="free_x")+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
com.sha.pwc.irri.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_sha_irri_boxplot.tiff",
       com.sha.pwc.irri.plot, device = "tiff",
       width = 10, height =5.5, 
       units= "in", dpi = 600)
#setwd('D:/Fina/INRAE_Project/microservices_fig/COM')

#################################################################################################################
# Separate between Bulk Soil and Rhizosphere Alpha Diversity

#_Bulk Soil Richness_________________________________________________________________________________________________________________

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
com.meta.bulk$DxD<-factor(com.meta.bulk$DxD)
com.meta.bulk$period<-factor(com.meta.bulk$period)
com.meta.bulk$period2<-factor(com.meta.bulk$period2)

com.meta.bulk.edit <- com.meta.bulk
#com.meta.bulk.edit$Date <- factor(com.meta.bulk.edit$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          #labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
#com.meta.bulk.edit$Date <- as.Date(com.meta.bulk.edit$Date)

com.meta.bulk.edit$Date <- factor(com.meta.bulk.edit$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("Apr", "Jun", "Jul5", "Jul20", "Sep"))
com.meta.bulk.edit <- com.meta.bulk.edit %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")

#stat_text.BS2.com <- data.frame(Date = as.Date("2022-05-01"), Richness = 120,Treatment="BIODYN", label="C *\nD x C *")
stat_text.BS2.com <- data.frame(Date = 0.5, Richness = 10,Treatment="BIODYN", label="C *\nD x C *")

COM.BS.Richness.plot <- ggplot(com.meta.bulk.edit , aes(x=Date, y=Richness)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control', 'Biodyn-drought', 'Confym-control', 
                             'Confym-drought', 'Conmin-control', 'Conmin-drought'))+
  labs(y="COMA Richness", subtitle = "C")+
  facet_wrap(~ Treatment)+
  scale_y_continuous(limits = c(0, 120))+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        #strip.text = element_blank(),
        axis.text.y = element_text(size = 18),
        #axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=18),
        plot.subtitle = element_text(size=20, face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 #geom_vline(data=filter(aob.meta.df.sub.ed2, Type=="Bulk Soil"),aes(xintercept = as.Date("2022-07-14")), linetype="dashed", colour="darkgrey") +
 #geom_vline(xintercept = as.Date("2022-07-14"), linetype="dashed", colour="darkgrey") +
 annotate(geom = "text", x = 3.5, y = 0, hjust = 0, size = 4, label = "Rewetting", color = "grey25")+
 geom_vline(xintercept = 3.4, linetype="dashed", colour="darkgrey") +
 geom_label(data = stat_text.BS2.com,label=stat_text.BS2.com$label,hjust=0, colour="black", size=4, fontface="bold")
 
COM.BS.Richness.plot

#_Rhizosphere Richness________________________________________________________________________________________________________________________________________

com.meta.rh <- com.meta.df.ed[119:190,]
view(com.meta.rh)
str(com.meta.rh)
com.meta.rh$SampleID<-factor(com.meta.rh$SampleID)
com.meta.rh$PlotID<-factor(com.meta.rh$PlotID)
com.meta.rh$Irrigation<-factor(com.meta.rh$Irrigation)
com.meta.rh$Block<-factor(com.meta.rh$Block)
com.meta.rh$x<-factor(com.meta.rh$x)
com.meta.rh$rep<-factor(com.meta.rh$rep)
com.meta.rh$rep2<-factor(com.meta.rh$rep2)
com.meta.rh$var2<-factor(com.meta.rh$var2)
com.meta.rh$var3<-factor(com.meta.rh$var3)
com.meta.rh$DxD<-factor(com.meta.rh$DxD)
com.meta.rh$period<-factor(com.meta.rh$period)
#com.meta.rh$period2<-factor(com.meta.rh$period2)

com.meta.rh.edit <- com.meta.rh
com.meta.rh.edit$Date <- factor(com.meta.rh.edit$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th"),
                          labels = c("Apr", "Jun", "Jul"))
com.meta.rh.edit <- com.meta.rh.edit %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN (D)", 
           `K` ="CONFYM (K)", 
           `M` ="CONMIN (M)")

com.stat_text.RS.rich <- data.frame(Date = 0.5, Richness = 10,Treatment="BIODYN", label="C *\nT *")

COM.RS.Richness.plot <- ggplot(com.meta.rh.edit , aes(x=Date, y=Richness)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control', 'Biodyn-drought', 'Confym-control', 
                             'Confym-drought', 'Conmin-control', 'Conmin-drought'))+
  labs(y="COMA Richness", subtitle = "I")+
  facet_wrap(~ Treatment)+
  scale_y_continuous(limits = c(0, 140))+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        #strip.text = element_blank(),
        strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18),
        axis.title.x =element_blank(),
        plot.subtitle = element_text(size=20, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 #annotate(geom = "text", x = as.Date("2022-07-15"), y = 0, hjust = 0, size = 4, label = "Rewetting", color = "grey25")+
 geom_label(data = com.stat_text.RS.rich,label=com.stat_text.RS.rich$label, hjust=0,colour="black", size=4, fontface="bold")

COM.RS.Richness.plot


#_Bulk Soil Shannon Index___________________________________________________________________________________________________________________

com.Sha.stat_text.BS2 <- data.frame(Date = 0.5, Shannon = 0.36,Treatment="BIODYN", label="C *\nD x C *")

COM.BS.Shannon.plot <- ggplot(com.meta.bulk.edit , aes(x=Date, y=Shannon)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control', 'Biodyn-drought', 'Confym-control', 
                             'Confym-drought', 'Conmin-control', 'Conmin-drought'))+
  labs(y="COMA Shannon", subtitle = "F")+
  facet_wrap(~ Treatment)+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18),
        plot.subtitle = element_text(size=20, face = "bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 #geom_vline(data=filter(aob.meta.df.sub.ed2, Type=="Bulk Soil"),aes(xintercept = as.Date("2022-07-14")), linetype="dashed", colour="darkgrey") +
 geom_vline(xintercept = 3.4, linetype="dashed", colour="darkgrey") +
 #geom_vline(xintercept = as.Date("2022-07-14"), linetype="dashed", colour="darkgrey") +
 #annotate(geom = "text", x = as.Date("2022-07-15"), y = 2, hjust = 0, size = 4, label = "Rewetting", color = "slategray")+
 annotate(geom = "text", x = 3.6, y = 0, hjust = 0, size = 4, label = "Rewetting", color = "gray25")+
 geom_label(data = com.Sha.stat_text.BS2,label=com.Sha.stat_text.BS2$label,hjust=0, colour="black", size=4, fontface="bold")
 
COM.BS.Shannon.plot

#_Rhizosphere Shannon Index_________________________________________________________________________________________________________________________

com.stat_text.RS.sha <- data.frame(Date = 0.5, Shannon = 0.36,Treatment="BIODYN", label="C *\nT **")

COM.RS.Shannon.plot <- ggplot(com.meta.rh.edit , aes(x=Date, y=Shannon)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control', 'Biodyn-drought', 'Confym-control', 
                             'Confym-drought', 'Conmin-control', 'Conmin-drought'))+
  labs(y="COMA Shannon", subtitle="L")+
  facet_wrap(~ Treatment)+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_blank(),
        #strip.text = element_text(size=18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.title.y = element_text(size=18),
        axis.title.x =element_blank(),
        plot.subtitle = element_text(size=20, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 #annotate(geom = "text", x = as.Date("2022-07-15"), y = 0, hjust = 0, size = 4, label = "Rewetting", color = "grey25")+
 geom_label(data = com.stat_text.RS.sha,label=com.stat_text.RS.sha$label, hjust=0,colour="black", size=4, fontface="bold")

COM.RS.Shannon.plot


#______________________________________________________________________________________________________________________________
# Inverse Simpson

# Inverse Simpson: plotting the significance across treatment
com.invsimp.pwc.plot <- ggplot(com.meta.df, aes(x=Irrigation, y=InvSimpson)) +
  geom_boxplot(aes(fill = Treatment))+
  theme_bw() +
  labs(y="Comammox Inverse Simpson")+
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
com.invsimp.pwc.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
com.xy.invsimp.bulk <- com.emm.invsimp.bulk %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8) # bulk soil
com.xy.invsimp.rh <- com.emm.invsimp.rh %>% 
  add_xy_position(x = "Irrigation", dodge = 0.8)# rhizosphere
# #combine two data frames and adding 'Type'
com.df.xy.invsimp.bulk <- as.data.frame(com.xy.invsimp.bulk)
com.df.xy.invsimp.rh <- as.data.frame(com.xy.invsimp.rh)
com.df.xy.invsimp.all <- rbind(com.df.xy.invsimp.bulk, com.df.xy.invsimp.rh) 
com.df.xy.invsimp.all$Type <-  c(rep("Bulk Soil", 30), rep("Rhizosphere", 18)) #adding 'Type'
# plotting the pairwise comparisons among treatments (emmeans results)
com.invsimp.pwc.plot2 <- com.invsimp.pwc.plot + 
  stat_pvalue_manual(com.df.xy.invsimp.all,label = "p.adj.signif", size=8, bracket.size = 0.6,bracket.nudge.y = -0.05,bracket.shorten = 0, color = "blue",tip.length = 0.01, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
com.invsimp.pwc.plot2
setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_min_invsimp_all.eps",
       com.min.invsimp.pwc.plot2, device = "eps",
       width = 14, height =5.8, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("COM_invsimp_all.tiff",
       com.invsimp.pwc.plot2, device = "tiff",
       width = 14, height =5.8, 
       units= "in", dpi = 600)

# inverse simpson between irrigations
com.invsimp.pwc.irri.plot <- ggplot(com.meta.df, aes(x=Date, y=InvSimpson)) +
  geom_boxplot(aes(group = var3, fill = Irrigation))+
  theme_bw() +
  labs(y="Comammox Inverse Simpson")+
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
com.invsimp.pwc.irri.plot

setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_invsimp_irri_boxplot.eps",
       com.invsimp.pwc.irri.plot, device = "eps",
       width = 10, height =5.5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("COM_invsimp_irri_boxplot.tiff",
       com.invsimp.pwc.irri.plot, device = "tiff",
       width = 10, height =5.5, 
       units= "in", dpi = 600)

###################################################################################
# Beta Diversity Analyses on Rarefied Data: COMAMMOX
# SEPARATE BETWEEN BULK SOIL AND RHIZOSPHERE
######################################################################################

# BULK SOIL

# 1. Calculating dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)

# Bray-Curtis - Bulk Soil :
com.rare.asv.df
com.asv.bulk <- com.rare.asv.df[,1:118]
com.asv.bulk1 <- com.asv.bulk[rowSums(com.asv.bulk)>0,]
sort(rowSums(com.asv.bulk1, na.rm = FALSE, dims = 1), decreasing = FALSE)
com.bulk_dist_bc <- vegdist(t(com.asv.bulk1), method = "bray")
# jaccard - Bulk Soil :
com.bulk_dist_jac <- vegdist(t(com.asv.bulk1), binary = TRUE, method = "jaccard")
# Weighted UniFrac (rarefied) - Bulk Soil:
com.physeq_bulk <- subset_samples(com.rare.min.physeq, Type=="BS")
com.physeq_bulk1 <- prune_taxa(taxa_sums(com.physeq_bulk)>0, com.physeq_bulk)
com.physeq_bulk1
sort(taxa_sums(com.physeq_bulk1), decreasing =F) #checking
com.bulk_dist_wUF <- UniFrac(com.physeq_bulk1, weighted=TRUE, normalized = TRUE)
com.bulk_dist_wUF
# Unweighted UniFrac (rarefied) -  Bulk Soil:
com.bulk_dist_uwUF <- UniFrac(com.physeq_bulk1, weighted=FALSE, normalized = TRUE)
com.bulk_dist_uwUF

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis - Bulk Soil:
com.bulk_pcoa_bc <- cmdscale(com.bulk_dist_bc, eig=T)
# Jaccard - Bulk Soil:
com.bulk_pcoa_jac <- cmdscale(com.bulk_dist_jac, eig=T)
# Weighted UniFrac - Bulk Soil:
com.bulk_pcoa_wUF <- cmdscale(com.bulk_dist_wUF, eig=T)
# Unweighted UniFrac - Bulk Soil:
com.bulk_pcoa.uwUF <- cmdscale(com.bulk_dist_uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis - Bulk Soil:
ax1.scores.com.BS <- com.bulk_pcoa_bc$points[,1]
ax2.scores.com.BS <- com.bulk_pcoa_bc$points[,2] 
# Jaccard - Bulk Soil:
ax1.scores.j.bulk <- com.bulk_pcoa_jac$points[,1]
ax2.scores.j.bulk <- com.bulk_pcoa_jac$points[,2]
# Weighted UniFrac - Bulk Soil:
ax1.scores.wUF.bulk <- com.bulk_pcoa_wUF$points[,1]
ax2.scores.wUF.bulk <- com.bulk_pcoa_wUF$points[,2]
# Unweighted UniFrac - Bulk Soil:
ax1.scores.uwUF.bulk <- com.bulk_pcoa.uwUF$points[,1]
ax2.scores.uwUF.bulk <- com.bulk_pcoa.uwUF$points[,2]

# 4. Envfit
env.com.bulk <- com.meta.df[1:118,c(13:19, 22,26:28)]
str(env.com.bulk)
env.com.bulk <- env.com.bulk %>% mutate_at(colnames(env.com.bulk), as.numeric)
# bray-curtis
set.seed(13)
env_fit.com.bc.bulk <- envfit(com.bulk_pcoa_bc, env.com.bulk, na.rm=TRUE)
env_fit.com.bc.bulk
# Jaccard 
set.seed(13)
env_fit.com.jac <- envfit(com.bulk_pcoa_jac, env.com.bulk, na.rm=TRUE)
# Weighted UniFrac
set.seed(13)
env_fit.com.wuF <- envfit(com.bulk_pcoa_wUF, env.com.bulk, na.rm=TRUE)
# UnWeighted UniFrac
set.seed(13)
env_fit.com.uwuF <- envfit(com.bulk_pcoa.uwUF, env.com.bulk, na.rm=TRUE)

# 5. calculate percent variance explained, then add to plot
com.meta.bulk <- com.meta.df[1:118,]
# Bray-curtis - Bulk Soil:
ax1.com.BS <- com.bulk_pcoa_bc$eig[1]/sum(com.bulk_pcoa_bc$eig)
ax2.com.BS <- com.bulk_pcoa_bc$eig[2]/sum(com.bulk_pcoa_bc$eig)
com.map.pcoa.bulk <- cbind(com.meta.bulk,ax1.scores.com.BS,ax2.scores.com.BS)
# Jaccard - Bulk Soil:
ax1.j.bulk <- com.bulk_pcoa_jac$eig[1]/sum(com.bulk_pcoa_jac$eig)
ax2.j.bulk <- com.bulk_pcoa_jac$eig[2]/sum(com.bulk_pcoa_jac$eig)
com.map.pcoa.j.bulk <- cbind(com.meta.bulk,ax1.scores.j.bulk,ax2.scores.j.bulk)
# Weighted UniFrac - Bulk Soil:
ax1.wUF.bulk <- com.bulk_pcoa_wUF$eig[1]/sum(com.bulk_pcoa_wUF$eig)
ax2.wUF.bulk <- com.bulk_pcoa_wUF$eig[2]/sum(com.bulk_pcoa_wUF$eig)
com.map.pcoa.wUF.bulk <- cbind(com.meta.bulk,ax1.scores.wUF.bulk,ax2.scores.wUF.bulk)
# Unweighted UniFrac - Bulk Soil:
ax1.uwUF.bulk <- com.bulk_pcoa.uwUF$eig[1]/sum(com.bulk_pcoa.uwUF$eig)
ax2.uwUF.bulk <- com.bulk_pcoa.uwUF$eig[2]/sum(com.bulk_pcoa.uwUF$eig)
com.map.pcoa.uwUF.bulk <- cbind(com.meta.bulk,ax1.scores.uwUF.bulk,ax2.scores.uwUF.bulk)

#################################################################################################

# RHIZOSPHERE

# Bray-Curtis - Rhizosphere :
com.asv.rh <- com.rare.asv.df[,119:190]
com.asv.rh1 <- com.asv.rh[rowSums(com.asv.rh)>0,]
sort(rowSums(com.asv.rh1, na.rm = FALSE, dims = 1), decreasing = FALSE)
dim(com.asv.rh1) 
com.rh_dist_bc <- vegdist(t(com.asv.rh1), method = "bray")
# jaccard - Rhizosphere :
com.rh_dist_jac <- vegdist(t(com.asv.rh1), binary = TRUE, method = "jaccard") 
# Weighted UniFrac (rarefied) - Rhizosphere :
com.physeq_rh <- subset_samples(com.rare.min.physeq, Type=="RS")
com.physeq_rh1 <- prune_taxa(taxa_sums(com.physeq_rh)>0, com.physeq_rh)
com.physeq_rh1
sort(taxa_sums(com.physeq_rh1), decreasing =F) #checking
com.rh_dist_wUF <- UniFrac(com.physeq_rh1, weighted=TRUE, normalized = TRUE)
com.rh_dist_wUF
# Unweighted UniFrac (rarefied) -  Rhizosphere :
com.rh_dist_uwUF <- UniFrac(com.physeq_rh1, weighted=FALSE, normalized = TRUE)
com.rh_dist_uwUF

# 2. CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis

# Bray-Curtis - Rhizosphere :
com.rh_pcoa_bc <- cmdscale(com.rh_dist_bc, eig=T)
# Jaccard - Rhizosphere :
com.rh_pcoa_jac <- cmdscale(com.rh_dist_jac, eig=T)
# Weighted UniFrac - Rhizosphere :
com.rh_pcoa_wUF <- cmdscale(com.rh_dist_wUF, eig=T)
# Unweighted UniFrac - Rhizosphere :
com.rh_pcoa.uwUF <- cmdscale(com.rh_dist_uwUF, eig=T)

# 3. scores of PC1 and PC2

# Bray-Curtis - Rhizosphere :
ax1.scores.com.RS <- com.rh_pcoa_bc$points[,1]
ax2.scores.com.RS <- com.rh_pcoa_bc$points[,2] 
# Jaccard - Rhizosphere :
ax1.scores.j.rh <- com.rh_pcoa_jac$points[,1]
ax2.scores.j.rh <- com.rh_pcoa_jac$points[,2]
# Weighted UniFrac - Rhizosphere :
ax1.scores.wUF.rh <- com.rh_pcoa_wUF$points[,1]
ax2.scores.wUF.rh <- com.rh_pcoa_wUF$points[,2]
# Unweighted UniFrac - Rhizosphere :
ax1.scores.uwUF.rh <- com.rh_pcoa.uwUF$points[,1]
ax2.scores.uwUF.rh <- com.rh_pcoa.uwUF$points[,2]

#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# 4. calculate percent variance explained, then add to plot
com.meta.rh <- com.meta.df[119:190,]
# Bray-curtis - Rhizosphere :
ax1.com.RS <- com.rh_pcoa_bc$eig[1]/sum(com.rh_pcoa_bc$eig)
ax2.com.RS <- com.rh_pcoa_bc$eig[2]/sum(com.rh_pcoa_bc$eig)
com.map.pcoa.rh <- cbind(com.meta.rh,ax1.scores.com.RS,ax2.scores.com.RS)
# Jaccard - Rhizosphere :
ax1.j.rh <- com.rh_pcoa_jac$eig[1]/sum(com.rh_pcoa_jac$eig)
ax2.j.rh <- com.rh_pcoa_jac$eig[2]/sum(com.rh_pcoa_jac$eig)
com.map.pcoa.j.rh <- cbind(com.meta.rh,ax1.scores.j.rh,ax2.scores.j.rh)
# Weighted UniFrac - Rhizosphere :
ax1.wUF.rh <- com.rh_pcoa_wUF$eig[1]/sum(com.rh_pcoa_wUF$eig)
ax2.wUF.rh <- com.rh_pcoa_wUF$eig[2]/sum(com.rh_pcoa_wUF$eig)
com.map.pcoa.wUF.rh <- cbind(com.meta.rh,ax1.scores.wUF.rh,ax2.scores.wUF.rh)
# Unweighted UniFrac - Rhizosphere :
ax1.uwUF.rh <- com.rh_pcoa.uwUF$eig[1]/sum(com.rh_pcoa.uwUF$eig)
ax2.uwUF.rh <- com.rh_pcoa.uwUF$eig[2]/sum(com.rh_pcoa.uwUF$eig)
com.map.pcoa.uwUF.rh <- cbind(com.meta.rh,ax1.scores.uwUF.rh,ax2.scores.uwUF.rh)

# 5. PCoA Plot 

#require("ggrepel")
library(ggrepel)
#install.packages("viridis")
library(viridis)

# A. Bray-Curtis - Bulk Soil :
set.seed(33)
A.bc <- as.list(env_fit.com.bc.bulk$vectors) #shortcutting ef$vectors
pvals.bc<-as.data.frame(A.bc$pvals) #creating the dataframe
#environment scores (vectors scaled by R2 values)
env.scores1.bc <- as.data.frame(scores(env_fit.com.bc.bulk, display="vectors"))
env.scores2.bc <- cbind(env.scores1.bc, pvals.bc)
env.scores3.bc <- cbind(env.scores2.bc,Variable=rownames(env.scores2.bc))
env.scores4.bc <- subset(env.scores3.bc,pvals.bc<0.05)
set.seed(33)
mult <-.65

COM.PCoA.BS <- ggplot(data = com.map.pcoa.bulk, aes(x=ax1.scores.com.BS, y=ax2.scores.com.BS, colour=Treatment))+
  theme_bw()+
  geom_point(aes(color = com.map.pcoa.bulk$Treatment, shape = com.map.pcoa.bulk$Irrigation), size = 4, stroke=2) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#E69F00","#E69F00","#009E73","#009E73","#FF618C","#FF618C","#FF618C","#FF618C",
                            "#E69F00","#E69F00","#009E73","#009E73","#009E73","#009E73","#E69F00","#E69F00",
                            "#FF618C","#FF618C","#FF618C","#FF618C","#009E73","#009E73","#E69F00","#E69F00")) +
 #scale_fill_manual(values = c("#E69F00","#009E73","#FF618C","#FF618C",
                             #"#E69F00","#009E73","#009E73","#E69F00",
                             #"#FF618C","#FF618C","#009E73","#E69F00")) +
  #geom_label(show.legend  = F,aes(label = Block))+
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.com.BS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.com.BS,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "E. Comammox")+
  theme(legend.position="none",
        #legend.title = element_blank(),
        #legend.text=element_text(size=12),
        #legend.spacing.x = unit(0.05, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))+
  geom_mark_ellipse(aes(fill = com.map.pcoa.bulk$PlotID,label = com.map.pcoa.bulk$Block),label.fontsize = 20, 
                    expand = 0, linewidth = NA, show.legend = FALSE)
COM.PCoA.BS
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("com.bray.plotid.tiff",
       com.pcoa_bulk.plot, device = "tiff",
       width = 12, height = 8, 
       units= "in", dpi = 600)

setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("com.bray.plotid.tiff",
     com.pcoa_bulk.plot, device = "tiff",
       width = 6, height =5, 
       units= "in", dpi = 600)


COM.PCoA.BS.fert <- ggplot(data = com.map.pcoa.bulk, aes(x=ax1.scores.com.BS, y=ax2.scores.com.BS, colour=Treatment))+
  theme_classic()+
  geom_point(aes(color = com.map.pcoa.bulk$Treatment, shape = com.map.pcoa.bulk$Irrigation), size = 4,stroke=1) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00")) +
  #geom_label(show.legend  = F,aes(label = Block))+
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.com.BS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.com.BS,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "E. Comammox")+
  theme(legend.position="none",
        #legend.title = element_blank(),
        #legend.text=element_text(size=12),
        #legend.spacing.x = unit(0.05, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face="bold"),
        plot.subtitle = element_text(size = 20, face="bold"),
        axis.text=element_text(size=18), 
        axis.title=element_text(size=19))+
  guides(colour=guide_legend(override.aes = list(size=4)))+
  stat_ellipse(aes(colour = com.map.pcoa.bulk$Treatment))
                    #expand = 0, linewidth = NA, show.legend = FALSE)
COM.PCoA.BS.fert



# B. Bray-Curtis - Rhizosphere :
COM.PCoA.RS <- ggplot(data = com.map.pcoa.rh, aes(x=ax1.scores.com.RS, y=ax2.scores.com.RS, colour=Treatment))+
  theme_bw()+
  geom_point(aes(color = com.map.pcoa.rh$Treatment, shape = com.map.pcoa.rh$Irrigation), size = 4, stroke=2) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  scale_fill_manual(values = c("#E69F00","#E69F00","#009E73","#009E73","#FF618C","#FF618C","#FF618C","#FF618C",
                            "#E69F00","#E69F00","#009E73","#009E73","#009E73","#009E73","#E69F00","#E69F00",
                            "#FF618C","#FF618C","#FF618C","#FF618C","#009E73","#009E73","#E69F00","#E69F00")) +
 #scale_fill_manual(values = c("#E69F00","#009E73","#FF618C","#FF618C",
                             #"#E69F00","#009E73","#009E73","#E69F00",
                             #"#FF618C","#FF618C","#009E73","#E69F00")) +
  #geom_label(show.legend  = F,aes(label = Block))+
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.com.RS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.com.RS,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "F. Comammox")+
  theme(legend.position="none",
        #legend.title = element_blank(),
        #legend.text=element_text(size=12),
        #legend.spacing.x = unit(0.05, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #plot.title = element_text(size = 40, face="bold"),
        plot.subtitle = element_text(size = 35, face="bold"),
        axis.text=element_text(size=25), 
        axis.title=element_text(size=30))+
  guides(colour=guide_legend(override.aes = list(size=4)))+
  geom_mark_ellipse(aes(fill = com.map.pcoa.rh$PlotID,label = com.map.pcoa.rh$Block),label.fontsize = 20, 
                    expand = 0, linewidth = NA, show.legend = FALSE)
COM.PCoA.RS


setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("com.bray.plotid.rh.tiff",
     com.pcoa_rh.plot, device = "tiff",
       width = 6, height =5, 
       units= "in", dpi = 600)

COM.PCoA.RS.fert <- ggplot(data = com.map.pcoa.rh, aes(x=ax1.scores.com.RS, y=ax2.scores.com.RS, colour=Treatment))+
  theme_classic()+
  geom_point(aes(color = com.map.pcoa.rh$Treatment, shape = com.map.pcoa.rh$Irrigation), size = 4,stroke=1) +
  scale_color_manual(values = c("#009E73","#FF618C","#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 1),
                     name = "Irrigation treatment",
                     labels = c("control", "drought")) + theme_classic() +
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00")) +
  #geom_label(show.legend  = F,aes(label = Block))+
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.com.RS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.com.RS,3)*100,"% var. explained", sep=""))+
  labs(subtitle = "F. Comammox")+
  theme(legend.position="none",
        #legend.title = element_blank(),
        #legend.text=element_text(size=12),
        #legend.spacing.x = unit(0.05, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face="bold"),
        plot.subtitle = element_text(size = 20, face="bold"),
        axis.text=element_text(size=18), 
        axis.title=element_text(size=19))+
  guides(colour=guide_legend(override.aes = list(size=4)))+
  stat_ellipse(aes(colour = com.map.pcoa.rh$Treatment))
                    #expand = 0, linewidth = NA, show.legend = FALSE)
COM.PCoA.RS.fert



#install.packages("patchwork")
library(patchwork)

com.bray.plot.envfit <- com.pcoa_bulk.plot |  com.pcoa_rh.plot
com.bray.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("com.bray.tiff",
       com.bray.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("com.bray.tiff",
       com.bray.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
ggsave("com.bray.envfit.tiff",
       com.bray.plot.envfit, device = "tiff",
       width = 16, height = 6, 
       units= "in", dpi = 600)

# A. Jaccard - Bulk Soil :
com.pcoa_bulk.jac <- ggplot(data = com.map.pcoa.j.bulk, aes(x=ax1.scores.j.bulk, y=ax2.scores.j.bulk, colour=Treatment))+
  theme_bw()+
  geom_point(data = com.map.pcoa.j.bulk, aes(x = ax1.scores.j.bulk, y = ax2.scores.j.bulk, shape=Irrigation),size=5, alpha= 0.8)+
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
com.pcoa_bulk.jac

# B. Jaccard - Rhizosphere :
com.pcoa_rh.jac <- ggplot(data = com.map.pcoa.j.rh, aes(x=ax1.scores.j.rh, y=ax2.scores.j.rh, colour=Treatment))+
  theme_bw()+
  geom_point(data = com.map.pcoa.j.rh, aes(x = ax1.scores.j.rh, y = ax2.scores.j.rh, shape=Irrigation),size=5, alpha= 0.8)+
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
com.pcoa_rh.jac

com.jac.plot <- com.pcoa_bulk.jac |  com.pcoa_rh.jac
com.jac.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("com.jac.tiff",
       com.jac.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("com.jac.tiff",
       com.jac.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)

# A. Weighted UniFrac - Bulk Soil :
set.seed(33)
A.wu <- as.list(env_fit.com.wuF$vectors) #shortcutting ef$vectors
pvals.wu<-as.data.frame(A.wu$pvals) #creating the dataframe
#environment scores (vectors scaled by R2 values)
env.scores1.wu <- as.data.frame(scores(env_fit.com.wuF, display="vectors"))
env.scores2.wu <- cbind(env.scores1.wu, pvals)
env.scores3.wu <- cbind(env.scores2.wu,Variable=rownames(env.scores2.wu))
env.scores4.wu <- subset(env.scores3.wu,pvals<0.05)
set.seed(33)
mult <-.25

com.pcoa_bulk.wUF <- ggplot(data = com.map.pcoa.wUF.bulk, aes(x=ax1.scores.wUF.bulk, y=ax2.scores.wUF.bulk, colour=Treatment))+
  theme_bw()+
  geom_point(data = com.map.pcoa.wUF.bulk, aes(x = ax1.scores.wUF.bulk, y = ax2.scores.wUF.bulk, shape=Irrigation),size=5, alpha= 0.6)+
  scale_color_viridis(discrete = T) +
  scale_x_continuous(name=paste("PCoA1:\n",round(ax1.wUF.bulk,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2:\n",round(ax2.wUF.bulk,3)*100,"% var. explained", sep=""))+
  labs(colour = "Treatment",  title = "A. Bulk Soil")+
  geom_segment(data=env.scores4.wu,
               aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), 
               arrow = arrow(length = unit(0.3, "cm")),
               colour = "grey",inherit.aes = FALSE)+
  geom_text_repel(data = env.scores4.wu,
                  aes(x = mult*Dim1, y = mult*Dim2, label = Variable),
                  size = 5,fontface="bold",
                  position=position_jitter(width=0.03,height=0.001), inherit.aes = FALSE)+
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
com.pcoa_bulk.wUF

# B. Weighted UniFrac - Rhizosphere :
com.pcoa_rh.wUF <- ggplot(data = com.map.pcoa.wUF.rh, aes(x=ax1.scores.wUF.rh, y=ax2.scores.wUF.rh, colour=Treatment))+
  theme_bw()+
  geom_point(data = com.map.pcoa.wUF.rh, aes(x = ax1.scores.wUF.rh, y = ax2.scores.wUF.rh, shape=Irrigation),size=5, alpha= 0.6)+
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
com.pcoa_rh.wUF

com.wUF.plot.envfit <- com.pcoa_bulk.wUF |  com.pcoa_rh.wUF
com.wUF.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("com.wUF.tiff",
       com.wUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("com.wUF.tiff",
       com.wUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
ggsave("com.wUF.envfit.tiff",
       com.wUF.plot.envfit, device = "tiff",
       width = 16, height = 6, 
       units= "in", dpi = 600)


# A. Unweighted UniFrac - Bulk Soil :
com.pcoa_bulk.uwUF <- ggplot(data = com.map.pcoa.uwUF.bulk, aes(x=ax1.scores.uwUF.bulk, y=ax2.scores.uwUF.bulk, colour=Treatment))+
  theme_bw()+
  geom_point(data = com.map.pcoa.uwUF.bulk, aes(x = ax1.scores.uwUF.bulk, y = ax2.scores.uwUF.bulk, shape=Irrigation),size=5, alpha= 0.8)+
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
com.pcoa_bulk.uwUF
#com.pcoa_bulk.uwUF.id <- com.pcoa_bulk.uwUF+geom_text_repel(aes(label = SampleID),size = 3, max.overlaps = Inf)
#com.pcoa_bulk.uwUF.id
setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_UnweightedUniFrac_bulk_id.tiff",
       com.pcoa_bulk.uwUF.id, device = "tiff",
       width = 8, height =6, 
       units= "in", dpi = 600)

# B. Unweighted UniFrac - Rhizosphere :
com.pcoa_rh.uwUF <- ggplot(data = com.map.pcoa.uwUF.rh, aes(x=ax1.scores.uwUF.rh, y=ax2.scores.uwUF.rh, colour=Treatment))+
  theme_bw()+
  geom_point(data = com.map.pcoa.uwUF.rh, aes(x = ax1.scores.uwUF.rh, y = ax2.scores.uwUF.rh, shape=Irrigation),size=5, alpha= 0.8)+
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
com.pcoa_rh.uwUF

com.uwUF.plot <- com.pcoa_bulk.uwUF |  com.pcoa_rh.uwUF
com.uwUF.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("com.uwUF.tiff",
       com.uwUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)
setwd('D:/Fina/INRAE_Project/microservices_fig/COM')
ggsave("com.uwUF.tiff",
       com.uwUF.plot, device = "tiff",
       width = 12, height = 5, 
       units= "in", dpi = 600)

############################################################################################
# PERMANOVA FOR BULK SOIL AND RHIZOSPHERE
############################################################################################

block.com=as.factor(com.meta.bulk$Block)
plot.com=as.factor(com.meta.bulk$PlotID)
TxI.com=as.factor(com.meta.bulk$x)
trt.com=as.factor(com.meta.bulk$Treatment)
irri.com=as.factor(com.meta.bulk$Irrigation)

## Betadisper for treatment
com.trt.mod <- betadisper(com.bulk_dist_bc, trt.com)
com.trt.mod
boxplot(com.trt.mod)
# Null hypothesis of no difference in dispersion between groups
set.seed(13)
#permutation-based test for multivariate homogeneity of group dispersion (variances)
permod.com.bs <- permutest(com.trt.mod, permutations = 999, pairwise = T)
permod.com.bs # there is significant differences in dispersion between groups
# the variances among groups are not homogeneous,
hsd.com.bs <- TukeyHSD(com.trt.mod) #which groups differ in relation to their variances
hsd.com.bs

# 1. Using adonis2 package with defined perm to restrict the permutation 
set.seed(133)
com.adonis.bulk.bc <- adonis2(com.bulk_dist_bc ~ Irrigation*Treatment*Date+Block, #strata=block.com, 
                              data=com.meta.bulk, 
                              permutations = 999) # significant
com.adonis.bulk.bc

# similar with below:
perm1.com = how(within = Within(type="free"), 
                plots = Plots(type = "none"),
                blocks = block.com,
                nperm = 999,
                observed = TRUE)
set.seed(133)
com.adonis.bulk.bc.perm1 <- adonis2(com.bulk_dist_bc ~ Irrigation, data=com.meta.bulk, 
                                    permutations = perm1.com)
com.adonis.bulk.bc.perm1

# another way to use how()

# Since our intent is to focus on the variation among treatments, 
# we need to restrict the permutations so that plots are permuted within each block, but plots are not permuted across blocks.
# these two ways are equivalent:
CTRL.t1.com <- how(within = Within(type = "free"),
                   plots = Plots(type = "none"),
                   blocks = block.com,
                   nperm = 999,
                   observed = TRUE)
# and
CTRL.t2.com <- how(within = Within(type = "free"),
                   plots = Plots(strata = block.com, type = "none"),
                   nperm = 999,
                   observed = TRUE)
#they specify that plots are to be freely permuted within blocks but that blocks are not allowed to permute
set.seed(132222) #132222
com.adonis.bulk.bc.CTRL.t2 <- adonis2(com.bulk_dist_bc ~ Irrigation, data=com.meta.bulk, 
                                      permutations = CTRL.t2.com)
com.adonis.bulk.bc.CTRL.t2

# 2. Using ANOSIM package and define the strata
set.seed(132)
com.bc.anosim <- anosim(com.bulk_dist_bc,
                        grouping = irri.com, permutations = 999, strata = block.com)
summary(com.bc.anosim) # SIGNIFICANT


# test the permanova for farming system
set.seed(133)
com.adonis.bulk <- adonis2(com.bulk_dist_bc ~ Irrigation*Treatment*Date+block.com, data=com.meta.bulk, 
                           permutation=999) # only treatment is significant
com.adonis.bulk
###################################################################################################

# B. Bray-Curtis - Rhizosphere : 
block.com.rh=as.factor(com.meta.rh$Block)
plot.com.rh=as.factor(com.meta.rh$PlotID)
TxI.com.rh=as.factor(com.meta.rh$x)
trt.com.rh=as.factor(com.meta.rh$Treatment)
irri.com.rh=as.factor(com.meta.rh$Irrigation)

# and
CTRL.t2.rh.com <- how(within = Within(type = "free"),
                      plots = Plots(strata = block.com.rh, type = "none"),
                      nperm = 999,
                      observed = TRUE)
#they specify that plots are to be freely permuted within blocks but that blocks are not allowed to permute
set.seed(13)
com.adonis.rh.bc.CTRL.t2 <- adonis2(com.rh_dist_bc ~ Irrigation, data=com.meta.rh, 
                                    permutations = CTRL.t2.rh.com)
com.adonis.rh.bc.CTRL.t2

# 2. Using ANOSIM package and define the strata
set.seed(13)
com.bc.anosim.rh <- anosim(com.rh_dist_wUF,
                        grouping = irri.com.rh, permutations = 999, strata = block.com.rh)
summary(com.bc.anosim.rh) # SIGNIFICANT

set.seed(13)
com.adonis.rh <- adonis2(com.rh_dist_bc ~ Irrigation*Treatment*Date+Block, data=com.meta.rh, 
                         permutation=999) # only treatment is significant
com.adonis.rh







set.seed(13)
com.adonis.rh.irri <- adonis2(com.rh_dist_bc ~ Irrigation, data=com.meta.rh, 
                              permutation=999,
                              method="bray", 
                              strata = NULL) # not significant
com.adonis.rh.irri

set.seed(13)
com.adonis.rh.irri2 <- adonis2(com.rh_dist_bc ~ Irrigation, data=com.meta.rh, 
                               permutation=999,
                               method="bray", 
                               strata = com.meta.rh$Treatment) # not significant
com.adonis.rh.irri2

set.seed(13)
com.adonis.rh.trt <- adonis2(com.rh_dist_bc ~ Treatment, data=com.meta.rh, 
                             permutation=999,
                             method="bray", 
                             strata = NULL) # treatment is significant ( p val = 0.001***)
com.adonis.rh.trt

set.seed(13)
com.adonis.rh.date <- adonis2(com.rh_dist_bc ~ Date, data=com.meta.rh, 
                              permutation=999,
                              method="bray", 
                              strata = NULL) # not significant
com.adonis.rh.date

# A. Jaccard - Bulk Soil : 
set.seed(13)
com.adonis.jac.bulk <- adonis2(com.bulk_dist_jac ~ Irrigation*Treatment*Date, data=com.meta.bulk, 
                               permutation=999,
                               method="jaccard", 
                               strata = NULL)
com.adonis.jac.bulk
# B. Jaccard - Rhizosphere : 
set.seed(13)
com.adonis.jac.rh <- adonis2(com.rh_dist_jac ~ Irrigation*Treatment*Date, data=com.meta.rh, 
                             permutation=999,
                             method="jaccard", 
                             strata = NULL)
com.adonis.jac.rh

# A. Weighted UniFrac - Bulk Soil : 
set.seed(13)
com.adonis.wuF.bulk <- adonis2(com.bulk_dist_wUF ~ Irrigation*Treatment*Date, data=com.meta.bulk, 
                               permutation=999, 
                               strata = NULL)
com.adonis.wuF.bulk

# B. Weighted UniFrac - Rhizosphere : 
set.seed(13)
com.adonis.wuF.rh <- adonis2(com.rh_dist_wUF ~ Irrigation*Treatment*Date, data=com.meta.rh, 
                             permutation=999, 
                             strata = NULL)
com.adonis.wuF.rh

# A. Unweighted UniFrac - Bulk Soil : 
set.seed(13)
com.adonis.uwuF.bulk <- adonis2(com.bulk_dist_uwUF ~ Irrigation*Treatment*Date, data=com.meta.bulk, 
                                permutation=999, 
                                strata = NULL)
com.adonis.uwuF.bulk
# B. Unweighted UniFrac - Rhizosphere : 
set.seed(13)
com.adonis.uwuF.rh <- adonis2(com.rh_dist_uwUF ~ Irrigation*Treatment*Date, data=com.meta.rh, 
                              permutation=999, 
                              strata = NULL)
com.adonis.uwuF.rh

########################################################################################
# Pairwise comparison analyses across treatments and between irrigation within date
########################################################################################
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)


# A. Pairwise Adonis Among Treatment (overall)

# 1. Bray-Curtis - Bulk Soil:
set.seed(13)
com.pw.bulk.trt_bc <- pairwiseAdonis::pairwise.adonis(com.bulk_dist_bc, 
                                                      com.meta.bulk$Treatment,
                                                  p.adjust.m = "BH")
com.pw.bulk.trt_bc # all pairwise comparisons are significant (p val =0.001**)

# 2. weighted UniFrac - Bulk Soil:
set.seed(13)
com.pw.bulk.trt_wUF <- pairwiseAdonis::pairwise.adonis(com.bulk_dist_wUF, 
                                                   com.meta.bulk$Treatment,
                                                   p.adjust.m = "BH")
com.pw.bulk.trt_wUF # all pairwise comparisons are significant (p val =0.001**)

# 3. Unweighted UniFrac - Bulk Soil:
set.seed(13)
com.pw.bulk.trt_uwUF <- pairwiseAdonis::pairwise.adonis(com.bulk_dist_uwUF, 
                                                    com.meta.bulk$Treatment,
                                                    p.adjust.m = "BH")
com.pw.bulk.trt_uwUF # all pairwise comparisons are significant (p val =0.001**)

# B. Pairwise Adonis Among Date

# 1. Bray-Curtis - Bulk Soil:
set.seed(13)
com.pw.bulk.dat_bc <- pairwiseAdonis::pairwise.adonis(com.bulk_dist_bc, 
                                                  com.meta.bulk$Date,
                                                  p.adjust.m = "BH")
com.pw.bulk.dat_bc # none are significant 

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
# COMAMMOX Community Composition
########################################################################################
# Phyloseq object of rarefied data and unrarefied data:
# 1. rarefied data
com.rare.min.physeq
phyloseq::tax_table(com.rare.min.physeq)
# merge taxa by species
com.sp <- tax_glom(com.rare.min.physeq, taxrank = "Species", NArm = F)
com.sp.ra <- transform_sample_counts(com.sp, function(x) x/sum(x))
sample_data(com.sp.ra)
#detach(package:plyr)
com.sp.df <- psmelt(com.sp.ra) %>%
  group_by(var2, Type, Date, Treatment, Irrigation, Clade, Species) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
com.abund.trt.subcla <- psmelt(com.cla.ra) %>%
  group_by(Type,Clade) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)




colours <- ColourPalleteMulti(com.sp.df, "Species")
colours
#install.packages("Polychrome")
library(Polychrome)
# build-in color palette
grad.col2 <- c("#990F0F", "#CC5151", "#E57E7E", "#FFB2B2",
               "#99540F", "#B26F2C", "#CC8E51", "#E5B17E", "#FFD8B2",
              "#6B990F", "#85B22C", "#A3CC51",
              "#C3E57E")    


           
#install.packages("colorBlindness")
library(colorBlindness)
displayAvailablePalette(color="white")
SteppedSequential5Steps
Brown2Blue12Steps
Blue2DarkOrange18Steps
str(com.sp.df)
#install.packages("ggh4x")
library(ggh4x)

com.sp.df$Type <- factor(com.sp.df$Type, levels = c("BS", "RS"),
                         labels = c("Bulk Soil", "Rhizosphere"))
com.sp.df$Treatment <- factor(com.sp.df$Treatment, levels = c("D", "K", "M"),
                              labels = c("Biodynamic", "Conventional", "Mineral"))
#com.sp.df$Sub_Clade[com.sp.df$Mean<0.001] <- "Other (less than 0.1%)"
#com.sp.df$Sub_Clade[is.na(com.sp.df$Sub_Clade)] <- "Other (less than 0.1%)"

legend <- "Comammox Taxa"
library(scales)
library(forcats)
library(dplyr)
df <- com.sp.df
df2 <- select(df, Clade, Species, Mean) %>%
  mutate(Clade=factor(Clade, levels=c("Clade A", "Clade B")),
         Species=fct_reorder(Species, 10*as.integer(Clade)))
df2$Species <- paste(df2$Clade,df2$Species, sep = "-")


# plotting
com.sp.plot <- ggplot(df2, aes(x=interaction(Date, Irrigation), y=Mean, fill=Species)) + 
  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(legend, values=grad.col2)+
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
com.sp.plot

setwd('/Users/arifinabintarti/Documents/France/Figures/COM/')
ggsave("COM_meanRA_barplot2.eps",
       com.sp.plot, device = "eps",
       width = 15, height =6, 
       units= "in", dpi = 600)
#setwd('D:/Fina/INRAE_Project/microservices_fig/COM')

#________________________________________________________________________________________________________________________________-
# Clade Level

com.cla <- tax_glom(com.rare.min.physeq, taxrank = "Clade", NArm = F)
com.cla.ra <- transform_sample_counts(com.cla, function(x) x/sum(x))
com.abund.trt.cla <- psmelt(com.cla.ra) %>%
  group_by(var2, Type, Date, Treatment, Irrigation, Clade) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

com.abund.trt.cla$Type <- factor(com.abund.trt.cla$Type, levels = c("BS", "RS"),
                         labels = c("Bulk Soil", "Rhizosphere"))
com.abund.trt.cla$Treatment <- factor(com.abund.trt.cla$Treatment, levels = c("D", "K", "M"),
                              labels = c("BIODYN", "CONFYM", "CONMIN"))
com.abund.trt.cla$Irrigation <- factor(com.abund.trt.cla$Irrigation, levels = c("Rainout", "Control"),
                         labels = c("Drought", "Control"))
com.abund.trt.cla$Date <- factor(com.abund.trt.cla$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("Apr", "Jun", "Jul5", "Jul20", "Sep"))

display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
colorblindFriendly=T)
brewer.pal(4, "Set2")

set.seed(13)
com.clade.plot <- ggplot(com.abund.trt.cla, aes(x=interaction(Date, Irrigation), y=Mean, fill=Clade)) + 
  geom_bar(aes(), stat="identity", position="fill") + 
  theme_bw()+
  #scale_fill_manual(legend, values=c("#1B9E77","#D95F02","#7570B3","#E7298A"))+
  scale_fill_manual(values=c("#66C2A5","#FC8D62"))+
  facet_nested(~Type+Treatment, 
               #nest_line = element_line(linetype = 1, linewidth = 0.5), 
               scales="free")+
               #resect = unit(5, "cm"))+
  theme(legend.direction = "vertical",legend.position="right") + 
  guides(fill=guide_legend(ncol=1))+
  labs(y= "Mean Relative Abundance", title = "C. Comammox")+
  theme(plot.title = element_text(size = 25, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y =element_text(size=18),
        legend.text=element_text(size = 16),
        legend.title = element_text(size=17, face="bold"),
        panel.grid = element_blank(), 
        panel.background = element_blank(),
        #strip.background = element_blank(),
        strip.background=element_rect(color="grey30", fill="white"),
        strip.text.x = element_text(size = 20),
        #panel.border = element_rect(colour = "black", fill = NA,linewidth= 0.2),
        panel.border=element_rect(color="grey30",fill = NA))+
  scale_y_continuous(expand = c(0,0))+ guides(x="axis_nested")
  #geom_vline(xintercept = 5, linetype="dotted", 
             #color = "blue", linewidth=1.5)
com.clade.plot





















