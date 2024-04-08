

### Mantel test separated between control and drought

# Project: Microservices, BIODIVERSA
# Author: Ari Fina Bintarti
# Date: 04/03/2024

# Data set: AOB, AOA, Comammox, Soil properties, NO3, NH4
# Note: This analysis was performed separately between control and drought

# 1. CONTROL

# Preparing the microtable class
# rarefied AOB ASV table of bulk soil
aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS") #subset bulk soil from phyloseq object
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1
aob.asv.tab <- as.data.frame(otu_table(aob.physeq_bulk1))
dim(aob.asv.tab)
# subset AOB CONTROL
aob.BS.cont_seq <- subset_samples(aob.physeq_bulk1, Irrigation=="Control")
aob.BS.cont_seq1 <- prune_taxa(taxa_sums(aob.BS.cont_seq)>0, aob.BS.cont_seq)
sort(rowSums(otu_table(aob.BS.cont_seq1), na.rm = FALSE, dims = 1), decreasing = F)
aob.BS.cont_df <- as.data.frame(otu_table(aob.BS.cont_seq1))
dim(aob.BS.cont_df) # 698 otus, 59 samples
view(aob.BS.cont_df)
# taxonomy table AOB CONTROL
aob.BS.cont_tax <- as.data.frame(phyloseq::tax_table(aob.BS.cont_seq1))
dim(aob.BS.cont_tax)
aob.BS.cont_tax$Group <- 'AOB'#adding a new column with AOB as value
view(aob.BS.cont_tax)
# environmental data
aob.asv.env <- aob.meta.df.sub[1:119,]
colnames(aob.asv.env)
dim(aob.asv.env) # 119 46
aob.asv.env <- aob.asv.env %>% 
  dplyr::rename("UniqueID"="SampleID",
                "AOB_richness" = "Richness",
                 "AOB_Shannon" = "Shannon")
# environmental data CONTROL
view(aob.asv.env)
aob.BS.cont_env <- aob.asv.env %>% filter(Irrigation == "Control")
colnames(aob.BS.cont_env)
view(aob.BS.cont_env)
dim(aob.BS.cont_env) # 59 46
#extract alpha diversity CONTROL
aob.BS.cont_alpha <- aob.BS.cont_env[,43:44]
head(aob.BS.cont_alpha)
aob.BS.cont_alpha <- rownames_to_column(aob.BS.cont_alpha, var = "SampleID")
dim(aob.BS.cont_alpha) # 59 3
# adding qPCR data
setwd('/Users/arifinabintarti/Documents/France/microservices/')
qpcr.BS <- read.csv("qPCR.BS.csv", row.names = 1)
# adding qPCR data CONTROL
qpcr.BS.control <- qpcr.BS %>% filter(irrigation == "control")
view(qpcr.BS.control)
colnames(qpcr.BS.control)
qpcr.BS.cont_logdws <- qpcr.BS.control[,c(29:32,38:39)]
dim(qpcr.BS.cont_logdws) # 60 6
qpcr.BS.cont_logdws.x <- rownames_to_column(qpcr.BS.cont_logdws, var = "SampleID")
aob.BS.cont_logdws <- qpcr.BS.cont_logdws.x %>% filter(SampleID != "S11")# filter out S11 from the metadata
aob.BS.cont_logdws <- column_to_rownames(aob.BS.cont_logdws, var = "SampleID")
dim(aob.BS.cont_logdws) # 59 6
view(aob.BS.cont_logdws)
# combine ONLY CONTROL
aob.cont_env.ed <- cbind(aob.BS.cont_env,aob.BS.cont_logdws)
head(aob.cont_env.ed)
dim(aob.cont_env.ed) # 59 52
str(aob.cont_env.ed)
# rename CONTROL environment
aob.cont_env.ed <- aob.cont_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "AOA_AOB_rat" = "AOA_AOB.arc.ratio",
         "ComA_ComB_rat" = "ComA_ComB.arc.ratio",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") # change column names
aob.cont_env.ed$UniqueID <- factor(aob.cont_env.ed$UniqueID)
aob.cont_env.ed$PlotID <- factor(aob.cont_env.ed$PlotID)
aob.cont_env.ed$Block <- factor(aob.cont_env.ed$Block)
aob.cont_env.ed$Irrigation <- factor(aob.cont_env.ed$Irrigation)
aob.cont_env.ed$Treatment <- factor(aob.cont_env.ed$Treatment)
aob.cont_env.ed$Type <- factor(aob.cont_env.ed$Type)
aob.cont_env.ed$rep <- factor(aob.cont_env.ed$rep)
aob.cont_env.ed$rep2 <- factor(aob.cont_env.ed$rep2)
aob.cont_env.ed$Date <- factor(aob.cont_env.ed$Date)
aob.cont_env.ed$x <- factor(aob.cont_env.ed$x)
aob.cont_env.ed$var3 <- factor(aob.cont_env.ed$var3)
aob.cont_env.ed$var2 <- factor(aob.cont_env.ed$var2)
# adding alpha data of AOA CONTROL
dim(aoa.BS.cont_alpha) # 60 2 still contain S11
view(aoa.BS.cont_alpha)
aoa.BS.cont_alpha <- rownames_to_column(aoa.BS.cont_alpha, var = "SampleID")
aoa.BS.cont_alpha.ed <- aoa.BS.cont_alpha %>% filter(SampleID != "S11") # Remove S11 to match the OTU table
dim(aoa.BS.cont_alpha.ed) #59 3
head(aoa.BS.cont_alpha.ed)
aoa.BS.cont_alpha.ed.x <- column_to_rownames(aoa.BS.cont_alpha.ed, var = "SampleID")
aob.cont_env.ed2 <- cbind(aob.cont_env.ed,aoa.BS.cont_alpha.ed.x)
dim(aob.cont_env.ed2) # 59 54
head(aob.cont_env.ed2)
# adding alpha data of COMA CONTROL
aob.cont_env.ed2 <- rownames_to_column(aob.cont_env.ed2, var = "SampleID")
head(com.BS.cont_alpha)
dim(com.BS.cont_alpha) # 60 3 still contain S11
com.BS.cont_alpha.ed <- com.BS.cont_alpha %>% filter(SampleID != "S11") # Remove S11 to match the OTU table
dim(com.BS.cont_alpha.ed) # 59 3
head(com.BS.cont_alpha.ed)
#com.BS.cont_alpha.x <- rownames_to_column(com.BS.cont_alpha, var = "SampleID")
aob.cont_env.ed3 <- aob.cont_env.ed2 %>%
 left_join(com.BS.cont_alpha.ed, by="SampleID") 
dim(aob.cont_env.ed3) #59 57
colnames(aob.cont_env.ed3)
view(aob.cont_env.ed3)
# relocate alpha data of AOB in CONTROL data
aob.cont_env.ed4 <- aob.cont_env.ed3 %>% relocate("AOB_richness", .after="AOA_Shannon")
aob.cont_env.ed5 <- aob.cont_env.ed4 %>% relocate("AOB_Shannon", .after="AOB_richness")
view(aob.cont_env.ed5) # there are nno NAs in the dataset
dim(aob.cont_env.ed5) #59 57
aob.cont_env.ed5 <- column_to_rownames(aob.cont_env.ed5, var = "SampleID")
# create  a microtable CONTROL
aob.cont_microdata <- microtable$new(sample_table = aob.cont_env.ed5, otu_table = aob.BS.cont_df, tax_table = aob.BS.cont_tax)
# calculate beta diversity in CONTROL
aob.cont_microdata$tidy_dataset()
set.seed(13)
aob.cont_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test in CONTROL
aob.t1.cont <- trans_env$new(dataset = aob.cont_microdata, env_cols = c(15:18,20:21,24,28:30,45:56))
set.seed(13)
aob.t1.cont$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")
aob.t1.cont$res_mantel
# extract a part of the results in CONTROL
aob.x1.cont <- data.frame(aob.t1.cont$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aob.x1.cont[aob.x1.cont=="All"] <- "AOB"
# rename columns in CONTROL
colnames(aob.x1.cont)  <- c("spec", "env", "r", "p.value","p.adj")
view(aob.x1.cont)

################################################################################
# rarefied AOA ASV table of bulk soil CONTROL

aoa.physeq_bulk <- subset_samples(aoa.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
aoa.physeq_bulk1 <- prune_taxa(taxa_sums(aoa.physeq_bulk)>0, aoa.physeq_bulk)
aoa.physeq_bulk1
aoa.asv.tab <- as.data.frame(otu_table(aoa.physeq_bulk1))
dim(aoa.asv.tab)
# subset AOA CONTROL
aoa.BS.cont_seq <- subset_samples(aoa.physeq_bulk1, Irrigation=="Control")
aoa.BS.cont_seq1 <- prune_taxa(taxa_sums(aoa.BS.cont_seq)>0, aoa.BS.cont_seq)
sort(rowSums(otu_table(aoa.BS.cont_seq1), na.rm = FALSE, dims = 1), decreasing = F)
aoa.BS.cont_df <- as.data.frame(otu_table(aoa.BS.cont_seq1))
dim(aoa.BS.cont_df) # 379 otus, 60 samples nothing is missing
view(aoa.BS.cont_df)
# taxonomy table AOA CONTROL
aoa.BS.cont_tax <- as.data.frame(phyloseq::tax_table(aoa.BS.cont_seq1))
dim(aoa.BS.cont_tax) # 379 6
aoa.BS.cont_tax$Group <- 'AOA'#adding a new column with AOB as value
view(aoa.BS.cont_tax)
# enviromental data
aoa.asv.env <- aoa.meta.df[1:120,]
colnames(aoa.asv.env)
aoa.asv.env <- aoa.asv.env %>% 
  dplyr::rename("UniqueID"="SampleID",
                "AOA_richness" = "Richness",
                "AOA_Shannon" = "Shannon")
# environmental data CONTROL
view(aoa.asv.env)
aoa.BS.cont_env <- aoa.asv.env %>% filter(Irrigation == "Control")
colnames(aoa.BS.cont_env)
view(aoa.BS.cont_env)
dim(aoa.BS.cont_env) # 60, 46
#extract alpha diversity CONTROL
aoa.BS.cont_alpha <- aoa.BS.cont_env[,43:44]
head(aoa.BS.cont_alpha)
aoa.BS.cont_alpha <- rownames_to_column(aoa.BS.cont_alpha, var = "SampleID")
dim(aoa.BS.cont_alpha) # 60 3
# adding qPCR data CONTROL
qpcr.BS.cont_logdws <- qpcr.BS.control[,c(29:32,38:39)]
head(qpcr.BS.cont_logdws)
aoa.BS.cont_logdws <- qpcr.BS.cont_logdws
dim(aoa.BS.cont_logdws) #60 6
view(aoa.BS.cont_logdws)
# combine ONLY CONTROL
aoa.cont_env.ed <- cbind(aoa.BS.cont_env,aoa.BS.cont_logdws)
head(aoa.cont_env.ed)
dim(aoa.cont_env.ed) # 60 52
str(aoa.cont_env.ed)
# rename CONTROL environment
aoa.cont_env.ed <- aoa.cont_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "AOA_AOB_rat" = "AOA_AOB.arc.ratio",
         "ComA_ComB_rat" = "ComA_ComB.arc.ratio",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") # change column names
aoa.cont_env.ed$UniqueID <- factor(aoa.cont_env.ed$UniqueID)
aoa.cont_env.ed$PlotID <- factor(aoa.cont_env.ed$PlotID)
aoa.cont_env.ed$Block <- factor(aoa.cont_env.ed$Block)
aoa.cont_env.ed$Irrigation <- factor(aoa.cont_env.ed$Irrigation)
aoa.cont_env.ed$Treatment <- factor(aoa.cont_env.ed$Treatment)
aoa.cont_env.ed$Type <- factor(aoa.cont_env.ed$Type)
aoa.cont_env.ed$rep <- factor(aoa.cont_env.ed$rep)
aoa.cont_env.ed$rep2 <- factor(aoa.cont_env.ed$rep2)
aoa.cont_env.ed$Date <- factor(aoa.cont_env.ed$Date)
aoa.cont_env.ed$x <- factor(aoa.cont_env.ed$x)
aoa.cont_env.ed$var3 <- factor(aoa.cont_env.ed$var3)
aoa.cont_env.ed$var2 <- factor(aoa.cont_env.ed$var2)
# adding alpha data of AOB of CONTROL
dim(aob.BS.cont_alpha) # 59 3, it does not contain S11
head(aob.BS.cont_alpha)
head(aoa.cont_env.ed)
aoa.cont_env.ed <- rownames_to_column(aoa.cont_env.ed, var = "SampleID")
aoa.cont_env.ed2 <- aoa.cont_env.ed %>%
 left_join(aob.BS.cont_alpha, by="SampleID") # USE left join for different length of two dataframes
colnames(aoa.cont_env.ed2)
dim(aoa.cont_env.ed2) # 60 55, BUT THERE IS A MISSING DATA ON AOB ALPHA (S11)
view(aoa.cont_env.ed2) 
# adding alpha data of COMA of CONTROL
dim(com.BS.cont_alpha) # 60 3
aoa.cont_env.ed3 <- aoa.cont_env.ed2 %>%
 left_join(com.BS.cont_alpha, by="SampleID") 
dim(aoa.cont_env.ed3) # 60 57
view(aoa.cont_env.ed3)
# relocate alpha data of AOA in CONTROL
colnames(aoa.cont_env.ed3)
aoa.cont_env.ed4 <- aoa.cont_env.ed3 %>% relocate("AOA_richness", .after="ComA_ComB_rat")
aoa.cont_env.ed5 <- aoa.cont_env.ed4 %>% relocate("AOA_Shannon", .after="AOA_richness")
head(aoa.cont_env.ed5)
view(aoa.cont_env.ed5)
dim(aoa.cont_env.ed5) # 60 57, with S11 missing AOB alpha
aoa.cont_env.ed5 <- column_to_rownames(aoa.cont_env.ed5, var = "SampleID")
#aoa.cont_env.ed5$AOB_richness = ifelse(aoa.cont_env.ed5$AOB_richness == "NA", NA, aoa.cont_env.ed5$AOB_richness)
#aoa.cont_env.ed5$AOB_Shannon = ifelse(aoa.cont_env.ed5$AOB_Shannon == "NA", NA, aoa.cont_env.ed5$AOB_Shannon)
# create  a microtable for CONTROL
aoa.cont_microdata <- microtable$new(sample_table = aoa.cont_env.ed5, otu_table = aoa.BS.cont_df, tax_table = aoa.BS.cont_tax)
# calculate beta diversity for CCONTROL
aoa.cont_microdata$tidy_dataset()
set.seed(13)
aoa.cont_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test for CONTROL
set.seed(13)
# set the complete_na=TRUE becuse there is NA, if not it would not working, removal is not an option because then i need to remove everything from the otu table
aoa.t1.cont <- trans_env$new(dataset = aoa.cont_microdata, env_cols = c(15:18,20:21,24,28:30,45:56),complete_na = TRUE) 
set.seed(13)
aoa.t1.cont$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")
aoa.t1.cont$res_mantel
# extract a part of the results CONTROL
aoa.x1.cont <- data.frame(aoa.t1.cont$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aoa.x1.cont[aoa.x1.cont=="All"] <- "AOA"
# rename columns CONTROL
colnames(aoa.x1.cont)  <- c("spec", "env", "r", "p.value","p.adj")
view(aoa.x1.cont)

################################################################################
# rarefied COMAMMOX ASV table of bulk soil CONTROL

com.physeq_bulk <- subset_samples(com.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
com.physeq_bulk1 <- prune_taxa(taxa_sums(com.physeq_bulk)>0, com.physeq_bulk)
com.physeq_bulk1
com.asv.tab <- as.data.frame(otu_table(com.physeq_bulk1))
dim(com.asv.tab) # 497, 118 samples (missing 2 values: SS26 & S52)
# subset AOA CONTROL
com.BS.cont_seq <- subset_samples(com.physeq_bulk1, Irrigation=="Control")
com.BS.cont_seq1 <- prune_taxa(taxa_sums(com.BS.cont_seq)>0, com.BS.cont_seq)
sort(rowSums(otu_table(com.BS.cont_seq1), na.rm = FALSE, dims = 1), decreasing = F)
com.BS.cont_df <- as.data.frame(otu_table(com.BS.cont_seq1))
dim(com.BS.cont_df) # 371 otus, 60 samples (COMPLETE becausse SS26 & S52 are nnot in the CONTROL)
view(com.BS.cont_df)
# taxonomy table COM CONTROL
com.BS.cont_tax <- as.data.frame(phyloseq::tax_table(com.BS.cont_seq1))
dim(com.BS.cont_tax) # 371 4
com.BS.cont_tax$Group <- 'COM'#adding a new column with AOB as value
view(com.BS.cont_tax)
# enviromental data
com.asv.env <- com.meta.df[1:118,]
colnames(com.asv.env)
com.asv.env <- com.asv.env %>% 
  dplyr::rename("UniqueID"="SampleID",
                "Coma_richness" = "Richness",
                "Coma_Shannon" = "Shannon")
dim(com.asv.env)
# environmental data CONTROL
view(com.asv.env)
com.BS.cont_env <- com.asv.env %>% filter(Irrigation == "Control")
colnames(com.BS.cont_env)
view(com.BS.cont_env)
dim(com.BS.cont_env) # 60 46 (complete)
#extract alpha diversity CONTROL
com.BS.cont_alpha <- com.BS.cont_env[,43:44]
head(com.BS.cont_alpha)
com.BS.cont_alpha <- rownames_to_column(com.BS.cont_alpha, var = "SampleID")
dim(com.BS.cont_alpha) # 60 3 (complete)
# adding qPCR data CONTROL
qpcr.BS.cont_logdws <- qpcr.BS.control[,c(29:32,38:39)]
head(qpcr.BS.cont_logdws)
#qpcr.BS.cont_logdws.x <- rownames_to_column(qpcr.BS.cont_logdws, var = "SampleID")
#com.BS.cont_logdws <- qpcr.BS.cont_logdws.x %>% filter(SampleID != "S26",
                                                  #SampleID !="S52")# filter out S11 from the metadata
#com.BS.cont_logdws <- column_to_rownames(com.BS.cont_logdws, var = "SampleID")
com.BS.cont_logdws <- qpcr.BS.cont_logdws
dim(com.BS.cont_logdws) # 60 6 (complete)
view(com.BS.cont_logdws)
# combine ONLY CONTROL
com.cont_env.ed <- cbind(com.BS.cont_env,com.BS.cont_logdws)
head(com.cont_env.ed)
dim(com.cont_env.ed) # 60 52
str(com.cont_env.ed)
# rename CONTROL environment
com.cont_env.ed <- com.cont_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "AOA_AOB_rat" = "AOA_AOB.arc.ratio",
         "ComA_ComB_rat" = "ComA_ComB.arc.ratio",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") # change column names

com.cont_env.ed$UniqueID <- factor(com.cont_env.ed$UniqueID)
com.cont_env.ed$PlotID <- factor(com.cont_env.ed$PlotID)
com.cont_env.ed$Block <- factor(com.cont_env.ed$Block)
com.cont_env.ed$Irrigation <- factor(com.cont_env.ed$Irrigation)
com.cont_env.ed$Treatment <- factor(com.cont_env.ed$Treatment)
com.cont_env.ed$Type <- factor(com.cont_env.ed$Type)
com.cont_env.ed$rep <- factor(com.cont_env.ed$rep)
com.cont_env.ed$rep2 <- factor(com.cont_env.ed$rep2)
com.cont_env.ed$Date <- factor(com.cont_env.ed$Date)
com.cont_env.ed$x <- factor(com.cont_env.ed$x)
com.cont_env.ed$var3 <- factor(com.cont_env.ed$var3)
com.cont_env.ed$var2 <- factor(com.cont_env.ed$var2)

# adding alpha data of AOA CONTROL
head(aoa.BS.cont_alpha)
dim(aoa.BS.cont_alpha) # 60 3 complete
aoa.BS.cont_alpha <- column_to_rownames(aoa.BS.cont_alpha, var = "SampleID")
com.cont_env.ed2 <- cbind(com.cont_env.ed,aoa.BS.cont_alpha)
head(com.cont_env.ed2) # 60 54
# adding alpha data of AOB CONTROL
com.cont_env.ed2 <- rownames_to_column(com.cont_env.ed2, var = "SampleID")
dim(aob.BS.cont_alpha) # 59 3 bbeause it does not contain S11,thus use left_join
com.cont_env.ed3 <- com.cont_env.ed2 %>%
 left_join(aob.BS.cont_alpha, by="SampleID") 
dim(com.cont_env.ed3) # 60 57, but the AOB alpha of S11 are missing
view(com.cont_env.ed3)
# relocate alpha data of COMA CONTROL
com.cont_env.ed4 <- com.cont_env.ed3 %>% relocate("Coma_richness", .after="AOB_Shannon")
com.cont_env.ed5 <- com.cont_env.ed4 %>% relocate("Coma_Shannon", .after="Coma_richness")
head(com.cont_env.ed5)
str(com.cont_env.ed5) # 60 57
com.cont_env.ed5 <- column_to_rownames(com.cont_env.ed5, var = "SampleID")

# create  a microtable CONTROL
com.cont_microdata <- microtable$new(sample_table = com.cont_env.ed5, otu_table = com.BS.cont_df, tax_table = com.BS.cont_tax)
# calculate beta diversity CONTROL
com.cont_microdata$tidy_dataset()
set.seed(13)
com.cont_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test CONTROL
set.seed(13)
# set the complete_na=TRUE becuse there is NA, if not it would not working, removal is not an option because then i need to remove everything from the otu table
com.t1.cont <- trans_env$new(dataset = com.cont_microdata, env_cols = c(15:18,20:21,24,28:30,45:56), complete_na = TRUE)
set.seed(13)
com.t1.cont$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")
com.t1.cont$res_mantel
# extract a part of the results CONTROL
com.x1.cont <- data.frame(com.t1.cont$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
com.x1.cont[com.x1.cont=="All"] <- "comammox"
# rename columns CONTROL
colnames(com.x1.cont)  <- c("spec", "env", "r", "p.value","p.adj")
view(com.x1.cont)

###############################################################################
# generate interval data CONTROL
aob.x1.cont %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
aoa.x1.cont %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
com.x1.cont %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# combine three tables CONTROL
plot_table.CONT <- rbind(aob.x1.cont, aoa.x1.cont,com.x1.cont)
plot_table.CONT <- plot_table.CONT %>% 
  mutate(spec= replace(spec, spec == "AOB", "AOB β diversity")) %>%
  mutate(spec= replace(spec, spec =="AOA","AOA β diversity")) %>%
  mutate(spec= replace(spec, spec=="comammox", "Comammox β diversity")) %>%# change column names
  mutate(env= replace(env, env =="Total_C","Total C")) %>%
  mutate(env= replace(env, env =="Total_N","Total N")) %>%
  mutate(env= replace(env, env =="AOA_abundance","AOA abundance")) %>%
  mutate(env= replace(env, env =="AOB_abundance","AOB abundance")) %>%
  mutate(env= replace(env, env =="Coma_A_abundance","Coma A abundance")) %>%
  mutate(env= replace(env, env =="Coma_B_abundance","Coma B abundance")) %>%
  mutate(env= replace(env, env =="AOA_AOB_rat","AOA/AOB")) %>%
  mutate(env= replace(env, env =="ComA_ComB_rat","Coma A/Coma B")) %>%
  mutate(env= replace(env, env =="AOA_richness","AOA richness")) %>%
  mutate(env= replace(env, env =="AOA_Shannon","AOA Shannon")) %>%
  mutate(env= replace(env, env =="AOB_richness","AOB richness")) %>%
  mutate(env= replace(env, env =="AOB_Shannon","AOB Shannon")) %>%
  mutate(env= replace(env, env =="Coma_richness","Coma richness")) %>%
  mutate(env= replace(env, env =="Coma_Shannon","Coma Shannon")) 
view(plot_table.CONT)
# plotting
#setwd('D:/Fina/INRAE_Project/microservices/')
setwd('/Users/arifinabintarti/Documents/France/microservices/')
aoa.env.cont <- aoa.t1.cont$data_env # the complete one
dim(aoa.env.cont) # 60 22
com.env.cont <- com.t1.cont$data_env # this one is also complete (include all samples and no missing data as a result of interpolation by MICE)
view(com.env.cont)
#write.csv(aob.env, file = "aob.env.mantel.csv")
all.env.cont <- aoa.env.cont
#all.env.no.alpha <- all.env[,c(1:4,6:11)]
colnames(all.env.cont)
str(all.env.cont)
#qpcr.BS <- read.csv("qPCR.BS.csv", row.names = 1)
#colnames(qpcr.BS)
#qpcr.BS.ed <- qpcr.BS[,c(11:12,14:15,19:20,22:23,27:33)]
#qpcr.BS.logdws <- qpcr.BS[,29:32]
#all.env.qpcr <- cbind(all.env.no.alpha,qpcr.BS.logdws)
#all.env.qpcr$Mg_mgkg <- as.numeric(all.env.qpcr$Mg_mgkg)
all.env.cont.ed <- all.env.cont %>% 
  dplyr::rename("AOA abundance" = "AOA_abundance",
         "AOB abundance" = "AOB_abundance",
         "Coma A abundance" = "Coma_A_abundance",
         "Coma B abundance" = "Coma_B_abundance",
         "AOA/AOB" ="AOA_AOB_rat",
         "Coma A/Coma B" ="ComA_ComB_rat",
         "Total C" = "Total_C",
         "Total N" = "Total_N",
         "AOA richness" = "AOA_richness",
         "AOA Shannon" = "AOA_Shannon",
         "AOB richness" = "AOB_richness",
         "AOB Shannon" = "AOB_Shannon",
         "Coma richness" = "Coma_richness",
         "Coma Shannon" = "Coma_Shannon") 
str(all.env.cont.ed)
head(all.env.cont.ed)
# combine all alpha diversity
#ao.alpha <- aoa.alpha %>%
 #left_join(aob.alpha, by="SampleID") %>%
 #left_join(com.alpha, by="SampleID")
#ao.alpha <- column_to_rownames(ao.alpha, var="SampleID")
#head(ao.alpha)
# combine
#all.env.qpcr.ed2 <- cbind(all.env.qpcr.ed,ao.alpha)

mantelplot.control <- quickcor(all.env.cont.ed, type = "upper", cor.test = TRUE,method="spearman",exact=FALSE) + #cor.test = TRUE #use="complete.obs"
  geom_square() + 
  geom_mark(sig.thres = 0.05, color = "black", size=2.7) +
  add_link(plot_table.CONT, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE,
           spec.label.hspace = 1.5, #0.75
           spec.label.vspace = -0.8) + #-0.3
           #spec.label.hspace = 0.75,
           #spec.label.vspace = -0.3) +
  #scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_size_manual(values = c(0.4, 1.1, 2.2)) +
 scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  #scale_size_area(max_size = 3) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))+
  theme(axis.ticks = element_blank())
 #geom_diag_label() +
  #remove_axis("y")
mantelplot.control
setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("mantel.control_everything.tiff",
       mantelplot.control, device = "tiff",
       width = 10, height = 10, 
       units= "in", dpi = 600)

# try with NA of AOB alpha
all.env.cont.ed2 <- rownames_to_column(all.env.cont.ed, var = "SampleID")
all.env.cont.ed3 <- all.env.cont.ed2 %>% 
  mutate(`AOB richness`=ifelse(SampleID=="S11",NA,`AOB richness`),
         `AOB Shannon`=ifelse(SampleID=="S11",NA,`AOB Shannon`))
view(all.env.cont.ed3)
all.env.cont.ed3 <- column_to_rownames(all.env.cont.ed3, var = "SampleID")

mantelplot.control.noAOBalph <- quickcor(all.env.cont.ed3, type = "upper", cor.test = TRUE,method="spearman",exact=FALSE,use="pairwise.complete.obs") + #cor.test = TRUE #use="complete.obs"
  geom_square(position = "identity") + 
  geom_mark(sig.thres = 0.05, color = "black", size=2.7) +
  add_link(plot_table.CONT, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE,
           spec.label.hspace = 1.5, #0.75
           spec.label.vspace = -0.8) + #-0.3
           #spec.label.hspace = 0.75,
           #spec.label.vspace = -0.3) +
  #scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_size_manual(values = c(0.4, 1.1, 2.2)) +
 scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  #scale_size_area(max_size = 3) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))+
  theme(axis.ticks = element_blank(),
        axis.text= element_text(size=15),
        plot.margin = margin(0, 0, 0, 0, "cm"))
 #geom_diag_label() +
  #remove_axis("y")
mantelplot.control.noAOBalph
ggsave("mantelplot.control.noAOBalph.tiff",
       mantelplot.control.noAOBalph, device = "tiff",
       width = 13, height = 10, 
       units= "in", dpi = 600)

###################################################################################################################

# 2. DROUGHT

# Preparing the microtable class
# rarefied AOB ASV table of bulk soil
aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS") #subset bulk soil from phyloseq object
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1
aob.asv.tab <- as.data.frame(otu_table(aob.physeq_bulk1))
dim(aob.asv.tab)
# subset AOB DROUGHT
aob.BS.dro_seq <- subset_samples(aob.physeq_bulk1, Irrigation=="Rainout")
aob.BS.dro_seq1 <- prune_taxa(taxa_sums(aob.BS.dro_seq)>0, aob.BS.dro_seq)
sort(rowSums(otu_table(aob.BS.dro_seq1), na.rm = FALSE, dims = 1), decreasing = F)
aob.BS.dro_df <- as.data.frame(otu_table(aob.BS.dro_seq1))
dim(aob.BS.dro_df) # 669 otus, 60 samples
view(aob.BS.dro_df)
# taxonomy table AOB DROUGHT
aob.BS.dro_tax <- as.data.frame(phyloseq::tax_table(aob.BS.dro_seq1))
dim(aob.BS.dro_tax) # 669 2
aob.BS.dro_tax$Group <- 'AOB'#adding a new column with AOB as value
view(aob.BS.dro_tax)
# environmental data
aob.asv.env <- aob.meta.df.sub[1:119,]
colnames(aob.asv.env)
dim(aob.asv.env) # 119 46
aob.asv.env <- aob.asv.env %>% 
  dplyr::rename("UniqueID"="SampleID",
                "AOB_richness" = "Richness",
                 "AOB_Shannon" = "Shannon")
# environmental data DROUGHT
view(aob.asv.env)
aob.BS.dro_env <- aob.asv.env %>% filter(Irrigation == "Rainout")
colnames(aob.BS.dro_env)
view(aob.BS.dro_env)
dim(aob.BS.dro_env) # 60 46
#extract alpha diversity DROUGHT
aob.BS.dro_alpha <- aob.BS.dro_env[,43:44]
head(aob.BS.dro_alpha)
aob.BS.dro_alpha <- rownames_to_column(aob.BS.dro_alpha, var = "SampleID")
dim(aob.BS.dro_alpha) # 60 3
# adding qPCR data
setwd('/Users/arifinabintarti/Documents/France/microservices/')
qpcr.BS <- read.csv("qPCR.BS.csv", row.names = 1)
# adding qPCR data DROUGHT
qpcr.BS.drought <- qpcr.BS %>% filter(irrigation == "rainout")
view(qpcr.BS.drought)
colnames(qpcr.BS.drought)
qpcr.BS.dro_logdws <- qpcr.BS.drought[,c(29:32,38:39)]
dim(qpcr.BS.dro_logdws) # 60 6
#qpcr.BS.dro_logdws.x <- rownames_to_column(qpcr.BS.dro_logdws, var = "SampleID")
#aob.BS.dro_logdws <- qpcr.BS.dro_logdws.x %>% filter(SampleID != "S11")# filter out S11 from the metadata
#aob.BS.dro_logdws <- column_to_rownames(aob.BS.dro_logdws, var = "SampleID")
aob.BS.dro_logdws <- qpcr.BS.dro_logdws
head(aob.BS.dro_logdws)
# combine ONLY CONTROL
aob.dro_env.ed <- cbind(aob.BS.dro_env,aob.BS.dro_logdws)
head(aob.dro_env.ed)
dim(aob.dro_env.ed) # 60 52
str(aob.dro_env.ed)
# rename CONTROL environment
aob.dro_env.ed <- aob.dro_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "AOA_AOB_rat" = "AOA_AOB.arc.ratio",
         "ComA_ComB_rat" = "ComA_ComB.arc.ratio",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") # change column names
aob.dro_env.ed$UniqueID <- factor(aob.dro_env.ed$UniqueID)
aob.dro_env.ed$PlotID <- factor(aob.dro_env.ed$PlotID)
aob.dro_env.ed$Block <- factor(aob.dro_env.ed$Block)
aob.dro_env.ed$Irrigation <- factor(aob.dro_env.ed$Irrigation)
aob.dro_env.ed$Treatment <- factor(aob.dro_env.ed$Treatment)
aob.dro_env.ed$Type <- factor(aob.dro_env.ed$Type)
aob.dro_env.ed$rep <- factor(aob.dro_env.ed$rep)
aob.dro_env.ed$rep2 <- factor(aob.dro_env.ed$rep2)
aob.dro_env.ed$Date <- factor(aob.dro_env.ed$Date)
aob.dro_env.ed$x <- factor(aob.dro_env.ed$x)
aob.dro_env.ed$var3 <- factor(aob.dro_env.ed$var3)
aob.dro_env.ed$var2 <- factor(aob.dro_env.ed$var2)
# adding alpha data of AOA DROUGHT
dim(aoa.BS.dro_alpha) # 60 3 complete
view(aoa.BS.dro_alpha)
aoa.BS.dro_alpha <- column_to_rownames(aoa.BS.dro_alpha, var = "SampleID")
aob.dro_env.ed2 <- cbind(aob.dro_env.ed,aoa.BS.dro_alpha)
dim(aob.dro_env.ed2) # 60 54
head(aob.dro_env.ed2)
# adding alpha data of COMA DROUGHT
aob.dro_env.ed2 <- rownames_to_column(aob.dro_env.ed2, var = "SampleID")
head(com.BS.dro_alpha)
dim(com.BS.dro_alpha) # 58 3 it does not contain S26 and S52 (missing)
aob.dro_env.ed3 <- aob.dro_env.ed2 %>%
 left_join(com.BS.dro_alpha, by="SampleID") 
dim(aob.dro_env.ed3) # 60 3 but it has NA in the S26 and S52 for Coma alpha
view(aob.dro_env.ed3)
# relocate alpha data of AOB DROUGHT data
aob.dro_env.ed4 <- aob.dro_env.ed3 %>% relocate("AOB_richness", .after="AOA_Shannon")
aob.dro_env.ed5 <- aob.dro_env.ed4 %>% relocate("AOB_Shannon", .after="AOB_richness")
view(aob.dro_env.ed5) # there are nno NAs in the dataset
dim(aob.dro_env.ed5) #60 57 contains missing values (NA)
aob.dro_env.ed5 <- column_to_rownames(aob.dro_env.ed5, var = "SampleID")
# create  a microtable DROUGHT
aob.dro_microdata <- microtable$new(sample_table = aob.dro_env.ed5, otu_table = aob.BS.dro_df, tax_table = aob.BS.dro_tax)
# calculate beta diversity in DROUGHT
aob.dro_microdata$tidy_dataset()
set.seed(13)
aob.dro_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test in DROUGHT
set.seed(13)
aob.t1.dro <- trans_env$new(dataset = aob.dro_microdata, env_cols = c(15:18,20:21,24,28:30,45:56),complete_na = TRUE) # interpolation for the missing values
set.seed(13)
aob.t1.dro$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")
aob.t1.dro$res_mantel
# extract a part of the results in DROUGHT
aob.x1.dro <- data.frame(aob.t1.dro$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aob.x1.dro[aob.x1.dro=="All"] <- "AOB"
# rename columns in DROUGHT
colnames(aob.x1.dro)  <- c("spec", "env", "r", "p.value","p.adj")
view(aob.x1.dro)

################################################################################
# rarefied AOA ASV table of bulk soil DROUGHT

aoa.physeq_bulk <- subset_samples(aoa.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
aoa.physeq_bulk1 <- prune_taxa(taxa_sums(aoa.physeq_bulk)>0, aoa.physeq_bulk)
aoa.physeq_bulk1
aoa.asv.tab <- as.data.frame(otu_table(aoa.physeq_bulk1))
dim(aoa.asv.tab)
# subset AOA DROUGHT
aoa.BS.dro_seq <- subset_samples(aoa.physeq_bulk1, Irrigation=="Rainout")
aoa.BS.dro_seq1 <- prune_taxa(taxa_sums(aoa.BS.dro_seq)>0, aoa.BS.dro_seq)
sort(rowSums(otu_table(aoa.BS.dro_seq1), na.rm = FALSE, dims = 1), decreasing = F)
aoa.BS.dro_df <- as.data.frame(otu_table(aoa.BS.dro_seq1))
dim(aoa.BS.dro_df) # 375 otus, 60 samples nothing is missing
view(aoa.BS.dro_df)
# taxonomy table AOA DROUGHT
aoa.BS.dro_tax <- as.data.frame(phyloseq::tax_table(aoa.BS.dro_seq1))
dim(aoa.BS.dro_tax) # 375 5
aoa.BS.dro_tax$Group <- 'AOA'#adding a new column with AOB as value
view(aoa.BS.dro_tax)
# enviromental data
aoa.asv.env <- aoa.meta.df[1:120,]
colnames(aoa.asv.env)
aoa.asv.env <- aoa.asv.env %>% 
  dplyr::rename("UniqueID"="SampleID",
                "AOA_richness" = "Richness",
                "AOA_Shannon" = "Shannon")
# environmental data DROUGHT
view(aoa.asv.env)
aoa.BS.dro_env <- aoa.asv.env %>% filter(Irrigation == "Rainout")
colnames(aoa.BS.dro_env)
view(aoa.BS.dro_env)
dim(aoa.BS.dro_env) # 60, 46
#extract alpha diversity DROUGHT
aoa.BS.dro_alpha <- aoa.BS.dro_env[,43:44]
head(aoa.BS.dro_alpha)
aoa.BS.dro_alpha <- rownames_to_column(aoa.BS.dro_alpha, var = "SampleID")
dim(aoa.BS.dro_alpha) # 60 3
# adding qPCR data DROUGHT
qpcr.BS.dro_logdws <- qpcr.BS.drought[,c(29:32,38:39)]
head(qpcr.BS.dro_logdws)
aoa.BS.dro_logdws <- qpcr.BS.dro_logdws
dim(aoa.BS.dro_logdws) #60 6
view(aoa.BS.dro_logdws)
# combine ONLY DROUGHT
aoa.dro_env.ed <- cbind(aoa.BS.dro_env,aoa.BS.dro_logdws)
head(aoa.dro_env.ed)
dim(aoa.dro_env.ed) # 60 52
str(aoa.dro_env.ed)
# rename DROUGHT environment
aoa.dro_env.ed <- aoa.dro_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "AOA_AOB_rat" = "AOA_AOB.arc.ratio",
         "ComA_ComB_rat" = "ComA_ComB.arc.ratio",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") # change column names
aoa.dro_env.ed$UniqueID <- factor(aoa.dro_env.ed$UniqueID)
aoa.dro_env.ed$PlotID <- factor(aoa.dro_env.ed$PlotID)
aoa.dro_env.ed$Block <- factor(aoa.dro_env.ed$Block)
aoa.dro_env.ed$Irrigation <- factor(aoa.dro_env.ed$Irrigation)
aoa.dro_env.ed$Treatment <- factor(aoa.dro_env.ed$Treatment)
aoa.dro_env.ed$Type <- factor(aoa.dro_env.ed$Type)
aoa.dro_env.ed$rep <- factor(aoa.dro_env.ed$rep)
aoa.dro_env.ed$rep2 <- factor(aoa.dro_env.ed$rep2)
aoa.dro_env.ed$Date <- factor(aoa.dro_env.ed$Date)
aoa.dro_env.ed$x <- factor(aoa.dro_env.ed$x)
aoa.dro_env.ed$var3 <- factor(aoa.dro_env.ed$var3)
aoa.dro_env.ed$var2 <- factor(aoa.dro_env.ed$var2)
# adding alpha data of AOB of DROUGHT
dim(aob.BS.dro_alpha) # 60 3 complete
head(aob.BS.dro_alpha)
head(aoa.dro_env.ed)
aoa.dro_env.ed <- rownames_to_column(aoa.dro_env.ed, var = "SampleID")
aoa.dro_env.ed2 <- aoa.dro_env.ed %>%
 left_join(aob.BS.dro_alpha, by="SampleID") 
colnames(aoa.dro_env.ed2)
dim(aoa.dro_env.ed2) # 60 55
view(aoa.dro_env.ed2) 
# adding alpha data of COMA of DROUGHT
dim(com.BS.dro_alpha) # 58 3
aoa.dro_env.ed3 <- aoa.dro_env.ed2 %>%
 left_join(com.BS.dro_alpha, by="SampleID") 
dim(aoa.dro_env.ed3) # 60 57 (it has missing values of comammox alpha)
view(aoa.dro_env.ed3)
# relocate alpha data of AOA in DROUGHT
colnames(aoa.dro_env.ed3)
aoa.dro_env.ed4 <- aoa.dro_env.ed3 %>% relocate("AOA_richness", .after="ComA_ComB_rat")
aoa.dro_env.ed5 <- aoa.dro_env.ed4 %>% relocate("AOA_Shannon", .after="AOA_richness")
head(aoa.dro_env.ed5)
view(aoa.dro_env.ed5)
dim(aoa.dro_env.ed5) # 60 57, with S26 & S52 missing Comammox alpha
aoa.dro_env.ed5 <- column_to_rownames(aoa.dro_env.ed5, var = "SampleID")
# create  a microtable for DROUGHT
aoa.dro_microdata <- microtable$new(sample_table = aoa.dro_env.ed5, otu_table = aoa.BS.dro_df, tax_table = aoa.BS.dro_tax)
# calculate beta diversity for DROUGHT
aoa.dro_microdata$tidy_dataset()
set.seed(13)
aoa.dro_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test for DROUGHT
set.seed(13)
# set the complete_na=TRUE becuse there is NA, if not it would not working, removal is not an option because then i need to remove everything from the otu table
aoa.t1.dro <- trans_env$new(dataset = aoa.dro_microdata, env_cols = c(15:18,20:21,24,28:30,45:56),complete_na = TRUE) 
set.seed(13)
aoa.t1.dro$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")
aoa.t1.dro$res_mantel
# extract a part of the results DROUGHT
aoa.x1.dro <- data.frame(aoa.t1.dro$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aoa.x1.dro[aoa.x1.dro=="All"] <- "AOA"
# rename columns DROUGHT
colnames(aoa.x1.dro) <- c("spec", "env", "r", "p.value","p.adj")
view(aoa.x1.dro)

################################################################################
# rarefied COMAMMOX ASV table of bulk soil DROUGHT

com.physeq_bulk <- subset_samples(com.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
com.physeq_bulk1 <- prune_taxa(taxa_sums(com.physeq_bulk)>0, com.physeq_bulk)
com.physeq_bulk1
com.asv.tab <- as.data.frame(otu_table(com.physeq_bulk1))
dim(com.asv.tab) # 497, 118 samples (missing 2 values: SS26 & S52)
# subset COM DROUGHT
com.BS.dro_seq <- subset_samples(com.physeq_bulk1, Irrigation=="Rainout")
com.BS.dro_seq1 <- prune_taxa(taxa_sums(com.BS.dro_seq)>0, com.BS.dro_seq)
sort(rowSums(otu_table(com.BS.dro_seq1), na.rm = FALSE, dims = 1), decreasing = F)
com.BS.dro_df <- as.data.frame(otu_table(com.BS.dro_seq1))
dim(com.BS.dro_df) # 385 otus, 58 samples (missing because S26 & S52 are removed)
view(com.BS.dro_df)
# taxonomy table COM DROUGHT
com.BS.dro_tax <- as.data.frame(phyloseq::tax_table(com.BS.dro_seq1))
dim(com.BS.dro_tax) # 385 4
com.BS.dro_tax$Group <- 'COM'#adding a new column with AOB as value
view(com.BS.dro_tax)
# enviromental data
com.asv.env <- com.meta.df[1:118,]
colnames(com.asv.env)
com.asv.env <- com.asv.env %>% 
  dplyr::rename("UniqueID"="SampleID",
                "Coma_richness" = "Richness",
                "Coma_Shannon" = "Shannon")
dim(com.asv.env)
# environmental data DROUGHT
view(com.asv.env)
com.BS.dro_env <- com.asv.env %>% filter(Irrigation == "Rainout")
colnames(com.BS.dro_env)
view(com.BS.dro_env)
dim(com.BS.dro_env) # 58 46 (S26 and S52 were removed from the metadata to match the otu table)
#extract alpha diversity DROUGHT
com.BS.dro_alpha <- com.BS.dro_env[,43:44]
head(com.BS.dro_alpha)
com.BS.dro_alpha <- rownames_to_column(com.BS.dro_alpha, var = "SampleID")
dim(com.BS.dro_alpha) # 58 3 
# adding qPCR data DROUGHT
qpcr.BS.dro_logdws <- qpcr.BS.drought[,c(29:32,38:39)]
dim(qpcr.BS.dro_logdws) # 60 6 still contain S26 and S52
qpcr.BS.dro_logdws.x <- rownames_to_column(qpcr.BS.dro_logdws, var = "SampleID")
com.BS.dro_logdws <- qpcr.BS.dro_logdws.x %>% filter(SampleID != "S26",
                                                  SampleID !="S52")# filter out S11 from the qPCR data
com.BS.dro_logdws <- column_to_rownames(com.BS.dro_logdws, var = "SampleID")
dim(com.BS.dro_logdws) # 58 6 
view(com.BS.dro_logdws)
# combine ONLY DROUGHT
com.dro_env.ed <- cbind(com.BS.dro_env,com.BS.dro_logdws)
head(com.dro_env.ed)
dim(com.dro_env.ed) # 58 52
str(com.dro_env.ed)
# rename DROUGHT environment
com.dro_env.ed <- com.dro_env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "AOA_AOB_rat" = "AOA_AOB.arc.ratio",
         "ComA_ComB_rat" = "ComA_ComB.arc.ratio",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") # change column names
com.dro_env.ed$UniqueID <- factor(com.dro_env.ed$UniqueID)
com.dro_env.ed$PlotID <- factor(com.dro_env.ed$PlotID)
com.dro_env.ed$Block <- factor(com.dro_env.ed$Block)
com.dro_env.ed$Irrigation <- factor(com.dro_env.ed$Irrigation)
com.dro_env.ed$Treatment <- factor(com.dro_env.ed$Treatment)
com.dro_env.ed$Type <- factor(com.dro_env.ed$Type)
com.dro_env.ed$rep <- factor(com.dro_env.ed$rep)
com.dro_env.ed$rep2 <- factor(com.dro_env.ed$rep2)
com.dro_env.ed$Date <- factor(com.dro_env.ed$Date)
com.dro_env.ed$x <- factor(com.dro_env.ed$x)
com.dro_env.ed$var3 <- factor(com.dro_env.ed$var3)
com.dro_env.ed$var2 <- factor(com.dro_env.ed$var2)
str(com.dro_env.ed)
# adding alpha data of AOA DROUGHT
head(aoa.BS.dro_alpha)
dim(aoa.BS.dro_alpha) # 60 2 complete
head(com.dro_env.ed) # 58 52 
aoa.BS.dro_alpha <- rownames_to_column(aoa.BS.dro_alpha, var = "SampleID")
com.dro_env.ed <- rownames_to_column(com.dro_env.ed, var = "SampleID")
com.dro_env.ed2 <- com.dro_env.ed %>%
 left_join(aoa.BS.dro_alpha, by="SampleID") 
dim(com.dro_env.ed2) #58 55 #the NA are automatically removed to match the comammox meta data
head(com.dro_env.ed2)
# adding alpha data of AOB DROUGHT
head(aob.BS.dro_alpha) # 60 3 
com.dro_env.ed3 <- com.dro_env.ed2 %>%
 left_join(aob.BS.dro_alpha, by="SampleID") 
dim(com.dro_env.ed3) # 58 57, #the NA are automatically removed to match the comammox meta data
view(com.dro_env.ed3)
# relocate alpha data of COMA DROUGHT
com.dro_env.ed4 <- com.dro_env.ed3 %>% relocate("Coma_richness", .after="AOB_Shannon")
com.dro_env.ed5 <- com.dro_env.ed4 %>% relocate("Coma_Shannon", .after="Coma_richness")
head(com.dro_env.ed5)
dim(com.dro_env.ed5) # 58 56
com.dro_env.ed5 <- column_to_rownames(com.dro_env.ed5, var = "SampleID")

# create  a microtable DROUGHT
com.dro_microdata <- microtable$new(sample_table = com.dro_env.ed5, otu_table = com.BS.dro_df, tax_table = com.BS.dro_tax)
# calculate beta diversity DROUGHT
com.dro_microdata$tidy_dataset()
set.seed(13)
com.dro_microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test DROUGHT
set.seed(13)
com.t1.dro <- trans_env$new(dataset = com.dro_microdata, env_cols = c(15:18,20:21,24,28:30,45:56))
set.seed(13)
com.t1.dro$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")
com.t1.dro$res_mantel
# extract a part of the results DROUGHT
com.x1.dro <- data.frame(com.t1.dro$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
com.x1.dro[com.x1.dro=="All"] <- "comammox"
# rename columns DROUGHT
colnames(com.x1.dro)  <- c("spec", "env", "r", "p.value","p.adj")
view(com.x1.dro)

###############################################################################
# generate interval data DROUGHT
aob.x1.dro %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
aoa.x1.dro %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
com.x1.dro %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# combine three tables DROUGHT
plot_table.DRO <- rbind(aob.x1.dro, aoa.x1.dro,com.x1.dro)
plot_table.DRO <- plot_table.DRO %>% 
  mutate(spec= replace(spec, spec == "AOB", "AOB β diversity")) %>%
  mutate(spec= replace(spec, spec =="AOA","AOA β diversity")) %>%
  mutate(spec= replace(spec, spec=="comammox", "Comammox β diversity")) %>%# change column names
  mutate(env= replace(env, env =="Total_C","Total C")) %>%
  mutate(env= replace(env, env =="Total_N","Total N")) %>%
  mutate(env= replace(env, env =="AOA_abundance","AOA abundance")) %>%
  mutate(env= replace(env, env =="AOB_abundance","AOB abundance")) %>%
  mutate(env= replace(env, env =="Coma_A_abundance","Coma A abundance")) %>%
  mutate(env= replace(env, env =="Coma_B_abundance","Coma B abundance")) %>%
  mutate(env= replace(env, env =="AOA_AOB_rat","AOA/AOB")) %>%
  mutate(env= replace(env, env =="ComA_ComB_rat","Coma A/Coma B")) %>%
  mutate(env= replace(env, env =="AOA_richness","AOA richness")) %>%
  mutate(env= replace(env, env =="AOA_Shannon","AOA Shannon")) %>%
  mutate(env= replace(env, env =="AOB_richness","AOB richness")) %>%
  mutate(env= replace(env, env =="AOB_Shannon","AOB Shannon")) %>%
  mutate(env= replace(env, env =="Coma_richness","Coma richness")) %>%
  mutate(env= replace(env, env =="Coma_Shannon","Coma Shannon")) 
view(plot_table.DRO)
# plotting
#setwd('D:/Fina/INRAE_Project/microservices/')
setwd('/Users/arifinabintarti/Documents/France/microservices/')
# extract env data from one of the group
aoa.env.dro <- aoa.t1.dro$data_env # the complete one
dim(aoa.env.dro) # 60 22
aob.env.dro <- aob.t1.dro$data_env 
view(aoa.env.dro)
view(aob.env.dro)
all.env.dro <- aoa.env.dro
colnames(all.env.dro)
str(all.env.dro)
all.env.dro.ed <- all.env.dro %>% 
  dplyr::rename("AOA abundance" = "AOA_abundance",
         "AOB abundance" = "AOB_abundance",
         "Coma A abundance" = "Coma_A_abundance",
         "Coma B abundance" = "Coma_B_abundance",
         "AOA/AOB" ="AOA_AOB_rat",
         "Coma A/Coma B" ="ComA_ComB_rat",
         "Total C" = "Total_C",
         "Total N" = "Total_N",
         "AOA richness" = "AOA_richness",
         "AOA Shannon" = "AOA_Shannon",
         "AOB richness" = "AOB_richness",
         "AOB Shannon" = "AOB_Shannon",
         "Coma richness" = "Coma_richness",
         "Coma Shannon" = "Coma_Shannon") 
str(all.env.dro.ed)
head(all.env.dro.ed)

# plot DROUGHT
mantelplot.drought <- quickcor(all.env.dro.ed, type = "upper", cor.test = TRUE,method="spearman",exact=FALSE) + #cor.test = TRUE #use="complete.obs"
  geom_square() + 
  geom_mark(sig.thres = 0.05, color = "black", size=2.7) +
  add_link(plot_table.DRO, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE,
           spec.label.hspace = 1.5, #0.75
           spec.label.vspace = -0.8) + #-0.3
           #spec.label.hspace = 0.75,
           #spec.label.vspace = -0.3) +
  #scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_size_manual(values = c(0.4, 1.1, 2.2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  #scale_size_area(max_size = 3) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))+
  theme(axis.ticks = element_blank())
 #geom_diag_label() +
  #remove_axis("y")
mantelplot.drought

setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures')
ggsave("mantel.drought_everything.tiff",
       mantelplot.drought, device = "tiff",
       width = 10, height = 10, 
       units= "in", dpi = 600)

# try with NA of COMAMMOX alpha
all.env.dro.ed2 <- rownames_to_column(all.env.dro.ed, var = "SampleID")
all.env.dro.ed3 <- all.env.dro.ed2 %>% 
  mutate(`Coma richness`=ifelse(SampleID=="S26",NA,`Coma richness`),
         `Coma Shannon`=ifelse(SampleID=="S26",NA,`Coma Shannon`),
         `Coma richness`=ifelse(SampleID=="S52",NA,`Coma richness`),
         `Coma Shannon`=ifelse(SampleID=="S52",NA,`Coma Shannon`))
view(all.env.dro.ed3)
all.env.dro.ed3 <- column_to_rownames(all.env.dro.ed3, var = "SampleID")
# plot
#install.packages("latex2exp")
library(latex2exp)
mantelplot.drought.noCOMalph <- quickcor(fixed.xy = F,all.env.dro.ed3, type = "upper", 
                                         cor.test = TRUE,
                                         #cluster = TRUE,
                                         method="spearman",
                                         exact=F,
                                         show.diag = F,
                                         use="pairwise.complete.obs") + #cor.test = TRUE #use="complete.obs"
  geom_square() + 
  geom_mark(sig.thres = 0.05, 
            color = "black", 
            size=2.7)+
            #digits = 2,
            #nudge_x = 0, 
            #nudge_y = 0,) +
  add_link(plot_table.DRO, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE,
           spec.label.hspace = 1.5, #0.75
           spec.label.vspace = -0.8) + #-0.3
           #spec.label.hspace = 0.75,
           #spec.label.vspace = -0.3) +
  #scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_size_manual(values = c(0.4, 1.1, 2.2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  #scale_size_area(max_size = 10) +
  guides(fill = guide_colorbar(title = "Spearman's r",order = 3),
         size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p-value", override.aes = list(size = 3), order = 1))+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=15),
        plot.margin = margin(0, 0, 0, 0, "cm"))
 #geom_diag_label() +
  #remove_axis("y")
mantelplot.drought.noCOMalph

ggsave("mantelplot.drought.noCOMalph2.tiff",
       mantelplot.drought.noCOMalph, device = "tiff",
       width = 12, height = 8.1, 
       units= "in", dpi = 600)
