library(dplyr)
library(vegan)
library(ggplot2)
library(devtools)
#install.packages("microeco")
library(microeco)
library(magrittr)
###devtools::install_github("hannet91/ggcor")
#devtools::install_github("houyunhuang/ggcor")
#.rs.restartR()
#install.packages("igraph")
library(ggcor)


# Preparing the microtable class

# rarefied AOB ASV table of bulk soil
aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS") #subset bulk soil from phyloseq object
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1
aob.asv.tab <- as.data.frame(otu_table(aob.physeq_bulk1))
dim(aob.asv.tab)

# taxonomy table AOB
aob.asv.tax <- as.data.frame(phyloseq::tax_table(aob.physeq_bulk1))
dim(aob.asv.tax)
aob.asv.tax$Group <- 'AOB'#adding a new column with AOB as value
view(aob.asv.tax)

# environmental data
aob.asv.env <- aob.meta.df.sub[1:119,]
colnames(aob.asv.env)
dim(aob.asv.env)
aob.asv.env <- aob.asv.env %>% 
  dplyr::rename("UniqueID"="SampleID",
                "AOB_richness" = "Richness",
                 "AOB_Shannon" = "Shannon")

#extract alpha diversity
aob.alpha <- aob.asv.env[,43:44]
head(aob.alpha)
aob.alpha <- rownames_to_column(aob.alpha, var = "SampleID")


# adding qPCR data
setwd('/Users/arifinabintarti/Documents/France/microservices/')
qpcr.BS <- read.csv("qPCR.BS.csv", row.names = 1)
view(qpcr.BS)
colnames(qpcr.BS)
#qpcr.BS.ed <- qpcr.BS[,c(11:12,14:15,19:20,22:23,27:33)]
qpcr.BS.logdws <- qpcr.BS[,c(29:32,38:39)]
dim(qpcr.BS.logdws)
qpcr.BS.logdws.x <- rownames_to_column(qpcr.BS.logdws, var = "SampleID")
qpcr.BS.logdws.aob <- qpcr.BS.logdws.x %>% filter(SampleID != "S11")# filter out S11 from the metadata
qpcr.BS.logdws.aob.x <- column_to_rownames(qpcr.BS.logdws.aob, var = "SampleID")
dim(qpcr.BS.logdws.aob.x)

# combine
aob.asv.env.ed <- cbind(aob.asv.env,qpcr.BS.logdws.aob.x)
head(aob.asv.env.ed)
dim(aob.asv.env.ed)
str(aob.asv.env.ed)

# rename
aob.asv.env.ed <- aob.asv.env.ed %>% 
  dplyr::rename(
         "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "Total_N_mineral" = "Nmin_tot",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") # change column names
# adding alpha data of AOA 
aoa.alpha.ed <- aoa.alpha %>% filter(SampleID != "S11")
head(aoa.alpha.ed)
aoa.alpha.ed.x <- column_to_rownames(aoa.alpha.ed, var = "SampleID")
aob.asv.env.ed2 <- cbind(aob.asv.env.ed,aoa.alpha.ed.x)
head(aob.asv.env.ed2)

# adding alpha data of COMA 
head(aob.asv.env.ed2)
aob.asv.env.ed2 <- rownames_to_column(aob.asv.env.ed2, var = "SampleID")
head(com.alpha)
#com.alpha.x <- rownames_to_column(com.alpha, var = "SampleID")
aob.asv.env.ed3 <- aob.asv.env.ed2 %>%
 left_join(com.alpha, by="SampleID") 
colnames(aob.asv.env.ed3)
view(aob.asv.env.ed3)

# relocate alpha data of AOB
aob.asv.env.ed4 <- aob.asv.env.ed3 %>% relocate("AOB_richness", .after="AOA_Shannon")
aob.asv.env.ed5 <- aob.asv.env.ed4 %>% relocate("AOB_Shannon", .after="AOB_richness")
head(aob.asv.env.ed5)
dim(aob.asv.env.ed5)
aob.asv.env.ed5 <- column_to_rownames(aob.asv.env.ed5, var = "SampleID")

#aob.asv.env.ed5 <- aob.asv.env.ed5 %>% 
  #dplyr::rename(
         #"AOA_richness" = "AOA richness",
         #"AOA_Shannon" = "AOA Shannon")

#aob.asv.env <- aob.asv.env %>% mutate_at(c('GWC_g_g', 'TS', 'NH4', 'NO3', 'Nmin_tot', 'C_tot', 'N_tot', 'pH', 'K_mgkg', 'Mg_mgkg', 'P_mgkg','AOB_Richness'), as.numeric)

# create  a microtable
aob.microdata <- microtable$new(sample_table = aob.asv.env.ed5, otu_table = aob.asv.tab, tax_table = aob.asv.tax)

# calculate beta diversity
aob.microdata$tidy_dataset()
aob.microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)

# perform mantel test
#install.packages("mice")
library(mice)
aob.t1 <- trans_env$new(dataset = aob.microdata, env_cols = c(15:18,20:21,24,28:30,45:56),complete_na = TRUE)
set.seed(13)
aob.t1$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")

print(aob.t1$res_mantel)
view(aob.t1$data_env)

# extract a part of the results 
aob.x1 <- data.frame(aob.t1$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aob.x1[aob.x1=="All"] <- "AOB"
# rename columns
colnames(aob.x1)  <- c("spec", "env", "r", "p.value","p.adj")

################################################################################
# rarefied AOA ASV table of bulk soil
aoa.physeq_bulk <- subset_samples(aoa.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
aoa.physeq_bulk1 <- prune_taxa(taxa_sums(aoa.physeq_bulk)>0, aoa.physeq_bulk)
aoa.physeq_bulk1
aoa.asv.tab <- as.data.frame(otu_table(aoa.physeq_bulk1))
dim(aoa.asv.tab)

# taxonomy table 
aoa.asv.tax <- as.data.frame(phyloseq::tax_table(aoa.physeq_bulk1))
dim(aoa.asv.tax)
aoa.asv.tax$Group <- 'AOA'#adding a new column with AOA as value

# enviromental data
aoa.asv.env <- aoa.meta.df[1:120,]
colnames(aoa.asv.env)
aoa.asv.env <- aoa.asv.env %>% 
  dplyr::rename("UniqueID"="SampleID",
                "AOA_richness" = "Richness",
                "AOA_Shannon" = "Shannon")

#extract alpha diversity
aoa.alpha <- aoa.asv.env[,43:44]
head(aoa.alpha)
aoa.alpha <- rownames_to_column(aoa.alpha, var = "SampleID")

# adding qPCR data
qpcr.BS <- read.csv("qPCR.BS.csv", row.names = 1)
colnames(qpcr.BS)
#qpcr.BS.ed <- qpcr.BS[,c(11:12,14:15,19:20,22:23,27:33)]
qpcr.BS.logdws <- qpcr.BS[,29:32]
head(qpcr.BS.logdws)

# combine
aoa.asv.env.ed <- cbind(aoa.asv.env,qpcr.BS.logdws)
dim(aoa.asv.env.ed)


aoa.asv.env.ed <- aoa.asv.env.ed %>% 
  dplyr::rename(#"AOA_Richness" = "Richness",
                #"AOA_Shannon" = "Shannon",
                #"AOA_InvSimpson" = "InvSimpson",
                "UniqueID"="SampleID",
                "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "Total_N_mineral" = "Nmin_tot",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") # change column names
head(aoa.asv.env.ed)
str(aoa.asv.env.ed)
#aoa.asv.env <- aoa.asv.env %>% mutate_at(c('GWC_g_g', 'TS', 'NH4', 'NO3', 'Nmin_tot', 'C_tot', 'N_tot', 'pH', 'K_mgkg', 'Mg_mgkg', 'P_mgkg','AOA_Richness'), as.numeric)

# adding alpha data of AOB 
dim(aob.alpha)
head(aob.alpha)
head(aoa.asv.env.ed)
aoa.asv.env.ed <- rownames_to_column(aoa.asv.env.ed, var = "SampleID")
aoa.asv.env.ed2 <- aoa.asv.env.ed %>%
 left_join(aob.alpha, by="SampleID") 
colnames(aoa.asv.env.ed2)

# adding alpha data of COMA 
head(com.alpha)
aoa.asv.env.ed3 <- aoa.asv.env.ed2 %>%
 left_join(com.alpha, by="SampleID") 
colnames(aoa.asv.env.ed3)
View(aoa.asv.env.ed3)

# relocate alpha data of AOA
aoa.asv.env.ed4 <- aoa.asv.env.ed3 %>% relocate("AOA_richness", .after="Coma_B_abundance")
aoa.asv.env.ed5 <- aoa.asv.env.ed4 %>% relocate("AOA_Shannon", .after="AOA_richness")
head(aoa.asv.env.ed5)
dim(aoa.asv.env.ed5)
aoa.asv.env.ed5 <- column_to_rownames(aoa.asv.env.ed5, var = "SampleID")
aoa.asv.env.ed5 <- aoa.asv.env.ed5 %>% 
  dplyr::rename(
         "AOB_richness" = "AOB richness",
         "AOB_Shannon" = "AOB Shannon",
         "Coma_richness" = "Coma richness",
         "Coma_Shannon" = "Coma Shannon")

# create  a microtable
aoa.microdata <- microtable$new(sample_table = aoa.asv.env.ed5, otu_table = aoa.asv.tab, tax_table = aoa.asv.tax)
# create  a microtable for CONTROL
#dim(aoa.cont_env.ed5)
#aoa.cont_env.ed5 <- rownames_to_column(aoa.cont_env.ed5, var = "SampleID")
#aoa.cont_env.ed6 <- aoa.cont_env.ed5 %>% filter(SampleID != "S11")# filter out S11 from the metadata
#view(aoa.cont_env.ed6)
#aoa.cont_env.ed6 <- column_to_rownames(aoa.cont_env.ed6, var = "SampleID")

#aoa.BS.cont_df.ed <- aoa.BS.cont_df[,-6]
#sort(rowSums(aoa.BS.cont_df.ed, na.rm = FALSE, dims = 1), decreasing = F)
#aoa.BS.cont_df.ed1 <- aoa.BS.cont_df.ed[rowSums(aoa.BS.cont_df.ed)>0,]
#head(aoa.BS.cont_df.ed1)
#aoa.BS.cont_df.ed1 <- rownames_to_column(aoa.BS.cont_df.ed1, var = "ASV")
#dim(aoa.BS.cont_df.ed1)
#aoa.BS.cont_df.ed1 <- column_to_rownames(aoa.BS.cont_df.ed1, var = "ASV")

#head(aoa.BS.cont_tax)
#aoa.BS.cont_tax <- rownames_to_column(aoa.BS.cont_tax, var = "ASV")
#aoa.BS.cont_tax2 <- aoa.BS.cont_tax %>% semi_join(aoa.BS.cont_df.ed1, by = "ASV")
#dim(aoa.BS.cont_tax2)
#aoa.BS.cont_tax2 <- column_to_rownames(aoa.BS.cont_tax2, var = "ASV")
#head(aoa.BS.cont_tax2)


# calculate beta diversity
aoa.microdata$tidy_dataset()
aoa.microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)

# perform mantel test
aoa.t1 <- trans_env$new(dataset = aoa.microdata, env_cols = c(14:17,19:20,23,27:29,44:53),complete_na = TRUE)
aoa.t1$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")
aoa.t1$res_mantel

# extract a part of the results 
aoa.x1 <- data.frame(aoa.t1$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
aoa.x1[aoa.x1=="All"] <- "AOA"
# rename columns
colnames(aoa.x1)  <- c("spec", "env", "r", "p.value","p.adj")

################################################################################
# rarefied comammox ASV table of bulk soil
com.physeq_bulk <- subset_samples(com.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
com.physeq_bulk1 <- prune_taxa(taxa_sums(com.physeq_bulk)>0, com.physeq_bulk)
com.physeq_bulk1
com.asv.tab <- as.data.frame(otu_table(com.physeq_bulk1))
head(com.asv.tab)

# taxonomy table 
com.asv.tax <- as.data.frame(phyloseq::tax_table(com.physeq_bulk1))
dim(com.asv.tax)
com.asv.tax$Group <- 'COM'#adding a new column with COM as value

# enviromental data
com.asv.env <- com.meta.df[1:118,]
colnames(com.asv.env)
com.asv.env <- com.asv.env %>% 
  dplyr::rename("UniqueID"="SampleID",
                "Coma_richness" = "Richness",
                "Coma_Shannon" = "Shannon")
dim(com.asv.env)

#extract alpha diversity
com.alpha <- com.asv.env[,43:44]
head(com.alpha)
com.alpha <- rownames_to_column(com.alpha, var = "SampleID")

# adding qPCR data
qpcr.BS <- read.csv("qPCR.BS.csv", row.names = 1)
view(qpcr.BS)
#qpcr.BS.ed <- qpcr.BS[,c(11:12,14:15,19:20,22:23,27:33)]
qpcr.BS.logdws <- qpcr.BS[,29:32]
head(qpcr.BS.logdws)
qpcr.BS.logdws.x <- rownames_to_column(qpcr.BS.logdws, var = "SampleID")
qpcr.BS.logdws.com <- qpcr.BS.logdws.x %>% filter(SampleID != "S26",
                                                  SampleID !="S52")# filter out from the metadata
qpcr.BS.logdws.com.x <- column_to_rownames(qpcr.BS.logdws.com, var = "SampleID")

# combine
com.asv.env.ed <- cbind(com.asv.env,qpcr.BS.logdws.com.x)
head(com.asv.env.ed)
com.asv.env.ed <- com.asv.env.ed %>% 
  dplyr::rename(#"COM_Richness" = "Richness",
                #"COM_Shannon" = "Shannon",
                #"COM_InvSimpson" = "InvSimpson",
                "AOA_abundance" = "AOA_logDWS",
         "AOB_abundance" = "AOB_logDWS",
         "Coma_A_abundance" = "ComA_logDWS",
         "Coma_B_abundance" = "ComB_logDWS",
         "Total_N_mineral" = "Nmin_tot",
         "Total_C" = "C_tot",
         "Total_N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") # change column names
#com.asv.env <- com.asv.env %>% mutate_at(c('GWC_g_g', 'TS', 'NH4', 'NO3', 'Nmin_tot', 'C_tot', 'N_tot', 'pH', 'K_mgkg', 'Mg_mgkg', 'P_mgkg','COM_Richness'), as.numeric)
str(com.asv.env.ed)
dim(com.asv.env.ed)

# adding alpha data of AOA 
dim(aoa.alpha)
head(aoa.alpha)
aoa.alpha.ed2 <- aoa.alpha %>% filter(SampleID != "S26",
                                      SampleID != "S52")
aoa.alpha.ed2.x <- column_to_rownames(aoa.alpha.ed2, var = "SampleID")
com.asv.env.ed2 <- cbind(com.asv.env.ed,aoa.alpha.ed2.x)
head(com.asv.env.ed2)

# adding alpha data of AOB
com.asv.env.ed2 <- rownames_to_column(com.asv.env.ed2, var = "SampleID")
head(aob.alpha)
com.asv.env.ed3 <- com.asv.env.ed2 %>%
 left_join(aob.alpha, by="SampleID") 
colnames(com.asv.env.ed3)

# relocate alpha data of COMA
com.asv.env.ed4 <- com.asv.env.ed3 %>% relocate("Coma_richness", .after="AOB Shannon")
com.asv.env.ed5 <- com.asv.env.ed4 %>% relocate("Coma_Shannon", .after="Coma_richness")
head(com.asv.env.ed5)
dim(com.asv.env.ed5)
com.asv.env.ed5 <- column_to_rownames(com.asv.env.ed5, var = "SampleID")
com.asv.env.ed5 <- com.asv.env.ed5 %>% 
  dplyr::rename(
         "AOB_richness" = "AOB richness",
         "AOB_Shannon" = "AOB Shannon")
# create  a microtable
com.microdata <- microtable$new(sample_table = com.asv.env.ed5, otu_table = com.asv.tab, tax_table = com.asv.tax)
# calculate beta diversity
com.microdata$tidy_dataset()
com.microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test
com.t1 <- trans_env$new(dataset = com.microdata, env_cols = c(14:17,19:20,23,27:29,44:53), complete_na = TRUE)
com.t1$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  permutations=999,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")
com.t1$res_mantel
# extract a part of the results 
com.x1 <- data.frame(com.t1$res_mantel) %>% .[, c(1, 2, 5, 6,7)]
com.x1[com.x1=="All"] <- "comammox"
# rename columns
colnames(com.x1)  <- c("spec", "env", "r", "p.value","p.adj")


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

# plotting
#setwd('D:/Fina/INRAE_Project/microservices/')
setwd('/Users/arifinabintarti/Documents/France/microservices/')
aoa.env <- aoa.t1$data_env
aob.env <- aob.t1$data_env
com.env <- com.t1$data_env
#write.csv(aob.env, file = "aob.env.mantel.csv")
#write.csv(aoa.env, file = "aoa.env.mantel.csv")
#write.csv(com.env, file = "com.env.mantel.csv")

all.env <- read.csv("env.all.mantel.csv", row.names = 1)
all.env.no.alpha <- all.env[,c(1:4,6:11)]
colnames(all.env.no.alpha)

qpcr.BS <- read.csv("qPCR.BS.csv", row.names = 1)
colnames(qpcr.BS)
#qpcr.BS.ed <- qpcr.BS[,c(11:12,14:15,19:20,22:23,27:33)]
qpcr.BS.logdws <- qpcr.BS[,29:32]
all.env.qpcr <- cbind(all.env.no.alpha,qpcr.BS.logdws)

all.env.qpcr$Mg_mgkg <- as.numeric(all.env.qpcr$Mg_mgkg)
all.env.qpcr.ed <- all.env.qpcr %>% 
  dplyr::rename(#"COM_Richness" = "Richness",
                #"COM_Shannon" = "Shannon",
                #"COM_InvSimpson" = "InvSimpson",
                "AOA abundance" = "AOA_logDWS",
         "AOB abundance" = "AOB_logDWS",
         "Coma A abundance" = "ComA_logDWS",
         "Coma B abundance" = "ComB_logDWS",
         "Total C" = "C_tot",
         "Total N" = "N_tot",
         "K" = "K_mgkg",
         "Mg" = "Mg_mgkg",
         "P" = "P_mgkg",
         "GWC" = "GWC_g_g") 
str(all.env.qpcr.ed)
colnames(all.env.qpcr.ed)
head(all.env.qpcr.ed)
# combine all alpha diversity
ao.alpha <- aoa.alpha %>%
 left_join(aob.alpha, by="SampleID") %>%
 left_join(com.alpha, by="SampleID")
ao.alpha <- column_to_rownames(ao.alpha, var="SampleID")
head(ao.alpha)
# combine
all.env.qpcr.ed2 <- cbind(all.env.qpcr.ed,ao.alpha)
View(all.env.qpcr.ed2)
str(all.env.qpcr.ed2)

mantelplot <- quickcor(all.env.qpcr.ed2, type = "upper", cor.test = TRUE,method="spearman",exact=FALSE, use="complete.obs") + #cor.test = TRUE
  geom_square() + 
  geom_mark(sig.thres = 0.05, color = "black", size=2.7) +
  add_link(plot_table, mapping = aes(colour = pd, size = rd),
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
mantelplot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures')

ggsave("mantel.env.tiff",
       mantelplot, device = "tiff",
       width = 20, height = 7, 
       units= "in", dpi = 600)
ggsave("mantel.env.noalpha.tiff",
       mantelplot, device = "tiff",
       width = 10, height = 7, 
       units= "in", dpi = 600)


