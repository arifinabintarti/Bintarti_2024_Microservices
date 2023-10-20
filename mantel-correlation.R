# Preparing the microtable class

# rarefied AOB ASV table of bulk soil
aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS") #subset bulk soil from phyloseq object
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1
aob.asv.tab <- as.data.frame(otu_table(aob.physeq_bulk1))
dim(aob.asv.tab)
# taxonomy table AOB
aob.asv.tax <- as.data.frame(tax_table(aob.physeq_bulk1))
dim(aob.asv.tax)
aob.asv.tax$Group <- 'AOB'#adding a new column with AOB as value
# enviromental data
aob.asv.env <- aob.meta.df.sub[1:119,]
str(aob.asv.env)
aob.asv.env <- aob.asv.env %>% 
  dplyr::rename("AOB_Richness" = "Richness",
         "AOB_Shannon" = "Shannon",
         "AOB_InvSimpson" = "InvSimpson") # change column names
aob.asv.env <- aob.asv.env %>% mutate_at(c('GWC_g_g', 'TS', 'NH4', 'NO3', 'Nmin_tot', 'C_tot', 'N_tot', 'pH', 'K_mgkg', 'Mg_mgkg', 'P_mgkg','AOB_Richness'), as.numeric)
# create  a microtable
aob.microdata <- microtable$new(sample_table = aob.asv.env, otu_table = aob.asv.tab, tax_table = aob.asv.tax)
# calculate beta diversity
aob.microdata$tidy_dataset()
aob.microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test
aob.t1 <- trans_env$new(dataset = aob.microdata, env_cols = c(13:19,22,26:28,41:42,44))
aob.t1$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")

aob.t1$res_mantel
# extract a part of the results 
aob.x1 <- data.frame(aob.t1$res_mantel) %>% .[, c(1, 2, 5, 7)]
aob.x1[aob.x1=="All"] <- "AOB"
# rename columns
colnames(aob.x1)  <- c("spec", "env", "r", "p.value")
################################################################################
# rarefied AOA ASV table of bulk soil
aoa.physeq_bulk <- subset_samples(aoa.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
aoa.physeq_bulk1 <- prune_taxa(taxa_sums(aoa.physeq_bulk)>0, aoa.physeq_bulk)
aoa.physeq_bulk1
aoa.asv.tab <- as.data.frame(otu_table(aoa.physeq_bulk1))
dim(aoa.asv.tab)
# taxonomy table 
aoa.asv.tax <- as.data.frame(tax_table(aoa.physeq_bulk1))
dim(aoa.asv.tax)
aoa.asv.tax$Group <- 'AOA'#adding a new column with AOA as value
# enviromental data
aoa.asv.env <- aoa.meta.df[1:120,]
str(aoa.asv.env)
aoa.asv.env <- aoa.asv.env %>% 
  dplyr::rename("AOA_Richness" = "Richness",
                "AOA_Shannon" = "Shannon",
                "AOA_InvSimpson" = "InvSimpson") # change column names
aoa.asv.env <- aoa.asv.env %>% mutate_at(c('GWC_g_g', 'TS', 'NH4', 'NO3', 'Nmin_tot', 'C_tot', 'N_tot', 'pH', 'K_mgkg', 'Mg_mgkg', 'P_mgkg','AOA_Richness'), as.numeric)
# create  a microtable
aoa.microdata <- microtable$new(sample_table = aoa.asv.env, otu_table = aoa.asv.tab, tax_table = aoa.asv.tax)
# calculate beta diversity
aoa.microdata$tidy_dataset()
aoa.microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test
aoa.t1 <- trans_env$new(dataset = aoa.microdata, env_cols = c(13:19,22,26:28,41:42,44))
aoa.t1$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")

# extract a part of the results 
aoa.x1 <- data.frame(aoa.t1$res_mantel) %>% .[, c(1, 2, 5, 7)]
aoa.x1[aoa.x1=="All"] <- "AOA"
# rename columns
colnames(aoa.x1)  <- c("spec", "env", "r", "p.value")
################################################################################
# rarefied comammox ASV table of bulk soil
com.physeq_bulk <- subset_samples(com.rare.min.physeq, Type=="BS") #subset bulk soil from phyloseq object
com.physeq_bulk1 <- prune_taxa(taxa_sums(com.physeq_bulk)>0, com.physeq_bulk)
com.physeq_bulk1
com.asv.tab <- as.data.frame(otu_table(com.physeq_bulk1))
dim(com.asv.tab)
# taxonomy table 
com.asv.tax <- as.data.frame(tax_table(com.physeq_bulk1))
dim(com.asv.tax)
com.asv.tax$Group <- 'COM'#adding a new column with COM as value
# enviromental data
com.asv.env <- com.meta.df[1:118,]
str(com.asv.env)
com.asv.env <- com.asv.env %>% 
  dplyr::rename("COM_Richness" = "Richness",
                "COM_Shannon" = "Shannon",
                "COM_InvSimpson" = "InvSimpson") # change column names
com.asv.env <- com.asv.env %>% mutate_at(c('GWC_g_g', 'TS', 'NH4', 'NO3', 'Nmin_tot', 'C_tot', 'N_tot', 'pH', 'K_mgkg', 'Mg_mgkg', 'P_mgkg','COM_Richness'), as.numeric)
# create  a microtable
com.microdata <- microtable$new(sample_table = com.asv.env, otu_table = com.asv.tab, tax_table = com.asv.tax)
# calculate beta diversity
com.microdata$tidy_dataset()
com.microdata$cal_betadiv(method = "bray",
                          unifrac = F,
                          binary = F)
# perform mantel test
com.t1 <- trans_env$new(dataset = com.microdata, env_cols = c(13:19,22,26:28,41:42,44))
com.t1$cal_mantel(use_measure = "bray", 
                  partial_mantel = F,
                  add_matrix = NULL,
                  method = "spearman",
                  p_adjust_method = "fdr")
com.t1$res_mantel
# extract a part of the results 
com.x1 <- data.frame(com.t1$res_mantel) %>% .[, c(1, 2, 5, 7)]
com.x1[com.x1=="All"] <- "comammox"
# rename columns
colnames(com.x1)  <- c("spec", "env", "r", "p.value")
###############################################################################
# generate interval data
aob.x1 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
aoa.x1 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
com.x1 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# combine three tables
plot_table <- rbind(aob.x1, aoa.x1,com.x1)
# plotting
all.env <- merge(x = aoa.asv.env, y = aob.t1[ , c("Client", "LO")], by = "Client", all.x=TRUE)




aoa.env <- aoa.t1$data_env
aoa.env <- rownames_to_column(aoa.env, var = "ID")
aob.env <- aob.t1$data_env
aob.env <- rownames_to_column(aob.env, var = "ID")
com.env <- com.t1$data_env
com.env <- rownames_to_column(com.env, var = "ID")

quickcor(aoa.t1$data_env, type = "upper") + 
  geom_square() + 
  #geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  add_link(plot_table, mapping = aes(colour = pd, size = rd),
           diag.label = TRUE) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288"))+
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
  #geom_diag_label() + remove_axis("x")




