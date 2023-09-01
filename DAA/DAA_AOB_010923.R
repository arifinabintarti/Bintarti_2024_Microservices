##################################################################################################
#DIFFERENTIAL ABUNDANCE ANALYSIS TEST
#################################################################################################

# GROUP & SEPARATE PHYLOSEQ OBJECT BY TYPE, DATE AND TREATMENT:

##### 1. BULK SOIL - RAREFIED #####

aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS")
aob.physeq_bulk
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1 # 937 taxa, 119 samples
###############################################################################
# Subset the data set BY date

# Date: 04-28-2022
aob04seq<- subset_samples(aob.physeq_bulk1, Date=="04-28-22")
aob04seq1 <- prune_taxa(taxa_sums(aob04seq)>0, aob04seq)
sort(rowSums(otu_table(aob04seq1), na.rm = FALSE, dims = 1), decreasing = F)
aob04seq1 #393 taxa 23 samples
################################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.01 % relative abundance across all samples
physeq.subset <- aob04seq1
physeq.subset #393 Taxa, 23 Samples
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0001)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 316 taxa, 23 samples remain in the data set after filtering

################################################################################
#Lets generate a prevalence table (number of samples each taxa occurs in) for each taxa.
prevalencedf = apply(X = otu_table(physeq.subset),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(physeq.subset))
prevalencedf[1:10,]
dim(prevalencedf)
# calculate prevalence
ps = physeq.subset
df_tmp <- psmelt(ps)
df_tmp$sample <- 0
df_tmp$sample[df_tmp$Abundance > 0] <- 1 #E: DON'T UNDERSTAND WHY THIS IS DONE
df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$x)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$x)))))
#attention il ya Sample et sample
for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$x)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$x == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$x == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)

# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt

tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75 # 63 taxa
####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
#install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset.75
str(tmp_T3s)
#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           "treatment"= as.character(tmp_T3s@sam_data$x))
# force control as intercept
#a[a == "Control"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(tmp_T3s))
# random effect
z <- as.factor(tmp_T3s@sam_data$SampleID)

# model with pairwise comparison ---------------------------------------------------------------------------------
glmT3s.sum.global = data.frame()
glmT3s.pairwise.global = data.frame()
#fam=genpois(link = "log")

for (i in 1:length(taxa_names(tmp_T3s))) {
  
  OTU = taxa_names(tmp_T3s)[i] 
  
  # response variable
  y = as.vector(tmp_T3s@otu_table[OTU,]@.Data)
  
  tryCatch({
    ### model
    glmT3s <- glmmTMB(y ~ -1+a + (1 | z), family='poisson',offset = o)
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum$cond))
    glmT3s.sum
    glmT3s.sum.global = rbind(glmT3s.sum.global,glmT3s.sum)
    ### multiple comparison
    glmT3s.pairwise = emmeans(glmT3s,pairwise~a)
    # select p value
    glmT3s.pairwise.sum = summary(glmT3s.pairwise)
    glmT3s.pairwise.sum = glmT3s.pairwise.sum[["contrasts"]]
    # extract summary
    tmp_df = glmT3s.pairwise.sum
    # keep only comparisons of interest
    tmp = unlist(strsplit(as.character(tmp_df$contrast)," - "))
    tmp_df[,"a"] <- tmp[seq(1,length(tmp),by=2)]
    tmp_df[,"b"] <- tmp[seq(2,length(tmp),by=2)]
    #tmp_df = tmp_df[grep("Ni",tmp_df$b), ]
    tmp_df = cbind("OTU"=OTU,tmp_df)
    # extract results in data frame
    glmT3s.pairwise.global = rbind(glmT3s.pairwise.global,tmp_df)
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(OTU,y,glmT3s,glmT3s.sum)
}

glmT3s.model.global = glmT3s.sum.global
glmT3s.pairwise.global = glmT3s.pairwise.global
glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")

setwd('D:/Fina/INRAE_Project/microservices/DAA/')
write.csv(glmT3s.pairwise.global, file = "AOB_date1_rare_pw_300823.csv")
aob04_subset_table_75 <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aob04_subset_table_75, file = "AOB_date1_rare_inpt_300823.csv")

################################################################################
### Subset the data set per irrigation-treatment-date (RAREFIED)

# Date: 04-28-2022
# 1. MINERAL
M04seq<- subset_samples(aob.physeq_bulk1, Date=="04-28-22" & Treatment=="M")
M04seq1 <- prune_taxa(taxa_sums(M04seq)>0, M04seq)
sort(rowSums(otu_table(M04seq1), na.rm = FALSE, dims = 1), decreasing = F)
# 2. BIODYNAMIC
D04seq<- subset_samples(aob.physeq_bulk1, Date=="04-28-22" & Treatment=="D")
D04seq1 <- prune_taxa(taxa_sums(D04seq)>0, D04seq)
# 3. CONVENTIONAL
K04seq<- subset_samples(aob.physeq_bulk1, Date=="04-28-22" & Treatment=="K")
K04seq1 <- prune_taxa(taxa_sums(K04seq)>0, K04seq)

# Date: 06-01-2022
# 1. MINERAL
M06seq<- subset_samples(aob.physeq_bulk1, Date=="06-01-22" & Treatment=="M")
M06seq1 <- prune_taxa(taxa_sums(M06seq)>0, M06seq)
# 2. BIODYNAMIC
D06seq<- subset_samples(aob.physeq_bulk1, Date=="06-01-22" & Treatment=="D")
D06seq1 <- prune_taxa(taxa_sums(D06seq)>0, D06seq)
# 3. CONVENTIONAL
K06seq<- subset_samples(aob.physeq_bulk1, Date=="06-01-22" & Treatment=="K")
K06seq1 <- prune_taxa(taxa_sums(K06seq)>0, K06seq)

# Date: 07-05-2022
# 1. MINERAL
M0705seq<- subset_samples(aob.physeq_bulk1, Date=="07-05-22" & Treatment=="M")
M0705seq1 <- prune_taxa(taxa_sums(M0705seq)>0, M0705seq)
# 2. BIODYNAMIC
D0705seq<- subset_samples(aob.physeq_bulk1, Date=="07-05-22" & Treatment=="D")
D0705seq1 <- prune_taxa(taxa_sums(D0705seq)>0, D0705seq)
# 3. CONVENTIONAL
K0705seq<- subset_samples(aob.physeq_bulk1, Date=="07-05-22" & Treatment=="K")
K0705seq1 <- prune_taxa(taxa_sums(K0705seq)>0, K0705seq)

# Date: 07-20-2022
# 1. MINERAL
M0720seq<- subset_samples(aob.physeq_bulk1, Date=="07-20-22" & Treatment=="M")
M0720seq1 <- prune_taxa(taxa_sums(M0720seq)>0, M0720seq)
# 2. BIODYNAMIC
D0720seq<- subset_samples(aob.physeq_bulk1, Date=="07-20-22" & Treatment=="D")
D0720seq1 <- prune_taxa(taxa_sums(D0720seq)>0, D0720seq)
# 3. CONVENTIONAL
K0720seq<- subset_samples(aob.physeq_bulk1, Date=="07-20-22" & Treatment=="K")
K0720seq1 <- prune_taxa(taxa_sums(K0720seq)>0, K0720seq)

# Date: 09-13-2022
# 1. MINERAL
M09seq<- subset_samples(aob.physeq_bulk1, Date=="09-13-22" & Treatment=="M")
M09seq1 <- prune_taxa(taxa_sums(M09seq)>0, M09seq)
# 2. BIODYNAMIC
D09seq<- subset_samples(aob.physeq_bulk1, Date=="09-13-22" & Treatment=="D")
D09seq1 <- prune_taxa(taxa_sums(D09seq)>0, D09seq)
# 3. CONVENTIONAL
K09seq<- subset_samples(aob.physeq_bulk1, Date=="09-13-22" & Treatment=="K")
K09seq1 <- prune_taxa(taxa_sums(K09seq)>0, K09seq)


##### 2. RHIZOSPHERE - RAREFIED #####

aob.physeq_rh <- subset_samples(aob.rare.1282.seq, Type=="RS")
aob.physeq_rh # 1222 taxa 72 samples
aob.physeq_rh1 <- prune_taxa(taxa_sums(aob.physeq_rh)>0, aob.physeq_rh)
aob.physeq_rh1 # 831 taxa

# Date: 04-28-2022
# 1. MINERAL
M04.rh.seq<- subset_samples(aob.physeq_rh1, Date=="04-28-22" & Treatment=="M")
M04.rh.seq1 <- prune_taxa(taxa_sums(M04.rh.seq)>0, M04.rh.seq)
sort(rowSums(otu_table(M04.rh.seq1), na.rm = FALSE, dims = 1), decreasing = F)
# 2. BIODYNAMIC
D04.rh.seq<- subset_samples(aob.physeq_rh1, Date=="04-28-22" & Treatment=="D")
D04.rh.seq1 <- prune_taxa(taxa_sums(D04.rh.seq)>0, D04.rh.seq)
# 3. CONVENTIONAL
K04.rh.seq<- subset_samples(aob.physeq_rh1, Date=="04-28-22" & Treatment=="K")
K04.rh.seq1 <- prune_taxa(taxa_sums(K04.rh.seq)>0, K04.rh.seq)

# Date: 06-01-2022
# 1. MINERAL
M06.rh.seq<- subset_samples(aob.physeq_rh1, Date=="06-01-22" & Treatment=="M")
M06.rh.seq1 <- prune_taxa(taxa_sums(M06.rh.seq)>0, M06.rh.seq)
# 2. BIODYNAMIC
D06.rh.seq<- subset_samples(aob.physeq_rh1, Date=="06-01-22" & Treatment=="D")
D06.rh.seq1 <- prune_taxa(taxa_sums(D06.rh.seq)>0, D06.rh.seq)
# 3. CONVENTIONAL
K06.rh.seq<- subset_samples(aob.physeq_rh1, Date=="06-01-22" & Treatment=="K")
K06.rh.seq1 <- prune_taxa(taxa_sums(K06.rh.seq)>0, K06.rh.seq)

# Date: 07-05-2022
# 1. MINERAL
M0705.rh.seq<- subset_samples(aob.physeq_rh1, Date=="07-05-22" & Treatment=="M")
M0705.rh.seq1 <- prune_taxa(taxa_sums(M0705.rh.seq)>0, M0705.rh.seq)
# 2. BIODYNAMIC
D0705.rh.seq<- subset_samples(aob.physeq_rh1, Date=="07-05-22" & Treatment=="D")
D0705.rh.seq1 <- prune_taxa(taxa_sums(D0705.rh.seq)>0, D0705.rh.seq)
# 3. CONVENTIONAL
K0705.rh.seq<- subset_samples(aob.physeq_rh1, Date=="07-05-22" & Treatment=="K")
K0705.rh.seq1 <- prune_taxa(taxa_sums(K0705.rh.seq)>0, K0705.rh.seq)
################################################################################
###############################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.02 % relative abundance across all samples
physeq.subset <- M04seq1
physeq.subset #219 Taxa, 8 Samples
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0002)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 182 taxa, 8 samples remain in the data set after filtering

################################################################################
#Lets generate a prevalence table (number of samples each taxa occurs in) for each taxa.
prevalencedf = apply(X = otu_table(physeq.subset),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(physeq.subset))
prevalencedf[1:10,]
dim(prevalencedf)
# calculate prevalence
ps = physeq.subset
df_tmp <- psmelt(ps)
df_tmp$sample <- 0
df_tmp$sample[df_tmp$Abundance > 0] <- 1 #E: DON'T UNDERSTAND WHY THIS IS DONE
df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$Irrigation)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$Irrigation)))))
#attention il ya Sample et sample
for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$Irrigation)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$Irrigation == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$Irrigation == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)

# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75 # 32 taxa
####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
#install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset.75
str(tmp_T3s)
#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           "treatment"= as.character(tmp_T3s@sam_data$Irrigation))
# force control as intercept
a[a == "Control"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(tmp_T3s))
# random effect
z <- as.factor(tmp_T3s@sam_data$SampleID)

# model with pairwise comparison ---------------------------------------------------------------------------------
glmT3s.sum.global = data.frame()
glmT3s.pairwise.global = data.frame()
#fam=genpois(link = "log")

for (i in 1:length(taxa_names(tmp_T3s))) {
  
  OTU = taxa_names(tmp_T3s)[i] 
  
  # response variable
  y = as.vector(tmp_T3s@otu_table[OTU,]@.Data)
  
  tryCatch({
    ### model
    glmT3s <- glmmTMB(y ~ a + (1 | z), family='poisson',offset = o)
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum$cond))
    glmT3s.sum
    glmT3s.sum.global = rbind(glmT3s.sum.global,glmT3s.sum)
    ### multiple comparison
    glmT3s.pairwise = emmeans(glmT3s,pairwise~a)
    # select p value
    glmT3s.pairwise.sum = summary(glmT3s.pairwise)
    glmT3s.pairwise.sum = glmT3s.pairwise.sum[["contrasts"]]
    # extract summary
    tmp_df = glmT3s.pairwise.sum
    # keep only comparisons of interest
    tmp = unlist(strsplit(as.character(tmp_df$contrast)," - "))
    tmp_df[,"a"] <- tmp[seq(1,length(tmp),by=2)]
    tmp_df[,"b"] <- tmp[seq(2,length(tmp),by=2)]
    #tmp_df = tmp_df[grep("Ni",tmp_df$b), ]
    tmp_df = cbind("OTU"=OTU,tmp_df)
    # extract results in data frame
    glmT3s.pairwise.global = rbind(glmT3s.pairwise.global,tmp_df)
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(OTU,y,glmT3s,glmT3s.sum)
}

glmT3s.model.global = glmT3s.sum.global
glmT3s.pairwise.global = glmT3s.pairwise.global
glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")

setwd('D:/Fina/INRAE_Project/microservices/DAA/')
write.csv(glmT3s.pairwise.global, file = "AOB_date1_rare_pw_300823.csv")
aob04_subset_table_75 <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aob04_subset_table_75, file = "AOB_date1_rare_inpt_300823.csv")

####################################################################################################################################

##### 1. BULK SOIL - NOT RAREFIED #####

aob.physeq # still contain rhizosphere!, 1338 taxa 192 samples
sort(rowSums(otu_table(aob.physeq), na.rm = FALSE, dims = 1), decreasing = FALSE) # nothing is zero
aob.bulk.rawseq <- subset_samples(aob.physeq, Type=="BS")
aob.bulk.rawseq1 <- prune_taxa(taxa_sums(aob.bulk.rawseq)>0,aob.bulk.rawseq)
aob.bulk.rawseq1 # 1008 taxa

# Date: 04-28-2022

# 1. MINERAL
M04rawseq<- subset_samples(aob.bulk.rawseq1, Date=="04-28-22" & Treatment=="M")
M04rawseq1 <- prune_taxa(taxa_sums(M04rawseq)>0, M04rawseq)
M04_table <- as.data.frame(otu_table(M04rawseq1))
M04_table
cond.aldx <- sample_data(M04rawseq1)$Irrigation
# 2. BIODYNAMIC
D04rawseq<- subset_samples(aob.bulk.rawseq1, Date=="04-28-22" & Treatment=="D")
D04rawseq1 <- prune_taxa(taxa_sums(D04rawseq)>0, D04rawseq)
D04_table <- as.data.frame(otu_table(D04rawseq1))
D04_table
cond.aldx.D <- sample_data(D04rawseq1)$Irrigation
# 3. CONVENTIONAL
K04rawseq<- subset_samples(aob.bulk.rawseq1, Date=="04-28-22" & Treatment=="K")
K04rawseq1 <- prune_taxa(taxa_sums(K04rawseq)>0, K04rawseq)
K04_table <- as.data.frame(otu_table(K04rawseq1))
K04_table
cond.aldx.K <- sample_data(K04rawseq1)$Irrigation

# Date: 06-01-2022

# 1. MINERAL
M06rawseq<- subset_samples(aob.bulk.rawseq1, Date=="06-01-22" & Treatment=="M")
M06rawseq1 <- prune_taxa(taxa_sums(M06rawseq)>0, M06rawseq)
M06_table <- as.data.frame(otu_table(M06rawseq1))
M06_table
cond.aldx.M06 <- sample_data(M06rawseq1)$Irrigation
# 2. BIODYNAMIC
D06rawseq<- subset_samples(aob.bulk.rawseq1, Date=="06-01-22" & Treatment=="D")
D06rawseq1 <- prune_taxa(taxa_sums(D06rawseq)>0, D06rawseq)
D06_table <- as.data.frame(otu_table(D06rawseq1))
D06_table
# 3. CONVENTIONAL
K06rawseq<- subset_samples(aob.bulk.rawseq1, Date=="06-01-22" & Treatment=="K")
K06rawseq1 <- prune_taxa(taxa_sums(K06rawseq)>0, K06rawseq)
K06_table <- as.data.frame(otu_table(K06rawseq1))
K06_table
sample_data(K06rawseq1)$Irrigation

# Date: 07-05-2022

# 1. MINERAL
M0705rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-05-22" & Treatment=="M")
M0705rawseq1 <- prune_taxa(taxa_sums(M0705rawseq)>0, M0705rawseq)
M0705_table <- as.data.frame(otu_table(M0705rawseq1))
M0705_table
cond.aldx.M0705 <- sample_data(M0705rawseq1)$Irrigation
# 2. BIODYNAMIC
D0705rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-05-22" & Treatment=="D")
D0705rawseq1 <- prune_taxa(taxa_sums(D0705rawseq)>0, D0705rawseq)
D0705_table <- as.data.frame(otu_table(D0705rawseq1))
D0705_table
# 3. CONVENTIONAL
K0705rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-05-22" & Treatment=="K")
K0705rawseq1 <- prune_taxa(taxa_sums(K0705rawseq)>0, K0705rawseq)
K0705_table <- as.data.frame(otu_table(K0705rawseq1))
K0705_table
sample_data(K0705rawseq1)$Irrigation

# Date: 07-20-2022

# 1. MINERAL
M0720rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-20-22" & Treatment=="M")
M0720rawseq1 <- prune_taxa(taxa_sums(M0720rawseq)>0, M0720rawseq)
M0720_table <- as.data.frame(otu_table(M0720rawseq1))
M0720_table
cond.aldx.M0720 <- sample_data(M0720rawseq1)$Irrigation
# 2. BIODYNAMIC
D0720rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-20-22" & Treatment=="D")
D0720rawseq1 <- prune_taxa(taxa_sums(D0720rawseq)>0, D0720rawseq)
D0720_table <- as.data.frame(otu_table(D0720rawseq1))
D0720_table
# 3. CONVENTIONAL
K0720rawseq<- subset_samples(aob.bulk.rawseq1, Date=="07-20-22" & Treatment=="K")
K0720rawseq1 <- prune_taxa(taxa_sums(K0720rawseq)>0, K0720rawseq)
K0720_table <- as.data.frame(otu_table(K0720rawseq1))
K0720_table
sample_data(K0720rawseq1)$Irrigation

# Date: 09-13-2022

# 1. MINERAL
M09rawseq<- subset_samples(aob.bulk.rawseq1, Date=="09-13-22" & Treatment=="M")
M09rawseq1 <- prune_taxa(taxa_sums(M09rawseq)>0, M09rawseq)
M09_table <- as.data.frame(otu_table(M09rawseq1))
M09_table
cond.aldx.M09 <- sample_data(M09rawseq1)$Irrigation
# 2. BIODYNAMIC
D09rawseq<- subset_samples(aob.bulk.rawseq1, Date=="09-13-22" & Treatment=="D")
D09rawseq1 <- prune_taxa(taxa_sums(D09rawseq)>0, D09rawseq)
D09_table <- as.data.frame(otu_table(D09rawseq1))
D09_table
# 3. CONVENTIONAL
K09rawseq<- subset_samples(aob.bulk.rawseq1, Date=="09-13-22" & Treatment=="K")
K09rawseq1 <- prune_taxa(taxa_sums(K09rawseq)>0, K09rawseq)
K09_table <- as.data.frame(otu_table(K09rawseq1))
K09_table
sample_data(K09rawseq1)$Irrigation

##### 2. RHIZOSPHERE - NOT RAREFIED #####

aob.rh.rawseq <- subset_samples(aob.physeq, Type=="RS")
aob.rh.rawseq1 <- prune_taxa(taxa_sums(aob.rh.rawseq)>0,aob.rh.rawseq)
aob.rh.rawseq1 # 940 taxa 72 samples

# Date: 04-28-2022

# 1. MINERAL
M04.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="04-28-22" & Treatment=="M")
M04.rh.rawseq1 <- prune_taxa(taxa_sums(M04.rh.rawseq)>0, M04.rh.rawseq)
M04.rh_table <- as.data.frame(otu_table(M04.rh.rawseq1))
M04.rh_table
# 2. BIODYNAMIC
D04.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="04-28-22" & Treatment=="D")
D04.rh.rawseq1 <- prune_taxa(taxa_sums(D04.rh.rawseq)>0, D04.rh.rawseq)
D04.rh_table <- as.data.frame(otu_table(D04.rh.rawseq1))
D04.rh_table
# 3. CONVENTIONAL
K04.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="04-28-22" & Treatment=="K")
K04.rh.rawseq1 <- prune_taxa(taxa_sums(K04.rh.rawseq)>0, K04.rh.rawseq)
K04.rh_table <- as.data.frame(otu_table(K04.rh.rawseq1))
K04.rh_table

# Date: 06-01-2022

# 1. MINERAL
M06.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="06-01-22" & Treatment=="M")
M06.rh.rawseq1 <- prune_taxa(taxa_sums(M06.rh.rawseq)>0, M06.rh.rawseq)
M06.rh_table <- as.data.frame(otu_table(M06.rh.rawseq1))
M06.rh_table
# 2. BIODYNAMIC
D06.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="06-01-22" & Treatment=="D")
D06.rh.rawseq1 <- prune_taxa(taxa_sums(D06.rh.rawseq)>0, D06.rh.rawseq)
D06.rh_table <- as.data.frame(otu_table(D06.rh.rawseq1))
D06.rh_table
# 3. CONVENTIONAL
K06.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="06-01-22" & Treatment=="K")
K06.rh.rawseq1 <- prune_taxa(taxa_sums(K06.rh.rawseq)>0, K06.rh.rawseq)
K06.rh_table <- as.data.frame(otu_table(K06.rh.rawseq1))
K06.rh_table

# Date: 07-05-2022

# 1. MINERAL
M0705.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="07-05-22" & Treatment=="M")
M0705.rh.rawseq1 <- prune_taxa(taxa_sums(M0705.rh.rawseq)>0, M0705.rh.rawseq)
M0705.rh_table <- as.data.frame(otu_table(M0705.rh.rawseq1))
M0705.rh_table
# 2. BIODYNAMIC
D0705.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="07-05-22" & Treatment=="D")
D0705.rh.rawseq1 <- prune_taxa(taxa_sums(D0705.rh.rawseq)>0, D0705.rh.rawseq)
D0705.rh_table <- as.data.frame(otu_table(D0705.rh.rawseq1))
D0705.rh_table
# 3. CONVENTIONAL
K0705.rh.rawseq<- subset_samples(aob.rh.rawseq1, Date=="07-05-22" & Treatment=="K")
K0705.rh.rawseq1 <- prune_taxa(taxa_sums(K0705.rh.rawseq)>0, K0705.rh.rawseq)
K0705.rh_table <- as.data.frame(otu_table(K0705.rh.rawseq1))
K0705.rh_table

###################################################################################################################################
###############################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.02 % relative abundance across all samples
physeq.subset <- M04rawseq1
physeq.subset #245 Taxa, 8 Samples
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0002)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 192 taxa, 8 samples remain in the data set after filtering

################################################################################
#Lets generate a prevalence table (number of samples each taxa occurs in) for each taxa.
prevalencedf = apply(X = otu_table(physeq.subset),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(physeq.subset))
prevalencedf[1:10,]
dim(prevalencedf)
# calculate prevalence
ps = physeq.subset
df_tmp <- psmelt(ps)
df_tmp$sample <- 0
df_tmp$sample[df_tmp$Abundance > 0] <- 1 #E: DON'T UNDERSTAND WHY THIS IS DONE
df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$Irrigation)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$Irrigation)))))
#attention il ya Sample et sample
for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$Irrigation)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$Irrigation == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$Irrigation == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)

# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75 # 36 taxa
####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
#install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset.75
str(tmp_T3s)
#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           "treatment"= as.character(tmp_T3s@sam_data$Irrigation))
# force control as intercept
a[a == "Control"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(tmp_T3s))
# random effect
z <- as.factor(tmp_T3s@sam_data$SampleID)

# model with pairwise comparison ---------------------------------------------------------------------------------
glmT3s.sum.global = data.frame()
glmT3s.pairwise.global = data.frame()
#fam=genpois(link = "log")

for (i in 1:length(taxa_names(tmp_T3s))) {
  
  OTU = taxa_names(tmp_T3s)[i] 
  
  # response variable
  y = as.vector(tmp_T3s@otu_table[OTU,]@.Data)
  
  tryCatch({
    ### model
    glmT3s <- glmmTMB(y ~ a + (1 | z), family='poisson',offset = o)
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum$cond))
    glmT3s.sum
    glmT3s.sum.global = rbind(glmT3s.sum.global,glmT3s.sum)
    ### multiple comparison
    glmT3s.pairwise = emmeans(glmT3s,pairwise~a)
    # select p value
    glmT3s.pairwise.sum = summary(glmT3s.pairwise)
    glmT3s.pairwise.sum = glmT3s.pairwise.sum[["contrasts"]]
    # extract summary
    tmp_df = glmT3s.pairwise.sum
    # keep only comparisons of interest
    tmp = unlist(strsplit(as.character(tmp_df$contrast)," - "))
    tmp_df[,"a"] <- tmp[seq(1,length(tmp),by=2)]
    tmp_df[,"b"] <- tmp[seq(2,length(tmp),by=2)]
    #tmp_df = tmp_df[grep("Ni",tmp_df$b), ]
    tmp_df = cbind("OTU"=OTU,tmp_df)
    # extract results in data frame
    glmT3s.pairwise.global = rbind(glmT3s.pairwise.global,tmp_df)
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(OTU,y,glmT3s,glmT3s.sum)
}

glmT3s.model.global = glmT3s.sum.global
glmT3s.pairwise.global = glmT3s.pairwise.global
glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")

setwd('D:/Fina/INRAE_Project/microservices/DAA/')
write.csv(glmT3s.pairwise.global, file = "AOB_date1_rare_pw_300823.csv")
aob04_subset_table_75 <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aob04_subset_table_75, file = "AOB_date1_rare_inpt_300823.csv")

























