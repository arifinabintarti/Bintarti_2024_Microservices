##################################################################################################
#DIFFERENTIAL ABUNDANCE ANALYSIS TEST AOA
#################################################################################################

# GROUP & SEPARATE PHYLOSEQ OBJECT BY TYPE, DATE AND TREATMENT:

##### 1. BULK SOIL - RAREFIED #####
# rarefied ASV table # 592 taxa
aoa.rare.min.physeq
# separate the bulk soil
aoa.physeq_bulk <- subset_samples(aoa.rare.min.physeq, Type=="BS")
aoa.physeq_bulk1 <- prune_taxa(taxa_sums(aoa.physeq_bulk)>0, aoa.physeq_bulk)
aoa.physeq_bulk1 # 487 taxa
################################################################################
# Subset the data set per irrigation-treatment-date

# Date: 04-28-2022
# 1. MINERAL
aoaM04seq<- subset_samples(aoa.physeq_bulk1, Date=="04-28-22" & Treatment=="M")
aoaM04seq1 <- prune_taxa(taxa_sums(aoaM04seq)>0, aoaM04seq)
sort(rowSums(otu_table(aoaM04seq1), na.rm = FALSE, dims = 1), decreasing = F)
# 2. BIODYNAMIC
aoaD04seq<- subset_samples(aoa.physeq_bulk1, Date=="04-28-22" & Treatment=="D")
aoaD04seq1 <- prune_taxa(taxa_sums(aoaD04seq)>0, aoaD04seq)
# 3. CONVENTIONAL
aoaK04seq<- subset_samples(aoa.physeq_bulk1, Date=="04-28-22" & Treatment=="K")
aoaK04seq1 <- prune_taxa(taxa_sums(aoaK04seq)>0, aoaK04seq)

# Date: 06-01-2022
# 1. MINERAL
aoaM06seq<- subset_samples(aoa.physeq_bulk1, Date=="06-01-22" & Treatment=="M")
aoaM06seq1 <- prune_taxa(taxa_sums(aoaM06seq)>0, aoaM06seq)
# 2. BIODYNAMIC
aoaD06seq<- subset_samples(aoa.physeq_bulk1, Date=="06-01-22" & Treatment=="D")
aoaD06seq1 <- prune_taxa(taxa_sums(aoaD06seq)>0, aoaD06seq)
# 3. CONVENTIONAL
aoaK06seq<- subset_samples(aoa.physeq_bulk1, Date=="06-01-22" & Treatment=="K")
aoaK06seq1 <- prune_taxa(taxa_sums(aoaK06seq)>0, aoaK06seq)
################################################################################
# 1. Mineral 04-28-2022 rarefied
# asv phyloseq object
aoaM04seq1 # 206 taxa
################################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.02 % relative abundance across all samples
physeq.subset <- aoaM04seq1
physeq.subset #206 Taxa, 8 Samples
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0002)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 159 taxa, 8 samples remain in the data set after filtering

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
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 50,])
physeq.subset.50 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps) # 99 taxa

tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 60,])
physeq.subset.60 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.60 # 46 taxa

tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75 # 46 taxa
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
    glmT3s <- glmer(y ~ a + (1 | z), family='poisson',offset = o)
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum))
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
write.csv(glmT3s.pairwise.global, file = "AOA_M04_rare_pw_300823.csv")
aoam04_subset_table_75 <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aoam04_subset_table_75, file = "AOA_M04_rare_inpt_300823.csv")
################################################################################

# Subset the data set BY date

# Date: 04-28-2022
aoa04seq<- subset_samples(aoa.physeq_bulk1, Date=="04-28-22")
aoa04seq1 <- prune_taxa(taxa_sums(aoa04seq)>0, aoa04seq)
sort(rowSums(otu_table(aoa04seq1), na.rm = FALSE, dims = 1), decreasing = F)
aoa04seq1 #292 taxa 24 samples
################################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.02 % relative abundance across all samples
physeq.subset <- aoa04seq1
physeq.subset #292 Taxa, 24 Samples
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0002)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 167 taxa, 24 samples remain in the data set after filtering

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
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 50,])
physeq.subset.50 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps) # 127 taxa

tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 60,])
physeq.subset.60 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.60 # 90 taxa

tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75 # 90 taxa
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
    glmT3s <- glmer(y ~ a + (1 | z), family='poisson',offset = o)
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum))
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
write.csv(glmT3s.pairwise.global, file = "AOA_date1_rare_pw_300823.csv")
aoa04_subset_table_75 <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aoa04_subset_table_75, file = "AOA_date1_rare_inpt_300823.csv")

glmT3s.pairwise.global.rare.ed <- subset(glmT3s.pairwise.global.rare,glmT3s.pairwise.global.rare$contrast=="cont.D - rain.D" | glmT3s.pairwise.global.rare$contrast=="cont.M - rain.M" | glmT3s.pairwise.global.rare$contrast=="cont.K - rain.K")






















################################################################################
##### 1. BULK SOIL - NOT RAREFIED #####

aoa.physeq # still contain rhizosphere!, 646 taxa 192 samples
aoa.bulk.rawseq <- subset_samples(aoa.physeq, Type=="BS")
aoa.bulk.rawseq1 <- prune_taxa(taxa_sums(aoa.bulk.rawseq)>0,aoa.bulk.rawseq)
aoa.bulk.rawseq1 # 531 taxa

# separate data by treatment and date

# Date: 04-28-2022
# 1. MINERAL
aoaM04rawseq<- subset_samples(aoa.bulk.rawseq1, Date=="04-28-22" & Treatment=="M")
aoaM04rawseq1 <- prune_taxa(taxa_sums(aoaM04rawseq)>0, aoaM04rawseq)
aoaM04_table <- as.data.frame(otu_table(aoaM04rawseq1))
# 2. BIODYNAMIC
aoaD04rawseq<- subset_samples(aoa.bulk.rawseq1, Date=="04-28-22" & Treatment=="D")
aoaD04rawseq1 <- prune_taxa(taxa_sums(aoaD04rawseq)>0, aoaD04rawseq)
aoaD04_table <- as.data.frame(otu_table(aoaD04rawseq1))
# 3. CONVENTIONAL
aoaK04rawseq<- subset_samples(aoa.bulk.rawseq1, Date=="04-28-22" & Treatment=="K")
aoaK04rawseq1 <- prune_taxa(taxa_sums(aoaK04rawseq)>0, aoaK04rawseq)
aoaK04_table <- as.data.frame(otu_table(aoaK04rawseq1))

# Date: 06-01-2022
# 1. MINERAL
aoaM06rawseq<- subset_samples(aoa.bulk.rawseq1, Date=="06-01-22" & Treatment=="M")
aoaM06rawseq1 <- prune_taxa(taxa_sums(aoaM06rawseq)>0, aoaM06rawseq)
aoaM06_table <- as.data.frame(otu_table(aoaM06rawseq1))
# 2. BIODYNAMIC
aoaD06rawseq<- subset_samples(aoa.bulk.rawseq1, Date=="06-01-22" & Treatment=="D")
aoaD06rawseq1 <- prune_taxa(taxa_sums(aoaD06rawseq)>0, aoaD06rawseq)
aoaD06_table <- as.data.frame(otu_table(aoaD06rawseq1))
# 3. CONVENTIONAL
aoaK06rawseq<- subset_samples(aoa.bulk.rawseq1, Date=="06-01-22" & Treatment=="K")
aoaK06rawseq1 <- prune_taxa(taxa_sums(aoaK06rawseq)>0, aoaK06rawseq)
aoaK06_table <- as.data.frame(otu_table(aoaK06rawseq1))

################################################################################
################################################################################
# 1. Mineral 04-28-2022 not rarefied
# asv phyloseq object
aoaM04rawseq1 # 221 taxa
################################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.02 % relative abundance across all samples
physeq.subset <- aoaM04rawseq1
physeq.subset #221 Taxa, 8 Samples
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0002)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 159 taxa, 8 samples remain in the data set after filtering

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
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 50,])
physeq.subset.50 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps) # 101 taxa

tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 60,])
physeq.subset.60 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.60 # 48 taxa

tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
physeq.subset.75 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.75 # 48 taxa
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
    glmT3s <- glmer(y ~ a + (1 | z), family='poisson',offset = o)
    glmT3s.sum = summary(glmT3s)$coefficients
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum))
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
write.csv(glmT3s.pairwise.global, file = "AOA_M04_unrare_pw_300823.csv")
aoam04_subset_table_75 <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aoam04_subset_table_75, file = "AOA_M04_unrare_inpt_300823.csv")
