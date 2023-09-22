### DAA on AOB when the data is glommed into species instead of ASV 
# rarefied phyloseq object
aob.rare.1282.seq
# bulk soil only
aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS")
aob.physeq_bulk
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1 # 937 taxa, 119 samples
# glom into species
set.seed(13)
aob.sp <- tax_glom(aob.rare.1282.seq, taxrank = "Species", NArm = F)
aob.sp.tab <- as.data.frame(otu_table(aob.sp))
aob.sp.tax <- as.data.frame(tax_table(aob.sp))
rownames(aob.sp.tab) <- aob.sp.tax$Species
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
#write.csv(aob.sp.tab, file = "AOB.sp.tab_220923.csv")

# subset the data by date
# 1. 04/28/22
aob.spT1 <- subset_samples(aob.sp, Date=="04-28-22")
aob.spT1
aob.spT1.1 <- prune_taxa(taxa_sums(aob.spT1)>0, aob.spT1)
sort(rowSums(otu_table(aob.spT1.1), na.rm = FALSE, dims = 1), decreasing = F)
aob.spT1.1 # 18 taxa

# 1. 06/01/22
aob.spT2 <- subset_samples(aob.sp, Date=="06-01-22")
aob.spT2
aob.spT2.1 <- prune_taxa(taxa_sums(aob.spT2)>0, aob.spT2)
sort(rowSums(otu_table(aob.spT2.1), na.rm = FALSE, dims = 1), decreasing = F)
aob.spT2.1 # 19 taxa



###############################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.01 % relative abundance across all samples
physeq.subset <- aob.spT2.1
physeq.subset 
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0001)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 15 TAXA

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
df_otu_prev_ttt <- data.frame(matrix(ncol=nlevels(as.factor(df_tmp$Treatment)),
                                     nrow=nlevels(as.factor(df_tmp$OTU)), 
                                     dimnames=list(levels(as.factor(df_tmp$OTU)),
                                                   levels(as.factor(df_tmp$Treatment)))))
#attention il ya Sample et sample
for (i in unique(df_tmp$OTU)) {
  for (j in unique(df_tmp$Treatment)) {
    df_otu_prev_ttt[i,j] <- sum(df_tmp$sample[df_tmp$OTU == i & df_tmp$Treatment == j],na.rm = T) / nrow(df_tmp[df_tmp$OTU == i & df_tmp$Treatment == j,]) *100
    print(paste(i,j,df_otu_prev_ttt[i,j]),sep="\t")
    #print(df_otu_prev_ttt[i,j])
  }
  
}

df_otu_prev_ttt$max_prev <- apply(df_otu_prev_ttt,MARGIN=1, FUN=max)

# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 60,])
physeq.subset.60 <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
rm(ps,df_prev,tmp_otu_F)
physeq.subset.60  # 6 taxa

####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
#install.packages("glmmTMB", type="source")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset.60
str(tmp_T3s)
#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           "treatment"= as.character(tmp_T3s@sam_data$Treatment))
# force control as intercept
a[a == "D"] <- "1a"
a = as.factor(a$treatment)
# offset
o = log(sample_sums(tmp_T3s)) # using unfiltered data
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
    glmT3s <- glmmTMB(y ~ a + (1|z) , family="poisson", offset = o)
    #glmT3s <- glm(y ~ a, family='poisson')
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

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
write.csv(glmT3s.pairwise.global, file = "AOB_sp_treatment_pwc_220923.csv")
aob.spT1.tab <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aob.spT1.tab, file = "AOB_sp_filt_tab_220923.csv")
##############################################################################################################################



# By irrigation
# subset the data by date
# 1. 04/28/22
aob.spM1 <- subset_samples(aob.sp, Date == "04-28-22" & Treatment == "M")
aob.spM1
aob.spM1.1 <- prune_taxa(taxa_sums(aob.spM1)>0, aob.spM1)
sort(rowSums(otu_table(aob.spM1.1), na.rm = FALSE, dims = 1), decreasing = F)
aob.spM1.1 

aob.spD1 <- subset_samples(aob.sp, Date == "04-28-22" & Treatment == "D")
aob.spD1
aob.spD1.1 <- prune_taxa(taxa_sums(aob.spD1)>0, aob.spD1)
sort(rowSums(otu_table(aob.spD1.1), na.rm = FALSE, dims = 1), decreasing = F)
aob.spD1.1 

aob.spK1 <- subset_samples(aob.sp, Date == "04-28-22" & Treatment == "K")
aob.spK1
aob.spK1.1 <- prune_taxa(taxa_sums(aob.spK1)>0, aob.spK1)
sort(rowSums(otu_table(aob.spK1.1), na.rm = FALSE, dims = 1), decreasing = F)
aob.spK1.1 

aob.spM2 <- subset_samples(aob.sp, Date == "06-01-22" & Treatment == "M")
aob.spM2
aob.spM2.1 <- prune_taxa(taxa_sums(aob.spM2)>0, aob.spM2)
sort(rowSums(otu_table(aob.spM2.1), na.rm = FALSE, dims = 1), decreasing = F)
aob.spM2.1 

###############################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.01 % relative abundance across all samples
physeq.subset <- aob.spM2.1
physeq.subset 
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0001)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 10 TAXA

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
physeq.subset.75  # 6 taxa

####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
#install.packages("glmmTMB", type="source")
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
o = log(sample_sums(tmp_T3s)) # using unfiltered data
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
    glmT3s <- glmmTMB(y ~ a + (1|z) , family="poisson", offset = o)
    #glmT3s <- glm(y ~ a, family='poisson')
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

setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
write.csv(glmT3s.pairwise.global, file = "AOB_sp_treatment_pwc_220923.csv")
aob.spT1.tab <- as.data.frame(otu_table(physeq.subset.75))
write.csv(aob.spT1.tab, file = "AOB_sp_filt_tab_220923.csv")



















