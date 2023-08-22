
#######################################################
# AOB
######################################################

### Subset the data set per treatment-date

# Rarefied Data

# Date: 04-28-2022

dat04seq <- subset_samples(aob.physeq_bulk1, Date=="04-28-22")
dat04seq1 <- prune_taxa(taxa_sums(dat04seq)>0, dat04seq)
sort(rowSums(otu_table(dat04seq1), na.rm = FALSE, dims = 1), decreasing = F)
dat04seq1 #393 taxa, 23 samples

################################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.02 % relative abundance across all samples
physeq.subset <- dat04seq1
physeq.subset #393 Taxa, 23 Samples
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0002)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 254 taxa, 23 samples remain in the data set after filtering

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

#setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out')
#write.csv(df_otu_prev_ttt, file = "df_otu_prev_ttt_bulk.csv")

# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
physeq.subset <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
physeq.subset# 63 taxa, 23 samples


####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset

str(tmp_T3s)

#  treatment
a = tibble("sample"= as.factor(tmp_T3s@sam_data$SampleID),
           "treatment"= as.character(tmp_T3s@sam_data$x))
# force control as intercept
#a[a == "Control"] <- "1a"
a = as.matrix(a$treatment)
# offset
o = log(sample_sums(tmp_T3s))
# random effect
z <- as.matrix(tmp_T3s@sam_data$SampleID)
#tmp_T3s@sam_data$block <- paste(c("b"),tmp_T3s@sam_data$block, sep="")
# x <- as.factor(tmp_T3s@sam_data$block)

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

glmT3s.model.global.rare = glmT3s.sum.global
glmT3s.pairwise.global.rare = glmT3s.pairwise.global
glmT3s.pairwise.global.rare$p.adjust <- p.adjust(glmT3s.pairwise.global.rare$p.value, method = "fdr")

glmT3s.pairwise.global.rare.ed <- subset(glmT3s.pairwise.global.rare,glmT3s.pairwise.global.rare$contrast=="cont.D - rain.D" | glmT3s.pairwise.global.rare$contrast=="cont.M - rain.M" | glmT3s.pairwise.global.rare$contrast=="cont.K - rain.K")
glmT3s.pairwise.global.rare.ed$p.adjust <- p.adjust(glmT3s.pairwise.global.rare.ed$p.value, method = "fdr")

#nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value < glmT3s.pairwise.global$p.adjust,])
#nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value > glmT3s.pairwise.global$p.adjust,])
#nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value == glmT3s.pairwise.global$p.adjust,])


## nb of pval <= 0.05 before and after filter
table(glmT3s.pairwise.global.rare.ed$p.value <= 0.05)
table(glmT3s.pairwise.global.rare.ed$p.adjust <= 0.05)

## nb of OTU with a pval <= 0.05 before and after filter
tmp_otu3s.rare.ed = unique(glmT3s.pairwise.global.rare.ed$OTU[glmT3s.pairwise.global.rare.ed$p.adjust <= 0.05])
glmT3s.pairwise.global.signif.rare.ed = glmT3s.pairwise.global.rare.ed[glmT3s.pairwise.global.rare.ed$p.adjust <=0.05,]

length(tmp_otu3s) 

######### OJO
# save(glmT3s.pairwise.global, glmT3s.model.global, file = "glm-T3s.RData")
# load("/home/mathilde/Bureau/STAT/R-objects/glm-T3s.RData")

################################################################################

# NOT Rarefied Data

# Date: 04-28-2022

dat04rawseq <- subset_samples(aob.bulk.rawseq1, Date=="04-28-22")
dat04rawseq1 <- prune_taxa(taxa_sums(dat04rawseq)>0, dat04rawseq)
sort(rowSums(otu_table(dat04rawseq1), na.rm = FALSE, dims = 1), decreasing = F)
dat04rawseq1 #424 taxa, 24 samples

################################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.02 % relative abundance across all samples
physeq.subset <- dat04rawseq1
physeq.subset #424 Taxa, 24 Samples
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0002)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 264 taxa, 24 samples remain in the data set after filtering

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
physeq.subset <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
physeq.subset# 69 taxa, 24 samples


####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset

str(tmp_T3s)

#  treatment
a = tibble("sample"= as.factor(tmp_T3s@sam_data$SampleID),
           "treatment"= as.character(tmp_T3s@sam_data$x))
# force control as intercept
#a[a == "Control"] <- "1a"
a = as.matrix(a$treatment)
# offset
o = log(sample_sums(tmp_T3s))
# random effect
z <- as.matrix(tmp_T3s@sam_data$SampleID)
#tmp_T3s@sam_data$block <- paste(c("b"),tmp_T3s@sam_data$block, sep="")
# x <- as.factor(tmp_T3s@sam_data$block)

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




#glmT3s.pairwise.global.ed <- subset(glmT3s.pairwise.global,glmT3s.pairwise.global$contrast=="cont.D - rain.D" | glmT3s.pairwise.global$contrast=="cont.M - rain.M" | glmT3s.pairwise.global$contrast=="cont.K - rain.K")
#glmT3s.pairwise.global.ed$p.adjust <- p.adjust(glmT3s.pairwise.global.ed$p.value, method = "fdr")

nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value < glmT3s.pairwise.global$p.adjust,])
nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value > glmT3s.pairwise.global$p.adjust,])
nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value == glmT3s.pairwise.global$p.adjust,])



## nb of pval <= 0.05 before and after filter
table(glmT3s.pairwise.global$p.value <= 0.05)
table(glmT3s.pairwise.global$p.adjust <= 0.05)

## nb of OTU with a pval <= 0.05 before and after filter
tmp_otu3s = unique(glmT3s.pairwise.global$OTU[glmT3s.pairwise.global$p.adjust <= 0.05])
glmT3s.pairwise.global.signif = glmT3s.pairwise.global[glmT3s.pairwise.global$p.adjust <=0.05,]

length(tmp_otu3s) 

######### OJO
# save(glmT3s.pairwise.global, glmT3s.model.global, file = "glm-T3s.RData")
# load("/home/mathilde/Bureau/STAT/R-objects/glm-T3s.RData")


#######################################################
# COMA
######################################################

### Subset the data set per treatment-date

# Rarefied Data

# Date: 04-28-2022

dat04seq <- subset_samples(aob.physeq_bulk1, Date=="04-28-22")
dat04seq1 <- prune_taxa(taxa_sums(dat04seq)>0, dat04seq)
sort(rowSums(otu_table(dat04seq1), na.rm = FALSE, dims = 1), decreasing = F)
dat04seq1 #393 taxa, 23 samples

################################################################################
# Filter low-abundant taxa
# keeping OTUs with at least 0.02 % relative abundance across all samples
physeq.subset <- dat04seq1
physeq.subset #393 Taxa, 23 Samples
data.obs <- as.data.frame(otu_table(physeq.subset))
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0002)
data.F=data.obs[keep.taxa.id,,drop=FALSE]
new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object
physeq.subset # 254 taxa, 23 samples remain in the data set after filtering

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

#setwd('D:/Fina/INRAE_Project/microservices/070623_AOB_out')
#write.csv(df_otu_prev_ttt, file = "df_otu_prev_ttt_bulk.csv")

# filter otu par prevalence
physeq.subset 
ps =  physeq.subset 
df_prev = df_otu_prev_ttt
tmp_otu_F = rownames(df_prev[df_prev$max_prev >= 75,])
physeq.subset <- prune_taxa(taxa_names(ps) %in% tmp_otu_F, ps)
physeq.subset# 63 taxa, 23 samples


####################################################
# DIFFERENTIAL ABUNDANCE
##################################################
install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

tmp_T3s <- physeq.subset

str(tmp_T3s)

#  treatment
a = tibble("sample"= as.factor(tmp_T3s@sam_data$SampleID),
           "treatment"= as.character(tmp_T3s@sam_data$x))
# force control as intercept
#a[a == "Control"] <- "1a"
a = as.matrix(a$treatment)
# offset
o = log(sample_sums(tmp_T3s))
# random effect
z <- as.matrix(tmp_T3s@sam_data$SampleID)
#tmp_T3s@sam_data$block <- paste(c("b"),tmp_T3s@sam_data$block, sep="")
# x <- as.factor(tmp_T3s@sam_data$block)

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

glmT3s.pairwise.global.ed <- subset(glmT3s.pairwise.global,glmT3s.pairwise.global$contrast=="cont.D - rain.D" | glmT3s.pairwise.global$contrast=="cont.M - rain.M" | glmT3s.pairwise.global$contrast=="cont.K - rain.K")

glmT3s.pairwise.global.ed$p.adjust <- p.adjust(glmT3s.pairwise.global.ed$p.value, method = "fdr")

nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value < glmT3s.pairwise.global$p.adjust,])
nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value > glmT3s.pairwise.global$p.adjust,])
nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value == glmT3s.pairwise.global$p.adjust,])



## nb of pval <= 0.05 before and after filter
table(glmT3s.pairwise.global$p.value <= 0.05)
table(glmT3s.pairwise.global$p.adjust <= 0.05)

## nb of OTU with a pval <= 0.05 before and after filter
tmp_otu3s = unique(glmT3s.pairwise.global$OTU[glmT3s.pairwise.global$p.adjust <= 0.05])
glmT3s.pairwise.global.signif = glmT3s.pairwise.global[glmT3s.pairwise.global$p.adjust <=0.05,]

length(tmp_otu3s) 



