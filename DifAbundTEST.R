##################################################################################################
#DIFFERENTIAL ABUNDANCE ANALYSIS TEST
#################################################################################################

# GROUP & SEPARATE PHYLOSEQ OBJECT BY TYPE, DATE AND TREATMENT:

# 1. BULK SOIL
aob.physeq_bulk <- subset_samples(aob.rare.1282.seq, Type=="BS")
aob.physeq_bulk1 <- prune_taxa(taxa_sums(aob.physeq_bulk)>0, aob.physeq_bulk)
aob.physeq_bulk1
# Date: 04-28-2022
# 1. MINERAL
M04seq<- subset_samples(aob.physeq_bulk1, Date=="04-28-22" & Treatment=="M")
M04seq1 <- prune_taxa(taxa_sums(M04seq)>0, M04seq)
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
# Date: 07-05-2022
# 1. MINERAL
M09seq<- subset_samples(aob.physeq_bulk1, Date=="09-13-22" & Treatment=="M")
M09seq1 <- prune_taxa(taxa_sums(M09seq)>0, M09seq)
# 2. BIODYNAMIC
D09seq<- subset_samples(aob.physeq_bulk1, Date=="09-13-22" & Treatment=="D")
D09seq1 <- prune_taxa(taxa_sums(D09seq)>0, D09seq)
# 3. CONVENTIONAL
K09seq<- subset_samples(aob.physeq_bulk1, Date=="09-13-22" & Treatment=="K")
K09seq1 <- prune_taxa(taxa_sums(K09seq)>0, K09seq)
################################################################################
# Filter low-abundant taxa
################################################################################
# 1. BULK SOIL : 04-28-22 : MINERAL
physeq.subset1 <- M04seq1
data.obs1 <- as.data.frame(otu_table(physeq.subset1))

### keeping OTUs with at least 0.01 % relative abundance across all samples
keep.taxa.id1=which((rowSums(data.obs1)/sum(data.obs1))>0.0001)
data.F1=data.obs1[keep.taxa.id1,,drop=FALSE]

new.otu1 <- as.matrix(data.F1) # convert it into a matrix.
new.otu1 <- otu_table(data.F1, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset1) <- new.otu # incorporate into phyloseq Object


physeq.subset1 # 197 taxa remain in the data set after filtering
















































































