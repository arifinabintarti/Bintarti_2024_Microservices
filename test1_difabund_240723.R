aob.physeq_bulk1

# filtering low abundance of ASV
physeq.subset <- aob.physeq_bulk1
data.obs <- as.data.frame(otu_table(physeq.subset))

### keeping OTUs with at least 0.01 % relative abundance across all samples
keep.taxa.id=which((rowSums(data.obs)/sum(data.obs))>0.0001)
data.F=data.obs[keep.taxa.id,,drop=FALSE]

new.otu <- as.matrix(data.F) # convert it into a matrix.
new.otu <- otu_table(data.F, taxa_are_rows = TRUE) # convert into phyloseq compatible file.
otu_table(physeq.subset) <- new.otu # incorporate into phyloseq Object


physeq.subset # 428 ASV, 191 samples
####################################################
# DIFFERENTIAL ABUNDANCE
##################################################

tmp_T3s <- physeq.subset

str(tmp_T3s)

#  treatment
a = tibble("sample"= tmp_T3s@sam_data$SampleID,
           #"date" = tmp_T3s@sam_data$Date,
           #"fertilization" = tmp_T3s@sam_data$Treatment,
           "treatment"= tmp_T3s@sam_data$Irrigation)
# force control as intercept
a[a == "Control"] <- "1a"
a = as.factor(a$treatment)
#a.ed = a$treatment
#a. <- as.factor(a.ed)
#b
#b <- as.factor(tmp_T3s@sam_data$Treatment)
#c
#c <- as.factor(tmp_T3s@sam_data$Date)
# offset
o = log(sample_sums(tmp_T3s))
# random effect
z <- as.factor(tmp_T3s@sam_data$SampleID)
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
    glmT3s <- glmer(y ~ a + (1 | z), family='poisson', offset = o)
    glmT3s.sum = summary(glmT3s)$coefficients
    #glmT3s.sum = summary(glmT3s)$coefficients$cond
    glmT3s.sum = tibble("OTU"= OTU,
                        "treatment"=rownames(glmT3s.sum),
                        as_tibble(glmT3s.sum))
    glmT3s.sum
    
    glmT3s.sum.global = rbind(glmT3s.sum.global,glmT3s.sum)
    
    ### multiple comparaison
    
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
    tmp_df = tmp_df[tmp_df$a == "1a" | tmp_df$b == "1a" ,]
    
    # adjust pval n=nb of contrasts
    #tmp_df[,"p.adjust"] <- p.adjust(tmp_df$p.value,"fdr",n=21)
    #tmp_df[,"p.adjust"] <- p.adjust(tmp_df$p.value,"bonferroni",n=21)
    
    tmp_df = cbind("OTU"=OTU,tmp_df)
    
    glmT3s.pairwise.global = rbind(glmT3s.pairwise.global,tmp_df)
    
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  rm(OTU,y,glmT3s,glmT3s.sum)
  
}

glmT3s.model.global = glmT3s.sum.global
glmT3s.pairwise.global = glmT3s.pairwise.global

 glmT3s.pairwise.global$p.adjust <- p.adjust(glmT3s.pairwise.global$p.value, method = "fdr")

nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value < glmT3s.pairwise.global$p.adjust,])
nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value > glmT3s.pairwise.global$p.adjust,])
nrow(glmT3s.pairwise.global[glmT3s.pairwise.global$p.value == glmT3s.pairwise.global$p.adjust,])
