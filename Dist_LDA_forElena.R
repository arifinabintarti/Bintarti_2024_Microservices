library(usedist)

# Run CAP on your Bray-Curtis distance
cap <- CAPdiscrim(bray.dist ~ TxI, data = metadata, m = 39, permutations = 9999, add = TRUE) 
# define the group
groups <- sample_data(step1.rarefied) # Load the sample data from a phyloseq object, but you can load it directly from your metadata
groups <- groups$TxI # Treatment x Irrigation 
# Calculate distances (simple euclidean distance) from the CAP 
dist_matrix <- dist(cap$x) # x = the positions of the sites provided by the discriminant analysis (LDA)
# calculate distance between groups
d.calcul <- dist_groups(dist_matrix, groups)
# select the group comparisons within the same farming system
tab.distance = as_tibble(d.calcul) 
tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.D and rain.D","Within cont.D","Within rain.D","Between cont.K and rain.K","Within cont.K","Within rain.K","Between cont.M and rain.M","Within cont.M","Within rain.M"))
#tab.distance.C <- subset(tab.distance, Label%in% c("Between cont.M and rain.M","Between cont.D and rain.D","Between cont.K and rain.K"))
tab.distance.C$Label <- factor(tab.distance.C$Label)
tab.distance.C$Group2 <- as.character(tab.distance.C$Group2) 
tab.distance.C.ed <- tab.distance.C %>%  #adding one column based on the fertilization
  mutate(Treatment = case_when(
    endsWith(Group2, "D") ~ "BIODYN (D)",
    endsWith(Group2, "K") ~ "CONFYM (K)",
    endsWith(Group2, "M") ~ "CONMIN (M)"))
tab.distance.C.ed2 <- tab.distance.C.ed %>% # re-order level of the Label for the plot
  mutate(Label = factor(Label, 
                        levels = c("Between cont.D and rain.D",
                                   "Within cont.D","Within rain.D",
                                   "Between cont.K and rain.K",
                                   "Within cont.K","Within rain.K",
                                   "Between cont.M and rain.M",
                                   "Within cont.M","Within rain.M")))

# Kruskal-Wallis Test
kruskal.test(Distance ~ Label,
             data = tab.distance.C.ed2) 
# Post Hoc Dunn Test
dunn <- dunnTest(Distance ~ Label,
                        data = tab.distance.C.ed2, method = "bh")
dunn

# The mean within-treatment (within-control & within-drought) distances were significantly lower than 
# between-treatment distances (between control-drought), indicating lower diversity and closer similarity within tretment, 
# and that there is a drought effect that lead to the dissimilarity between the two group of treatment.


# Make the significance letter
CLD.all <- cldList(P.adj ~ Comparison, data=dunn$res)
CLD.all
CLD.all.ord <- CLD.all[c(1,4,7,2,5,8,3,6,9),] # re-order
# Plot
sumData.all <- ddply(tab.distance.C.ed2, "Label", summarise,
                     Max = max(Distance),
                     N    = length(Distance),
                     Mean = mean(Distance),
                     Sd   = sd(Distance),
                     Se   = Sd / sqrt(N))
rownames(CLD.all.ord) = sumData.all$Label
CLD.all.ed <- rownames_to_column(CLD.all.ord, var="Label") # to have the same label
CLD.all.ed

CLD.all.ed2 <- CLD.all.ed %>%  #adding one column based on the farming system
  mutate(Treatment = case_when(
    endsWith(Group, "D") ~ "BIODYN (D)",
    endsWith(Group, "K") ~ "CONFYM (K)",
    endsWith(Group, "M") ~ "CONMIN (M)"))

tab.distance.C.ed2$Treatment <- as.factor(tab.distance.C.ed2$Treatment)
# Plot
plot <- ggplot(tab.distance.C.ed2, aes(x = Label, y = Distance)) + 
  geom_boxplot(width=0.7,lwd=0.7,aes(fill=Label)) + 
  scale_fill_manual(values = c("#009E73","#6DC7AF","#DAF1EB","#FF618C","#FF8EAC","#FFE8EE","#E69F00","#EDBA48","#FBF1DA"))+
  scale_x_discrete(labels = c("Between","within \n control","within \n drought","Between","within \n control","within \n drought","Between","within \n control","within \n drought"))+
  geom_text(data=CLD.all.ed2,aes(x=Label, y = sumData.all$Max + 1, label=Letter), vjust=0, size=6) +
  ylab("Distances") +
  facet_wrap(~ Treatment,scales = "free_x")+
  ylim(0,15)+
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_text(size=13, face='bold'),
        plot.title = element_text(hjust = 0.5, size = 15, face='bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=15), 
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        legend.text=element_text(size=13),
        legend.spacing.x = unit(0.05, 'cm'))
plot

