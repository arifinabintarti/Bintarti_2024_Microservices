# ################################################################################
# Microservices Project: Percentage of Drought-responsive ASVs
################################################################################
# Date : 26 March 2024
# Author : Ari Fina Bintarti


# 1. AOB-BS
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
rr.AOB <- read.csv("AOB_RR_Bulk_p.val_270224.csv", row.names = 1)
names(rr.AOB)=str_sub(names(rr.AOB),4)
View(rr.AOB)

results_list <- list()
# Calculate the percentage of increased and decreased OTUs for each treatment
treatments <- colnames(rr.AOB)  # get the number of treatments
# please be careful with the divisor (length(taxa_names(physeq.subset))), as it is the dominant OTUs
for (treatment in treatments) {
  increased <- sum(rr.AOB[[treatment]] > 0) / length(unique(glmT3s.pairwise.global.AOB$OTU))
  decreased <- sum(rr.AOB[[treatment]] < 0) / length(unique(glmT3s.pairwise.global.AOB$OTU))
  results_list <- c(
    results_list,
    list(data.frame(response = "positive", treatment = treatment, percentage = increased)),
    list(data.frame(response = "negative", treatment = treatment, percentage = -decreased))
  )
}
# Combine the results into a single dataframe and sort the treatments according to the desired order
result.aob.BS <- do.call(rbind, results_list)
view(result.aob.BS)
str(result.aob.BS)
# adding column 'fertilization'
result.aob.BS.ed <- result.aob.BS %>%
 mutate(fertilization = map_chr(treatment, ~unlist(strsplit(., ""))[4])) %>%
 mutate(fertilization = case_when(fertilization == "D" ~ "BIODYN",
 fertilization == "K" ~ "CONFYM",
 fertilization == "M" ~ "CONMIN"))
view(result.aob.BS.ed)
# adding column 'date'
result.aob.BS.ed2 <- result.aob.BS.ed %>%
 mutate(date = map_chr(treatment, ~unlist(strsplit(., ""))[9])) %>%
 mutate(date = case_when(date == "8" ~ "Apr",
 date == "1" ~ "Jun",
 date == "5" ~ "Jul5",
 date == "0" ~ "Jul20",
 date == "3" ~ "Sep"))
view(result.aob.BS.ed2)
str(result.aob.BS.ed2)
result.aob.BS.ed2$treatment <- factor(result.aob.BS.ed2$treatment)
result.aob.BS.ed2$fertilization <- factor(result.aob.BS.ed2$fertilization)
# change factor level for date
result.aob.BS.ed2$date <- factor(result.aob.BS.ed2$date,levels = c("Apr", "Jun", "Jul5", "Jul20", "Sep"))
## Make the graph regarding the % of OTUs that are postiively or negatively changing in the treatments
abund_counts_response.AOB <- ggplot(result.aob.BS.ed2, aes(x = date, y = percentage*100, fill = response)) +
  geom_bar(stat = "identity")+geom_hline(yintercept =0,color="white")+ scale_fill_manual(values=c("#FF6666","#6666FF"))+geom_hline(yintercept=0, color="black")+
  facet_wrap(~ fertilization)+
  #labs(title="Bulk Soil")+#subtitle="A. AOB")+
  labs(title="A. AOB")+
  theme_bw() + theme(legend.position="none",
                     plot.title = element_text(size=25, face="bold"),
                     #plot.subtitle = element_text(size=20, face="bold"),
                     #axis.text.x = element_text(angle=90, size = 16, hjust=1, vjust=0.5),
                     strip.text.x= element_text(size = 20),
                     strip.background = element_rect(fill="white"),
                     axis.text.y = element_text(size = 16),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.y = element_text(size=18),
                     axis.title.x = element_blank(),
                     legend.text = element_text(size=16),
                     legend.title = element_text(size=18),
                     panel.grid = element_blank())+ylab("Drought-affected ASVs (%)")+
  scale_y_continuous(limits = c(-10,10),breaks = c(-10, -5, 0, 5, 10))
abund_counts_response.AOB

# 2. AOB-RS
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
rr.AOB.RS <- read.csv("AOB_RR_Rhizo_p.val_270224.csv", row.names = 1)
View(rr.AOB.RS)

results_list <- list()
# Calculate the percentage of increased and decreased OTUs for each treatment
treatments <- colnames(rr.AOB.RS)  # get the number of treatments
# please be careful with the divisor (length(taxa_names(physeq.subset))), as it is the dominant OTUs
for (treatment in treatments) {
  increased <- sum(rr.AOB.RS[[treatment]] > 0) / length(unique(glmT3s.pairwise.global.AOB.RS$OTU))
  decreased <- sum(rr.AOB.RS[[treatment]] < 0) / length(unique(glmT3s.pairwise.global.AOB.RS$OTU))
  results_list <- c(
    results_list,
    list(data.frame(response = "positive", treatment = treatment, percentage = increased)),
    list(data.frame(response = "negative", treatment = treatment, percentage = -decreased))
  )
}
# Combine the results into a single dataframe and sort the treatments according to the desired order
result.aob.RS <- do.call(rbind, results_list)
view(result.aob.RS)
str(result.aob.RS)
# adding column 'fertilization'
result.aob.RS.ed <- result.aob.RS %>%
 mutate(fertilization = map_chr(treatment, ~unlist(strsplit(., ""))[4])) %>%
 mutate(fertilization = case_when(fertilization == "D" ~ "BIODYN",
 fertilization == "K" ~ "CONFYM",
 fertilization == "M" ~ "CONMIN"))
view(result.aob.RS.ed)
# adding column 'date'
result.aob.RS.ed2 <- result.aob.RS.ed %>%
 mutate(date = map_chr(treatment, ~unlist(strsplit(., ""))[7])) %>%
 mutate(date = case_when(date == "4" ~ "Apr",
 date == "6" ~ "Jun",
 date == "7" ~ "Jul"))
 #date == "0" ~ "Jul20"))
 #date == "3" ~ "Sep"))
view(result.aob.RS.ed2)
str(result.aob.RS.ed2)
result.aob.RS.ed2$treatment <- factor(result.aob.RS.ed2$treatment)
result.aob.RS.ed2$fertilization <- factor(result.aob.RS.ed2$fertilization)
# change factor level for date
result.aob.RS.ed2$date <- factor(result.aob.RS.ed2$date,levels = c("Apr", "Jun", "Jul"))
## Make the graph regarding the % of OTUs that are postiively or negatively changing in the treatments
abund_counts_response.AOB.RS <- ggplot(result.aob.RS.ed2, aes(x = date, y = percentage*100, fill = response)) +
  geom_bar(stat = "identity")+geom_hline(yintercept =0,color="white")+ scale_fill_manual(values=c("#FF6666","#6666FF"))+geom_hline(yintercept=0, color="black")+
  facet_wrap(~ fertilization)+
  labs(title="Rhizosphere")+#subtitle="A. AOB")+
  theme_bw() + theme(legend.position="none",
                     plot.title = element_text(size=25, face="bold"),
                     #plot.subtitle = element_text(size=20, face="bold"),
                     #axis.text.x = element_text(angle=90, size = 16, hjust=1, vjust=0.5),
                     strip.text.x= element_text(size = 20),
                     strip.background = element_rect(fill="white"),
                     axis.text.y = element_text(size = 16),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     #axis.title.y = element_text(size=18),
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     legend.text = element_text(size=16),
                     legend.title = element_text(size=18),
                     panel.grid = element_blank())+ylab("AOB\nDrought-affected ASVs (%)")+
  scale_y_continuous(limits = c(-10,10),breaks = c(-10, -5, 0, 5, 10))
abund_counts_response.AOB.RS

#__________________________________________________________________________________________________________________________________________________
# 3. AOA-BS
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
rr.AOA.BS <- read.csv("AOA_RR_Bulk_p.val_280224.csv", row.names = 1)
View(rr.AOA.BS)

results_list <- list()
# Calculate the percentage of increased and decreased OTUs for each treatment
treatments <- colnames(rr.AOA.BS)  # get the number of treatments
# please be careful with the divisor (length(taxa_names(physeq.subset))), as it is the dominant OTUs
for (treatment in treatments) {
  increased <- sum(rr.AOA.BS[[treatment]] > 0) / length(unique(AOA.glmT3s.pairwise.global.ALL$OTU))
  decreased <- sum(rr.AOA.BS[[treatment]] < 0) / length(unique(AOA.glmT3s.pairwise.global.ALL$OTU))
  results_list <- c(
    results_list,
    list(data.frame(response = "positive", treatment = treatment, percentage = increased)),
    list(data.frame(response = "negative", treatment = treatment, percentage = -decreased))
  )
}
# Combine the results into a single dataframe and sort the treatments according to the desired order
result.aoa.BS <- do.call(rbind, results_list)
view(result.aoa.BS)
str(result.aoa.BS)
# adding column 'fertilization'
result.aoa.BS.ed <- result.aoa.BS %>%
 mutate(fertilization = map_chr(treatment, ~unlist(strsplit(., ""))[4])) %>%
 mutate(fertilization = case_when(fertilization == "D" ~ "BIODYN",
 fertilization == "K" ~ "CONFYM",
 fertilization == "M" ~ "CONMIN"))
view(result.aoa.BS.ed)
# adding column 'date'
result.aoa.BS.ed2 <- result.aoa.BS.ed %>%
 mutate(date = map_chr(treatment, ~unlist(strsplit(., ""))[9])) %>%
 mutate(date = case_when(date == "8" ~ "Apr",
 date == "1" ~ "Jun",
 date == "5" ~ "Jul5",
 date == "0" ~ "Jul20",
 date == "3" ~ "Sep"))
view(result.aoa.BS.ed2)
str(result.aoa.BS.ed2)
result.aoa.BS.ed2$treatment <- factor(result.aoa.BS.ed2$treatment)
result.aoa.BS.ed2$fertilization <- factor(result.aoa.BS.ed2$fertilization)
# change factor level for date
result.aoa.BS.ed2$date <- factor(result.aoa.BS.ed2$date,levels = c("Apr", "Jun", "Jul5", "Jul20", "Sep"))
## Make the graph regarding the % of OTUs that are postiively or negatively changing in the treatments
abund_counts_response.AOA <- ggplot(result.aoa.BS.ed2, aes(x = date, y = percentage*100, fill = response)) +
  geom_bar(stat = "identity")+geom_hline(yintercept =0,color="white")+ scale_fill_manual(values=c("#FF6666","#6666FF"))+geom_hline(yintercept=0, color="black")+
  facet_wrap(~ fertilization)+
  labs(title="B. AOA")+
  theme_bw() + theme(legend.position="none",
                     plot.title = element_text(size=25, face="bold"),
                     #plot.subtitle = element_text(size=20, face="bold"),
                     #axis.text.x = element_text(angle=90, size = 16, hjust=1, vjust=0.5),
                     strip.text= element_blank(),
                     #strip.background = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_text(size = 16),
                     axis.title.y = element_text(size=18),
                     axis.title.x = element_blank(),
                     legend.text = element_text(size=16),
                     legend.title = element_text(size=18),
                     panel.grid = element_blank())+ylab("Drought-affected ASVs (%)")+
  scale_y_continuous(limits = c(-10,10),breaks = c(-10, -5, 0, 5, 10))
abund_counts_response.AOA

# 4. AOA-RS
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
rr.AOA.RS <- read.csv("AOA_RR_Rhizo_p.val_280224.csv", row.names = 1)
View(rr.AOA.RS)

results_list <- list()
# Calculate the percentage of increased and decreased OTUs for each treatment
treatments <- colnames(rr.AOA.RS)  # get the number of treatments
# please be careful with the divisor (length(taxa_names(physeq.subset))), as it is the dominant OTUs
for (treatment in treatments) {
  increased <- sum(rr.AOA.RS[[treatment]] > 0) / length(unique(glmT3s.pairwise.global.AOA.RS$OTU))
  decreased <- sum(rr.AOA.RS[[treatment]] < 0) / length(unique(glmT3s.pairwise.global.AOA.RS$OTU))
  results_list <- c(
    results_list,
    list(data.frame(response = "positive", treatment = treatment, percentage = increased)),
    list(data.frame(response = "negative", treatment = treatment, percentage = -decreased))
  )
}
# Combine the results into a single dataframe and sort the treatments according to the desired order
result.aoa.RS <- do.call(rbind, results_list)
view(result.aoa.RS)
str(result.aoa.RS)
# adding column 'fertilization'
result.aoa.RS.ed <- result.aoa.RS %>%
 mutate(fertilization = map_chr(treatment, ~unlist(strsplit(., ""))[4])) %>%
 mutate(fertilization = case_when(fertilization == "D" ~ "BIODYN",
 fertilization == "K" ~ "CONFYM",
 fertilization == "M" ~ "CONMIN"))
view(result.aoa.RS.ed)
# adding column 'date'
result.aoa.RS.ed2 <- result.aoa.RS.ed %>%
 mutate(date = map_chr(treatment, ~unlist(strsplit(., ""))[7])) %>%
 mutate(date = case_when(date == "4" ~ "Apr",
 date == "6" ~ "Jun",
 date == "7" ~ "Jul"))
 #date == "0" ~ "Jul20"))
 #date == "3" ~ "Sep"))
view(result.aoa.RS.ed2)
str(result.aoa.RS.ed2)
result.aoa.RS.ed2$treatment <- factor(result.aoa.RS.ed2$treatment)
result.aoa.RS.ed2$fertilization <- factor(result.aoa.RS.ed2$fertilization)
# change factor level for date
result.aoa.RS.ed2$date <- factor(result.aoa.RS.ed2$date,levels = c("Apr", "Jun", "Jul"))
## Make the graph regarding the % of OTUs that are postiively or negatively changing in the treatments
abund_counts_response.AOA.RS <- ggplot(result.aoa.RS.ed2, aes(x = date, y = percentage*100, fill = response)) +
  geom_bar(stat = "identity")+geom_hline(yintercept =0,color="white")+ scale_fill_manual(values=c("#FF6666","#6666FF"))+geom_hline(yintercept=0, color="black")+
  facet_wrap(~ fertilization)+
  #labs(title="Rhizosphere")+#subtitle="A. AOB")+
  theme_bw() + theme(legend.position="none",
                     plot.title = element_text(size=25, face="bold"),
                     #plot.subtitle = element_text(size=20, face="bold"),
                     #axis.text.x = element_text(angle=90, size = 16, hjust=1, vjust=0.5),
                     strip.text= element_blank(),
                     #strip.text.x= element_text(size = 20),
                     strip.background = element_rect(fill="white"),
                     axis.text.y = element_text(size = 16),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     #axis.title.y = element_text(size=18),
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     legend.text = element_text(size=16),
                     legend.title = element_text(size=18),
                     panel.grid = element_blank())+ylab("AOA\nDrought-affected ASVs (%)")+
  scale_y_continuous(limits = c(-10,10),breaks = c(-10, -5, 0, 5, 10))
abund_counts_response.AOA.RS

#__________________________________________________________________________________________________________________________________________________
# 4. COMA-BS
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
rr.COM.BS <- read.csv("COM_RR_Bulk_p.val_280224.csv", row.names = 1)
View(rr.COM.BS)

results_list <- list()
# Calculate the percentage of increased and decreased OTUs for each treatment
treatments <- colnames(rr.COM.BS)  # get the number of treatments
# please be careful with the divisor (length(taxa_names(physeq.subset))), as it is the dominant OTUs
for (treatment in treatments) {
  increased <- sum(rr.COM.BS[[treatment]] > 0) / length(unique(glmT3s.pairwise.global.COM.BS$OTU))
  decreased <- sum(rr.COM.BS[[treatment]] < 0) / length(unique(glmT3s.pairwise.global.COM.BS$OTU))
  results_list <- c(
    results_list,
    list(data.frame(response = "positive", treatment = treatment, percentage = increased)),
    list(data.frame(response = "negative", treatment = treatment, percentage = -decreased))
  )
}
# Combine the results into a single dataframe and sort the treatments according to the desired order
result.com.BS <- do.call(rbind, results_list)
view(result.com.BS)
str(result.com.BS)
# adding column 'fertilization'
result.com.BS.ed <- result.com.BS %>%
 mutate(fertilization = map_chr(treatment, ~unlist(strsplit(., ""))[4])) %>%
 mutate(fertilization = case_when(fertilization == "D" ~ "BIODYN",
 fertilization == "K" ~ "CONFYM",
 fertilization == "M" ~ "CONMIN"))
view(result.com.BS.ed)
# adding column 'date'
result.com.BS.ed2 <- result.com.BS.ed %>%
 mutate(date = map_chr(treatment, ~unlist(strsplit(., ""))[9])) %>%
 mutate(date = case_when(date == "8" ~ "Apr",
 date == "1" ~ "Jun",
 date == "5" ~ "Jul5",
 date == "0" ~ "Jul20",
 date == "3" ~ "Sep"))
view(result.com.BS.ed2)
str(result.com.BS.ed2)
result.com.BS.ed2$treatment <- factor(result.com.BS.ed2$treatment)
result.com.BS.ed2$fertilization <- factor(result.com.BS.ed2$fertilization)
# change factor level for date
result.com.BS.ed2$date <- factor(result.com.BS.ed2$date,levels = c("Apr", "Jun", "Jul5", "Jul20", "Sep"))
## Make the graph regarding the % of OTUs that are postiively or negatively changing in the treatments
abund_counts_response.COM <- ggplot(result.com.BS.ed2, aes(x = date, y = percentage*100, fill = response)) +
  geom_bar(stat = "identity")+geom_hline(yintercept =0,color="white")+ scale_fill_manual(values=c("#FF6666","#6666FF"))+geom_hline(yintercept=0, color="black")+
  facet_wrap(~ fertilization)+
  labs(title="C. Comammox")+
  theme_bw() + theme(legend.position="none",
                     plot.title = element_text(size=25, face="bold"),
                     #plot.subtitle = element_text(size=20, face="bold"),
                     axis.text.x = element_text(angle=90, size = 16, hjust=1, vjust=0.5),
                     strip.text= element_blank(),
                     #strip.background = element_blank(),
                     axis.text.y = element_text(size = 16),
                     axis.title.y = element_text(size=18),
                     axis.title.x = element_blank(),
                     legend.text = element_text(size=16),
                     legend.title = element_text(size=18),
                     panel.grid = element_blank())+ylab("Drought-affected ASVs (%)")+
  scale_y_continuous(limits = c(-10,10),breaks = c(-10, -5, 0, 5, 10))
abund_counts_response.COM

# 6. Comammox-RS
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
rr.COM.RS <- read.csv("COM_RR_Rhizo_p.val_280224.csv", row.names = 1)
View(rr.COM.RS)

results_list <- list()
# Calculate the percentage of increased and decreased OTUs for each treatment
treatments <- colnames(rr.COM.RS)  # get the number of treatments
# please be careful with the divisor (length(taxa_names(physeq.subset))), as it is the dominant OTUs
for (treatment in treatments) {
  increased <- sum(rr.COM.RS[[treatment]] > 0) / length(unique(glmT3s.pairwise.global.COM.RS$OTU))
  decreased <- sum(rr.COM.RS[[treatment]] < 0) / length(unique(glmT3s.pairwise.global.COM.RS$OTU))
  results_list <- c(
    results_list,
    list(data.frame(response = "positive", treatment = treatment, percentage = increased)),
    list(data.frame(response = "negative", treatment = treatment, percentage = -decreased))
  )
}
# Combine the results into a single dataframe and sort the treatments according to the desired order
result.com.RS <- do.call(rbind, results_list)
view(result.com.RS)
str(result.com.RS)
# adding column 'fertilization'
result.com.RS.ed <- result.com.RS %>%
 mutate(fertilization = map_chr(treatment, ~unlist(strsplit(., ""))[4])) %>%
 mutate(fertilization = case_when(fertilization == "D" ~ "BIODYN",
 fertilization == "K" ~ "CONFYM",
 fertilization == "M" ~ "CONMIN"))
view(result.com.RS.ed)
# adding column 'date'
result.com.RS.ed2 <- result.com.RS.ed %>%
 mutate(date = map_chr(treatment, ~unlist(strsplit(., ""))[7])) %>%
 mutate(date = case_when(date == "4" ~ "Apr",
 date == "6" ~ "Jun",
 date == "7" ~ "Jul"))
 #date == "0" ~ "Jul20"))
 #date == "3" ~ "Sep"))
view(result.com.RS.ed2)
str(result.com.RS.ed2)
result.com.RS.ed2$treatment <- factor(result.com.RS.ed2$treatment)
result.com.RS.ed2$fertilization <- factor(result.com.RS.ed2$fertilization)
# change factor level for date
result.com.RS.ed2$date <- factor(result.com.RS.ed2$date,levels = c("Apr", "Jun", "Jul"))
## Make the graph regarding the % of OTUs that are postiively or negatively changing in the treatments
abund_counts_response.COM.RS <- ggplot(result.com.RS.ed2, aes(x = date, y = percentage*100, fill = response)) +
  geom_bar(stat = "identity")+geom_hline(yintercept =0,color="white")+ scale_fill_manual(values=c("#FF6666","#6666FF"))+geom_hline(yintercept=0, color="black")+
  facet_wrap(~ fertilization)+
  #labs(subtitle="C. Comammox")+
  theme_bw() + theme(legend.position="none",
                     plot.title = element_text(size=25, face="bold"),
                     #plot.subtitle = element_text(size=20, face="bold"),
                     axis.text.x = element_text(angle=90, size = 16, hjust=1, vjust=0.5),
                     strip.text= element_blank(),
                     #strip.text.x= element_text(size = 20),
                     strip.background = element_rect(fill="white"),
                     axis.text.y = element_text(size = 16),
                     #axis.title.y = element_text(size=18),
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     legend.text = element_text(size=16),
                     legend.title = element_text(size=18),
                     panel.grid = element_blank())+ylab("Comammox\nDrought-affected ASVs (%)")+
  scale_y_continuous(limits = c(-10,10),breaks = c(-10, -5, 0, 5, 10))
abund_counts_response.COM.RS

#####################################################################################################################################################
# All figures

percentage.ASV.plot <- ((abund_counts_response.AOB / abund_counts_response.AOA / abund_counts_response.COM) | (abund_counts_response.AOB.RS / abund_counts_response.AOA.RS / abund_counts_response.COM.RS)) +
 plot_layout(guides = "collect") & 
 theme(legend.position = 'bottom')
percentage.ASV.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Fig.8dpi600.tiff",
       percentage.ASV.plot, device = "tiff",bg="white",
       width = 10, height = 12, 
       units= "in", dpi = 600) #compression="lzw")

percentage.ASV.BS.plot <- (abund_counts_response.AOB / abund_counts_response.AOA / abund_counts_response.COM)+
 plot_layout(guides = "collect") & 
 theme(legend.position = 'bottom')
percentage.ASV.BS.plot
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("Fig.perc.ASV.BS.tiff",
       percentage.ASV.BS.plot, device = "tiff",bg="white",
       width = 6, height = 12, 
       units= "in", dpi = 600,compression="lzw")


