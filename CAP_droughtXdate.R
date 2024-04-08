### CAP ANALYSIS FOR DROUGHT X SAMPLING DATE ###
# Author: Ari Fina Bintarti
# Date: 07/03/2024

library(parallel)
library(BiodiversityR) # ALWAYS LOAD FROM THE R CONSOLE!!!!
library(ggforce)

# I. BULK SOIL 

# 1. AOA

# run Bray-Curtis beta diversity on bulk soil
aoa.asv.bulk <- aoarare.asv.df[,1:120]
aoa.asv.bulk1 <- aoa.asv.bulk[rowSums(aoa.asv.bulk)>0,]
sort(rowSums(aoa.asv.bulk1, na.rm = FALSE, dims = 1), decreasing = FALSE)
aoa.bulk_dist_bc <- vegdist(t(aoa.asv.bulk1), method = "bray")
aoa.bulk_dist_bc
# metadata
aoa.meta.bulk <- aoa.meta.df[1:120,]
aoa.meta.bulk$period2 <- factor(aoa.meta.bulk$period, levels = c("Control.D", "Control.R11W", "Control.R1W", 
                                                              "Rainout.D", "Rainout.R11W",  "Rainout.R1W"),
                          labels = c("Drought", "R11", "R1", "Drought", "R11", "R1"))
#aoa.meta.bulk.ed <- aoa.meta.bulk[,-30:-45]
str(aoa.meta.bulk)
aoa.meta.bulk$SampleID<-factor(aoa.meta.bulk$SampleID)
aoa.meta.bulk$PlotID<-factor(aoa.meta.bulk$PlotID)
aoa.meta.bulk$Irrigation<-factor(aoa.meta.bulk$Irrigation)
aoa.meta.bulk$Block<-factor(aoa.meta.bulk$Block)
aoa.meta.bulk$x<-factor(aoa.meta.bulk$x)
aoa.meta.bulk$rep<-factor(aoa.meta.bulk$rep)
aoa.meta.bulk$rep2<-factor(aoa.meta.bulk$rep2)
aoa.meta.bulk$DxD<-factor(aoa.meta.bulk$DxD)
aoa.meta.bulk$period<-factor(aoa.meta.bulk$period)
aoa.meta.bulk$period2<-factor(aoa.meta.bulk$period2)
aoa.meta.bulk$var2<-factor(aoa.meta.bulk$var2)
aoa.meta.bulk$var3<-factor(aoa.meta.bulk$var3)

# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(aoa.bulk_dist_bc ~ period, data = aoa.meta.bulk, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 32 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.BS.DxPer <- CAPdiscrim(aoa.bulk_dist_bc ~ period, data = aoa.meta.bulk, m = 32, permutations = 999, add = TRUE) # % 
success <- cbind(data.frame(aoa.cap.BS.DxPer$group), data.frame(aoa.cap.BS.DxPer$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.BS.DxPer$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1.BS.DxPer <- paste("CAP1 (", round((100/sum(aoa.cap.BS.DxPer$lda.other$svd^2) * aoa.cap.BS.DxPer$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aoa.cap2.BS.DxPer <- paste("CAP2 (", round((100/sum(aoa.cap.BS.DxPer$lda.other$svd^2) * aoa.cap.BS.DxPer$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with ggplot2

aoa.cap.BS.DxPer.plot <- ggplot(as.data.frame(aoa.cap.BS.DxPer$x), aes(x = aoa.cap.BS.DxPer$x[,1], y = aoa.cap.BS.DxPer$x[,2])) +
  geom_point(aes(color = aoa.meta.bulk$Irrigation, shape = aoa.meta.bulk$period2), size = 4) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aoa.cap1.BS.DxPer) + ylab(aoa.cap2.BS.DxPer) +
  scale_color_manual(values = c("#CC6677","#117733"),
                     name = "Irrigation treatment",
                     labels = c("Control", "Drought")) +
  scale_shape_manual(values = c(15, 16, 3),
                     name = "Sampling Date",
                     labels = c("Drought", "Rewetting 1 week","Rewetting 11 weeks")) + theme_classic() +
  scale_fill_manual(values = c("#CC6677","#CC6677","#CC6677", 
                               "#117733","#117733","#117733")) +
  geom_mark_ellipse(aes(fill = aoa.meta.bulk$period), 
                   expand = 0, linewidth = NA, show.legend = F)  +
  #labs(title = "AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.position = "right",
        legend.title = element_text(size=14),
        legend.text = element_text(size=14)) 
  #annotate("text",x=-4.4,y=-3.6,hjust = 0, size = 5,
  #label= "Control (drought): 47.2%\nControl (R1W): 0%\nControl (R11W): 50%\nDrought (drought): 55.6%\nDrought (R1W): 8.3%\nDrought (R11W): 25%")
aoa.cap.BS.DxPer.plot


#setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOA_CAP_BS_DxPer.tiff",
       aoa.cap.BS.DxPer.plot, device = "tiff",
       width = 12, height =9, 
       units= "in", dpi = 600)

#___________________________________________________________________________________________________________________

# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:110) {
  cap <- CAPdiscrim(aoa.bulk_dist_bc ~ DxD, data = aoa.meta.bulk, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 100 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.BS.DxD <- CAPdiscrim(aoa.bulk_dist_bc ~ DxD, data = aoa.meta.bulk, m = 100, permutations = 999, add = TRUE) # 19.2% 
success <- cbind(data.frame(aoa.cap.BS.DxD$group), data.frame(aoa.cap.BS.DxD$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.BS.DxD$PCoA)
success <- success[order(success$source), ]
success

aoa.cap1.BS.DxD <- paste("CAP1 (", round((100/sum(aoa.cap.BS.DxD$lda.other$svd^2) * aoa.cap.BS.DxD$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aoa.cap2.BS.DxD <- paste("CAP2 (", round((100/sum(aoa.cap.BS.DxD$lda.other$svd^2) * aoa.cap.BS.DxD$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with ggplot2

aoa.cap.BS.DxD.plot <- ggplot(as.data.frame(aoa.cap.BS.DxD$x), aes(x = aoa.cap.BS.DxD$x[,1], y = aoa.cap.BS.DxD$x[,2])) +
  geom_point(aes(color = aoa.meta.bulk$Irrigation, shape = aoa.meta.bulk$Date), size = 4) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aoa.cap1.BS.DxD) + ylab(aoa.cap2.BS.DxD) +
  scale_color_manual(values = c("#CC6677","#117733"),
                     name = "Irrigation treatment",
                     labels = c("Control", "Drought")) +
  scale_shape_manual(values = c(15, 3, 16, 8, 17),
                     name = "Sampling Date",
                     labels = c("2022-04-28", "2022-06-01","2022-07-05","2022-07-20","2022-09-13")) + theme_classic() +
  scale_fill_manual(values = c("#CC6677","#CC6677","#CC6677","#CC6677","#CC6677", 
                               "#117733","#117733","#117733","#117733","#117733")) +
  geom_mark_ellipse(aes(fill = aoa.meta.bulk$DxD), 
                   expand = 0, linewidth = NA, show.legend = F)  +
  #labs(title = "AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.position = "right",
        legend.title = element_text(size=14),
        legend.text = element_text(size=14)) +
  annotate("text",x=-4.4,y=-3.6,hjust = 0, size = 5,
  label= "Overall reclassification rate: 78.8%")
aoa.cap.BS.DxD.plot

#setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOA_CAP_BS_DxD.tiff",
       aoa.cap.BS.DxD.plot, device = "tiff",
       width = 12, height =9, 
       units= "in", dpi = 600)
############################################################################

# 2. AOB

# run Bray-Curtis beta diversity on bulk soil
aob.asv.rare1k <- as.data.frame(otu_table(aob.rare.1282.seq))
aob.asv.bulk <- aob.asv.rare1k[,1:119]
aob.asv.bulk1 <- aob.asv.bulk[rowSums(aob.asv.bulk)>0,]
aob.bulk_dist_bc <- vegdist(t(aob.asv.bulk1), method = "bray")
aob.bulk_dist_bc
# metadata
aob.meta.bulk <- aob.meta.df.sub[1:119,]
aob.meta.bulk$period2 <- factor(aob.meta.bulk$period, levels = c("Control.D", "Control.R11W", "Control.R1W", 
                                                              "Rainout.D", "Rainout.R11W",  "Rainout.R1W"),
                          labels = c("Drought", "R11", "R1", "Drought", "R11", "R1"))
#aoa.meta.bulk.ed <- aoa.meta.bulk[,-30:-45]
str(aob.meta.bulk)
aob.meta.bulk$SampleID<-factor(aob.meta.bulk$SampleID)
aob.meta.bulk$PlotID<-factor(aob.meta.bulk$PlotID)
aob.meta.bulk$Irrigation<-factor(aob.meta.bulk$Irrigation)
aob.meta.bulk$Block<-factor(aob.meta.bulk$Block)
aob.meta.bulk$x<-factor(aob.meta.bulk$x)
aob.meta.bulk$rep<-factor(aob.meta.bulk$rep)
aob.meta.bulk$rep2<-factor(aob.meta.bulk$rep2)
aob.meta.bulk$DxD<-factor(aob.meta.bulk$DxD)
aob.meta.bulk$period<-factor(aob.meta.bulk$period)
aob.meta.bulk$period2<-factor(aob.meta.bulk$period2)
aob.meta.bulk$var2<-factor(aob.meta.bulk$var2)
aob.meta.bulk$var3<-factor(aob.meta.bulk$var3)

# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:110) {
  cap <- CAPdiscrim(aob.bulk_dist_bc ~ period, data = aob.meta.bulk, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 42 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.BS.DxPer <- CAPdiscrim(aob.bulk_dist_bc ~ period, data = aob.meta.bulk, m = 42, permutations = 999, add = TRUE) # 29.41% 
success <- cbind(data.frame(aob.cap.BS.DxPer$group), data.frame(aob.cap.BS.DxPer$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.BS.DxPer$PCoA)
success <- success[order(success$source), ]
success

aob.cap1.BS.DxPer <- paste("CAP1 (", round((100/sum(aob.cap.BS.DxPer$lda.other$svd^2) * aob.cap.BS.DxPer$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aob.cap2.BS.DxPer <- paste("CAP2 (", round((100/sum(aob.cap.BS.DxPer$lda.other$svd^2) * aob.cap.BS.DxPer$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with ggplot2

aob.cap.BS.DxPer.plot <- ggplot(as.data.frame(aob.cap.BS.DxPer$x), aes(x = aob.cap.BS.DxPer$x[,1], y = aob.cap.BS.DxPer$x[,2])) +
  geom_point(aes(color = aob.meta.bulk$Irrigation, shape = aob.meta.bulk$period2), size = 4) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aob.cap1.BS.DxPer) + ylab(aob.cap2.BS.DxPer) +
  scale_color_manual(values = c("#CC6677","#117733"),
                     name = "Irrigation treatment",
                     labels = c("Control", "Drought")) +
  scale_shape_manual(values = c(15, 16, 3),
                     name = "Sampling Date",
                     labels = c("Drought", "Rewetting 1 week","Rewetting 11 weeks")) + theme_classic() +
  scale_fill_manual(values = c("#CC6677","#CC6677","#CC6677", 
                               "#117733","#117733","#117733")) +
  geom_mark_ellipse(aes(fill = aob.meta.bulk$period), 
                   expand = 0, linewidth = NA, show.legend = F)  +
  #labs(title = "AOA")+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.position = "right",
        legend.title = element_text(size=14),
        legend.text = element_text(size=14)) 
  #annotate("text",x=-4.4,y=-3.6,hjust = 0, size = 5,
  #label= "Control (drought): 37.14%\nControl (R1W): 8.33%\nControl (R11W): 41.67%\nDrought (drought): 33.33%\nDrought (R1W): 16.67%\nDrought (R11W): 16.67%")
aob.cap.BS.DxPer.plot
#setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOB_CAP_BS_DxPer.tiff",
       aob.cap.BS.DxPer.plot, device = "tiff",
       width = 12, height =9, 
       units= "in", dpi = 600)

#___________________________________________________________________________________________________________________

# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:110) {
  cap <- CAPdiscrim(aob.bulk_dist_bc ~ DxD, data = aob.meta.bulk, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 28 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.BS.DxD <- CAPdiscrim(aob.bulk_dist_bc ~ DxD, data = aob.meta.bulk, m = 28, permutations = 999, add = TRUE) # 18.48% 
success <- cbind(data.frame(aob.cap.BS.DxD$group), data.frame(aob.cap.BS.DxD$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.BS.DxD$PCoA)
success <- success[order(success$source), ]
success

aob.cap1.BS.DxD <- paste("CAP1 (", round((100/sum(aob.cap.BS.DxD$lda.other$svd^2) * aob.cap.BS.DxD$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
aob.cap2.BS.DxD <- paste("CAP2 (", round((100/sum(aob.cap.BS.DxD$lda.other$svd^2) * aob.cap.BS.DxD$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with ggplot2

aob.cap.BS.DxD.plot <- ggplot(as.data.frame(aob.cap.BS.DxD$x), aes(x = aob.cap.BS.DxD$x[,1], y = aob.cap.BS.DxD$x[,2])) +
  geom_point(aes(color = aob.meta.bulk$Irrigation, shape = aob.meta.bulk$Date), size = 4) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(aob.cap1.BS.DxD) + ylab(aob.cap2.BS.DxD) +
  scale_color_manual(values = c("#CC6677","#117733"),
                     name = "Irrigation treatment",
                     labels = c("Control", "Drought")) +
  scale_shape_manual(values = c(15, 3, 16, 8, 17),
                     name = "Sampling Date",
                     labels = c("2022-04-28", "2022-06-01","2022-07-05","2022-07-20","2022-09-13")) + theme_classic() +
  scale_fill_manual(values = c("#CC6677","#CC6677","#CC6677","#CC6677","#CC6677", 
                               "#117733","#117733","#117733","#117733","#117733")) +
  geom_mark_ellipse(aes(fill = aob.meta.bulk$DxD), 
                   expand = 0, linewidth = NA, show.legend = F)  +
  #labs(title = "AOB")+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.position = "right",
        legend.title = element_text(size=14),
        legend.text = element_text(size=14)) 
  #annotate("text",x=-4.4,y=-3.6,hjust = 0, size = 5,
  #label= "Overall reclassification rate: 78.8%")
aob.cap.BS.DxD.plot

#setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOB_CAP_BS_DxD.tiff",
       aob.cap.BS.DxD.plot, device = "tiff",
       width = 12, height =9, 
       units= "in", dpi = 600)
############################################################################################################################

# 3. COMAMMOX

# run Bray-Curtis beta diversity on bulk soil
com.rare.asv.df
com.asv.bulk <- com.rare.asv.df[,1:118]
com.asv.bulk1 <- com.asv.bulk[rowSums(com.asv.bulk)>0,]
sort(rowSums(com.asv.bulk1, na.rm = FALSE, dims = 1), decreasing = FALSE)
com.bulk_dist_bc <- vegdist(t(com.asv.bulk1), method = "bray")
com.bulk_dist_bc
# metadata
com.meta.bulk <- com.meta.df[1:118,]
com.meta.bulk$period2 <- factor(com.meta.bulk$period, levels = c("Control.D", "Control.R11W", "Control.R1W", 
                                                              "Rainout.D", "Rainout.R11W",  "Rainout.R1W"),
                          labels = c("Drought", "R11", "R1", "Drought", "R11", "R1"))
str(com.meta.bulk)
com.meta.bulk$SampleID<-factor(com.meta.bulk$SampleID)
com.meta.bulk$PlotID<-factor(com.meta.bulk$PlotID)
com.meta.bulk$Irrigation<-factor(com.meta.bulk$Irrigation)
com.meta.bulk$Block<-factor(com.meta.bulk$Block)
com.meta.bulk$x<-factor(com.meta.bulk$x)
com.meta.bulk$rep<-factor(com.meta.bulk$rep)
com.meta.bulk$rep2<-factor(com.meta.bulk$rep2)
com.meta.bulk$DxD<-factor(com.meta.bulk$DxD)
com.meta.bulk$period<-factor(com.meta.bulk$period)
com.meta.bulk$period2<-factor(com.meta.bulk$period2)
com.meta.bulk$var2<-factor(com.meta.bulk$var2)
com.meta.bulk$var3<-factor(com.meta.bulk$var3)

# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(com.bulk_dist_bc ~ period, data = com.meta.bulk, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 17 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.BS.DxPer <- CAPdiscrim(com.bulk_dist_bc ~ period, data = com.meta.bulk, m = 17, permutations = 999, add = TRUE) # 33.05% 
success <- cbind(data.frame(com.cap.BS.DxPer$group), data.frame(com.cap.BS.DxPer$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.BS.DxPer$PCoA)
success <- success[order(success$source), ]
success

com.cap1.BS.DxPer <- paste("CAP1 (", round((100/sum(com.cap.BS.DxPer$lda.other$svd^2) * com.cap.BS.DxPer$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
com.cap2.BS.DxPer <- paste("CAP2 (", round((100/sum(com.cap.BS.DxPer$lda.other$svd^2) * com.cap.BS.DxPer$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with ggplot2

com.cap.BS.DxPer.plot <- ggplot(as.data.frame(com.cap.BS.DxPer$x), aes(x = com.cap.BS.DxPer$x[,1], y = com.cap.BS.DxPer$x[,2])) +
  geom_point(aes(color = com.meta.bulk$Irrigation, shape = com.meta.bulk$period2), size = 4) +
  #geom_text(label=com.meta.bulk.ed$PlotID)+
  xlab(com.cap1.BS.DxPer) + ylab(com.cap2.BS.DxPer) +
  scale_color_manual(values = c("#CC6677","#117733"),
                     name = "Irrigation treatment",
                     labels = c("Control", "Drought")) +
  scale_shape_manual(values = c(15, 16, 3),
                     name = "Sampling Date",
                     labels = c("Drought", "Rewetting 1 week","Rewetting 11 weeks")) + theme_classic() +
  scale_fill_manual(values = c("#CC6677","#CC6677","#CC6677", 
                               "#117733","#117733","#117733")) +
  geom_mark_ellipse(aes(fill = com.meta.bulk$period), 
                   expand = 0, linewidth = NA, show.legend = F)  +
  #labs(title = "COMAMMOX")+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.position = "right",
        legend.title = element_text(size=14),
        legend.text = element_text(size=14)) 
  #annotate("text",x=-4.4,y=-3.6,hjust = 0, size = 5,
  #label= "Control (drought): 37.14%\nControl (R1W): 8.33%\nControl (R11W): 41.67%\nDrought (drought): 33.33%\nDrought (R1W): 16.67%\nDrought (R11W): 16.67%")
com.cap.BS.DxPer.plot
#setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("COMA_CAP_BS_DxPer.tiff",
       com.cap.BS.DxPer.plot, device = "tiff",
       width = 12, height =9, 
       units= "in", dpi = 600)

#___________________________________________________________________________________________________________________

# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:110) {
  cap <- CAPdiscrim(com.bulk_dist_bc ~ DxD, data = com.meta.bulk, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 100 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.BS.DxD <- CAPdiscrim(com.bulk_dist_bc ~ DxD, data = com.meta.bulk, m = 100, permutations = 999, add = TRUE) # 7.6% 
success <- cbind(data.frame(com.cap.BS.DxD$group), data.frame(com.cap.BS.DxD$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.BS.DxD$PCoA)
success <- success[order(success$source), ]
success

com.cap1.BS.DxD <- paste("CAP1 (", round((100/sum(com.cap.BS.DxD$lda.other$svd^2) * com.cap.BS.DxD$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
com.cap2.BS.DxD <- paste("CAP2 (", round((100/sum(com.cap.BS.DxD$lda.other$svd^2) *com.cap.BS.DxD$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")

# Plot with ggplot2
com.cap.BS.DxD.plot <- ggplot(as.data.frame(com.cap.BS.DxD$x), aes(x = com.cap.BS.DxD$x[,1], y = com.cap.BS.DxD$x[,2])) +
  geom_point(aes(color = com.meta.bulk$Irrigation, shape = com.meta.bulk$Date), size = 4) +
  #geom_text(label=aoa.meta.bulk.ed$PlotID)+
  xlab(com.cap1.BS.DxD) + ylab(com.cap2.BS.DxD) +
  scale_color_manual(values = c("#CC6677","#117733"),
                     name = "Irrigation treatment",
                     labels = c("Control", "Drought")) +
  scale_shape_manual(values = c(15, 3, 16, 8, 17),
                     name = "Sampling Date",
                     labels = c("2022-04-28", "2022-06-01","2022-07-05","2022-07-20","2022-09-13")) + theme_classic() +
  scale_fill_manual(values = c("#CC6677","#CC6677","#CC6677","#CC6677","#CC6677", 
                               "#117733","#117733","#117733","#117733","#117733")) +
  geom_mark_ellipse(aes(fill = com.meta.bulk$DxD), 
                   expand = 0, linewidth = NA, show.legend = F)  +
  #labs(title = "COMAMMOX")+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(legend.position = "right",
        legend.title = element_text(size=14),
        legend.text = element_text(size=14)) 
  #annotate("text",x=-4.4,y=-3.6,hjust = 0, size = 5,
  #label= "Overall reclassification rate: 78.8%")
com.cap.BS.DxD.plot

#setwd('D:/Fina/INRAE_Project/microservices_fig/')
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("COMA_CAP_BS_DxD.tiff",
       com.cap.BS.DxD.plot, device = "tiff",
       width = 12, height =9, 
       units= "in", dpi = 600)

