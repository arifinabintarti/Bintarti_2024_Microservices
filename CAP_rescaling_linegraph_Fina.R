library(BiodiversityR) 
### 1. AOA

### 1 A. AOA Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
aoa.bulk_dist_bc <- vegdist(t(aoa.asv.bulk1), method = "bray")
aoa.bulk_dist_bc
# metadata
aoa.meta.bulk
aoa.meta.bulk.ed <- aoa.meta.bulk[,-30:-45]
aoa.meta.bulk.ed$x <- as.factor(aoa.meta.bulk.ed$x)
aoa.meta.bulk.ed$Block <- as.factor(aoa.meta.bulk.ed$Block)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aoa.bulk_dist_bc ~ Irrigation, data = aoa.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 24 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.bulk.irri <- CAPdiscrim(aoa.bulk_dist_bc ~ Irrigation, data = aoa.meta.bulk.ed, m = 24, permutations = 9999, add = TRUE) # 79.17% 
aoa.cap.bulk.irri 

success <- cbind(data.frame(aoa.cap.bulk.irri$group), data.frame(aoa.cap.bulk.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.bulk.irri$PCoA)
success <- success[order(success$source), ]
success
#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "AOA.cap.sucess.irri.BS.csv")

aoa.cap1.bulk <- paste("CAP1 (", round((100/sum(aoa.cap.bulk.irri$lda.other$svd^2) * aoa.cap.bulk.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
aoa.cap2.bulk <- paste("CAP2 (", round((100/sum(aoa.cap.bulk.irri$lda.other$svd^2) * aoa.cap.bulk.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(aoa.cap.bulk.irri$x), aes(x = aoa.cap.bulk.irri$x[,1])) +
  geom_density(aes(fill = aoa.meta.bulk.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(aoa.cap1.bulk) 

# point plot ==============================

ggplot(as.data.frame(aoa.cap.bulk.irri$x), aes(x = aoa.cap.bulk.irri$x[,1], y = aoa.cap.bulk.irri$x[,2])) +
  geom_point(aes(shape = factor(aoa.meta.bulk.ed$Irrigation), color = aoa.meta.bulk.ed$Treatment), size = 2) +
  xlab(aoa.cap1.bulk) + ylab(aoa.cap2.bulk)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

cap.B.points <- aoa.cap.bulk.irri$x[,1]

cap.points.B.df <- cap.B.points %>% as.data.frame() %>%
  group_by(aoa.meta.bulk.ed$Irrigation, aoa.meta.bulk.ed$Treatment, aoa.meta.bulk.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'aoa.meta.bulk.ed$Irrigation', Treatment = 'aoa.meta.bulk.ed$Treatment', Date = 'aoa.meta.bulk.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

cap.points.B.df <- mutate(cap.points.B.df, x = paste(Irrigation, Treatment)) 
str(cap.points.B.df)

cap.points.B.df$Date <- factor(cap.points.B.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
cap.points.B.df$Date <- as.Date(cap.points.B.df$Date)
cap.points.B.df$Date
cap.points.B.df$x <- as.factor(cap.points.B.df$x)

# point point of mean cap ==============================

#cap.points.B <- ggplot(cap.points.B.df, aes(x = Date , y = mean, color = Treatment, shape = Irrigation)) +
  #theme_classic() +
  #geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd), 
                  #position=position_jitter(width=1.5))+
  #xlab("") +
  #ylab("Mean of CAP1") +
  #scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     #name = "Cropping system",
                     #labels = c("BIODYN", "CONFYM", "CONMIN")) +
  #theme(legend.position = "bottom")  +
  #annotate(geom = "text", x = as.Date("2022-05-01"), y = -3.4, hjust = 0, size = 3, label = "100%", color = "black") +
  #annotate(geom = "text", x = as.Date("2022-06-04"), y = -2.3, hjust = 0, size = 3, label = "100%", color = "black") +
  #annotate(geom = "text", x = as.Date("2022-07-04"), y = -4.2, hjust = 0, size = 3, label = "100%", color = "black") +
  #annotate(geom = "text", x = as.Date("2022-07-12"), y = -2.3, hjust = 0, size = 3, label = "100%", color = "black") +
  #annotate(geom = "text", x = as.Date("2022-09-06"), y = -2.7, hjust = 0, size = 3, label = "100%", color = "black") +
  #annotate(geom = "text", x = as.Date("2022-05-02"), y = 3, hjust = 0, size = 3, label = "100%", color = "black") +
  #annotate(geom = "text", x = as.Date("2022-06-04"), y = 2.3, hjust = 0, size = 3, label = "100%", color = "black") +
  #annotate(geom = "text", x = as.Date("2022-06-27"), y = 3.6, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  #annotate(geom = "text", x = as.Date("2022-07-07"), y = 3.8, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  #annotate(geom = "text", x = as.Date("2022-07-08"), y = 0.8, hjust = 0, size = 3, label = "75%", color = "#E69F00") +
  #annotate(geom = "text", x = as.Date("2022-07-17"), y = 1.3, hjust = 0, size = 3, label = "100%", color = "black") +
  #annotate(geom = "text", x = as.Date("2022-09-04"), y = 2.6, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  #annotate(geom = "text", x = as.Date("2022-09-17"), y = 2, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  #annotate(geom = "text", x = as.Date("2022-09-06"), y = 3.5, hjust = 0, size = 3, label = "100%", color = "#E69F00") +
  #annotate(geom = "text", x = as.Date("2022-04-28"), y =-5.5, hjust = 0, size = 3, label = "Reclassification rate drought-induced: 97%", color = "black") +
  #annotate(geom = "text", x = as.Date("2022-04-28"), y =-6.5, hjust = 0, size = 3, label = "Reclassification rate control: 100%", color = "black")
#cap.points.B

# line point of mean cap ==============================

AOA.irri.cap.line <- ggplot(cap.points.B.df, aes(x = Date , y = mean, color = Treatment, group = x)) +
  geom_point(size = 1.5) + theme_classic() +
  geom_line(aes(linetype = Irrigation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Treatment), linetype=0, alpha=0.05) +
  xlab("") +
  ylab("Mean of CAP1") +
  labs(title = "A. Bulk Soil")+
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00"),
                    name = "Cropping system",
                    labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_linetype_manual(name = "Irrigation treatment",
                        labels = c("control", "drought-induced"),
                        values = c("solid", "dashed")) +
  scale_x_date(date_breaks = "1 month",date_labels = "%B")+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 13, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 13, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(size = 20, face="bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.background = element_blank()) +
  theme(legend.position = "none", legend.box="horizontal") +
  #theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = -1.3, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = -2, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = -0.98, hjust = 0, size = 3, label = "75%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-05-30"), y = -1.1, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-05-30"), y = -1.4, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-05-30"), y = -0.3, hjust = 0, size = 3, label = "50%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-07-03"), y = -1.2, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = -1.8, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = -0.55, hjust = 0, size = 3, label = "50%", color = "#E69F00") +
 
  annotate(geom = "text", x = as.Date("2022-07-18"), y = -0.5, hjust = 0, size = 3, label = "75%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-07-18"), y = -1.8, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-07-18"), y = -0.1, hjust = 0, size = 3, label = "25%", color = "#E69F00") +
 
  annotate(geom = "text", x = as.Date("2022-09-13"), y = -1.9, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-09-13"), y = -1, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-09-13"), y = -1.3, hjust = 0, size = 3, label = "75%", color = "#E69F00") +
 
  annotate(geom = "text", x = as.Date("2022-04-24"), y = 0.9, hjust = 0, size = 3, label = "75%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = 1.45, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = 0.3, hjust = 0, size = 3, label = "50%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-05-30"), y = 2, hjust = 0, size = 3, label = "75%", color = "black") +
  
  annotate(geom = "text", x = as.Date("2022-07-03"), y = 0.4, hjust = 0, size = 3, label = "50%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = 2, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = 1.2, hjust = 0, size = 3, label = "75%", color = "#E69F00") +
 
  annotate(geom = "text", x = as.Date("2022-07-18"), y = 1.8, hjust = 0, size = 3, label = "100%", color = "black") +
 
  annotate(geom = "text", x = as.Date("2022-09-13"), y = 0.6, hjust = 0, size = 3, label = "75%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-09-13"), y = 1.9, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-09-13"), y = 1.3, hjust = 0, size = 3, label = "50%", color = "#E69F00") +

  annotate(geom = "text", x = as.Date("2022-04-28"), y =-5, hjust = 0, size = 4, label = "Reclassification rate drought-induced: 76.67%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-04-28"), y =-5.4, hjust = 0, size = 4, label = "Reclassification rate control: 81.67%", color = "black") +
  geom_vline(xintercept = as.Date("2022-07-14"), linetype="dashed", color = "black", linewidth = 0.5) +
  annotate(geom = "text", x = as.Date("2022-07-10"), y = -5, hjust = 0, size = 4, label = "Rewetting", color = "black")

AOA.irri.cap.line

setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("aoa.resscalling.irri.tiff",
       AOA.irri.cap.line, device = "tiff",
       width = 9.5, height = 6, 
       units= "in", dpi = 600)

### 1 B. AOA Rhizosphere

# run Bray-Curtis beta diversity on bulk soil
aoa.rh_dist_bc <- vegdist(t(aoa.asv.rh1), method = "bray")
aoa.rh_dist_bc 
# metadata
aoa.meta.rh <- aoa.meta.df[121:192,]
aoa.meta.rh.ed <- aoa.meta.rh[,c(-14:-29,-42:-45)]
aoa.meta.rh.ed$x <- as.factor(aoa.meta.rh.ed$x)
aoa.meta.rh.ed$Block <- as.factor(aoa.meta.rh.ed$Block)
#aoa.meta.rh.ed$Date <- factor(aoa.meta.rh.ed$Date)
str(aoa.meta.rh.ed)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aoa.rh_dist_bc ~ Irrigation, data = aoa.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 22 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.rh.irri <- CAPdiscrim(aoa.rh_dist_bc ~ Irrigation, data = aoa.meta.rh.ed, m = 22, permutations = 9999, add = TRUE) # 86.11% 
aoa.cap.rh.irri 

success <- cbind(data.frame(aoa.cap.rh.irri$group), data.frame(aoa.cap.rh.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.rh.irri$PCoA)
success <- success[order(success$source), ]
success
#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "AOA.cap.sucess.irri.RS.csv")

aoa.cap1.rh <- paste("CAP1 (", round((100/sum(aoa.cap.rh.irri$lda.other$svd^2) * aoa.cap.rh.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
aoa.cap2.rh <- paste("CAP2 (", round((100/sum(aoa.cap.rh.irri$lda.other$svd^2) * aoa.cap.rh.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(aoa.cap.rh.irri$x), aes(x = aoa.cap.rh.irri$x[,1])) +
  geom_density(aes(fill =  aoa.meta.rh.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(aoa.cap1.rh) 

# point plot ==============================

ggplot(as.data.frame(aoa.cap.rh.irri$x), aes(x = aoa.cap.rh.irri$x[,1], y = aoa.cap.rh.irri$x[,2])) +
  geom_point(aes(shape = factor(aoa.meta.rh.ed$Irrigation), color = aoa.meta.rh.ed$Treatment), size = 2) +
  xlab(aoa.cap1.rh) + ylab(aoa.cap2.rh)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

AOA.cap.rh.points <- aoa.cap.rh.irri$x[,1]

AOA.cap.points.rh.df <- AOA.cap.rh.points %>% as.data.frame() %>%
  group_by(aoa.meta.rh.ed$Irrigation, aoa.meta.rh.ed$Treatment, aoa.meta.rh.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'aoa.meta.rh.ed$Irrigation', Treatment = 'aoa.meta.rh.ed$Treatment', Date = 'aoa.meta.rh.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

AOA.cap.points.rh.df <- mutate(AOA.cap.points.rh.df, x = paste(Irrigation, Treatment)) 
str(AOA.cap.points.rh.df)

AOA.cap.points.rh.df$Date <- factor(AOA.cap.points.rh.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05"))
AOA.cap.points.rh.df$Date <- as.Date(AOA.cap.points.rh.df$Date)
AOA.cap.points.rh.df$Date
AOA.cap.points.rh.df$x <- as.factor(AOA.cap.points.rh.df$x)


# line point of mean cap ==============================

AOA.irri.cap.line.rh <- ggplot(AOA.cap.points.rh.df, aes(x = Date , y = mean, color = Treatment, group = x)) +
  geom_point(size = 1.5) + 
  theme_classic() +
  geom_line(aes(linetype = Irrigation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Treatment), linetype=0, alpha=0.05) +
  xlab("") +
  ylab("Mean of CAP1") +
  labs(title = "B. Rhizosphere")+
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00"),
                    name = "Cropping system",
                    labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_linetype_manual(name = "Irrigation treatment",
                        labels = c("control", "drought-induced"),
                        values = c("solid", "dashed")) +
  scale_x_date(date_breaks = "1 month",date_labels = "%B")+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 13, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 13, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(size = 20, face="bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.background = element_blank()) +
  theme(legend.position = "none", legend.box="horizontal") +
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = -1.5, hjust = 0, size = 3, label = "75%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = -2.2, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = -1, hjust = 0, size = 3, label = "100%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-05-30"), y = -1, hjust = 0, size = 3, label = "75%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-05-30"), y = -1.5, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-05-30"), y = -2, hjust = 0, size = 3, label = "100%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-07-03"), y = -2.3, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = -1.1, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = -2.1, hjust = 0, size = 3, label = "100%", color = "#E69F00") +
 
  annotate(geom = "text", x = as.Date("2022-04-24"), y = 0.7, hjust = 0, size = 3, label = "75%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = 1.1, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-04-24"), y = 2, hjust = 0, size = 3, label = "100%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-05-30"), y = 2.9, hjust = 0, size = 3, label = "100%", color = "black") +
  
  annotate(geom = "text", x = as.Date("2022-07-03"), y = 1.2, hjust = 0, size = 3, label = "75%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = 2, hjust = 0, size = 3, label = "50%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = 1, hjust = 0, size = 3, label = "75%", color = "#E69F00") +

  annotate(geom = "text", x = as.Date("2022-04-28"), y =-5, hjust = 0, size = 4, label = "Reclassification rate drought-induced: 83.33%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-04-28"), y =-5.4, hjust = 0, size = 4, label = "Reclassification rate control: 88.88%", color = "black") 
  #coord_fixed()
  #geom_vline(xintercept = as.Date("2022-07-14"), linetype="dashed", color = "black", linewidth = 0.5) +
  #annotate(geom = "text", x = as.Date("2022-07-10"), y = -5, hjust = 0, size = 4, label = "Rewetting", color = "black")

AOA.irri.cap.line.rh

setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("aoa.resscalling.irri.rh.tiff",
       AOA.irri.cap.line.rh, device = "tiff",
       width = 9.5, height = 6, 
       units= "in", dpi = 600)

# Combine between AOA Bulk Soil and Rhizosphere
library(patchwork)
library(gridExtra)
#AOA.irri.cap.line+ AOA.irri.cap.line.rh +
 #plot_layout(design ="1#
 #2#") 
#AOA.CAP.line.irri.ALL <- AOA.irri.cap.line / AOA.irri.cap.line.rh
#AOA.CAP.line.irri.ALL
# layout
#layout <- matrix(c(1, 1,
                   #0, 2), ncol = 2, byrow = T)
#grid.arrange(AOA.irri.cap.line.rh,AOA.irri.cap.line ,
             #layout_matrix = layout) 

# combine
cap.points.B.df2 <- cap.points.B.df %>%
  mutate(Type = "Bulk Soil")
AOA.cap.points.rh.df2 <- AOA.cap.points.rh.df %>%
  mutate(Type = "Rhizossphere")
AOA.cap.irri.all <- merge(cap.points.B.df2,AOA.cap.points.rh.df2,all = T)


# plot
AOA.cap.irri.all.line <- ggplot(AOA.cap.irri.all, aes(x = Date , y = mean, color = Treatment, group = x)) +
  geom_point(size = 1.5) +
  theme_bw() +
  #theme_classic() +
  geom_line(aes(linetype = Irrigation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Treatment), linetype=0, alpha=0.05) +
  xlab("") +
  ylab("Mean of CAP1") +
  facet_wrap(~ Type, strip.position="right", nrow = 2)+
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00"),
                    name = "Cropping system",
                    labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_linetype_manual(name = "Treatment",
                        labels = c("control", "drought"),
                        values = c("solid", "dashed")) +
  scale_x_date(date_breaks = "1 month",date_labels = "%B")+
  theme(axis.text.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 18, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(size = 20, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.background = element_blank()) +
  theme(legend.position = "bottom", legend.box="horizontal")+
  geom_vline(data=filter(AOA.cap.irri.all, Type=="Bulk Soil"), aes(xintercept = as.Date("2022-07-14")), linetype="dashed",color = "black", linewidth = 0.5) +
  geom_text(data=filter(AOA.cap.irri.all, Type=="Bulk Soil"), aes(x = as.Date("2022-07-10"), y = -5), hjust = 0, size = 7, label = "Rewetting", color = "black")
AOA.cap.irri.all.line

library(cowplot)
AOA.cap.irri.all.lineplot <- ggdraw(AOA.cap.irri.all.line) + draw_plot_label(x = c(0.08,0.08), y = c(.67,.24), 
                       label =  c("Reclassification rate drought-induced: 76.67%\nReclassification rate control: 81.67%",  "Reclassification rate drought-induced: 83.33%\nReclassification rate control: 88.88%"),
                       hjust = 0, size = 14)
AOA.cap.irri.all.lineplot 
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOA.cap.irri.all.lineplot.tiff",
       AOA.cap.irri.all.lineplot, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 600)

### 2. AOB

### 2 A. AOB Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
aob.bulk_dist_bc <- vegdist(t(aob.asv.bulk1), method = "bray")
aob.bulk_dist_bc
# metadata
# metadata
aob.meta.bulk
aob.meta.bulk.ed <- aob.meta.bulk[,-30:-45]
aob.meta.bulk.ed$x <- as.factor(aob.meta.bulk.ed$x)
aob.meta.bulk.ed$Block <- as.factor(aob.meta.bulk.ed$Block)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(aob.bulk_dist_bc ~ Irrigation, data = aob.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 68 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.bulk.irri <- CAPdiscrim(aob.bulk_dist_bc ~ Irrigation, data = aob.meta.bulk.ed, m = 68, permutations = 9999, add = TRUE) # 69.747% 
aob.cap.bulk.irri  #Significance of this percentage was 0.00020002 
#Control (n=59) correct: 66.1016949152542 percent
#Rainout (n=60) correct: 73.3333333333333 percent

success <- cbind(data.frame(aoa.cap.bulk.irri$group), data.frame(aoa.cap.bulk.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.bulk.irri$PCoA)
success <- success[order(success$source), ]
success
#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "AOB.cap.sucess.irri.BS.csv")

aob.cap1.bulk <- paste("CAP1 (", round((100/sum(aob.cap.bulk.irri$lda.other$svd^2) * aob.cap.bulk.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
aob.cap2.bulk <- paste("CAP2 (", round((100/sum(aob.cap.bulk.irri$lda.other$svd^2) * aob.cap.bulk.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(aob.cap.bulk.irri$x), aes(x = aob.cap.bulk.irri$x[,1])) +
  geom_density(aes(fill = aob.meta.bulk.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(aob.cap1.bulk) 

# point plot ==============================

ggplot(as.data.frame(aob.cap.bulk.irri$x), aes(x = aob.cap.bulk.irri$x[,1], y = aob.cap.bulk.irri$x[,2])) +
  geom_point(aes(shape = factor(aob.meta.bulk.ed$Irrigation), color = aob.meta.bulk.ed$Treatment), size = 2) +
  xlab(aob.cap1.bulk) + ylab(aob.cap2.bulk)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

AOB.cap.BS.points <- aob.cap.bulk.irri$x[,1]

AOB.cap.BS.points.df <- AOB.cap.BS.points %>% as.data.frame() %>%
  group_by(aob.meta.bulk.ed$Irrigation, aob.meta.bulk.ed$Treatment, aob.meta.bulk.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'aob.meta.bulk.ed$Irrigation', Treatment = 'aob.meta.bulk.ed$Treatment', Date = 'aob.meta.bulk.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

AOB.cap.BS.points.df <- mutate(AOB.cap.BS.points.df, x = paste(Irrigation, Treatment)) 
str(AOB.cap.BS.points.df)

AOB.cap.BS.points.df$Date <- factor(AOB.cap.BS.points.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
AOB.cap.BS.points.df$Date <- as.Date(AOB.cap.BS.points.df$Date)
AOB.cap.BS.points.df$Date
AOB.cap.BS.points.df$x <- as.factor(AOB.cap.BS.points.df$x)

### 2 B. AOB Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
aob.rh_dist_bc <- vegdist(t(aob.asv.rh1), method = "bray")
aob.rh_dist_bc
# metadata
aob.meta.rh
aob.meta.rh.ed <- aob.meta.rh[,c(-14:-29,-42:-45)]
aob.meta.rh.ed$x <- as.factor(aob.meta.rh.ed$x)
str(aob.meta.rh.ed)
aob.meta.rh.ed$Date <- factor(aob.meta.rh.ed$Date)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aob.rh_dist_bc ~ Irrigation, data = aob.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 14 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.rh.irri <- CAPdiscrim(aob.rh_dist_bc ~ Irrigation, data = aob.meta.rh.ed, m = 14 , permutations = 9999, add = TRUE) # 69.44 % 
aob.cap.rh.irri #Significance of this percentage was 0.00260026 
#Control (n=36) correct: 72.2222222222222 percent
#Rainout (n=36) correct: 66.6666666666667 percent

success <- cbind(data.frame(aob.cap.rh.irri$group), data.frame(aob.cap.rh.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.rh.irri$PCoA)
success <- success[order(success$source), ]
success

#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "AOB.cap.sucess.irri.RS.csv")

aob.cap1.rh <- paste("CAP1 (", round((100/sum(aob.cap.rh.irri$lda.other$svd^2) * aob.cap.rh.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
aob.cap2.rh <- paste("CAP2 (", round((100/sum(aob.cap.rh.irri$lda.other$svd^2) * aob.cap.rh.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(aob.cap.rh.irri$x), aes(x = aob.cap.rh.irri$x[,1])) +
  geom_density(aes(fill = aob.meta.rh.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(aob.cap1.rh) 

# point plot ==============================

ggplot(as.data.frame(aob.cap.rh.irri$x), aes(x = aob.cap.rh.irri$x[,1], y = aob.cap.rh.irri$x[,2])) +
  geom_point(aes(shape = factor(aob.meta.rh.ed$Irrigation), color = aob.meta.rh.ed$Treatment), size = 2) +
  xlab(aob.cap1.rh) + ylab(aob.cap2.rh)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

AOB.cap.RS.points <- aob.cap.rh.irri$x[,1]

AOB.cap.RS.points.df <- AOB.cap.RS.points %>% as.data.frame() %>%
  group_by(aob.meta.rh.ed$Irrigation, aob.meta.rh.ed$Treatment, aob.meta.rh.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'aob.meta.rh.ed$Irrigation', Treatment = 'aob.meta.rh.ed$Treatment', Date = 'aob.meta.rh.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

AOB.cap.RS.points.df <- mutate(AOB.cap.RS.points.df, x = paste(Irrigation, Treatment)) 
str(AOB.cap.RS.points.df)

AOB.cap.RS.points.df$Date <- factor(AOB.cap.RS.points.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05"))
AOB.cap.RS.points.df$Date <- as.Date(AOB.cap.RS.points.df$Date)
AOB.cap.RS.points.df$Date
AOB.cap.RS.points.df$x <- as.factor(AOB.cap.RS.points.df$x)

# combine
AOB.cap.BS.points.df2 <- AOB.cap.BS.points.df %>%
  mutate(Type = "Bulk Soil")
AOB.cap.RS.points.df2 <- AOB.cap.RS.points.df %>%
  mutate(Type = "Rhizossphere")
AOB.cap.irri.all <- merge(AOB.cap.BS.points.df2,AOB.cap.RS.points.df2,all = T)

# plot
AOB.cap.irri.all.line <- ggplot(AOB.cap.irri.all, aes(x = Date , y = mean, color = Treatment, group = x)) +
  geom_point(size = 1.5) +
  theme_bw() +
  #theme_classic() +
  geom_line(aes(linetype = Irrigation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Treatment), linetype=0, alpha=0.05) +
  xlab("") +
  ylab("Mean of CAP1") +
  facet_wrap(~ Type, strip.position="right", nrow = 2)+
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00"),
                    name = "Cropping system",
                    labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_linetype_manual(name = "Treatment",
                        labels = c("control", "drought"),
                        values = c("solid", "dashed")) +
  scale_x_date(date_breaks = "1 month",date_labels = "%B")+
  theme(axis.text.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 18, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(size = 20, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.background = element_blank()) +
  theme(legend.position = "bottom", legend.box="horizontal")+
  geom_vline(data=filter(AOB.cap.irri.all, Type=="Bulk Soil"), aes(xintercept = as.Date("2022-07-14")), linetype="dashed",color = "black", linewidth = 0.5) +
  geom_text(data=filter(AOB.cap.irri.all, Type=="Bulk Soil"), aes(x = as.Date("2022-07-10"), y = -4), hjust = 0, size = 7, label = "Rewetting", color = "black")
AOB.cap.irri.all.line

library(cowplot)
AOB.cap.irri.all.lineplot <- ggdraw(AOB.cap.irri.all.line) + draw_plot_label(x = c(0.08,0.08), y = c(.67,.25), 
                       label =  c("Reclassification rate drought-induced: 73.33%\nReclassification rate control: 66.1%",  "Reclassification rate drought-induced: 66.67%\nReclassification rate control: 72.22%"),
                       hjust = 0, size = 14)
AOB.cap.irri.all.lineplot 
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOB.cap.irri.all.lineplot.tiff",
       AOB.cap.irri.all.lineplot, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 600)

### 3. COMAMMOX

### 3 A. COMA Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
com.bulk_dist_bc <- vegdist(t(com.asv.bulk1), method = "bray")
com.bulk_dist_bc
# metadata
# metadata
com.meta.bulk
com.meta.bulk.ed <- com.meta.bulk[,-30:-45]
com.meta.bulk.ed$x <- as.factor(com.meta.bulk.ed$x)
com.meta.bulk.ed$Block <- as.factor(com.meta.bulk.ed$Block)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(com.bulk_dist_bc ~ Irrigation, data = com.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 25 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.bulk.irri <- CAPdiscrim(com.bulk_dist_bc ~ Irrigation, data = com.meta.bulk.ed, m = 25, permutations = 9999, add = TRUE) # 72.03% 
com.cap.bulk.irri  #Significance of this percentage was 0.00010001 
#Control (n=60) correct: 76.6666666666667 percent
#Rainout (n=58) correct: 67.2413793103448 percent

success <- cbind(data.frame(com.cap.bulk.irri$group), data.frame(com.cap.bulk.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.bulk.irri$PCoA)
success <- success[order(success$source), ]
success
#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "COM.cap.sucess.irri.BS.csv")

com.cap1.bulk <- paste("CAP1 (", round((100/sum(com.cap.bulk.irri$lda.other$svd^2) * com.cap.bulk.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
com.cap2.bulk <- paste("CAP2 (", round((100/sum(com.cap.bulk.irri$lda.other$svd^2) * com.cap.bulk.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(com.cap.bulk.irri$x), aes(x = com.cap.bulk.irri$x[,1])) +
  geom_density(aes(fill = com.meta.bulk.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(com.cap1.bulk) 

# point plot ==============================

ggplot(as.data.frame(com.cap.bulk.irri$x), aes(x = com.cap.bulk.irri$x[,1], y = com.cap.bulk.irri$x[,2])) +
  geom_point(aes(shape = factor(com.meta.bulk.ed$Irrigation), color = com.meta.bulk.ed$Treatment), size = 2) +
  xlab(com.cap1.bulk) + ylab(com.cap2.bulk)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

COM.cap.BS.points <- com.cap.bulk.irri$x[,1]

COM.cap.BS.points.df <- COM.cap.BS.points %>% as.data.frame() %>%
  group_by(com.meta.bulk.ed$Irrigation, com.meta.bulk.ed$Treatment, com.meta.bulk.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'com.meta.bulk.ed$Irrigation', Treatment = 'com.meta.bulk.ed$Treatment', Date = 'com.meta.bulk.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

COM.cap.BS.points.df <- mutate(COM.cap.BS.points.df, x = paste(Irrigation, Treatment)) 
str(COM.cap.BS.points.df)

COM.cap.BS.points.df$Date <- factor(COM.cap.BS.points.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
COM.cap.BS.points.df$Date <- as.Date(COM.cap.BS.points.df$Date)
COM.cap.BS.points.df$Date
COM.cap.BS.points.df$x <- as.factor(COM.cap.BS.points.df$x)

### 3 B. COMAMMOX Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
com.rh_dist_bc <- vegdist(t(com.asv.rh1), method = "bray")
com.rh_dist_bc
# metadata
com.meta.rh
com.meta.rh.ed <- com.meta.rh[,c(-14:-29,-42:-45)]
com.meta.rh.ed$x <- as.factor(com.meta.rh.ed$x)
com.meta.rh.ed$Date <- factor(com.meta.rh.ed$Date)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(com.rh_dist_bc ~ Irrigation, data = com.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 7 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.rh.irri <- CAPdiscrim(com.rh_dist_bc ~ Irrigation, data = com.meta.rh.ed, m = 7 , permutations = 9999, add = TRUE) #  77.78% 
com.cap.rh.irri #Significance of this percentage was 0.00010001  
#Control (n=36) correct: 77.78 percent
#Rainout (n=36) correct: 77.78 percent

success <- cbind(data.frame(com.cap.rh.irri$group), data.frame(com.cap.rh.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.rh.irri$PCoA)
success <- success[order(success$source), ]
success

#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "COM.cap.sucess.irri.RS.csv")

com.cap1.rh <- paste("CAP1 (", round((100/sum(com.cap.rh.irri$lda.other$svd^2) * com.cap.rh.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
com.cap2.rh <- paste("CAP2 (", round((100/sum(com.cap.rh.irri$lda.other$svd^2) * com.cap.rh.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(com.cap.rh.irri$x), aes(x = com.cap.rh.irri$x[,1])) +
  geom_density(aes(fill = com.meta.rh.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(com.cap1.rh) 

# point plot ==============================

ggplot(as.data.frame(com.cap.rh.irri$x), aes(x = com.cap.rh.irri$x[,1], y = com.cap.rh.irri$x[,2])) +
  geom_point(aes(shape = factor(com.meta.rh.ed$Irrigation), color = com.meta.rh.ed$Treatment), size = 2) +
  xlab(com.cap1.rh) + ylab(com.cap2.rh)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

COM.cap.RS.points <- com.cap.rh.irri$x[,1]

COM.cap.RS.points.df <- COM.cap.RS.points %>% as.data.frame() %>%
  group_by(com.meta.rh.ed$Irrigation, com.meta.rh.ed$Treatment, com.meta.rh.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'com.meta.rh.ed$Irrigation', Treatment = 'com.meta.rh.ed$Treatment', Date = 'com.meta.rh.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

COM.cap.RS.points.df <- mutate(COM.cap.RS.points.df, x = paste(Irrigation, Treatment)) 
str(COM.cap.RS.points.df)

COM.cap.RS.points.df$Date <- factor(COM.cap.RS.points.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05"))
COM.cap.RS.points.df$Date <- as.Date(COM.cap.RS.points.df$Date)
COM.cap.RS.points.df$Date
COM.cap.RS.points.df$x <- as.factor(COM.cap.RS.points.df$x)

# combine
COM.cap.BS.points.df2 <- COM.cap.BS.points.df %>%
  mutate(Type = "Bulk Soil")
COM.cap.RS.points.df2 <- COM.cap.RS.points.df %>%
  mutate(Type = "Rhizossphere")
COM.cap.irri.all <- merge(COM.cap.BS.points.df2,COM.cap.RS.points.df2,all = T)

# plot
COM.cap.irri.all.line <- ggplot(COM.cap.irri.all, aes(x = Date , y = mean, color = Treatment, group = x)) +
  geom_point(size = 1.5) +
  theme_bw() +
  #theme_classic() +
  geom_line(aes(linetype = Irrigation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Treatment), linetype=0, alpha=0.05) +
  xlab("") +
  ylab("Mean of CAP1") +
  facet_wrap(~ Type, strip.position="right", nrow = 2)+
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00"),
                    name = "Cropping system",
                    labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_linetype_manual(name = "Treatment",
                        labels = c("control", "drought"),
                        values = c("solid", "dashed")) +
  scale_x_date(date_breaks = "1 month",date_labels = "%B")+
  theme(axis.text.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 18, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(size = 20, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.background = element_blank()) +
  theme(legend.position = "bottom", legend.box="horizontal")+
  geom_vline(data=filter(COM.cap.irri.all, Type=="Bulk Soil"), aes(xintercept = as.Date("2022-07-14")), linetype="dashed",color = "black", linewidth = 0.5) +
  geom_text(data=filter(COM.cap.irri.all, Type=="Bulk Soil"), aes(x = as.Date("2022-07-10"), y = -4), hjust = 0, size = 7, label = "Rewetting", color = "black")
COM.cap.irri.all.line

library(cowplot)
COM.cap.irri.all.lineplot <- ggdraw(COM.cap.irri.all.line) + draw_plot_label(x = c(0.08,0.08), y = c(.67,.25), 
                       label =  c("Reclassification rate drought-induced: 67.24%\nReclassification rate control: 76.67%",  "Reclassification rate drought-induced: 77.78%\nReclassification rate control: 77.78%"),
                       hjust = 0, size = 14)
COM.cap.irri.all.lineplot 
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("COM.cap.irri.all.lineplot.tiff",
       COM.cap.irri.all.lineplot, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 600)

##################################################################################################################################################################################################################
##################################################################################################################################################################################################################

### Rescaling analysis of Drought x Cropping System
### 1. AOA

### 1 A. AOA Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
aoa.bulk_dist_bc <- vegdist(t(aoa.asv.bulk1), method = "bray")
aoa.bulk_dist_bc
# metadata
aoa.meta.bulk
aoa.meta.bulk.ed <- aoa.meta.bulk[,-30:-45]
aoa.meta.bulk.ed$x <- as.factor(aoa.meta.bulk.ed$x)
aoa.meta.bulk.ed$Block <- as.factor(aoa.meta.bulk.ed$Block)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aoa.bulk_dist_bc ~ x, data = aoa.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 44 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.bulk.x <- CAPdiscrim(aoa.bulk_dist_bc ~ x, data = aoa.meta.bulk.ed, m = 44, permutations = 9999, add = TRUE) # 94.167% 
aoa.cap.bulk.x #Significance of this percentage was 0.00010001 
#cont.D (n=20) correct: 90 percent
#cont.K (n=20) correct: 100 percent
#cont.M (n=20) correct: 95 percent
#rain.D (n=20) correct: 95 percent
#rain.K (n=20) correct: 100 percent
#rain.M (n=20) correct: 85 percent

success <- cbind(data.frame(aoa.cap.bulk.x$group), data.frame(aoa.cap.bulk.x$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.bulk.irri$PCoA)
success <- success[order(success$source), ]
success
#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "AOA.cap.sucess.x.BS.csv")

aoa.cap1.bulk <- paste("CAP1 (", round((100/sum(aoa.cap.bulk.x$lda.other$svd^2) * aoa.cap.bulk.x$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
aoa.cap2.bulk <- paste("CAP2 (", round((100/sum(aoa.cap.bulk.x$lda.other$svd^2) * aoa.cap.bulk.x$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(aoa.cap.bulk.x$x), aes(x = aoa.cap.bulk.x$x[,1])) +
  geom_density(aes(fill = aoa.meta.bulk.ed$x), alpha = 0.4) + 
  theme_classic() +
  xlab(aoa.cap1.bulk) 

# point plot ==============================

ggplot(as.data.frame(aoa.cap.bulk.x$x), aes(x = aoa.cap.bulk.x$x[,1], y = aoa.cap.bulk.x$x[,2])) +
  geom_point(aes(shape = factor(aoa.meta.bulk.ed$x), color = aoa.meta.bulk.ed$Treatment), size = 2) +
  xlab(aoa.cap1.bulk) + ylab(aoa.cap2.bulk)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

AOA.cap.BS.points <- aoa.cap.bulk.x$x[,1]

AOA.cap.BS.points.df <- AOA.cap.BS.points %>% as.data.frame() %>%
  group_by(aoa.meta.bulk.ed$Irrigation, aoa.meta.bulk.ed$Treatment, aoa.meta.bulk.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'aoa.meta.bulk.ed$Irrigation', Treatment = 'aoa.meta.bulk.ed$Treatment', Date = 'aoa.meta.bulk.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

AOA.cap.BS.points.df <- mutate(AOA.cap.BS.points.df, x = paste(Irrigation, Treatment)) 
str(AOA.cap.BS.points.df)

AOA.cap.BS.points.df$Date <- factor(AOA.cap.BS.points.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
AOA.cap.BS.points.df$Date <- as.Date(AOA.cap.BS.points.df$Date)
AOA.cap.BS.points.df$Date
AOA.cap.BS.points.df$x <- as.factor(AOA.cap.BS.points.df$x)

### 1 B. AOA Rhizosphere

# run Bray-Curtis beta diversity on bulk soil
aoa.rh_dist_bc <- vegdist(t(aoa.asv.rh1), method = "bray")
aoa.rh_dist_bc 
# metadata
aoa.meta.rh <- aoa.meta.df[121:192,]
aoa.meta.rh.ed <- aoa.meta.rh[,c(-14:-29,-42:-45)]
aoa.meta.rh.ed$x <- as.factor(aoa.meta.rh.ed$x)
aoa.meta.rh.ed$Block <- as.factor(aoa.meta.rh.ed$Block)
aoa.meta.rh.ed$Date <- factor(aoa.meta.rh.ed$Date)
str(aoa.meta.rh.ed)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aoa.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aoa.rh_dist_bc ~ x, data = aoa.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 38  PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.rh.x <- CAPdiscrim(aoa.rh_dist_bc ~ x, data = aoa.meta.rh.ed, m = 38, permutations = 9999, add = TRUE) # 90.278% 
aoa.cap.rh.x #Overall classification success (m=38) : 90.2777777777778 percent
#Significance of this percentage was 0.00010001
#cont.D (n=12) correct: 100 percent
#cont.K (n=12) correct: 91.6666666666667 percent
#cont.M (n=12) correct: 100 percent
#rain.D (n=12) correct: 83.3333333333333 percent
#rain.K (n=12) correct: 75 percent
#rain.M (n=12) correct: 91.6666666666667 percent

success <- cbind(data.frame(aoa.cap.rh.x$group), data.frame(aoa.cap.rh.x$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.rh.x$PCoA)
success <- success[order(success$source), ]
success
#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "AOA.cap.sucess.x.RS.csv")

aoa.cap1.rh <- paste("CAP1 (", round((100/sum(aoa.cap.rh.x$lda.other$svd^2) * aoa.cap.rh.x$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
aoa.cap2.rh <- paste("CAP2 (", round((100/sum(aoa.cap.rh.x$lda.other$svd^2) * aoa.cap.rh.x$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(aoa.cap.rh.x$x), aes(x = aoa.cap.rh.x$x[,1])) +
  geom_density(aes(fill =  aoa.meta.rh.ed$x), alpha = 0.4) + 
  theme_classic() +
  xlab(aoa.cap1.rh) 

# point plot ==============================

ggplot(as.data.frame(aoa.cap.rh.x$x), aes(x = aoa.cap.rh.x$x[,1], y = aoa.cap.rh.x$x[,2])) +
  geom_point(aes(shape = factor(aoa.meta.rh.ed$x), color = aoa.meta.rh.ed$Treatment), size = 2) +
  xlab(aoa.cap1.rh) + ylab(aoa.cap2.rh)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

AOA.cap.rh.points <- aoa.cap.rh.x$x[,1]

AOA.cap.points.rh.df <- AOA.cap.rh.points %>% as.data.frame() %>%
  group_by(aoa.meta.rh.ed$Irrigation, aoa.meta.rh.ed$Treatment, aoa.meta.rh.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'aoa.meta.rh.ed$Irrigation', Treatment = 'aoa.meta.rh.ed$Treatment', Date = 'aoa.meta.rh.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

AOA.cap.points.rh.df <- mutate(AOA.cap.points.rh.df, x = paste(Irrigation, Treatment)) 
str(AOA.cap.points.rh.df)

AOA.cap.points.rh.df$Date <- factor(AOA.cap.points.rh.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05"))
AOA.cap.points.rh.df$Date <- as.Date(AOA.cap.points.rh.df$Date)
AOA.cap.points.rh.df$Date
AOA.cap.points.rh.df$x <- as.factor(AOA.cap.points.rh.df$x)


# line point of mean cap ==============================

# combine
AOA.cap.BS.points.df2 <- AOA.cap.BS.points.df %>%
  mutate(Type = "Bulk Soil")
AOA.cap.points.rh.df2 <- AOA.cap.points.rh.df %>%
  mutate(Type = "Rhizossphere")
AOA.cap.x.all <- merge(AOA.cap.BS.points.df2,AOA.cap.points.rh.df2,all = T)


# plot
AOA.cap.x.all.line <- ggplot(AOA.cap.x.all, aes(x = Date , y = mean, color = Treatment, group = x)) +
  geom_point(size = 1.5) +
  theme_bw() +
  #theme_classic() +
  geom_line(aes(linetype = Irrigation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Treatment), linetype=0, alpha=0.05) +
  xlab("") +
  ylab("Mean of CAP1") +
  facet_wrap(~ Type, strip.position="right", nrow = 2)+
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00"),
                    name = "Cropping system",
                    labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_linetype_manual(name = "Treatment",
                        labels = c("control", "drought"),
                        values = c("solid", "dashed")) +
  scale_x_date(date_breaks = "1 month",date_labels = "%B")+
  theme(axis.text.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 18, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(size = 20, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.background = element_blank()) +
  theme(legend.position = "bottom", legend.box="horizontal")+
  geom_vline(data=filter(AOA.cap.x.all, Type=="Bulk Soil"), aes(xintercept = as.Date("2022-07-14")), linetype="dashed",color = "black", linewidth = 0.5) 
  geom_text(data=filter(AOA.cap.x.all, Type=="Bulk Soil"), aes(x = as.Date("2022-07-10"), y = -5), hjust = 0, size = 7, label = "Rewetting", color = "black")
AOA.cap.x.all.line

library(cowplot)
AOA.cap.x.all.lineplot <- ggdraw(AOA.cap.x.all.line) + draw_plot_label(x = c(0.08,0.08), y = c(.67,.24), 
                       label =  c("Reclassification rate drought-induced: 76.67%\nReclassification rate control: 81.67%",  "Reclassification rate drought-induced: 83.33%\nReclassification rate control: 88.88%"),
                       hjust = 0, size = 14)
AOA.cap.x.all.lineplot 

setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOA.cap.irri.all.lineplot.tiff",
       AOA.cap.irri.all.lineplot, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 600)

### 2. AOB

### 2 A. AOB Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
aob.bulk_dist_bc <- vegdist(t(aob.asv.bulk1), method = "bray")
aob.bulk_dist_bc
# metadata
# metadata
aob.meta.bulk
aob.meta.bulk.ed <- aob.meta.bulk[,-30:-45]
aob.meta.bulk.ed$x <- as.factor(aob.meta.bulk.ed$x)
aob.meta.bulk.ed$Block <- as.factor(aob.meta.bulk.ed$Block)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:100) {
  cap <- CAPdiscrim(aob.bulk_dist_bc ~ Irrigation, data = aob.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 68 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.bulk.irri <- CAPdiscrim(aob.bulk_dist_bc ~ Irrigation, data = aob.meta.bulk.ed, m = 68, permutations = 9999, add = TRUE) # 69.747% 
aob.cap.bulk.irri  #Significance of this percentage was 0.00020002 
#Control (n=59) correct: 66.1016949152542 percent
#Rainout (n=60) correct: 73.3333333333333 percent

success <- cbind(data.frame(aoa.cap.bulk.irri$group), data.frame(aoa.cap.bulk.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.bulk.irri$PCoA)
success <- success[order(success$source), ]
success
#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "AOB.cap.sucess.irri.BS.csv")

aob.cap1.bulk <- paste("CAP1 (", round((100/sum(aob.cap.bulk.irri$lda.other$svd^2) * aob.cap.bulk.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
aob.cap2.bulk <- paste("CAP2 (", round((100/sum(aob.cap.bulk.irri$lda.other$svd^2) * aob.cap.bulk.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(aob.cap.bulk.irri$x), aes(x = aob.cap.bulk.irri$x[,1])) +
  geom_density(aes(fill = aob.meta.bulk.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(aob.cap1.bulk) 

# point plot ==============================

ggplot(as.data.frame(aob.cap.bulk.irri$x), aes(x = aob.cap.bulk.irri$x[,1], y = aob.cap.bulk.irri$x[,2])) +
  geom_point(aes(shape = factor(aob.meta.bulk.ed$Irrigation), color = aob.meta.bulk.ed$Treatment), size = 2) +
  xlab(aob.cap1.bulk) + ylab(aob.cap2.bulk)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

AOB.cap.BS.points <- aob.cap.bulk.irri$x[,1]

AOB.cap.BS.points.df <- AOB.cap.BS.points %>% as.data.frame() %>%
  group_by(aob.meta.bulk.ed$Irrigation, aob.meta.bulk.ed$Treatment, aob.meta.bulk.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'aob.meta.bulk.ed$Irrigation', Treatment = 'aob.meta.bulk.ed$Treatment', Date = 'aob.meta.bulk.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

AOB.cap.BS.points.df <- mutate(AOB.cap.BS.points.df, x = paste(Irrigation, Treatment)) 
str(AOB.cap.BS.points.df)

AOB.cap.BS.points.df$Date <- factor(AOB.cap.BS.points.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
AOB.cap.BS.points.df$Date <- as.Date(AOB.cap.BS.points.df$Date)
AOB.cap.BS.points.df$Date
AOB.cap.BS.points.df$x <- as.factor(AOB.cap.BS.points.df$x)

### 2 B. AOB Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
aob.rh_dist_bc <- vegdist(t(aob.asv.rh1), method = "bray")
aob.rh_dist_bc
# metadata
aob.meta.rh
aob.meta.rh.ed <- aob.meta.rh[,c(-14:-29,-42:-45)]
aob.meta.rh.ed$x <- as.factor(aob.meta.rh.ed$x)
str(aob.meta.rh.ed)
aob.meta.rh.ed$Date <- factor(aob.meta.rh.ed$Date)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(aob.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(aob.rh_dist_bc ~ Irrigation, data = aob.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 14 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aob.cap.rh.irri <- CAPdiscrim(aob.rh_dist_bc ~ Irrigation, data = aob.meta.rh.ed, m = 14 , permutations = 9999, add = TRUE) # 69.44 % 
aob.cap.rh.irri #Significance of this percentage was 0.00260026 
#Control (n=36) correct: 72.2222222222222 percent
#Rainout (n=36) correct: 66.6666666666667 percent

success <- cbind(data.frame(aob.cap.rh.irri$group), data.frame(aob.cap.rh.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aob.cap.rh.irri$PCoA)
success <- success[order(success$source), ]
success

#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "AOB.cap.sucess.irri.RS.csv")

aob.cap1.rh <- paste("CAP1 (", round((100/sum(aob.cap.rh.irri$lda.other$svd^2) * aob.cap.rh.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
aob.cap2.rh <- paste("CAP2 (", round((100/sum(aob.cap.rh.irri$lda.other$svd^2) * aob.cap.rh.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(aob.cap.rh.irri$x), aes(x = aob.cap.rh.irri$x[,1])) +
  geom_density(aes(fill = aob.meta.rh.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(aob.cap1.rh) 

# point plot ==============================

ggplot(as.data.frame(aob.cap.rh.irri$x), aes(x = aob.cap.rh.irri$x[,1], y = aob.cap.rh.irri$x[,2])) +
  geom_point(aes(shape = factor(aob.meta.rh.ed$Irrigation), color = aob.meta.rh.ed$Treatment), size = 2) +
  xlab(aob.cap1.rh) + ylab(aob.cap2.rh)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

AOB.cap.RS.points <- aob.cap.rh.irri$x[,1]

AOB.cap.RS.points.df <- AOB.cap.RS.points %>% as.data.frame() %>%
  group_by(aob.meta.rh.ed$Irrigation, aob.meta.rh.ed$Treatment, aob.meta.rh.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'aob.meta.rh.ed$Irrigation', Treatment = 'aob.meta.rh.ed$Treatment', Date = 'aob.meta.rh.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

AOB.cap.RS.points.df <- mutate(AOB.cap.RS.points.df, x = paste(Irrigation, Treatment)) 
str(AOB.cap.RS.points.df)

AOB.cap.RS.points.df$Date <- factor(AOB.cap.RS.points.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05"))
AOB.cap.RS.points.df$Date <- as.Date(AOB.cap.RS.points.df$Date)
AOB.cap.RS.points.df$Date
AOB.cap.RS.points.df$x <- as.factor(AOB.cap.RS.points.df$x)

# combine
AOB.cap.BS.points.df2 <- AOB.cap.BS.points.df %>%
  mutate(Type = "Bulk Soil")
AOB.cap.RS.points.df2 <- AOB.cap.RS.points.df %>%
  mutate(Type = "Rhizossphere")
AOB.cap.irri.all <- merge(AOB.cap.BS.points.df2,AOB.cap.RS.points.df2,all = T)

# plot
AOB.cap.irri.all.line <- ggplot(AOB.cap.irri.all, aes(x = Date , y = mean, color = Treatment, group = x)) +
  geom_point(size = 1.5) +
  theme_bw() +
  #theme_classic() +
  geom_line(aes(linetype = Irrigation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Treatment), linetype=0, alpha=0.05) +
  xlab("") +
  ylab("Mean of CAP1") +
  facet_wrap(~ Type, strip.position="right", nrow = 2)+
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00"),
                    name = "Cropping system",
                    labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_linetype_manual(name = "Treatment",
                        labels = c("control", "drought"),
                        values = c("solid", "dashed")) +
  scale_x_date(date_breaks = "1 month",date_labels = "%B")+
  theme(axis.text.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 18, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(size = 20, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.background = element_blank()) +
  theme(legend.position = "bottom", legend.box="horizontal")+
  geom_vline(data=filter(AOB.cap.irri.all, Type=="Bulk Soil"), aes(xintercept = as.Date("2022-07-14")), linetype="dashed",color = "black", linewidth = 0.5) +
  geom_text(data=filter(AOB.cap.irri.all, Type=="Bulk Soil"), aes(x = as.Date("2022-07-10"), y = -4), hjust = 0, size = 7, label = "Rewetting", color = "black")
AOB.cap.irri.all.line

library(cowplot)
AOB.cap.irri.all.lineplot <- ggdraw(AOB.cap.irri.all.line) + draw_plot_label(x = c(0.08,0.08), y = c(.67,.25), 
                       label =  c("Reclassification rate drought-induced: 73.33%\nReclassification rate control: 66.1%",  "Reclassification rate drought-induced: 66.67%\nReclassification rate control: 72.22%"),
                       hjust = 0, size = 14)
AOB.cap.irri.all.lineplot 
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOB.cap.irri.all.lineplot.tiff",
       AOB.cap.irri.all.lineplot, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 600)

### 3. COMAMMOX

### 3 A. COMA Bulk Soil

# run Bray-Curtis beta diversity on bulk soil
com.bulk_dist_bc <- vegdist(t(com.asv.bulk1), method = "bray")
com.bulk_dist_bc
# metadata
# metadata
com.meta.bulk
com.meta.bulk.ed <- com.meta.bulk[,-30:-45]
com.meta.bulk.ed$x <- as.factor(com.meta.bulk.ed$x)
com.meta.bulk.ed$Block <- as.factor(com.meta.bulk.ed$Block)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.bulk_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(com.bulk_dist_bc ~ Irrigation, data = com.meta.bulk.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 25 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.bulk.irri <- CAPdiscrim(com.bulk_dist_bc ~ Irrigation, data = com.meta.bulk.ed, m = 25, permutations = 9999, add = TRUE) # 72.03% 
com.cap.bulk.irri  #Significance of this percentage was 0.00010001 
#Control (n=60) correct: 76.6666666666667 percent
#Rainout (n=58) correct: 67.2413793103448 percent

success <- cbind(data.frame(com.cap.bulk.irri$group), data.frame(com.cap.bulk.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.bulk.irri$PCoA)
success <- success[order(success$source), ]
success
#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "COM.cap.sucess.irri.BS.csv")

com.cap1.bulk <- paste("CAP1 (", round((100/sum(com.cap.bulk.irri$lda.other$svd^2) * com.cap.bulk.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
com.cap2.bulk <- paste("CAP2 (", round((100/sum(com.cap.bulk.irri$lda.other$svd^2) * com.cap.bulk.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(com.cap.bulk.irri$x), aes(x = com.cap.bulk.irri$x[,1])) +
  geom_density(aes(fill = com.meta.bulk.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(com.cap1.bulk) 

# point plot ==============================

ggplot(as.data.frame(com.cap.bulk.irri$x), aes(x = com.cap.bulk.irri$x[,1], y = com.cap.bulk.irri$x[,2])) +
  geom_point(aes(shape = factor(com.meta.bulk.ed$Irrigation), color = com.meta.bulk.ed$Treatment), size = 2) +
  xlab(com.cap1.bulk) + ylab(com.cap2.bulk)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

COM.cap.BS.points <- com.cap.bulk.irri$x[,1]

COM.cap.BS.points.df <- COM.cap.BS.points %>% as.data.frame() %>%
  group_by(com.meta.bulk.ed$Irrigation, com.meta.bulk.ed$Treatment, com.meta.bulk.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'com.meta.bulk.ed$Irrigation', Treatment = 'com.meta.bulk.ed$Treatment', Date = 'com.meta.bulk.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

COM.cap.BS.points.df <- mutate(COM.cap.BS.points.df, x = paste(Irrigation, Treatment)) 
str(COM.cap.BS.points.df)

COM.cap.BS.points.df$Date <- factor(COM.cap.BS.points.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
COM.cap.BS.points.df$Date <- as.Date(COM.cap.BS.points.df$Date)
COM.cap.BS.points.df$Date
COM.cap.BS.points.df$x <- as.factor(COM.cap.BS.points.df$x)

### 3 B. COMAMMOX Rhizosphere

# run Bray-Curtis beta diversity on rhizosphere
com.rh_dist_bc <- vegdist(t(com.asv.rh1), method = "bray")
com.rh_dist_bc
# metadata
com.meta.rh
com.meta.rh.ed <- com.meta.rh[,c(-14:-29,-42:-45)]
com.meta.rh.ed$x <- as.factor(com.meta.rh.ed$x)
com.meta.rh.ed$Date <- factor(com.meta.rh.ed$Date)
# run CAP on increasing numbers of PCoA axes to check how many axes need to be included in the model (diagnostics).
nc <- nrow(as.matrix(com.rh_dist_bc))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
set.seed(133)
for (i in 1:50) {
  cap <- CAPdiscrim(com.rh_dist_bc ~ Irrigation, data = com.meta.rh.ed, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

#here for example you would choose 7 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
com.cap.rh.irri <- CAPdiscrim(com.rh_dist_bc ~ Irrigation, data = com.meta.rh.ed, m = 7 , permutations = 9999, add = TRUE) #  77.78% 
com.cap.rh.irri #Significance of this percentage was 0.00010001  
#Control (n=36) correct: 77.78 percent
#Rainout (n=36) correct: 77.78 percent

success <- cbind(data.frame(com.cap.rh.irri$group), data.frame(com.cap.rh.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(com.cap.rh.irri$PCoA)
success <- success[order(success$source), ]
success

#setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
#write.csv(success, file = "COM.cap.sucess.irri.RS.csv")

com.cap1.rh <- paste("CAP1 (", round((100/sum(com.cap.rh.irri$lda.other$svd^2) * com.cap.rh.irri$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
com.cap2.rh <- paste("CAP2 (", round((100/sum(com.cap.rh.irri$lda.other$svd^2) * com.cap.rh.irri$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(com.cap.rh.irri$x), aes(x = com.cap.rh.irri$x[,1])) +
  geom_density(aes(fill = com.meta.rh.ed$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(com.cap1.rh) 

# point plot ==============================

ggplot(as.data.frame(com.cap.rh.irri$x), aes(x = com.cap.rh.irri$x[,1], y = com.cap.rh.irri$x[,2])) +
  geom_point(aes(shape = factor(com.meta.rh.ed$Irrigation), color = com.meta.rh.ed$Treatment), size = 2) +
  xlab(com.cap1.rh) + ylab(com.cap2.rh)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

COM.cap.RS.points <- com.cap.rh.irri$x[,1]

COM.cap.RS.points.df <- COM.cap.RS.points %>% as.data.frame() %>%
  group_by(com.meta.rh.ed$Irrigation, com.meta.rh.ed$Treatment, com.meta.rh.ed$Date) %>%
  dplyr::rename(c(value = ., Irrigation = 'com.meta.rh.ed$Irrigation', Treatment = 'com.meta.rh.ed$Treatment', Date = 'com.meta.rh.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

COM.cap.RS.points.df <- mutate(COM.cap.RS.points.df, x = paste(Irrigation, Treatment)) 
str(COM.cap.RS.points.df)

COM.cap.RS.points.df$Date <- factor(COM.cap.RS.points.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05"))
COM.cap.RS.points.df$Date <- as.Date(COM.cap.RS.points.df$Date)
COM.cap.RS.points.df$Date
COM.cap.RS.points.df$x <- as.factor(COM.cap.RS.points.df$x)

# combine
COM.cap.BS.points.df2 <- COM.cap.BS.points.df %>%
  mutate(Type = "Bulk Soil")
COM.cap.RS.points.df2 <- COM.cap.RS.points.df %>%
  mutate(Type = "Rhizossphere")
COM.cap.irri.all <- merge(COM.cap.BS.points.df2,COM.cap.RS.points.df2,all = T)

# plot
COM.cap.irri.all.line <- ggplot(COM.cap.irri.all, aes(x = Date , y = mean, color = Treatment, group = x)) +
  geom_point(size = 1.5) +
  theme_bw() +
  #theme_classic() +
  geom_line(aes(linetype = Irrigation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Treatment), linetype=0, alpha=0.05) +
  xlab("") +
  ylab("Mean of CAP1") +
  facet_wrap(~ Type, strip.position="right", nrow = 2)+
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00"),
                    name = "Cropping system",
                    labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_linetype_manual(name = "Treatment",
                        labels = c("control", "drought"),
                        values = c("solid", "dashed")) +
  scale_x_date(date_breaks = "1 month",date_labels = "%B")+
  theme(axis.text.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 18, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(size = 20, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.background = element_blank()) +
  theme(legend.position = "bottom", legend.box="horizontal")+
  geom_vline(data=filter(COM.cap.irri.all, Type=="Bulk Soil"), aes(xintercept = as.Date("2022-07-14")), linetype="dashed",color = "black", linewidth = 0.5) +
  geom_text(data=filter(COM.cap.irri.all, Type=="Bulk Soil"), aes(x = as.Date("2022-07-10"), y = -4), hjust = 0, size = 7, label = "Rewetting", color = "black")
COM.cap.irri.all.line

library(cowplot)
COM.cap.irri.all.lineplot <- ggdraw(COM.cap.irri.all.line) + draw_plot_label(x = c(0.08,0.08), y = c(.67,.25), 
                       label =  c("Reclassification rate drought-induced: 67.24%\nReclassification rate control: 76.67%",  "Reclassification rate drought-induced: 77.78%\nReclassification rate control: 77.78%"),
                       hjust = 0, size = 14)
COM.cap.irri.all.lineplot 
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("COM.cap.irri.all.lineplot.tiff",
       COM.cap.irri.all.lineplot, device = "tiff",
       width = 10, height =8, 
       units= "in", dpi = 600)


