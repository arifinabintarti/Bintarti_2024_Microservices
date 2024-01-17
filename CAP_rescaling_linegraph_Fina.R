### 1. AOA

### 1 A. Bulk Soil

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

#here for example you would choose 23 PCoA axes

# run the final CAP by including PCoA axes showing the highest reclassification rate
set.seed(13)
aoa.cap.bulk.irri <- CAPdiscrim(aoa.bulk_dist_bc ~ Irrigation, data = aoa.meta.bulk.ed, m = 23, permutations = 9999, add = TRUE) # 79.17% 
aoa.cap.bulk.irri 

success <- cbind(data.frame(aoa.cap.bulk.irri$group), data.frame(aoa.cap.bulk.irri$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(aoa.cap.bulk.irri$PCoA)
success <- success[order(success$source), ]
success
setwd('/Users/arifinabintarti/Documents/France/microservices/CAP_success/')
write.csv(success, file = "AOA.cap.sucess.irri.BS.csv")

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
  rename(c(value = ., Irrigation = 'aoa.meta.bulk.ed$Irrigation', Treatment = 'aoa.meta.bulk.ed$Treatment', Date =  'aoa.meta.bulk.ed$Date')) %>%
  summarise(mean = mean(value),
            sd = sd(value))

cap.points.B.df <- mutate(cap.points.B.df, x = paste(Irrigation, Treatment)) 
str(cap.points.B.df)

cap.points.B.df$Date <- factor(cap.points.B.df$Date, levels = c("Apr 28th", "Jun 1st", "Jul 5th", "Jul 20th", "Sept 13th"),
                          labels = c("2022-04-28", "2022-06-01", "2022-07-05", "2022-07-20", "2022-09-13"))
cap.points.B.df$Date <- as.Date(cap.points.B.df$Date)
cap.points.B.df$Date

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

cap.line.B <- ggplot(cap.points.B.df, aes(x = Date , y = mean, color = Treatment, group = x)) +
  geom_point(size = 1.5) + theme_classic() +
  geom_line(aes(linetype = Irrigation)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Treatment), linetype=0, alpha=0.05) +
  xlab("") +
  ylab("Mean of CAP1") +
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
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.background = element_blank()) +
  theme(legend.position = "bottom", legend.box="horizontal") +
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
  annotate(geom = "text", x = as.Date("2022-04-25"), y = -1.3, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-04-25"), y = -2, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-04-25"), y = -0.98, hjust = 0, size = 3, label = "75%", color = "#E69F00") +
  
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
 
  annotate(geom = "text", x = as.Date("2022-04-25"), y = 0.9, hjust = 0, size = 3, label = "75%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-04-25"), y = 1.45, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-04-25"), y = 0.3, hjust = 0, size = 3, label = "50%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-05-30"), y = 2, hjust = 0, size = 3, label = "75%", color = "black") +
  
  annotate(geom = "text", x = as.Date("2022-07-03"), y = 0.4, hjust = 0, size = 3, label = "50%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = 2, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-07-03"), y = 1.2, hjust = 0, size = 3, label = "75%", color = "#E69F00") +
 
  annotate(geom = "text", x = as.Date("2022-07-18"), y = 1.8, hjust = 0, size = 3, label = "100%", color = "black") +
 
  annotate(geom = "text", x = as.Date("2022-09-13"), y = 0.6, hjust = 0, size = 3, label = "75%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-09-13"), y = 1.9, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-09-13"), y = 1.3, hjust = 0, size = 3, label = "50%", color = "#E69F00") +

  annotate(geom = "text", x = as.Date("2022-04-28"), y =-5, hjust = 0, size = 3.5, label = "Reclassification rate drought-induced: 76.67%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-04-28"), y =-5.4, hjust = 0, size = 3.5, label = "Reclassification rate control: 81.67%", color = "black") +
  geom_vline(xintercept = as.Date("2022-07-14"), linetype="dashed", color = "black", linewidth = 0.5) +
  annotate(geom = "text", x = as.Date("2022-07-10"), y = -5, hjust = 0, size = 3.5, label = "Rewetting", color = "black")

cap.line.B
