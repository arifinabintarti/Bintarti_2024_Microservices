nc <- nrow(as.matrix(ISS.bray.B.BS))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:29) {
  cap <- CAPdiscrim(ISS.bray.B.BS ~ Irrigation, data = design.B.BS, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

ISS.cap.B <- CAPdiscrim(ISS.bray.B.BS ~ Irrigation, data = design.B.BS, m = 10, permutations = 9999, add = TRUE) #try with different

success <- cbind(data.frame(ISS.cap.B$group), data.frame(ISS.cap.B$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(ISS.cap.B$PCoA)
success <- success[order(success$source), ]
success

ISS.cap1.B <- paste("CAP1 (", round((100/sum(ISS.cap.B$lda.other$svd^2) * ISS.cap.B$lda.other$svd^2)[1],
                                    digits = 1), "%)", sep = "")
ISS.cap2.B <- paste("CAP2 (", round((100/sum(ISS.cap.B$lda.other$svd^2) * ISS.cap.B$lda.other$svd^2)[2],
                                    digits = 1), "%)", sep = "")

# density plot ==============================

ggplot(as.data.frame(ISS.cap.B$x), aes(x = ISS.cap.B$x[,1])) +
  geom_density(aes(fill = design.B.BS$Irrigation), alpha = 0.4) + 
  theme_classic() +
  xlab(ISS.cap1.B) 

# point plot ==============================

ggplot(as.data.frame(ISS.cap.B$x), aes(x = ISS.cap.B$x[,1], y = ISS.cap.B$x[,2])) +
  geom_point(aes(shape = factor(design.B.BS$Irrigation), color = design.B.BS$Treatment), size = 2) +
  xlab(ISS.cap1.B) + ylab(ISS.cap2.B)  +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) + theme_classic()

# mean cap ==============================

cap.B.points <- ISS.cap.B$x[,1]

cap.points.B.df <- cap.B.points %>% as.data.frame() %>%
  group_by(design.B.BS$Irrigation, design.B.BS$Treatment, design.B.BS$Date) %>%
  rename(c(value = ., Irrigation = `design.B.BS$Irrigation`, Treatment =  `design.B.BS$Treatment`, Date =  `design.B.BS$Date`)) %>%
  summarise(mean = mean(value),
            sd = sd(value))

cap.points.B.df <- mutate(cap.points.B.df, x = paste(Irrigation, Treatment)) 

# point point of mean cap ==============================

cap.points.B <- ggplot(cap.points.B.df, aes(x = as.Date(Date) , y = mean, color = Treatment, shape = Irrigation)) +
  theme_classic() +
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd), 
                  position=position_jitter(width=1.5)) +
  xlab("") +
  ylab("Mean of CAP1") +
  scale_color_manual(values = c("#009E73", "#FF618C",  "#E69F00"),
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  theme(legend.position = "bottom")  +
  annotate(geom = "text", x = as.Date("2022-05-01"), y = -3.4, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-06-04"), y = -2.3, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-07-04"), y = -4.2, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-07-12"), y = -2.3, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-09-06"), y = -2.7, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-05-02"), y = 3, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-06-04"), y = 2.3, hjust = 0, size = 3, label = "100%", color = "black") +
  
  annotate(geom = "text", x = as.Date("2022-06-27"), y = 3.6, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-07-07"), y = 3.8, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-07-08"), y = 0.8, hjust = 0, size = 3, label = "75%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-07-17"), y = 1.3, hjust = 0, size = 3, label = "100%", color = "black") +
  
  annotate(geom = "text", x = as.Date("2022-09-04"), y = 2.6, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-09-17"), y = 2, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-09-06"), y = 3.5, hjust = 0, size = 3, label = "100%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-04-28"), y =-5.5, hjust = 0, size = 3, label = "Reclassification rate drought-induced: 97%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-04-28"), y =-6.5, hjust = 0, size = 3, label = "Reclassification rate control: 100%", color = "black")
cap.points.B

# line point of mean cap ==============================

cap.line.B <- ggplot(cap.points.B.df, aes(x = as.Date(Date) , y = mean, color = Treatment, group = x)) +
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
  annotate(geom = "text", x = as.Date("2022-04-25"), y = -3.5, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-05-31"), y = -2.1, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-07-04"), y = -3.8, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-07-17"), y = -2, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-09-10"), y = -2.7, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-04-25"), y = 3.4, hjust = 0, size = 3, label = "100%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-05-31"), y = 2, hjust = 0, size = 3, label = "100%", color = "black") +
  
  annotate(geom = "text", x = as.Date("2022-06-28"), y = 3.6, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-07-04"), y = 4.2, hjust = 0, size = 3, label = "100%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-07-02"), y = 0.8, hjust = 0, size = 3, label = "75%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-07-17"), y = 1.8, hjust = 0, size = 3, label = "100%", color = "black") +
  
  annotate(geom = "text", x = as.Date("2022-09-06"), y = 1, hjust = 0, size = 3, label = "100%", color = "#FF618C") +
  annotate(geom = "text", x = as.Date("2022-09-15"), y = 1.4, hjust = 0, size = 3, label = "75%", color = "#009E73") +
  annotate(geom = "text", x = as.Date("2022-09-10"), y = 3.3, hjust = 0, size = 3, label = "100%", color = "#E69F00") +
  
  annotate(geom = "text", x = as.Date("2022-04-28"), y =-5, hjust = 0, size = 3.5, label = "Reclassification rate drought-induced: 97%", color = "black") +
  annotate(geom = "text", x = as.Date("2022-04-28"), y =-6, hjust = 0, size = 3.5, label = "Reclassification rate control: 100%", color = "black") +
  
  geom_vline(xintercept = as.Date("2022-07-14"), linetype="dashed", color = "black", linewidth = 0.5) +
  annotate(geom = "text", x = as.Date("2022-07-01"), y = -5, hjust = 0, size = 3.5, label = "Rewetting", color = "black")

cap.line.B
