# Plotting

qPCR.RS$x
qPCR.RS.ed <- qPCR.RS %>%
  mutate(x = factor(x,levels = c("cont.D","rain.D","cont.K","rain.K","cont.M","rain.M")))
label <- c(`D` ="BIODYN", 
           `K` ="CONFYM", 
           `M` ="CONMIN")
qPCR.RS.ed$sampling.date <- factor(qPCR.RS.ed$sampling.date, levels = c("Apr 28th", "Jun 1st", "Jul 5th"),
                          labels = c("Apr", "Jun", "Jul"))

# 1. Rhizosphere - AOA

aoa.cop.rh.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=AOA_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  #ylim(0,7e+08)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('AOA'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('AOA')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("B. Rhizosphere")
aoa.cop.rh.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/')
ggsave("AOA_copies.Rhizos.tiff",
       aoa.cop.rh.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)


# 2. Rhizosphere - AOB

aob.cop.rh.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=AOB_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  #ylim(0,7e+08)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('AOB'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('AOB')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        #strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("B. Rhizosphere")
aob.cop.rh.plot


# 3. Rhizosphere - Comammox A

comA.cop.rh.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=ComA_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  #ylim(0,7e+08)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox A'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('Comammox A')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        #strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("B. Rhizosphere")
comA.cop.rh.plot

# 4. Rhizosphere- Comammox B

comB.cop.rh.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=ComB_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  #ylim(0,7e+08)+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('Comammox B')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  scale_alpha_manual(values = c(1, 0.5),
                     labels=c('Control', 'Drought'),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        #strip.text = element_text(size=15),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("B. Rhizosphere")
comB.cop.rh.plot

# 4. Rhizosphere - 16S

tot.cop.rh.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=Tot_nbc_per_ngDNA)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  ylim(0,6e+05)+
  ylab(bquote('16S'~(copies~ng^-1~DNA)))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab('16S')+
  #scale_fill_manual(values = c("#009E73","#FF618C","#E69F00"))+
  #scale_alpha_manual(values = c(1, 0.5),
  #labels=c('Control', 'Drought'),
  #guide = guide_legend(override.aes = list(fill = "black"))) +
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("B. Rhizosphere")
tot.cop.rh.plot

# adding xy position for the pairwise comparisons among treatments (emmeans results)
tot.ng.dws.rh.emm.rstat <- qPCR.RS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(Tot_nbc_per_ngDNA  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.tot.rh)
tot.ng.dws.rh.emm.rstat 
# add x y position
tot.ng.dws.rh.xy <- tot.ng.dws.rh.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
tot.cop.rh.plot2 <- tot.cop.rh.plot + 
  stat_pvalue_manual(tot.ng.dws.rh.xy,
                     #step.increase = 1,
                     label = "p = {scales::pvalue(p.adj)}",size=3.5, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1.2, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
tot.cop.rh.plot3 <- tot.cop.rh.plot2 +ylim(0,3.2e+10)

setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("16Scopies_RS.tiff",
       tot.cop.rh.plot, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)

copies16S.All <-  tot.cop.pwc.plot3 / tot.cop.rh.plot 
copies16S.All
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("16S.All.tiff",
       copies16S.All, device = "tiff",
       width = 9, height =7, 
       units= "in", dpi = 600)

# 5. Rhizosphere - AOA/16S RATIO

AOA_16S.rat.rh.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=AOA_16S_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('AOA/16S (%)')+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("B. Rhizosphere")
AOA_16S.rat.rh.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
AOA_16S.rat.rh.emm.rstat <- qPCR.RS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(AOA_16S_ratio_percent  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.aoa_16S.rh)
AOA_16S.rat.rh.emm.rstat 
# add x y position
AOA_16S.rat.rh.emm.xy <- AOA_16S.rat.rh.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
AOA_16S.rat.rh.emm.pwc.plot2 <- AOA_16S.rat.rh.plot + 
  stat_pvalue_manual(AOA_16S.rat.rh.emm.xy,
                     #step.increase = 1,
                     #label = "p.adj.signif",size=3.5,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
AOA_16S.rat.rh.plot3 <- AOA_16S.rat.rh.emm.pwc.plot2 +   ylim(0,3.5)
AOA_16S.rat.rh.plot3 
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("AOA_16S_ratio_RS.tiff",
       AOA_16S.rat.rh.plot3, device = "tiff",
       width = 11, height =6, 
       units= "in", dpi = 600)

AOA_16S.All <-  AOA_16S.rat.plot3  / AOA_16S.rat.rh.plot3 
AOA_16S.All
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("AOA_16S.All.tiff",
       AOA_16S.All, device = "tiff",
       width = 9, height =7, 
       units= "in", dpi = 600)



#6. Rhizosphere - AOB/16S RATIO

AOB_16S.stat_text.RS <- data.frame(sampling.date = 0.5, AOB_16S_ratio_percent = 0.4 , fertilization="D", label="D x C x T *")

AOB_16S.rat.RS.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=AOB_16S_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_y_continuous(limits = c(0, 5))+
  #ylim(0,5)+
  labs(title = "Rhizosphere", subtitle = "B")+
  ylab('AOB/16S (%)')+
  #ylab(bquote('Comammmox B abundance'~(copies~g^-1~dry~soil)))+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control', 'Biodyn-drought', 'Confym-control', 
                             'Confym-drought', 'Conmin-control', 'Conmin-drought'))+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        #strip.text = element_blank(),
        axis.text.y = element_text(size = 18),
        #axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_markdown(size=19),
        plot.title = element_text(size=25, face="bold"),
        plot.subtitle = element_text(size=20, face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 geom_label(data = AOB_16S.stat_text.RS,label=AOB_16S.stat_text.RS$label,hjust=0, colour="black", size=4, fontface="bold")
AOB_16S.rat.RS.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
aob_16S_percent_rat_rh_emm <- qPCR.RS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(AOB_16S_ratio_percent ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model =  t.aob.RS.arc.aov2)
aob_16S_percent_rat_rh_emm 
aob_16S_percent_rh.xy <- aob_16S_percent_rat_rh_emm %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
AOB_16S.rat.RS.plot2 <- AOB_16S.rat.RS.plot + 
  stat_pvalue_manual(aob_16S_percent_rh.xy,x = "sampling.date", y.position = 4.8,
                     label = "p.adj.signif",size=5,
                     tip.length = 0.01, hide.ns = F)
AOB_16S.rat.RS.plot2

setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("AOB_16S.All.tiff",
       AOB_16S.All, device = "tiff",
       width = 9, height =7, 
       units= "in", dpi = 600)



# 7. Rhizosphere - ComammoxA/16S RATIO

#comA_16S.stat_text.RS <- data.frame(sampling.date = 0.5, ComA_16S_ratio_percent = 0.16 , fertilization="D", label="")

comA_16S.rat.RS.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=ComA_16S_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_y_continuous(limits = c(0, 1.2))+
  #ylim(0,5)+
  labs(subtitle = "D")+
  ylab('Comammox A/16S (%)')+
  #ylab(bquote('Comammmox B abundance'~(copies~g^-1~dry~soil)))+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control', 'Biodyn-drought', 'Confym-control', 
                             'Confym-drought', 'Conmin-control', 'Conmin-drought'))+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        #strip.text = element_text(size=18),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 18),
        #axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_markdown(size=19),
        #plot.title = element_text(size=25, face="bold"),
        plot.subtitle = element_text(size=20, face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
comA_16S.rat.RS.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("ComA_16S.All.tiff",
       ComA_16S.All, device = "tiff",
       width = 9, height =7, 
       units= "in", dpi = 600)



# 7. Rhizosphere - ComammoxB/16S RATIO

comB_16S.stat_text.RS <- data.frame(sampling.date = 0.5, ComB_16S_ratio_percent = 0.015 , fertilization="D", label="C *")

comB_16S.rat.RS.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=ComB_16S_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_classic() +
  scale_y_continuous(limits = c(0, 0.22))+
  #ylim(0,5)+
  labs(subtitle = "F")+
  ylab('Comammox B/16S (%)')+
  #ylab(bquote('Comammmox B abundance'~(copies~g^-1~dry~soil)))+
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control', 'Biodyn-drought', 'Confym-control', 
                             'Confym-drought', 'Conmin-control', 'Conmin-drought'))+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        #strip.text = element_text(size=18),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_markdown(size=19),
        #plot.title = element_text(size=25, face="bold"),
        plot.subtitle = element_text(size=20, face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 geom_label(data = comB_16S.stat_text.RS,label=comB_16S.stat_text.RS$label,hjust=0, colour="black", size=4, fontface="bold")
comB_16S.rat.RS.plot

setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("ComB_16S.All.tiff",
       ComB_16S.All, device = "tiff",
       width = 9, height =7, 
       units= "in", dpi = 600)




# 9. Rhizosphere - AOA/AOB RATIO
AOB_AOB.stat_text.RS <- data.frame(sampling.date = 0.5, AOA_AOB_ratio_percent = 40 , fertilization="D", label="C *\nD x T *")
AOA_AOB.rat.rh.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=AOA_AOB_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('Biodyn-control', 'Biodyn-drought', 'Confym-control', 
                             'Confym-drought', 'Conmin-control', 'Conmin-drought'))+
  labs(title="B. Rhizosphere")+
  ylab('AOA/AOB (%)')+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.position = "none",
        legend.title = element_text(size=15, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_text(size=18),
        #strip.text = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 16,angle = 45, hjust = 1),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_markdown(size=19),
        plot.title = element_text(size=25, face="bold"),
        plot.subtitle = element_text(size=20, face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
 geom_label(data = AOB_AOB.stat_text.RS,label=AOB_AOB.stat_text.RS$label,hjust=0, colour="black", size=4, fontface="bold")
AOA_AOB.rat.rh.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
AOA_AOB_rh_emm <- qPCR.RS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(AOA_AOB_ratio_percent ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model =  AOA_AOB.RS.arc.aov2)
AOA_AOB_rh_emm
AOA_AOB_rh.xy <- AOA_AOB_rh_emm %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
AOA_AOB.rat.rh.plot2 <- AOA_AOB.rat.rh.plot + 
  stat_pvalue_manual(AOA_AOB_rh.xy,x = "sampling.date", y.position = 165,
                     label = "p.adj.signif",size=5,
                     tip.length = 0.01, hide.ns = F)
AOA_AOB.rat.rh.plot2 

AOA_AOB.All <- (AOA_AOB.rat.plot2  / AOA_AOB.rat.rh.plot2)+
 plot_layout(guides = "collect") & 
 theme(legend.position = 'bottom',legend.title = element_blank())
AOA_AOB.All
setwd('/Users/arifinabintarti/Documents/France/Figures/')
ggsave("AOA_AOB.All.tiff",
       AOA_AOB.All, device = "tiff",
       width = 10, height =12, 
       units= "in", dpi = 300, compression="lzw")

# 9. Rhizosphere - comA/comB RATIO

ComA_ComB.rat.rh.plot <- ggplot(qPCR.RS.ed, aes(x=sampling.date, y=ComA_ComB_ratio_percent)) +
  geom_boxplot(aes(group = var3, fill = x))+
  theme_bw() +
  scale_fill_manual(values = c("#009E73","#DAF1EB","#FF618C","#FFE8EE","#E69F00","#FBF1DA"),
                    labels=c('control (D)', 'drought (D)', 'control (K)', 
                             'drought (K)', 'control (M)', 'drought (M)'))+
  #labs(fill='Farming system', alpha= 'Drought')+
  #ylab(bquote('Comammox B'~italic(amoA)~'gene'~(copies~g^-1~dry~soil)))+
  ylab('Coma A/Coma B (%)')+
  facet_wrap(~ fertilization,scales="free_x", labeller = as_labeller(label))+
  theme(legend.title = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(size = 20, face='bold'),
        legend.text = element_text(size=15),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x =element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none", alpha="none")+ggtitle("B. Rhizosphere")
ComA_ComB.rat.rh.plot
# adding xy position for the pairwise comparisons among treatments (emmeans results)
ComA_ComB.rat.rh.emm.rstat <- qPCR.RS.ed %>%
  group_by(sampling.date, fertilization) %>%
  emmeans_test(ComA_ComB_ratio_percent  ~ irrigation, 
               p.adjust.method = "BH", 
               conf.level = 0.95, model = t.ComA_ComB_percent.rh)
ComA_ComB.rat.rh.emm.rstat
# add x y position
ComA_ComB.rat.rh.emm.xy <- ComA_ComB.rat.rh.emm.rstat  %>% 
  add_xy_position(x = "sampling.date", dodge = 0.8) # bulk soil
# plotting the pairwise comparisons among treatments (emmeans results)
ComA_ComB.rat.rh.plot2 <- ComA_ComB.rat.rh.plot + 
  stat_pvalue_manual(ComA_ComB.rat.rh.emm.xy ,
                     #step.increase = 1,
                     #label = "p.adj.signif",size=3.5,
                     label = "p = {scales::pvalue(p.adj)}",size=3, 
                     bracket.size = 0.6,#bracket.nudge.y = -0.05,
                     bracket.shorten = 1, color = "black",
                     tip.length = 0.005, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
ComA_ComB.rat.rh.plot3 <- ComA_ComB.rat.rh.plot2 +   ylim(0,1650)
ComA_ComB.rat.rh.plot3 

ComA_ComB.All <-  comA_comB.rat.plot3  / ComA_ComB.rat.rh.plot3 
ComA_ComB.All
setwd('D:/Fina/INRAE_Project/microservices_fig/qPCR')
ggsave("ComA_ComB.All.tiff",
       ComA_ComB.All, device = "tiff",
       width = 9, height =7, 
       units= "in", dpi = 600)
