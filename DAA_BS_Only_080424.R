## Heatmap DAA for Bulk Soil only ##

################################################################################
###compile 3 genes in one heatmap Bulk Soil Only###
################################################################################

#### Bulk Soil #####

#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/log2fold/')
#rr.comp <- read.csv("3genes.bulk.RR.csv", row.names = 1)
rr.comp.BS <- read.csv("RR_3Genes_BS_NoRhizo_080424.csv", row.names = 1)
head(rr.comp.BS)
#names(rr.comp)=str_sub(names(rr.comp),4)
view(rr.comp.BS)
#Set annotation
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
#ann.comp <- read.csv("3genes.anno.csv", row.names = 1)
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB')
ann.comp.BS <- read.csv("Anno_3Genes_BS_NoRhizo_080424.csv", row.names = 1)
view(ann.comp.BS)
#order rownames
rr.comp.ord.BS <- rr.comp.BS[rownames(ann.comp.BS), ]
view(rr.comp.ord.BS)
#remove the character before "_"
#rownames(rr.comp.ord) <- sub('.*_', '', rownames(rr.comp.ord))


#relative abund of bulk soil


# RA AOB
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/')
ann.AOB.BS <- read.csv("AOB.anno_prev80_BSOnly_080424.csv",row.names = 1)
# calculate relative abundance 
aob.asv.ra <- transform_sample_counts(aob.rare.1282.seq, function(x) x/sum(x))
aob.asv.ra
#aob.asv.ra.melt <- psmelt(aob.asv.ra)
aob.asv.ra.melt <- psmelt(aob.asv.ra) %>%
  group_by(OTU) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
View(aob.asv.ra.melt)
aob.asv.ra.melt$ra.perc <- aob.asv.ra.melt$Mean*100
aob.asv.ra.melt <- column_to_rownames(aob.asv.ra.melt, var = "OTU")
rownames(aob.asv.ra.melt)
# make asv and ra
aob.ra.BS <- ann.AOB.BS
rownames(aob.ra.BS)
aob.ra.BS$RA <- aob.asv.ra.melt$ra.perc[match(row.names(aob.ra.BS), row.names(aob.asv.ra.melt))]
view(aob.ra.BS)
rownames(aob.ra.BS) <- sub("ASV_", "AOB_ASV ", rownames(aob.ra.BS))

# RA AOA
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB/')
ann.AOA.BS <- read.csv("AOA_anno_prev80_BSOnly_080424.csv", row.names = 1)
# calculate relative abundance 
aoa.asv.ra <- transform_sample_counts(aoa.rare.min.physeq, function(x) x/sum(x))
aoa.asv.ra
aoa.asv.ra.melt <- psmelt(aoa.asv.ra) %>%
  group_by(OTU) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
aoa.asv.ra.melt$ra.perc <- aoa.asv.ra.melt$Mean*100
aoa.asv.ra.melt <- column_to_rownames(aoa.asv.ra.melt, var = "OTU")
# make asv and ra
aoa.ra.BS <- ann.AOA.BS
rownames(aoa.ra.BS)
aoa.ra.BS$RA <- aoa.asv.ra.melt$ra.perc[match(row.names(aoa.ra.BS), row.names(aoa.asv.ra.melt))]
view(aoa.ra.BS)
rownames(aoa.ra.BS) <- sub("ASV_", "AOA_ASV ", rownames(aoa.ra.BS))

# # RA COMAMMOX
setwd('/Users/arifinabintarti/Documents/France/microservices/DAA/glmmTMB')
ann.COM.BS <- read.csv("COM_anno_prev80_BSOnly_080424.csv", row.names = 1)
# calculate relative abundance 
com.asv.ra <- transform_sample_counts(com.rare.min.physeq, function(x) x/sum(x))
com.asv.ra
com.asv.ra.melt <- psmelt(com.asv.ra) %>%
  group_by(OTU) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
com.asv.ra.melt$ra.perc <- com.asv.ra.melt$Mean*100
view(com.asv.ra.melt)
com.asv.ra.melt <- column_to_rownames(com.asv.ra.melt, var = "OTU")
# make asv and ra
com.ra.BS <- ann.COM.BS
rownames(com.ra.BS)
com.ra.BS$RA <- com.asv.ra.melt$ra.perc[match(row.names(com.ra.BS), row.names(com.asv.ra.melt))]
view(com.ra.BS)
rownames(com.ra.BS) <- sub("ASV_", "COM_ASV ", rownames(com.ra.BS))

# Join everything
RA.comp.BS <- rbind(aob.ra.BS, aoa.ra.BS, com.ra.BS)
RA.comp.BS <- rownames_to_column(RA.comp.BS, var="ASV")

RA.comp.BS.df <- as.data.frame(RA.comp.BS)
view(RA.comp.BS.df)
RA.comp.BS.df2 <- RA.comp.BS.df[,-2]
view(RA.comp.BS.df2)
RA.comp.BS.df2 <- column_to_rownames(RA.comp.BS.df2, var="ASV")
str(RA.comp.BS.df2)

library(ComplexHeatmap)
# SET LEGEND
lgd1.BS <- Legend(labels = c("Nitrosolobus multiformis",
                              "Nitrosomonas communis",
                              "Nitrosospira sp."),
                              #"Nitrosospira sp",
                              #"Nitrosospira sp"),
               legend_gp = gpar(fill=c("#33A02C","#FF7F00","#6A3D9A")), #"#6A3D9A","#6A3D9A")),
               title= "AOB", labels_gp = gpar(fontsize=18),title_gp = gpar(fontsize = 20,fontface='bold'))

lgd2.BS <- Legend(labels = c("NS-Beta",
                          "NS-Delta",
                          #"NS-Delta",
                          "NS-Gamma",
                          #"NS-Gamma",
                          #"NS-Gamma",
                          "NT-Alpha"),
               legend_gp = gpar(fill=c("#DD701E","#E1823A","#E9A672","#7570B3")),
               title= "AOA",labels_gp = gpar(fontsize=18),title_gp = gpar(fontsize = 20,fontface='bold'))

lgd3.BS <- Legend(labels = c("Clade A-Nitrospira sp.",
                                 "Clade B-Nitrospira sp."),
                                 #"Clade B-Nitrospira-sp.LM-bin98",
                                 #"Clade B-Nitrospira-sp.LPPL-bin249",
                                 #"Clade B-Nitrospira-sp.Smid-bin44"),
               legend_gp = gpar(fill=c("#66C2A5","#FC8D62")),#"#FC8D62","#FC8D62","#FC8D62")),
               title= "COMAMMOX",labels_gp = gpar(fontsize=18),title_gp = gpar(fontsize = 20, fontface='bold'))
pd.BS = packLegend(lgd,lgd1.BS, lgd2.BS, lgd3.BS, direction = "vertical")
draw(pd.BS)
              
col.comp.ord.BS <- list("Taxonomy"=c("Nitrosolobus-multiformis-Nl1_2667636517"="#33A02C",
                              "Nitrosomonas-communis-Nm44_2676397764"="#FF7F00",
                              #"Nitrosomonas-europaea-ATCC-19718_637427314"="#FF7F00",
                              "Nitrosospira-sp-17Nsp14_2671457573"="#6A3D9A",
                              "Nitrosospira-sp_2630434854"="#6A3D9A",
                              "Nitrosospira-sp_2636913388"="#6A3D9A",
                                 "NS-Beta-1.3_OTU1_1-EU885554"="#DD701E",
                                 "NS-Delta-1.Incertae_sedis.2_OTU1_1-EU885561"="#E1823A",
                                 "NS-Delta-1.Incertae_sedis.2_OTU2_1-EU885632"="#E1823A",
                                 "NS-Gamma-1.2_OTU1_1-EU671146"="#E9A672",
                                 "NS-Gamma-1.Incertae_sedis_OTU1_1-EU651089"="#E9A672",
                                 "NS-Gamma-2.3.1_OTU2_1-KC469632"="#E9A672",
                                 "NT-Alpha-1.1.2.2-JN179533"="#7570B3",
                                 "Clade A-Nitrospira-sp.CTRL-LIN-TMP1"="#66C2A5",
                                 "Clade B-Nitrospira-sp.GGF-bin22"="#FC8D62",
                                 "Clade B-Nitrospira-sp.LM-bin98"="#FC8D62",
                                 "Clade B-Nitrospira-sp.LPPL-bin249"="#FC8D62",
                                 "Clade B-Nitrospira-sp.Smid-bin44"="#FC8D62"))
col_level.BS <- factor(ann.comp.BS$Taxonomy, levels = c("Nitrosolobus-multiformis-Nl1_2667636517",
                                                  "Nitrosomonas-communis-Nm44_2676397764",
                                                  #"Nitrosomonas-europaea-ATCC-19718_637427314",
                                                  "Nitrosospira-sp-17Nsp14_2671457573",
                                                  "Nitrosospira-sp_2630434854",
                                                  "Nitrosospira-sp_2636913388",
                                                  "NS-Beta-1.3_OTU1_1-EU885554",
                                                  "NS-Delta-1.Incertae_sedis.2_OTU1_1-EU885561",
                                                  "NS-Delta-1.Incertae_sedis.2_OTU2_1-EU885632",
                                                  "NS-Gamma-1.2_OTU1_1-EU671146",
                                                  "NS-Gamma-1.Incertae_sedis_OTU1_1-EU651089",
                                                  "NS-Gamma-2.3.1_OTU2_1-KC469632",
                                                  "NT-Alpha-1.1.2.2-JN179533",
                                                  "Clade A-Nitrospira-sp.CTRL-LIN-TMP1",
                                 "Clade B-Nitrospira-sp.GGF-bin22",
                                 "Clade B-Nitrospira-sp.LM-bin98",
                                 "Clade B-Nitrospira-sp.LPPL-bin249",
                                 "Clade B-Nitrospira-sp.Smid-bin44"))

tax_level.BS=levels(col_level.BS)
str(tax_level.BS)
tax_level.BS

colAnn.comp.BS <- rowAnnotation(annotation_name_gp= gpar(fontsize = 15),df=ann.comp.BS,
                             col=col.comp.ord.BS,
                             show_legend =F,
                             annotation_legend_param = list(Taxonomy = list(
                             title="Taxonomy",
                             #annotation_name_gp= gpar(fontsize = 20)
                             #ncol=1,
                             at = tax_level.BS)),
                             annotation_width=unit(c(1, 4), "cm"), 
                             gap=unit(1, "mm"))
colAnn.comp.BS

bar.ann.comp.BS <- rowAnnotation(annotation_name_gp= gpar(fontsize = 15),'Relative\nAbundance (%)' = anno_barplot(RA.comp.BS.df2,
                                                               axis_param=list(gp=gpar(fontsize = 15)),
                                                  gp = gpar(fill = c("#33A02C", "#33A02C","#33A02C","#33A02C","#33A02C","#33A02C","#33A02C","#33A02C",
                                                                     "#FF7F00","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A",
                                                                     "#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A","#6A3D9A",
                                                                     
                                                                     "#DD701E","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A","#E1823A","#E9A672","#E9A672","#E9A672","#E9A672","#E9A672","#E9A672",
                                                                     "#7570B3","#7570B3",
                                                                     
                                                                     "#66C2A5","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62",
                                                                     "#FC8D62","#FC8D62","#FC8D62","#FC8D62","#FC8D62")),
                                                  ylim=c(0,0.18),
                                                  extend = 100,
                                                  width  = unit(4, "cm"),
                                                  height = unit(6, "cm")))
#setwd('D:/Fina/INRAE_Project/microservices/DAA/glmmTMB/')
ann.fert.BS <- read.csv("BulkSoil.anno.ed_080424.csv", row.names = 1)
colours.fert.BS <- list("Fertilization"=c("CONMIN"="#E69F00",
                             "BIODYN"="#009E73",
                             "CONFYM"="#FF618C"))
#label = as.vector(c("group1","group2","group3"))
colFert.Ann.BS <- columnAnnotation(df=ann.fert.BS,
                                col=colours.fert.BS,
                                show_legend =F,
                                #labels=label,
                                show_annotation_name =F,
                                annotation_width=unit(c(1, 4), "cm"),
                                annotation_name_gp= gpar(fontsize = 20),
                                gap=unit(1, "mm"))

library(colorRamp2)
col_fun = colorRamp2(c(10, 0, -10), c("blue", "white", "red"))
lgd = Legend(title="Log2-ratio",col_fun = col_fun, labels_gp = gpar(fontsize=18),
             direction = "horizontal", 
             title_position = "topcenter",
             title_gp = gpar(fontsize=20,fontface='bold'))
#tax = Heatmap(as.matrix(ann), cluster_rows  = F)
#row_split = rep("AOB", 17)
row_split = rep("AOB", 30)
#row_split[18:22] = "AOA"
row_split[31:46] = "AOA"
#row_split[23:28] = "COMAMMOX"
row_split[47:63] = "COMAMMOX"
row_split.fa = factor(row_split, levels = c("AOB", "AOA", "COMAMMOX"))
col_split = rep("CONMIN", 5)
col_split[6:10] = "BIODYN"
col_split[11:15] = "CONFYM"
col_split.fa = factor(col_split, levels = c("BIODYN", "CONFYM", "CONMIN"))
comp.BS.HM <- Heatmap(as.matrix(rr.comp.ord.BS),
                     name = "Log2-ratio",
                     #column_title = "Bulk Soil",
                     cluster_rows  = F,
                     cluster_row_slices=F,
                     column_order=colnames(as.matrix(rr.comp.ord.BS)),
                     #column_order = order(colnames(as.matrix(rr.comp.ord))),
                     row_split = row_split.fa, 
                     column_split = col_split.fa, 
                     left_annotation = bar.ann.comp.BS,
                     right_annotation = colAnn.comp.BS,
                     #bottom_annotation = colFert.Ann.BS,
                     show_column_dend = F,
                     show_row_dend = F,
                     row_gap = unit(0.4, "cm"),
                     column_gap = unit(0, "cm"),
                     row_names_gp = gpar(fontsize = 12),
                      row_title_gp = gpar(fontsize = 25),
                      column_title_gp = gpar(fontsize = 25),
                      column_names_gp = gpar(fontsize=20),
                     border_gp = gpar(col = "black", lty = 1),
                     show_heatmap_legend = F,
                    col= col_fun)
comp.BS.HM

comp.BS.HM2 <- draw(comp.BS.HM,
                   ht_gap = unit(1, "cm"),
                    heatmap_legend_list=pd.BS,
                    align_heatmap_legend="heatmap_top")
comp.BS.HM2

setwd('/Users/arifinabintarti/Documents/France/Figures/')
png("Fig.comp.BS.HM.tiff",width=12.6,height=11,units="in",res=300)
comp.BS.HM2
dev.off()


