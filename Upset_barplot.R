library(file2meco)
library(microeco)
library(magrittr)
library(ggplot2)
library(MicrobiotaProcess)
library(UpSetR)
library(paletteer)

# to create the upsetR plot
# replace Genotype with your variable or variables and sets with the levels of your variable
upsetda <- get_upset(physeq, factorNames="Genotype")
upset(upsetda, sets=c("Nox1","Nia1Nia2","GSNOR1","Ahb1","Col0","No"), sets.bar.color = "#56B4E9",
      order.by = "freq", keep.order=T,empty.intersections = "on",nintersects = 21)


# the code from here is to create the barplot with taxonomy

# I have my data in a phyloseq object called physeq
dataset<-phyloseq2meco(physeq)

# I order the levels according to what I want
dataset$sample_table$Genotype %<>% factor(., levels = c("No","Col0","Ahb1","GSNOR1","Nia1Nia2","Nox1"))


# from here, it will do several things, creating a similar plot (that will not be used because I like more the other one)
tmp <- dataset$merge_samples(use_group = "Genotype")

t1 <- trans_venn$new(dataset = tmp)

# only show some sets with large intersection numbers
t1$data_summary <-t1$data_summary[order(-t1$data_summary$Counts),]

#This number needs to be adjusted depending on how many intersections are shown
t1$data_summary %<>% .[.[, 1] > 28, ]
t1$data_summary

# I will not use the plot but it will be used to know that the order of the intersects is correct
p1<-t1$plot_bar()
p1


t2 <- t1$trans_comm(use_frequency = F)
# t2 is a new microtable class, each part is considered a sample
class(t2)
t2$cal_abund()

#sorting
t2$taxa_abund$Rank2<-t2$taxa_abund$Rank2[,row.names(t1$data_summary)]
# transform and plot
t3 <- trans_abund$new(dataset = t2, taxrank = "Rank2",ntaxa = 12)

p2<-t3$plot_bar(bar_type = "full", others_color = "grey70",legend_text_italic = F,xtext_angle = 45, color_values=paletteer_d("RColorBrewer::Paired"),
                order_x = row.names(t1$data_summary  ))

p2

# check that the order of the intersections is the same as in the graphs above. 
# If everything is correct, remove the xtext_angle option, save the images and merge them on Inkscape