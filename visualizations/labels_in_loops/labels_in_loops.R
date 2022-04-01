

setwd("~/Documents/UNLV/Year4/nmf3D/visualizations/labels_in_loops/")

library(tidyverse)


chromHMM <- read_tsv("../../get_contact_labels/chromHMM_coordonly.bed", col_names = c("chr","pos1", "pos2", "state"))
chromHMM$bp <- chromHMM$pos2 - chromHMM$pos1

hic1 <- read_tsv("../../get_contact_labels/intersect_loops_states1.bed",
                col_names = c("chr_1","pos1_1", "pos2_1", "chr_2","pos1_2", "pos2_2", "score",
                              "chr_3","pos1_3", "pos2_3","state"))

hic2 <- read_tsv("../../get_contact_labels/intersect_loops_states2.bed",
                 col_names = c("chr_1","pos1_1", "pos2_1", "chr_2","pos1_2", "pos2_2", "score",
                               "chr_3","pos1_3", "pos2_3","state"))

hic <- rbind(hic1,hic2)

hic$label <- (paste0(hic$chr_3,hic$pos1_3,hic$pos2_3,hic$state))
hic$bp <- hic$pos2_3 - hic$pos1_3


label_counts <- table((paste0(hic$chr_3,hic$pos1_3,hic$pos2_3,hic$state)))
dups <- names(label_counts)[label_counts >1]
hic$dup <- ifelse(hic$label %in% dups, T, F)


x <- hic[!duplicated(hic$label), ]

h1_x <- x[grepl("_h1",x$chr_1),]
chromHMM_h1 <-  chromHMM[grepl("_h1",chromHMM$chr),]
  
endo_x <- x[grepl("_endo",x$chr_1),]
chromHMM_endo <-  chromHMM[grepl("_endo",chromHMM$chr),]


# table(hic$state)
# table(chromHMM$state)
# table(x$state)


# hic %>% group_by(state) %>% 
#   dplyr::summarise(median(bp),sum(bp),mean(bp),n = n()) 

sum_x_h1 <- h1_x %>% group_by(state) %>% 
  dplyr::summarise(median(bp),sum(bp),mean(bp),n = n()) 

sum_chrom_h1 <- chromHMM_h1 %>% group_by(state) %>% 
  summarise(sum(bp),n = n()) 

sum_x_endo <- endo_x %>% group_by(state) %>% 
  dplyr::summarise(median(bp),sum(bp),mean(bp),n = n()) 

sum_chrom_endo <- chromHMM_endo %>% group_by(state) %>% 
  summarise(sum(bp),n = n()) 


dev.off()

par(mfrow = c(1:2))

h1_z <- data.frame(table(h1_x$state) / table(chromHMM_h1$state))
endo_z <- data.frame(table(endo_x$state) / table(chromHMM_endo$state))

h1_z <- cbind.data.frame(cell = rep("H1",nrow(h1_z)), h1_z, base = sum_x_h1$`sum(bp)` / sum_chrom_h1$`sum(bp)`)
endo_z <- cbind.data.frame(cell = rep("Endoderm",nrow(endo_z)), endo_z, base = sum_x_endo$`sum(bp)` / sum_chrom_endo$`sum(bp)`)

z <- rbind.data.frame(h1_z,endo_z)



colnames(z) <- c("Cell","State", "Percent of ChromHMM labelled regions in contact loops",
                 "Percent of ChromHMM labelled bases in contact loops")
z <- z[,-4]
zz <- gather(z, key = "type",value = "percent",-State, -Cell)

dev.off()


# ADD n= underneath x axis
# zz <- zz %>% mutate(State = factor(State),
#                     newstate = factor(State, labels = paste0(levels(State), "\n", "n=",sum_x$n )))
  
png("labels_in_loops.png",width = 1000,height = 1000, res = 200)
ggplot(zz, aes(x = State, y = percent, fill = Cell)) +
  geom_bar(stat="identity", position = position_dodge(),) +
 # geom_text(aes(label=c(rep('',18),paste0("n=",sum_x$n))),position=position_dodge(width=0.9),
        #    vjust=-0.25, size = 4) +
  labs(title ="% ChromHMM Labels present in Micro-C 5kb loops",
       x = "") +
  facet_wrap(~Cell, nrow = 2) +
  theme_minimal() + ylim(c(0,0.155)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  guides(fill="none")
dev.off()
  