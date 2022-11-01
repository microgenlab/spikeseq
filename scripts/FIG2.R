
library(stringr)
library(dplyr)
library(ggplot2)
library(ggalluvial)

########################################################################
###################FAST-S METHOD#####################################
########################################################################

#get genome info

n <- read.csv("./genomes/nextclade_genomes.tsv", sep = "\t", header = T)
n <- select(n, seqName, clade)
colnames(n)[1] <- c("id")
n$id <- gsub("barcode61-CUY23-004816\t", "barcode61-CUY23-004816", n$id)
n$id <- gsub("barcode61-CUY42-007289\t", "barcode61-CUY42-007289", n$id)


l <- read.csv("./genomes/lineage_report.csv", header = T)
l <- l[,1:2]
colnames(l)[1] <- c("id")
l$id <- gsub("barcode61-CUY23-004816\t", "barcode61-CUY23-004816", l$id)
l$id <- gsub("barcode61-CUY42-007289\t", "barcode61-CUY42-007289", l$id)


p <- read.csv("./genomes/report.tsv", sep = "\t", header = T)
p$pct <-as.numeric(p$length_query/p$length_reference*100)
colnames(p)[1] <- c("id")
p <- select(p, id, pct)


gen <- left_join(n, l, by = "id")
head(gen)
dim(gen)
colnames(gen) <- c("id", "genome_clade", "genome_lineage")
gen$genome_lineage[26] = "AY.20"
gen$genome_lineage[43] = "BA.1"


sp <- as.data.frame(str_split_fixed(gen$id, "-", 2))
gen <- as.data.frame(cbind(sp[,2], sp[,1], gen))
colnames(gen) <- c("sample", "barcode", "id", "genome_clade", "genome_lineage")
gen$id <- NULL


fst <- read.table("./inputs/fastS_tab.tsv", sep = "\t", header = T) 
head(fst)
dim(fst)

tab <- left_join(gen, fst, by ="sample")
tab[is.na(tab)] <- 0
tab$depth <- as.numeric(tab$depth)
mean(tab$depth)
sd(tab$depth)
min(tab$depth)
max(tab$depth)

mean(tab$pct)
sd(tab$pct)
min(tab$pct)
max(tab$pct)

write.table(tab, "./outputs/metadata_single_samples.tsv", row.names = F, quote = F, sep = "\t")
unique(tab$lineage_pangolin)


tab <- tab[-29,]
tab <- tab %>%
  mutate(spike_lineage_pangolin_group = case_when(
    endsWith(lineage_pangolin, "P.2") ~ "P.2",
    endsWith(lineage_pangolin, "P.6") ~ "P.6",
    endsWith(lineage_pangolin, "C.37.1") ~ "C.37",
    endsWith(lineage_pangolin, "P.1.10") ~ "P.1",
    endsWith(lineage_pangolin, "B.1.1") ~ "B.1.1",
    endsWith(lineage_pangolin, "B.1") ~ "B.1",
    endsWith(lineage_pangolin, "B.1.351") ~ "B.1.351",
    endsWith(lineage_pangolin, "B.1.1.7") ~ "B.1.1.7",
    endsWith(lineage_pangolin, "AY.26") ~ "AY",
    endsWith(lineage_pangolin, "AY.62") ~ "AY",
    endsWith(lineage_pangolin, "B.1.621") ~ "B.1.621",
    endsWith(lineage_pangolin, "B.1.621.1") ~ "B.1.621",
    endsWith(lineage_pangolin, "AY.30") ~ "AY",
    endsWith(lineage_pangolin, "C.37") ~ "C.37",
    endsWith(lineage_pangolin, "AY.20") ~ "AY",
    endsWith(lineage_pangolin, "AY.30") ~ "AY",
    endsWith(lineage_pangolin, "AY.48") ~ "AY",
    endsWith(lineage_pangolin, "AY.4") ~ "AY", 
    endsWith(lineage_pangolin, "AY.20") ~ "AY",
    endsWith(lineage_pangolin, "AY.25.1") ~ "AY",
    endsWith(lineage_pangolin, "BA.1") ~ "BA.1",
    endsWith(lineage_pangolin, "BA.1.1") ~ "BA.1",
    endsWith(lineage_pangolin, "Unassigned") ~ "Unassigned"
  ))
head(tab)

#define genome lineage_group
tab <- tab %>%
  mutate(genome_lineage_group = case_when(
    endsWith(genome_lineage, "P.2") ~ "P.2",
    endsWith(genome_lineage, "P.6") ~ "P.6",
    endsWith(genome_lineage, "C.37.1") ~ "C.37",
    endsWith(genome_lineage, "P.1") ~ "P.1",
    endsWith(genome_lineage, "B.1.617.2") ~ "B.1.617.2",
    endsWith(genome_lineage, "B.1") ~ "B.1",
    endsWith(genome_lineage, "B.1.1") ~ "B.1.1",
    endsWith(genome_lineage, "B.1.351") ~ "B.1.351",
    endsWith(genome_lineage, "B.1.1.7") ~ "B.1.1.7",
    endsWith(genome_lineage, "AY.26") ~ "AY",
    endsWith(genome_lineage, "B.1.621") ~ "B.1.621",
    endsWith(genome_lineage, "B.1.621.1") ~ "B.1.621",
    endsWith(genome_lineage, "AY.30") ~ "AY",
    endsWith(genome_lineage, "C.37") ~ "C.37",
    endsWith(genome_lineage, "AY.20") ~ "AY",
    endsWith(genome_lineage, "AY.25.1") ~ "AY",
    endsWith(genome_lineage, "AY.99.2") ~ "AY",
    endsWith(genome_lineage, "AY.48") ~ "AY",
    endsWith(genome_lineage, "AY.43") ~ "AY", 
    endsWith(genome_lineage, "AY.122") ~ "AY",
    endsWith(genome_lineage, "BA.1") ~ "BA.1",
    endsWith(genome_lineage, "BA.1.1") ~ "BA.1",
    endsWith(genome_lineage, "Unassigned") ~ "Unassigned"))
head(tab)

tab$status_lineage_group <- tab$genome_lineage_group == tab$spike_lineage_pangolin_group
tab$status_lineage_group <-gsub("FALSE", "Incorrect", tab$status_lineage_group)
tab$status_lineage_group <-gsub("TRUE", "Correct", tab$status_lineage_group)


#group lineage groups
list <- tab %>%
  group_by(genome_clade, genome_lineage, clade, lineage_pangolin, spike_lineage_pangolin_group, genome_lineage_group, status_lineage_group) %>%
  tally() %>%
  ungroup()
list <- as.data.frame(list)
head(list)

#plot matching and mismatching lineage groups
lin <- ggplot(list, aes(axis1 = genome_lineage_group, axis2= spike_lineage_pangolin_group, y = n)) +
  geom_alluvium(aes(fill = status_lineage_group), aes.bind=TRUE, width = 1/12) +
  geom_stratum(width = 1/4, fill = "white", color = "black") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size =5) +
  scale_x_discrete(limits = c("Genome-based", "S gene-based"),
                   expand = c(.05, .05)) +
  scale_fill_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  labs(y = "") + guides(fill=guide_legend(title="Status")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 18, face = "bold"), 
        plot.title = element_text(size = 18)) + guides(fill=guide_legend(title="Status")) +
  ggtitle("Fast-S method", subtitle = "PANGO lineage group")
lin



tab$status_lineage_group <-gsub("TRUE", "Correct", tab$status_lineage_group)
tab$status_lineage_group <-gsub("FALSE", "Incorrect", tab$status_lineage_group)


tab$depth <- as.numeric(as.character(tab$depth))
tab$depth <- round(tab$depth, 1)
tab$pct <- as.numeric(as.character(tab$pct))
tab$pct <- round(tab$pct, 1)

library(ggsignif)
difc <- ggplot(tab, aes(status_lineage_group, pct, fill = status_lineage_group)) + 
  geom_boxplot(alpha = 0.5) + 
  theme(panel.background = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position="none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 18, face = "bold")) +
  xlab("") + ylab("Gene completeness (%)") + scale_fill_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  guides(fill=guide_legend(title="Status")) + ggtitle("") +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T) 
difc

difd <- ggplot(tab, aes(status_lineage_group, depth, fill = status_lineage_group)) + 
  geom_boxplot(alpha = 0.5) + 
  theme(panel.background = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position="none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 18, face = "bold")) +
  xlab("") + ylab("Average sequencing depth (X)") + scale_fill_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  guides(fill=guide_legend(title="Status")) + ggtitle("") +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T) 
difd


library(ggpubr)
figa <- ggarrange(lin, legend ="bottom")
figb <- ggarrange(difd, difc, ncol = 2)

fig3 <- ggarrange(figa, figb, nrow = 2, labels = c("A", "B"))
fig3


head(tab)
inc <- tab[which(tab$status_lineage_group == "Incorrect"),]
cor <- tab[which(tab$status_lineage_group == "Correct"),]

mean(cor$depth)
sd(cor$depth)
min(cor$depth)
max(cor$depth)

mean(inc$depth)
sd(inc$depth)
min(inc$depth)
max(inc$depth)

mean(cor$pct)
sd(cor$pct)
min(cor$pct)
max(cor$pct)

mean(inc$pct)
sd(inc$pct)
min(inc$pct)
max(inc$pct)


########################################
###Supplementary figures S3 and S4
########################################
library(dplyr)
library(ggplot2)

merged <- read.table("inputs/fast_sampling.tsv", sep = "\t", header = T)
head(merged)

merged2 <- merged[-which(merged$seq == "general"),]
head(merged2)

#Omicron samples
bc01 <- merged2[which(merged2$barcode == "barcode01"),]
bc13 <- merged2[which(merged2$barcode == "barcode13"),]
bc25 <- merged2[which(merged2$barcode == "barcode25"),]
bc37 <- merged2[which(merged2$barcode == "barcode37"),]
bc49 <- merged2[which(merged2$barcode == "barcode49"),]


omicron2 <- bc01
omicron3 <- bc13 
omicron4 <- bc25 
omicron5 <- bc37 
omicron6 <- bc49 


lin_omicron <- as.data.frame(rbind(omicron2, omicron3, omicron4, omicron5, omicron6))
unique(lin_omicron$pangolin)

lin_omicron <- lin_omicron %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "BA.1") ~ "Correct",
    endsWith(pangolin, "BA.1.1") ~ "Correct",
    endsWith(pangolin, "BA.1.20") ~ "Correct",
    endsWith(pangolin, "Unassigned") ~ "Incorrect",))
head(lin_omicron)

lin_omicron$sampling <- as.numeric(as.character(lin_omicron$sampling))
lin_omicron$depth <- as.numeric(as.character(lin_omicron$depth))
lin_omicron$completeness <- as.numeric(as.character(lin_omicron$completeness))


omicron1 <- ggplot(lin_omicron, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10), 
                     legend.title = element_text(size = 18), 
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("Omicron")
omicron1 

omicron2 <-  ggplot(lin_omicron, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.1, size = 3) +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10), 
                     legend.title = element_text(size = 18), 
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
omicron2


merged3 <- merged[-which(merged$seq == "omicron"),]
merged3$seq <- NULL

#alfa
bc18 <- merged3[which(merged3$barcode == "barcode18"),]
bc06 <- merged3[which(merged3$barcode == "barcode06"),]
bc30 <- merged3[which(merged3$barcode == "barcode30"),]

#beta
bc07 <- merged3[which(merged3$barcode == "barcode07"),]
bc90 <- merged3[which(merged3$barcode == "barcode90"),]
bc66 <- merged3[which(merged3$barcode == "barcode66"),]
bc78 <- merged3[which(merged3$barcode == "barcode78"),]
bc19 <- merged3[which(merged3$barcode == "barcode19"),]

#Gamma
bc67 <- merged3[which(merged3$barcode == "barcode67"),]
bc43 <- merged3[which(merged3$barcode == "barcode43"),]
bc55 <- merged3[which(merged3$barcode == "barcode55"),]
bc79 <- merged3[which(merged3$barcode == "barcode79"),]
bc31 <- merged3[which(merged3$barcode == "barcode31"),]

#Delta
bc08 <- merged3[which(merged3$barcode == "barcode08"),]
bc20 <- merged3[which(merged3$barcode == "barcode20"),]
bc32 <- merged3[which(merged3$barcode == "barcode32"),]
bc44 <- merged3[which(merged3$barcode == "barcode44"),]
bc63 <- merged3[which(merged3$barcode == "barcode63"),]

#P.6
bc22 <- merged3[which(merged3$barcode == "barcode22"),]
bc58 <- merged3[which(merged3$barcode == "barcode58"),]
bc70 <- merged3[which(merged3$barcode == "barcode70"),]
bc34 <- merged3[which(merged3$barcode == "barcode34"),]
bc46 <- merged3[which(merged3$barcode == "barcode46"),]

#Mu
bc56 <- merged3[which(merged3$barcode == "barcode56"),]
bc80 <- merged3[which(merged3$barcode == "barcode80"),]
bc68 <- merged3[which(merged3$barcode == "barcode68"),]
bc54 <- merged3[which(merged3$barcode == "barcode54"),] 

#Lambda
bc57 <- merged3[which(merged3$barcode == "barcode57"),]
bc45 <- merged3[which(merged3$barcode == "barcode45"),]
bc33 <- merged3[which(merged3$barcode == "barcode33"),]

#P.2
bc10 <- merged3[which(merged3$barcode == "barcode10"),]
bc81 <- merged3[which(merged3$barcode == "barcode81"),]
bc69 <- merged3[which(merged3$barcode == "barcode69"),]



##############
###ALFA#######
##############


alfa1 <- bc18
alfa2 <- bc06
alfa3 <- bc30 

lin_alfa <- as.data.frame(rbind(alfa1, alfa2, alfa3))
unique(lin_alfa$pangolin)


lin_alfa$Pango_lineage_assigment <- ifelse(grepl("B.1.1.7", lin_alfa$pangolin, ignore.case = T), "Correct",
                                           ifelse(grepl("", lin_alfa$pangolin, ignore.case = T), "Incorrect", "Incorrect"))


colnames(lin_alfa)[7] <- c("PANGO_lineage")
head(lin_alfa)

lin_alfa[which(lin_alfa$Pango_lineage == "Correct"),]
unique(lin_alfa$pangolin)
lin_alfa$sampling <- as.numeric(as.character(lin_alfa$sampling))
lin_alfa$depth <- as.numeric(as.character(lin_alfa$depth))
lin_alfa$completeness <- as.numeric(as.character(lin_alfa$completeness))

a <- ggplot(lin_alfa, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("Alfa") + labs(color='PANGO lineage')
a 

aa <-  ggplot(lin_alfa, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("") +
  labs(color='PANGO lineage')
aa


##############
###BETA#######
##############


beta1 <- bc07
beta2 <- bc90
beta3 <- bc66 
beta4 <- bc78
beta5 <- bc19

lin_beta <- as.data.frame(rbind(beta1, beta2, beta3, beta4, beta5))
unique(lin_beta$pangolin)
head(lin_beta)

lin_beta <- lin_beta %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "B.1.351") ~ "Correct",
    endsWith(pangolin, "B.1") ~ "Incorrect"))
head(lin_beta)


lin_beta$sampling <- as.numeric(as.character(lin_beta$sampling))
lin_beta$depth <- as.numeric(as.character(lin_beta$depth))
lin_beta$completeness <- as.numeric(as.character(lin_beta$completeness))

b <- ggplot(lin_beta, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("Beta") + labs(color='PANGO lineage')
b 

bb <-  ggplot(lin_beta, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
bb

##############
###GAMMA######
##############


gamma1 <- bc67
gamma2 <- bc43
gamma3 <- bc55
gamma4 <- bc79
gamma6 <- bc31

lin_gamma <- as.data.frame(rbind(gamma1, gamma2, gamma3, gamma4, gamma6))
unique(lin_gamma$pangolin)


lin_gamma <- lin_gamma %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "P.1.10") ~ "Correct"))
head(lin_gamma)


lin_gamma$sampling <- as.numeric(as.character(lin_gamma$sampling))
lin_gamma$depth <- as.numeric(as.character(lin_gamma$depth))
lin_gamma$completeness <- as.numeric(as.character(lin_gamma$completeness))


c <- ggplot(lin_gamma, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("Gamma") + labs(color='PANGO lineage')
c

cc <- ggplot(lin_gamma, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
cc

##############
###DELTA######
##############

delta1 <- bc08
delta2 <- bc20
delta3 <- bc32
delta4 <- bc44
delta5 <- bc63


lin_delta <- as.data.frame(rbind(delta1, delta2, delta3, delta4, delta5))
unique(lin_delta$pangolin)

lin_delta <- lin_delta %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "AY.26") ~ "Correct",
    endsWith(pangolin, "AY.4") ~ "Correct",
    endsWith(pangolin, "AY.48") ~ "Correct",
    endsWith(pangolin, "AY.20") ~ "Correct",
    endsWith(pangolin, "AY.80") ~ "Correct",
    endsWith(pangolin, "B.1.621") ~ "Incorrect"))
head(lin_delta)



lin_delta$sampling <- as.numeric(as.character(lin_delta$sampling))
lin_delta$depth <- as.numeric(as.character(lin_delta$depth))
lin_delta$completeness <- as.numeric(as.character(lin_delta$completeness))

d <- ggplot(lin_delta, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("Delta") + labs(color='PANGO lineage')
d 

dd <- ggplot(lin_delta, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
dd


#############
###P.6#######
#############

p1 <- bc22
p2 <- bc58
p3 <- bc70
p4 <- bc34
p5 <- bc46

lin_p <- as.data.frame(rbind(p1, p2, p3, p4))
unique(lin_p$pangolin)
head(lin_p)
unique(lin_p)


lin_p <- lin_p %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "P.6") ~ "Correct",
    endsWith(pangolin, "B.1.1") ~ "Incorrect"))
head(lin_p)

lin_p$sampling <- as.numeric(as.character(lin_p$sampling))
lin_p$depth <- as.numeric(as.character(lin_p$depth))
lin_p$completeness <- as.numeric(as.character(lin_p$completeness))

e <- ggplot(lin_p, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("P.6") + labs(color='PANGO lineage')
e

ee <- ggplot(lin_p, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
ee

#############
###MU########
#############

mu1 <- bc56
mu2 <- bc80
mu3 <- bc68
mu4 <- bc54

lin_mu <- as.data.frame(rbind(mu1, mu2, mu3, mu4))
unique(lin_mu$pangolin)

lin_mu <- lin_mu %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "B.1.621") ~ "Correct",
    endsWith(pangolin, "AY.30") ~ "Incorrect"))
head(lin_mu)

lin_mu$sampling <- as.numeric(as.character(lin_mu$sampling))
lin_mu$depth <- as.numeric(as.character(lin_mu$depth))
lin_mu$completeness <- as.numeric(as.character(lin_mu$completeness))



f <- ggplot(lin_mu, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("Mu") + labs(color='PANGO lineage')
f

ff <- ggplot(lin_mu, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
ff


#############
###LAMBDA####
#############

lambda1 <- bc57
lambda2 <- bc45
lambda3 <- bc33


lin_lambda <- as.data.frame(rbind(lambda1, lambda2, lambda3))
unique(lin_lambda$pangolin)

lin_lambda <- lin_lambda %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "C.37") ~ "Correct",
    endsWith(pangolin, "C.37.1") ~ "Correct"))
head(lin_lambda)

lin_lambda$sampling <- as.numeric(as.character(lin_lambda$sampling))
lin_lambda$depth <- as.numeric(as.character(lin_lambda$depth))
lin_lambda$completeness <- as.numeric(as.character(lin_lambda$completeness))


h <- ggplot(lin_lambda, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("Lambda") + labs(color='PANGO lineage')
h

hh <- ggplot(lin_lambda, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
hh

#############
###P.2#######
#############

p21 <- bc10
p22 <- bc81
p23 <- bc69


lin_p2 <- as.data.frame(rbind(p21, p22, p23))
unique(lin_p2$pangolin)

lin_p2 <- lin_p2 %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "P.2") ~ "Correct"))
head(lin_p2)


lin_p2$sampling <- as.numeric(as.character(lin_p2$sampling))
lin_p2$depth <- as.numeric(as.character(lin_p2$depth))
lin_p2$completeness <- as.numeric(as.character(lin_p2$completeness))

i <- ggplot(lin_p2, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("P.2") + labs(color='PANGO lineage')
i


ii <- ggplot(lin_p2, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage group")) +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
ii


######################
##Supplementary 2
######################

library(ggpubr)
library(grid)
library(gridExtra)

f1a <- ggarrange(a, b, c, d, omicron1,  ncol = 5, nrow = 1, common.legend = T,
                 align = "hv",
                 font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f1a <- annotate_figure(f1a, left = textGrob("Average sequencing depth (X)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                       bottom = textGrob("", gp = gpar(cex = 1.3)))
f1a 


f1b <- ggarrange(aa, bb, cc, dd, omicron2, ncol = 5, nrow = 1, legend = "none",
                 align = "hv",
                 font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "bottom"))

f1b
f1b <- annotate_figure(f1b, left = textGrob("% S gene completeness", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                       bottom = textGrob("Sampled reads", gp = gpar(cex = 1.3)))
f1b


S3 <- ggarrange(f1a, f1b, ncol = 1, nrow = 2, common.legend = T,
                align = "hv",
                labels = c("A", "B"),
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
S3

dir.create("supplementary3")
png('./supplementary3/voc_depth_completness.png', res = 600, height = 25, width = 40, units = 'cm')
S3
dev.off()


#####################
##Supplementary 3
#####################

f1a <- ggarrange(f, h, i, e,  ncol = 4, nrow = 1, common.legend = T,
                 align = "hv",
                 font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f1a <- annotate_figure(f1a, left = textGrob("Average sequencing depth (X)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                       bottom = textGrob("Sampled reads", gp = gpar(cex = 1.3)))
f1a 

f1b <- ggarrange(ff, hh, ii, ee, ncol = 4, nrow = 1, common.legend = T,
                 align = "hv",
                 font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "bottom"))

f1b

f1b <- annotate_figure(f1b, left = textGrob("% S gene completeness", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                       bottom = textGrob("Sampled reads", gp = gpar(cex = 1.3)))
f1b


S4 <- ggarrange(f1a, f1b, ncol = 1, nrow = 2, common.legend = T,
                align = "hv",
                labels = c("A", "B"),
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
S4
dir.create("supplementary4")
png('./supplementary4/novocvoi_depth_completness.png', res = 600, height = 25, width = 40, units = 'cm')
S4
dev.off()


############################################################################################################################################################
head(lin_omicron)
lin_omicron$seq <- NULL
head(lin_alfa)
head(lin_beta)
head(lin_gamma)
head(lin_delta)
head(lin_mu)
head(lin_lambda)
head(lin_p2)
head(lin_p)


tab <- as.data.frame(rbind(lin_alfa, lin_beta, lin_gamma, lin_delta, lin_omicron, lin_mu, lin_lambda, lin_p2, lin_p))

tab$completeness <- as.numeric(as.character(tab$completeness))
tab$completeness[is.na(tab$completeness)] <- 0
tab$depth <- as.numeric(as.character(tab$depth))
tab$depth[is.na(tab$depth)] <- 0
head(tab)
colnames(tab)[7] <- c("status")

failed <- tab[which(tab$status == "Incorrect"),]
mean_faield_depth <- as.data.frame(mean(failed$depth))
sd_failed_depth <- as.data.frame(sd(failed$depth))
failed_stats_depth <- as.data.frame(cbind(mean_faield_depth, sd_failed_depth))
failed_stats_depth$maximum <- max(failed$depth)
failed_stats_depth$minimum <- min(failed$depth)
rownames(failed_stats_depth) <- c("Incorrect")
colnames(failed_stats_depth) <- c("mean", "SD", "max", "min")
failed_stats_depth <- round(failed_stats_depth, 2)


mean_failed_completness <- as.data.frame(mean(failed$completeness))
sd_failed_completeness <- as.data.frame(sd(failed$completeness))
failed_stats_completeness <- as.data.frame(cbind(mean_failed_completness, sd_failed_completeness))
failed_stats_completeness$maximum <- max(failed$completeness)
failed_stats_completeness$minimum <- min(failed$completeness)
rownames(failed_stats_completeness) <- c("Incorrect")
colnames(failed_stats_completeness) <- c("mean", "SD", "max", "min")
failed_stats_completeness <- round(failed_stats_completeness, 2)


pass <- tab[which(tab$status == "Correct"),]
mean_pass_depth <- as.data.frame(mean(pass$depth))
sd_pass_depth <- as.data.frame(sd(pass$depth))
pass_stats_depth <- as.data.frame(cbind(mean_pass_depth, sd_pass_depth))
pass_stats_depth$maximum <- max(pass$depth)
pass_stats_depth$minimum <- min(pass$depth)
rownames(pass_stats_depth) <- c("Correct")
colnames(pass_stats_depth) <- c("mean", "SD", "max", "min")
pass_stats_depth <- round(pass_stats_depth, 2)

mean_pass_completness <- as.data.frame(mean(pass$completeness))
sd_pass_completeness <- as.data.frame(sd(pass$completeness))
pass_stats_completeness <- as.data.frame(cbind(mean_pass_completness, sd_pass_completeness))
pass_stats_completeness$maximum <- max(pass$completeness)
pass_stats_completeness$minimum <- min(pass$completeness)
rownames(pass_stats_completeness) <- c("Correct")
colnames(pass_stats_completeness) <- c("mean", "SD", "max", "min")
pass_stats_completeness <- round(pass_stats_completeness, 2)

summary_pass_depth <- as.data.frame(rbind(pass_stats_depth, failed_stats_depth))
summary_completeness <- as.data.frame(rbind(pass_stats_completeness, failed_stats_completeness))

library(formattable)
tab_depth <- formattable(summary_pass_depth,  
                         align = c("c",rep("r", NCOL(summary_pass_depth))),
                         list(`Indicator Name` = formatter("span", style = ~ style(color = "grey", font.weight = "bold"))))


tab_completeness <- formattable(summary_completeness, 
                                align = c("c",rep("r", NCOL(summary_completeness))),
                                list(`Indicator Name` = formatter("span", style = ~ style(color = "grey", font.weight = "bold"))))

depth_g <- ggplot(tab, aes(status, depth, fill = status)) + 
  geom_boxplot(alpha = 0.5) + 
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)) +
  guides(fill=guide_legend(title="Status"))+
  xlab("") + ylab("Average depth (X)") + scale_fill_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T) 
depth_g

comp_g <- ggplot(tab, aes(status, completeness, fill = status)) + 
  geom_boxplot(alpha = 0.5) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)) +
  guides(fill=guide_legend(title="Status")) +
  xlab("") + ylab("Completeness (%)") + scale_fill_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T)
comp_g

stats <- ggarrange(depth_g, comp_g, ncol = 2, nrow = 1, common.legend = T, font.label = list(size = 18, color = "black", face = "bold", family = NULL, position = "top"))
stats


cor <- ggplot(tab, aes(x=depth, y=completeness, color = status)) + geom_point(alpha = 0.01, size = 5) +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     legend.position = "none",
                     panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(size = 16, face = "bold"),
                     axis.text.y = element_text(size = 16, face = "bold"),
                     axis.title = element_text(size = 18, face = "bold"),
                     axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0, 500, by = 50)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  geom_vline(xintercept = 243, linetype = "dashed")+
  geom_hline(yintercept = 77, linetype = "dashed") +
  annotate("text", x=290, y=78, label= "mean completeness") +
  annotate("text", x=240, y=55, label= "mean depth", angle = 90) +
  xlab("Average depth (X)") + ylab("Completeness (%)") +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red"))
cor


fig4 <- ggarrange(stats, cor, ncol = 1, nrow = 2, common.legend = F,
                  labels = c("C", "D"), font.label = list(size = 14, color = "black", face = "bold", family = NULL, position = "top"))
fig4


FIG2f <- ggarrange(fig3, fig4, ncol = 2)
FIG2f

dir.create("figure2")
png('./figure2/figure2.png', res = 600, height = 35, width = 40, units = 'cm')
FIG2f 
dev.off()

