library(stringr)
library(dplyr)
library(ggplot2)
library(ggalluvial)

########################################################################
###################STANDARD-S METHOD####################################
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
colnames(gen) <- c("id", "genome_clade", "genome_lineage")
gen$genome_lineage[26] = "AY.20"
gen$genome_lineage[43] = "BA.1"


sp <- as.data.frame(str_split_fixed(gen$id, "-", 2))
gen <- as.data.frame(cbind(sp[,2], sp[,1], gen))
colnames(gen) <- c("sample", "barcode", "id", "genome_clade", "genome_lineage")
gen$id <- NULL

std <- read.table("./inputs/stdS_tab.tsv", sep = "\t", header = T)

tab <- left_join(gen, std, by ="sample")
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

dir.create("outputs")
write.table(tab, "./outputs/metadata_standard_samples.tsv", row.names = F, quote = F, sep = "\t")

unique(tab$lineage_pangolin)
tab <- tab %>%
  mutate(spike_lineage_set = case_when(
    endsWith(lineage_pangolin, "P.2") ~ "P.2",
    endsWith(lineage_pangolin, "P.6") ~ "P.6",
    endsWith(lineage_pangolin, "C.37.1") ~ "C.37",
    endsWith(lineage_pangolin, "P.1.10") ~ "P.1",
    endsWith(lineage_pangolin, "B.1.1") ~ "B.1.1",
    endsWith(lineage_pangolin, "B.1") ~ "B.1",
    endsWith(lineage_pangolin, "B.1.351") ~ "B.1.351",
    endsWith(lineage_pangolin, "B.1.1.7") ~ "B.1.1.7",
    endsWith(lineage_pangolin, "AY.26") ~ "AY",
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



unique(tab$genome_lineage)
tab <- tab %>%
  mutate(genome_lineage_set = case_when(
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

tab$status_lineage_set <- tab$genome_lineage_set == tab$spike_lineage_set
tab$status_lineage <- tab$genome_lineage == tab$lineage_pangolin
tab$status_clade <- tab$genome_clade == tab$clade

tab$status_lineage_set <-gsub("FALSE", "Incorrect", tab$status_lineage_set)
tab$status_lineage_set <-gsub("TRUE", "Correct", tab$status_lineage_set)
tab$status_lineage <-gsub("FALSE", "Incorrect", tab$status_lineage)
tab$status_lineage <-gsub("TRUE", "Correct", tab$status_lineage)
tab$status_clade <- gsub("TRUE", "Correct", tab$status_clade)
tab$status_clade <- gsub("FALSE", "Incorrect", tab$status_clade)

#lineage group
inc <- tab[which(tab$status_lineage_set == "Incorrect"),]
s1 <-as.data.frame(table(inc$genome_lineage))
pct <- 100-(sum(s1$Freq/44*100))
cor <- tab[which(tab$status_lineage_set == "Correct"),]
mean(cor$depth)
sd(cor$depth)
min(cor$depth)
max(cor$depth)

mean(cor$pct)
sd(cor$pct)
min(cor$pct)
max(cor$pct)

#Exact lineage
inc <- tab[which(tab$status_lineage == "Incorrect"),]
s2 <- as.data.frame(table(inc$genome_lineage))
pct <- 100-(sum(s2$Freq/44*100))
table(inc$genome_lineage)

#clade
inc <- tab[which(tab$status_clade == "Incorrect"),]
s3 <- as.data.frame(table(inc$genome_lineage))
pct <- 100-(sum(s3$Freq/44*100))


list <- tab %>%
  group_by(genome_clade, genome_lineage, clade, lineage_pangolin, spike_lineage_set, genome_lineage_set, status_lineage_set) %>%
  tally() %>%
  ungroup()
list <- as.data.frame(list)
head(list)

lin <- ggplot(list, aes(axis1 = genome_lineage_set, axis2= spike_lineage_set, y = n)) +
  geom_alluvium(aes(fill = status_lineage_set), aes.bind=TRUE, width = 1/12) +
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
  ggtitle("Standard-S method", subtitle = "PANGO lineage group")
lin

library(ggsignif)
difc <- ggplot(tab, aes(status_lineage_set, pct, fill = status_lineage_set)) + 
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
  xlab("") + ylab("Completeness (%)") + scale_fill_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  guides(fill=guide_legend(title="Status")) + ggtitle("") +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T, test = "wilcox.test") 
difc

c <- tab[which(tab$status_lineage_set == "Correct"),]
mean(c$pct)
sd(c$pct)
max(c$pct)
min(c$pct)

i <- tab[which(tab$status_lineage_set == "Incorrect"),]
mean(i$pct)
sd(i$pct)
max(i$pct)
min(i$pct)

difd <- ggplot(tab, aes(status_lineage_set, as.numeric(depth), fill = status_lineage_set)) + 
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
  xlab("") + ylab("Average depth (X)") + scale_fill_manual(values = c("Correct" = "orange", "Incorrect" = "red")) +
  guides(fill=guide_legend(title="Status")) + ggtitle("") +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T, test = "wilcox.test") 
difd


library(ggpubr)

figa <- ggarrange(lin, legend ="bottom")
figb <- ggarrange(difd, difc, ncol = 2)

fig1 <- ggarrange(figa, figb, nrow = 2, labels = c("A", "B"))
fig1


########################################
###Supplementary figures S1 and S2
########################################
library(dplyr)
library(ggplot2)


merged <- read.table("./inputs/standard_sampling.tsv", sep = "\t", header = T)
head(merged)

#Get omicron samples obtained with the standard amplicon sequencing protocol
merged2 <- merged[-which(merged$seq == "general"),]

#Omicron samples
bc36 <- merged2[which(merged2$barcode == "barcode36"),]
bc48 <- merged2[which(merged2$barcode == "barcode48"),]
bc60 <- merged2[which(merged2$barcode == "barcode60"),]
bc72 <- merged2[which(merged2$barcode == "barcode72"),]
bc84 <- merged2[which(merged2$barcode == "barcode84"),]


omicron2 <- bc36
omicron3 <- bc48 
omicron4 <- bc60 
omicron5 <- bc72 
omicron6 <- bc84 


lin_omicron <- as.data.frame(rbind(omicron2, omicron3, omicron4, omicron5, omicron6))
unique(lin_omicron$pangolin)

lin_omicron <- lin_omicron %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "BA.1") ~ "Correct",
    endsWith(pangolin, "BA.1.1") ~ "Correct",
    endsWith(pangolin, "BA.2") ~ "Incorrect"))
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
bc08 <- merged3[which(merged3$barcode == "barcode08"),]
bc20 <- merged3[which(merged3$barcode == "barcode20"),]
bc32 <- merged3[which(merged3$barcode == "barcode32"),]

#beta
bc09 <- merged3[which(merged3$barcode == "barcode09"),]
bc21 <- merged3[which(merged3$barcode == "barcode21"),]
bc68 <- merged3[which(merged3$barcode == "barcode68"),]
bc80 <- merged3[which(merged3$barcode == "barcode80"),]
bc92 <- merged3[which(merged3$barcode == "barcode92"),]

#Gamma
bc33 <- merged3[which(merged3$barcode == "barcode33"),]
bc45 <- merged3[which(merged3$barcode == "barcode45"),]
bc57 <- merged3[which(merged3$barcode == "barcode57"),]
bc69 <- merged3[which(merged3$barcode == "barcode69"),]
bc81 <- merged3[which(merged3$barcode == "barcode81"),]

#Delta
bc10 <- merged3[which(merged3$barcode == "barcode10"),]
bc22 <- merged3[which(merged3$barcode == "barcode22"),]
bc34 <- merged3[which(merged3$barcode == "barcode34"),]
bc46 <- merged3[which(merged3$barcode == "barcode46"),]
bc56 <- merged3[which(merged3$barcode == "barcode56"),]

#P.6
bc24 <- merged3[which(merged3$barcode == "barcode24"),]
bc36 <- merged3[which(merged3$barcode == "barcode36"),]
bc48 <- merged3[which(merged3$barcode == "barcode48"),]
bc60 <- merged3[which(merged3$barcode == "barcode60"),]
bc72 <- merged3[which(merged3$barcode == "barcode72"),]

#Mu
bc58 <- merged3[which(merged3$barcode == "barcode58"),]
bc70 <- merged3[which(merged3$barcode == "barcode70"),]
bc82 <- merged3[which(merged3$barcode == "barcode82"),]
bc93 <- merged3[which(merged3$barcode == "barcode93"),] 

#Lambda
bc35 <- merged3[which(merged3$barcode == "barcode35"),]
bc47 <- merged3[which(merged3$barcode == "barcode47"),]
bc59 <- merged3[which(merged3$barcode == "barcode59"),]

#P.2
bc12 <- merged3[which(merged3$barcode == "barcode12"),]
bc71 <- merged3[which(merged3$barcode == "barcode71"),]
bc83 <- merged3[which(merged3$barcode == "barcode83"),]


##############
###ALFA#######
##############


alfa1 <- bc08
alfa2 <- bc20
alfa3 <- bc32 

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

beta1 <- bc09
beta2 <- bc21
beta3 <- bc68 
beta4 <- bc80
beta5 <- bc92

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

gamma1 <- bc33
gamma2 <- bc45
gamma3 <- bc57
gamma4 <- bc69
gamma6 <- bc81

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

delta1 <- bc10
delta2 <- bc22
delta3 <- bc34
delta4 <- bc46
delta5 <- bc56

lin_delta <- as.data.frame(rbind(delta1, delta2, delta3, delta4, delta5))
unique(lin_delta$pangolin)

lin_delta <- lin_delta %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "AY.26") ~ "Correct",
    endsWith(pangolin, "AY.4") ~ "Correct",
    endsWith(pangolin, "AY.48") ~ "Correct",
    endsWith(pangolin, "AY.20") ~ "Correct",
    endsWith(pangolin, "AY.80") ~ "Correct",
    endsWith(pangolin, "AY.30") ~ "Correct"))
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

##############
##P.6#########
##############


p1 <- bc24
p2 <- bc36
p3 <- bc48
p4 <- bc60
p5 <- bc72

lin_p <- as.data.frame(rbind(p1, p2, p3, p4))
unique(lin_p$pangolin)
head(lin_p)
unique(lin_p)


lin_p <- lin_p %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "B.1.1.28") ~ "Correct",
    endsWith(pangolin, "P.6") ~ "Correct",
    endsWith(pangolin, "B.1.1") ~ "Warning"))
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

##############
###MU#########
##############

mu1 <- bc58
mu2 <- bc70
mu3 <- bc82
mu4 <- bc93

lin_mu <- as.data.frame(rbind(mu1, mu2, mu3, mu4))
unique(lin_mu$pangolin)

lin_mu <- lin_mu %>%
  mutate(PANGO_lineage = case_when(
    endsWith(pangolin, "B.1.621") ~ "Correct",
    endsWith(pangolin, "B.1.621.1") ~ "Correct"))
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

##############
###LAMBDA#####
##############


lambda1 <- bc35
lambda2 <- bc47
lambda3 <- bc59


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


##############
###P.2########
##############


p21 <- bc12
p22 <- bc71
p23 <- bc83


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
##Supplementary 1
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


S1 <- ggarrange(f1a, f1b, ncol = 1, nrow = 2, common.legend = T,
                align = "hv",
                labels = c("A", "B"),
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
S1

dir.create("supplementary1")
png('./supplementary1/voc_depth_completness.png', res = 600, height = 25, width = 40, units = 'cm')
S1
dev.off()


#####################
##Supplementary 2
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


S2 <- ggarrange(f1a, f1b, ncol = 1, nrow = 2, common.legend = T,
                align = "hv",
                labels = c("A", "B"),
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
S2

dir.create("supplementary2")
png('./supplementary2/novocvoi_depth_completness.png', res = 600, height = 25, width = 40, units = 'cm')
S2
dev.off()


############################################################################################################################################################
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


mean_failed_completness <- as.data.frame(mean(failed$depth))
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
  geom_vline(xintercept = 311, linetype = "dashed")+
  geom_hline(yintercept = 87, linetype = "dashed") +
  annotate("text", x=370, y=89, label= "mean completeness") +
  annotate("text", x=305, y=75, label= "mean depth", angle = 90) +
  xlab("Average depth (X)") + ylab("Completeness (%)") +
  scale_color_manual(values = c("Correct" = "orange", "Incorrect" = "red"))
cor


fig2 <- ggarrange(stats, cor, ncol = 1, nrow = 2, common.legend = F,
                  labels = c("C", "D"), font.label = list(size = 14, color = "black", face = "bold", family = NULL, position = "top"))
fig2


FIG1f <- ggarrange(fig1, fig2, ncol = 2)
FIG1f 

dir.create("figure1")

png('./figure1/figure1.png', res = 600, height = 35, width = 40, units = 'cm')
FIG1f 
dev.off()





