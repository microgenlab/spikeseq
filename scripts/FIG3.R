library(dplyr)
library(ggplot2)
library(circlize)


dir.create("./figure3")

c <- read.table("./inputs/omicron_gisaid_Spike.tsv", sep = "\t", header = T)


############BA.5 

s <- c[grepl("BA.5", c$genome_lineage),]
dim(s)

png('./figure3/BA5.png', res = 600, height = 20, width = 20, units = 'cm')

m <- select(s, genome_lineage, spike_lineage_group)
dim(m)
head(m)

m$spike_lineage_group <- paste("S", m$spike_lineage_group)
row.names(m) <- paste0("Fila", 1:161)
par(cex = 0.5, mar = c(0, 0, 0, 0))
chordDiagram(m,annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 1), col = "black")
}, bg.border = NA)


circos.clear()
dev.off()

############BA.4

png('./figure3/BA4.png', res = 600, height = 25, width = 35, units = 'cm')
s <- c[grepl("BA.4", c$genome_lineage),]
m <- select(s, genome_lineage, spike_lineage_group)
dim(m)
m$spike_lineage_group<- paste("S", m$spike_lineage_group)
row.names(m) <- paste0("Fila", 1:57)
par(cex = 1, mar = c(0, 0, 0, 0))
chordDiagram(m,annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 1), col = "black")
}, bg.border = NA)

circos.clear()
dev.off()

############BA.2 

s <- c[grepl("BA.2", c$genome_lineage),]

png('./figure3/BA2.png', res = 600, height = 25, width = 35, units = 'cm')

m <- select(s, genome_lineage, spike_lineage_group)
dim(m)
m$spike_lineage_group<- paste("S", m$spike_lineage_group)
row.names(m) <- paste0("Fila", 1:2097)
par(cex = 1, mar = c(0, 0, 0, 0))
chordDiagram(m,annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 1), col = "black")
}, bg.border = NA)

circos.clear()
dev.off()


############BA.1

s <- c[grepl("BA.1", c$genome_lineage),]
png('./figure3/BA1.png', res = 600, height = 25, width = 35, units = 'cm')

m <- select(s, genome_lineage, spike_lineage_group)
dim(m)
m$spike_lineage_group<- paste("S", m$spike_lineage_group)
row.names(m) <- paste0("Fila", 1:2634)
par(cex = 1, mar = c(0, 0, 0, 0))
chordDiagram(m,annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 1), col = "black")
}, bg.border = NA)

circos.clear()
dev.off()












