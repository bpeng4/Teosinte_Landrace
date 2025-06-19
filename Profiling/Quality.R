setwd("/work/benson/bpeng4/76Seeds_Peptides/data")
suppressMessages({
library("tidyverse")
library("dada2")
library("gridExtra")
})
no_of_cores = 36
metadata = "SampleMeta_Majumder_Peptide.csv"
metadata_df = read.csv(metadata)
fnFs <- metadata_df$fq1
fnRs <- metadata_df$fq2

sample.names <- metadata_df$sample

fPlot <- dada2::plotQualityProfile(fnFs, aggregate = TRUE) +
  ggtitle("Forward") +
  geom_hline(yintercept =  30, colour = "blue")
rPlot <- dada2::plotQualityProfile(fnRs, aggregate = TRUE) +
  ggtitle("Reverse") +
  geom_hline(yintercept =  30, colour = "blue")

Quality<-grid.arrange(fPlot, rPlot, nrow = 1)
ggsave("Quality.png", plot = Quality, width = 14, height = 9)


