setwd("/work/benson/bpeng4/Maize_Teo_WholeSeed/")
#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions
#library("devtools")
#devtools::install_github("benjjneb/dada2")
#getwd()
#install.packages("path/to/dada2",
#                 repos = NULL,
#                 type = "source",
#                 dependencies = c("Depends", "Suggests","Imports"))
#install.packages("dada2", lib = "/work/benson/bpeng4/TEO_BMR_TOTAL")
#list.files()
#suppressMessages({
library("tidyverse")
library("dada2")
library("gridExtra")
#})
no_of_cores = 16
metadata = "Part1_Meta.csv"
metadata_df = read.csv(metadata)
metadata_df = metadata_df[,-1]
fnFs <- metadata_df$fq1
fnRs <- metadata_df$fq2

sample.names <- metadata_df$sample

#fnFs
#length(fnFs)
#fnRs
#length(fnRs)
#sample.names

#Test Quality
#dada2::plotQualityProfile(fnFs[1:20])
#dada2::plotQualityProfile(fnRs[1:20])

#fPlot <- dada2::plotQualityProfile(fnFs, aggregate = TRUE) +
#  ggtitle("Forward") +
#  geom_hline(yintercept =  30, colour = "blue")
#rPlot <- dada2::plotQualityProfile(fnRs, aggregate = TRUE) +
#  ggtitle("Reverse") +
#  geom_hline(yintercept =  30, colour = "blue")

#Quality<-grid.arrange(fPlot, rPlot, nrow = 1)
#ggsave("QualityPart1.png", plot = Quality, width = 1500, height = 1000)

filt_path <- paste0(getwd(), '/part1filtered') # don't change this 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

FORWARD_TRUNC <- 250 # determine from quality plots
REVERSE_TRUNC <- 150 # determine from quality plots

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(FORWARD_TRUNC,REVERSE_TRUNC), 
                     trimLeft=c(20, 20), maxEE=c(2,2), 
                     multithread=no_of_cores,
                     matchIDs=TRUE, compress=TRUE, 
                     verbose=TRUE)

derepFs <- derepFastq(filtFs, n = 1e+06, verbose = TRUE)
derepRs <- derepFastq(filtRs, n = 1e+06, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

errF <- learnErrors(filtFs, verbose=TRUE, multithread=no_of_cores)
errR <- learnErrors(filtRs, verbose=TRUE, multithread=no_of_cores)

dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=no_of_cores, 
               verbose=TRUE)
dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread=no_of_cores, 
               verbose=TRUE)

cat("dada-class: object describing DADA2 denoising results", sep="\n")

cat(paste(length(dadaFs[[1]]$denoised), 
          "sequence variants were inferred from", 
          length(derepFs[[1]]$uniques), 
          "input unique sequences."), sep="\n")

cat(paste("Key parameters: OMEGA_A =", dadaFs[[1]]$opts$OMEGA_A, 
          "OMEGA_C =", 
          dadaFs[[1]]$opts$OMEGA_C, "BAND_SIZE =", 
          dadaFs[[1]]$opts$BAND_SIZE), sep="\n")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

# barplot(table(nchar(getSequences(seqtab))))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=no_of_cores, verbose=TRUE)

#View(seqtab.nochim)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x)) # getUniques() gets abundance of unique sequences
# Calculate the number of reads at different steps
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

track.long <- as.data.frame(track, 
                            row.names = row.names(track)) %>% 
  gather(., key = steps, value = counts, input:nonchim, factor_key = TRUE)

Datatracking<-ggplot(track.long, aes(x = steps, y = counts, color = steps)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.3)) +
  scale_y_continuous(labels = scales::comma)
ggsave("Datatrackingpart1.png", plot = Datatracking, width = 6, height = 4, dpi = 300)

dir.create("./intermediate/", showWarnings = FALSE)
save(metadata_df, seqtab.nochim, file = "./intermediate/part1.rda")