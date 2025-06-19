#
rm(list = ls())
gc()
suppressMessages({
  library("tidyverse")
  library("phyloseq")
  library("rstatix")
  library("vegan")
  library("picante")
  library("ggpubr")
})
no_of_cores = 16
setwd( "/work/benson/bpeng4/76Seeds_Peptides/data")
load("intermediate/Peptide_ps.rda")
ls()
ps

ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("Mitochondria"))
ps.clean
#Filter Out Taxa Exist in At Least 25% samples
ps.clean.p0 <- filter_taxa(ps.clean, function (x) {sum(x > 0) >= 112}, prune=TRUE)
ps.clean.p0

#Complete Metadata
sample_meta = sample_data(ps.clean.p0)
sample_meta$Category[907:930]<-sample_meta$Sample_number[907:930]
sample_meta$Abbrv[907:930]<-sample_meta$Sample_number[907:930]
sample_meta$reps <- as.character(sample_meta$reps)
sam_data(ps.clean.p0)<-sample_meta
sam_data(ps)<-sample_meta

#Check Total Reads to find the Thresholds for Rarafication
OTU<-otu_table(ps.clean.p0) |>as.data.frame()
OTU$Total<- rowSums(OTU)

#Rarefication
ps.rare<-rarefy_even_depth(ps.clean.p0, sample.size = 14239,
                           rngseed = 111, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
#Change to relative abundance
ps.clean.re <- transform_sample_counts(ps.rare, function(x) x / sum(x))
ps.sample.re<-subset_samples(ps.clean.re, !ps.clean.re@sam_data$Category %in% c("FBB0","FBB16"))

#Plot Stack Bar chart for Baseline relative abundance
ps.re.FBB<-subset_samples(ps.clean.re, ps.clean.re@sam_data$Abbrv %in% c("FBB0","FBB16"))

#Set up mycolor
mycolor<-c(
  '#F781BF', '#CCEBC5', '#FFFF99', '#FF7F00', '#8DD3C7', '#7FC97F', '#999999',  
  '#CAB2D6', '#FDBF6F', '#D95F02', '#666666', '#A6761D', '#B2DF8A',  
  '#6A3D9A', '#FB9A99', '#E41A1C', '#FFFFB3', '#FDB462', '#F0027F', '#D9D9D9', '#4DAF4A', 
  '#FFFF33', '#B15928', '#FDC086', '#80B1D3', '#A6CEE3', '#BC80BD', '#1B9E77', '#E31A1C',  
  '#FCCDE5', '#386CB0', '#377EB8', '#984EA3', '#FFED6F', '#66A61E', '#E6AB02',  
  '#BF5B17', '#1F78B4', '#A65628', '#B3DE69', '#7570B3', '#E7298A', '#33A02C','#FB8072',
  '#f7754f', '#ce9032', '#97a431', '#32b166', '#35ad9c', '#38a9c5', '#3ca2f4','#a48cf4', 
  '#f45cf2', '#f66bad')

phyloseq::plot_bar(ps.re.FBB, fill = "Family") + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~sample, scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  scale_fill_manual(values=mycolor)

####Plotting Microbiome Family Level Profile For the average of each Category
##Subject1
#Subset for Each Subject
ps.S1<-subset_samples(ps.rare, ps.rare@sam_data$Sub == "S1")
#Calculate the average for each category
ps.S1.avg <- merge_samples(ps.S1, "Category")
#Convert to RA
ps.S1.avg.re <- transform_sample_counts(ps.S1.avg , function(x) x / sum(x))
#Plot
phyloseq::plot_bar(ps.S1.avg.re, fill = "Family") + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~sample, scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  scale_fill_manual(values=mycolor)+
  ggtitle("Subject One Family Level Relative Abundance Profile")


#####Alpha Diversity
## alpha diversity should be calculated before filtering on abundance and prevalence
#Remove NCs
ps.nanc<-subset_samples(ps, !ps@sam_data$Category %in% c("FBB0", "FBB16") )
tree = phyloseq::phy_tree(ps.nanc)
samp = data.frame(phyloseq::otu_table(ps.nanc))

##Alpha Diversity Among the Categories
adiv <- data.frame(
  phyloseq::estimate_richness(ps.nanc, measures = c( "Shannon", "Chao1", "Simpson", "InvSimpson", "Fisher")),
  "PD" = picante::pd(samp, tree, include.root=FALSE)[,1],
  dplyr::select(as_tibble(phyloseq::sample_data(ps.nanc)), sample, Lines, Genotype.Accesion.IDs, Abbrv, Category, Subject, reps)) %>%
  dplyr::select(-se.chao1)
#Glance the median values
adiv %>%
  group_by( Category ) %>%
  dplyr::summarise(median_Chao1 = median(Chao1),
                   median_Shannon = median(Shannon),
                   median_Simpson = median(Simpson),
                   median_InvSimpson = median(InvSimpson),
                   median_Fisher = median(Fisher),
                   median_PD = median(PD))
#Plot
adivglm<-glm(Chao1 ~ Category + Subject + reps, data = adiv)
Anova(adivglm, type = "III")


# Plotting
adiv[adiv$Category %in% c("Landrace","Teosinte"),] %>%
  pivot_longer(cols = c("Shannon", "Chao1", "Simpson", "PD"), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Shannon", "Chao1", "Simpson", "PD"))) %>%
  ggplot(aes(x = Category, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Category), height = 0, width = .2) +
  labs(x = "", y = "Alpha Diversity Measures") +
  facet_wrap(~ metric, scales = "free", ncol = 5) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 14)) +
  scale_x_discrete(labels = c("Landrace",  "Teosinte"))+
  stat_anova_test()



  aov()

  Anova(adivglm, type = "III")



##Alpha Diversity Within the Modern Maize
ps.maize<-subset_samples(ps.nanc,ps.nanc@sam_data$Category=="Modern maize")
tree = phyloseq::phy_tree(ps.maize)
samp = data.frame(phyloseq::otu_table(ps.maize))
adiv.maize <- data.frame(
  phyloseq::estimate_richness(ps.maize, measures = c("Shannon", "Chao1", "Simpson", "InvSimpson", "Fisher")),
  "PD" = picante::pd(samp, tree, include.root=FALSE)[,1],
  dplyr::select(as_tibble(phyloseq::sample_data(ps.maize)), sample, Lines, Genotype.Accesion.IDs, Abbrv, Category, Subject, reps)) %>%
  dplyr::select(-se.chao1)


#Plot
adiv_maizeglm<-glm(PD ~ Abbrv + Subject + reps, data = adiv.maize)
Anova(adiv_maizeglm, type = "III")



adiv.maize %>%
  pivot_longer(cols = c("Shannon", "Chao1", "Simpson", "PD"), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c( "Shannon", "Chao1", "Simpson", "PD"))) %>%
  ggplot(aes(x = Abbrv, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Abbrv), height = 0, width = .2) +
  labs(x = "", y = "Alpha Diversity Measures") +
  facet_wrap(~ metric, scales = "free", ncol = 1) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 14))


##Alpha Diversity Within the Landrace
ps.landrace<-subset_samples(ps.nanc,ps.nanc@sam_data$Category=="Landrace")
tree = phyloseq::phy_tree(ps.landrace)
samp = data.frame(phyloseq::otu_table(ps.landrace))
adiv.landrace <- data.frame(
  phyloseq::estimate_richness(ps.landrace, measures = c( "Shannon", "Chao1", "Simpson", "InvSimpson", "Fisher")),
  "PD" = picante::pd(samp, tree, include.root=FALSE)[,1],
  dplyr::select(as_tibble(phyloseq::sample_data(ps.landrace)), sample, Lines, Genotype.Accesion.IDs, Abbrv, Category, Sub, reps)) %>%
  dplyr::select(-se.chao1)
#Plot
adiv_landraceglm<-glm(Chao1 ~ Abbrv + Sub + reps, data = adiv.landrace)
Anova(adiv_landraceglm, type = "III")

adiv.landrace %>%
  pivot_longer(cols = c( "Shannon", "Chao1", "Simpson", "PD"), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Shannon", "Chao1", "Simpson", "PD"))) %>%
  ggplot(aes(x = Abbrv, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Abbrv), height = 0, width = .2) +
  labs(x = "", y = "Alpha Diversity Measures") +
  facet_wrap(~ metric, scales = "free", ncol = 1) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 14)) +
  stat_anova_test()


##Alpha Diversity Within Teosinte
ps.teosinte<-subset_samples(ps.nanc,ps.nanc@sam_data$Category=="Teosinte")
tree = phyloseq::phy_tree(ps.teosinte)
samp = data.frame(phyloseq::otu_table(ps.teosinte))
adiv.teosinte <- data.frame(
  phyloseq::estimate_richness(ps.teosinte, measures = c( "Shannon", "Chao1", "Simpson", "InvSimpson", "Fisher")),
  "PD" = picante::pd(samp, tree, include.root=FALSE)[,1],
  dplyr::select(as_tibble(phyloseq::sample_data(ps.teosinte)), sample, Lines, Genotype.Accesion.IDs, Abbrv, Category, Sub, reps)) %>%
  dplyr::select(-se.chao1)
#Plot
adiv.teosinte %>%
  pivot_longer(cols = c( "Shannon", "Chao1", "Simpson", "PD"), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Shannon", "Chao1", "Simpson", "PD"))) %>%
  ggplot(aes(x = Abbrv, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Abbrv), height = 0, width = .2) +
  labs(x = "", y = "Alpha Diversity Measures") +
  facet_wrap(~ metric, scales = "free", ncol = 1) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 14)) +
  stat_anova_test()


##########Beta Diversity
ps.beta= ps.sample.re

WholeSeed_otu<-otu_table(ps.beta) |> as.data.frame()

#Save the otu&meta files for permanova analysis
write.csv(WholeSeed_otu, file = "/work/benson/bpeng4/Maize_Teo_WholeSeed_16S/WholeSeed_otu.csv")


###Principle Component Analysis
ordBC <- ordinate(ps.beta, "PCoA", "bray")
ordJC <- ordinate(ps.beta, "PCoA", "jaccard")
ordUF <- ordinate(ps.beta, "PCoA", "unifrac")
ordwUF <- ordinate(ps.beta, "PCoA", "wunifrac")
smpID <- sample_data(ps.beta)$sample

# Keep first 2 vectors (latent variables, PCs) of each distance matrix
df <- rbind(data.frame(ordBC$vectors[,1:2], sample = smpID, method = 'BC'),
            data.frame(ordJC$vectors[,1:2], sample = smpID,method = 'Jaccard'),
            data.frame(ordUF$vectors[,1:2], sample = smpID,method = 'unifrac'),
            data.frame(ordwUF$vectors[,1:2], sample = smpID,method = 'wunifrac'))
# add sample_data info
df <- merge(df, data.frame(sample_data(ps.beta)), by = 'sample')

#Plotting PCoA
ggplot(data = df, aes(Axis.1,Axis.2, color = Subject, shape = Category )) + 
  geom_point() + 
  stat_ellipse() + 
  facet_wrap(~method,scales = 'free')

#For Each Microbiome
ggplot(data = df[df$Subject=="S1",], aes(Axis.1,Axis.2, color  = Category )) + 
  geom_point() + 
  stat_ellipse() + 
  facet_wrap(~method,scales = 'free')


###PERMANOVA
ps_clr <- microbiome::transform(ps.rare, "clr")
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 


dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$Status)
dispr


D_BC <- phyloseq::distance(ps.beta, "bray")
D_JC <- phyloseq::distance(ps.beta, "jaccard")
D_UF <- phyloseq::distance(ps.beta,  "unifrac")
D_wUF <- phyloseq::distance(ps.beta, "wunifrac")

# Assuming D_BC, D_JC, D_UF, and D_wUF are already defined

dist_list = list("Bray Curtis" = D_BC, 'Jaccard' = D_JC, 'Unifrac' = D_UF, 'weighted Unifrac' = D_wUF)

# Initialize an empty list to store results
res_variance_list <- list()

# Loop over each distance metric
for (i in 1:length(dist_list)) {
  # Get unique values for "Sub", "Category", and "Abbrv"
  unique_subs <- unique(sample_meta$Sub)
  unique_categories <- unique(sample_meta$Category)
  unique_abbrv <- unique(sample_meta$Abbrv)
  
  # Initialize a list for each iteration of the outer loop
  res_variance_list[[i]] <- list()
  
  # Loop through each combination of Sub, Category, and Abbrv
  for (sub in unique_subs) {
    for (category in unique_categories) {
      for (abbrv in unique_abbrv) {
        formula_str <- paste0(
          "dist_list[[i]] ~ phyloseq::sample_data(ps.beta)$Sub == '", sub, 
          "' & phyloseq::sample_data(ps.beta)$Category == '", category, 
          "' & phyloseq::sample_data(ps.beta)$Abbrv == '", abbrv, "'"
        )
        
        res_variance <- vegan::adonis2(as.formula(formula_str))
        
        # Store the result in the list for this iteration
        res_variance_list[[i]][[paste(sub, category, abbrv)]] <- res_variance
      }
    }
  }
}

# Combine all data frames from res_variance_list into one data frame
combined_df <- do.call(rbind, lapply(res_variance_list, function(x) do.call(rbind, x)))

# View the results
combined_df

write.csv(combined_df, file = "PERMANOVA_results.csv", row.names = TRUE)

# Agglomerate to phylum-level and rename
ps.clean.phylum = phyloseq::tax_glom(ps.clean.p0, taxrank = rank_names(ps.clean.p0)[2])
phyloseq::taxa_names(ps.clean.phylum) = phyloseq::tax_table(ps.clean.phylum)[,"Phylum"]
otu_table(ps.clean.phylum)[1:5,1:5]
##Melt and plot
#psmelt_phylum<-phyloseq::psmelt(ps.clean.phylum)

#Change the Category variable to a factor
ps.clean.phylum.NANC@sam_data$Category<-as.factor(ps.clean.phylum.NANC@sam_data$Category)
#Plot with p_values at the phylum level
phyloseq::psmelt(ps.clean.phylum.NANC) %>%
  ggplot(data = ., aes(x = Abbrv, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Category), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free", nrow = 2) +
  ggtitle("Phylum absolute abundance comparison among samples") +
  stat_compare_means(aes(group=Category), method = "kruskal.test")



####
ps.re.fam = phyloseq::tax_glom(ps.clean.re, taxrank = rank_names(ps.clean.re)[5])
psmelt_fam<-phyloseq::psmelt(ps.re.fam)

counts <- df %>%
  mutate(TotalCount = rowSums(.[2:ncol(df)]))

counts <- counts %>%
  mutate(across(starts_with("D"), ~./TotalCount, .names = "RelativeAbundance_{.col}"))

RA<-counts[180:356]
RA$ID<-df$ID
#write.csv(RA,"GenusTableRA1.csv")

RA<-read.csv("GenusTableRA2.csv")


df1 <- RA %>%
  group_by(ID2) %>%
  summarise(across(-c(ID),mean,.names="{.col}"))

meta$Rep<-NULL  
RA2<-merge(df1,meta)
RA3<-unique(RA2)

RA4<-RA3[-c(1)]
RAT<-as.data.frame(t(RA4[-c(178:180)]))

RAT$AVE<-rowMeans(RAT)
RAT2<-subset(RAT,AVE>=.003)
RAT2$AVE<-NULL
RA5<-as.data.frame(t(RAT2))
RA5$ID2<-RA3$ID2

RA6 <- RA5 %>%
  mutate(Subject = sub("^(.*)-.*$", "\\1", ID2),
         Type = sub("^.*-(.*)$", "\\1", ID2))

RA7<-subset(RA6,Type!="FBB0")
RA7<-subset(RA7,Type!="FBB16")

RA8<-unique(RA7)

calculate_log_fold_change <- function(subject) {
  S <- subset(RA8, Subject == subject)
  S_means <- as.data.frame(t(aggregate(S[1:35], list(S$Type), FUN = mean)))
  names(S_means) <- c("B73", "Teosinte")
  S_means <- S_means[-c(1), ]
  S_means <- transform(S_means, B73 = as.numeric(B73), Teosinte = as.numeric(Teosinte))
  S_means$LogFoldChange <- log2((S_means$Teosinte + 0.0001) / (S_means$B73 + 0.0001))
  S_Log<-as.data.frame(S_means$LogFoldChange)
  names(S_Log)<-paste0(subject)
  rownames(S_Log)<-rownames(S_means)
  return(S_Log)
}

Subjects<-unique(RA8$Subject)
Subjects<-c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S11","S12","S13")
X<-calculate_log_fold_change("S1")

for (i in Subjects){
  X<-cbind(X,calculate_log_fold_change(i))
}

X<-X[-c(1)]
#X[X==0]<-NA

XM<-as.matrix(X)
XM<-XM[rownames(XM) != "sub",]
XM<-XM[rownames(XM) != "Pedigree",]
write.csv(XM,"LogFoldChangeTeoVB73.csv")
library(RColorBrewer)
heatmap(XM,scale="column", Colv = NA, col = pal)

pal<-colorpanel(50,"darkblue","white","violetred1")
heatmap.2(XM, scale="col", col = pal, Colv = NA, trace = "none", margins=c(5,11),
          dendrogram = 'row', key.par=list(mar=c(5,4,5,2)),denscol = 'black')
mtext("Microbiome Relative Abundance 
      Comparsion between 
      Teosinte and Maize", 
      side = 3, 
      line = -4, 
      cex = 1.8)



######
#Visualize the Over-dispersion
otu_tab <- t(ps.rare@otu_table@.Data)
theme_sets <- theme_bw() +  theme(panel.border = element_blank(),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "black"))
ggplot(data=NULL, aes(x=rowMeans(otu_tab), y=apply(otu_tab, 1, sd))) + 
  geom_abline(slope = 1, intercept = 0) + geom_point() + 
  labs(title='Over-dispersion', x='Mean', y='SD') + theme_sets


#######
#Heatmap
ps.rm.na<-subset_samples(ps.clean.p0, !ps.clean.p0@sam_data$Category %in% c("FBB0","FBB16") )
ps.clean.S1<-subset_samples(ps.rm.na, ps.rm.na@sam_data$Subject=="S1")
ps.clean.S9<-subset_samples(ps.rm.na, ps.rm.na@sam_data$Subject=="S9")

#Subject One
#Subset Bacteria
ps.S1.Bac<-subset_taxa(ps.clean.S1, Family=="Bacteroidaceae")
ps.S1.Pre<-subset_taxa(ps.clean.S1, Family=="Prevotellaceae")
ps.S1.Tan<-subset_taxa(ps.clean.S1, Family=="Tannerellaceae")
ps.S1.Ahallii<-subset_taxa(ps.clean.S1, Genus=="[Eubacterium] hallii group")
ps.S1.L.Anaerostipes<-subset_taxa(ps.clean.S1, Genus=="Anaerostipes")
ps.S1.L.Dorea<-subset_taxa(ps.clean.S1, Genus=="Dorea")
ps.S1.R.Ruminococcus<-subset_taxa(ps.clean.S1, Genus=="Ruminococcus")
ps.S1.R.Faeca<-subset_taxa(ps.clean.S1, Genus=="Faecalibacterium")
ps.S1.L.Blautia<-subset_taxa(ps.clean.S1, Genus=="Blautia")
ps.S1.Bifidobacterium<-subset_taxa(ps.clean.S1, Genus=="Bifidobacterium")
ps.S1.Multiple<-subset_taxa(ps.clean.S1, Genus %in% c("[Eubacterium] hallii group",
                                                      "Anaerostipes","Dorea",
                                                      "Ruminococcus","Faecalibacterium",
                                                      "Blautia","Bifidobacterium",
                                                      "Streptococcus"))

#glommate to Genues Level
ps.S1.Bac = phyloseq::tax_glom(ps.S1.Bac, taxrank = rank_names(ps.S1.Bac)[7])
ps.S1.Pre = phyloseq::tax_glom(ps.S1.Pre, taxrank = rank_names(ps.S1.Pre)[7])
ps.S1.Tan = phyloseq::tax_glom(ps.S1.Tan, taxrank = rank_names(ps.S1.Tan)[7])
ps.S1.Multiple = phyloseq::tax_glom(ps.S1.Multiple, taxrank = rank_names(ps.S1.Multiple)[6])
#Plotting
library(microbiomeutilities)
library(viridis)
library(RColorBrewer)

heat.sample <- plot_taxa_heatmap(ps.S1.Multiple,
                                 subset.top = 50,
                                 VariableA = "Category",
                                 heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                 transformation = "log10")

keyi#Subject Nine
#Subset Bacteria
ps.S9.Bac<-subset_taxa(ps.clean.S9, Family=="Bacteroidaceae")
ps.S9.Aci<-subset_taxa(ps.clean.S9, Family=="Acidaminococcaceae")
ps.S9.Rik<-subset_taxa(ps.clean.S9, Family=="Rikenellaceae")
ps.S9.Tan<-subset_taxa(ps.clean.S9, Family=="Tannerellaceae")

#glommate to Genues Level
ps.S9.Bac = phyloseq::tax_glom(ps.S9.Bac, taxrank = rank_names(ps.S9.Bac)[7])
ps.S9.Rik = phyloseq::tax_glom(ps.S9.Rik, taxrank = rank_names(ps.S9.Rik)[7])
ps.S9.Tan = phyloseq::tax_glom(ps.S9.Tan, taxrank = rank_names(ps.S9.Tan)[7])
#Plotting
library(microbiomeutilities)
library(viridis)
library(RColorBrewer)

heat.sample <- plot_taxa_heatmap(ps.S9.Tan,
                                 subset.top = 20,
                                 VariableA = "Category",
                                 heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                 transformation = "log10")


#######
###Permanova Calculation and Visualization
##############
#glommate to Genus Level
ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])
ps.genus <- subset_samples(ps.genus, !ps.genus@sam_data$Category %in% c("FBB0", "FBB16","Modern maize") )
#####Write the Function for the Permanova Test 
ps.filter<-function(ps.genus){
  #Extract OTU/ASV table and metadata:
  otu_table <- as.data.frame(otu_table(ps.genus))
  metadata <- data.frame(sample_data(ps.genus))  # Force conversion
  taxa_table<- as.data.frame(tax_table(ps.genus))
  
  #Use the adonis function (PERMANOVA) from vegan to analyze treatment(Sample) effects
  metadata$Category <- as.factor(metadata$Category)  # Ensure it's a factor
  adonis_results <- adonis2(otu_table ~ Category, data = metadata, permutations = 999, method = "bray")
  
  #Extract p-values
  p_value <- adonis_results$`Pr(>F)`[1]  # Extract p-value for treatment
  
  #test each taxon individually:
  p_values <- apply(otu_table, 2, function(x) {
    df <- data.frame(x = x, Category = metadata$Category)
    fit <- aov(x ~ Category, data = df)
    summary(fit)[[1]][["Pr(>F)"]][1]  # Extract p-value
  })
  
  #Adjust for multiple testing (FDR correction):
  p_values_adj <- p.adjust(p_values, method = "fdr")
  
  #Mark taxa as significant if p < 0.05
  signif_taxa <- names(p_values_adj[p_values_adj < 0.05])
  
  #Filter for only the significant taxa
  ps.genus.filtered <- prune_taxa(taxa_names(ps.genus) %in% signif_taxa, ps.genus)
  
  return(ps.genus.filtered)
}
#####
#####Write the Function for Plotting Heatmap
heatmapsample<-function(ps.genus){
  # Extract OTU (abundance) table and metadata from phyloseq
  otu_table_df <- as.data.frame(otu_table(ps.genus))
  meta_df <- data.frame(sample_data(ps.genus))
  
  # Add category information to the OTU table
  otu_table_df$Sample <- meta_df$Genotype.Accesion.IDs
  otu_table_df$Category <- meta_df$Category
  
  # Aggregate by sample (mean abundance per sample)
  otu_table_avg <- otu_table_df %>%
    group_by(Category, Sample) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  
  # Convert back to matrix and ensure row names are correct
  otu_table_avg <- otu_table_avg[, !colnames(otu_table_avg) %in% "Category"]
  otu_table_avg <- column_to_rownames(otu_table_avg, var = "Sample") %>%
    as.matrix()
  
  otu_table_avg_phylo <- otu_table(otu_table_avg, taxa_are_rows = FALSE) 
  
  
  # Create a new sample metadata table
  meta_avg_phylo <- meta_df %>% distinct(Genotype.Accesion.IDs, .keep_all = TRUE)
  rownames(meta_avg_phylo) <- meta_avg_phylo$Genotype.Accesion.IDs
  meta_avg_phylo <- sample_data(meta_avg_phylo)
  
  # Ensure tax_table matches the taxa in otu_table_avg_phylo
  common_taxa <- intersect(rownames(tax_table(ps.genus)), colnames(otu_table_avg_phylo))
  tax_table_filtered <- tax_table(ps.genus)[common_taxa, ] 
  
  
  # Create a new phyloseq object with matching OTU and tax tables
  ps.genus.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)
  
  # Generate heatmap of top 50 taxa, now grouped by category
  heat.sample <- plot_taxa_heatmap(ps.genus.avg,
                                   subset.top = 50,
                                   VariableA =  c("Category"),
                                   heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                   transformation = "log10")
  heat.sample
}

###Heatmap for significant taxa calculated by Permanova for each Subject 
#Subset for each Subject
ps.S1<-subset_samples(ps.genus,ps.genus@sam_data$Subject=="S1")

ps.S1.filtered <- ps.filter(ps.S1)
heatmapsample(ps.S1.filtered)






