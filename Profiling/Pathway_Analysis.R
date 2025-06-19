##Set Up Pathway Data Frame
Pathway <- read.delim("/Users/bopeng/Documents/76Seeds_Peptides/Cache/path.tsv", header = TRUE, sep = "\t", fill = TRUE, quote = "") 
Pathway <- t(Pathway)
colnames(Pathway) <- as.character(Pathway[2, ])
Pathwaylist <- Pathway[1:2,]
Pathway <- Pathway[-c(1:2), ]  # Remove the first row since it is now the header
Pathway <- as.data.frame(Pathway)
write.csv(Pathway, file = "/Users/bopeng/Documents/76Seeds_Peptides/Cache/Pathway_CombinedName.csv")

##Average the Values by Line Names
# Extract the middle label (e.g., L27, L28) as a new column
Pathway <- read.csv(file = "/Users/bopeng/Documents/76Seeds_Peptides/Cache/Pathway_SeperatedName.csv")
PathwayS1 <- Pathway[Pathway$Subject == "S1",]
# Average values by group
library(dplyr)
Pathway_avg <- PathwayS1 %>%
  group_by(Abbrv) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
# Clean the data frame
Pathway_avg <- as.data.frame(Pathway_avg)
rownames(Pathway_avg) <- Pathway_avg[,1]
Pathway_avg <- Pathway_avg[,-c(1:3)]  
write.csv(Pathway_avg, file = "/Users/bopeng/Documents/76Seeds_Peptides/Cache/Pathway_Avg.csv")

############
###Select Significant Pathways
data <- read.csv("/Users/bopeng/Documents/76Seeds_Peptides/Cache/Pathway_Avg.csv", row.names = 1)
data <- as.data.frame(data)
data$Category <- c(rep("Landrace",24), rep("Teosinte", 25))
Category <- data$Category
expr_data <- data[, !colnames(data) %in% "Category"]

##T-test with FDR to Select
#Initialize a vector to store p-values
p_values <- c()

#Loop over each patheay column and perform t-test
for (pathway in colnames(expr_data)) {
  group1 <- expr_data[Category=="Landrace", pathway]
  group2 <- expr_data[Category=="Teosinte", pathway]
  test_result <- t.test(group1, group2)
  p_values[pathway] <- test_result$p.value
  
#Adjust p-values
  p_adj <- p.adjust(p_values, method = "fdr")
}

#Filter the significant pathways
signif_pathways <- names(p_adj[p_adj < 0.05])
#Subset the data to only significant pathways
TTest_data <- expr_data[, signif_pathways]



##Random Forest + PLS_DA
#install.packages("randomForest")
#remotes::install_github("mixOmicsTeam/mixOmics")
library(randomForest)
library(mixOmics)

#Random_Forest
#Prepare data
expr_data <- read.csv("/Users/bopeng/Documents/76Seeds_Peptides/Pathway/Pathway_Avg.csv", row.names = 1)
expr_data <- as.data.frame(expr_data)
Category <- factor(c(rep("Landrace", 24), rep("Teosinte", 25)))

#Fit Random Forest Model
rf_model <- randomForest(x= expr_data, y=Category, importance = TRUE, ntree =1000)

#View variable importance
importance <- importance(rf_model)
importance_df <- data.frame(Pathway = rownames(importance),
                            MeanDecreaseGini = importance[, "MeanDecreaseGini"])
#Sort by importance
importance_df <- importance_df[order(-importance_df$MeanDecreaseGini),]
#Select Top Pathways
top_pathways <- importance_df$Pathway[1:50]
Random_Forest_Data <- expr_data[, top_pathways]


#PLS-DA
# Transpose data: samples in rows, variables in columns
#Prepare data
expr_data <- read.csv("/Users/bopeng/Documents/76Seeds_Peptides/Pathway/Pathway_Avg.csv", row.names = 1)
expr_data <- as.data.frame(expr_data)
Category <- factor(c(rep("Landrace", 24), rep("Teosinte", 25)))
X <- expr_data
Y <- Category

# Run PLS-DA
plsda_model <- plsda(X, Y, ncomp = 2)

# Plot sample scores
plotIndiv(plsda_model, comp = 1:2, group = Y, legend = TRUE)

# Plot Variable Importance (VIP)
vip <- selectVar(plsda_model, comp = 1)$value
vip_df <- data.frame(Pathway = rownames(vip), VIP = vip[,1])

# View top VIP pathways
vip_df <- vip_df[order(-vip_df$VIP), ]
top_pathways <- vip_df$Pathway[1:50]
PLS_DA_Data <- expr_data[,top_pathways]

###Heatmap
Heatmap <- function(filtered_data){
  data <- t(filtered_data)
  data <- as.matrix(data)
  mode(data) <- "numeric"
  #install.packages("pheatmap")  # if not already installed
  library(pheatmap)
  pheatmap(data, 
           scale = "row",       # optional: standardize each row
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           fontsize_row = 8,
           fontsize_col = 10,
           main = "Pathway Abundance Heatmap")
}

#Heatmap(TTest_data)
Heatmap(Random_Forest_Data)
Heatmap(PLS_DA_Data)
