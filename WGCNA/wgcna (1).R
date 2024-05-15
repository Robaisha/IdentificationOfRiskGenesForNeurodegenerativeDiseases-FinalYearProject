#Step 1: Loading WGCNA

install.packages("BiocManager")
BiocManager::install("GO.db", force = TRUE)
install.packages("WGCNA")
BiocManager::install("WGCNA", force= TRUE)
library(GO.db)
#Loading WGCNA
library(WGCNA)
#Setting string not as factor
options(stringsAsFactors = FALSE)
#Enable multithread
enableWGCNAThreads()

#Step 2: Preparing Dataset
# Reading the raw data (rows are the samples and columns the genes)
expressiondata = read.csv("D:/FYP/output_matrix.txt", sep = "\t")

# Create a new format expression data - remove gene name column (ID_REF)
expression1 = as.data.frame(expressiondata[, -c(1)])
expression1 = t(expression1)
expression = as.data.frame(expression1[0:401, 0:5000 ])


# Column 1 - gene names
colnames(expression) = expressiondata$ID_REF
rownames(expression) = names(expression)[-c(1)]

# Group data in a dendrogram to check for outliers
sampleTree = hclust(dist(expression), method = "average")

pdf(file = "D:/FYP/sampleClustering.pdf", width = 12, height = 9)

par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line showing the cut-off
abline(h = 130, col = "red")

dev.off()

# Determine clusters below the line
clust = cutreeStatic(sampleTree, cutHeight = 31000, minSize = 10)

table(clust)

# Assuming cluster 1 contains the samples.
keepSamples = (clust == 1)
expression = expression[keepSamples, ]

# Now, 'expression' contains only the data from the selected cluster.
# Display the dimensions of the filtered expression data
nGenes = ncol(expression)
nSamples = nrow(expression)

print(nGenes)
print(nSamples)

# Read trait data from the provided file
traitData <- read.csv("D:/FYP/output_trait.txt", sep = "\t", check.names = FALSE)


# Extract sample names (column names) from the expression data
Samples <- rownames(expression)

print(head(traitData))

# Sample names should match column names in traitData starting from the second column
sampleNames <- colnames(traitData)[-1]  # Excluding 'ID_REF'

sampleIDs <- traitData[1, -1]  # Excludes 'ID_REF' which is not a sample ID

# Extracting corresponding tissue types
tissueTypes <- traitData[1, -1]  # Excludes the first column which is just a label, not data

# Constructing the 'datTraits' data frame
datTraits <- data.frame(TissueType = as.character(tissueTypes), stringsAsFactors = FALSE)

# Converting 'TissueType' to a factor for WGCNA analysis
datTraits$TissueType <- as.factor(datTraits$TissueType)

# Checking the constructed 'datTraits'
print(head(datTraits))

library(RColorBrewer)
distinctTissueTypes <- unique(datTraits$TissueType)
colorPalette <- brewer.pal(n = length(distinctTissueTypes), name = "Set1")

# Map tissue types to colors
tissueColors <- setNames(colorPalette, distinctTissueTypes)

# Apply the color mapping to your data frame
traitColors <- sapply(datTraits$TissueType, function(t) tissueColors[t])

# Now, 'traitColors' contains the color coding for each sample based on its tissue type
print(traitColors)

pdf("D:/FYP/SampleDendrogramAndTraitHeatmap77.pdf", width = 12, height = 8)

# Regrouping samples using hierarchical clustering
sampleTree2 <- hclust(dist(expression), method = "average")

print(head(traitColors))

# Plotting the dendrogram and trait heatmap
plotDendroAndColors(sampleTree2, colors = traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

if (length(traitColors) < length(sampleTree2$labels)) {
  traitColors <- rep(traitColors, length.out = length(sampleTree2$labels))
} else if (length(traitColors) > length(sampleTree2$labels)) {
  traitColors <- traitColors[1:length(sampleTree2$labels)]
}
# Closing the PDF file
dev.off()

#Step 3: Creating the Network
powers <- c(1:10)
sft = pickSoftThreshold(expression, powerVector = powers, networkType = "unsigned")

pdf("D:/FYP/networkAnalysisOutput.pdf")

# Plotting Scale Free Topology Fit Indices
par(mfrow = c(1, 2))
cex1 = 0.9

# Scale Independence Plot
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

# Mean Connectivity Plot
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

dev.off()

# Continue with the analysis after plotting
softPower = 6
adjacency = adjacency(expression, power = softPower, type = "unsigned")
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

#Step 4: Grouping Genes in Modules

hc <- hclust(as.dist(dissTOM), method = "average")

# Cut the dendrogram to obtain modules
min_module_size <- 30
dynamicColors <- cutreeDynamic(dendro = hc, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = min_module_size)

# Convert module colors to a named vector
module_colors <- labels2colors(dynamicColors)

MEList = moduleEigengenes(expression, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate the module dissimilarity based on eigengenes correlation
MEDiss = 1 - cor(MEs)

# Cluster the module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

pdf("D:/FYP/clustering_of_module_eigengenes.pdf")

# Plot the result
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
# Add a cut-off line
MEDissThres = 0.6
abline(h = MEDissThres, col = "red")

# Close the plotting device
dev.off()

library(dynamicTreeCut)

print(merge)
class(merge)

merge_result <- mergeCloseModules(expression, dynamicColors, cutHeight = 1, verbose = 3)

# Access merged module colors and eigengenes from the result
mergedColors <- merge_result$colors
mergedMEs <- merge_result$newMEs


geneTree <- hclust(as.dist(dissTOM), method = "average")

pdf("D:/FYP/GENEDENDRO.pdf", width = 9, height = 6)

# Plot the gene dendrogram with dynamic and merged module colors
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Close the PDF device
dev.off()

#Step 5: Associating Modules and Phenotypes

# Convert tissue names to a factor and then to numeric to facilitate correlation analysis
datTraits$TissueType <- as.numeric(factor(datTraits$TissueType))

# Defining the number of genes and samples
nGenes <- ncol(expression)
nSamples <- nrow(expression)

# Recalculating MEs with label colors
MEs0 <- moduleEigengenes(expression, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs0)

dim(MEs)  # Check dimensions of MEs
dim(datTraits)  # Check dimensions of datTraits


# Compute correlations and p-values between module eigengenes and the trait (now numeric)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

pdf("D:/FYP/ModuleTraitRelationships.pdf", width = 8, height = 4)

# Prepare text for the heatmap: combine correlation values and their p-values
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Plotting the heatmap
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-trait relationships")

dev.off()

TissueTypess <- as.numeric(factor(datTraits$TissueType))
names(TissueTypess) = "TissueType"

#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(expression, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(expression, TissueTypess, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <- "GS.TissueType"

# Example, assuming you have a specific tissue type name available
tissueTypeName <- "SpecificTissueType"
names(geneTraitSignificance) <- paste("GS", tissueTypeName, sep=".")

names(GSPvalue) <- "p.GS.TissueType"

# If you had a specific name for the tissue type analyzed, you would use:
tissueTypeNameGSP <- "SpecificTissueType" 
names(GSPvalue) <- paste("p.GS", tissueTypeNameGSP, sep=".")

pdf("D:/FYP/Module membership vs. gene significance.pdf", width = 7, height = 7)
module = "pink" 
column = match(module, modNames)
moduleGenes = module_colors==module

#sizeGrWindow(7, 7)
par(mfrow = c(1,1))


moduleMembershipMMturquoise <- geneModuleMembership$MMturquoise

# Extract the gene significance values
geneSignificanceSpecificTissue <- geneTraitSignificance$GS.SpecificTissueType

# Check lengths to ensure they match
length(moduleMembershipMMturquoise)  # Check length of moduleMembershipMMturquoise
length(geneSignificanceSpecificTissue) 

# Plot the relationship
plot(moduleMembershipMMturquoise,
     geneSignificanceSpecificTissue,
     xlab = "Module Membership in MM4",
     ylab = "Gene significance for Specific Tissue Type",
     main = "Module membership (MM4) vs. gene significance\n",
     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'pink')

dev.off()

write.csv(dynamicColors, file = "D:/FYP/gene_module_assignments.csv", row.names = TRUE)

#HUB GENES
# Ensure that moduleLabels are factorized if not already
moduleLabels <- labels2colors(dynamicColors)
moduleLabels <- factor(moduleLabels)

# Recalculate MEs if necessary
MEs <- moduleEigengenes(expression, colors = moduleLabels)$eigengenes

# Ensure adjacency matrix is correct
adjacency <- adjacency(expression, power = softPower, type = "unsigned")

# Calculate total connectivity for each gene
totalConnectivity <- rowSums(adjacency)

# Assuming moduleLabels is a vector of module assignments for each gene
modules <- unique(moduleLabels)

# Initialize a list to store connectivity data for each module
moduleConnectivity <- list()

# Loop through each module and calculate intramodular connectivity
for(module in modules) {
  moduleGenes <- which(moduleLabels == module)
  moduleAdjacency <- adjacency[moduleGenes, moduleGenes]
  moduleConnectivity[[module]] <- rowSums(moduleAdjacency)
}

# Optionally, combine all module connectivity data into a single data frame
connectivityData <- do.call(rbind, lapply(names(moduleConnectivity), function(x) data.frame(GeneID = names(moduleConnectivity[[x]]), Module = x, Connectivity = moduleConnectivity[[x]])))

# Step 7: Identify top hub genes in each module
# We have the 'connectivityData' dataframe with GeneID, Module, and Connectivity from the previous steps.

# Convert 'Module' to factor to handle categorical data efficiently
connectivityData$Module <- as.factor(connectivityData$Module)

# Sorting the data to find top hub genes within each module
# Assuming we want the top 5 hub genes from each module
topHubGenes <- by(data = connectivityData, connectivityData$Module, function(x) {
  topGenes <- x[order(-x$Connectivity), ][1:5, ]
  return(topGenes)
})

# Optional: Convert the list to a data frame for easier viewing and further analysis
topHubGenesDF <- do.call(rbind, topHubGenes)
row.names(topHubGenesDF) <- NULL  # Clean up row names for a clearer presentation

# Print out the top hub genes for review
print(topHubGenesDF)

# Step 8: Output the results to a CSV file
write.csv(topHubGenesDF, "D:/FYP/top_hub_genes.csv", row.names = FALSE)

# Read the hub genes output data
hubGenesData <- read.csv("D:/FYP/top_hub_genes.csv", header = TRUE)

# Ensure GeneIDs match (assuming 'GeneID' columns exist in both data frames)
# This step assumes there's a way to map or validate IDs between the two sources
hubGenesData$ID_REF <- expressiondata$ID_REF[match(hubGenesData$GeneID, expressiondata$ID_REF)]

# Save the corrected data
write.csv(hubGenesData, "D:/FYP/corrected_top_hub_genes.csv", row.names = FALSE)
