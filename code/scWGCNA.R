######The construct_metacells function outputs a new Seurat object containing the aggregated expression profiles.
library(scWGCNA)

adstroke$metacell_group <- paste0(
  as.character(adstroke$celltype), '_',
  as.character(adstroke$condition)
)

genes.keep<-VariableFeatures(adstroke)


seurat_odc<-sce.big  ######

seurat_list <- list()
for(group in unique(seurat_odc$metacell_group)){
  print(group)
  cur_seurat <- subset(seurat_odc, metacell_group == group)
  cur_seurat <- cur_seurat[genes.keep,]
  cur_metacell_seurat <- scWGCNA::construct_metacells(
    cur_seurat, name=group,
    k=100, reduction='umap',
    assay='RNA', slot='data'
  )
  cur_metacell_seurat$Condition <- as.character(unique(cur_seurat$Condition))
  cur_metacell_seurat$odc_group <- as.character(unique(cur_seurat$odc_group))
  seurat_list[[group]] <- cur_metacell_seurat
}



metacell_seurat <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)])

metacell_seurat <- ScaleData(metacell_seurat, features = rownames(metacell_seurat))
metacell_seurat <- RunPCA(metacell_seurat, features=rownames(metacell_seurat))
metacell_seurat <- RunHarmony(metacell_seurat, group.by='orig.ident', dims=1:15,project.dim = F)
metacell_seurat <- RunUMAP(metacell_seurat, reduction='harmony', dims=1:15)



DimPlot(metacell_seurat, group.by='odc_group', reduction='umap', label=TRUE)

library(tidyverse)
library(WGCNA)
enableWGCNAThreads(nThreads = 64)

# how many groups are there
nclusters <- length(unique(metacell_seurat$odc_group))
genes.use <- rownames(metacell_seurat)
targets <- metacell_seurat@meta.data
group <- as.factor(metacell_seurat$Condition)

# format the expression matrix for WGCNA
datExpr <- as.data.frame(GetAssayData(metacell_seurat, assay='RNA', slot='data')[genes.use,])
datExpr <- as.data.frame(t(datExpr))


# only keep good genes:
datExpr <- datExpr[,goodGenes(datExpr)]

# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,30, by=2))


# Call the network topology analysis function for each set in turn
powerTable = list(
  data = pickSoftThreshold(
    datExpr,
    powerVector=powers,
    verbose = 100,
    networkType="signed",
    corFnc="bicor"
  )[[2]]
)


# Plot the results:
pdf("./1_Power.pdf", height=10, width=18)

colors = c("blue", "red","black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "mean connectivity",
             "Max connectivity");

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (col in 1:length(plotCols)){
  ylim[1, col] = min(ylim[1, col], powerTable$data[, plotCols[col]], na.rm = TRUE)
  ylim[2, col] = max(ylim[2, col], powerTable$data[, plotCols[col]], na.rm = TRUE)
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2))
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7

for (col in 1:length(plotCols)){
  plot(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
       xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
       main = colNames[col])
  addGrid()
  
  if (col==1){
    text(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
         labels=powers,cex=cex1,col=colors[1])
  } else
    text(powerTable$data[,1], powerTable$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[1])
  if (col==1){
    legend("bottomright", legend = 'Metacells', col = colors, pch = 20) 
  } else
    legend("topright", legend = 'Metacells', col = colors, pch = 20) 
}
dev.off()


####选择一个softpower


softPower=10

nSets = 1
setLabels = 'ODC'
shortLabels = setLabels

multiExpr <- list()
multiExpr[['ODC']] <- list(data=datExpr)

checkSets(multiExpr) # check data size


# construct network
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              maxBlockSize = 30000, ## This should be set to a smaller size if the user has limited RAM
                              randomSeed = 12345,
                              corType = "pearson",
                              power = softPower,
                              consensusQuantile = 0.3,
                              networkType = "signed",
                              TOMType = "unsigned",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0.8,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 4,
                              pamStage=FALSE,
                              detectCutHeight = 0.995, minModuleSize = 50,
                              mergeCutHeight = 0.2,
                              saveConsensusTOMs = TRUE,
                              consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")



consMEs = net$multiMEs

moduleLabels = net$colors
moduleColors = as.character(moduleLabels)
consTree = net$dendrograms[[1]]

# module eigengenes
MEs=moduleEigengenes(multiExpr[[1]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs=orderMEs(MEs)
meInfo<-data.frame(rownames(datExpr), MEs)
colnames(meInfo)[1]= "SampleID"

# intramodular connectivity
KMEs<-signedKME(datExpr, MEs,outputColumnName = "kME",corFnc = "bicor")

# compile into a module metadata table
geneInfo=as.data.frame(cbind(colnames(datExpr),moduleColors, KMEs))

# how many modules did we get?
nmodules <- length(unique(moduleColors))

# merged gene symbol column
colnames(geneInfo)[1]= "GeneSymbol"
colnames(geneInfo)[2]= "Initially.Assigned.Module.Color"

plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = paste0("ODC lineage gene dendrogram and module colors"))
