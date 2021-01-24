
library(WGCNA)


top_abundance<-function(table,cutoff){
  talbe_rmean<-rowMeans(table)
  newtable<-table[order(rowMeans(table),decreasing = T),]
  keeprow<-ceiling(nrow(table)*cutoff)
  return(newtable[1:keeprow,])
}
int <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}
remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}
```{r}
enableWGCNAThreads(20)

otuXsample <- otu<-read.delim("pathway_int.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,header = T)


dim(otuXsample)



rownames(otuXsample)
min(otuXsample)
min(otuXsample)
colnames(otu)
colnames(otuXsample)[1:5]
dim(otuXsample)
otu_ana<-top_abundance(otuXsample,0.75)
dim(otu_ana)

sampleXotu <- as.data.frame(t(otu_ana))
dim(sampleXotu)

dim(sampleXotu)
#行是样本，列是菌群
str(sampleXotu)


rownames(sampleXotu)[1:5]


sampleXotu_f8<-sampleXotu

powers <- c(c(1:10), seq(from = 12, to=30, by=2))
powers

sft <- pickSoftThreshold(sampleXotu_f8, powerVector = powers, verbose = 5,networkType="signed")
sft
softPower <- sft$powerEstimate
softPower
# Scale-free topology fit index as a function of the soft-thresholding power
png("power.png")
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.85, col="red")

# this line corresponds to using an R^2 cut-off of h 左图：Soft Threshold (power)表示权重，纵坐标表示连接度k与p(k)的相关性。右图：Soft Threshold (power)表示权重，纵坐标表示平均连接度。一般要求k与p(k)的相关性达到0.85时的power作为β值
#综合考虑这些指标选择合适的软阈powers值，使R2可能大但也要保证连通度不能太低。看到网上有教程中提到，一般选择R2大于0.85时的power值，那么在该示例中则需选择powers=8。


abline(h=0.85,col="red")
abline(h=0.5,col="red")
abline(h=0.75,col="red")
dev.off()

pdf(file ="power.pdf",width =4,height =5.5 )
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")



abline(h=0.85,col="red")
dev.off()

getwd()





# Mean connectivity as a function of the soft-thresholding power
png("connectivity.png")
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()


pdf(file ="connectivity.pdf",width =4,height =5.5 )
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()
#save.image("1.Rdata")



softPower<-16
# check whether network approach scale free  检验选定的β值下记忆网络是否逼近 scale free 可以看出k与p(k)成负相关(相关性系数0.9),说明选择的β值能够建立基因无尺度网络。
k <- softConnectivity(sampleXotu_f8,power=softPower)
png("hist_k.png")
hist(k)
dev.off()
png("scale_free_test.png")
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()
#save.image("1.Rdata")
# obtain adjacency matrix  默认就是unsigned
adjacency = adjacency(sampleXotu_f8, power = softPower, type = "signed",corOptions = list(use = "p",method="spearman"))
#adjacency = adjacency(sampleXotu_f8, power = softPower, type = "signed")

# convert adj matrix into TOM
TOM = TOMsimilarity(adjacency)

rownames(TOM) <- colnames(TOM) <- rownames(adjacency)
rownames(TOM)
colnames(TOM)
# calculate dist between OTUs
dissTOM = 1 - TOM

## use TOM to cluster OTU
otuTree = hclust(as.dist(dissTOM), method = "average")
#save.image("2.Rdata")
# Plot the resulting clustering tree (dendrogram)
png("cluster_tree.png")
plot(otuTree, xlab="", sub="", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = F, hang = 0.04)
dev.off()

# 6.  use dynamic tree to cluster OTU into module .对刚才得到的拓扑矩阵使用相异度dissimilarity between genes对基因进行聚类，然后使用动态剪切法对树进行剪切成不同的模块（模块最小基因数为30）
minModuleSize = 5
dynamicMods = cutreeDynamic(dendro = otuTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)

write.table(table(dynamicColors), "module_color_number.txt", quote = F, sep = "\t")

otu2colors <- data.frame(OTUID = otuTree$labels, color=dynamicColors, stringsAsFactors = F)
write.table(otu2colors, "otu_colour_group.txt", quote = F, sep = "\t")

# 7. plot TOM 可以随机选择400个基因画拓扑重叠热图，图中行和列都表示单个基因，深黄色和红色表示高度的拓扑重叠

png("TOM.png")
TOMplot(dissTOM^softPower,
        otuTree,
        dynamicColors,
        main = "Network heatmap plot")
dev.off()
#save.image("2.Rdata")
# 8. calculate module eigengenes 计算每个模块的特征向量基因，为某一特定模块第一主成分基因E。代表了该模块内基因表达的整体水平
#取幂指数使dissTOM发生变化，从而使中等强度的连接在热图中更加可见



MEList = moduleEigengenes(sampleXotu_f8, colors = dynamicColors)
### 每一列是一个本征基因eigengenes 丰度代表模块基因表达水平
MEs = MEList$eigengenes

MEList$varExplained
write.table(MEList$varExplained, "pc1_varexplain.txt", quote = F, sep = "\t")
write.table(MEs, "old_eigenesgens_abu.txt", quote = F, sep = "\t")
# distance between module eigengens  

MEDiss = 1 - cor(MEs,method="spearman")
#MEDiss = 1 - cor(MEs)
# Cluster module eigengenes

METree = hclust(as.dist(MEDiss), method = "average")

png("Eigengene_adjacency_heatmap.png")
# Plot the result 特征向量基因临近热图
plotEigengeneNetworks(MEs,
                      "Eigengene adjacency heatmap",
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE,
                      xLabelsAngle = 90)
dev.off()




####### 将相关性系数大于0.8的模块合并掉，即相异性系数小于0.2:(本次合并掉2个模块)
png("module_tree.png")
plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")

# plot the cut line
abline(h=0.3, col = "red")
abline(h=0.25, col = "yellow")
abline(h=0.2, col = "blue")
dev.off()
# mrege modules
mergeModules = mergeCloseModules(sampleXotu_f8, dynamicColors, cutHeight = 0.25, verbose = 3)

# merged modules' colors #### 新模块的颜色
mergedColors = mergeModules$colors
length(table(mergedColors))
table(mergedColors)
#table(dynamicMods)

as.data.frame(table(mergedColors))
write.table(as.data.frame(table(mergedColors)), "merge_color.txt", quote = F, sep = "\t")

otu2colors$mergedcolor <- mergedColors
write.table(otu2colors[,c("OTUID","mergedcolor")], "merge_otu_color.txt", quote = F, sep = "\t")

# eigensgens for merged modules 
mergedMEs = mergeModules$newMEs

write.table(mergedMEs, "merge_eigenesgens_abu.txt", quote = F, sep = "\t")
write.table(cor(mergedMEs), "cor_mergedME.txt", quote = F, sep = "\t")



png("new_tree.png")

plotDendroAndColors(otuTree, mergedColors,
                    c("Module Color"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



pdf(file ="new_tree.pdf",width =7,height =6 )
plotDendroAndColors(otuTree, mergedColors,
                    c("Module Color"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()





KIM = intramodularConnectivity(adjacency, otu2colors$mergedcolor, scaleByMax= FALSE)
otu_mean<-as.data.frame(rowMeans(otu_ana))
KIM$mean_abundance<-rowMeans(otu_ana)

write.csv(KIM,"gene_connectivity.csv")




HubGenes <-chooseTopHubInEachModule(sampleXotu_f8,mergedColors,power=softPower,type="signed")

write.table (HubGenes,file = "HubGenes_of_each_module.xls",quote=F,sep='\t')


dim(mergedMEs)
dim(sampleXotu_f8)

out_color<-otu2colors[,c("OTUID","mergedcolor")]
module_mean<-matrix(NA,nrow(mergedMEs),ncol(mergedMEs))
rownames(module_mean)<-rownames(sampleXotu_f8)
colnames(module_mean)<-unique(out_color[,2])



otu_ana_data<-as.data.frame(t(otu_ana))


for (sample_id in rownames(module_mean)){
     for (color in colnames(module_mean)){
       module_mean[sample_id,color]=mean(as.numeric(otu_ana_data[sample_id,out_color[,1][out_color[,2]==color]]))

     }
}
write.table(module_mean, "merge_mean_abu.txt", quote = F, sep = "\t")



