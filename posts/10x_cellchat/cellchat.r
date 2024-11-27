#! /data/users/jlmao/micromamba/envs/singlecell/bin/Rscript

library(Seurat)
library(CellChat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(showtext)
# setwd("/storage/jlmao/test/sc_Data/FW2024034/cellchat")

plan("multicore",workers=6)
options(future.globals.maxSize = 8 * 1024^3)

data <- readRDS("../group_s_expr.RDS")
Idents(data)<- data$celltype
data <- subset(data,idents =c("erythroid", "platelet"),invert = TRUE)
data$celltype <- data@active.ident

Idents(data) <- data$celltype
data$celltype <- Idents(data)
data$samples <- data$samplename
cellchat <- createCellChat(object = data, group.by="celltype",datatype = "RNA")

groupSize <- as.numeric(table(cellchat@idents))
groupSize
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB)# , search = "Secreted Signaling") 
# > ?subsetDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)


cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- smoothData(cellchat, adj = PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
cellchat <- readRDS("cc_S_add_AM.RDS")
groupSize <- as.numeric(table(cellchat@idents))
pdf("interaction_network.pdf")
par(xpd=T)
showtext_auto()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
showtext_auto(FALSE)
png("interaction_network.png",width=2000,height = 1200,units = "px")
par(mfrow = c(1,2),xpd=T)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
mat <- cellchat@net$weight
pdf("interaction_circle.pdf",width = 20,height = 12)
par(mfrow = c(1,2),xpd = T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  plotname <- rownames(mat)[i]
  if (plotname == "γδT/ILC"){
    plotname <- "ILC"
  }
  # pdf(paste0(plotname,"_interaction_circle.pdf"))
  
  # par(xpd = T)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  # dev.off()
}
dev.off()
png("interaction_circle.png",width=2000,height = 1200,units = "px")
par(mfrow=c(4,4),xpd = T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  plotname <- rownames(mat)[i]
  if (plotname == "γδT/ILC"){
    plotname <- "ILC"
  }
  # pdf(paste0(plotname,"_interaction_circle.pdf"))
  
  # par(xpd = T)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  # dev.off()
}
dev.off()

pathways.show <- c("GDF")

vertex.receiver = seq(1,length(levels(cellchat@idents))) # a numeric vector.
for (pathways.show in unique(cellchat@netP[["pathways"]])) {
  p1 <-
    netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  pdf(paste0(pathways.show, "_Hierarchy.pdf"))
  par(xpd = T)
  print(p1)
  dev.off()
  png(paste0(pathways.show, "_Hierarchy.png"))
  par(xpd = T)
  print(p1)
  dev.off()  
  
  p2 <-
    netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  pdf(paste0(pathways.show, "_HeatMap.pdf"))
  par(xpd = T)
  print(p2)
  dev.off()
  png(paste0(pathways.show, "_HeatMap.png"))
  par(xpd = T)
  print(p2)
  dev.off()
  
  p3 <- netAnalysis_contribution(cellchat, signaling = pathways.show)
  pdf(paste0(pathways.show, "_netAnalysis_contribution.pdf"))
  par(xpd = T)
  print(p3)
  dev.off()
  png(paste0(pathways.show, "_netAnalysis_contribution.png"))
  par(xpd = T)
  print(p3)
  dev.off()
  
}
vertex.receiver = seq(1, length(levels(cellchat@idents))) # a numeric vector.

dir.create("hierarchy_part")
setwd("hierarchy_part")
for (pathways.show in unique(cellchat@netP[["pathways"]])) {
p <- netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = c(1:9),layout = "hierarchy",arrow.size = 1,arrow.width = 1.5)
q <- netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = c(10:17),layout = "hierarchy",arrow.size = 1,arrow.width = 1.5)
pdf(paste0(pathways.show, "_hierarchy_1.pdf"),width = 14,height = 9)
par(xpd = T)
print(p)
dev.off()
png(paste0(pathways.show, "_hierarchy_1.png"),width = 900,height = 700)
par(xpd = T)
print(p)
dev.off()
pdf(paste0(pathways.show, "_hierarchy_2.pdf"),width = 14,height = 9)
par(xpd = T)
print(q)
dev.off()
png(paste0(pathways.show, "_hierarchy_2.png"),width = 900,height = 700)
par(xpd = T)
print(q)
dev.off()
}
setwd("../")
for (pathways.show in unique(cellchat@netP[["pathways"]])) {
  pairLR <-
    extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  nums <- nrow(pairLR)
  #LR.show <- pairLR.CXCL[1, ] # show one ligand-receptor pair
  # Hierarchy plot
 # vertex.receiver = seq(1, 15) # a numeric vector
  for (i in 1:nums) {
    LR.show <- pairLR[i, ]
    p1 <- netVisual_individual(
      cellchat,
      signaling = pathways.show,
      pairLR.use = LR.show,
      vertex.receiver = vertex.receiver
    )
    pdf(paste0(pathways.show,".",pairLR[i,1],"_LR_hierarchy.pdf"))
    par(xpd = T)
    print(p1)
    dev.off()
    png(paste0(pathways.show,".",pairLR[i,1],"_LR_hierarchy.png"))
    par(xpd = T)
    print(p1)
    dev.off()
    

    p2 <- netVisual_individual(
      cellchat,
      signaling = pathways.show,
      pairLR.use = LR.show,
      layout = "circle"
    )
    p3 <- netVisual_individual(
      cellchat,
      signaling = pathways.show,
      pairLR.use = LR.show,
      layout = "chord"
    )
    pdf(paste0(pathways.show,".",pairLR[i,1],"_LR_circle.pdf"))
    par(xpd = T)
    print(p2)
    dev.off()
    pdf(paste0(pathways.show,".",pairLR[i,1],"_LR_chord.pdf"))
    par(xpd = T)
    print(p3)
    dev.off()
    png(paste0(pathways.show,".",pairLR[i,1],"_LR_circle.png"))
    par(xpd = T)
    print(p2)
    dev.off()
    png(paste0(pathways.show,".",pairLR[i,1],"_LR_chord.png"))
    par(xpd = T)
    print(p3)
    dev.off()
  }
}
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver_a <- list(part_1=seq(1,length(levels(cellchat@idents))))#,part_2=seq(6,10),part_3=seq(11,15))
ori <- getwd()
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  # for (j in 1:3){
  #   setwd(paste0(ori,"/part",j))
    netVisual(cellchat, signaling = pathways.show.all[i])# ,layout = "hierarchy")# , vertex.receiver = seq(1,2), 
  # }  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  # pdf(paste0(pathways.show.all[i], "_L-R_contribution.pdf"))
  # netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  # dev.off()
}
## all
setwd(ori)
for(i in 1:length(levels(cellchat@idents))){
  p <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:length(levels(cellchat@idents))), remove.isolate = FALSE)
  from <- levels(cellchat@idents)[i]
  if (from == "γδT/ILC"){
    from <- "ILC"
  }
  pdf(paste0(from,"_to_all_bubble.pdf"))
  par(xpd = T)
  print(p)
  dev.off()
  png(paste0(from,"_to_all_bubble.png"))
  par(xpd = T)
  print(p)
  dev.off()
}
## signal
# netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:17), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
for (pathways.show in unique(cellchat@netP[["pathways"]])) {
  for (i in 1:length(levels(cellchat@idents))) {
    p <-
      netVisual_bubble(
        cellchat,
        sources.use = i,
        targets.use = c(1:length(levels(cellchat@idents))),
        signaling = pathways.show,
        remove.isolate = FALSE
      )
    from <- levels(cellchat@idents)[i]
    if (from == "γδT/ILC") {
      from <- "ILC"
    }
    pdf(paste0(from, "to_all,",pathways.show,"_bubble.pdf"))
    par(xpd = T)
    print(p)
    dev.off()
    png(paste0(from, "to_all,",pathways.show,"_bubble.png"))
    par(xpd = T)
    print(p)
    dev.off()
  }
}
## user select
pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show.all)
# pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
for (i in 1:length(levels(cellchat@idents))) {
  from <- levels(cellchat@idents)[i]
  p <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:length(levels(cellchat@idents))), pairLR.use = pairLR.use, remove.isolate = TRUE)
}
# set the order of interacting cell pairs on x-axis
# (4) Default: first sort cell pairs based on the appearance of sources in levels(object@idents), and then based on the appearance of targets in levels(object@idents)
# (5) sort cell pairs based on the targets.use defined by users
# netVisual_bubble(cellchat, targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.target = T)
# # (6) sort cell pairs based on the sources.use defined by users
# netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T)
# # (7) sort cell pairs based on the sources.use and then targets.use defined by users
# netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T)
# # (8) sort cell pairs based on the targets.use and then sources.use defined by users
# netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
for (i in c(1:length(levels(cellchat@idents)))){ 
  from <- levels(cellchat@idents)[i]
  if (from == "γδT/ILC") {
    from <- "ILC"
  }
  p <- netVisual_chord_gene(cellchat, sources.use = i, targets.use = c(1:length(levels(cellchat@idents))), lab.cex = 0.5,show.legend = F)
  pdf(paste0(from, "_to_all_net_chord.pdf"))
  par(xpd = T)
  print(p)
  dev.off()
  png(paste0(from, "_to_all_net_chord.png"))
  par(xpd = T)
  print(p)
  dev.off()
  
}
# netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
for (i in 1:length(levels(cellchat@idents))){ 
  from <- levels(cellchat@idents)[i]
  if (from == "γδT/ILC") {
    from <- "ILC"
  }
  target_use = c(1:length(levels(cellchat@idents)))
  p <- netVisual_chord_gene(cellchat, sources.use = i, targets.use = target_use,slot.name = "netP", show.legend = F)
  pdf(paste0(from, "_to_all_netP_chord.pdf"))
  par(xpd = T)
  print(p)
  dev.off()
  png(paste0(from, "_to_all_netP_chord.png"))
  par(xpd = T)
  print(p)
  dev.off()
  
}
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
## 信号基因表达分布
for (pathways.show in unique(cellchat@netP[["pathways"]])) {
  p <-
    plotGeneExpression(
      cellchat,
      signaling = pathways.show,
      enriched.only = TRUE,
      type = "violin"
    )
  pdf(paste0(pathways.show, "_gene_enrich_violin.pdf"))
  par(xpd = T)
  print(p)
  dev.off()
  png(paste0(pathways.show, "_gene_enrich_violin.png"))
  par(xpd = T)
  print(p)
  dev.off()
}
for (pathways.show in unique(cellchat@netP[["pathways"]])) {
  p <-
    plotGeneExpression(
      cellchat,
      signaling = pathways.show,
      enriched.only = FALSE,
      type = "violin"
    )
  pdf(paste0(pathways.show, "_gene_all_violin.pdf"))
  par(xpd = T)
  print(p)
  dev.off()
  png(paste0(pathways.show, "_gene_all_violin.png"))
  par(xpd = T)
  print(p)
  dev.off()
}
# plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)




# 计算和可视化网络中心性分数
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
for (pathways.show in unique(cellchat@netP[["pathways"]])) {
  pdf(paste0(pathways.show, "_signalingRole_heatmap.pdf"))
  par(xpd = T)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show)
  dev.off()
  png(paste0(pathways.show, "_signalingRole_heatmap.png"))
  par(xpd = T)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show)
  dev.off()
}

# 可视化 2D 空间中的主要发送者（sources）和接收者（targets）
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)

ggsave(filename = "all_signalingRole_scatter.pdf",width = 8,height = 8,plot = gg1)
ggsave(filename = "all_signalingRole_scatter.png",width = 300,height = 300,units = "px",plot = gg1)
for (pathways.show in unique(cellchat@netP[["pathways"]])) {
  # gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
  gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
#> Signaling role analysis on the cell-cell communication network from user's input
  # gg_a <- gg1 + gg2
  pdf(paste0(pathways.show, "_signalingRole_scatter.pdf"))
  par(xpd = T)
  print(gg2)
  dev.off()
  png(paste0(pathways.show, "_signalingRole_scatter.png"))
  par(xpd = T)
  print(gg2)
  dev.off()
}
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
for (pathways.show in unique(cellchat@netP[["pathways"]])) {
ht1 <- netAnalysis_signalingRole_heatmap(cellchat,signaling = pathways.show, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,signaling = pathways.show, pattern = "incoming")
ht <- ht1 + ht2
pdf(paste0(pathways.show, "_signalingRole_cell2cell.pdf"))
par(xpd = T)
print(ht)
dev.off()
png(paste0(pathways.show, "_signalingRole_cell2cell.png"))
par(xpd = T)
print(ht)
dev.off()
}



# ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))

library(NMF)
library(ggalluvial)

q <- selectK(cellchat, pattern = "outgoing")
png("selectK_outgoing.png")
print(q)
dev.off()
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat,pattern = "outgoing",k=2)#, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function
q1 <- selectK(cellchat, pattern = "incoming")
png("selectK_incoming.png")
print(q1)
dev.off()
nPatterns = 7
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 2)
netAnalysis_dot(cellchat, pattern = "incoming")


netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat,type = "functional",umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural",umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
saveRDS(cellchat,"cellchat_grouprful.RDS")

