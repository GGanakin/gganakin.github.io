library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(SingleR)
library(showtext)
library(reshape2)
library(openxlsx)
library(loupeR)



args <- commandArgs(TRUE)

allRDS <- args[1]
useRDS <- args[2]
res <- as.numeric(args[3])
ori_path <- getwd()
data <- readRDS(allRDS)

all_count <- table(Idents(data),data$sample.ident)

all_df <- dcast(as.data.frame(all_count), Var1 ~ Var2)
all_col_sums <- colSums(all_df[, -1])


colors_dutch <- c(
  '#FFC312',
  '#C4E538',
  '#12CBC4',
  '#FDA7DF',
  '#ED4C67',
  '#F79F1F',
  '#A3CB38',
  '#1289A7',
  '#D980FA',
  '#B53471',
  '#EE5A24',
  '#009432',
  '#0652DD',
  '#9980FA',
  '#833471',
  '#EA2027',
  '#006266',
  '#1B1464',
  '#5758BB',
  '#6F1E51',
  '#40407a',
  '#706fd3',
  '#f7f1e3',
  '#34ace0',
  '#33d9b2',
  '#2c2c54',
  '#474787',
  '#aaa69d',
  '#227093',
  '#218c74',
  '#ff5252',
  '#ff793f',
  '#d1ccc0',
  '#ffb142',
  '#ffda79',
  '#b33939',
  '#cd6133',
  '#84817a',
  '#cc8e35',
  '#ccae62',
  "#E64B35",
  "#4DBBD5",
  "#00A087",
  "#3C5488",
  "#F39B7F",
  "#8491B4",
  "#91D1C2",
  "#DC0000",
  "#7E6148",
  "#B09C85"
)

## T
subset_t <- readRDS(useRDS)
subset_t <-
  FindVariableFeatures(subset_t, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(subset_t)
subset_t <- ScaleData(subset_t, features = scale.genes)

subset_t <- RunPCA(subset_t, features = VariableFeatures(subset_t))

ElbowPlot(subset_t, ndims = 20, reduction = "pca")
pc.num = 1:30
##细胞聚类
subset_t <- FindNeighbors(subset_t, dims = pc.num)
subset_t <- FindClusters(subset_t, resolution = res)
table(subset_t@meta.data$seurat_clusters)

##非线性降维
#tSNE
subset_t <- RunTSNE(subset_t, dims = pc.num)
embed_tsne <- Embeddings(subset_t, 'tsne')
write.csv(embed_tsne, 'embed_tsne.csv')


#UMAP
subset_t <- RunUMAP(subset_t, dims = pc.num)

embed_umap <- Embeddings(subset_t, 'umap')
write.csv(embed_umap, 'embed_umap.csv')

diff.wilcox <- FindAllMarkers(subset_t, only.pos = TRUE)

write.csv(diff.wilcox, "all_marker.csv", row.names = F)

## plot


p <-
  DimPlot(subset_t,
          split.by = "groupName",
          label = T,
          pt.size = 1.5) + scale_color_manual(values = colors_dutch) + theme(panel.spacing = unit(0.6, "cm", data = NULL))
showtext_auto()
ggsave(
  plot = p,
  filename = "umap_group.pdf",
  width = 37,
  height = 7
)
showtext_auto(FALSE)
ggsave(
  plot = p,
  filename = "umap_group.png",
  width = 37,
  height = 7,
  type = "cairo-png"
)
p <-
  DimPlot(subset_t, label = T, pt.size = 2) + scale_color_manual(values = colors_dutch) + theme(panel.spacing = unit(0.6, "cm", data = NULL))
ggsave(
  plot = p,
  filename = "umap_all.pdf",
  width = 9,
  height = 7
)
showtext_auto(FALSE)
ggsave(
  plot = p,
  filename = "umap_all.png",
  width = 9,
  height = 7,
  type = "cairo-png"
)
p <-
  DimPlot(subset_t,
          reduction = "tsne",
          label = T,
          pt.size = 2) + scale_color_manual(values = colors_dutch) + theme(panel.spacing = unit(0.6, "cm", data = NULL))
showtext_auto()
ggsave(
  plot = p,
  filename = "tSNE_all.pdf",
  width = 9,
  height = 7
)
showtext_auto(FALSE)
ggsave(
  plot = p,
  filename = "tSNE_all.png",
  width = 9,
  height = 7,
  type = "cairo-png"
)

q <-
  DimPlot(
    subset_t,
    reduction = "tsne",
    split.by = "groupName",
    label = T,
    pt.size = 2
  ) + scale_color_manual(values = colors_dutch) + theme(panel.spacing = unit(0.6, "cm", data = NULL))

showtext_auto()
ggsave(
  plot = q,
  filename = "tSNE_group.pdf",
  width = 37,
  height = 7
)
showtext_auto(FALSE)
ggsave(
  plot = q,
  filename = "tSNE_group.png",
  width = 37,
  height = 7,
  type = "cairo-png"
)

##cell count in this part
cell_count <- table(Idents(subset_t),subset_t$sample.ident)
df <- dcast(as.data.frame(cell_count), Var1 ~ Var2)
colnames(df)[1] <- "cluster"
# df <- df[,-1]
df_transformed <- df
df_transformed[, -1] <- lapply(df[, -1], function(x) x / sum(x))
df_transformed_all <- df
df_transformed_all[, -1] <- lapply(seq_along(all_col_sums), function(i) {
  df[, i + 1] / all_col_sums[i]
})
# df_transformed_all[,-1] <- lapply(df[, -1], function(x, divisor) x / divisor, divisor = all_col_sums)

wb <- createWorkbook()
addWorksheet(wb, "count")
writeData(wb, "count", cell_count, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb, "pct")
writeData(wb, "pct", df_transformed, colNames = TRUE, rowNames = FALSE)
addWorksheet(wb,"pct_all")
writeData(wb, "pct_all", df_transformed_all, colNames = TRUE, rowNames = FALSE)
saveWorkbook(wb, "cell_count.xlsx", overwrite = TRUE)
saveRDS(subset_t,"re_expr.RDS")
create_loupe_from_seurat(subset_t,executable_path ="/storage/jlmao/software/bin/louper-linux-x64" )