cat("library packages...\n")
library(Seurat)
library(SeuratObject)
library(CiteFuse)
library(cowplot)
library(dplyr)
library(future)
library(ggplot2)
library(JXYlightHippo)
library(SingleCellExperiment)
library(SingleR)
library(SummarizedExperiment)
library(tidyr)
library(tidyverse)
##修改点
dir_outs_name <- "BGI_PBMC2_test_QC_outs/"
save_dir <- paste0("/database/jiangxinyu/Result/",dir_outs_name)
dir.create(paste0(save_dir))
cat("define Split_HTODemux_quantile...\n")
Split_HTODemux_quantile <- function(mat,sample_num=6){

  Coassay <- CreateSeuratObject(counts = mat[["Gene Expression"]])
  ADT_HTO <- mat$`Antibody Capture`

  HTO <- ADT_HTO[1:sample_num,]
  Coassay[["HTO"]] <- CreateAssayObject(counts = HTO)

  ADT <- ADT_HTO[(sample_num+1):143,]
  Coassay[["ADT"]] <- CreateAssayObject(counts = ADT)

  Coassay <- NormalizeData(Coassay, assay = "HTO", normalization.method = "CLR")
  positive.quantile <- c(0.95,0.96,0.97,0.98,0.99)
  i=1
  singlets_rate <- NULL
  best.quantile <- NULL
  ##取Coassay的测试子集
  cat("取Coassay的测试子集\n")
  Coassay@active.assay <- "HTO"
  Coassay_subset <- Coassay[,sample(1:ncol(Coassay),round(ncol(Coassay)*0.6),replace =  F)]


  for (quantile in positive.quantile) {
    #check HTODemux function
    Coassay_subset <- HTODemux(Coassay_subset, assay = "HTO", positive.quantile = quantile)
    Coassay_table <-table(Coassay_subset$HTO_classification.global)
    print(Coassay_table)
    singlets_rate <- c(singlets_rate,as.double(100*Coassay_table[3]/sum(Coassay_table)))
  }
  singlets_rate
  library(ggplot2)
  ggplot()+
    geom_line(aes(x=positive.quantile, y=singlets_rate))
  ##找到最优quantile值
  cat("找到最优quantile值\n")
  best.quantile <- positive.quantile[which.max(singlets_rate)]
  print(paste0("best.quantile:",best.quantile))

  Coassay <- HTODemux(Coassay, assay = "HTO", positive.quantile = best.quantile)
  cat("拆分结果：\n")
  Coassay_table <-table(Coassay$HTO_classification.global)
  print(Coassay_table)

  siglets_rate <- as.double(100*Coassay_table[3]/sum(Coassay_table))
  doublets_rate <- as.double(100*Coassay_table[1]/sum(Coassay_table))
  negtive_rate <- as.double(100*Coassay_table[2]/sum(Coassay_table))
  print(paste0("siglets_rate:",siglets_rate))
  print(paste0("doublets_rate:",doublets_rate))
  print(paste0("negtive_rate:",negtive_rate))

  return(Coassay)
}
cat("Done!")
time1 <- Sys.time()
cat("Step1:read data...\n")
##修改点
Data_dir = '/database/jiangxinyu/cellranger/cellranger_outs/BGI_test_202205/outs/filtered_feature_bc_matrix'
sample_num <- 6
Data_mat <- Read10X(Data_dir)

cat("Step2:split data of Seurat...\n")
Coassay_Data <- Split_HTODemux_quantile(Data_mat,sample_num)

cat("Step3:cluster and annotation...\n")
Singlet_Coassay_Data<- Coassay_Data[,which(Coassay_Data@meta.data[["HTO_classification.global"]] == "Singlet")]

Singlet_Coassay_Data@active.assay <- "RNA"

# Select the top 1000 most variable features
Singlet_Coassay_Data <- FindVariableFeatures(Singlet_Coassay_Data)

# Scaling RNA data, we only scale the variable features here for efficiency
Singlet_Coassay_Data <- ScaleData(Singlet_Coassay_Data, features = VariableFeatures(Singlet_Coassay_Data))

# Run PCA
Singlet_Coassay_Data <- RunPCA(Singlet_Coassay_Data,features = VariableFeatures(Singlet_Coassay_Data))

# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
Singlet_Coassay_Data <- FindNeighbors(Singlet_Coassay_Data, reduction = "pca", dims = 1:20)
Singlet_Coassay_Data <- FindClusters(Singlet_Coassay_Data, resolution = 0.2, verbose = FALSE)
Singlet_Coassay_Data <- RunUMAP(Singlet_Coassay_Data,dims = 1:20)
##修改点
pdf(paste0(save_dir,"Seurat_Dimplot_umap.pdf"))
print(DimPlot(Singlet_Coassay_Data,reduction = "umap"))
dev.off()

cat("SingleR annotation:\n")

load("/home/jiangxinyu/R/R_script/biomamba_single-sample/SingR annotation lab/BlueprintEncode_bpe.se_human.RData")
load("/home/jiangxinyu/R/R_script/biomamba_single-sample/SingR annotation lab/HumanPrimaryCellAtlas_hpca.se_human.RData")
Data_for_SingleR <- GetAssayData(Singlet_Coassay_Data, slot="data")  ###提取RNA的转录表达数据
clusters <- Singlet_Coassay_Data@meta.data$seurat_clusters
###把scRNA数据中的seurat_clusters提取出来，注意是因子类型
cellpred <- SingleR(test = Data_for_SingleR, ref = list(bpe.se,hpca.se), labels = list(bpe.se$label.main, hpca.se$label.main))
Singlet_Coassay_Data@meta.data$celltype <- cellpred$labels
##修改点
pdf(paste0(save_dir,"Seurat_celltype_umap.pdf"))
print(DimPlot(Singlet_Coassay_Data, group.by = c("celltype"),reduction = "umap",label = T))
dev.off()
##合并同类细胞类型
cat("Step4:merge same celltype and calculate the proportion of celltype...\n")
celltype_Data <- Singlet_Coassay_Data@meta.data[["celltype"]]
celltype_Data[celltype_Data=="B_cell"] <- "B-cells"
celltype_Data[celltype_Data=="NK_cell"] <- "NK cells"
celltype_Data[celltype_Data=="Monocyte"] <- "Monocytes"
cat("merge celltype number:\n")
table(celltype_Data)
Singlet_Coassay_Data@meta.data[["celltype_merge"]] <- celltype_Data

pdf(paste0(save_dir,"Seurat_celltype_merge_umap.pdf"))
print(DimPlot(Singlet_Coassay_Data, group.by = c("celltype_merge"),reduction = "umap",label = T))
dev.off()
celltype_merge <- Singlet_Coassay_Data@meta.data[["celltype_merge"]]
cat("the proportion of celltype:\n")
table(celltype_merge)/sum(table(celltype_merge))
cat("save seurat data...\n")
##修改点
saveRDS(Singlet_Coassay_Data,paste0(save_dir,"Singlet_Coassay_BGI_test_PBMC2.rds"))


cat("Step5:split data of CiteFuse...\n")


CiteFuse_Data <- NULL
CiteFuse_Data$RNA <- Data_mat[["Gene Expression"]]
ADT_HTO <- Data_mat$`Antibody Capture`
CiteFuse_Data$HTO <- ADT_HTO[1:sample_num,]
CiteFuse_Data$ADT <- ADT_HTO[(sample_num+1):143,]

Sce_CiteFuse_Data <- preprocessing(CiteFuse_Data)

Sce_CiteFuse_Data <- normaliseExprs(sce = Sce_CiteFuse_Data,
                                      altExp_name = "HTO",
                                      transform = "log")
cat("CiteFuse split crossSampleDoublets:\n")
Sce_CiteFuse_Data <- crossSampleDoublets(Sce_CiteFuse_Data)

table(Sce_CiteFuse_Data$doubletClassify_between_label)
table(Sce_CiteFuse_Data$doubletClassify_between_class)
data_class <- table(Sce_CiteFuse_Data$doubletClassify_between_class)
siglets_rate <- as.double(100*data_class[3]/sum(data_class))
doublets_rate <- as.double(100*data_class[1]/sum(data_class))
negtive_rate <- as.double(100*data_class[2]/sum(data_class))

print(paste0("siglets_rate:",siglets_rate))
print(paste0("doublets_rate:",doublets_rate))
print(paste0("negtive_rate:",negtive_rate))

cat("CiteFuse split withinSampleDoublets:\n")
Sce_CiteFuse_Data <- withinSampleDoublets(Sce_CiteFuse_Data,minPts = 10)
table(Sce_CiteFuse_Data$doubletClassify_within_label)
table(Sce_CiteFuse_Data$doubletClassify_within_class)

cat("save sce data...\n")
##修改点
saveRDS(Sce_CiteFuse_Data,paste0(save_dir,"Sce_CiteFuse_BGI_test_PBMC2.rds"))

cat("Step6:CiteFuse result to cluster and annotation\n")

Coassay_Data@meta.data[["doubletClassify_between_class"]] <- Sce_CiteFuse_Data$doubletClassify_between_class
Coassay_Data@meta.data[["doubletClassify_within_class"]] <- Sce_CiteFuse_Data$doubletClassify_within_class

Singlet_CiteFuse_Data<- Coassay_Data[,which(Coassay_Data@meta.data[["doubletClassify_between_class"]] == "Singlet" & Coassay_Data@meta.data[["doubletClassify_within_class"]] == "Singlet")]

Singlet_CiteFuse_Data@active.assay <- "RNA"

# Select the top 1000 most variable features
Singlet_CiteFuse_Data <- FindVariableFeatures(Singlet_CiteFuse_Data)

# Scaling RNA data, we only scale the variable features here for efficiency
Singlet_CiteFuse_Data <- ScaleData(Singlet_CiteFuse_Data, features = VariableFeatures(Singlet_CiteFuse_Data))

# Run PCA
Singlet_CiteFuse_Data <- RunPCA(Singlet_CiteFuse_Data,features = VariableFeatures(Singlet_CiteFuse_Data))

# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
Singlet_CiteFuse_Data <- FindNeighbors(Singlet_CiteFuse_Data, reduction = "pca", dims = 1:20)
Singlet_CiteFuse_Data <- FindClusters(Singlet_CiteFuse_Data, resolution = 0.2, verbose = FALSE)
Singlet_CiteFuse_Data <- RunUMAP(Singlet_CiteFuse_Data,dims = 1:20)
##修改点
pdf(paste0(save_dir,"CiteFuse_Dimplot_umap.pdf"))
print(DimPlot(Singlet_CiteFuse_Data,reduction = "umap"))
dev.off()

cat("SingleR annotation:\n")

load("/home/jiangxinyu/R/R_script/biomamba_single-sample/SingR annotation lab/BlueprintEncode_bpe.se_human.RData")
load("/home/jiangxinyu/R/R_script/biomamba_single-sample/SingR annotation lab/HumanPrimaryCellAtlas_hpca.se_human.RData")
Data_for_SingleR <- GetAssayData(Singlet_CiteFuse_Data, slot="data")  ###提取RNA的转录表达数据
clusters <- Singlet_CiteFuse_Data@meta.data$seurat_clusters
###把scRNA数据中的seurat_clusters提取出来，注意是因子类型
cellpred <- SingleR(test = Data_for_SingleR, ref = list(bpe.se,hpca.se), labels = list(bpe.se$label.main, hpca.se$label.main))
Singlet_CiteFuse_Data@meta.data$celltype <- cellpred$labels
##修改点
pdf(paste0(save_dir,"CiteFuse_celltype_umap.pdf"))
print(DimPlot(Singlet_CiteFuse_Data, group.by = c("celltype"),reduction = "umap",label = T))
dev.off()
##合并同类细胞类型
cat("Step4:merge same celltype and calculate the proportion of celltype...\n")
celltype_Data <- Singlet_CiteFuse_Data@meta.data[["celltype"]]
celltype_Data[celltype_Data=="B_cell"] <- "B-cells"
celltype_Data[celltype_Data=="NK_cell"] <- "NK cells"
celltype_Data[celltype_Data=="Monocyte"] <- "Monocytes"
cat("merge celltype number:\n")
table(celltype_Data)
Singlet_CiteFuse_Data@meta.data[["celltype_merge"]] <- celltype_Data

pdf(paste0(save_dir,"CiteFuse_celltype_merge_umap.pdf"))
print(DimPlot(Singlet_CiteFuse_Data, group.by = c("celltype_merge"),reduction = "umap",label = T))
dev.off()

celltype_merge <- Singlet_CiteFuse_Data@meta.data[["celltype_merge"]]
cat("the proportion of celltype:\n")
table(celltype_merge)/sum(table(celltype_merge))
cat("save seurat data...\n")
##修改点
saveRDS(Singlet_CiteFuse_Data,paste0(save_dir,"Singlet_CiteFuse_BGI_test_PBMC2.rds"))


time2 <- Sys.time()
print(paste0("using times:",time2-time1))
