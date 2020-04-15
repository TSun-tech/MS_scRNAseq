#attach analysis packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

#read in two datasets seperately
#gcb: all inflamattory cells subset from GSE118257
#mp:all cells from GSE124335
gcb<-readRDS("/media/tsun/Data/Tsun/small_seminar_room_backup/data/MPI_HS_GCB_DATA/GCB_HS_snRNA_progressed_inflammatory cells.rds")
mp<-readRDS("/media/tsun/Data/Tsun/external_data/Marco_Prinz/Marco_Prinz_2019_HS_microglia_scRNAseq/MP_HS_snRNA_processed.rds")

#check objects
gcb
mp

#remove RM Lesion from gcb dataset sue to low cell number
Idents(gcb)<-"Lesion"
gcb<-subset(gcb, idents = "RM", invert= T)

#create new meta.data tags to prepare data merging
gcb@meta.data$group_for_merge<-gcb@meta.data$new_condition
gcb@meta.data$data<-"gcb"
mp@meta.data$data<-"mp"
mp@meta.data$group_for_merge<-mp@meta.data$Condition

#organise group_for_merge based on lesion type "Control/NAWM" or "MS"
for (i in 1:nrow(mp@meta.data)){
    mp@meta.data$group_for_merge[[i]]<-strsplit(mp@meta.data$group_for_merge[[i]], "_")[[1]][1]
}
for (i in 1:nrow(mp@meta.data)){
     mp@meta.data$group_for_merge[[i]]<-gsub("Healthy", "Control/NAWM", mp@meta.data$group_for_merge[[i]])
}
for (k in 1:nrow(gcb@meta.data)){
     gcb@meta.data$group_for_merge[[k]]<-gsub("Lesion", "MS", gcb@meta.data$group_for_merge[[k]])
}

#check if meta.data tags are correct
unique(gcb@meta.data$group_for_merge)
unique(mp@meta.data$group_for_merge)

#merge two dataswets
seu<-merge(gcb, y = mp, add.cell.ids = c("gcb","mp"), project = "Berghoff")

#check merged object and clean meta.data
seu
length(colnames(seu@meta.data))
colnames(seu@meta.data)

seu@meta.data<-seu@meta.data[,-c(1,2, 8:10,15:18,21,22,23:25,27:31)]

#Seurat scRNA-seq analysis pipeline
#use group CCA for data interptation

#data filter
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT.")

Idents(seu)<-"group_for_merge"
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "Condition")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "Condition")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15)

#CCA
seu.list<-SplitObject(seu, split.by = "data")

for (i in 1:length(seu.list)) {
    seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
    seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}

seu.anchors <- FindIntegrationAnchors(object.list = seu.list, dims = 1:30)

seu<-IntegrateData(anchorset = seu.anchors, dims = 1:30)

seu<-ScaleData(seu)
seu<-RunPCA(seu, npcs = 50)

ElbowPlot(seu)

#select the first 20 PC for emedding
seu<-RunUMAP(object = seu, reduction = "pca", dims = 1:20)

#check merged data
DimPlot(seu, reduction = "umap", group.by = "group_for_merge")
DimPlot(seu, reduction = "umap", group.by = "Condition")
DimPlot(seu, reduction = "umap", group.by = "data")

#Use RNA assay to plot gene expression levels on two dimensional plot
DefaultAssay(seu)<-"RNA"
genes<-c("DHCR24","PLP1","APOE","ABCA1","NR1H3","NR1H2","TREM2", "CX3CR1", "ITGAM", "P2RY12", "CD68", "AIF1", "LAMP2",
        "CD3G","CD3D","CD3E")
for (k in genes){
print(FeaturePlot(seu, features = k))
    }

#use intergrated assay for cluster detection
DefaultAssay(seu)<-"integrated"

seu<-FindNeighbors(seu, dims = 1:20)

#use relatively low resolution (0.1)
seu<-FindClusters(seu, resolution =0.1)

#visualize cluster distribution
DimPlot(seu, reduction = "umap", label = T)

#compute cell proportion in each cluster
sample_per_cluster<-as.data.frame(table(seu@meta.data$group_for_merge,seu@meta.data$integrated_snn_res.0.1))
sample_proportion<-table(seu@meta.data$integrated_snn_res.0.1,seu@meta.data$group_for_merge)
colnames(sample_per_cluster)<-c("Group", "Cluster", "Freq")

#visualize with barplot
ggplot(sample_per_cluster, aes(y= Freq, x= Cluster, 
                                  fill= Group)) +  geom_bar( stat="identity", position="fill")+
theme(axis.text.x = element_text(angle = 90, 
                                 hjust = 1, vjust = 0.5))+
ylab("Contribution of microglia cells from healthy or MS patient in cluster")+
scale_x_discrete(labels= rownames(sample_proportion))+ coord_cartesian(expand = FALSE)+
scale_fill_manual(values=c("lightsteelblue4", "indianred2"))



#example script of calculating certain genes avarage expression, proportion and positive rate presented in each cluster
#example gene: DHCR24

position<-match("DHCR24", rownames(seu))
position

DefaultAssay(seu)<-"RNA"
seu@meta.data$DHCR24_exp<-"NA"

seu@meta.data$DHCR24_exp[which(GetAssayData(object = seu, slot = "counts")[325,]>0)]<-"Pos"
seu@meta.data$DHCR24_exp[which(GetAssayData(object = seu, slot = "counts")[325,]==0)]<-"Neg"

table(seu@meta.data$group_for_merge,seu@meta.data$DHCR24_exp)

seu@meta.data$DHCR24_exp_condition<-"NA"
for (i in 1:nrow(seu@meta.data)){
    seu@meta.data$DHCR24_exp_condition[[i]]<-paste0(seu@meta.data$group_for_merge[[i]],
                                                    "_DHCR24_",
                                                   seu@meta.data$DHCR24_exp[[i]])
}
seu@meta.data$cluster_condition<-"NA"
for (i in 1:nrow(seu@meta.data)){
    seu@meta.data$cluster_condition[[i]]<-paste0(seu@meta.data$seurat_clusters[[i]],
                                                    "_",
                                                   seu@meta.data$group_for_merge[[i]])
}

unique(seu@meta.data$DHCR24_exp_condition)
unique(seu@meta.data$cluster_condition)

sample_per_cluster<-as.data.frame(table(seu@meta.data$DHCR24_exp,seu@meta.data$cluster_condition))
sample_proportion<-table(seu@meta.data$cluster_condition,seu@meta.data$group_for_merge)

colnames(sample_per_cluster)<-c("Group", "Cluster", "Freq")

ggplot(sample_per_cluster, aes(y= Freq, x= Cluster, 
                                  fill= Group)) +  geom_bar( stat="identity", position="fill")+
theme(axis.text.x = element_text(angle = 90, 
                                 hjust = 1, vjust = 0.5))+
ylab("Percentage of DHCR24 positive and negative cells")+
scale_x_discrete(labels= rownames(sample_proportion))+ coord_cartesian(expand = FALSE)+
scale_fill_manual(values=c("grey90","indianred2"))

#calculate average expression from positive cells
Idents(seu)<-"DHCR24_exp"
seu_dhcr24_pos<-subset(seu, idents = "Pos")

Idents(seu_dhcr24_pos)<-"cluster_condition"
dhcr<-AverageExpression(seu_dhcr24_pos)


