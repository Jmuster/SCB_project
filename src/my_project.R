#Setting environment
library(Seurat)
library(dplyr)

#Setting global variables
working_directory='F:/Documents/WORK/EPFL/Master/2019-2020/Spring semester/Single cell biology/Project'
data_path='./data/MOUSE_BRAIN_DATASET_3_COUNTS.tsv'
project_name='Single cell biology'
min_cells=3
min_features=200
max_features=2500
max_percent_MT=3
top_var=2000

setwd(working_directory)

#Importing data
raw_data=read.table(data_path,sep='\t',header = T,as.is = T, row.names=1)
crm=CreateSeuratObject(raw_data, project= project_name, min.cells = min_cells, min.features = min_features)

#Calculating perc mitochondrial RNA
crm[['percents.MT']]=PercentageFeatureSet(crm, pattern="^mt-")

#plotting features prior to selection
jpeg("./data/results/pre_selection_plot.jpeg")
print(VlnPlot(crm, features = c('nCount_RNA','nFeature_RNA','percents.MT'), ncol = 3))
dev.off()

#Filtering Data based on UMIs and percentage MT RNA
crm=subset(crm, subset = nFeature_RNA>min_features & nFeature_RNA<max_features & percents.MT<max_percent_MT)

#Normalizing the data
crm=NormalizeData(crm)

#Selecting only "top_var" most variable features
crm=FindVariableFeatures(crm, selection.method = "vst", nfeatures = top_var)

#plotting feature variance
jpeg("./data/results/features_variance.jpeg")
print(VariableFeaturePlot(crm))
dev.off()

#Data scaling and pca
all.genes= rownames(crm)
crm <- ScaleData(crm, features = all.genes)
crm<-RunPCA(crm, features = VariableFeatures(object = crm))

#pca vizualisation
print(crm[["pca"]], dims = 1:5 , nfeatures =5)
VizDimLoadings(crm, dims= 1:2, reduction = "pca")
DimPlot(crm, reduction ="pca")

#dimensionality
crm = JackStraw(crm, num.replicate = 100)
crm = ScoreJackStraw(crm, dims = 1:15)

#Visualtisation
JackStrawPlot(crm, dims = 1:15)

#Elbow
ElbowPlot(crm)
# Most of the variance of the dataset is included in the 10 first PCs, EVTL 15

#Clustering
crm = FindNeighbors(crm, dims = 1:10)
crm = FindClusters(crm, resolution = 0.5)
head(Idents(crm),5)

#UMAP
crm = RunUMAP(crm,dims = 1:10)
DimPlot(crm, reduction = "umap")

saveRDS(crm , file = "./data/results/scb_project.rds")

working_directory='F:/Documents/WORK/EPFL/Master/2019-2020/Spring semester/Single cell biology/Project'
setwd(working_directory)
readRDS("./data/results/scb_project.rds")