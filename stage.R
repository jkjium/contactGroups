####################################################
# pin cg_ver
library(Matrix)
library(tidyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(ape)
library(reticulate)
library(DoubletFinder)
library(glue)
library(readr)
library(doParallel)
library(WGCNA)
library(gridExtra)
library(ROCR)
library(GO.db)
library(goseq)
library(conflicted)
#BiocManager::install("goseq")
#BiocManager::install("GO.db")

setwd('C:\\Users\\kjia\\workspace\\coral\\stage.acropora_raw')

rm(list=ls())
source('coral_ppl.R')

#load('ad3456.clean.sym_split.rd')
load('ad3456.clean.sym_split.filtered.rd') # after remove small clusters
load('wgcna.bwnet.ad3456.split.rd')

obj <- ad3456.clean

# 06.ad2all.sym_split.summary.vec2.tsv # alias for heatmap
# 06.ad2all.sym_split.summary.r.vec2.tsv # alias for dotplot

#####################################################

setwd('C:\\Users\\kjia\\workspace\\foxd\\stage')
rm(list=ls())
source('coral_ppl.R')

###################################################
# foxd vlnplot & barplot ad, at, ch
# kjia@DESKTOP-L0MMU09 ~/workspace/foxd/stage 2024-12-10 17:04:27

std_obj <-readRDS('javier_clytia_converted_to_seurat_byjake_241210.rds')
std_obj$cell_type_family <- std_obj$annos
std_obj$cell_type <- std_obj$annosSub

t <- GetAssayData(std_obj, layer = "counts")

sum(colSums(t) == 0)
sum(rowSums(t) == 0)
zero_counts_per_cell <- colSums(t)  # Zero counts per cell
zero_counts_per_gene <- rowSums(t)  # Zero counts per gene
print(length(zero_counts_per_cell))  # Cells with zero counts
print(length(zero_counts_per_gene))  # Genes with zero counts


std_obj <- std_seurat_ppl(std_obj)
std_obj <- NormalizeData(std_obj , normalization.method = "LogNormalize", scale.factor = 10000)
std_obj <- FindVariableFeatures(std_obj, selection.method = "vst", nfeatures = 3000)
std_obj <- ScaleData(std_obj, verbose = F, vars.to.regress = "nCount_RNA")
std_obj <- RunPCA(std_obj, npcs = 20, verbose = FALSE)
std_obj <- RunTSNE(std_obj, dims = 1:20)


pref='ch'
genes_file <- paste0(pref,'.genes.tsv')
barcodes_file <- paste0(pref,'.barcodes.tsv')
mtx_file <- paste0(pref,'.counts.mtx')
counts <- ReadMtx(mtx_file, cells=barcodes_file, features=genes_file, feature.column = 1)
s_obj <- CreateSeuratObject(counts = counts)

std_obj <- std_seurat_ppl(s_obj)

clustering_data <- read.table(paste0(pref,'.all_assignments.tsv'), header = TRUE, sep = "\t", row.names = 1)
for (col_name in colnames(clustering_data)) {
  std_obj@meta.data[[col_name]] <- clustering_data[[col_name]][match(rownames(std_obj@meta.data), rownames(clustering_data))]
}
std_obj <- RunTSNE(std_obj, dims = 1:20)
std_obj <- RunUMAP(std_obj, dims = 1:20, min.dist = 0.1, umap.method = "umap-learn", metric = "cosine", verbose = FALSE)
save(std_obj, file='01.cn.clustering2.rd')

DimPlot(std_obj, reduction = "tsne", group.by = "annos", pt.size = 2, alpha = 0.7, label=T)
DimPlot(std_obj, reduction = "umap", group.by = "annos", pt.size = 2, alpha = 0.7, label=T)



gene_expression_values <- as.vector(GetAssayData(std_obj, layer = "data"))

num_expressed_genes_raw <- sum(rowSums(GetAssayData(std_obj, layer = "counts") > 3) > 0)
print(num_expressed_genes_raw)

num_expressed_genes_raw <- sum(rowSums(GetAssayData(ad3456.clean, layer = "data") > 3) > 0)
print(num_expressed_genes_raw)

ggplot(data.frame(expression = gene_expression_values), aes(x = expression)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Gene Expression Across All Genes",
       x = "Normalized Expression",
       y = "Frequency")


# ad
load('ad3456.clean.sym_split.filtered.rd')
std_obj <- ad3456.clean
write.table(std_obj$merged_clusters, file = "ad3456.clean.merged_clusters.tsv", sep = "\t", quote = FALSE, col.names = FALSE)
#$ dos2unix ad3456.clean.merged_clusters.tsv
# $ awk -F '\t' 'FNR==NR{a[$1]=$2;next} {if(a[$2]!="") print $1"\t"a[$2]}' 05.ad2all_merged_clusters_summary_0.40.tsv ad3456.clean.merged_clusters.tsv > 01.ad3456.samap_cluster.tsv
# output: 01.ad3456.samap_cluster.tsv (need add {index\tsamap_family} as header)

clustering_data <- read.table('01.ad3456.samap_cluster.tsv', header = TRUE, sep = "\t", row.names = 1)
col_name <- 'samap_family'
std_obj@meta.data[[col_name]] <- clustering_data[[col_name]][match(rownames(std_obj@meta.data), rownames(clustering_data))]
Idents(std_obj) <- std_obj$samap_family
save(std_obj, file='01.ad3456.clean.samap_family.rd')

load('01.ad3456.clean.samap_family.rd')
pref <- "ad" 
gene_set <- c('adig-s0019.g96', 'adig-s0020.g92', 'adig-s0031.g180', 'adig-s0032.g47', 'adig-s0032.g48', 'adig-s0042.g184', 'adig-s0125.g82') # orthofinder
gene_set <- c('adig-s0022.g99', 'adig-s0042.g184') # prost
cname = 'samap_family'
Idents(std_obj) <- std_obj$samap_family

# at
load('at345.clean.sym_split_v2.rd')
std_obj <- at345.clean
# save assignments for matching samap annotations
write.table(std_obj$merged_clusters, file = "at345.clean.merged_clusters.tsv", sep = "\t", quote = FALSE, col.names = FALSE)
#$ dos2unix ad3456.clean.merged_clusters.tsv
# $ awk -F '\t' 'FNR==NR{a[$1]=$2;next} {if(a[$2]!="") print $1"\t"a[$2]}' 05.at2all_merged_clusters_summary_0.40.tsv at345.clean.merged_clusters.tsv > 01.at345.samap_cluster.tsv
# output: 01.at345.samap_cluster.tsv (need add {index\tsamap_family} as header)

clustering_data <- read.table('01.at345.samap_cluster.tsv', header = TRUE, sep = "\t", row.names = 1)
col_name <- 'samap_family'
std_obj@meta.data[[col_name]] <- clustering_data[[col_name]][match(rownames(std_obj@meta.data), rownames(clustering_data))]
Idents(std_obj) <- std_obj$samap_family
save(std_obj, file='01.at345.clean.samap_family.rd')

load('01.at345.clean.samap_family.rd')
pref <- "at" 
gene_set <- c('aten-s0004.g46', 'aten-s0017.g115', 'aten-s0019.g48', 'aten-s0028.g98', 'aten-s0041.g37', 'aten-s0133.g54', 'aten-s0133.g55')
gene_set <- c('aten-s0061.g34') # prost
cname = 'samap_family'
Idents(std_obj) <- std_obj$samap_family

rm(list=ls())
source('coral_ppl.R')

load('xe.seurat.info3.rd')
pref='xe'
gene_set <- c('Xesp-000527', 'Xesp-003050', 'Xesp-006291', 'Xesp-012614', 'Xesp-018653', 'Xesp-022000', 'Xesp-022001', 'Xesp-002944')

load('sp.seurat.info3.rd')
pref <- "sp" 
gene_set <- c('Spis-XP-022780413-1', 'Spis-XP-022781399-1', 'Spis-XP-022781401-1', 'Spis-XP-022789920-1', 'Spis-XP-022799810-1', 'Spis-XP-022800136-1', 'Spis-XP-022800147-1', 'Spis-XP-022810126-1')

load('hy.seurat.info3.rd')
pref <- "hy" 
gene_set <- c('Hvul-g24126-1', 'Hvul-g30219-1') # orthofinder

# nt -----------
load('nt.seurat.info2.rd')
load('alison.Robj')
std_obj <- AllData
#std_obj <- std_seurat_ppl(std_obj)
pref <- "nt" 
gene_set <- c('FoxA', 'FoxB', 'FOXD1-like-1', 'FOXL1-like-1', 'FOXL1-like-2', 'FOXL2', 'FXC1B-like-1') # orthofinder
std_obj$cell_type_family <- std_obj$IDs
std_obj$cell_type <- std_obj$ID.separate
save(std_obj, file='01.nt.info2.rd')


cname = 'cell_type_family'
Idents(std_obj) <- std_obj$cell_type_family

cname = 'cell_type'
Idents(std_obj) <- std_obj$cell_type

cname = 'metacell'
Idents(std_obj) <- std_obj$metacell


#### func ##################################
# foxd vlnPlt
for (gene in gene_set) {
  gg <- VlnPlot(std_obj, features = gene, ncol=1, pt.size = 0, group.by = cname) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + NoLegend()
  outname = paste0(pref,'_vlnplot_orthofinder_',cname,'_',gene,'_','.pdf')
  print(outname)
  ggsave(file=outname, plot=gg, width = 16, height = 10)  
}

# barplt
average_expression <- AverageExpression(std_obj, features = gene_set, layer="counts")
for (gene in gene_set) {
  gene_expression <- average_expression$RNA[gene, ]
  gene_df <- data.frame(Cluster = names(gene_expression),Expression = as.numeric(gene_expression))
  
  gg <- ggplot(gene_df, aes(x = Cluster, y = Expression, fill = 'red'))  +
    geom_bar(stat = "identity") +
    ggtitle(gene) +
    labs(x = "Cluster", y = "Average raw counts") +
    theme_minimal() +
    NoLegend()
  
  gg <- gg + theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1))
  
  outname = paste0(pref,'_barplot_orthofinder_',cname,'_',gene,'_','.pdf')
  print(outname)
  ggsave(file=outname, plot=gg, width = 16, height = 10) # for samap_family
  #ggsave(file=outname, plot=gg, width = 32, height = 8) # metacell
  #ggsave(file=outname, plot=gg, width = 24, height = 8) # nt cell type
  #ggsave(file=outname, plot=gg, width = 16, height = 8)
}

###################################################
load('xe.seurat.info3.rd')
pref='xe'
gene_set <- c('Xesp-000527', 'Xesp-003050', 'Xesp-006291', 'Xesp-012614', 'Xesp-018653', 'Xesp-022000', 'Xesp-022001', 'Xesp-002944')

load('sp.seurat.info3.rd')
pref <- "sp" 
gene_set <- c('Spis-XP-022780413-1', 'Spis-XP-022781399-1', 'Spis-XP-022781401-1', 'Spis-XP-022789920-1', 'Spis-XP-022799810-1', 'Spis-XP-022800136-1', 'Spis-XP-022800147-1', 'Spis-XP-022810126-1')

load('hy.seurat.info3.rd')
pref <- "hy" 
gene_set <- c('Hvul-g24126-1', 'Hvul-g30219-1') # orthofinder


# nt -----------
load('nt.seurat.info2.rd')
load('alison.Robj')
std_obj <- AllData
std_obj <- std_seurat_ppl(std_obj)
pref <- "nt" 
gene_set <- c('FoxA', 'FoxB', 'FOXD1-like-1', 'FOXL1-like-1', 'FOXL1-like-2', 'FOXL2', 'FXC1B-like-1') # orthofinder
std_obj$cell_type_family <- std_obj$IDs
std_obj$cell_type <- std_obj$ID.separate



#####################################################
setwd('C:\\Users\\kjia\\workspace\\foxd\\stage')
t<-readRDS('Hydra_Seurat_Whole_Transcriptome.rds')
t1<-UpdateSeuratObject(t)
DimPlot(t1, reduction = "tsne", group.by = "cluster.short", pt.size = 2, alpha = 0.7)
pref <- "nt" 
gene_set <- c('FoxA', 'FoxB', 'FOXD1-like-1', 'FOXL1-like-1', 'FOXL1-like-2', 'FOXL2', 'FXC1B-like-1') # orthofinder


###
#std_obj <- RunPCA(std_obj, npcs = 20, verbose = FALSE)
std_obj <- RunTSNE(std_obj, dims = 1:20)
#save(std_obj, file=paste0(pref,'.seurat.info2.umap.rd'))


cname = 'cell_type_family'
Idents(std_obj) <- std_obj$cell_type_family

cname = 'cell_type'
Idents(std_obj) <- std_obj$cell_type

cname = 'metacell'
Idents(std_obj) <- std_obj$metacell

#cl_df <- std_obj[[cname]]
#Idents(std_obj) <- setNames(cl_df$cell_type_family, rownames(cl_df))

# Umap
gg <- DimPlot(std_obj, reduction = "tsne", alpha = 0.7, label=T, group.by=cname, pt.size = 1, label.size = 4) 
#gg <- DimPlot(std_obj, reduction = "umap", alpha = 0.7, label=T, group.by=cname, pt.size = 1, label.size = 3)+ NoLegend()
gg
outname = paste0(pref,'_FOX_dimplot_orthofinder_',cname,'.pdf')
outname
ggsave(file=outname, plot=gg, width = 20, height = 12)


# violinPlot
gg <- VlnPlot(std_obj, features = gene_set, ncol=1, pt.size = 0, group.by = cname)
gg <- gg+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
#gg
outname = paste0(pref,'_FOX_vlnplot_orthofinder_',cname,'.pdf')
outname
ggsave(file=outname, plot=gg, width = 40, height = 8)#24)


# FeaturePlot
#gg <- FeaturePlot(std_obj,features=gene_set, alpha = 0.5, pt.size = 2, blend = F,reduction = "umap",cols =c("lightgrey","red"), label=T, label.size = 3) # for nt
gg <- FeaturePlot(std_obj,features=gene_set, alpha = 0.5, pt.size = 2, blend = F,reduction = "tsne",cols =c("lightgrey","red"), label=T, label.size = 3)
outname = paste0(pref,'_FOX_featureplot_orthofinder_',cname,'.pdf')
outname
ggsave(file=outname, plot=gg, width = 24, height = 12)#24)
gg

### alternative 
coor_df <- as.data.frame(std_obj[['tsne']]@cell.embeddings)
coor_df$cluster <- std_obj$cell_type_family
centroids <- coor_df %>%  group_by(cluster) %>% summarize(tSNE_1 = mean(tSNE_1), tSNE_2 = mean(tSNE_2))
gg + geom_text(data = centroids, aes(x=tSNE_1, y=tSNE_2, label = cluster), 
              check_overlap = TRUE, 
              size = 3, 
              color = "black")

t<-readRDS('Hydra_Seurat_Whole_Transcriptome.rds')
t1<-UpdateSeuratObject(t)

#rm(t)
Hvul_sc_UMI_counts.RDS


######################################################
# foxd init seurat obj for sp,xe,hy 1205
# hy.barcodes.tsv
# hy.clustering2.tsv
# hy.counts.mtx
# hy.genes.tsv
# 
# sp.barcodes.tsv
# sp.clustering2.tsv
# sp.counts.mtx
# sp.genes.tsv
# 
# xe.barcodes.tsv
# xe.clustering2.tsv
# xe.counts.mtx
# xe.genes.tsv

setwd('C:\\Users\\kjia\\workspace\\samap_coral\\stage3.sandbox')
source('coral_ppl.R')


genes_file <- paste0(pref,'.genes.tsv')
barcodes_file <- paste0(pref,'.barcodes.tsv')
mtx_file <- paste0(pref,'.counts.mtx')

counts <- ReadMtx(mtx_file, cells=barcodes_file, features=genes_file, feature.column = 1)
s_obj <- CreateSeuratObject(counts = counts)

std_obj <- std_seurat_ppl(s_obj)

# append metacell info
#rm(list=ls())
#source('coral_ppl.R')
#load('hy.seurat.info2.tsne.rd')
#pref='hy'
#load('sp.seurat.info2.tsne.rd')
#pref='sp'
#load('xe.seurat.info2.tsne.rd')
#pref='xe'

clustering_data <- read.table(paste0(pref,'.clustering3.tsv'), header = TRUE, sep = "\t", row.names = 1)
for (col_name in colnames(clustering_data)) {
  std_obj@meta.data[[col_name]] <- clustering_data[[col_name]][match(rownames(std_obj@meta.data), rownames(clustering_data))]
}

save(std_obj, file=paste0(pref,'.seurat.info3.rd'))


# ad
std_obj <- ad3456.clean
pref <- "ad" 
gene_set <- c('adig-s0019.g96', 'adig-s0020.g92', 'adig-s0031.g180', 'adig-s0032.g47', 'adig-s0032.g48', 'adig-s0042.g184', 'adig-s0125.g82') # orthofinder
gene_set <- c('adig-s0022.g99', 'adig-s0042.g184') # prost

# at
load('at345_nosym.mc.RData')
std_obj <- at345_nosym.mc
pref <- "at" 
gene_set <- c('aten-s0004.g46', 'aten-s0017.g115', 'aten-s0019.g48', 'aten-s0028.g98', 'aten-s0041.g37', 'aten-s0133.g54', 'aten-s0133.g55')
gene_set <- c('aten-s0061.g34')

# --------
source("coral_ppl.R")
# hy
load('hy.seurat.info2.rd')
#hy_obj <- std_obj
std_obj <- hy_obj
pref <- "hy" 
gene_set <- c('Hvul-g24126-1', 'Hvul-g30219-1') # orthofinder
gene_set <- c('Hvul-g5578-1', 'Hvul-g15144-1') # prost
gene_set <- c('Hvul-g5578-1', 'Hvul-g15144-1','Hvul-g24126-1', 'Hvul-g30219-1') 

proc_foxgene(std_obj, pref, 'orthofinder')


# nt
load('nt.seurat.info2.rd')
pref <- "nt" 
gene_set <- c('FoxA', 'FoxB', 'FOXD1-like-1', 'FOXL1-like-1', 'FOXL1-like-2', 'FOXL2', 'FXC1B-like-1') # orthofinder
gene_set <- c('FOXD3-like-1', 'FOXD1-like-1', 'FOXO3-like-1', 'FoxO2-like	', 'FoxO1-like', 'FoxQ2a') #prost
proc_foxgene(std_obj, pref, 'orthofinder')


# sp
load('sp.seurat.info2.rd')
pref <- "sp" 
gene_set <- c('Spis-XP-022780413-1', 'Spis-XP-022781399-1', 'Spis-XP-022781401-1', 'Spis-XP-022789920-1', 'Spis-XP-022799810-1', 'Spis-XP-022800136-1', 'Spis-XP-022800147-1', 'Spis-XP-022810126-1')
gene_set <- c('Spis-XP-022792313-1', 'Spis-XP-022793008-1', 'Spis4284-1') # prost
proc_foxgene(std_obj, pref, 'orthofinder')

# xe
load('xe.seurat.info2.rd')
pref <- "xe" 
gene_set <- c('Xesp-000527', 'Xesp-003050', 'Xesp-006291', 'Xesp-012614', 'Xesp-018653', 'Xesp-022000', 'Xesp-022001', 'Xesp-002944')
proc_foxgene(std_obj, pref, 'orthofinder')



## testing
# violin
expression_data <- FetchData(std_obj, vars = gene_set)
total_expression <- rowSums(expression_data)
std_obj$fox_expr <- total_expression
gg <- VlnPlot(std_obj, features = c("fox_expr"), ncol = 1, pt.size = 0, group.by = "merged_clusters")
gg <- gg + ggtitle(paste0('FOX gene expression (', pref,')')) + theme(legend.position = "none")
gg
outname = paste0(pref,'_FOX_violin_orthofinder.pdf')
#outname = paste0(pref,'_FOX_violin_prost.pdf')
outname
ggsave(file=outname, plot=gg, width = 15, height = 8)


# bar plots
# v1 <- VlnPlot(ad3456, features = c("nFeature_RNA"), ncol = 1, pt.size = 0, y.max = 3000, group.by = "sample")
#VlnPlot(std_obj, features = gene_set, ncol = 1, pt.size = 0, group.by = "cell_type_family")

expression_data <- FetchData(std_obj, vars = gene_set)
cluster_labels <- std_obj$merged_clusters
#cluster_labels <- std_obj$cell_type_family
total_expression_by_cluster <- aggregate(expression_data, by = list(cluster_labels), FUN = sum)
colnames(total_expression_by_cluster) <- c("Cluster", gene_set)
#total_expression_long <- reshape2::melt(total_expression_by_cluster, id.vars = "Cluster", variable.name = "Gene", value.name = "TotalExpression")
ndf <- data.frame(Cluster=total_expression_by_cluster$Cluster, TotalExpression=rowSums(total_expression_by_cluster[, -1]))

Fox_family<-paste0(length(gene_set), ' genes')

gplot<-ggplot(ndf, aes(x = Cluster, y = TotalExpression, fill=Fox_family)) +
#ggplot(total_expression_long, aes(x = Cluster, y = TotalExpression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Total Expression of Selected Genes in Different Clusters",
       x = "Cluster",
       y = "Total Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14))
gplot
outname = paste0(pref,'_FOX_total_expression_orthofinder.pdf')
#outname = paste0(pref,'_FOX_total_expression_prost.pdf')
outname
ggsave(file=outname, plot=gplot, width = 15, height = 8)

######################################################
# acropora cell type tree 20241218
# 
library(tidyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(ape)
library(reticulate)
library(DoubletFinder)
library(glue)
library(readr)

setwd('C:\\Users\\kjia\\workspace\\coral\\stage.acropora_raw')
rm(list=ls())
source('coral_ppl.R')

# load seurat objects {ad,at,sp,xe,nt}
objs <- list()
load('01.ad3456.clean.samap_family.rd')
objs[['Acropora.digitifera']] <- std_obj
load('01.at345.clean.samap_family.rd')
objs[['Acropora.tenuis']] <- std_obj
load('01.nt.seurat.info2.umap.rd')
objs[['Nematostella.vectensis']] <- std_obj
load('01.sp.seurat.info3.rd')
objs[['Stylophora.pistillata']] <- std_obj
load('01.xe.seurat.info3.rd')
objs[['Xenia.sp.']] <- std_obj

# integrate sub objects
install.packages("future")
options(future.globals.maxSize = 10 * 1024^3)
sct_obj <- lapply(objs, function(x){ SCTransform(x, return.only.var.genes = FALSE)})


for(n in colnames(ortho_genes)) {write(x = rownames(sct_obj[[n]]), file = paste0("011.sct.",n,".genes.txt"))}
#kjia@DESKTOP-L0MMU09 ~/workspace/coral/stage.acropora_raw 2024-12-18 16:58:04
#$ ls 011* > 012.seurat.genes.stub
#$ python proc_coral_samap.py filter_1to1 012.seurat.genes.stub 02.ortho4.all.tsv 025.ortho.filtered.tsv
#2024-12-18 16:58:59|2732|0|INFO|ignore absent libraries
#2024-12-18 16:58:59|2732|0|INFO|load 134672 genes
#2024-12-18 16:59:14|2732|15|INFO|save 4538 filtered 1to1 orthologs to 025.ortho.filtered.tsv
# add header to 025.ortho.filtered.tsv

ortho_genes <- read.delim('025.ortho.filtered.tsv', header = TRUE, sep = "\t")

# subset seurat by 1to1 orthologs
sub_objs <- list()
for(n in colnames(ortho_genes)) {
  sub_objs[[n]] <- subset(sct_obj[[n]]$SCT, features = ortho_genes[[n]])
  sub_objs[[n]]$organism_name <- n
}
# check # of genes in each obj
for(n in colnames(ortho_genes)){print(n);print(dim(sub_objs[[n]]))}

# append organism name to clustering assignments
sub_objs[['Acropora.digitifera']]$samap_family <- paste0('Acropora.digitifera - ', sub_objs[['Acropora.digitifera']]$samap_family)
sub_objs[['Acropora.tenuis']]$samap_family <- paste0('Acropora.tenuis - ', sub_objs[['Acropora.tenuis']]$samap_family)
sub_objs[['Nematostella.vectensis']]$cell_type_family <- paste0('Nematostella.vectensis - ', sub_objs[['Nematostella.vectensis']]$cell_type_family)
sub_objs[['Stylophora.pistillata']]$cell_type_family <- paste0('Stylophora.pistillata - ', sub_objs[['Stylophora.pistillata']]$cell_type_family)
sub_objs[['Xenia.sp.']]$cell_type_family <- paste0('Xenia.sp. - ', sub_objs[['Xenia.sp.']]$cell_type_family)




# calculate average expression
avg_expr <- list()
avg_expr[['Acropora.digitifera']] <- AverageExpression(sub_objs[['Acropora.digitifera']], group.by = "samap_family")
avg_expr[['Acropora.tenuis']] <- AverageExpression(sub_objs[['Acropora.tenuis']], group.by = "samap_family")
avg_expr[['Nematostella.vectensis']] <- AverageExpression(sub_objs[['Nematostella.vectensis']], group.by = "cell_type_family")
avg_expr[['Stylophora.pistillata']] <- AverageExpression(sub_objs[['Stylophora.pistillata']], group.by = "cell_type_family")
avg_expr[['Xenia.sp.']] <- AverageExpression(sub_objs[['Xenia.sp.']], group.by = "cell_type_family")

# append organism name 
for(n in colnames(ortho_genes)) { 
  colnames(avg_expr[[n]]$RNA) <- paste0(n,'|',colnames(avg_expr[[n]]$RNA)) 
  avg_expr[[n]]$RNA <- avg_expr[[n]]$RNA[ortho_genes[[n]],]
}

# combine matrices and update gene names
cmat <- avg_expr[[1]]$RNA
crname <- rownames(avg_expr[[1]]$RNA)
for (i in 2:length(avg_expr)) {
  cmat <- cbind(cmat, avg_expr[[i]]$RNA)
  crname <- paste(crname, rownames(avg_expr[[i]]$RNA), sep = ",")
}
rownames(cmat) <- crname
save(cmat, file='06.combnied_avg_mat.5.rd')

dmat <- dist(t(cmat), method = "euclidean")
ct_tree <- nj(dmat)
pdf("07.cell_type_tree.pdf", width = 20, height = 36) 
plot(ct_tree, main = "Phylogenetic Tree of Cell types")
dev.off()
#### debug
avg <- AverageExpression(std_obj, group.by = "samap_family")





######################################################
# acropora distribution of symbiodinium genes
#   transfer merged clusters to full count table
#   get full rdata of ad
#   subset cells by referencing no_sym data
#   umap with color diff by symbiodinium
library(tidyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(ape)
library(reticulate)
library(DoubletFinder)
library(glue)
library(readr)

setwd('C:\\Users\\kjia\\workspace\\coral\\stage.acropora_raw')
source('coral_ppl.R')
############################### pre-production analysis ###############################
# preprocess 10X data
# 1. running standard seurat pipeline (normalization, variable genes, PCA, neighbor, clustering, umap)
# 2. find doublets
# 3. pre-QC for each dataset
ad323 <- fn_prep("./ad323", "ad323", 0.064) 	#   8202 ./ad323/barcodes.tsv
ad423 <- fn_prep("./ad423", "ad423", 0.080) 	#   9727 ./ad423/barcodes.tsv
ad523 <- fn_prep("./ad523", "ad523", 0.080) 	#  11133 ./ad523/barcodes.tsv
ad623 <- fn_prep("./ad623", "ad623", 0.080) 	#   9967 ./ad623/barcodes.tsv

# remove doublets
# .nd: no doublets
ad3.nd <- subset(ad323, subset=(DF.classifications_0.25_0.01_505=="Singlet"))
ad4.nd <- subset(ad423, subset=(DF.classifications_0.25_0.29_748=="Singlet"))
ad5.nd <- subset(ad523, subset=(DF.classifications_0.25_0.01_854=="Singlet"))
ad6.nd <- subset(ad623, subset=(DF.classifications_0.25_0.03_761=="Singlet"))

# remove ribosomal genes
## run ppl: > Identify ribosomal genes using blast
## remove ".t1", change "_" to "-" before loading
ad.ribo.genes <- readLines('02.ad.ribo.gene.txt')
ad3.nd.nr <- subset(ad3.nd, features = setdiff(rownames(ad3.nd) , ad.ribo.genes))
ad4.nd.nr <- subset(ad4.nd, features = setdiff(rownames(ad4.nd) , ad.ribo.genes))
ad5.nd.nr <- subset(ad5.nd, features = setdiff(rownames(ad5.nd) , ad.ribo.genes))
ad6.nd.nr <- subset(ad6.nd, features = setdiff(rownames(ad6.nd) , ad.ribo.genes))

ad3c <- std_seurat_ppl(ad3.nd.nr)
ad4c <- std_seurat_ppl(ad4.nd.nr)
ad5c <- std_seurat_ppl(ad5.nd.nr)
ad6c <- std_seurat_ppl(ad6.nd.nr)

# merge dataset and re-run seurat std ppl
ad3456 <- merge(ad3c, y = c(ad4c, ad5c, ad6c), add.cell.ids = c("ad3", "ad4", "ad5", "ad6"), project = "ad3456c")
ad3456[["RNA"]]<-JoinLayers(ad3456[["RNA"]])

ad3456_new <- CreateSeuratObject(counts = ad3456@assays$RNA$counts)
ad3456.std <- std_seurat_ppl(ad3456_new)
save(ad3456.std, file="ad3456_full_std.rd")


# remove small clusters
load('ad3456.clean.sym_split.rd')
obj <- ad3456.clean
Idents(obj) <- obj$sym_split
t=table(Idents(obj))
clusters_to_keep<-names(t[t>=10])

#t[t<10]
# g22.apo g64.sym g37.sym g24.sym g14.sym 
# 4       4       2       4       3
#clusters_to_remove <- c('g22.apo', 'g64.sym', 'g37.sym', 'g24.sym', 'g14.sym')

#######################################
# sync cluster annotation files
#######################################

subset_obj <- subset(obj, subset = sym_split %in% clusters_to_keep)
ad3456.clean <- subset_obj
save(ad3456.clean, file="ad3456.clean.sym_split.filtered.rd")

# Subset full ad3456.std using the cell IDs from # ad3456.clean.sym_split.filtered.rd
ad3456.clean.full <- ad3456.std[, colnames(ad3456.clean)]
ad3456.clean.full@meta.data <-ad3456.clean@meta.data
save(ad3456.clean.full, file="ad3456.clean.full.sym_split.filtered.rd")

# 


###
# symbiont genes QC
ad3456.std$symbiont_umis <- count_symbiont_umis_per_cell(ad3456.std) 
ad3456.std$symbiont_genes <- count_symbiont_genes_per_cell(ad3456.std)





#######################################################
# WGCNA without small clusters 20241203
# WGCNA without low expressed genes 20241125
# ft WGCNA 20241111
setwd('C:\\Users\\kjia\\workspace\\coral\\stage.acropora_raw')
rm(list=ls())
source('coral_ppl.R')
# load('ad3456.clean.sym_split.rd') # without removing small clusters
load('ad3456.clean.sym_split.filtered.rd')
obj <- ad3456.clean
Idents(obj) <- obj$sym_split

# 1201: use all genes
s=0.8
#w_mat <- wgcna_ppl_1(obj, ngene=as.integer(nrow(obj)*s), expr_filter=20)
w_mat <- wgcna_ppl_1(obj, ngene='all', expr_filter=20) # 7858 / 21038
adj_mat <- adjacency(w_mat, power = 8, type = "signed")
tom_mat <- TOMsimilarity(adj_mat)
dis_tom <- 1-tom_mat
geneTree = hclust(as.dist(dis_tom), method = "average")

# plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dis_tom,
                                 deepSplit = 4, pamRespectsDendro = FALSE,
                                 minClusterSize = 15)
dynamicColors = labels2colors(dynamicMods)
# plotDendroAndColors(geneTree, dynamicMods, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

MEList <- moduleEigengenes(as.matrix(w_mat), colors = dynamicColors)
mes <- MEList$eigengenes

# determine 
METree = hclust(as.dist(1-cor(mes)), method = "average")
plot(METree, main = "Clustering of module eigengenes",  xlab = "", sub = "")
mes_dis_thr = 0.2
abline(h=mes_dis_thr, col = "red")

# old bwnet
merged_list <- mergeCloseModules(as.matrix(w_mat), dynamicColors, cutHeight = mes_dis_thr, verbose = 3)
merged_colors <- merged_list$colors

# overall heatmap
# append cluster alias
# old mes
merged_mes <- merged_list$newMEs
colnames(merged_mes) <- sub("ME", "", colnames(merged_mes))
merged_mes_tr <- t(merged_mes)
col_df <- data.frame(idx = colnames(merged_mes_tr))
#xn_map <- read.csv("05.ad2all_sym_split_summary_0.40_fmt.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("06.ad2all.sym_split.summary.nj.vec2.tsv", header = TRUE, sep = "\t")
cluster_stub <- left_join(col_df, xn_map, by="idx")
colnames(merged_mes_tr) <- cluster_stub$alias
save(merged_mes_tr,file="wgcan.all_genes.heatmap.rd")
pheatmap(seriation_mat(merged_mes_tr), cluster_cols = F, cluster_rows = F, scale = 'none', treeheight_row = 0, treeheight_col = 0)



# wgcna dotplot
## goseq params
geneModuleMembership = as.data.frame(cor(x = w_mat, y = merged_mes, use = "p"))
goseq_params <- list(
  gene.vector = NULL, # named vecotor {aten-0001s.g1: 0/1}
  gene2go = NULL,
  gene.lengths = NULL,
  top_n = 150,
  name = ""
)
# load gene.lengths
gl<-read.csv('ad.gene.length.tsv', sep='\t', header=FALSE)
gene.lengths <- setNames(gl$V2,  gsub("_", "-", gl$V1))

# load gene2go
df <- read.csv('ad.gene2go.tsv', sep='\t', header=FALSE)
df$V1 <- gsub("_", "-", df$V1)
gene2go <- split(df$V2, df$V1)

# assemble params and run goseq
goseq_params$gene.lengths <- gene.lengths
goseq_params$gene2go <- gene2go


# load cluster orders, gene alias, cluster alias
ordered_clusters <- read.table(file = "ad.sym_split_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
#xn_map <- read.csv("05.ad2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")

#$ awk -F '\t' 'FNR==NR{a[$1]=$2;next} {if(a[$1]!="")printf "%s - (%s) %s - %s\n", $1,a[$1],$1,$2}' 05.ad2all_sym_split.cluster_size.tsv 05.ad2all_sym_split_summary_0.40.tsv|awk -F ' - ' 'BEGIN{printf "idx\talias\n"}{printf "%s\t",$1; for(i=3;i<=NF-1;i++) printf "%s - ", $i; printf "%s\n", $2}' > 06.ad2all.sym_split.summary.r.vec2.tsv
xn_map <- read.csv("06.ad2all.sym_split.summary.r.vec2.tsv", header = TRUE, sep = "\t")
source('coral_ppl.R')
wgcna_dotplots_withmaps(obj, merged_list, geneModuleMembership, yn_map, xn_map, 150, ordered_clusters, goseq_params, 'ad3456')


# debug
## percentage of expression
bwnet <- merged_list
me_wt <- geneModuleMembership
module_names <- unique(bwnet$colors)
i = 'black'
genes_to_plot_all <- me_wt[row.names(me_wt)[bwnet$colors == i], i, drop = F]
dim(genes_to_plot_all)
gene_expression <- GetAssayData(obj, layer = "data") 
module_expression <- gene_expression[row.names(genes_to_plot_all), ]
cluster_ids <- Idents(obj)
cluster<-'g4.apo'

cells_in_cluster <- names(cluster_ids[cluster_ids == cluster])
length(cells_in_cluster)
cluster_expr_data <- module_expression[, cells_in_cluster]
percent_in_cluster <- (rowSums(cluster_expr_data > 0) / length(cells_in_cluster)) * 100


## goseq + plot
bwnet <- merged_list
me_wt <- geneModuleMembership
module_names <- unique(bwnet$colors)
i = 'tan'
genes_to_plot_all <- me_wt[row.names(me_wt)[bwnet$colors == i], i, drop = F]
genes_to_plot <- genes_to_plot_all
all.genes <- names(goseq_params$gene.length)
genes_to_plot_ordered <- row.names(genes_to_plot)[order(-genes_to_plot[[i]])]
de.genes <- genes_to_plot_ordered
gene.vector <- as.integer(all.genes %in% de.genes)
names(gene.vector) <- all.genes	
goseq_params$gene.vector <- gene.vector
params <- goseq_params
pwf <- nullp(params$gene.vector, bias.data=params$gene.lengths, plot.fit = FALSE)
res <- goseq(pwf = pwf, gene2cat = params$gene2go, method = "Hypergeometric")
goseq_ret <- res

t <- goseq_ret %>% dplyr::filter(!is.na(term), over_represented_pvalue < 0.05, numDEInCat > 2) %>% arrange(over_represented_pvalue)
#t$term <- factor(t$term, levels = t$term)

df<-head(t, params$top_n) %>% arrange(-over_represented_pvalue)

go_plot<- ggplot(df, aes(x=numDEInCat, y=term, colour=over_represented_pvalue, size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="#.DEInCat", y="GO term", colour="p value", size="Count")+
  scale_y_discrete(limits = unique(df$term))

go_plot <- go_plot + ggtitle(i)
go_plot

ggsave(filename = "output/tt.pdf",
       plot = go_plot, 
       path = ".",
       width = 13, # 18
       height = (length(genes_to_plot_ordered) * (1/6)) + 20, # +2
       #height = 8, 
       units = "in", 
       limitsize = F)
# debug
mes <- t(mes)

merge_6000 = mergeCloseModules(as.matrix(w_mat), dynamicColors, cutHeight = MEDissThres_6000, verbose = 3)

outprefix <- 'ad3456'
#rownames(mes) <- sub("ME", "", rownames(mes))
# load cluster alias 
# awk -F '\t' 'FNR==NR{a[$1]=$1"("$2")";next} {if(a[$1]!=""){printf "%s - %s - %s\n", $1,$2, a[$1]}}' 05.ad2all_sym_split.cluster_size.tsv 05.ad2all_sym_split_summary_0.40_fmt.tsv > t
# awk -F ' - ' 'BEGIN{printf "idx\talias\n"} {printf "%s\t%s", $1,$NF;for(i=3;i<NF;i++){printf " - %s",$i} printf "\n";}' t > 06.ad2all.sym_split.summary.vec2.tsv
# df_cluster_alias <- read.csv("06.ad2all.sym_split.summary.vec2.tsv", header = T, sep = "\t")
# load manually ajusted category file
# category information is for splitting the heatmap
# df_cluster_category <- read.csv("05.ad2all_sym_split.group.tsv", header = F, sep = "\t")
#colnames(df_cluster_category) <- c('cat', 'idx', 'alias')
#wgcna_subgroup_dotplots_withmaps (mes, df_cluster_category, df_cluster_alias, outprefix)

# wgcna dotplot
# load cluster orders, gene alias, cluster alias
ordered_clusters <- read.table(file = "at.sym_split_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, goseq_params, 'at345')




#######################################################
# sub WGCNA split heatmap 
rm(list=ls())
source('coral_ppl.R')
load('ad3456.clean.sym_split.rd')
load('wgcna.bwnet.ad3456.split.rd')

obj <- ad3456.clean
# set activate clustering
Idents(obj) <- obj$sym_split

mes <- t(bwnet$MEs)
rownames(mes) <- sub("ME", "", rownames(mes))

# load cluster alias 
# awk -F '\t' 'FNR==NR{a[$1]=$1"("$2")";next} {if(a[$1]!=""){printf "%s - %s - %s\n", $1,$2, a[$1]}}' 05.ad2all_sym_split.cluster_size.tsv 05.ad2all_sym_split_summary_0.40_fmt.tsv > t
# awk -F ' - ' 'BEGIN{printf "idx\talias\n"} {printf "%s\t%s", $1,$NF;for(i=3;i<NF;i++){printf " - %s",$i} printf "\n";}' t > 06.ad2all.sym_split.summary.vec2.tsv
df_cluster_alias <- read.csv("06.ad2all.sym_split.summary.vec2.tsv", header = T, sep = "\t")

# load manually ajusted category file
# category information is for splitting the heatmap
df_cluster_category <- read.csv("05.ad2all_sym_split.group.tsv", header = F, sep = "\t")
colnames(df_cluster_category) <- c('cat', 'idx', 'alias')

outprefix <- 'ad3456'
wgcna_subgroup_dotplots_withmaps (mes, df_cluster_category, df_cluster_alias, outprefix)



######################################################
# nj neighbor joning with alias
rm(list=ls())
source('coral_ppl.R')
load('ad3456.clean.sym_split.rd')
load('wgcna.bwnet.ad3456.split.rd')
obj <- ad3456.clean
# set activate clustering
Idents(obj) <- obj$merged_clusters

## save cluster size
df_cluster_size <- as.data.frame(table(Idents(obj)))
write.table(df_cluster_size, file = "05.ad3456.merged.cluster_size.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#$ awk -F '\t' 'FNR==NR{a[$1]=$1" - "$1"("$2")";next} {if(a[$1]!=""){printf "%s - %s\n",a[$1], $2} else {print $0}}' 05.ad3456.merged.cluster_size.tsv 05.ad2all_merged_clusters_summary_0.40.tsv |awk -F ' - ' '{if(NR==1) print $0; else {printf "%s\t%s", $1,$2; for(i=NF-1;i>=3;i--){printf " - %s", $i} printf "\n"} }'> 06.ad2all.merged_clusters.summary.vec2.tsv
# output: 06.ad2all.merged_clusters.summary.vec2.tsv

njmet <- as.matrix(read.table("ad3456.wmat.merged_clusters.tsv", sep='\t', header = TRUE))
# append cluster alias
df_njidx <- data.frame(idx=colnames(njmet))
cluster_alias <- read.csv("06.ad2all.merged_clusters.summary.vec2.tsv", header = TRUE, sep = "\t")
cluster_stub <- left_join(df_njidx, cluster_alias, by="idx")
colnames(njmet) <- cluster_stub$alias

#obj_nj <- nj(as.dist(as.matrix(read.table("ad3456.wmat.merged_clusters.tsv", sep='\t', header = TRUE))))
obj_nj <- nj(njmet)
plot(obj_nj, cex = 0.8, no.margin = TRUE)







#######################################################
# WGCNA pheatmap with ordered cluster 
rm(list=ls())
source('coral_ppl.R')
load('ad3456.clean.sym_split.rd')
load('wgcna.bwnet.ad3456.split.rd')
obj <- ad3456.clean
levels(obj)
# set activate clustering
Idents(obj) <- obj$sym_split

mes <- t(bwnet$MEs)
rownames(mes) <- sub("ME", "", rownames(mes))

# exclude small clusters
cluster_counts <- table(Idents(obj))
small_cluster_ids <- names(cluster_counts[cluster_counts < 10])
mes_filered <- mes[, !(colnames(mes) %in% small_cluster_ids)]

# append cluster alias
col_df <- data.frame(idx = colnames(mes_filered))
xn_map <- read.csv("05.ad2all_sym_split_summary_0.40_fmt.tsv", header = TRUE, sep = "\t")
cluster_stub <- left_join(col_df, xn_map, by="idx")
colnames(mes_filered) <- cluster_stub$alias

# seriation and heatmap
out <- seriation_mat(mes_filered)
pheatmap(out, cluster_cols = F, cluster_rows = F, scale = 'none', treeheight_row = 0, treeheight_col = 0)


###############################################
# samap heatmap
raw_table <- read.table('04.heatmap.ad2all.merged_clusters_cell_type_family.blast.tsv', row.names=1, header=TRUE, sep='\t')
dfm <- as.matrix(raw_table)
source('coral_ppl.R')
out <- seriation_mat(dfm)
pheatmap(out, cluster_cols = F, cluster_rows = F, scale = 'none', treeheight_row = 0, treeheight_col = 0)



# x label: nj ordered cluster alias
# $ awk -F '\t' '{print $2}' 05.ad2all_sym_split_summary_0.40.tsv|awk -F ' - ' '{printf "%s",$NF; for(i=NF-1;i>0;i--) printf " - %s", $i; printf "\n"}' |awk '{if(NR==1) print "idx\t"$1; else print $1"\t"$0}' > 05.ad2all_sym_split_summary_0.40_fmt.tsv
xn_map <- read.csv("05.ad2all_sym_split_summary_0.40_fmt.tsv", header = TRUE, sep = "\t")
ordered_clusters <- read.table(file = "ad.sym_split_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
ordered_clusters_df <- tibble(idx = ordered_clusters)
cluster_stub <- left_join(ordered_clusters_df, xn_map, by="idx")

mes_co <- mes[, cluster_stub$idx]
colnames(mes_co) <- cluster_stub$alias

out <- seriation_mat(mes_co)
pheatmap(out, cluster_cols = F, cluster_rows = F, scale = 'none', treeheight_row = 0, treeheight_col = 0)

# misc
## save cluster size
df_cluster_size <- as.data.frame(table(Idents(obj)))
write.table(df_cluster_size, file = "05.ad2all_sym_split.cluster_size.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#### Jake's seriation function #################
order_genes <- function (x, object, cell.type.order, avg.exp.matrix)
{
  x_average_matrix <- avg.exp.matrix$RNA[x,]
  x_average_matrix_max <- colnames(x_average_matrix)[apply(x_average_matrix,1,which.max)]
  names(x_average_matrix_max) <- x
  x_average_matrix_max_ordered <- x_average_matrix_max[order(match(x_average_matrix_max, cell.type.order))]
  return(names(x_average_matrix_max_ordered))
}



########################################################
# 241022 ad comparing two clusters
load('ad3456.clean.sym_split.rd')
obj <- ad3456.clean

## cascade sym_split  ---------------------------------------------------

Idents(obj) <- obj$sym_split
ordered_clusters <- read.table(file = "ad.sym_split_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
obj@active.ident <- factor(Idents(obj), levels = ordered_clusters)

# dotplot params
dotplot_params <- list(
  thr.p_val = 0.01,
  thr.avg_log2FC = 1,
  gene.alias = NULL,
  cluster.alias = NULL
)

yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
dotplot_params$gene.alias<-yn_map
xn_map <- read.csv("05.ad2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
dotplot_params$cluster.alias <-xn_map


## goseq params
goseq_params <- list(
  gene.vector = NULL, # named vecotor {aten-0001s.g1: 0/1}
  gene2go = NULL,
  gene.lengths = NULL,
  top_n = 20
)
# load gene.lengths
gl<-read.csv('ad.gene.length.tsv', sep='\t', header=FALSE)
gene.lengths <- setNames(gl$V2,  gsub("_", "-", gl$V1))

# load gene2go
df <- read.csv('ad.gene2go.tsv', sep='\t', header=FALSE)
df$V1 <- gsub("_", "-", df$V1)
gene2go <- split(df$V2, df$V1)

# assemble params and run goseq
goseq_params$gene.lengths <- gene.lengths
goseq_params$gene2go <- gene2go

# function call
# 26 g12 g12.apo 333 g12.sym 590 0.57 g12 gastrodermis (4) - g12
# 27 g25 g25.apo 233 g25.sym 191 -0.20 g25 gastrodermis (4) - g25
# 28 g15 g15.apo 237 g15.sym 613 0.95 g15 gastrodermis (5) - g15
# 29 g21 g21.apo 209 g21.sym 347 0.51 g21 gastrodermis (5) - g21
# 30 g24 g24.apo 421 g24.sym 4 -4.44 g24 gastrodermis (4) - g24
# 31 g46 g46.apo 58 g46.sym 107 0.60 g46 gastrodermis (4) - g46
# 32 g22 g22.apo 4 g22.sym 498 4.60 g22 alga-hosting_cells - gastrodermis (2) - g22

source('coral_ppl.R')
dotplot_pair_de_goseq(obj, 'g12.sym', 'g46.sym', dotplot_params, goseq_params, 'ad3456')
dotplot_pair_de_goseq(obj, 'g46.sym', 'g12.sym', dotplot_params, goseq_params, 'ad3456')

dotplot_pair_de_goseq(obj, 'g12.apo', 'g46.apo', dotplot_params, goseq_params, 'ad3456')
dotplot_pair_de_goseq(obj, 'g46.apo', 'g12.apo', dotplot_params, goseq_params, 'ad3456')


dotplot_pair_de_goseq(obj, 'g12.sym', 'g22.sym', dotplot_params, goseq_params, 'ad3456')
dotplot_pair_de_goseq(obj, 'g22.sym', 'g12.sym', dotplot_params, goseq_params, 'ad3456')


# generate sym_split clusters dotplot #####################################
source('coral_ppl.R')
generate_cluster_dotplots_withmaps(obj, yn_map, xn_map, ordered_clusters, goseq_params, 'ad3456')


# wgcna split heatmap
# obj<-ad3456.clean
w_mat <- wgcna_ppl_1(obj) # 12 for ad
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 10, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
save(bwnet, file='wgcna.bwnet.ad3456.split.rd')

# pheatmap
mes <- t(bwnet$MEs)
rownames(mes) <- sub("ME", "", rownames(mes))
mes_order_stub<-unique(colnames(mes)[apply(mes,1,which.max)])
pheatmap(mes[, mes_order_stub], cluster_cols = F, cluster_rows = F, scale = 'none', treeheight_row = 0, treeheight_col = 0)

cluster.order<-sort(levels(obj))
avg_max <- colnames(mes)[apply(mes,1,which.max)]
mes_order_stub <- avg_max[order(match(avg_max, cluster.order))]
#############################
order_genes <- function (x, object, cell.type.order, avg.exp.matrix)
{
  x_average_matrix <- avg.exp.matrix$RNA[x,]
  x_average_matrix_max <- colnames(x_average_matrix)[apply(x_average_matrix,1,which.max)]
  names(x_average_matrix_max) <- x
  x_average_matrix_max_ordered <- x_average_matrix_max[order(match(x_average_matrix_max, cell.type.order))]
  return(names(x_average_matrix_max_ordered))
}
##############################

## cascade merged  ---------------------------------------------------
# 26 g12 g12.apo 333 g12.sym 590 0.57 g12 gastrodermis (4) - g12
# 27 g25 g25.apo 233 g25.sym 191 -0.20 g25 gastrodermis (4) - g25
# 28 g15 g15.apo 237 g15.sym 613 0.95 g15 gastrodermis (5) - g15
# 29 g21 g21.apo 209 g21.sym 347 0.51 g21 gastrodermis (5) - g21
# 30 g24 g24.apo 421 g24.sym 4 -4.44 g24 gastrodermis (4) - g24
# 31 g46 g46.apo 58 g46.sym 107 0.60 g46 gastrodermis (4) - g46
# 32 g22 g22.apo 4 g22.sym 498 4.60 g22 alga-hosting_cells - gastrodermis (2) - g22

Idents(obj) <- obj$merged_clusters
ordered_clusters <- read.table(file = "ad.merged_clusters_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
obj@active.ident <- factor(Idents(obj), levels = ordered_clusters)

# dotplot params
dotplot_params <- list(
  thr.p_val = 0.01,
  thr.avg_log2FC = 1,
  gene.alias = NULL,
  cluster.alias = NULL
)

yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
dotplot_params$gene.alias<-yn_map
xn_map <- read.csv("05.ad2all_merged_clusters_summary_0.40.tsv", header = TRUE, sep = "\t")
dotplot_params$cluster.alias <-xn_map


## goseq params
goseq_params <- list(
  gene.vector = NULL, # named vecotor {aten-0001s.g1: 0/1}
  gene2go = NULL,
  gene.lengths = NULL,
  top_n = 20
)
# load gene.lengths
gl<-read.csv('ad.gene.length.tsv', sep='\t', header=FALSE)
gene.lengths <- setNames(gl$V2,  gsub("_", "-", gl$V1))

# load gene2go
df <- read.csv('ad.gene2go.tsv', sep='\t', header=FALSE)
df$V1 <- gsub("_", "-", df$V1)
gene2go <- split(df$V2, df$V1)

# assemble params and run goseq
goseq_params$gene.lengths <- gene.lengths
goseq_params$gene2go <- gene2go

# function call
# 32 g22 g22.apo 4 g22.sym 498 4.60 g22 alga-hosting_cells - gastrodermis (2) - g22
# 28 g15 g15.apo 237 g15.sym 613 0.95 g15 gastrodermis (5) - g15
# 31 g46 g46.apo 58 g46.sym 107 0.60 g46 gastrodermis (4) - g46
# 26 g12 g12.apo 333 g12.sym 590 0.57 g12 gastrodermis (4) - g12
# 29 g21 g21.apo 209 g21.sym 347 0.51 g21 gastrodermis (5) - g21
# 27 g25 g25.apo 233 g25.sym 191 -0.20 g25 gastrodermis (4) - g25
# 30 g24 g24.apo 421 g24.sym 4 -4.44 g24 gastrodermis (4) - g24

# 26 g12 g12.apo 333 g12.sym 590 0.57 g12 gastrodermis (4) - g12
# 27 g25 g25.apo 233 g25.sym 191 -0.20 g25 gastrodermis (4) - g25
# 28 g15 g15.apo 237 g15.sym 613 0.95 g15 gastrodermis (5) - g15
# 29 g21 g21.apo 209 g21.sym 347 0.51 g21 gastrodermis (5) - g21
# 30 g24 g24.apo 421 g24.sym 4 -4.44 g24 gastrodermis (4) - g24
# 31 g46 g46.apo 58 g46.sym 107 0.60 g46 gastrodermis (4) - g46
# 32 g22 g22.apo 4 g22.sym 498 4.60 g22 alga-hosting_cells - gastrodermis (2) - g22

source('coral_ppl.R')
## cascade.merged
pair_de_goseq(obj, 'g12', 'g24', dotplot_params, goseq_params, 'ad3456')
pair_de_goseq(obj, 'g24', 'g12', dotplot_params, goseq_params, 'ad3456')

pair_de_goseq(obj, 'g22', 'g46', dotplot_params, goseq_params, 'ad3456')
pair_de_goseq(obj, 'g46', 'g22', dotplot_params, goseq_params, 'ad3456')

pair_de_goseq(obj, 'g12', 'g22', dotplot_params, goseq_params, 'ad3456')
pair_de_goseq(obj, 'g22', 'g12', dotplot_params, goseq_params, 'ad3456')




# generate merged wgcna dotplot #####################################
Idents(obj) <- obj$sym_split
sum(unique(Idents(obj)) != unique(obj$sym_split))

w_mat <- wgcna_ppl_1(obj)
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 12, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
#load('wgcna.bwnet.at345.rd') # yield bwnet
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))


# load cluster orders, gene alias, cluster alias
ordered_clusters <- read.table(file = "at.sym_split_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, goseq_params, 'at345')

# merged nj 
# get x_ordered neighbor joining
obj_nj <- nj(as.dist(as.matrix(read.table("ad3456.wmat.merged_clusters.tsv", sep='\t', header = TRUE))))

branch_colors <- rep("gray", length(obj_nj$edge.length))
branch_colors[50:64] <- "red"

branch_thicknesses <- rep(1, length(obj_nj$edge.length))
branch_thicknesses[50:64] <- 6

plot(obj_nj, 'u',  edge.color = branch_colors, lwd = branch_thicknesses, cex = 0.9, no.margin = TRUE)





ordered_clusters <- ordered_tips_apenj(obj_nj)
#ordered_clusters <- read.table(file = "ad.sym_split_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
obj@active.ident <- factor(Idents(obj), levels = ordered_clusters)





#######################################
# 241021 ad marker gene between two clusters
Idents(obj) <- obj$merged_clusters
fold.difference = 2
adj.p.value = 0.01
de.genes.df <- FindMarkers(obj = obj, ident.1 = "g25", ident.2 = "g22", print.bar = F, min.pct = 0.25, assay = "RNA", verbose = F)
de.genes <- rownames(de.genes.df[de.genes.df$p_val_adj < 0.05 & de.genes.df$avg_log2FC > 0.25, ])

# goseq
goseq_params <- list(
  gene.vector = NULL, # named vecotor {aten-0001s.g1: 0/1}
  gene2go = NULL,
  gene.lengths = NULL,
  top_n = 20
)

# load gene.lengths
gl<-read.csv('ad.gene.length.tsv', sep='\t', header=FALSE)
gene.lengths <- setNames(gl$V2,  gsub("_", "-", gl$V1))

# load gene2go
df <- read.csv('ad.gene2go.tsv', sep='\t', header=FALSE)
df$V1 <- gsub("_", "-", df$V1)
gene2go <- split(df$V2, df$V1)

# assemble params and run goseq
goseq_params$gene2go <- gene2go
goseq_params$gene.lengths <- gene.lengths
all.genes <- names(gene.lengths)
gene.vector <- as.integer(all.genes %in% de.genes)
names(gene.vector) <- all.genes			
goseq_params$gene.vector <- gene.vector
go_plot <- go_enrichment(goseq_params)	
print(go_plot)
# goseq_clusters_dotplots(obj, goseq_params, 'ad3456')



### debug
goseq_params$gene.vector
goseq_params$gene2go
goseq_params$gene.lengths
goseq_params$top_n


gene.alias<-yn_map
gene.alias$gene<-gsub("_","-", gene.alias$geneID)
colnames(de_alias)[colnames(de_alias)=="geneID"] <- "gene"
de_alias<-merge(de_m_filtered, gene.alias, by="gene", all.x = TRUE)

g2g <- read.table("ad.gene2go.tsv", header = FALSE, sep = "\t", col.names = c("gene", "go"))
g2g$gene <- gsub("_", "-", g2g$gene)
de_alias<-merge(de_alias, g2g, by="gene", all.x = TRUE)






####################################
# goseq ad3456 sym_split
load('ad3456.clean.sym_split.rd')

# no merge (seruat clusters) ppl
obj <- ad3456.clean
table(obj$sym_split)

# goseq on sym_split
Idents(obj) <- obj$sym_split

# goseq
goseq_params <- list(
  gene.vector = NULL, # named vecotor {aten-0001s.g1: 0/1}
  gene2go = NULL,
  gene.lengths = NULL,
  top_n = 20
)

# load gene.lengths
gl<-read.csv('ad.gene.length.tsv', sep='\t', header=FALSE)
gene.lengths <- setNames(gl$V2,  gsub("_", "-", gl$V1))

# load gene2go
df <- read.csv('ad.gene2go.tsv', sep='\t', header=FALSE)
df$V1 <- gsub("_", "-", df$V1)
gene2go <- split(df$V2, df$V1)

# assemble params and run goseq
goseq_params$gene2go <- gene2go
goseq_params$gene.lengths <- gene.lengths
# goseq_clusters_dotplots(obj, goseq_params, 'ad3456')

# get x_ordered neighbor joining
obj_nj <- nj(as.dist(as.matrix(read.table("ad3456.wmat.sym_split.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
#ordered_clusters <- read.table(file = "ad.sym_split_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
obj@active.ident <- factor(Idents(obj), levels = ordered_clusters)

yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.ad2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")

source('coral_ppl.R')
# generate_cluster_dotplots_withmaps(obj, yn_map, xn_map, ordered_clusters, goseq_params, 'ad3456')
generate_cluster_dotplots_withmaps(obj, yn_map, xn_map, ordered_clusters, goseq_params, 'ad3456')

#generate_cluster_dotplots_withmaps(obj, yn_map, xn_map, sort(ordered_clusters), 'ad3456')











# get x_ordered neighbor joining
obj_nj <- nj(as.dist(as.matrix(read.table("ad3456.wmat.merged_clusters.tsv", sep='\t', header = TRUE))))
branch_colors <- rainbow(nrow(obj_nj$edge))
plot(obj_nj,'u', edge.color = branch_colors, cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
write.table(ordered_clusters , file = "ad.merged_clusters_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)

 


####################################################
# 10/13/2024 goseq clear
# output.klog2::

library(Matrix)
library(GO.db)
library(goseq)
library(tidyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(ape)
library(reticulate)
library(DoubletFinder)
library(glue)
library(readr)
library(doParallel)
library(WGCNA)
library(gridExtra)
library(GO.db)
library(goseq)

setwd('C:\\Users\\kjia\\workspace\\coral\\stage.acropora_raw')


# at generate cluster dotplot #####################################
rm(list=ls())
load('at345.clean.sym_split.rd')
source('coral_ppl.R')

obj<-at345.clean
Idents(obj) <- obj$merged_clusters

# goseq
goseq_params <- list(
  gene.vector = NULL, # named vecotor {aten-0001s.g1: 0/1}
  gene2go = NULL,
  gene.lengths = NULL,
  top_n = 20
)

# load gene.lengths
gl<-read.csv('at.gene.length.tsv', sep='\t', header=FALSE)
gene.lengths <- setNames(gl$V2,  gsub("_", "-", gl$V1))

# load gene2go
df <- read.csv('at.gene2go.tsv', sep='\t', header=FALSE)
df$V1 <- gsub("_", "-", df$V1)
gene2go <- split(df$V2, df$V1)

goseq_params$gene2go <- gene2go
goseq_params$gene.lengths <- gene.lengths

goseq_clusters_dotplots(obj, goseq_params, 'at345')


# generate wgcna dotplot #####################################
Idents(obj) <- obj$sym_split
sum(unique(Idents(obj)) != unique(obj$sym_split))

w_mat <- wgcna_ppl_1(obj)
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 12, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
#load('wgcna.bwnet.at345.rd') # yield bwnet
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))


# load cluster orders, gene alias, cluster alias
ordered_clusters <- read.table(file = "at.sym_split_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, goseq_params, 'at345')




# ad generate cluster dotplot #####################################
rm(list=ls())
source('coral_ppl.R')
load('ad3456.clean.sym_split.rd')


obj<-ad3456.clean
Idents(obj) <- obj$merged_clusters

# goseq
goseq_params <- list(
  gene.vector = NULL, # named vecotor {aten-0001s.g1: 0/1}
  gene2go = NULL,
  gene.lengths = NULL,
  top_n = 20
)

# load gene.lengths
gl<-read.csv('ad.gene.length.tsv', sep='\t', header=FALSE)
gene.lengths <- setNames(gl$V2,  gsub("_", "-", gl$V1))

# load gene2go
df <- read.csv('ad.gene2go.tsv', sep='\t', header=FALSE)
df$V1 <- gsub("_", "-", df$V1)
gene2go <- split(df$V2, df$V1)

# assemble params and run goseq
goseq_params$gene2go <- gene2go
goseq_params$gene.lengths <- gene.lengths
goseq_clusters_dotplots(obj, goseq_params, 'ad3456')


# generate wgcna dotplot #####################################
Idents(obj) <- obj$sym_split
sum(unique(Idents(obj)) != unique(obj$sym_split))

w_mat <- wgcna_ppl_1(obj)
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 10, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))

# load cluster orders, gene alias, cluster alias
ordered_clusters <- read.table(file = "ad.sym_split_tree_order.tsv", sep = "\t", header=FALSE, stringsAsFactors = FALSE)[,1]
yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.ad2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, goseq_params, 'ad3456')


# debug ########################

Idents(obj) <- obj$sym_split
obj_nj <- nj(as.dist(as.matrix(read.table("ad3456.wmat.sym_split.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
write.table(ordered_clusters , file = "ad.sym_split_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)



levels(obj)
de.genes.df <- FindMarkers(obj, ident.1="g14", only.pos = T)
all.genes <- names(goseq_params$gene.length)
de.genes <- rownames(de.genes.df[de.genes.df$p_val_adj < 0.05 & de.genes.df$avg_log2FC > 0.25, ])
gene.vector <- as.integer(all.genes %in% de.genes)
names(gene.vector) <- all.genes		
goseq_params$gene.vector <- gene.vector
params<-goseq_params
pwf <- nullp(params$gene.vector, bias.data=params$gene.lengths, plot.fit = FALSE)








# ppl without merging #####################################
rm(list=ls())
source('coral_ppl.R')
load('at345.clean.sym_split.rd')

# generate new clusters
obj <- at345.clean
# merged_clusters
default_clusters <- obj$merged_clusters
obj$merged_split <- ifelse(obj$sym_label == 'apo',paste0(default_clusters, '.apo'), 
                        ifelse(obj$sym_label == 'sym', paste0(default_clusters, '.sym'),
                               as.character(default_clusters)))
# seurat_clusters
default_clusters <- obj$seurat_clusters
obj$seurat_split <- ifelse(obj$sym_label == 'apo',paste0(default_clusters, '.apo'), 
                           ifelse(obj$sym_label == 'sym', paste0(default_clusters, '.sym'),
                                  as.character(default_clusters)))

# export obj for python
export_mtx(obj, c("merged_clusters", "merged_split", "seurat_clusters", "seurat_split"), 'at345')
at345.clean <- obj
save(at345.clean, file='at345.clean.sym_split_v2.rd')


# h5ad -> sam -> samap -> {cluster annotation & alias} {pc_projected average expression wmat}

# get cluster order according to wmat nj
Idents(obj) <- obj$merged_clusters
obj_nj <- nj(as.dist(as.matrix(read.table("at345.wmat.merged_clusters.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)


# generate cluster dotplot #####################################
# goseq
goseq_params <- list(
  gene.vector = NULL, # named vecotor {aten-0001s.g1: 0/1}
  gene2go = NULL,
  gene.lengths = NULL,
  top_n = 20
)

# load gene.lengths
gl<-read.csv('at.gene.length.tsv', sep='\t', header=FALSE)
gene.lengths <- setNames(gl$V2,  gsub("_", "-", gl$V1))

# load gene2go
df <- read.csv('at.gene2go.tsv', sep='\t', header=FALSE)
df$V1 <- gsub("_", "-", df$V1)
gene2go <- split(df$V2, df$V1)

goseq_params$gene2go <- gene2go
goseq_params$gene.lengths <- gene.lengths


# x,y annotation
yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_merged_clusters_summary_0.40.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmaps(obj, yn_map, xn_map, ordered_clusters, goseq_params, 'at345')







sapply(goseq_params, head)

# go_enrichment() unit test
markers <- FindMarkers(at345.clean, ident.1 = "g12", only.pos = T)
de.genes <- rownames(markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.25, ])
all.genes <- gsub("_", "-", gl$V1)
gene.vector <- as.integer(all.genes %in% de.genes)
names(gene.vector) <- all.genes
goseq_params$gene.vector <- gene.vector
# print(go_enrichment(goseq_params))

go_plot <- go_enrichment(goseq_params)
print(go_plot)

sapply(goseq_params, head)






  



# neighbor joining for label order
## calculating wmat for ad
## s=SAM()
## s.load_data('01.at345.sam.h5ad')
## pmat = us._calc_distance_matrix(s, 'sym_split', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('at345.wmat.sym_split.tsv', sep='\t', index=False, header=True)







# no merging
##################
# restore original seurat clusters
# Idents(obj) <- obj$seurat_clusters
# obj <- merge_clusters3j(obj, thresh_corr = 0.8)

default_clusters <- Idents(obj)
obj$sym_split <- ifelse(obj$sym_label == 'apo',paste0(default_clusters, '.apo'), 
                        ifelse(obj$sym_label == 'sym', paste0(default_clusters, '.sym'),
                               as.character(default_clusters)))










##############
 # 12 for at too
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 12, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
save(bwnet, file='wgcna.bwnet.at345.rd')

## wgcna dotplot
## ranking genes for module membership
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))


Idents(obj) <- obj$sym_split
obj_nj <- nj(as.dist(as.matrix(read.table("at345.wmat.sym_split.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
write.table(ordered_clusters , file = "at.sym_split_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
#at345.clean@active.ident <- factor(Idents(at345.clean), levels = ordered_clusters)

unique(Idents(obj))==unique(obj$sym_split)
yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, 'at345')
##############














####################################################
# 10/05/2024 goseq
# output.klog2::

# BiocManager::install("org.Dm.eg.db")
# BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
# library(org.Dm.eg.db)  # or another organism-specific package
# library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

BiocManager::install("org.Hs.eg.db")  # for human annotation (use relevant annotation package)
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(GO.db)
library(goseq)
#library(org.Hs.eg.db)  # or another organism-specific package
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# check supported organisms
supportedOrganisms()

gl<-read.csv('at.gene.length.tsv', sep='\t', header=FALSE)
gene.lengths <- setNames(gl$V2, gl$V1)

# get all.gene.vector
all.genes <- gsub("_", "-", gl$V1)


markers <- FindMarkers(at345.clean, ident.1 = "g12", only.pos = T)
de.genes <- rownames(markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.25, ])

# Create a binary vector: 1 if the gene is DE, 0 if not
gene.vector <- as.integer(all.genes %in% de.genes)
names(gene.vector) <- all.genes


# load gene2go
df <- read.csv('at.gene2go.tsv', sep='\t', header=FALSE)
df$V1 <- gsub("_", "-", df$V1)
gene2go <- split(df$V2, df$V1)


pwf <- nullp(gene.vector, bias.data=gene.lengths, plot.fit = FALSE)
res <- goseq(pwf = pwf, gene2cat = gene2go, method = "Hypergeometric")
res <- res[!is.na(res$term),]

go_plot<- res %>% top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
             geom_point() +
             expand_limits(x=0) +
             labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

print(go_plot)
outfolder="output/"
outname <- paste0(outfolder, outprefix, "_c_", i, "_dotplot.pdf")

ggsave(filename = outname,
       plot = go_plot, 
       path = outfolder,
       width = 18,
       height = (length(genes_for_dotplot) * (1/6)) + 2, 
       units = "in", 
       limitsize = F)


####################################################
# 09/27/2024 re-run clear ppl at no merging
# output.klog2::- coral acropora project clear ppl

# at
load("at345.clean.sym_split.rd") # at345.clean
at345.nm <- at345.clean
obj<-at345.nm
Idents(obj) <- obj$seurat_clusters
Idents(obj) <- paste0("g", as.character(Idents(obj)))

default_clusters <- Idents(obj)
obj$sym_split <- ifelse(obj$sym_label == 'apo',paste0(default_clusters, '.apo'), 
                          ifelse(obj$sym_label == 'sym', paste0(default_clusters, '.sym'),
                          as.character(default_clusters)))


# save results
table_df <- as.data.frame(table(obj$sym_split))
write.table(table_df, file = "at345.nm.split.cell_num_per_cluster.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
save(at345.clean, file='at345.clean.sym_split.rd')


# run samap merged_clusters vs {cell_type_family, cell_type}
# export data for samp ppl
library(Matrix)
counts<-at345.clean[["RNA"]]$counts
writeMM(obj = at345.clean[["RNA"]]$counts, file="at345.counts.mtx")
write(x = rownames(counts), file = "at345.row.txt")
write(x = colnames(counts), file = "at345.col.txt")
# Save cluster information
mdata <- at345.clean@meta.data
write.table(mdata[, c("merged_clusters", "sym_split")], file = "at345.clusters.tsv", sep = "\t", row.names = TRUE, quote = FALSE)


# python ppl for samap running
# check gene name format consistency before initialize h5 object
# maps.blast/at_to_xx.txt vs at345.row.txt
## output.klog2::## at samap procedure



############################### production analysis ###############################

# neighbor joining for label order
## calculating wmat for ad
## s=SAM()
## s.load_data('01.at345.sam.h5ad')
## pmat = us._calc_distance_matrix(s, 'sym_split', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('at345.wmat.sym_split.tsv', sep='\t', index=False, header=True)

## pmat = us._calc_distance_matrix(s, 'merged_clusters', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('at345.wmat.merged_clusters.tsv', sep='\t', index=False, header=True)
Idents(at345.clean) <- at345.clean$merged_clusters
obj_nj <- nj(as.dist(as.matrix(read.table("at345.wmat.merged_clusters.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)

write.table(ordered_clusters , file = "at.merged_clusters_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
at345.clean@active.ident <- factor(Idents(at345.clean), levels = ordered_clusters)


# dotplot with alias for all merged clusters
sum(unique(Idents(at345.clean))!=unique(at345.clean$merged_clusters))
# Idents(at345.clean) <- at345.clean$merged_clusters
# yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
# xn_map <- read.csv("05.at2all_merged_clusters_alias_0.40.tsv", header = TRUE, sep = "\t")
# generate_cluster_dotplots_withmaps(at345.clean, yn_map, xn_map, ordered_clusters, 'at345')

# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-24 11:27:44
# $ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.at2all.merged_clusters_cell_type_family.blast.tsv 0.4 ad_ 05.at2all_merged_clusters


yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_merged_clusters_summary_0.40.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmaps(at345.clean, yn_map, xn_map, ordered_clusters, 'at345')


yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_merged_clusters_alias_0.40.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmaps(at345.clean, yn_map, xn_map, ordered_clusters, 'at345')


# wgcna dotplot with alias for sym_split clusters
## load samap result to calculate sym_split cell_type_cluster mapping

# sm = load_samap('02.adhy.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'hy', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adnv.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'nv', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adsp.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'sp', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adxe.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'xe', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adnt.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'nt', 'IDs', sankey_cutoff=0.15

##
# len(np.unique(sm.sams['at'].adata.obs['sym_split']))=100
#
# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-24 17:57:31
#  $ ls 03.at*sym_split_cell_type_family.tsv|awk 'BEGIN {printf "python proc_coral_samap.py combine_scorefiles \""} {printf "%s ", $1} END{printf "\" 100 04.heatmap.at2all.sym_split_cell_type_family.blast.tsv"}' |sh


# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-24 17:57:38
# $ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.at2all.sym_split_cell_type_family.blast.tsv 0.4 at_ 05.at2all_sym_split


library(gridExtra)
obj<-at345.clean
w_mat <- wgcna_ppl_1(obj) # 12 for at too
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 12, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
save(bwnet, file='wgcna.bwnet.at345.rd')

## wgcna dotplot
## ranking genes for module membership
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))

# neighbor joining for label order
## calculating wmat for ad
## s=SAM()
## s.load_data('01.at345.sam.h5ad')
## pmat = us._calc_distance_matrix(s, 'sym_split', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('at345.wmat.sym_split.tsv', sep='\t', index=False, header=True)

Idents(obj) <- obj$sym_split
obj_nj <- nj(as.dist(as.matrix(read.table("at345.wmat.sym_split.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
write.table(ordered_clusters , file = "at.sym_split_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
#at345.clean@active.ident <- factor(Idents(at345.clean), levels = ordered_clusters)

unique(Idents(obj))==unique(obj$sym_split)
yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, 'at345')



Idents(obj) <- obj$sym_split
obj_nj <- nj(as.dist(as.matrix(read.table("at345.wmat.sym_split.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
write.table(ordered_clusters , file = "at.sym_split_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
#at345.clean@active.ident <- factor(Idents(at345.clean), levels = ordered_clusters)

unique(Idents(obj))==unique(obj$sym_split)
yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_sym_split_alias_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, 'at345')








####################################################
# 09/20/2024 re-run clear ppl
# add final std_seurat_ppl (clustering) after junk clusters removed.
# output.klog2::- coral acropora project clear ppl

# ad
load("ad3456_nosym.mc.RData")

# remove additional ribosomal genes
# confirm gene name format are consistent first
# "adig-s0001.g1" %in% rownames(ad3456_nosym.mc)
ad.ribo.genes <- readLines('02.ad.ribo.gene.txt')
ad3456.nosym.nr <- subset(ad3456_nosym.mc, features = setdiff(rownames(ad3456_nosym.mc) , ad.ribo.genes))

# remove junk clusters
# Idents(ad3456.nosym.nr)<-ad3456.nosym.nr$merged_clusters
ad3456.clean_0 <- subset(ad3456.nosym.nr, cells = setdiff(Cells(ad3456.nosym.nr), WhichCells(ad3456.nosym.nr, ident = c('g4', 'g12'))))

# re-initialize a new seurat object from clean data
ad3456.clean_1 <- CreateSeuratObject(counts = ad3456.clean_0@assays$RNA$counts)
ad3456.clean_1$orig.ident <- ad3456_nosym.mc$orig.ident


# run std_seurat_ppl
ad3456.std <- std_seurat_ppl(ad3456.clean_1)
ad3456.std <- RunUMAP(ad3456.std, dims = 1:40, min.dist = 0.1, umap.method = "umap-learn", metric = "cosine", verbose = FALSE)

obj<-ad3456.std
# merge clusters
Idents(obj) <- paste0("g", as.character(Idents(obj)))
table(Idents(obj))
obj <- merge_clusters3j(obj)

# remove small clusters can't merge
obj_m1 <- subset(obj, cells = setdiff(Cells(obj), WhichCells(obj, ident = c('g68', 'g69'))))
obj_m1 <- std_seurat_ppl(obj_m1)
obj_m1 <- RunUMAP(obj_m1, dims = 1:40, min.dist = 0.1, umap.method = "umap-learn", metric = "cosine", verbose = FALSE)
Idents(obj_m1) <- paste0("g", as.character(Idents(obj_m1)))
ad3456.clean <- merge_clusters3j(obj_m1)

# add sym/apo split label
ad3456.clean$sym_label <- ifelse(ad3456.clean$orig.ident %in% c('ad323', 'ad423'), 'apo', 
                                 ifelse(ad3456.clean$orig.ident %in% c('ad523', 'ad623'), 'sym', NA)) 

default_clusters <- Idents(ad3456.clean)
ad3456.clean$sym_split <- ifelse(ad3456.clean$sym_label == 'apo',paste0(default_clusters, '.apo'), 
                                 ifelse(ad3456.clean$sym_label == 'sym', paste0(default_clusters, '.sym'),
                                 as.character(default_clusters)))

# save results
table_df <- as.data.frame(table(ad3456.clean$sym_split))
write.table(table_df, file = "ad3456.split.cell_num_per_cluster.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
save(ad3456.clean, file='ad3456.clean.sym_split.rd')


# run samap merged_clusters vs {cell_type_family, cell_type}
# export data for samp ppl
library(Matrix)
counts<-ad3456.clean[["RNA"]]$counts
writeMM(obj = ad3456.clean[["RNA"]]$counts, file="ad3456.counts.mtx")
write(x = rownames(counts), file = "ad3456.row.txt")
write(x = colnames(counts), file = "ad3456.col.txt")

# Save cluster information
mdata <- ad3456.clean@meta.data
write.table(mdata[, c("merged_clusters", "sym_split")], file = "ad3456.clusters.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
# python ppl for samap procedure



############################### production analysis ###############################

# neighbor joining for label order
## calculating wmat for ad
## s=SAM()
## s.load_data('01.ad3456.sam.h5ad')
## pmat = us._calc_distance_matrix(s, 'sym_split', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('ad3456.wmat.sym_split.tsv', sep='\t', index=False, header=True)

## pmat = us._calc_distance_matrix(s, 'merged_clusters', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('ad3456.wmat.merged_clusters.tsv', sep='\t', index=False, header=True)

adnj <- nj(as.dist(as.matrix(read.table("ad3456.wmat.merged_clusters.tsv", sep='\t', header = TRUE))))
plot(adnj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(adnj)
write.table(ordered_cells , file = "ad.merged_clusters_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
ad3456.clean@active.ident <- factor(Idents(ad3456.clean), levels = ordered_clusters)


# dotplot with alias for all merged clusters
unique(Idents(ad3456.clean))==unique(ad3456.clean$merged_clusters)
Idents(ad3456.clean) <- ad3456.clean$merged_clusters
yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.ad2all_merged_clusters_alias_0.40.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmaps(ad3456.clean, yn_map, xn_map, ordered_clusters, 'ad3456')

# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-24 11:27:44
# $ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.ad2all.merged_clusters_cell_type_family.blast.tsv 0.4 ad_ 05.ad2all_merged_clusters
# 2024-09-24 11:27:46|14090|0|INFO|ignore absent libraries
# 2024-09-24 11:27:46|14090|0|INFO|save alias to 05.ad2all_merged_clusters_alias_0.40.tsv
# 2024-09-24 11:27:46|14090|0|INFO|save summary to 05.ad2all_merged_clusters_desc_0.40.tsv

yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.ad2all_merged_clusters_summary_0.40.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmaps(ad3456.clean, yn_map, xn_map, ordered_clusters, 'ad3456')

# wgcna dotplot with alias for sym_split clusters
## load samap result to calculate sym_split cell_type_cluster mapping

# sm = load_samap('02.adhy.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'hy', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adnv.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'nv', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adsp.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'sp', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adxe.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'xe', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adnt.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'nt', 'IDs', sankey_cutoff=0.15

##
# len(np.unique(sm.sams['ad'].adata.obs['sym_split']))=133
#
# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-24 17:57:31
#  $ ls 03*sym_split_cell_type_family.tsv|awk 'BEGIN {printf "python proc_coral_samap.py combine_scorefiles \""} {printf "%s ", $1} END{printf "\" 133 04.heatmap.ad2all.sym_split_cell_type_family.blast.tsv"}' |sh
# 2024-09-24 17:57:38|14274|0|INFO|ignore absent libraries
# 2024-09-24 17:57:38|14274|0|INFO|5 files loaded.
# 2024-09-24 17:57:38|14274|0|INFO|save all scores (133, 53) to 04.heatmap.ad2all.sym_split_cell_type_family.blast.tsv

# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-24 17:57:38
# $ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.ad2all.sym_split_cell_type_family.blast.tsv 0.4 ad_ 05.ad2all_sym_split
# 2024-09-24 18:00:57|14275|0|INFO|ignore absent libraries
# 2024-09-24 18:00:57|14275|0|INFO|save alias to 05.ad2all_sym_split_alias_0.40.tsv
# 2024-09-24 18:00:57|14275|0|INFO|save summary to 05.ad2all_sym_split_summary_0.40.tsv


# load('ad3456.clean.sym_split.rd')
#load('ad3456.clean.sym_split.rd')
# obj<- ad3456.clean

library(gridExtra)
obj<-ad3456.clean
w_mat <- wgcna_ppl_1(obj) # 12 for ad
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 12, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
save(bwnet, file='wgcna.bwnet.ad3456.rd')

## wgcna dotplot
## ranking genes for module membership
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))

# neighbor joining for label order
## calculating wmat for ad
## s=SAM()
## s.load_data('01.ad3456.sam.h5ad')
## pmat = us._calc_distance_matrix(s, 'sym_split', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('ad3456.wmat.sym_split.tsv', sep='\t', index=False, header=True)

# load('ad3456.clean.sym_split.rd')
# obj<- ad3456.clean

Idents(obj) <- obj$sym_split
obj_nj <- nj(as.dist(as.matrix(read.table("ad3456.wmat.sym_split.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
write.table(ordered_cells , file = "ad.sym_split_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)

unique(Idents(obj)) == unique(obj$sym_split)
yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.ad2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, 'ad3456')

xn_map <- read.csv("05.ad2all_sym_split_alias_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, 'ad3456')



Idents(obj) <- obj$merged_clusters
obj_nj <- nj(as.dist(as.matrix(read.table("ad3456.wmat.merged_clusters.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
write.table(ordered_cells , file = "ad.sym_split_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)

unique(Idents(obj)) == unique(obj$merged_clusters)
yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.ad2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, 'ad3456')
# ------------------------------------------------------------------------------------------------------------

# at
load("at345_nosym.mc.RData")

# remove additional ribosomal genes
# manually change "_" to "-" in 02.at.ribo.gene.txt
# "aten-s0890.g1" %in% rownames(at345_nosym.mc)
at.ribo.genes <- readLines('02.at.ribo.gene.txt')
at345.nosym.nr <- subset(at345_nosym.mc, features = setdiff(rownames(at345_nosym.mc) , at.ribo.genes))

# remove junk clusters
# Idents(at345.nosym.nr)<-at345.nosym.nr$merged_clusters
at345.clean_0 <- subset(at345.nosym.nr, cells = setdiff(Cells(at345.nosym.nr), WhichCells(at345.nosym.nr, ident = 'g16')))

# re-initialize a new seurat object from clean data
at345.clean_1 <- CreateSeuratObject(counts = at345.clean_0@assays$RNA$counts)
at345.clean_1$orig.ident <- at345_nosym.mc$orig.ident


# run std_seurat_ppl
at345.std <- std_seurat_ppl(at345.clean_1)
at345.std <- RunUMAP(at345.std, dims = 1:40, min.dist = 0.1, umap.method = "umap-learn", metric = "cosine", verbose = FALSE)

obj<-at345.std
# merge clusters
Idents(obj) <- paste0("g", as.character(Idents(obj)))
table(Idents(obj))
obj <- merge_clusters3j(obj)

# remove small clusters can't merge
cn<-table(obj$merged_clusters)
cn[cn<10]
# g53 g54 g55 
# 7   3   2  
obj_m1 <- subset(obj, cells = setdiff(Cells(obj), WhichCells(obj, ident = c('g53', 'g54', 'g55'))))

obj_new <- CreateSeuratObject(counts = obj_m1@assays$RNA$counts)
obj_new$orig.ident <- at345_nosym.mc$orig.ident

obj_new <- std_seurat_ppl(obj_new)
obj_new <- RunUMAP(obj_new, dims = 1:40, min.dist = 0.1, umap.method = "umap-learn", metric = "cosine", verbose = FALSE)
Idents(obj_new) <- paste0("g", as.character(Idents(obj_new)))
at345.clean <- merge_clusters3j(obj_new)

cn<-table(at345.clean$merged_clusters)
cn[cn<10]

# add sym/apo cluster assignments based on merged_clusters
at345.clean$sym_label <- ifelse(at345.clean$orig.ident == 'at3', 'apo', 
                                ifelse(at345.clean$orig.ident %in% c('at4', 'at5'), 'sym', NA)) 

default_clusters <- Idents(at345.clean)
at345.clean$sym_split <- ifelse(at345.clean$sym_label == 'apo',paste0(default_clusters, '.apo'), 
                                ifelse(at345.clean$sym_label == 'sym', paste0(default_clusters, '.sym'),
                                       as.character(default_clusters)))


# save results
table_df <- as.data.frame(table(at345.clean$sym_split))
write.table(table_df, file = "at345.split.cell_num_per_cluster.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
save(at345.clean, file='at345.clean.sym_split.rd')


# run samap merged_clusters vs {cell_type_family, cell_type}
# export data for samp ppl
library(Matrix)
counts<-at345.clean[["RNA"]]$counts
writeMM(obj = at345.clean[["RNA"]]$counts, file="at345.counts.mtx")
write(x = rownames(counts), file = "at345.row.txt")
write(x = colnames(counts), file = "at345.col.txt")
# Save cluster information
mdata <- at345.clean@meta.data
write.table(mdata[, c("merged_clusters", "sym_split")], file = "at345.clusters.tsv", sep = "\t", row.names = TRUE, quote = FALSE)


# python ppl for samap running
# check gene name format consistency before initialize h5 object
# maps.blast/at_to_xx.txt vs at345.row.txt
## output.klog2::## at samap procedure



############################### production analysis ###############################

# neighbor joining for label order
## calculating wmat for ad
## s=SAM()
## s.load_data('01.at345.sam.h5ad')
## pmat = us._calc_distance_matrix(s, 'sym_split', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('at345.wmat.sym_split.tsv', sep='\t', index=False, header=True)

## pmat = us._calc_distance_matrix(s, 'merged_clusters', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('at345.wmat.merged_clusters.tsv', sep='\t', index=False, header=True)
Idents(at345.clean) <- at345.clean$merged_clusters
obj_nj <- nj(as.dist(as.matrix(read.table("at345.wmat.merged_clusters.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)

write.table(ordered_clusters , file = "at.merged_clusters_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
at345.clean@active.ident <- factor(Idents(at345.clean), levels = ordered_clusters)


# dotplot with alias for all merged clusters
unique(Idents(at345.clean))==unique(at345.clean$merged_clusters)
# Idents(at345.clean) <- at345.clean$merged_clusters
# yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
# xn_map <- read.csv("05.at2all_merged_clusters_alias_0.40.tsv", header = TRUE, sep = "\t")
# generate_cluster_dotplots_withmaps(at345.clean, yn_map, xn_map, ordered_clusters, 'at345')

# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-24 11:27:44
# $ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.at2all.merged_clusters_cell_type_family.blast.tsv 0.4 ad_ 05.at2all_merged_clusters


yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_merged_clusters_summary_0.40.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmaps(at345.clean, yn_map, xn_map, ordered_clusters, 'at345')


yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_merged_clusters_alias_0.40.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmaps(at345.clean, yn_map, xn_map, ordered_clusters, 'at345')


# wgcna dotplot with alias for sym_split clusters
## load samap result to calculate sym_split cell_type_cluster mapping

# sm = load_samap('02.adhy.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'hy', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adnv.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'nv', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adsp.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'sp', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adxe.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'xe', 'cell_type_family', sankey_cutoff=0.15)
# sm = load_samap('02.adnt.samap.blast.merged_clusters_cell_type_family.pkl')
# mt = pc.samap_alignment_and_sankey(sm, 'ad', 'sym_split', 'nt', 'IDs', sankey_cutoff=0.15

##
# len(np.unique(sm.sams['at'].adata.obs['sym_split']))=100
#
# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-24 17:57:31
#  $ ls 03.at*sym_split_cell_type_family.tsv|awk 'BEGIN {printf "python proc_coral_samap.py combine_scorefiles \""} {printf "%s ", $1} END{printf "\" 100 04.heatmap.at2all.sym_split_cell_type_family.blast.tsv"}' |sh


# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-24 17:57:38
# $ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.at2all.sym_split_cell_type_family.blast.tsv 0.4 at_ 05.at2all_sym_split


library(gridExtra)
obj<-at345.clean
w_mat <- wgcna_ppl_1(obj) # 12 for at too
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 12, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
save(bwnet, file='wgcna.bwnet.at345.rd')

## wgcna dotplot
## ranking genes for module membership
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))

# neighbor joining for label order
## calculating wmat for ad
## s=SAM()
## s.load_data('01.at345.sam.h5ad')
## pmat = us._calc_distance_matrix(s, 'sym_split', fn_cluster_vector=us._mean_expression_wpca)  # fn_cluster_vector=us._mean_expression
## pmat.to_csv('at345.wmat.sym_split.tsv', sep='\t', index=False, header=True)

Idents(obj) <- obj$sym_split
obj_nj <- nj(as.dist(as.matrix(read.table("at345.wmat.sym_split.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
write.table(ordered_clusters , file = "at.sym_split_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
#at345.clean@active.ident <- factor(Idents(at345.clean), levels = ordered_clusters)

unique(Idents(obj))==unique(obj$sym_split)
yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_sym_split_summary_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, 'at345')



Idents(obj) <- obj$sym_split
obj_nj <- nj(as.dist(as.matrix(read.table("at345.wmat.sym_split.tsv", sep='\t', header = TRUE))))
plot(obj_nj,'u', cex = 0.7)
ordered_clusters <- ordered_tips_apenj(obj_nj)
write.table(ordered_clusters , file = "at.sym_split_tree_order.tsv", sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
#at345.clean@active.ident <- factor(Idents(at345.clean), levels = ordered_clusters)

unique(Idents(obj))==unique(obj$sym_split)
yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map <- read.csv("05.at2all_sym_split_alias_0.40.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmaps(obj, bwnet, me_wt, yn_map, xn_map, 300, ordered_clusters, 'at345')


####################################################
# 09/17/2024
# adjust seurat cluster dotplot - append prost results and cluster mappings
source('coral_ppl.R')

# ad
load("ad3456.clean.sym_split.RData")
# unique(Idents(ad3456))==unique(ad3456$merged_clusters)

#kjia@DESKTOP-L0MMU09 ~/workspace/coral/stage.acropora_raw 2024-09-17 15:54:29
#$ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.ad2all.merged_clusters_cell_type_family.blast.tsv 0.4 ad_ 05.ad2all_merged_clusters_alias_04.tsv
#2024-09-17 15:56:00|10862|2|INFO|ignore absent libraries
#2024-09-17 15:56:01|10862|0|INFO|save alias to 05.ad2all_merged_clusters_alias_04.tsv

yn_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t") 
xn_map<- read.csv("05.ad2all_merged_clusters_alias_04.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmaps(obj, yn_map, xn_map, 'ad3456')


# at
load("at345.clean.sym_split.RData")
unique(Idents(obj))==unique(obj$merged_clusters)
Idents(obj) <- obj$merged_clusters

#kjia@DESKTOP-L0MMU09 ~/workspace/coral/stage.acropora_raw 2024-09-17 15:54:29
#$ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.at2all.merged_clusters_cell_type_family.blast.tsv 0.4 at_ 05.at2all_merged_clusters_alias_04.tsv
#2024-09-17 15:56:00|10862|2|INFO|ignore absent libraries
#2024-09-17 15:56:01|10862|0|INFO|save alias to 05.at2all_merged_clusters_alias_04.tsv

yn_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t") 
xn_map<- read.csv("05.at2all_merged_clusters_alias_04.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmaps(obj, yn_map, xn_map, 'at345')

####################################################
# 09/16/2024
# adjust wgcna dotplot - append prost results

# at345 ---------------------------------------------
load('at345.clean.sym_split.RData')
at345<-obj

load('wgcna.bwnet.at345.RData')
load('wgcna.me_wt.at345.RData')

Idents(obj) <- obj$sym_split
name_map <- read.csv("at.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map<- read.csv("04.at2all_cluster_alias_04.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmap(bwnet, me_wt, name_map, xn_map, 300, 'at345')



# ad3456
load('ad3456.clean.sym_split.RData')
ad3456<-obj

load('wgcna.bwnet.ad3456.RData')
load('wgcna.me_wt.ad3456.RData')

Idents(obj) <- obj$sym_split
name_map <- read.csv("ad.gene_alias.tsv", header = TRUE, sep = "\t")
xn_map<- read.csv("04.ad2all_cluster_alias_04.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmap(bwnet, me_wt, name_map, xn_map, 300, 'ad3456')



####################################################
# 09/11/2024
# adjust wgcna dotplot
source('coral_ppl.R')

# ad3456 ---------------------------------------------
load('ad3456.clean.sym_split.RData')
ad3456<-obj

# determine optimal power
library(gridExtra)
w_mat <- wgcna_ppl_1(obj) # choose 10 for ad, 12 for at
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 10, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
save(bwnet, file='wgcna.bwnet.ad3456.RData')

# dotplot
# ranking genes for module membership
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))
save(me_wt, file='wgcna.me_wt.ad3456.RData')

# kjia@DESKTOP-L0MMU09 ~/workspace/src/contactGroups 2024-09-12 18:23:56
# $ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.ad2all.sym_split_cell_type_family.blast.tsv 0.4 ad_ 04.ad2all_cluster_alias_04.tsv
# 2024-09-12 18:25:23|2487|0|INFO|ignore absent libraries
# 2024-09-12 18:25:23|2487|0|INFO|save alias to 04.ad2all_cluster_alias_04.tsv

Idents(obj) <- obj$sym_split
name_map <- read.csv("adig_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
xn_map<- read.csv("04.ad2all_cluster_alias_04.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmap(bwnet, me_wt, name_map, xn_map, 300, 'ad3456')



# at345 ---------------------------------------------
load('at345.clean.sym_split.RData')
at345<-obj

# determine optimal power
library(gridExtra)
w_mat <- wgcna_ppl_1(obj) # choose 10 for ad, 12 for at
bwnet <- blockwiseModules(w_mat, maxBlockSize = 6000, TOMType = "signed", power = 12, mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 0, verbose = 3)
save(bwnet, file='wgcna.bwnet.at345.RData')

# dotplot
# ranking genes for module membership
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))
save(me_wt, file='wgcna.me_wt.at345.RData')

# kjia@DESKTOP-L0MMU09 ~/workspace/coral/samap 2024-09-12 19:10:03
# $ python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.at2all.sym_split_cell_type_family.blast.tsv 0.4 at_ 04.at2all_cluster_alias_04.tsv
# 2024-09-12 19:10:16|2639|0|INFO|ignore absent libraries
# 2024-09-12 19:10:16|2639|0|INFO|save alias to 04.at2all_cluster_alias_04.tsv


Idents(obj) <- obj$sym_split
name_map <- read.csv("aten_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
xn_map<- read.csv("04.at2all_cluster_alias_04.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmap(bwnet, me_wt, name_map, xn_map, 300, 'at345')





# pheatmap -----------------------------------------------
library(pheatmap)
infile="04.heatmap.at2all.merged_clusters_cell_type_family.blast.tsv"
infile="04.heatmap.ad2all.merged_clusters_cell_type_family.blast.tsv"
infile="04.heatmap.at2all.sym_split_cell_type_family.blast.tsv"
infile="04.heatmap.ad2all.sym_split_cell_type_family.blast.tsv"

h=pheatmap(t(read.table(infile, row.names=1, header=TRUE, sep='\t')),cluster_rows=F, cluster_cols=T)


####################################################
# 202409

# restore project
# rstudio::New -> project-> [path]

install.packages("tidyr")
install.packages("Seurat")
install.packages("tibble")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("ape")
install.packages("reticulate")
install.packages("glue")
install.packages("readr")
install.packages("WGCNA")

install.packages('devtools')
install.packages("C:\\cygwin64\\home\\kjia\\workspace\\library\\repository\\presto")
install.packages("C:\\cygwin64\\home\\kjia\\workspace\\library\\repository\\DoubletFinder")

install.packages("BiocManager")
BiocManager::install("WGCNA")

library(tidyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(ape)

library(DoubletFinder)
library(glue)
library(readr)
library("WGCNA")
library(gridExtra)
library(Matrix)
# create r-reticulate using python3.8
# (base) C:\Users\kjia>conda create --name r-reticulate python=3.8
# conda activate r-reticulate
# pip install numpy scipy pandas leidenalg

library(reticulate)
reticulate::py_config()
py_module_available('leidenalg')


source("coral_ppl.R")
load("at345_nosym.mc.RData")
name_map <- read.csv("aten_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
#generate_cluster_dotplots_withmap(at345_nosym.mc, name_map, "at345.mc")
# remove junk clusters
at345.clean <- subset(at345_nosym.mc, cells = setdiff(Cells(at345_nosym.mc), WhichCells(at345_nosym.mc, ident = 'g16')))

# add apo/sym label
at345.clean$sym_label <- ifelse(at345.clean$orig.ident == 'at3', 'apo', 
                      ifelse(at345.clean$orig.ident %in% c('at4', 'at5'), 'sym', NA)) 
obj <- at345.clean

default_clusters <- Idents(obj)
obj$sym_split <- ifelse(obj$sym_label == 'apo',paste0(default_clusters, '.apo'),
                        ifelse(obj$sym_label == 'sym', paste0(default_clusters, '.sym'),
                               as.character(default_clusters)))

save(obj, file='at345.clean.sym_split.RData')
at345.split <- obj

table_df <- as.data.frame(table(obj$sym_split))
write.table(table_df, file = "at345.split.cell_num_per_cluster.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

library(Matrix)
counts<-at345.split[["RNA"]]$counts
writeMM(obj = at345.split[["RNA"]]$counts, file="at345.counts.mtx")
write(x = rownames(counts), file = "at345.row.txt")
write(x = colnames(counts), file = "at345.col.txt")

# Save cluster information if needed by your SAMap analysis

mdata <- at345.split@meta.data
write.table(mdata[, c("merged_clusters", "sym_split")], file = "at345.clusters.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

# ad
load("ad3456_nosym.mc.RData")
name_map <- read.csv("adig_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
#generate_cluster_dotplots_withmap(ad3456_nosym.mc, name_map, "ad3456.mc")
# remove junk clusters
ad3456.clean <- subset(ad3456_nosym.mc, cells = setdiff(Cells(ad3456_nosym.mc), WhichCells(ad3456_nosym.mc, ident = c('g4', 'g12'))))

# add apo/sym label
ad3456.clean$sym_label <- ifelse(ad3456.clean$orig.ident %in% c('ad323', 'ad423'), 'apo', 
                                ifelse(ad3456.clean$orig.ident %in% c('ad523', 'ad623'), 'sym', NA)) 


# split by apo/sym
obj <- ad3456.clean
default_clusters <- Idents(obj)
obj$sym_split <- ifelse(obj$sym_label == 'apo',paste0(default_clusters, '.apo'),
                      ifelse(obj$sym_label == 'sym', paste0(default_clusters, '.sym'),
                       as.character(default_clusters)))

table_df <- as.data.frame(table(obj$sym_split))
write.table(table_df, file = "ad3456.split.cell_num_per_cluster.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
save(obj, file='ad3456.clean.sym_split.RData')

library(Matrix)
ad3456.split <- obj
counts<-ad3456.split[["RNA"]]$counts
writeMM(obj = ad3456.split[["RNA"]]$counts, file="ad3456.counts.mtx")
write(x = rownames(counts), file = "ad3456.row.txt")
write(x = colnames(counts), file = "ad3456.col.txt")

# Save cluster information if needed by your SAMap analysis
mdata <- ad3456.split@meta.data
write.table(mdata[, c("merged_clusters", "sym_split")], file = "ad3456.clusters.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

# use python to calculate pmat
# neighbor joining
# par("mar")
# par("mar=c(1,1,1,1))
# if figure margins too large
plot(nj(as.dist(as.matrix(read.table("at345.wmat.tsv", sep='\t', header = TRUE)))),'u', cex = 0.7)
plot(nj(as.dist(as.matrix(read.table("at345.wmat.merged_clusters.tsv", sep='\t', header = TRUE)))),'u', cex = 0.7)

# ---------------------------------------------------------------
# skip as function---------------------------------------------------------------
# wgcna normal ppl
g6000_seu <- FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = 6000, binning.method = "equal_frequency")
g6000_vec <- VariableFeatures(g6000_seu)
averages <- AverageExpression(obj, verbose = F) %>% .[[DefaultAssay(obj)]] %>% log1p()
w_mat <- t(as.matrix(averages)[g6000_vec,]) # wgcna requires rows as samples and columns as geneIDs

# optimizing power
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
sft <- pickSoftThreshold(w_mat, powerVector = powers, verbose = 5, networkType = "signed")

# visualization to pick power
library(gridExtra)
a1 <- ggplot(sft$fitIndices, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.9, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft$fitIndices, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
# choose 10 for at345.clean
# ---------------------------------------------------------------
# skip as function end----------------------------------------------------------


# determine optimal power
w_mat <- wgcna_ppl_1(obj) # choose 10 for ad, 12 for at
save(w_mat, file='w_mat.at345.RData')

# network construction
bwnet <- blockwiseModules(w_mat,
                          maxBlockSize = 6000,
                          TOMType = "signed",
                          power = 10,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 0,
                          verbose = 3)

save(bwnet, file='bwnet.at345.RData')
# module eigen genes
# bwnet$MEs
# get number of genes for each module
# table(bwnet$colors)


plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


plotDendroAndColors(bwnet$dendrograms[[1]], bwnet$colors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# pheatmap
mes <- t(bwnet$MEs)
rownames(mes) <- sub("ME", "", rownames(mes))
mes_order_stub<-unique(colnames(mes)[apply(mes,1,which.max)])
pheatmap(mes[, mes_order_stub], cluster_cols = F, cluster_rows = F, scale = 'none', treeheight_row = 0, treeheight_col = 0)

cluster.order<-sort(levels(obj))
avg_max <- colnames(mes)[apply(mes,1,which.max)]
mes_order_stub <- avg_max[order(match(avg_max, cluster.order))]
#############################
order_genes <- function (x, object, cell.type.order, avg.exp.matrix)
{
  x_average_matrix <- avg.exp.matrix$RNA[x,]
  x_average_matrix_max <- colnames(x_average_matrix)[apply(x_average_matrix,1,which.max)]
  names(x_average_matrix_max) <- x
  x_average_matrix_max_ordered <- x_average_matrix_max[order(match(x_average_matrix_max, cell.type.order))]
  return(names(x_average_matrix_max_ordered))
}
##############################



# dotplot
# ranking genes for module membership
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))

#assigned previously
#name_map <- read.csv("aten_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmap(bwnet, me_wt, name_map, 300, 'ad3456')


# for splited object ------------------------------------------
Idents(obj) <- obj$sym_split
w_mat <- wgcna_ppl_1(obj) # choose 12 for at, 10 for ad.split


# network construction
bwnet <- blockwiseModules(w_mat,
                          maxBlockSize = 6000,
                          TOMType = "signed",
                          power = 10,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 0,
                          verbose = 3)

table(bwnet$colors)
table_df <- as.data.frame(table(obj$sym_split))
write.table(table_df, file = "ad3456.split.cell_num_per_cluster.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


plotDendroAndColors(bwnet$dendrograms[[1]], bwnet$colors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# pheatmap
mes <- t(bwnet$MEs)
rownames(mes) <- sub("ME", "", rownames(mes))
### order_genes
mes_order_stub<-unique(colnames(mes)[apply(mes,1,which.max)])
pheatmap(mes[, mes_order_stub], cluster_cols = F, cluster_rows = F, scale = 'none', treeheight_row = 0, treeheight_col = 0)


# dotplot by ranking genes for module membership
mes <- bwnet$MEs
colnames(mes) <- sub("ME", "", colnames(mes))
me_wt = as.data.frame(cor(x = w_mat, y = mes, use="p"))

#name_map <- read.csv("aten_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
wgcna_dotplots_withmap(bwnet, me_wt, name_map, 300, 'ad3456.split')












load("ad3456_nosym.mc.RData")
name_map <- read.csv("adig_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmap(ad3456_nosym.mc, name_map, "ad3456.mc")

# ---------------------------------------------------------------
# extract no sym sample
load("ad3456.std.RData")

ns_genes <- grep("adig-", rownames(ad3456.std))
ad3456.ns <- ad3456.std[ns_genes,]
#ad3456.ns <- std_seurat_ppl(ad3456.ns)

name_map <- read.csv("adig_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmap(ad3456.ns, name_map)

# -------------------------------------------------
load("at345.std.RData")

ns_genes <- grep("aten-", rownames(at345.std))
at345.ns <- at345.std[ns_genes,]
#at345.ns <- std_seurat_ppl(at345.ns)

name_map <- read.csv("aten_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmap(at345.ns, name_map)



####################################################
# distribtuion of symbiodinium related genes

obj<-ad3456.std
gene_list <- read.csv('ad_symbiodinium.list', header=F, sep='\t')[[1]]
gene_list <- grep("^Sym", rownames(obj), value = TRUE, ignore.case = TRUE)
obj$sdlike <- ifelse(rownames(obj) %in% gene_list,"s","n")
DimPlot(obj, reduction = "umap", group.by = "sdlike") + scale_color_manual(values = c("s" = "red", "n" = "lightgray"))

obj<-at345.std
gene_list <- read.csv('at_symbiodinium.list', header=F, sep='\t')[[1]]
gene_list <- grep("^Sym", rownames(obj), value = TRUE, ignore.case = TRUE)
obj$sdlike <- ifelse(rownames(obj) %in% gene_list,"s","n")
DimPlot(obj, reduction = "umap", group.by = "sdlike") + scale_color_manual(values = c("s" = "red", "n" = "lightgray"))
####################################################
# locating gene annotation redundancy
adig_emapper_orthotable <- read_tsv("problem/aten.prot.t1_longestprotein.emapper.orthologs.tsv", skip = 4) # redundant
adig_protein2gene_lookup <- read_tsv("problem/aten_proteinID_to_geneID_lookup.tsv") # non redundant
adig_transcript2gene_lookup <- read_tsv("problem/aten_transcriptID_to_geneID_lookup.tsv") # non redundant

human_genes <- read_tsv("problem/human_proteinIDS_geneNAMES.txt", show_col_types = F)
fly_genes <- read_tsv("problem/fbgn_fbtr_fbpp_expanded_fb_2022_01.tsv", show_col_types = F, skip = 4)
fly_genes <- fly_genes %>% select(polypeptide_ID, gene_ID, gene_symbol, gene_fullname)
goterms_table <- read_tsv("goterms_lookup.txt", col_names = c("go_id", "term"))


adig_human_orthos <- adig_emapper_orthotable %>%
  filter(species == "Homo sapiens(9606)") %>%
  separate_rows(orthologs, sep = ",") %>%
  left_join(human_genes, by = c("orthologs" = "Protein_stable_ID")) 

adig_human_orthos$Gene_name[is.na(adig_human_orthos$Gene_name)] <- adig_human_orthos$orthologs[is.na(adig_human_orthos$Gene_name)]

adig_human_orthos <- adig_human_orthos %>%
  group_by(query) %>%
  mutate(human_ortholog_symbols = paste(Gene_name, collapse = ","), human_proteinIDs = paste(orthologs, collapse = ","), human_geneIDs = paste(Gene_stable_ID, collapse = ",")) %>%
  rename(proteinID = query, human_orth_type = orth_type) %>%
  left_join(adig_protein2gene_lookup) %>%
  ungroup() %>%
  select(proteinID, geneID, human_orth_type, human_proteinIDs, human_geneIDs, human_ortholog_symbols) %>%
  unique()


adig_dm_orthos <- adig_emapper_orthotable %>%
  filter(species == "Drosophila melanogaster(7227)") %>%
  separate_rows(orthologs, sep = ",") %>%
  left_join(fly_genes, by = c("orthologs" = "polypeptide_ID")) 

adig_dm_orthos$Gene_name[is.na(adig_dm_orthos$gene_symbol)] <- adig_dm_orthos$polypeptide_ID[is.na(adig_dm_orthos$gene_symbol)]

adig_dm_orthos <- adig_dm_orthos %>%
  group_by(query) %>%
  mutate(fly_ortholog_symbols = paste(gene_symbol, collapse = ","), 
         fly_ortholog_names = paste(gene_fullname, collapse = ","), 
         fly_proteinIDs = paste(orthologs, collapse = ","), 
         fly_geneIDs = paste(gene_ID, collapse = ",")) %>%
  rename(proteinID = query, fly_orth_type = orth_type) %>%
  left_join(adig_protein2gene_lookup) %>%
  ungroup() %>%
  select(proteinID, geneID, fly_orth_type, fly_proteinIDs, fly_geneIDs, fly_ortholog_symbols, fly_ortholog_names) %>%
  unique()

## it's here!
adig_names <- tibble(geneID = unique(adig_transcript2gene_lookup$geneID)) %>%
  left_join(adig_human_orthos) %>%
  left_join(adig_dm_orthos)

write.table(unique(adig_transcript2gene_lookup$geneID), "ugid.tsv", row.names=F, col.names=F, sep="\t", quote=F)
write.table(adig_human_orthos, "hs.tsv", row.names=F, col.names=F, sep="\t", quote=F)
write.table(adig_dm_orthos, "dm.tsv", row.names=F, col.names=F, sep="\t", quote=F)
write.table(adig_names, "ret.tsv", row.names=F, col.names=F, sep="\t", quote=F)


## kj's demo
t_names <- tibble(
  geneID = c("geneA", "geneB", "geneC")
)
t_hs <- tibble(
  geneID = c("geneA", "geneA", "geneB", "geneD"),
  human_ortholog = c("hGene1", "hGene2", "hGene3", "hGene4")
)
# ret <- t_names %>% left_join(t_hs, by = "geneID")
ret <- t_names %>% left_join(t_hs)

#> ret
# A tibble: 4 × 2
#geneID human_ortholog
#<chr>  <chr>         
#  1 geneA  hGene1        
#2 geneA  hGene2        
#3 geneB  hGene3        
#4 geneC  NA   



adig_names <- adig_names %>%
  add_column(final_emapper_name = adig_names$geneID)

adig_names$final_emapper_name[!is.na(adig_names$fly_ortholog_symbols)] <- paste0(adig_names$geneID[!is.na(adig_names$fly_ortholog_symbols)],
                                                                                 " ",
                                                                                 adig_names$fly_orth_type[!is.na(adig_names$fly_ortholog_symbols)],
                                                                                 " dm-",
                                                                                 adig_names$fly_ortholog_symbols[!is.na(adig_names$fly_ortholog_symbols)])

adig_names$final_emapper_name[!is.na(adig_names$human_ortholog_symbols)] <- paste0(adig_names$geneID[!is.na(adig_names$human_ortholog_symbols)],
                                                                                   " ",
                                                                                   adig_names$human_orth_type[!is.na(adig_names$human_ortholog_symbols)],
                                                                                   " hs-",
                                                                                   adig_names$human_ortholog_symbols[!is.na(adig_names$human_ortholog_symbols)])


#########################################################
samples<-c("ad3c","ad4c","ad5c", "ad6c")
samples<-c("at3c","at4c","at5c")


symbiont_gene_table <- tibble(sample = samples) %>%
  mutate(sampleID = gsub("_trin", "", sample)) %>%
  rowwise() %>%
  mutate(num_sym_genes = length(grep("symb", rownames(get(sample)))) + length(grep("comp", rownames(get(sample))))) %>%
  mutate(num_sym_cells = count_symbiont_cells(object = get(sample))) %>%
  mutate(num_sym_umis = count_symbiont_umis(object = get(sample))) %>%
  ungroup()

symbiont_gene_table

ggplot(symbiont_gene_table, aes(x = num_sym_umis, fill = sampleID)) + geom_histogram()
ggplot(symbiont_gene_table, aes(x = num_sym_genes, fill = sampleID)) + geom_histogram()
ggplot(symbiont_gene_table, aes(x = num_sym_cells, fill = sampleID)) +  geom_histogram()


#################################################
ad3456.std$symbiont_umis <- count_symbiont_umis_per_cell(ad3456.std) ## return a vector len = number of cells
ad3456.std$symbiont_genes <- count_symbiont_genes_per_cell(ad3456.std)

adig_genes <- grep("adig-", rownames(ad3456.std))
ad3456_nosym <- ad3456.std[adig_genes,]
ad3456_nosym <- std_seurat_ppl(ad3456_nosym)
ad3456_nosym[['clusters_nosym_pcs40_res2']] <- Idents(ad3456_nosym)

dp <- DimPlot(ad3456_nosym, reduction = "umap", label = T, repel = T, group.by = 'clusters_nosym_pcs40_res2')
ggsave(dp, filename = "ad3456_nosym_umap_sample.pdf", units = "in", width = 8, height = 6)

ad3456_nosym.mc <- merge_clusters3j(ad3456_nosym)

# ad3456_nosym.mc$symbiont_umis <- ad3456.std$symbiont_umis
# ad3456_nosym.mc$symbiont_genes <- ad3456.std$symbiont_genes


# Count the number of cells starting with "at"
apo <- sum(grepl("^ad3", colnames(ad3456_nosym.mc)))+sum(grepl("^ad4", colnames(ad3456_nosym.mc)))
sym <- sum(grepl("^ad5", colnames(ad3456_nosym.mc)))+sum(grepl("^ad6", colnames(ad3456_nosym.mc)))

ad3456_nosym.mc$aposym <- c(rep("apo", apo), rep("sym", sym))
# problem!
# ad3456_nosym.mc$aposym <- c(rep("apo", apo), rep("sym", sym))
# VlnPlot(ad3456_nosym.mc, "symbiont_umis", split.by = "aposym", split.plot = T)

ad3456_nosym.mc_sym_tib <- tibble(ad3456_nosym_cluster = as.character(Idents(ad3456_nosym.mc)), sym = ad3456_nosym.mc@meta.data$symbiont_umis > 0) %>%
  group_by(ad3456_nosym_cluster) %>%
  summarize(sym_num = sum(sym == TRUE)) %>%
  ungroup() %>%
  mutate(sym_prop = sym_num/sum(sym_num))

ad3456_nosym.mc_sym_tib


ad3456_nosym.mc_sym_sum <- tibble(ad3456_nosym_cluster = as.character(levels(ad3456_nosym.mc)), cell_prop = (as.numeric(table(Idents(ad3456_nosym.mc))) / ncol(ad3456_nosym.mc))) %>%
  left_join(ad3456_nosym.mc_sym_tib) %>%
  select(ad3456_nosym_cluster, cell_prop, sym_prop) %>%
  gather(key = "type", value = "proportion", 2:3)


# ad3456_nosym.mc_sym_sum$ad3456_nosym_cluster <- ordered(ad3456_nosym.mc_sym_sum$ad3456_nosym_cluster)



gg <- ggplot(ad3456_nosym.mc_sym_sum, aes(x = ad3456_nosym_cluster, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position=position_dodge(0.6), width = .6) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(size = 6))

gg

#ggsave2(plot = gg, filename = paste0(output_directory, "symbiont_enrichment_plots/ad3456_nosym_cell_and_symbiont_proportions_barchart.pdf"), units = "in", width = 10, height = 6


# pending
ad3456_nosym_avgexp <- AverageExpression(ad3456_nosym.mc, return.seurat = F)
nj_matrix <- t(log1p(ad3456_nosym_avgexp$RNA))
nj_matrix_dist <- as.dist(cor(as.matrix(nj_matrix)))
nj_matrix_tree <- nj(nj_matrix_dist)

plot.phylo(midpoint.root(nj_matrix_tree))
plot.phylo(nj_matrix_tree)



##########################################

# at345.std$symbiont_umis <- count_symbiont_umis_per_cell(at345.std) ## return a vector len = number of cells
# at345.std$symbiont_genes <- count_symbiont_genes_per_cell(at345.std)

at_genes <- grep("aten-", rownames(at345.std))
at345_nosym <- at345.std[at_genes,]
at345_nosym <- std_seurat_ppl(at345_nosym)

at345_nosym[['clusters_nosym_pcs40_res2']] <- Idents(at345_nosym)

dp <- DimPlot(at345_nosym, reduction = "umap", label = T, repel = T, group.by = 'clusters_nosym_pcs40_res2')
ggsave(dp, filename = "at345_nosym_umap_sample.pdf", units = "in", width = 8, height = 6)

at345_nosym.mc <- merge_clusters3j(at345_nosym)

# at345_nosym.mc$symbiont_umis <- at345.std$symbiont_umis
# at345_nosym.mc$symbiont_genes <- at345.std$symbiont_genes


# Count the number of cells starting with "at"
# apo <- sum(grepl("^at3", colnames(at345_nosym.mc)))
# sym <- sum(grepl("^at4", colnames(at345_nosym.mc)))+sum(grepl("^at5", colnames(at345_nosym.mc)))
# at345_nosym.mc$aposym <- c(rep("apo", apo), rep("sym", sym))


at345_nosym.mc_sym_tib <- tibble(at345_nosym_cluster = as.character(Idents(at345_nosym.mc)), sym = at345_nosym.mc@meta.data$symbiont_umis > 0) %>%
  group_by(at345_nosym_cluster) %>%
  summarize(sym_num = sum(sym == TRUE)) %>%
  ungroup() %>%
  mutate(sym_prop = sym_num/sum(sym_num))

ad3456_nosym.mc_sym_tib


at345_nosym.mc_sym_sum <- tibble(at345_nosym_cluster = as.character(levels(at345_nosym.mc)), cell_prop = (as.numeric(table(Idents(at345_nosym.mc))) / ncol(at345_nosym.mc))) %>%
  left_join(at345_nosym.mc_sym_tib) %>%
  select(at345_nosym_cluster, cell_prop, sym_prop) %>%
  gather(key = "type", value = "proportion", 2:3)


gg <- ggplot(at345_nosym.mc_sym_sum, aes(x = at345_nosym_cluster, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position=position_dodge(0.6), width = .6) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))

gg

# output: at345_nosym.mc
##########################################
setRepositories()
renv::install("impute")

renv::install("BiocManager")
BiocManager::install("GO.db")

renv::install("WGCNA")
renv::install("cowplot")
###########################################
# dotplot
generate_cluster_dotplots(ad3456_nosym.mc)
generate_cluster_dotplots(at345_nosym.mc)

# violin plot
VlnPlot(at345_nosym.mc, features = c("nCount_RNA", "nFeature_RNA"), group.by = "merged_clusters", pt.size = 0.1, ncol = 1)
VlnPlot(ad3456_nosym.mc, features = c("nCount_RNA", "nFeature_RNA"), group.by = "merged_clusters", pt.size = 0.1, ncol = 1)


