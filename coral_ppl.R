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

# renv::install('devtools')
# renv::install("C:\\Users\\kjia\\workspace\\library\\repository\\presto")
# renv::install("C:\\Users\\kjia\\workspace\\library\\repository\\DoubletFinder")

# x: list of genes
# cell.type.order = levels(arcs)
# avg.exp.matrix
# order_genes <- function (x, cell.type.order, avg.exp.matrix){
# 	x_average_matrix <- avg.exp.matrix[, x]
# 	x_average_matrix_max <- colnames(x_average_matrix)[apply(x_average_matrix,1,which.max)] #finds the index of the maximum value in each row
# 	names(x_average_matrix_max) <- x
# 	x_average_matrix_max_ordered <- x_average_matrix_max[order(match(x_average_matrix_max, cell.type.order))]
# 	return(names(x_average_matrix_max_ordered))
# }
# mes_order_stub<-colnames(mes)[apply(mes,1,which.max)]


# obj: seurat object (after removing doublets, ribosomal genes, junk clusters)
# output: visualization of optimal power. User it in wgcna_ppl_2()
# w_mat <- t(as.matrix(averages)[g6000_vec,]) 
wgcna_ppl_1 <- function(obj){
	# wgcna normal ppl
	g6000_seu <- FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = 6000, binning.method = "equal_frequency")
	g6000_vec <- VariableFeatures(g6000_seu)
	averages <- AverageExpression(obj, verbose = F) %>% .[[DefaultAssay(obj)]] %>% log1p()
	w_mat <- t(as.matrix(averages)[g6000_vec,]) # wgcna requires rows as samples and columns as geneIDs

	# optimizing power
	powers <- c(c(1:10), seq(from = 12, to=30, by=2))
	sft <- pickSoftThreshold(w_mat, powerVector = powers, verbose = 5, networkType = "signed")

	# visualization to pick power
	#library(gridExtra)
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
	
	return(w_mat)
}


# name_map <- read.csv("aten_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
# bwnet: output of wgcna blockwiseModules() 
# me_wt: correlation matrix between module_color and (6000) genes
wgcna_dotplots_withmaps <- function(obj, bwnet, me_wt, name_map, xn_map, max_gene_num, x_ordered, outprefix){
	count=0
	if(setequal(levels(obj),x_ordered)==FALSE){
		message(glue("ordered label is different from levels(obj)"))
		return(FALSE)
	}	
	obj@active.ident <- factor(Idents(obj), levels = x_ordered)
	module_names <- unique(bwnet$colors)
	for(i in module_names){
		#genes_to_plot <- me_wt[row.names(me_wt)[bwnet$colors == "floralwhite"], "floralwhite", drop = F]
		#genes_to_plot_ordered <- row.names(genes_to_plot)[order(-genes_to_plot$floralwhite)]

		genes_to_plot <- me_wt[row.names(me_wt)[bwnet$colors == i], i, drop = F]
		genes_to_plot_ordered <- row.names(genes_to_plot)[order(-genes_to_plot[[i]])]
		
		if(length(genes_to_plot_ordered) > max_gene_num){
			message(glue("Trim gene number to {max_gene_num}."))
			genes_to_plot_ordered <- genes_to_plot_ordered[1:max_gene_num]
		}
		
		if(length(genes_to_plot_ordered) > 0){
		  #Idents(obj) <- factor(Idents(obj), levels = sort(levels(obj)))
		  gg <- Seurat::DotPlot(object = obj, features = rev(genes_to_plot_ordered), cols = "RdYlBu") + coord_flip() +
			theme(axis.text.x = element_text(hjust = 1, angle = 75))

		  gb<-ggplot_build(gg)

		  # alter gene(y) labels
		  ystr <- gb$layout$panel_params[[1]]$y$get_labels()
		  ystr <- gsub("-","_",ystr)
		  ydf <-data.frame(geneID=ystr)
		  mapped_names <- left_join(ydf, name_map, by = "geneID") %>% distinct(geneID, .keep_all=T)

		  # alter cluster(x) labels
		  xstr <- gb$layout$panel_params[[1]]$x$get_labels()
		  xdf <-data.frame(idx=xstr)
		  #xn_map<- read.csv("out", header = TRUE, sep = "\t")
		  x_names <- left_join(xdf, xn_map, by = "idx") %>% distinct(idx, .keep_all=T)

		  #mapped_names$final_emapper_name
		  # gg<-gg+scale_x_discrete(labels = mapped_names$final_emapper_name)+ ggtitle(i)
		  gg<-gg + scale_x_discrete(labels = mapped_names$gene_alias) + scale_y_discrete(labels=x_names$alias) + ggtitle(i)
			
		  outname <- paste0("output/", outprefix, "_m_", i, "_wgcna_dotplot.pdf")

		  ggsave(filename = outname,
				 plot = gg, 
				 path = ".",
				 width = 40, #30, 
				 height = (length(genes_to_plot_ordered) * (1/6)) + 20, # +2, 
				 units = "in", 
				 limitsize = F)
		  count=count+1
		  message(glue("Save to {outname}."))
		}
	}
	ncluster = length(module_names)
	message(glue("{count} / {ncluster} dotplot saved."))
}


# order nj tips according to tree structure
ordered_tips_apenj<-function(tree){
	is_tip <- tree$edge[,2] <= length(tree$tip.label)
	return(tree$tip.label[tree$edge[is_tip, 2]])
}

# yn_map: gene alias
# xn_map: cluster alias
# x_ordered: cluster ordered by neighbor joining; output from function ordered_tips_apenj()
generate_cluster_dotplots_withmaps <- function(obj, yn_map, xn_map, x_ordered, outprefix){
	count=0
	# check cluster IDs in ordered list is the same set of IDs from levels(obj)
	if(setequal(levels(obj),x_ordered)==FALSE){
		message(glue("ordered label is different from levels(obj)"))
		return(FALSE)
	}
	# for(i in levels(obj)){
	for(i in x_ordered){
		cluster_markers_tibble <- as_tibble(rownames_to_column(FindMarkers(obj, ident.1=i, only.pos = T), var = "gene"))
		
		genes_for_dotplot <- cluster_markers_tibble %>%
			filter(avg_log2FC > 0, p_val_adj < 0.05) %>% ## filter rows
			dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>% ## sort filter results first by p_val_adj ascending, then avg_log2FC decending
			dplyr::slice_head(n = 150) %>%  ## select top 150 rows
			pull(gene) ## get gene column value of previous result as a vector
		
		
		if(length(genes_for_dotplot) > 0){
			gg <- Seurat::DotPlot(object = obj, features = rev(genes_for_dotplot), cols = "RdYlBu") + coord_flip() +
			theme(axis.text.x = element_text(hjust = 1, angle = 75))
			
			gb<-ggplot_build(gg)
			# alter gene(y) labels
			ystr <- gb$layout$panel_params[[1]]$y$get_labels()
			ystr <- gsub("-","_",ystr)
			ydf <-data.frame(geneID=ystr)
			y_names <- left_join(ydf, yn_map, by = "geneID") %>% distinct(geneID, .keep_all=T)
			
			# alter cluster(x) labels
			xstr <- gb$layout$panel_params[[1]]$x$get_labels()
			xdf <-data.frame(idx=xstr)
			x_names <- left_join(xdf, xn_map, by = "idx") %>% distinct(idx, .keep_all=T)
			
			#y_names$final_emapper_name
			gg<-gg + scale_x_discrete(labels = y_names$gene_alias) + scale_y_discrete(labels=x_names$alias) + ggtitle(sprintf("%s (%d)", i, sum(Idents(ad3456.clean) == i)))
			
			outname <- paste0("output/", outprefix, "_c_", i, "_dotplot.pdf")
			
			ggsave(filename = outname,
				plot = gg, 
				path = ".",
				width = 30, # 18
				height = (length(genes_for_dotplot) * (1/6)) + 20, # +2
				units = "in", 
				limitsize = F)
			count=count+1
			message(glue("Save to {outname}."))
		}
	}
	ncluster = length(levels(obj))
	message(glue("{count} / {ncluster} dotplot saved."))
}

# name_map <- read.csv("aten_emapper_names_table_final.tsv", header = TRUE, sep = "\t")
generate_cluster_dotplots_withmap <- function(obj, name_map, outprefix){
	count=0
	for(i in levels(obj)){
		cluster_markers_tibble <- as_tibble(rownames_to_column(FindMarkers(obj, ident.1=i, only.pos = T), var = "gene"))

		genes_for_dotplot <- cluster_markers_tibble %>%
		  filter(avg_log2FC > 0, p_val_adj < 0.05) %>% ## filter rows
		  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>% ## sort filter results first by p_val_adj ascending, then avg_log2FC decending
		  dplyr::slice_head(n = 150) %>%  ## select top 150 rows
		  pull(gene) ## get gene column value of previous result as a vector


		if(length(genes_for_dotplot) > 0){
		  gg <- Seurat::DotPlot(object = obj, features = rev(genes_for_dotplot), cols = "RdYlBu") + coord_flip() +
			theme(axis.text.x = element_text(hjust = 1, angle = 45))

		  gb<-ggplot_build(gg)
		  ystr <- gb$layout$panel_params[[1]]$y$get_labels()
		  ystr <- gsub("-","_",ystr)
		  ydf <-data.frame(geneID=ystr)
		  mapped_names <- left_join(ydf, name_map, by = "geneID") %>% distinct(geneID, .keep_all=T)
		  #mapped_names$final_emapper_name
		  gg<-gg+scale_x_discrete(labels = mapped_names$final_emapper_name)
			
		  #outname <- paste0("dotplots/", outprefix, "_c_", i, "_dotplot.pdf")
		  outname <- paste0("output/", outprefix, "_c_", i, "_dotplot.jpg")

		  ggsave(filename = outname,
				 plot = gg, 
				 path = ".",
				 width = 18, 
				 height = (length(genes_for_dotplot) * (1/6)) + 2, 
				 units = "in", 
				 limitsize = F)
		  count=count+1
		  message(glue("Save to {outname}."))
		}
	}
	ncluster = length(levels(obj))
	message(glue("{count} / {ncluster} dotplot saved."))
}



generate_cluster_dotplots <- function(obj, outprefix){
	count=0
	for(i in levels(obj)){
		cluster_markers_tibble <- as_tibble(rownames_to_column(FindMarkers(obj, ident.1=i, only.pos = T), var = "gene"))

		genes_for_dotplot <- cluster_markers_tibble %>%
		  filter(avg_log2FC > 0, p_val_adj < 0.05) %>% ## filter rows
		  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>% ## sort filter results first by p_val_adj ascending, then avg_log2FC decending
		  dplyr::slice_head(n = 150) %>%  ## select top 150 rows
		  pull(gene) ## get gene column value of previous result as a vector


		if(length(genes_for_dotplot) > 0){
		  gg <- Seurat::DotPlot(object = obj, features = rev(genes_for_dotplot), cols = "RdYlBu") + coord_flip() +
			theme(axis.text.x = element_text(hjust = 1, angle = 45))
		  
		  outname <- paste0("dotplots/", outprefix, "_c_", i, "_dotplot.pdf")
		  ggsave(filename = outname,
				 plot = gg, 
				 path = ".",
				 width = 18, 
				 height = (length(genes_for_dotplot) * (1/6)) + 2, 
				 units = "in", 
				 limitsize = F)
		  count=count+1
		  message(glue("Save to {outname}."))
		}
	}
	ncluster = length(levels(obj))
	message(glue("{count} / {ncluster} dotplot saved."))
}


top_expr_genes <- function(obj, ntop){
	ids <- paste0("g", as.character(Idents(obj)))
	avg_exp <- AverageExpression(obj, return.seurat = FALSE)

	top_genes <- list()
	for (cluster in unique(ids)) {
	  sorted_genes <- sort(avg_exp$RNA[, cluster], decreasing = TRUE)
	  top_genes[[cluster]] <- names(sorted_genes)[1:ntop]
	}
	top_genes_flat<-unique(unlist(top_genes))
}


# wrapper for loading std doublet ppl
fn_prep <- function(data_path, name, r){
	d <- Read10X(data.dir = data_path, gene.column=1)
	d <- CreateSeuratObject(counts=d, project=name, min.cells=3, min.features=0)
	d <- std_seurat_ppl(d)
	d <- find_double(d, ratio=r)
return(d)
}

# wrapper for standard seurat ppl
# NormalizeData: calculate 1. scale each cell expression sum to be 10,000: 10,000 * count/sum. 2. take log(s+1); all positive decimal
# Log-Transformation: Applies a log transformation to stabilize the variance and make the data more normally distributed.
# FindVariableFeatures: apply on normalized data. calculate the difference between expected variance and observed variance. expected variance is calculated using regression between mean expression and mean variance of each gene.
# ScaleData: transform normalized data into z scores; has negative values

std_seurat_ppl <- function(obj){
	obj <- NormalizeData(obj , normalization.method = "LogNormalize", scale.factor = 10000) %>%
	FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
	ScaleData(verbose = F, vars.to.regress = "nCount_RNA") %>%
	RunPCA(verbose = F) %>%
	FindNeighbors(dims = 1:40, k.param = 20, verbose = FALSE, nn.method = "annoy", annoy.metric = "cosine") 
	message("Find clusters")
	obj <- FindClusters(obj , resolution = 2, algorithm = 4, verbose = FALSE)  
	message("skip UMap, run it seperately if needed")
	# obj <- RunUMAP(obj, dims = 1:40, min.dist = 0.1, umap.method = "umap-learn", metric = "cosine", verbose = FALSE)
	return(obj)
}

# wapper for running DoubleFinder
# ratio table: https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled
# Multiplet Rate	# of Cells Loaded # of Cells Recovered
# 0.004	825	500
# 0.008	1650	1000
# 0.016	3300	2000
# 0.024	4950	3000
# 0.032	6600	4000
# 0.040	8250	5000
# 0.048	9900	6000
# 0.056	11550	7000
# 0.064	13200	8000
# 0.072	14850	9000
# 0.080	16500	10000
find_double <- function(obj, pcs=1:20, pN = 0.25, ratio=0.075){
	sweep.obj<-paramSweep(obj, PCs = pcs, sct = FALSE)
	sweep.stat <- summarizeSweep(sweep.obj, GT=FALSE)
	bc.metric <- find.pK(sweep.stat)
	pK <- bc.metric %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
	pK <- as.numeric(as.character(pK[[1]]))
	message(glue("Find optimal pK: {pK}, identifying with ratio: {ratio}"))
	# to visualize pK values
	# ggplot(bc.metric, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()
	
	homotypic.prop <- modelHomotypic(obj@meta.data$seurat_clusters)
	nExp_poi <- round(ratio*nrow(obj@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	
	obj <- doubletFinder(obj, PCs = pcs, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
	return(obj)
}

# Jake's merging algorithm
# need to be refactored
# source("my_merge_clusters.R")
# obj<-at34
# Idents(obj) <- paste0("g", as.character(Idents(obj)))
# obj <- merge_clusters3j(obj)
merge_clusters3j <- function(obj, thresh_corr = 0.7, num.genes = 20, fold.difference = 2, adj.p.value = 0.05, assay = DefaultAssay(obj), genes = VariableFeatures(obj), min.cells = 25) {
	#Idents(obj) <- paste0("g", as.character(Idents(obj)))
	message(glue("Merging clusters with {min.cells} or fewer cells into nearest neighbor cluster with {thresh_corr} correlation threshold."))
	# first loop lumps clusters with 10 or fewer cells into nearest neighbor clusters (only if necessary)
	# while(min(table(obj@active.ident)) <= min.cells){
	count=0
	done_merging = FALSE
	while(!done_merging){
		merged_clusters <- as.character(obj@active.ident)
		averages <- AverageExpression(obj, verbose = F, assays = assay) %>% .[[assay]] %>% log1p()
		#averages <- log1p(averages[[assay]][genes,])
		
		corrs <- cor(as.matrix(averages))
		corrs[lower.tri(corrs,diag=TRUE)] <- NA  #convert lower triangle and diagonal to NAs
		corrs <- as.data.frame(as.table(corrs))  #convert to 3-column table
		corrs <- na.omit(corrs)  #Get rid of the NAs from lower triangle (duplicates)
		corrs <- corrs %>% filter(Freq > thresh_corr)
		corrs <- corrs[order(-abs(corrs$Freq)),]  #Sort by highest correlation (whether +ve or -ve)	

		## checking all candidates
		done_merging = TRUE
		for(i in unique(c(corrs$Var1,corrs$Var2))){
			n_cells <-length(which(obj@active.ident == i))
			if( n_cells <= min.cells) {
				message(glue("checking cluster {i}"))
				# found small cluster; calcualte correlations
				## find top correlated partner
				index <- which(corrs == i, arr.ind = T)[1,] ## find all rows contain "i" and (arr.ind) return the first index;## jk: now it's "g"+"i"				
				## index contains a list of (row_id, column_id)
				## index[1]: the first occurrence of "i" 
				closest_cluster <- corrs[index[1],] ## get the top correlated cluster pair 
			
				new_cluster = paste(as.character(closest_cluster[1,1]),	as.character(closest_cluster[1,2]),	sep = "-")
				merged_clusters[which(merged_clusters == closest_cluster[1,1])] <- new_cluster
				merged_clusters[which(merged_clusters == closest_cluster[1,2])] <- new_cluster
				
				message(glue("Merging clusters {i}[{n_cells}] to {new_cluster} corr: {round(closest_cluster[1,3],2)} "))
				count = count + 1
				done_merging=FALSE
				break # merge one pair per loop
			}
		}
		obj$merged_clusters <- merged_clusters
		Idents(obj) <- merged_clusters	
	}
	message(glue("Finish merging {count} small clusters.\n\n"))

	
	#LOOP 2 iteratively merges clusters that do not satisfy threshold criteria
	message(glue("Merging nearest-neighbor clusters with less than {as.character(num.genes)} genes with log fold difference of {as.character(fold.difference)}"))
	keep_merging = TRUE
	tested <- list() ## cache FindMarkers() results to avoid re-calculating	
	while (keep_merging){
		merged_clusters <- as.character(obj@active.ident)
		averages <- AverageExpression(obj, verbose = F, assays = assay) %>% .[[assay]] %>% log1p()
		#averages <- log1p(averages[[assay]][genes,])
		
		corrs <- cor(as.matrix(averages))
		corrs[lower.tri(corrs,diag=TRUE)] <- NA  #convert lower triangle and diagonal to NAs
		corrs <- as.data.frame(as.table(corrs))  #Turn into a 3-column table
		corrs=na.omit(corrs)  #Get rid of the junk we flagged above
		# corrs <- corrs %>% filter(Freq > thresh_corr)
		corrs <- corrs %>% filter(Freq > 0.8)
		corrs <- corrs[order(-abs(corrs$Freq)),]    #Sort by highest correlation (whether +ve or -ve)
		
		# examine criteria for merging
		progress=0
		keep_merging = FALSE
		cat(glue("searching {nrow(corrs)} pairs: "))
		for (i in 1:nrow(corrs)){
			evl_key <- paste(corrs$Var1[i], corrs$Var2[i])
			if(!(evl_key %in% tested)){
				de_m <- FindMarkers(obj = obj, ident.1 = corrs[i,1], ident.2 = corrs[i,2], print.bar = F, min.pct = 0.25, assay = assay, verbose = F)
				de_m <- tibble::rownames_to_column(de_m, var = "gene")
				de_m_filtered <- de_m %>% filter(abs(avg_log2FC) > log2(fold.difference), p_val_adj < adj.p.value)
				
				progress=progress+1
				#cat(glue("..{round(progress/nrow(corrs), 2)*100}%"))
				message(glue("[{round(progress/nrow(corrs), 2)*100}%]  Evaluating clusters {corrs[i,1]} and {corrs[i,2]} , DE genes = {nrow(de_m_filtered)}"))
				tested <- c(tested, evl_key)
				if((nrow(de_m_filtered)) < num.genes){
					new_cluster = paste(corrs[i,1], corrs[i,2], sep = ".")
					message("\n", glue("merging {corrs[i,1]} + {corrs[i,2]} corr: {round(corrs[i,3],2)}\n\n"))
					merged_clusters[which(merged_clusters == corrs[i,1])] <- new_cluster
					merged_clusters[which(merged_clusters == corrs[i,2])] <- new_cluster
					keep_merging = TRUE
					count=count+1
					break # merge one pair per loop
				}
			}
		}
		obj$merged_clusters <- merged_clusters
		Idents(obj) <- merged_clusters
	}
	message(glue("\n\nFinished merging {count} clusters by DE."))
  return(obj)
}



# redo Arun's algorithm
# source("my_merge_clusters.R")
# .obj<-at34
# Idents(.obj) <- paste0("g", as.character(Idents(.obj)))
# .obj<-merge_clusters(.obj, .thresh_n_de = 30)
merge_clusters <- function(.obj, .merged.idents.col="merged.idents", .thresh_corr = 0.8, .thresh_n_de = 30, .thresh_log2FC = 1) {
	time_start <- Sys.time()
	counter <- 1
	vf <- VariableFeatures(.obj)
	repeat {
		message(glue("iteration {counter}"))
		n_clusters_old <- length(unique(Idents(.obj)))
		.obj[[.merged.idents.col]] <- Idents(.obj)

		# calculate cluster expression measures {average expression}
		avg_exp <- AverageExpression(.obj, 
                             assays = DefaultAssay(.obj), 
                             slot = "data", ## by default
                             #features = VariableFeatures(.obj), 
							 features = vf, 
                             return.seurat = FALSE, 
                             verbose = FALSE) %>% # return a R list obj (similar to dictionary in python)
					.[[DefaultAssay(.obj)]] %>% # return "RNA" ## meaning get the "RNA" field from the returned list obj, which is a matrix with averaged expression 3000 genes x 52 clusters
					log1p()
					
		corrs <- cor(as.matrix(avg_exp))
		corrs[lower.tri(corrs, diag = TRUE)] <- NA
		corrs <- as.data.frame(as.table(corrs)) # convert to key1 key2 value flat format (Var1 Var2 Freq)
		corrs <- na.omit(corrs) # remove NA
		corrs <- dplyr::rename(corrs, corr = Freq, clust1 = Var1, clust2 = Var2) # rename header names
		corrs <- corrs %>% dplyr::filter(corr > .thresh_corr) %>% arrange(desc(corr)) 
		
		# iterate top correlated group pairs and merge one pair at a time 
		for(i in seq_len(nrow(corrs))){
			c1 <- as.character(corrs$clust1[i])
			c2 <- as.character(corrs$clust2[i])
			de_genes <- FindMarkers(.obj, ident.1 = c1, ident.2 = c2, assay = "RNA", features = VariableFeatures(.obj), base = 2, verbose = FALSE) %>% 
						filter(p_val_adj < 0.05 , abs(avg_log2FC) > .thresh_log2FC) # p_val_adj, avg_log2FC returned from FindMarkers()
						
			n_de <- nrow(de_genes)
			message(glue("  Evaluating clusters {c1} and {c2} ", "(corr = {round(corrs$corr[i], 2)}, DE genes = {n_de})"))
			if(n_de < .thresh_n_de){
				cn <- paste0(c1, ".", c2)
				current.idents <- as.character(Idents(.obj))
				current.idents[current.idents==c1] <- cn
				current.idents[current.idents==c2] <- cn
				.obj[[.merged.idents.col]] <- factor(current.idents, levels = sort(unique(current.idents)))
				.obj <- SetIdent(.obj, value = .merged.idents.col)	# value can be a slot name	
				
				message(glue("  Merged clusters {c1} and {c2} ", "(corr = {round(corrs$corr[i], 2)}, DE genes = {n_de})"))	
				break
			}
		}	
	
		# check for stop criteria
		n_clusters_new <- length(unique(Idents(.obj)))
		if(isTRUE(all.equal(n_clusters_old, n_clusters_new))){
			message("  No clusters merged")
			break
		} else counter <- counter + 1	
		
	} # end of repeat
	message("Finished merging clusters.")


	# store merged identities in new slot with renaming as size-sorted identies
	# .idents.new
	# clusters_by_freq <- names(sort(table(.obj[[.merged.idents.col, drop = TRUE]]), decreasing = TRUE))
	## drop=TRUE: return a single column vector
	# .obj[[.idents.new]] <- factor(.obj[[.merged.idents.col, drop = TRUE]], levels = clusters_by_freq) %>% as.numeric() %>% as.factor()
	# .obj <- SetIdent(.obj, value = .idents.new)
	# if (.remove.old) .obj[[.merged.idents.col]] <- NULL
	# message(glue("Renamed {n_clusters} clusters by size.", n_clusters = length(unique(.obj[[.idents.new, drop = TRUE]]))))
	
	time_elapsed <- Sys.time() - time_start
	message(glue("Elapsed time: {round(time_elapsed[[1]], 2)} ", "{attr(time_elapsed, \"units\")}."))
	
	return(.obj)	
}