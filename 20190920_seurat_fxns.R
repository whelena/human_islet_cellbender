# author: Helena Winata
# 20th Sept 2019
# seurat functions used to analyse human islet data

library(ggplot2)
library(dplyr)
library(sctransform)
library(future)
library(Seurat)
library(cowplot)
library(stringr)

###############################################################################
########################## CREATE SEURAT OBJECT =====########################## 

create_so_list <- function(scater.path, add.meta, sample.id, min.cells = 3, min.feat = 1000, 
                           save.rds = NULL, save.qc = NULL, SCT=F, regr.mt = T, regr.rp = T){
  # sample.id is a list of labels for each file in scater.path. This function take scater objects located in
  # scater.path and create seurat objects for each of them. Option to regress out mitochondrial and ribosomal
  # genes can be selected usin regr.mt and regr.rp. Count data is normalized using standard or SCTransform 
  # workflow. Returns a list of seurat objects.
  scater_objects <- list.files(scater.path, full.names = T)
  sample_list <- vector("list", length = length(sample.id))
  
  for (i in 1: length(scater_objects)){
    data <- readRDS(scater_objects[i])
    so_obj <- CreateSeuratObject(counts = SummarizedExperiment::assay(data,"counts"), min.cells = min.cells, min.features = min.feat)
    so_obj@meta.data$sample_id <- matrix(sample.id[i], nrow = ncol(so_obj))
    
    # add metadata from scater object to seurat object
    if(length(add.meta)>0){
      for(j in 1:length(add.meta)){
        so_obj@meta.data[, add.meta[j]]<-as.factor(colData(data)[colnames(so_obj),add.meta[j]])
        colnames(so_obj@meta.data)[ncol(so_obj@meta.data)]<-add.meta[j]
      }
      so_obj@meta.data$condition_id <- gsub("^R..._", "", so_obj@meta.data$condition)
    }
    
    # get percentage of mitochondrial and ribosomla subunit genes
    if (regr.mt==T){so_obj <- PercentageFeatureSet(so_obj, pattern = "^MT-", col.name = "percent.mt")}
    if (regr.rp==T){so_obj <- PercentageFeatureSet(so_obj, pattern = "^RP(S|L)", col.name = "percent.rp")}
    
    if (SCT==T){
      so_obj <- SCTransform(so_obj, vars.to.regress = c("percent.mt", "percent.rp"))
    }
    else{
      so_obj <- NormalizeData(so_obj)
      so_obj <- FindVariableFeatures(so_obj, selection.method = 'vst', nfeatures = 2000)
    }
    
    sample_list[[i]]<- so_obj
    
    if(is.null(save.qc) == F){
      # quality control plots - optional
      vln <- VlnPlot(so_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
                     ncol = 2, pt.size = 0.1)
      ggsave(plot = vln, file=paste0(save.qc, "qc_vln_plot_", sample.id[i], ".pdf"), width = 8.5, height = 11)
      
      plot1 <- FeatureScatter(so_obj, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1) + NoLegend()
      plot2 <- FeatureScatter(so_obj, feature1 = "nCount_RNA", feature2 = "percent.rp", pt.size = 0.1) + NoLegend()
      plot3 <- FeatureScatter(so_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1) + NoLegend()
      all_plt <- CombinePlots(plots = list(plot1, plot2, plot3))
      ggsave(plot = all_plt, file=paste0(save.qc, "qc_scatter_plot_", sample.id[i], ".pdf"), width = 8.5, height = 11)
    }
    
  }
  if(is.null(save.rds) == F){
    saveRDS(sample_list, file = paste0(save.rds, "_so_list.rds"))
  }
  return(sample_list) 
}
  
  
create_so_int <- function(so_list, save.as = getwd(), SCT=F){
  # integrates all the seurat objects listed in so_list following a standard or SCT workflow. 
  # Returns a single integrated seurat object
  if (SCT==T){
    sample_feat <- SelectIntegrationFeatures(object.list = so_list, nfeatures = 3000)
    sample_list <- PrepSCTIntegration(object.list = so_list, anchor.features = sample_feat)
    sample_anchors <- FindIntegrationAnchors(object.list = sample_list, normalization.method = "SCT",
                                             anchor.features = sample_feat)
    so_int <- IntegrateData(anchorset = sample_anchors, normalization.method = "SCT")
    saveRDS(so_int, file = paste0(save.as, "_int_SCT_seurat_object.rds"))
  }
  else{
    sample_anchors <- FindIntegrationAnchors(object.list = so_list, dims = 1:30)
    so_int <- IntegrateData(anchorset = sample_anchors, dims = 1:30)
    so_int <- ScaleData(so_int, vars.to.regress = c("percent.mt", "percent.rp"))
    saveRDS(so_int, file = paste0(save.as, "_int_seurat_object.rds"))
  }
  return(so_int)
}  
  
#############################################################################
########################## SEURAT ANALYSIS =====########################## 

downstream_analysis <- function(object, marker.genes = NULL, save.rds = getwd(), save.as = getwd(),
                                assay='integrated', cls.res=0.6, SCT=F, reduc = 'umap',
                                plt.metadata = c("seurat_clusters", "sample_id", "condition_id")){
  # Runs dimensional reduction (PCA, UMAP and tSNE), conduct unsupervised clustering. 
  # Plot dimensional reductions and feature plots
  DefaultAssay(object) <- assay
  object <- RunPCA(object, npcs = 30)
  object <- RunUMAP(object, reduction = "pca", dims = 1:30)
  object <- RunTSNE(object, dims = 1:30)
  object <- FindNeighbors(object, dims = 1:30)
  object <- FindClusters(object, reduction.type = "pca", dims.use = 1:30,
                         resolution = cls.res, print.output = 0, save.SNN = T)
  object@meta.data$cluster_num<- object@meta.data$seurat_clusters

  saveRDS(object, file = paste0(save.rds, "_int_processed_seurat_object.rds"))

  # plotting
  pca <- PCAPlot(object, dims = c(1, 2))
  # ggsave(pca, file=paste(save.as,"_PCAplot.eps",sep=""), width = 8.5, height = 11)
  ggsave(pca, file=paste(save.as,"_PCAplot.pdf",sep=""), width = 8.5, height = 11)

  for(i in 1:length(plt.metadata)){
    Idents(object) <- plt.metadata[i]
    if (reduc == 'umap'){
      UMAPPlot(object, label = T, repel = T)
      # ggsave(file=paste0(save.as, "_UMAPplot_",plt.metadata[i],".eps"), width = 8.5, height = 11)
      ggsave(file=paste0(save.as, "_UMAPplot_",plt.metadata[i],".pdf"), width = 8.5, height = 11)
    }
    if (reduc == 'tsne'){
      TSNEPlot(object, label = T, repel = T)
      # ggsave(file=paste0(save.as, "_TSNEplot_", plt.metadata[i], ".eps"), width = 8.5, height = 11)
      ggsave(file=paste0(save.as, "_TSNEplot_", plt.metadata[i], ".pdf"), width = 8.5, height = 11)
    }
    table_sum <- table(object$seurat_clusters, Idents(object))
    write.csv(table_sum, file = paste0(save.as, "_table_summary_",plt.metadata[i],".csv"))
  }
  Idents(object) <- "seurat_clusters"
  
  if(length(marker.genes)>0){
    if (SCT==T){DefaultAssay(object) <- 'SCT'} else{DefaultAssay(object) <- 'RNA'}

    FeaturePlot(object, features = marker.genes, reduction = reduc)
    # ggsave(file=paste0(save.as,"_feature_plots_marker.eps"), width = 8.5, height = 11)
    ggsave(file=paste0(save.as,"_feature_plots_marker.pdf"), width = 8.5, height = 11)

    VlnPlot(object, features = marker.genes, pt.size = 0, ncol = 1)
    # ggsave(file=paste0(save.as,"_vlnplot_marker.eps"), width = 11, height = 25)
    ggsave(file=paste0(save.as,"_vlnplot_marker.pdf"), width = 11, height = 25)
  }

  VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, ncol = 1)
  # ggsave(file=paste0(save.as,"_vlnplot_feat_count.eps"), width = 8.5, height = 15)
  ggsave(file=paste0(save.as,"_vlnplot_feat_count.pdf"), width = 8.5, height = 15)
  return(object)

}

Find_DEG <- function(object, by_ident = NULL, assay = "SCT",
                     save.as = getwd(), hm = F, top20 = F){
  # Find differentially expressed genes between identities specified by by.ident using the SCT or RNA assay
  DefaultAssay(object) <- assay
  if(is.null(by_ident) == F){Idents(object) <- by_ident}
  DEG <- FindAllMarkers(object)
  DEG <- DEG[(order(DEG[,"cluster"],-DEG[,"avg_logFC"])),]
  write.csv(DEG, file = paste0(save.as,"_highly_expressed_gene_cluster.csv"))
  print("writing csv")
  DEG <- DEG[which(DEG$p_val_adj<=0.05 & DEG$avg_logFC>=0),]
  #DEG <- DEG[which(DEG$avg_logFC>=0),]
  write.csv(DEG, file = paste0(save.as,"_highly_expressed_gene_cluster_sig_pos.csv"))
  
  if (top20 ==T){
    top20 <- DEG %>% group_by(cluster) %>% top_n(20, avg_logFC)
    top20 <- top20[order(top20$cluster,-top20$avg_logFC),]
    write.csv(top20, file = paste0(save.as,"_highly_expressed_gene_cluster_top20.csv"))
  }
  
  if (hm == T){
    print("making heatmap")
    top10 <- DEG %>% group_by(cluster) %>% top_n(10, avg_logFC)
    DoHeatmap(object, features = top10$gene, label = TRUE, assay = assay )
    ggsave(file=paste0(save.as,"_heatmap_top10_high_genes_cluster.eps"), width = 11, height = 11)
    ggsave(file=paste0(save.as,"_heatmap_top10_high_genes_cluster.pdf"), width = 11, height = 11)
  }
  
  # if (vlnplt == T){
  #   print("making violin plot")
    # top_deg <- DEG %>% group_by(cluster) %>% top_n(1, avg_logFC)
    # for (i in 1:length(top_deg)){
    #   DoHeatmap(object, features = top_deg$gene[i:i+6], label = TRUE, assay = assay )
    #   ggsave(file=paste0(save.as,"_heatmap_top10_high_genes_cluster.eps"), width = 11, height = 18)
    #   ggsave(file=paste0(save.as,"_heatmap_top10_high_genes_cluster.pdf"), width = 11, height = 18)
    #   i = i+7
    # }
  # }
  
  return(DEG)
}



#############################################################################
########################## CONDITION ANALYSIS =====########################## 

recluster <- function(object, cluster=NULL, save.as = getwd(), save.rds = getwd(),
                      assay="integrated", res=0.6, inv=F, SCT=F, reduc = 'umap',
                      plt.metadata = c("seurat_clusters", "sample_id", "condition"),
                      marker.genes = NULL){
  # subset and reclusters object; saves reclustered object; basically the same as downstream_analysis function
  # potentially redundant
  DefaultAssay(object) <- assay
  if(length(cluster)==0){obj_subset <- object}
  else {obj_subset <- subset(object, idents = cluster, invert = inv)}
  obj_subset <- RunPCA(obj_subset, npcs = 30)
  obj_subset <- RunTSNE(obj_subset, dims = 1:30)
  obj_subset <- RunUMAP(obj_subset, dims = 1:30)
  obj_subset <- FindNeighbors(obj_subset, dims = 1:30)
  obj_subset <- FindClusters(obj_subset, resolution = res)
  # obj_subset <- FindClusters(obj_subset, reduction.type = "pca", dims.use = 1:30, 
  #                           resolution = res, print.output = 0, save.SNN = T)
  obj_subset@meta.data$new_cluster_num <- obj_subset@active.ident
  saveRDS(obj_subset, file = paste0(save.rds, "_seurat_object.rds"))
  
  print("plotting")
  for(i in 1:length(plt.metadata)){
    Idents(object) <- plt.metadata[i]
    if (reduc == 'umap'){
      UMAPPlot(object, label = T, repel = T)
      # ggsave(file=paste0(save.as, "_UMAPplot_",plt.metadata[i],".eps"), width = 8.5, height = 11)
      ggsave(file=paste0(save.as, "_UMAPplot_",plt.metadata[i],".pdf"), width = 8.5, height = 11)
    }
    if (reduc == 'tsne'){
      TSNEPlot(object, label = T, repel = T)
      # ggsave(file=paste0(save.as, "_TSNEplot_", plt.metadata[i], ".eps"), width = 8.5, height = 11)
      ggsave(file=paste0(save.as, "_TSNEplot_", plt.metadata[i], ".pdf"), width = 8.5, height = 11)
    }
    table_sum <- table(object$seurat_clusters, Idents(object))
    write.csv(table_sum, file = paste0(save.as, "_table_summary_",plt.metadata[i],".csv"))
  }
  Idents(object) <- "seurat_clusters"
  
  if(length(marker.genes)>0){
    if (SCT==T){DefaultAssay(object) <- 'SCT'}
    else{DefaultAssay(object) <- 'RNA'}
    FeaturePlot(obj_subset, features = marker.genes, reduction = reduc)
    # ggsave(file=paste0(save.as,"_feature_plots_marker.eps"), width = 8.5, height = 11)
    ggsave(file=paste0(save.as,"_feature_plots_marker.pdf"), width = 8.5, height = 11)
    
    VlnPlot(object, features = marker.genes, pt.size = 0, ncol = 1)
    # ggsave(file=paste0(save.as,"_vlnplot_marker.eps"), width = 11, height = 25)
    ggsave(file=paste0(save.as,"_vlnplot_marker.pdf"), width = 11, height = 25)
  }
  VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, ncol = 1)
  # ggsave(file=paste0(save.as,"_vlnplot_feat_count.eps"), width = 11, height = 25)
  ggsave(file=paste0(save.as,"_vlnplot_feat_count.pdf"), width = 11, height = 25)
  return(obj_subset)
}

conditionDEG <- function(obj, cluster, ident, DEG.ident = "condition_id", save.as = getwd(), assay="RNA"){
  # Find DEG between each (experimental) condition for every cluster specified from the object
  DefaultAssay(obj) <- assay
  listDF <- vector(mode = "list", length = length(cluster))
  for (i in 1:length(cluster)) {
    Idents(obj) <- ident
    cls <- subset(obj, idents = cluster[i])
    Idents(cls) <- DEG.ident
    listDF[[i]] <- FindAllMarkers(cls)
    listDF[[i]]['cluster_num'] = rep(cluster[i], dim(listDF[[i]])[1])
    listDF[[i]] <-listDF[[i]][(order(listDF[[i]][,"cluster"],-listDF[[i]][,"avg_logFC"])),]
  }
  df <- do.call("rbind", listDF)
  write.csv(df, file = paste0(save.as, "_conditions_DEG.csv"))
  
  df_pos <- df[which(df$p_val_adj<=0.05 & df$avg_logFC>0),]
  write.csv(df_pos, file = paste0(save.as, "_conditions_DEG_sig_pos.csv"))
  
  df_neg <- df[which(df$p_val_adj<=0.05 & df$avg_logFC<=0),]
  write.csv(df_neg, file = paste0(save.as, "_conditions_DEG_sig_neg.csv"))
  return(df)
}
  
# change_levels <- function(object, sort.by){
#   # format plot ordering and stuff, plot_id is not saved in the
#   # seurat object
#   Idents(object) <- sort.by
#   levels(object) <- sort(levels(object))
#   return(object)
# }

plot_cond_DEG <- function(obj, DEG, cond, save.as = getwd(), 
                          col_scheme = c("cyan", "purple4"), no.rp = F, no.mt = F, add_genes = NULL,
                          assay="SCT", by.ident = "plot_id"){
  # Finds unique DEG from condition_DEG function and plot each condition for a cell type 
  DefaultAssay(obj) <- assay
  ori_ident <- Idents(obj)
  Idents(obj) <- by.ident
  levels(obj) <- sort(levels(obj))
  #list_deg <- list()
  for (i in 1:length(cond)){
    if (no.rp == T){DEG <- DEG[!grepl("^RP(S|L)", DEG$gene),]}
    if (no.mt == T){DEG <- DEG[!grepl("^MT)", DEG$gene),]}
    cond_DEG <- unique(DEG[which(DEG$cluster == cond[i]),]$gene)
    #list_deg[[i]] <- cond_DEG
    
    DotPlot(obj, features = cond_DEG, cols = col_scheme) + 
      FontSize(x.text = 10, y.text = 10) + RotatedAxis()
    if (length(cond_DEG) >= 15) {
      ggsave(file=paste0(save.as , "_", cond[i], "_dot_plot.pdf"), width = length(cond_DEG), height = 8.5)
      ggsave(file=paste0(save.as , "_", cond[i], "_dot_plot.eps"), width = length(cond_DEG), height = 8.5)
    } else {
      ggsave(file=paste0(save.as , "_", cond[i], "_dot_plot.pdf"), width = 11, height = 8.5)
      ggsave(file=paste0(save.as , "_", cond[i], "_dot_plot.eps"), width = 11, height = 8.5)
    }
  }
  Idents(obj) <- ori_ident
  write.table(DEG, file=paste0(save.as, "unique_cond_DEG.csv"))
  return(DEG)
}


# Scatter plot with highlighted genes
highlight_genes <- function(object, name=NULL, ident = NULL, save.as=getwd(), assay="SCT", col=c("red"), top.gene=10){
  plot_list <- list()
  for(i in 1:length(object)){
    obj <- object[[i]]
    if(is.null(ident) == F) {Idents(obj) <- ident}
    avg <- AverageExpression(obj, verbose = F, assay='SCT',return.seurat = T)$SCT
    pts <- as.data.frame(test@data)
    pts$gene <- row.names(pts)
    pts$diff <- pts[,2] - pts[,1]
    pts <- pts[order(-pts$diff),]
    
    plt <- CellScatter(avg, cell1 = colnames(avg)[1], cell2 = colnames(avg)[2]) + 
      geom_point(colour='grey') + ggtitle(name[i])

    if (length(col)== length(object)){
      plt <- plt <-LabelPoints(plot = plt, points=pts$gene[1:top.gene], repel=T) +
        geom_point(data=pts[1:top.gene,], colour=col[i])
    }
    else{plt <-LabelPoints(plot = plt, points=pts$gene[1:top.gene], repel=T) +
      geom_point(data=pts[1:top.gene,], colour=col[1])
    print("number of colour is not equal to number of objects")
    }
    plot_list[[i]] <- plt
  }
  pdf(file = paste0(save.as, "_gene_scatter_plt.pdf"), width = 8.5, height = 11)
  for (j in 1:length(plot_list)){
    print(plot_list[[j]])
  }
  dev.off()
}

cluster_DEG <- function(obj, ident = "cell_name", cls.list = NULL,
                        cond = NULL, save.as = "cluster_analysis/"){
  #Find DEG between each cluster, separated by condition
  Idents(obj) <- ident
  for (i in 1:length(cls.list)){
    #df_DEG <- vector(mode = "list", length = length(cls.list))
    for (j in 1:length(cond)){
      so_sub <- subset(obj, idents = as.vector(unlist(cls.list[i])), subset = condition_id == cond[j])
      DEG <- FindAllMarkers(so_sub, assay = "SCT")
      DEG["condition"] <- rep(cond[j], dim(DEG)[1])
      DEG <- DEG[(order(DEG[,"cluster"], -DEG[,"avg_logFC"])), ]
      #assign(paste0(names(cls.list)[i], cond[j], "_obj"), so_sub)
      write.csv(DEG, file = paste0(save.as, names(cls.list)[i],"_", cond[j], "_per_cond_DEG.csv"))
      
      DEG <- DEG[which(DEG$p_val_adj<=0.05 & DEG$p_val_adj>0 & DEG$avg_logFC>=0),]
      write.csv(DEG, file = paste0(save.as, names(cls.list)[i],"_", cond[j], "_per_cond_sig_pos.csv"))
    }}}

# cluster_DEG <- function(obj, ident = "cell_name", cls.list = list(delta = d_cluster),
#                         cond = c("Basal"), save.as = "cluster_analysis/"){
#   #Find DEG between each cluster, separated by condition
#   Idents(obj) <- ident
#   for (i in 1:length(cls.list)){
#     df_DEG <- vector(mode = "list", length = length(cls.list))
#     for (j in 1:length(cond)){
#       so_sub <- subset(obj, idents = as.vector(unlist(cls.list[i])), subset = condition_id == cond[j])
#       df_DEG[[j]] <- FindAllMarkers(so_sub, assay = "SCT")
#       df_DEG[[j]]["condition"] <- rep(cond[j], dim(df_DEG[[j]])[1])
#       df_DEG[[j]] <- df_DEG[[j]][(order(df_DEG[[j]][,"cluster"], -df_DEG[[j]][,"avg_logFC"])), ]
#       assign(paste0(names(cls.list)[i], cond[j], "_obj"), so_sub)
#     }
#     DEG <- do.call("rbind", df_DEG)
#     write.csv(DEG, file = paste0(save.as, names(cls.list)[i], "_per_cond_DEG.csv"))
#     
#     DEG <- DEG[which(DEG$p_val_adj<=0.05 & DEG$p_val_adj>0 & DEG$avg_logFC>=0),]
#     write.csv(DEG, file = paste0(save.as, names(cls.list)[i], "_per_cond_sig_pos.csv"))
#   }
#   return(df_DEG)
# }

