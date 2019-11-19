# master script for human pancreatic islet scRNA seq analysis
# Helena Winata 10/15/2019

#server_address <- "W:/LynnLabNew/"
server_address <- "~/cfrifs02/LynnLabNew/"

folder_address <- paste0(server_address, "Helena/20191015_human_islet_cellbender/")
setwd(paste0(folder_address, "analysis/"))

# specify h5 path (CellBender output) and set up file paths
h5_path <- paste0(folder_address, "10x_raw_aggr_files/")
scater_obj_path <- paste0(folder_address,"scater_obj_aggr/")
seurat_obj_path <- paste0(folder_address,"seurat_objects/")

# input stuff
name <- "hislet_aggr_cb" # name of dataset
marker_genes <- c("INS", "GCG", "SST", "PPY", "CPA1", "KRT19", "PECAM1", "PDGFRB", "LYZ")

#sample_id <- c("R253", "R282", "R317")  # length must be the same as number of files to be integrated
add_meta <- c("library_id","total_counts","pct_counts_MT","condition") # additional metadata from scater object

# Run Scater for each h5 file (see scater code for annotation)
raw_folder <- list.files(path = h5_path)
keyword_list <- gsub("_cb_aggr", "", raw_folder)
for(i in length(raw_h5)){
  keyword <- keyword_list[i]
  raw_folder_address <- paste0(h5_path, raw_folder[i])
  saveRDS(raw_folder_address,file="10x_raw_folder_address.rds")
  saveRDS(keyword,file = "keyword.rds")
  
  scater_file_address <- paste0(scater_obj_path, keyword,"_")
  saveRDS(scater_file_address, file = "scater_file_address.rds")

  #run Scater script=====
  source(paste0(server_address, 
                "Helena/10x_automated_pipeline/hislet_analysis/20191510_scater_hislet_h5.R"))
  print("finished scater")
}


# load pipeline functions
source(paste0(server_address, "Helena/10x_automated_pipeline/hislet_analysis/20190920_seurat_fxns.R"))
donor.id <- sub("_r.*", "", keyword_list)
sample_list <- create_so_list(scater_obj_path, add.meta = add_meta, sample.id = keyword_list, SCT=T,
                              save.rds = paste0(seurat_obj_path, name),
                              save.qc = seurat_obj_path)

so_int <- create_so_int(sample_list, save.as = paste0(seurat_obj_path, "hislet_aggr_cb"), SCT=T)

hislet <- downstream_analysis(so_int, marker.genes = marker_genes, SCT=T,
                              save.rds = paste0(seurat_obj_path, "hislet_aggr_cb"),
                              save.as = "seurat/hislet_aggr_cb")

# remove unrepresentative clusters (most cells from same donor/condition)
hislet <- subset(hislet, idents = c("14"), invert=T)
cell_names <- list("B1", "a1", "a2","B2", "B3", "a3", "a3", "a4", "INS/GCG1",
                   "Pericytes1", "d", "B4", "Ductal", "INS/SST", "Endothetial", 
                   "Pericytes2", "INS/GCG2","Acinar" , "PP", "INS/GCG3", "a5", "Immune")
names(cell_names) <- levels(hislet)
hislet <- RenameIdents(hislet, cell_names)
hislet[["cell_name"]] <- Idents(hislet)
saveRDS(hislet, file = paste0(seurat_obj_path, "hislet_cb_filtered_named_seurat_object.rds"))
UMAPPlot(hislet, label = T, repel = T)
ggsave(file="seurat/hislet_cb_UMAPplot_cellname.pdf", width = 11, height = 11)

# subset INS, GCG, SST cells
endo <- subset(hislet, idents = c(paste0("a", 1:5), paste0("B", 1:4), "d", "INS/SST",
                                  paste0("INS/GCG", 1:3))) 
endo <- downstream_analysis(endo, marker.genes = c("INS", "GCG", "SST"), SCT=T,
                            save.rds = paste0(seurat_obj_path, "endo"),
                            save.as = "seurat/endo/endo")
endo$endo_cls <- endo$cluster_num

# check out polyhormonal cells
# multi_endo <- subset(endo, idents = c("11", "17", "21"))
# multi_endo <- downstream_analysis(multi_endo, marker.genes = c("INS","GCG", "SST"), SCT=T,
#                             save.rds = paste0(seurat_obj_path, "multi_endo"),
#                             save.as = "seurat/endo/multi_endo")

# pick clusters to remove from endo, invert=T means remove clusters specified in idents
# invert=F means remove clusters NOT specified in idents
endo <- subset(endo, idents = c(13:17), invert = T)
endo_names <- c("B1", "B2", "B3", "a1", "a2", "B4", "a3", "a4", "a5", "INS/GCG",
                "d", "INS/SST", "B5")
names(endo_names) <- levels(endo)
endo <- RenameIdents(endo, endo_names)
endo[["cell_name"]] <- Idents(endo)
endo@meta.data$plot_id <- paste0(endo@meta.data$cell_name, "_", endo@meta.data$condition_id)
saveRDS(endo, file = paste0(seurat_obj_path, "endo_filtered_named_seurat_object.rds"))

endo_DEG <- Find_DEG(endo, by_ident = "cell_name", save.as = "seurat/endo", hm = T, 
                     top20 = T, assay = "SCT")

###########################################
loc <- "condition_analysis/"
# ploting colour_scheme
col_scheme <- c("cyan", "purple4")
heatmap_col <- c("darkturquoise", "white", "magenta")
cond <- c("Basal", "Negative", "Positive")

############# Alpha cell clusters =====
a_cluster <- c(paste0("a", 1:5), "INS/GCG")
saveRDS(a_cluster, file = "a_cls.rds")
# Find DEG between conditon in each cell type
a.plot <- subset(endo, idents = a_cluster)

# Find DEG between clusters and between conditions
a_DEG <- Find_DEG(a.plot, by_ident = "cell_name", assay="SCT", save.as="seurat/endo/alpha")
a_DEG_cond <- Find_DEG(a.plot, by_ident = "condition_id", assay = "SCT", save.as = paste0(loc,"alpha"))
DotPlot(a.plot, features = a_DEG_cond$gene, cols = col_scheme)
ggsave(file=paste0(loc, "alpha_cond_dot_plot.pdf"), width = 8.5, height = 11)
# ggsave(file=paste0(loc, "alpha_cond_dot_plot.eps"), width = 11, height = 8.5)

# Find DEG between each condition for each cluster
Idents(a.plot) <- "plot_id"
levels(a.plot) <- sort(levels(a.plot))
a_cond_DEG <- conditionDEG(endo, a_cluster, ident = "cell_name", assay = "SCT",
                           save.as = paste0(loc, "alpha"))
a_plt_DEG <- plot_cond_DEG(a.plot, a_cond_DEG, cond, save.as = "condition_analysis/alpha", no.rp = T)

############# Beta cell clusters =====
b_cluster <- c(paste0("B", 1:5), "INS/GCG", "INS/SST")
saveRDS(b_cluster, file = "b_cls.rds")
# Find DEG between conditon in each cell type
b.plot <- subset(endo, idents = b_cluster)

# Find DEG between clusters and between conditions
b_DEG <- Find_DEG(b.plot, by_ident = "cell_name", assay="SCT", save.as="seurat/endo/beta")
b_DEG_cond <- Find_DEG(b.plot, by_ident = "condition_id", assay = "SCT", save.as = paste0(loc, "beta"))
DotPlot(b.plot, features = b_DEG$gene, cols = col_scheme) + FontSize(x.text = 8)
ggsave(file=paste0(loc, "beta_cond_dot_plot.pdf"), width = 11, height = 8.5)
# ggsave(file=paste0(loc, "beta_cond_dot_plot.eps"), width = 11, height = 8.5)

# Find DEG between each condition for each cluster
Idents(b.plot) <- "plot_id"
levels(b.plot) <- sort(levels(b.plot))
b_cond_DEG <- conditionDEG(endo, b_cluster, ident = "cell_name", assay = "SCT",
                           save.as = paste0(loc, "beta"))
b_plt_DEG <- plot_cond_DEG(b.plot, b_cond_DEG, cond, save.as = "condition_analysis/beta", no.rp = T)

############# delta cell clusters =====
d_cluster <- c("d", "INS/SST")
saveRDS(d_cluster, file = "d_cls.rds")
# Find DEG between conditon in each cell type
d.plot <- subset(endo, idents = d_cluster)

# Find DEG between clusters and between conditions
d_DEG <- Find_DEG(d.plot, by_ident = "cell_name", assay="SCT", save.as="seurat/endo_cells/delta")
d_DEG_cond <- Find_DEG(d.plot, by_ident = "condition_id", save.as = paste0(loc, "delta"), assay = "SCT")

DotPlot(d.plot, features = d_DEG_cond$gene, cols = col_scheme) + FontSize(x.text = 8)
ggsave(file=paste0(loc, "delta_cond_dot_plot.pdf"), width = 11, height = 8.5)
# ggsave(file=paste0(loc, "delta_cond_dot_plot.eps"), width = 11, height = 8.5)

# Find DEG between each condition for each cluster
Idents(d.plot) <- "plot_id"
levels(d.plot) <- sort(levels(d.plot))
d_cond_DEG <- conditionDEG(endo, d_cluster, ident = "cell_name", assay = "SCT",
                           save.as = paste0(loc, "delta"))
d_plt_DEG <- plot_cond_DEG(d.plot, d_cond_DEG, cond, save.as = "condition_analysis/delta", no.rp = T)


#############################################################################
########################## CLUSTER ANALYSIS =====############################
cell_list <- list( alpha = a_cluster, beta = b_cluster, delta = d_cluster)
loc2 <- "cluster_analysis/"

clusterDEG <- cluster_DEG(endo, cls.list = cell_list, cond = cond, save.as = loc2)


