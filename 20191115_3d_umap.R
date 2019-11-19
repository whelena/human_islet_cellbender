library(Seurat)
library(plotly)

#server_address <- "~/cfrifs02/LynnLabNew/"
server_address <- "W:/LynnLabNew/"

folder_address <- paste0(server_address, "Helena/20191015_human_islet_cellbender/")
setwd(paste0(folder_address, "analysis/"))

endo <- readRDS(file="../seurat_objects/endo_8_filtered_named_seurat_object.rds")
#hislet <- readRDS(file="../seurat_objects/hislet_cb_filtered_named_seurat_object.rds")

endo_3d <- RunUMAP(endo, reduction = "pca", dims = 1:30, n.components = 3L)

# exctract umap info into variables
umap_1 <- endo_3d[["umap"]]@cell.embeddings[,1]
umap_2 <- endo_3d[["umap"]]@cell.embeddings[,2]
umap_3 <- endo_3d[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = endo_3d, reduction = "umap")

# Prepare a dataframe for cell plotting
plotting.data <- FetchData(object = endo_3d, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters", "cell_name"))

# Make a column of row name identities (these will be your cell/barcode names)
plotting.data$label <- plotting.data$cell_name

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plotting.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~seurat_clusters, 
        # colors = c("lightseagreen",
        #            "gray50",
        #            "darkgreen",
        #            "red4",
        #            "red",
        #            "turquoise4",
        #            "black",
        #            "yellow4",
        #            "royalblue1",
        #            "lightcyan3",
        #            "peachpuff3",
        #            "khaki3",
        #            "gray20",
        #            "orange2",
        #            "royalblue4",
        #            "yellow3",
        #            "gray80",
        #            "darkorchid1",
        #            "lawngreen",
        #            "plum2",
        #            "darkmagenta"),
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 2, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
