#server_address <- "W:/LynnLabNew/"
server_address <- "~/cfrifs02/LynnLabNew/"

folder_address <- paste0(server_address, "Helena/20190916_human_islet_refined/")
setwd(paste0(folder_address, "analysis/"))
source(paste0(server_address, "Helena/10x_automated_pipeline/hislet_analysis/20190920_seurat_fxns.R"))

endo <- readRDS(file = "../seurat_objects/endo_named_seurat_object.rds")
col_scheme <- c("cyan", "purple4")
heatmap_col <- c("darkturquoise", "white", "magenta") 
a_cluster <- readRDS(file = "a_cls.rds")
b_cluster <- readRDS(file = "b_cls.rds")
d_cluster <- readRDS(file = "d_cls.rds")

a.plot <- subset(endo, idents = a_cluster)
b.plot <- subset(endo, idents = b_cluster)
d.plot <- subset(endo, idents = d_cluster)

a_only <- subset(a.plot, idents = "INS/GCG", invert = T)
b_only <- subset(b.plot, idents = c("INS/GCG", "INS/SST"), invert = T)

############## calculate average expressions =====
b.only.avg <- AverageExpression(b_only, return.seurat = T, verbose = F)
b.no.sst.avg <- AverageExpression(b.no.sst, return.seurat = T, verbose = F)

#w separated by conditions
b.plot <- change_levels(b.plot, sort.by = "plot_id")
b.averages <- AverageExpression(b.plot, return.seurat = T)


############## Alpha cell maturity markers =====
mature_a_marker <- c("GCG", "TTR", "LOXL4", "MAFB", "ARX", "IRX1", "IRX2")
DotPlot(a_only, features = mature_a_marker[1:2], cols = col_scheme, assay = "SCT") + 
  FontSize(x.text = 10, y.text = 10) + RotatedAxis()
ggsave(file=paste0(loc, "alpha_mature_dot_plot_high.pdf"), width = 8.5, height = 11)
ggsave(file=paste0(loc, "alpha_mature_dot_plot_high.eps"), width = 8.5, height = 11)

DotPlot(a_only, features = mature_a_marker[3:7], cols = col_scheme, assay = "SCT", dot.scale = 10) + 
  FontSize(x.text = 10, y.text = 10) + RotatedAxis()
ggsave(file=paste0(loc, "alpha_mature_dot_plot_low.pdf"), width = 8.5, height = 11)
ggsave(file=paste0(loc, "alpha_mature_dot_plot_low.eps"), width = 8.5, height = 11)

############## Beta cell maturity markers =====
mature_b_marker <- c("INS", "IAPP", "CHGA", "PCSK1", "PAX6", "NEUROD1",
                     "MAFA", "FOXO1", "PDX1", "NKX2-2", "NKX6-1", "UCN3")
# "G6PC2", "ERO1LB", "GLUT2"
immature_b_marker <- c("MAFB", "PRDM16", "SOX4", "GATA4", "FEV" )
DotPlot(b.plot, features = mature_b_marker, cols = col_scheme, assay = "SCT") + 
  FontSize(x.text = 10, y.text = 10) + RotatedAxis()
ggsave(file=paste0(loc, "beta_mature_dot_plot.pdf"), width = 15, height = 8.5)
ggsave(file=paste0(loc, "beta_mature_dot_plot_1.eps"), width = 15, height = 8.5)

DotPlot(b.plot, features = immature_b_marker, cols = col_scheme, assay = "SCT") + 
  FontSize(x.text = 10, y.text = 10) + RotatedAxis()
ggsave(file=paste0(loc, "beta_immature_dot_plot.pdf"), width = 8.5, height = 11)
ggsave(file=paste0(loc, "beta_immature_dot_plot.eps"), width = 8.5, height = 11)

# revisions
DotPlot(b.plot, features = c("SOX4", "MAFB"), cols = col_scheme, assay = "SCT") + 
  FontSize(x.text = 10, y.text = 10) + RotatedAxis()
ggsave(file=paste0(loc, "beta_sox4_mafb_dot_plot.pdf"), width = 8.5, height = 11)
DoHeatmap(b.averages, features = c("SOX4", "MAFB"), draw.lines = F) +
  scale_fill_gradientn(colors = heatmap_col) + theme(legend.position="bottom")
ggsave(file=paste0(loc, "beta_sox4_mafb_heatmap.pdf"), width = 11, height = 8.5)

DotPlot(b.plot, features = c("CPE", "PCSK1", "PCSK2", "CHGA", "IAPP", "INS", "NEUROD1", "PAX6"), 
        cols = col_scheme, assay = "SCT") + FontSize(x.text = 10, y.text = 10) + RotatedAxis()
ggsave(file=paste0(loc, "beta_high_mature_dot_plot.pdf"), width = 15, height = 8.5)
DoHeatmap(b.averages, features = c("CPE", "PCSK1", "PCSK2", "CHGA", "IAPP", "INS", "NEUROD1", "PAX6"),
          draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom")
ggsave(file=paste0(loc, "beta_high_mature_heatmap.pdf"), width = 11, height = 8.5)

DotPlot(b.plot, features = c("NKX6-1", "NKX2-2", "PDX1"), cols = col_scheme, 
        assay = "SCT") + FontSize(x.text = 10, y.text = 10) + RotatedAxis()
ggsave(file=paste0(loc, "beta_low_TF_dot_plot.pdf"), width = 11, height = 8.5)

DoHeatmap(b.averages, features = c("NKX6-1", "NKX2-2", "PDX1"), draw.lines = F) +
  scale_fill_gradientn(colors = heatmap_col) + theme(legend.position="bottom")
ggsave(file=paste0(loc, "beta_low_TF_heatmap.pdf"), width = 11, height = 8.5)


############## Delta cell maturity markers =====
mature_d_marker <- c("SST", "RBP4", "HHEX")
DotPlot(d.plot, features = mature_d_marker, cols = col_scheme, assay = "SCT") + 
  FontSize(x.text = 10, y.text = 10) + RotatedAxis()
ggsave(file=paste0(loc, "delta_mature_dot_plot.pdf"), width = 8.5, height = 11)
ggsave(file=paste0(loc, "delta_mature_dot_plot.eps"), width = 8.5, height = 11)



################ Scatter Plots =====
# Ins vs GCG
FeatureScatter(test, cells = WhichCells(endo, idents = c("a1", "b1", "INS/GCG")),
               feature1 = "INS", feature2 = "GCG") + coord_cartesian(xlim = c(0, 8), ylim = c(0, 8))
ggsave(file=paste0(loc, "ins_vs_gcg_scatter_plt.pdf"), width = 11, height = 8.5)
ggsave(file=paste0(loc, "ins_vs_gcg_scatter_plt.pdf"), width = 11, height = 8.5)

# Ins vs SST
FeatureScatter(test, cells = WhichCells(endo, idents = c("d", "b1", "INS/SST")),
               feature1 = "INS", feature2 = "SST") + coord_cartesian(xlim = c(0, 8), ylim = c(0, 8))
ggsave(file=paste0(loc, "ins_vs_sst_scatter_plt.pdf"), width = 11, height = 8.5)
ggsave(file=paste0(loc, "ins_vs_sst_scatter_plt.pdf"), width = 11, height = 8.5)


# Scatter plot
highlight_genes(object=c(b.cells, a.cells, b2.cells, d.cells), ident = "cell_name",
                name=c("B vs. INS/GCG", "a vs. INS/GCG", "B vs. INS/SST", "d vs. INS/SST"), 
                col=c("red", "blue", "red", "darkgreen"), save.as=loc)

################# VLN PLOTS ====
for_beta <- c("UCN3", 'MAFA', "MAFB", "ERO1B", "PDX1", "NKX6-1", "NKX2-2", "INS", "IAPP", "CPE",
              "PCSK1", "PCSK2")
b.no.sst <- subset(b.plot, idents = "INS/SST", invert = T)
VlnPlot(b.no.sst, features = for_beta, assay='SCT', ncol = 1, pt.size = 0) + 
  FontSize(x.text=6, y.text=6, x.title=8, y.title=8, main = 8)
ggsave(file=paste0("seurat/endo_cells/beta_vlnplot.pdf"), width = 8.5, height = 35)
ggsave(file=paste0("seurat/endo_cells/beta_vlnplot.eps"), width = 8.5, height = 35)


DotPlot(b.no.sst, features = for_beta, assay = 'SCT', cols=col_scheme, dot.scale=10) + RotatedAxis()
ggsave(file=paste0("seurat/endo_cells/beta_dotplot.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/beta_dotplot.eps"), width = 11, height = 8.5)

DoHeatmap(b.only.avg, features = for_beta, assay='SCT', draw.lines = F) +
  scale_fill_gradientn(colors = heatmap_col) + theme(legend.position="bottom")
ggsave(file=paste0("seurat/endo_cells/beta_only_heatmap2.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/beta_only_heatmap2.eps"), width = 11, height = 8.5)

DoHeatmap(b.only.avg, features = c("INS", "NKX6-1", "PDX1", "GCK", "SERCA2", "UCN3", "G6PC2", "ERO1LB",
                                   "ATP2A1", "ATP2A2"),
          assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom")
ggsave(file=paste0("seurat/endo_cells/beta_only_heatmap.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/beta_only_heatmap.eps"), width = 11, height = 8.5)


tca.genes <- c("CS", "ACO2", "IDH3G", "IDH3B", "IDH3A", "OGDH", "DLST", "SUCLG1",
               "SUCLG2", "SUCLA2", "SDHB", "SDHA", "SDHD", "SDHC", "FH", "MDH1",
               "MDH2", "ME1", "ME2", "ME3", "PDHA1")
ox.phos.genes <- c("PPA1", "MT-ATP6", "MT-ATP8", "ATP4A", "ATP5A1","ATP5B", "ATP5C1", "ATP5D", "ATP5E",
                   "ATP5F1", "ATP5F1C", "ATP5G1", "ATP5G2", "ATP5G3", "ATP5H", "ATP5I", "ATP5J", "ATP5J2",
                   "ATP5L", "ATP5MC1", "ATP5MF", "ATP5O", "ATP5PB", "ATP6V1A", "COX6A1", "COX6A2", "COX6B",
                   "COX7A2", "COX8A", "NDUFA13", "NDUFB2", "NDUFB3", "NDUFC2", "NDUFS2")
fig.genes <- c("HK1", "HK2", "LDHA", "SLC16A1", "ACOT7", "DNMT3A", "RFX6", "ADCYAP1", "CCKAR",
               paste0("ADCY", 1:9), "PLCG1", "PLCG2", "PLCL2", "PLCB1", "PLCB3", "PLCB4", "PLCD1",
               "PLCD3", "CREB1", "ITPR3", "CAMK1", "CAMK4", "CAMK1G", "CAMK2G", "CAMK1D", "CAMK2D",
               "CAMK2B", "SYT1")
ca.dep <- c("IAPP", "IER3", "RGS16", "FOS", "NR4A2", "NR4A1", "BTG2", "GADD45B", "TXNIP", "ZNF331")
saveRDS(tca.genes, file = "TCA_cycle_genes.rds")
saveRDS(ox.phos.genes, file = "oxphos_genes.rds")
saveRDS(ca.dep, file = "ca_dependent_genes.rds")

DoHeatmap(b.only.avg, features = tca.genes, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/beta_tca_heatmap.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/beta_tca_heatmap.eps"), width = 11, height = 8.5)

DoHeatmap(b.only.avg, features = ca.dep, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/beta_ca_dep_heatmap.pdf"), width = 11, height = 8.5)

DoHeatmap(b.only.avg, features = ox.phos.genes, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/beta_oxphos_heatmap.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/beta_oxphos_heatmap.eps"), width = 11, height = 8.5)
# not found COX6B, ATP5O, ATP5L, ATP5J2, ATP5J, ATP5I, ATP5H, ATP5G3, ATP5G2, ATP5G1, ATP5F1, ATP5E, ATP5D, ATP5C1, ATP5B, ATP5A1

DoHeatmap(b.only.avg, features = fig.genes, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/beta_fig_heatmap.pdf"), width = 11, height = 11)
ggsave(file=paste0("seurat/endo_cells/beta_fig_heatmap.eps"), width = 11, height = 11)

DoHeatmap(b.only.avg, features = c("INS", "IAPP", "UCN3", "MAFA", "MAFB", "ERO1B", "PDX1", "NKX6-1",
                                   "NKX2-2", "CPE", "PCSK1", "PCSK2", "PCSK1N"), 
          assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/beta_markers_fin.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/beta_markers_fin.eps"), width = 11, height = 8.5)

DoHeatmap(b.only.avg, features = c("INS", "NKX6-1", "PDX1", "GCK", "ATP2A2", "UCN3", "G6PC2", "ERO1B"), 
          assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/beta_hub_virgin_hm.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/beta_hub_virgin_hm.eps"), width = 11, height = 8.5)

a.only.avg <- AverageExpression(a_only, return.seurat = T, verbose = F)
DoHeatmap(a.only.avg, features = c("GCG", "ARX", "LOXL4", "IRX2", "TTR", "POU3F4", "GLRB"), assay='SCT', draw.lines = F) + 
  scale_fill_gradientn(colors = heatmap_col) + theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/alpha_marker_fin.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/alpha_marker_fin.eps"), width = 11, height = 8.5)

DoHeatmap(a.only.avg, features = tca.genes, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/alpha_tca_heatmap.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/alpha_tca_heatmap.eps"), width = 11, height = 8.5)

DoHeatmap(a.only.avg, features = ox.phos.genes, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/alpha_oxphos_heatmap.pdf"), width = 11, height = 8.5)
ggsave(file=paste0("seurat/endo_cells/alpha_oxphos_heatmap.eps"), width = 11, height = 8.5)

# only B1, B3, B5
#sub.b.so <- subset(b.plot, idents = c("B1", "B3", "B5"))
Idents(b.plot) <- "plot_id"
levels(b.plot) <- sort(levels(b.plot))
b.cond.avg <- AverageExpression(b.plot, return.seurat = T, verbose = F)
sub.b.deg <- subset(b_cond_DEG, (cluster_num == "B1" | cluster_num == "B2" | cluster_num == "B3" |
                                   cluster_num == "B4" | cluster_num == "B5") & cluster != "Negative")
#sub.b.unique <- arrange(sub.b.deg, cluster) %>%distinct(gene)
sub.b.unique <- sub.b.deg %>% arrange(cluster) %>% distinct(gene)
sub.b.cond.avg <- subset(b.cond.avg, subset = orig.ident == "B1" | orig.ident == "B2" |
                           orig.ident == "B3" | orig.ident == "B4" | orig.ident == "B5")
DoHeatmap(sub.b.cond.avg, features = sub.b.unique$gene, assay='SCT', draw.lines = F) + 
  scale_fill_gradientn(colors = heatmap_col) + theme(legend.position="bottom") + 
  FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file=paste0(loc, "beta_only_cond_deg.pdf"), width = 11, height = 11)
ggsave(file=paste0(loc, "beta_only_cond_deg.eps"), width = 11, height = 11)

# DoHeatmap(b.cond.avg, features = ca.dep, assay='SCT', draw.lines = F) + 
#   scale_fill_gradientn(colors = heatmap_col) + theme(legend.position="bottom") + 
#   FontSize(x.text = 10, y.text = 12) + RotatedAxis()
# ggsave(file=paste0(loc, "beta_135_heatmap.pdf"), width = 11, height = 8.5)
# ggsave(file=paste0(loc, "beta_135_heatmap.eps"), width = 11, height = 8.5)

VlnPlot(sub.b.so, features =  c("IER3", "FOS", "ZNF331", "BTG2"), assay='SCT', pt.size = 0, ncol = 1)
ggsave(file=paste0(loc, "beta_135_vlnplot.pdf"), width = 8.5, height = 12)
ggsave(file=paste0(loc, "beta_135_vlnplot.eps"), width = 8.5, height = 12)

DoHeatmap(sub.b.cond.avg, features = ca.dep, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/beta_ca_dep_heatmap.pdf"), width = 11, height = 11)
ggsave(file=paste0("seurat/endo_cells/beta_ca_dep_heatmap.eps"), width = 11, height = 11)

# Check housekeeping genes
housekeeping <- c("RPL13A", "YWHAZ", "TUBB", "GAPDH", "C1ORF43", "CHMP2A", "EMCJ", "GPI",
                  "PSMB2", "PSMB4", "RABJA", "REEP5", "SNRPD3", "VCP", "VPS29") # "ACTB"
VlnPlot(endo, features = housekeeping, assay='SCT', group.by = "cell_name",
        pt.size = 0, ncol = 1)
ggsave(file=paste0("seurat/endo_cells/endo_vlnplot_housekeeping2.pdf"), width = 8.5, height = 30)

DoHeatmap(b.avg, group.by = "condition_id", features = housekeeping, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file=paste0("seurat/endo_cells/endo_hm_housekeeping.pdf"), width = 8.5, height = 11)

# Visualizing sample_id and condition_id for each clusters
vln <- VlnPlot(endo, features = c("nFeature_RNA", "nCount_RNA"), assay = 'SCT', group.by = "plot_id",
               pt.size = 0, ncol = 1)
ggsave(vln, file=paste0("seurat/endo_cells/endo_vlnplot_feat_count_1.pdf"), width = 15, height = 11)

# plot bargrapg of samples in each clusters
Idents(endo) <- "plot_id"
levels(endo) <- sort(levels(endo))
total.cells <- as.data.frame(table(Idents(endo)))
# ggplot(data = total.cells, aes(x = Var1, y = Freq)) + geom_bar(stat="identity") + RotatedAxis()
# ggsave(file=paste0("seurat/endo_cells/endo_bar_feat_count.pdf"), width = 15, height = 11)

ggplot(endo@meta.data, aes(x=cell_name, fill=sample_id)) + geom_bar() + RotatedAxis() # + scale_fill_manual(values = c("Green","tomato"))
ggsave(file=paste0("seurat/endo/endo_bar_sample_count.pdf"), width = 8.5, height = 11)

boxplt <- vln + geom_boxplot()


VlnPlot(endo, features = c("CREB1", "ATF4", "CREB2", "CREB3", "CREB5", "CREB3L1", "CREB3L2", "CREB3L3", "CREB3L4", 
                           "NFATC1", "NFATC2", "NFATC3", "NFATC4", "NFATC5", "NFAT5"), 
        assay='SCT', group.by = "cell_name", pt.size = 0, ncol = 1)
ggsave(file=paste0("seurat/endo_cells/endo_vlnplot_some_genes.pdf"), width = 8.5, height = 25)


VlnPlot(endo, features = c("SCN1B", "SCN2A", "SCN3A", "SCN3B", "SCNM1", "SCN9A", "KCNJ8", "KCNJ11", "ABCC8", "ABCC9", 
                           "CACNA1A", "CACNA1C", "CACNA1D", "CACNA1H", "CACNA2D1", "CACNA2D2", "CACNB2"), 
        assay='SCT', group.by = "cell_name", pt.size = 0, ncol = 1)
ggsave(file=paste0("seurat/endo_cells/endo_vlnplot_ion_channel.pdf"), width = 8.5, height = 30)

VlnPlot(endo, features = c("SCN1B", "ABCC8", "CACNA1A", "CACNA1D", "CACNA2D1"), 
       assay='SCT', group.by = "cell_name", pt.size = 0, ncol = 1)
ggsave(file=paste0("seurat/endo_cells/endo_vlnplot_ion_channel_1.pdf"), width = 8.5, height = 20)

ggsave(FeaturePlot(endo, features = "SCN9A", reduction = 'tsne'), file="seurat/endo_cells/endo_SCN9A_feat_plot.pdf", width = 11, height = 8.5)

endo.avg <- AverageExpression(endo, assay = "SCT", return.seurat = T, verbose = F)
DoHeatmap(endo.avg, features = SCN9A, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12)
ggsave(file=paste0("seurat/endo_cells/endo_SCN9A_heatmap.pdf"), width = 11, height = 11)

b.only.cond.avg <- AverageExpression(b_only, assays = "SCT", return.seurat = T, verbose = F)
a.only.cond.avg <- AverageExpression(a_only, assays = "SCT", return.seurat = T, verbose = F)
d.cond.avg <- AverageExpression(d.plot, assays = "SCT", return.seurat = T, verbose = F)
DoHeatmap(b.only.cond.avg, features = "SCN9A", assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file=paste0(loc, "beta_only_SCN9A_heatmap.pdf"), width = 11, height = 11)

DoHeatmap(a.only.cond.avg, features = "SCN9A", assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file=paste0(loc, "alpha_only_SCN9A_heatmap.pdf"), width = 11, height = 11)

DoHeatmap(d.cond.avg, features = "SCN9A", assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file=paste0(loc, "delta_SCN9A_heatmap.pdf"), width = 11, height = 11)

DotPlot(b_only, features = "SCN9A", cols = col_scheme)
ggsave(file=paste0(loc, "beta_SCN9A_dotplot.pdf"), width = 8.5, height = 11)

DotPlot(a_only, features = "SCN9A", cols = col_scheme)
ggsave(file=paste0(loc, "alpha_SCN9A_dotplot.pdf"), width = 8.5, height = 11)

glycolysis <- c("SLC2A1", "SLC2A2", "SLC2A3", "SLC2A4", "SLC2A5", "HK1", "HK2", "HK3", "GCK", "GPI", "PFKM",
                "PFKL", "PFKP", "ALDOA", "ALDOB", "ALDOC", "GAPDH", "GAPDHS", "PGK1", "PGK2", "PGAM1", "PGAM2",
                "ENO1", "ENO2", "ENO3", "PKM2", "PKLR", "MPC1", "MPC2", "PDHA1", "PDHA2", "PDHB", "PDHX", "DLAT",
                "DLD")
py_lac_conv <- c("LDHA", "LDHB", "LDHC", "LDHAL6B")
make_glu <- c("PCK1", "PC", "G6PC", "G6PC2", "FBP1", "FBP2")

DoHeatmap(b.only.avg, features = glycolysis, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file="seurat/endo_cells/beta_glycolysis_heatmap.pdf", width = 11, height = 11)

DoHeatmap(b.only.avg, features = py_lac_conv, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file="seurat/endo_cells/beta_pyruvate_conv_heatmap.pdf", width = 11, height = 11)

DoHeatmap(b.only.avg, features = make_glu, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file="seurat/endo_cells/beta_gluconeogenesis_heatmap.pdf", width = 11, height = 11)

DoHeatmap(a.only.avg, features = glycolysis, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file="seurat/endo_cells/alpha_glycolysis_heatmap.pdf", width = 11, height = 11)

DoHeatmap(a.only.avg, features = py_lac_conv, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file="seurat/endo_cells/alpha_pyruvate_conv_heatmap.pdf", width = 11, height = 11)

DoHeatmap(a.only.avg, features = make_glu, assay='SCT', draw.lines = F) + scale_fill_gradientn(colors = heatmap_col) + 
  theme(legend.position="bottom") + FontSize(x.text = 12, y.text = 12) + RotatedAxis()
ggsave(file="seurat/endo_cells/alpha_gluconeogenesis_heatmap.pdf", width = 11, height = 11)

FeaturePlot(endo, features = "PCSK1")
ggsave(file="seurat/endo_cells/endo_PCSK1_umap.pdf", width = 11, height = 11)

DimPlot()

b_pos <- c("NEAT1", "MT-ND3", "MT-ND6", "GADD45B", "FOS", "IER3")
b_pos2 <- c("NR4A1", "NR4A2", "NPAS4", "ZNF331", "C2CD4B", "RGS16", "BTG2")
b_basal <- c("HIST1C", "MT2A", "MT1X", "HIST1E", "ID1", "NPY", "DDIT3", "HIST1H4C")

DotPlot(b_only, features = c(b_pos, b_pos2), cols = col_scheme, dot.scale = 10) + RotatedAxis()
ggsave(file="beta_only_pos_dotplot.pdf", width = 11, height = 8.5)

