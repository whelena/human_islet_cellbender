library(scater)
library(DropletUtils)
library(Matrix)
library(Seurat)
library(stringr)
library(SingleCellExperiment)
library(mvoutlier)
library(limma)
library(ggplot2)
library(ggpmisc)

#load the dataset
raw_folder_address <- readRDS(file="10x_raw_folder_address.rds")
keyword <- readRDS(file = "keyword.rds")
scater_file_address<-readRDS(file = "scater_file_address.rds")

#convert 10x output folder into a scater object, by inputting the folder that contains the matrix.mtx file and two .tsv file.
raw_h5 <- list.files(path = raw_folder_address, pattern = "filtered.h5")
data <- Read10X_h5(paste0(raw_folder_address, raw_h5))
so_obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 1000)
scater_object <- as.SingleCellExperiment(so_obj)

# raw_h5 <- list.files(path = raw_folder_address, pattern = "filtered.h5")
# scater_object_h5 <- read10xCounts(samples = paste0(raw_folder_address, raw_h5), type = "auto")

scater_object <- calculateQCMetrics(scater_object)
#keep_cells <- scater_object$total_counts > 1000
#keep_cells <- colSums(counts(scater_object))>=1000
#scater_object<-scater_object[,keep_cells]
scater_object_colnames <- colnames(scater_object)
colnames(scater_object) <- paste(keyword,colnames(scater_object),sep="_")
rownames(scater_object)<-make.unique(rownames(scater_object),sep="_")
#extract the library id attached to cell names
#offset to the library_id is calculated based on the length of the keyword, as 
library_id <- as.integer(str_extract(scater_object_colnames,"[0-9]+"))

# load csv files in the folder.
files <- list.files(path=raw_folder_address, pattern="*.csv", full.names=F, recursive=FALSE)
#for each csv files, load the information and add it to the colData of the scater_object
for(i in files){
  name<-sub("^([^.]*).csv", "\\1",i)
  file<-as.character(read.csv(paste(raw_folder_address,i,sep=""))[,1])
  colData(scater_object)[,name]<-file[library_id]
}
# colData(scater_object)$condition <- sub("_r[1-2]", "", keyword)
# colData(scater_object)$library_id <- keyword
# colData(scater_object)$sample_id <- sub("_r[1-2].*", "", keyword)
# keep_feature <- rowSums(counts(scater_object) > 0)> 0
# keep_feature <- nexprs(scater_object, byrow=T) > 0
# scater_object_filtered <- scater_object[keep_feature, ]
isSpike(scater_object, "MT") <- grepl("^MT-", rownames(scater_object))

#calculate automatic QC metrics to filter out ones with an Median Absolute Deviation greater than 3
scater_object <- scater::calculateQCMetrics(scater_object,  use_spikes = T, exprs_value = 'counts')
#colData(scater_object)$total_features_by_counts<-colData(scater_object)$total_features_by_counts

#Current version of calculateQCMetrics does not generate automatic filter for total_counts and total_features_by_counts, use isOutlier to create these filter_ fields 
scater_object$filter_on_total_counts <- isOutlier(scater_object$total_counts, nmads=3, 
                                                           type="both", log=T)
scater_object$filter_on_total_features <- isOutlier(scater_object$total_features_by_counts, nmads=3, 
                                                             type="lower", log=F)
scater_object$filter_on_pct_counts_MT<- isOutlier(scater_object$pct_counts_MT, nmads=3, 
                                                           type="higher", log=F)

#setup manual filters of 3MAD on genes and counts
scater_object$use <- (!scater_object$filter_on_total_features & !scater_object$filter_on_total_counts&!scater_object$filter_on_pct_counts_MT)
rowData(scater_object)$use <- row.names(scater_object)
#summary(scater_object$use)
#output use filter
scater_filter_use <- summary(scater_object$use)

#run automatic filters based on default parameters without total_features_feature_control
parameters<-c("pct_counts_in_top_100_features","total_features_by_counts","log10_total_counts_endogenous","log10_total_counts_feature_control","pct_counts_feature_control")
scater_object<-normalize(scater_object)
#scater_object<-runPCA(scater_object,exprs_values = "counts",detect_outliers = T,use_coldata=T)

#scater_filter_outlier <- summary(scater_object$outlier)
#scater_object_filters<-cbind(scater_object$outlier,!scater_object$use)

#merge both manual and automatic filters, and remove genes that are expressed by less than 3 cells
#merge_manual_automatic_filters <- apply(counts(scater_object[ , colData(scater_object)$use & !colData(scater_object)$outlier]), 1, function(x) length(x[x > 1]) >= 3)
#add the merged filters to the $use column dataset
#rowData(scater_object)$use <- scater_object$use #merge_manual_automatic_filters
#output filtered genes by manual and auto outlier
#scater_filter_merged <- summary(scater_object$use & !scater_object$outlier)

#saves the above file
#saveRDS(scater_object, file=paste(scater_file_address,"scater_object.rds",sep=""))
#apply the filter with the $use settings
scater_object_qc <- scater_object[rowData(scater_object)$use,
                                           colData(scater_object)$use]

saveRDS(scater_object_qc, file=paste(scater_file_address,"scater_object_qc.rds",sep=""))
