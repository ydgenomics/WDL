### Date: 251011
### /home/stereonote/miniconda3/envs/r_env/bin/R

library(Seurat)
library(dplyr)
library(ggplot2)
library(CHOIR)
library(clustree)


library(optparse)
option_list <- list(
  make_option(c("-f", "--file_paths"),
              type = "character", default = "/data/input/Files/P-Ref/Gossypium_hirsutum/SC/01_dataget/Gh/CCRI12.hr.rds|/data/input/Files/P-Ref/Gossypium_hirsutum/SC/01_dataget/Gh/CCRI12gl.hr.rds",
              help = "Comma-separated RDS file paths"),
  make_option(c("-k", "--cluster_key"),
              type = "character", default = "leiden_res_0.50|leiden_res_0.50",
              help = "Meta.data column name for subsetting [default %default]"),
  make_option(c("-v", "--cluster_value"),
              type = "character", default = "all|all",
              help = "Comma-separated cluster values to keep [default %default]"),
  make_option(c("-p", "--plot_keys"),
              type = "character", default = "biosample|sample|leiden_res_0.50",
              help = "Comma-separated keys for UMAP plotting [default %default]"),
  make_option(c("-r", "--r_value"),
              type = "numeric", default = 0.5,
              help = "Clustering resolution [default %default]"),
  make_option(c("-n", "--name"),
              type = "character", default = "Gh",
              help = "Prefix for output files [default %default]")
)

## 2. 解析 --------------------------------------------------------------
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## 3. 取出并拆分 --------------------------------------------------------
file_paths  <- unlist(strsplit(opt$file_paths,  split = "|", fixed = TRUE))
cluster_key <- unlist(strsplit(opt$cluster_key,    split = "|", fixed = TRUE))
cluster_value <- unlist(strsplit(opt$cluster_value, split = "|", fixed = TRUE))
plot_keys   <- unlist(strsplit(opt$plot_keys,   split = "|", fixed = TRUE))
r_value     <- opt$r_value
name        <- opt$name



print(file_paths)
print(cluster_key)
print(cluster_value)
plot_keys <- c(plot_keys,paste0("RNA_snn_res.",as.character(r_value)))
plot_keys <- c(plot_keys,'CHOIR_clusters_0.05')
print(plot_keys)
print(r_value)
print(name)


seurat_pipeline <- function(seu, r_value, plot_keys, prefix){
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, nfeatures = 3000)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu, features = VariableFeatures(object = seu), verbose = FALSE)
    seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
    seu <- FindClusters(seu, resolution = as.numeric(r_value)) # RNA_snn_res.0.5
    seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
    print("--- CHOIR ---")
    seu <- seu |>
      buildTree() |>
      pruneTree()
    seu

    head(seu@meta.data)
    seu <- runCHOIRumap(seu, reduction = "P0_reduction")
    # plotCHOIR(seu,
    #           accuracy_scores = FALSE,
    #           plot_nearest = FALSE)
    # dev.off()
    # result <- getRecords(object, key = "CHOIR")
    # seu$leiden_res_choir <- seu$CHOIR_clusters_0.05
    # seu@meta.data[[paste0("leiden_res_r",as.character(r_value))]] <- seu@meta.data[[paste0("RNA_snn_res.",as.character(r_value))]]
    # clustree(object@meta.data, prefix="leiden_res_")
    prefix <- sub("\\.rds$", "", basename(prefix))
    pdf(paste0(prefix, "_umap.pdf"), width = 10, height = 8)
    for (plot_key in plot_keys){
        p1 <- DimPlot(seu, reduction = "umap", group.by = plot_key, shuffle = TRUE, label = TRUE)
        print(p1)
    }
    dev.off()
    saveRDS(seu, paste0(prefix,".rds"))
    return(seu)
}

merged_data <- readRDS(file_paths[[1]]); DefaultAssay(merged_data) <- "RNA"
if (!unlist(strsplit(cluster_value[[1]], split = ","))[1] == "all"){
    print(paste0("Dealing: ", cluster_key[[1]], " & ", cluster_value[[1]]))
    merged_data <- subset(merged_data, subset = !!sym(cluster_key[[1]]) == unlist(strsplit(cluster_value[[1]], split = ",")))
}
merged_data <- seurat_pipeline(merged_data, r_value, plot_keys, prefix=basename(file_paths[[1]]))

if (length(file_paths) > 1) {
    for (i in 2:length(file_paths)) {
        temp_data <- readRDS(file_paths[[i]]); DefaultAssay(temp_data) <- "RNA"
        if (!unlist(strsplit(cluster_value[[i]], split = ","))[1] == "all"){
            print(paste0("Dealing: ", file_paths[[i]]))
            temp_data <- subset(temp_data, subset = !!sym(cluster_key[[i]]) == unlist(strsplit(cluster_value[[i]], split = ",")))
        }
        temp_data <- seurat_pipeline(temp_data, r_value, plot_keys, prefix=basename(file_paths[[i]]))
        merged_data <- merge(merged_data, temp_data)
    }
}

merged_data <- NormalizeData(merged_data)
merged_data <- FindVariableFeatures(merged_data, nfeatures = 3000)
merged_data <- ScaleData(merged_data)
merged_data <- RunPCA(merged_data, features = VariableFeatures(object = seu), verbose = FALSE)
merged_data <- FindNeighbors(merged_data, reduction = "pca", dims = 1:30)
# merged_data <- FindClusters(merged_data, resolution = as.numeric(r_value)) # RNA_snn_res.0.5
merged_data <- RunUMAP(merged_data, dims = 1:20, verbose = FALSE)
pdf(paste0(name, "_umap.pdf"), width = 10, height = 8)
for (plot_key in plot_keys){
    p1 <- DimPlot(seu, reduction = "umap", group.by = plot_key, shuffle = TRUE, label = TRUE)
    print(p1)
}
dev.off()

saveRDS(merged_data, paste0(name,".rds"))