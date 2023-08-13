# This script takes in a processed Seurat object and uses hdWGCNA to 
# to prepare metacells for downstream analysis.

# This script uses optparse to parse arguments for specifying fragment files, output directory, 
# and optional parameters for data processing.

# Required arguments:
# rds_file: Path to RDS file containing processed Seurat object.
# output_dir: Path to directory to output final object and intermediate files.
# group_by: Comma-separated list of variables to use for grouping cells. Do not include white space between the commas.

# Optional arguments
# nearest_neighbors: number of nearest neighbors to use for creating metacells, by default will use 25
# genes: either a float between 0 and 1 or the number of variable genes to use, by default will use 5% of genes
# reduction: dimensionality reduction to use for KNN, by default will look for "harmony" in the reductions
# name: name to use for the WGCNA object. Default is the base filename of the rds_file.
# seed: random seed to use for reproducibility, by default will use 12345

# The script will:
# 1. Load in the rds file with the processed Seurat object
# 2. Plot the umap with the first grouping variable
# 3. Set-up the object for WGCNA using either variable genes or genes expressed in at least x% of cells
# 4. Construct metacells with x nearest neighbors and grouped by x. Note that the ident group is the first grouping variable.
# 5. Normalize the metacell expression matrix
# 6. Save the metacell stats
# 7. Save the metacell object
# 8. Save the object for the next step

# Usage:
# Rscript --vanilla 1_prepMetacells.R -r <rds_file> -o <output_dir> -b <grouping_vars> -g <genes> -n <nearest_neighbors> -u <reduction> -s <seed>

# Argument parsing
suppressMessages(library(optparse))
option_list <- list(
    make_option(c("-r", "--rds_file"), type="character", default=NULL, help="Path to RDS file containing processed Seurat object."),
    make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Path to directory to output final object and intermediate files."),
    make_option(c("-b", "--group_by"), type="character", default=NULL, 
        help="Comma-separated list of variables to use for grouping cells. Do not include white space between the commas."),
    make_option(c("-n", "--nearest_neighbors"), type="integer", default=25, help="Number of nearest neighbors to use for creating metacells."),
    make_option(c("-g", "--genes"), type="numeric", default=0.05, help="Either a float between 0 and 1 or the number of variable genes to use."),
    make_option(c("-u", "--reduction"), type="character", default="harmony", help="Dimensionality reduction to use for KNN."),
    make_option(c("-a", "--name"), type="character", default=NULL, help="Name to use for the WGCNA object. Default is the base filename of the rds_file."),
    make_option(c("-s", "--seed"), type="integer", default=12345, help="Random seed to use for reproducibility.")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
print(opt)

# Grab the arguments
cat("Parsing arguments")
rds_file <- opt$rds_file
output_dir <- opt$output_dir
group_by <- strsplit(opt$group_by, ",")[[1]]
nearest_neighbors <- as.numeric(opt$nearest_neighbors)
genes <- as.numeric(opt$genes)
reduction <- opt$reduction
name <- opt$name
seed <- opt$seed
cat("\n")

# Grab a name for this experiment
if (is.null(name)) {
    name <- sprintf("%s-genes_%s-neighbors", genes, nearest_neighbors)
}
out_prefix <- file.path(output_dir, name)
cat(sprintf("Using name %s\n", name))
cat(sprintf("Grouping metacells by %s\n", paste(group_by, collapse=", ")))
cat(sprintf("Using %s nearest neighbors\n", nearest_neighbors))
cat(sprintf("Using %s genes\n", genes))
cat(sprintf("Using %s reduction\n", reduction))
cat(sprintf("Using %s seed\n", seed))
cat(sprintf("Output prefix is %s\n", out_prefix))
cat("\n")

# Ensure all required arguments are provided
if(is.null(rds_file) || is.null(output_dir) || is.null(group_by)) {
    stop("Missing required arguments.")
}

# # Make the out_dirname if it doesn't exist using an R command
cat(sprintf("Creating output directory %s\n", output_dir))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("\n")

# Imports
cat("Loading libraries\n")
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratData))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(WGCNA))
suppressMessages(library(hdWGCNA))
theme_set(theme_cowplot())
options(repr.plot.width=6, repr.plot.height=6)
set.seed(opt$seed)
cat("\n")

# Actually read it
cat(sprintf("Reading RDS file %s\n", rds_file))
adata <- readRDS(rds_file)
cat("\n")

# Plot and save umap with grouping 1
cat("Plotting umap with grouping 1\n")
png(file.path(output_dir, "umap.png"), width=6, height=6, units="in", res=300)
DimPlot(adata, group.by=group_by[1], reduction="umap.wnn", label=FALSE) + umap_theme() + ggtitle(sprintf("%s", group_by[1]))
dev.off()
cat("\n")

# Set-up the object
DefaultAssay(adata) <- "RNA"

# Set-up a Seurat object for WGCNA
cat("Setting up Seurat object for WGCNA\n")
if (is.double(genes)) {
    print(sprintf("Using genes expressed in at least %s%% of cells", genes * 100))
    adata <- SetupForWGCNA(
        adata,
        gene_select = "fraction", # the gene selection approach
        fraction = genes, # fraction of cells that a gene needs to be expressed in order to be included
        wgcna_name = name # the name of the hdWGCNA experiment
    )
} else if (genes >= 1) {
    print("Using variable genes")
    adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
    adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = genes)
    adata <- SetupForWGCNA(
        adata,
        gene_select = "variable", # the gene selection approach
        wgcna_name = name # the name of the hdWGCNA experiment
    )
}
cat(sprintf("Using %s genes\n", length(get(name, adata@misc)$wgcna_genes)))
cat("\n")

# construct metacells n each group
cat(sprintf("Constructing metacells with %s nearest neighbors and grouped by %s\n", nearest_neighbors, paste(group_by, collapse=", ")))
adata <- MetacellsByGroups(
  seurat_obj = adata,
  group.by = group_by, # specify the columns in adata@meta.data to group by
  reduction = reduction, # select the dimensionality reduction to perform KNN on
  k = nearest_neighbors, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = group_by[1], # set the Idents of the metacell seurat object
  assay="RNA",  # default
  slot="counts",  # default
)
cat("\n")

# normalize metacell expression matrix
adata <- NormalizeMetacells(adata)
cat("\n")

# Save the object for the next step!
cat(sprintf("Saving object to %s.rds\n", out_prefix))
saveRDS(adata, file=sprintf('%s.rds', out_prefix))
cat(sprintf("Saved object to %s.rds\n", out_prefix))
cat("\n")

# Save the metacell stats
write.csv(
    x=get(name, adata@misc)$wgcna_params$metacell_stats, 
    file=sprintf('%s_metacell_stats.csv', out_prefix),
    quote=FALSE, 
    row.names=FALSE
)
cat(sprintf("Saved metacell stats to %s_metacell_stats.csv\n", out_prefix))
cat("\n")

# Save the metacell object
cat(sprintf("Saving metacell object to %s_metacells.rds\n", out_prefix))
saveRDS(get(name, adata@misc)$wgcna_metacell_obj, file=sprintf('%s_metacells.rds', out_prefix))
cat(sprintf("Saved metacell object to %s_metacells.rds\n", out_prefix))
cat("\n")
