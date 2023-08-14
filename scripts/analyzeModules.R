# This script takes in a Seurat object that has already been run through
# 3_constructNetwork.R and uses hdWGCNA analyze the co-expression modules in the network.

# This script uses optparse to parse arguments for specifying fragment files, output directory, 
# and optional parameters for data processing.

# Required arguments:
# rds_file: Path to RDS file containing processed Seurat object.
# output_dir: Path to directory to output final object and intermediate files.

# Optional arguments
# harmonize_by: Comma-separated list of metadata fields to perform harmonization of module eigengenes with, no white space between commas!
# group_by: column in seurat_obj@meta.data containing grouping info, ie clusters or celltypes
# groups: Comma-separated list name of the group(s) in harmonize_by to use for kME calculation
# name: Name to use for the WGCNA object. Default is the base filename of the rds_file.
# seed: Random seed to use for reproducibility.

# The script will:
# 1. Read in the Seurat object
# 2. Compute module eigengenes
# 3. Compute kMEs
# 4. Write out the module eigengenes and kMEs to a file
# 5. Write out the Seurat object that is ready for downstream analysis

# Usage:
# Rscript 4_analyzeModules.R -r <rds_file> -o <output_dir> -b <harmonize_by> -g <harmonize_by> -a <name> -s <seed>

# Argument parsing
suppressMessages(library(optparse))
option_list <- list(
    make_option(c("-r", "--rds_file"), type="character", default=NULL, help="Path to RDS file containing processed Seurat object."),
    make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Path to directory to output final object and intermediate files."),
    make_option(c("-q", "--harmonize_by"), type="character", default=NULL, 
        help="Comma-separated list of variables to perform harmonization of module eigengenes with, no white space between commas!"),
    make_option(c("-b", "--group_by"), type="character", default=NULL, 
        help="Column in meta.data containing grouping info, ie clusters or celltypes"),
    make_option(c("-g", "--groups"), type="character", default=NULL, 
        help="Comma-separated list of groups conatined within the harmonize_by parameter to find soft-thresholds for."),
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
harmonize_by <- opt$harmonize_by
group_by <- opt$group_by
groups <- opt$groups
name <- opt$name
seed <- opt$seed
cat("\n")

# Grab a name from the base filename of the rds and strip the file extension
if (is.null(name)) {
    name <- basename(rds_file)
    name <- substr(name, 1, nchar(name) - 4)
}
out_prefix <- file.path(output_dir, name)
cat(sprintf("Using name %s\n", name))
cat(sprintf("Using output prefix %s\n", out_prefix))
if (!is.null(harmonize_by)) {
    harmonize_by <- strsplit(harmonize_by, ",")[[1]]
    cat(sprintf("Harmonizing by %s\n", paste0(harmonize_by, collapse=", ")))
} else {
    cat("Not harmonizing\n")
}
if(is.null(group_by)) {
    cat("Using all cells during kME calculation\n")
    if(!is.null(groups)) {
        stop("Cannot include groups argument if not using group_by for kME calculation")
    }
} else {
    cat(sprintf("Using only selected cells for kME calculation using %s column\n", group_by))
    if (is.null(groups)) {
        stop("Must provide groups if using group_by for kME calculation")
    }
    groups <- strsplit(groups, ",")[[1]]
    cat(sprintf("Using groups %s\n", paste0(groups, collapse=", ")))
}
cat("\n")

# Ensure all required arguments are provided
if (is.null(rds_file) | is.null(output_dir)) {
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

# Set-up the object
DefaultAssay(adata) <- "RNA"
adata <- SetActiveWGCNA(adata, wgcna_name=name)

# Scaling data
cat("Scaling data\n")
if (is.null(VariableFeatures(adata))) {
    adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
}
adata <- ScaleData(adata, features=VariableFeatures(adata))
cat("\n")

# Analyze these modules
cat("Computing module eigengenes\n")
adata <- ModuleEigengenes(
    adata, 
    group.by.vars=harmonize_by,
    verbose=TRUE,
    assay="RNA",
    wgcna_name=name
)
if (!is.null(harmonize_by)) {
    hMEs <- GetMEs(adata, harmonized=TRUE)
    write.table(hMEs, sprintf("%s_hMEs.tsv", out_prefix), sep="\t")
    cat("Computing harmonized kMEs\n")
    adata <- ModuleConnectivity(
        adata,
        group.by = group_by,
        group_name = groups,
        assay="RNA",
        slot="data",
        harmonized=TRUE,
        wgcna_name=name
    )
} else {
    cat("Computing kMEs")
    adata <- ModuleConnectivity(
        adata,
        group.by = group_by,
        group_name = groups,
        assay="RNA",
        slot="data",
        harmonized=FALSE,
        wgcna_name=name
    )
}
MEs <- GetMEs(adata, harmonized=FALSE)
write.table(MEs, sprintf("%s_MEs.tsv", out_prefix), sep="\t")
cat("\n")

# Overwrite the WGCNA object with the module fun included
cat(sprintf("Saving object to %s.rds\n", out_prefix))
saveRDS(adata, file=sprintf('%s.rds', out_prefix))
cat(sprintf("Saved object to %s.rds\n", out_prefix))
cat("\n")
