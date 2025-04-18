# This script takes in a Seurat object that has already been run through
# 2_testSoftPowers.R and uses hdWGCNA to run a WGCNA

# This script uses optparse to parse arguments for specifying fragment files, output directory, 
# and optional parameters for data processing.

# Required arguments:
# rds_file: Path to RDS file containing processed Seurat object.
# output_dir: Path to directory to output final object and intermediate files.
# power: integer power to use for network construction

# Optional arguments
# name: Name to use for the WGCNA object. Default is the base filename of the rds_file.
# seed: Random seed to use for reproducibility.

# The script will:
# 1. Read in the Seurat object that has been run through 2_testSoftPowers.R
# 2. Construct a co-expression network using hdWGCNA
# 3. Save the WGCNA object as an RDS file
# 
# Usage:
# Rscript 3_runWGCNA.R -r <rds_file> -o <output_dir> -p <power> [-a <name>] [-s <seed>]

# Argument parsing
suppressMessages(library(optparse))
option_list <- list(
    make_option(c("-r", "--rds_file"), type="character", default=NULL, help="Path to RDS file containing processed Seurat object."),
    make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Path to directory to output final object and intermediate files."),
    make_option(c("-p", "--power"), type="integer", default=NULL, help="Soft power to use for network construction."),
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
power <- opt$power
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
cat("\n")

# Ensure all required arguments are provided
if(is.null(rds_file) || is.null(output_dir) || is.null(power)) {
    stop("Missing required arguments.")
}

# Check if output directory exists, if it does already exist, throw an error
if (dir.exists(output_dir)) {
    stop(sprintf("Output directory %s already exists. Exiting", output_dir))
} else {
    cat(sprintf("Creating output directory %s\n", output_dir))
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    cat("\n")
}

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

# construct co-expression network:
cat("Constructing co-expression network\n")
setwd(output_dir)
print(paste0("Working directory to throw the TOM object: ", getwd()))
adata <- ConstructNetwork(
  adata, 
  soft_power=power,
  use_metacells=TRUE,
  setDatExpr=FALSE,
  tom_name = name # name of the topoligical overlap matrix written to disk
)
cat("\n")

# Overwrite the WGCNA object after adding the network
cat(sprintf("Saving object to %s.rds\n", out_prefix))
saveRDS(adata, file=sprintf('%s.rds', out_prefix))
cat(sprintf("Saved object to %s.rds\n", out_prefix))
cat("\n")

# plot the dendrogram of the genes in modules
options(repr.plot.width=12, repr.plot.height=12)
png(sprintf("%s_moduleDendrogram.png", out_prefix), widt=600, height=600)
PlotDendrogram(adata, main=sprintf('%s Dendrogram', name))
dev.off()

# Save the module sizes to a dataframe
write.table(t(data.frame(table(get(name, adata@misc)$wgcna_net$colors))), sprintf("%s_module_sizes.tsv", out_prefix), sep="\t")
cat("\n")
