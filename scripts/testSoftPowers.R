# This script takes in a Seurat object that has already been run through
# 1_prepMetacells.R and uses hdWGCNA to find a soft power threshold to use

# This script uses optparse to parse arguments for specifying fragment files, output directory, 
# and optional parameters for data processing.

# Required arguments:
# rds_file: Path to RDS file containing processed Seurat object.
# output_dir: Path to directory to output final object and intermediate files.
# group_by: Column of the metadata with groups in it. Do not include white space between the commas.
#           Column must have been used in making metacells
# groups: Comma-separated list of groups conatined within the group_by parameter to find soft-thresholds for.

# Optional arguments
# name: Name to use for the WGCNA object. Default is the base filename of the rds_file.
# seed: Random seed to use for reproducibility.

# The script will:
# 1. Read in the Seurat object that must have been processed in 1_prepMetacells.R
# 2. Set the expression matrix based on the group_by and groups parameters
# 3. Run testSoftPowers to find a soft power threshold
# 4. Plot the results
# 5. Save the soft power threshold results
# 6. Overwrite the WGCNA object with the soft power threshold included

# Usage:
# Rscript 2_testSoftPowers.R -r <rds_file> -o <output_dir> -b <group_by> -g <groups> -s <seed>

# Argument parsing
suppressMessages(library(optparse))
option_list <- list(
    make_option(c("-r", "--rds_file"), type="character", default=NULL, help="Path to RDS file containing processed Seurat object."),
    make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Path to directory to output final object and intermediate files."),
    make_option(c("-b", "--group_by"), type="character", default=NULL, help="A variable in the metadata that contains groups. By default, all cells will be used"),
    make_option(c("-g", "--groups"), type="character", default=NULL, 
        help="Comma-separated list of groups conatined within the group_by parameter to find soft-thresholds for. Cannot contain whitespace between commas!"),
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
group_by <- opt$group_by
groups <- strsplit(opt$groups, ",")[[1]]
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
if(is.null(rds_file) || is.null(output_dir) || is.null(group_by) || is.null(groups)) {
    stop("Missing required arguments.")
}

# Ensure that groups is set if groupby is set
if(!is.null(group_by) && is.null(groups)) {
    stop("If group_by is set, groups must be set.")
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
SetActiveWGCNA(adata) <- name

# transpose the matrix, taking care of the
cat("Setting the expression matrix for running testSoftPowers\n")
adata <- SetDatExpr(
    adata,
    group.by=group_by,
    group_name=groups,
    assay="RNA",
    use_metacells=TRUE,
    wgcna_name=name,
    slot="data"
)
head(get(name, adata@misc)$datExpr)[1:5, c("LINC01409", "LINC01128", "SAMD11")]
cat(sprintf("Dimensions of the matrix are %s\n", dim(get(name, adata@misc)$datExpr)))
cat("\n")

# Test different soft powers:
cat("Running testSoftPowers\n")
adata <- TestSoftPowers(
  adata,
  use_metacells = TRUE,  # this is the default, I'm just being explicit
  setDatExpr = FALSE  # set this to FALSE since we did this above
)
cat("\n")

# plot the results:
cat("Plotting the soft power threshold results\n")
options(repr.plot.width=12, repr.plot.height=12)
png(sprintf("%s_softThreshold.png", out_prefix), widt=600, height=600)
plot_list <- PlotSoftPowers(adata)
wrap_plots(plot_list, ncol=2)
dev.off()
cat("\n")

# Save the soft power threshold results
power_table <- GetPowerTable(adata)
power <- power_table$Power[which(power_table$SFT.R.sq > 0.85)[1]]
cat(sprintf("Automatic soft power %s\n", power))
write.table(power_table, sprintf("%s_powerTable.tsv", out_prefix), sep="\t")

# Overwrite the WGCNA object with the soft power threshold included
saveRDS(adata, file=sprintf('%s.rds', out_prefix))
