print("Loading packages...")

# Conversion libraries and Seurat
library(SeuratDisk)
library(SeuratData)
library(Seurat)
library(Signac)

# Plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# Co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# Using the cowplot theme for ggplot
theme_set(theme_cowplot())

# Set random seed for reproducibility
set.seed(12345)

# Set parameters
ASSAY <- "RNA"  # assay to use
NN <- 10  # number of nearest neighbors for meta cell construction
TARGET_CELLS = 2000  # max cells to include for one group
#GROUPS <- c("celltypes", "sample") # grouping for metacells, first one will be used as idents#
GROUPS <- c("celltypes")
#GROUPS <- c("Gene")

# Get arguments from command line
args = commandArgs(trailingOnly=TRUE)
RDS <- args[1] # "/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/auxiliary_data/snrna/adrenal_Parse_10x_integrated.rds"
NAME <- args[2]  # "mouse_adrenal"
NAME <- paste0(NAME, "_", NN, "_", TARGET_CELLS)
OUT <- args[3] # "/cellar/users/aklie/projects/igvf/topic_grn_links/grn_inference/hdwgcna/results/mouse_adrenal/"
OUT <- file.path(OUT, NAME)

# Read in the R object
print(paste0("Loading ", RDS, "..."))
seurat_obj <- readRDS(RDS)
DefaultAssay(seurat_obj) <- ASSAY
print(paste0("Matrix is ", dim(seurat_obj)[1], " genes x ", dim(seurat_obj)[2], " cells..."))

# Set-up a Seurat object for WGCNA
print(paste0("Setting up for WGCNA..."))
#seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "variable",
#    gene_select = "fraction", # the gene selection approach
#    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = NAME # the name of the hdWGCNA experiment
)

# Construct metacells n each group
print(paste0("Creating metacells with ", NN, " neighbors..."))
seurat_obj <- MetacellsByGroups(
  seurat_obj=seurat_obj,
  group.by=GROUPS, # specify the columns in adata@meta.data to group by
  k=NN, # nearest-neighbors parameter
  max_shared=10, # maximum number of shared cells between two metacells
  ident.group = GROUPS[1], # set the Idents of the metacell seurat object
  assay=ASSAY,
  min_cells=10,
  slot="counts",
  target_metacells=TARGET_CELLS,
)

# Normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

# Transpose the matrix
print(paste0("Setting up metacell expression matrix..."))
seurat_obj <- SetDatExpr(
    seurat_obj, 
    assay="RNA", 
    use_metacells=TRUE, 
    wgcna_name=NAME, 
    slot="data"
)

# Test different soft powers
print(paste0("Testing soft powers..."))
seurat_obj <- TestSoftPowers(
  seurat_obj,
  use_metacells=TRUE,  # this is the default, I'm just being explicit
  setDatExpr=FALSE  # set this to FALSE since we did this above
)

# plot the results:
options(repr.plot.width=12, repr.plot.height=12)
png(sprintf("%s_softThreshold.png", OUT), widt=600, height=600)
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)
dev.off()

# Get the power table, can also access with head(get(NAME, seurat_obj@misc)$wgcna_powerTable)
power_table <- GetPowerTable(seurat_obj)
write.table(power_table, sprintf("%s_powerTable.tsv", OUT), sep="\t")

# Save the fully prepped Seurat object to be used in all the other notebooks
print(paste0("Saving to ", paste0(OUT, "_prelim_hdWGCNA.rds")))
saveRDS(seurat_obj, file = paste0(OUT, "_prelim_hdWGCNA.rds"))
