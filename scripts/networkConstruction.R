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

# Set params
ASSAY <- "RNA"

# Get arguments from command line
args = commandArgs(trailingOnly=TRUE)
RDS <- args[1] # "/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/auxiliary_data/snrna/adrenal_Parse_10x_integrated.rds"
NAME <- args[2]  # "mouse_adrenal"
OUT <- args[3] # "/cellar/users/aklie/projects/igvf/topic_grn_links/grn_inference/hdwgcna/results/mouse_adrenal/"
OUT <- file.path(OUT, NAME)
POWER <- as.integer(args[4]) # 10

# Read in the R object
print(paste0("Loading ", RDS, "..."))
seurat_obj <- readRDS(RDS)
DefaultAssay(seurat_obj) <- ASSAY

# Construct co-expression network:
print(paste0("Constructing networks..."))
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power=POWER,
  use_metacells=TRUE,
  setDatExpr=FALSE,
  tom_name=NAME # name of the topoligical overlap matrix written to disk
)

# Plot the dendrogram of the genes in modules
options(repr.plot.width=12, repr.plot.height=12)
png(sprintf("%s_moduleDendrogram.png", OUT), widt=600, height=600)
PlotDendrogram(seurat_obj, main=sprintf('%s Dendrogram', NAME))
dev.off()

# Store the TOM object for later
#TOM <- GetTOM(seurat_obj)
#write.table(TOM, sprintf("%s_TOM.tsv", OUT), sep="\t")
     
# Compute all MEs in the full single-cell dataset
print(paste0("Calculating MEs..."))
seurat_obj <- ModuleEigengenes(seurat_obj, assay=ASSAY, verbose=FALSE)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
write.table(MEs, sprintf("%s_MEs.tsv", OUT), sep="\t")

# Compute eigengene-based connectivity (kME):
print(paste0("Calculating kME..."))
seurat_obj <- ModuleConnectivity(
    seurat_obj,
    assay=ASSAY,
    slot="data",
    harmonized=FALSE
)

# Rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = sprintf("%s-M", NAME)
)

# Plot genes ranked by kME for each module
options(repr.plot.width=12, repr.plot.height=12)
png(sprintf("%s_topGenesPerModule.png", OUT), widt=1200, height=1200)
p <- PlotKMEs(seurat_obj, ncol=5)
dev.off()

# Get the module assignment table:
modules <- GetModules(seurat_obj)
write.table(modules, sprintf("%s_modules.tsv", OUT), sep="\t")

# Save the fully processed Seurat object to be used in all the other notebooks
print(paste0("Saving to ", paste0(OUT, "_hdWGCNA.rds")))
saveRDS(seurat_obj, file=paste0(OUT, "_hdWGCNA.rds"))
