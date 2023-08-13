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
NN <- 25  # number of nearest neighbors for meta cell construction
TARGET_CELLS = 1000  # max cells to include for one group
GROUPS <- c("celltypes", "sample") # grouping for metacells, first one will be used as idents

# Get arguments from command line
args = commandArgs(trailingOnly=TRUE)
RDS <- args[1] # "/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/auxiliary_data/snrna/adrenal_Parse_10x_integrated.rds"
NAME <- args[2]  # "mouse_adrenal"
OUT <- args[3] # "/cellar/users/aklie/projects/igvf/topic_grn_links/grn_inference/hdwgcna/results/mouse_adrenal"

# Read in the R object
print(paste0("Loading ", RDS, " ..."))
seurat_obj <- readRDS(RDS)
DefaultAssay(seurat_obj) <- ASSAY

# Set-up a Seurat object for WGCNA
print(paste0("Setting up for WGCNA..."))
seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = NAME # the name of the hdWGCNA experiment
)

# Construct metacells n each group
print(paste0("Creating metacells with ", NN, " ..."))
seurat_obj <- MetacellsByGroups(
  seurat_obj=seurat_obj,
  group.by=GROUPS, # specify the columns in adata@meta.data to group by
  k=NN, # nearest-neighbors parameter
  max_shared=10, # maximum number of shared cells between two metacells
  ident.group = GROUPS[1], # set the Idents of the metacell seurat object
  assay=ASSAY,
  slot="counts",
  target_metacells=TARGET_CELLS,
)

# Transpose the matrix
print(paste0("Setting up expression data..."))
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

# Plot the results:
options(repr.plot.width=12, repr.plot.height=12)
png(sprintf("%s_softThreshold.png", OUT), widt=600, height=600)
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)
dev.off()

# Get the power table, can also access with head(get(NAME, seurat_obj@misc)$wgcna_powerTable)
power_table <- GetPowerTable(seurat_obj)
power <- power_table$Power[which(power_table$SFT.R.sq > 0.85)][1]
print(paste0("Constructing networks with power ", power, "..."))

# Construct co-expression network:
print(paste0("Constructing networks..."))
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power=power,
  use_metacells=TRUE,
  setDatExpr=FALSE,
  tom_name=NAME # name of the topoligical overlap matrix written to disk
)

# plot the dendrogram of the genes in modules
options(repr.plot.width=12, repr.plot.height=12)
png(sprintf("%s_moduleDendrogram.png", OUT), widt=600, height=600)
PlotDendrogram(seurat_obj, main=sprintf('%s Dendrogram', NAME))
dev.off()

# store the TOM object for later
TOM <- GetTOM(seurat_obj)
write.table(TOM, sprintf("%s_TOM.tsv", OUT), sep="\t")
     
# Compute all MEs in the full single-cell dataset
print(paste0("Calculating MEs..."))
seurat_obj <- ModuleEigengenes(seurat_obj, assay=ASSAY, verbose=FALSE)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
write.table(MEs, sprintf("%s_MEs.tsv", OUT), sep="\t")


# compute eigengene-based connectivity (kME):
print(paste0("Calculating kME..."))
seurat_obj <- ModuleConnectivity(
    seurat_obj,
    assay=ASSAY,
    slot="data",
    harmonized=FALSE
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = sprintf("%s-M", NAME)
)

# plot genes ranked by kME for each module
options(repr.plot.width=12, repr.plot.height=12)
png(sprintf("%s_topGenesPerModule.png", OUT), widt=1200, height=1200)
p <- PlotKMEs(seurat_obj, ncol=5)
dev.off()

# get the module assignment table:
modules <- GetModules(seurat_obj)
write.table(modules, sprintf("%s_modules.tsv", OUT), sep="\t")
     
# Save the fully prepped Seurat object to be used in all the other notebooks
print(paste0("Saving to ", paste0(OUT, "_hdWGCNA.rds")))
saveRDS(seurat_obj, file=paste0(OUT, "_hdWGCNA.rds"))
