# https://themilolab.github.io/SPATA2/
# # Install packages -----------------
# install.packages("devtools")
# 
# if (!base::requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# }
# 
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'EBImage'))
# 
# install.packages("Seurat")
# 
# devtools::install_github(repo = "kueckelj/confuns")
# devtools::install_github(repo = "theMILOlab/SPATA2")
# 
# # if you want to use monocle3 related wrappers
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')
# # Default installation misses hdf5r which is required to read the filtered feature matrix
# install.packages("hdf5r")
# install.packages('pracma')
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
# install.packages("viridis")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GSVA")
# install.packages("viridisLite")
# install.packages("msigdbr")
# devtools::install_github('xzhoulab/SPARK')
# 
# Declare constants ---------
library(SPATA2)
library(magrittr)
library(tidyverse)
library(monocle3)
library(viridisLite)

# MILO Lab samples:
#directory_10X = "C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/242_T"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/243_T"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/248_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/251_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/255_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/259_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/260_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/262_T"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/275_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/296_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/304_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/334_T"

#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__56_A1" #Need to reverse slashes when copied in Windows
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/50D" #Need to reverse slashes when copied in Windows
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/16_1" #Need to reverse slashes when copied in Windows
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/45B"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/47"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/A1 - 15A"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/B1(missing currently) - 19"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/C1 - 33B"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/D1_-_35"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__60_B1"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__61_C1"
directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__63_D1"


sample_name <- substr(directory_10X, 67, nchar(directory_10X))
#sample_name <- substr(directory_10X, 55, nchar(directory_10X))
subdir='SPATA2_figures'
sig_trajectory_subdir='signature_trajectories'
custom_fig_subdir='custom_figures'
directory_genesets="C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/all_SPATA2_external1.csv"
external_geneset=read.csv(directory_genesets)
#Remove empty columns from external_geneset
empty_columns <- colSums(is.na(external_geneset) | external_geneset == "") == nrow(external_geneset)
external_geneset=external_geneset[, !empty_columns]

res=100 #Resolution for figures
# Create output directory ------------
dir.create(file.path(directory_10X, subdir), showWarnings = FALSE) #create base figure directory
setwd(file.path(directory_10X, subdir))
dir.create(sig_trajectory_subdir, showWarnings = FALSE) #create subdirectory for iterative sig plotting
dir.create(custom_fig_subdir, showWarnings = FALSE) #create subdirectory for iterative sig plotting

# Create directory to read/write spata objects
spata_obj_dir=paste(getwd(),"spata_obj.RDS",sep="/")
spata_obj_dir2=paste(getwd(),"spata_obj2.RDS",sep="/")
#spata_obj_dir = "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/50D/SPATA2_figures/spata_obj.RDS"
# # Create spata object -----
# #Note, this also runs various preprocessing, DR, and clustering techniques too.
# spata_obj <-
#   initiateSpataObject_10X(
#     directory_10X = paste(directory_10X,"/outs",sep=""), # the directory from which to load the data
#     sample_name = sample_name,
#     directory_spata = spata_obj_dir
#   )
# 
# 
# Or alternatively, load premade spata obj --------
spata_obj <- loadSpataObject(spata_obj_dir, verbose = TRUE, update = TRUE)
# Adjust directory instructions for spata object ---------
spata_obj <-
  setSpataDir(
    object = spata_obj,
    #to = "spata_object",
    dir = spata_obj_dir, # the directory under which to store the object
  )
# Load external gene sets -------
#Load in a .csv where each row is a gene, and each column is a gene set name
sigs_to_add=colnames(external_geneset)
for (sig in sigs_to_add) {
  print(sig)
  genes=dplyr::pull(external_geneset, sig)
  genes <- genes[genes != "" & genes %in% c(getGenes(spata_obj))]
  spata_obj <- addGeneSet(spata_obj,class_name="external",gs_name=sig,genes=genes[genes != ""],overwrite=TRUE,check_genes = TRUE)
}

spata_obj <- discardGeneSets(spata_obj,gs_names="external_SpOT.ALLGENES")

# Denoise data ------
#This function constructs and uses a neural network to denoise expression levels spatially.
#Output is a new expression matrix called denoised
spata_obj <-
  runAutoencoderDenoising(
    object = spata_obj,
    activation = "selu",
    bottleneck = 56,
    epochs = 20,
    layers = c(128, 64, 32),
    dropout = 0.1,
    set_as_active = TRUE
  )
#Preview denoising output
spata_obj <- setActiveExpressionMatrix(object = spata_obj, mtr_name = "scaled")
plotSurface(object = spata_obj,
            color_by = "VEGFA",
            smooth = FALSE,
            pt_size = 1.8,
            display_image = TRUE) +
  ggplot2::labs(title = "Scaled (not denoised)") +
  legendNone()
ggplot2::ggsave('before_denoise.png')
spata_obj <- setActiveExpressionMatrix(object = spata_obj, mtr_name = "denoised")

plotSurface(object = spata_obj,
            color_by = "VEGFA",
            smooth = FALSE,
            pt_size = 1.8,
            display_image=TRUE) +
  ggplot2::labs(title = "Denoised")+
  legendNone()

ggplot2::ggsave('after_denoise.png')



# Infer CNV (doesn't work) -------
# directory_cnv_folder=file.path(directory_10X, subdir,"CNV")
# dir.create(directory_cnv_folder, showWarnings = FALSE)
# 
# spata_obj <-
#   runCnvAnalysis(
#     object = spata_obj,
#     directory_cnv_folder = directory_cnv_folder, # example directory
#     ref_annotation = cnv_ref$annotation,
#     ref_mtr = cnv_ref$mtr,
#     ref_regions = cnv_ref$regions
#   )
# 
# Create user-specified segmentations and trajectories------
#spata_obj <- createSegmentation(object = spata_obj)
spata_obj <- createSpatialTrajectories(object = spata_obj)

#spata_obj <- createTrajectories(object = spata_obj_subset)

# segment the "good quality area" as well as the "bad quality area"
# spata_obj <- createSegmentation(object = spata_obj)
# # visualize the segmentation
# plotSegmentation(object = spata_obj,
#                  segment_subset = "good_quality",
#                  encircle = TRUE,
#                  clrp_adjust = c("good_quality" = "forestgreen")
# )
#
# # name the segment you want to keep
# spata_obj_subset <- subsetBySegment_CountMtr(object = spata_obj, segment_name = "good_quality")
#
# plotSurface(object = spata_obj_subset, color_by = "nCount_Spatial", smooth_span = 0.1)
# # Save SpataObj----
# saveSpataObject(object = spata_obj)#Save SpataObj
# #saveCorrespondingCDS(cds = sample_cds, object = spata_obj) 
# 
# 
###---- 

# Generate misc figures ----
#for colour spectra (pt_clrsp), use "inferno" for regular colour scale; "Inferno" for reverse colour scale
# color_by='SOX4'
color_by='external_SpOT.upreg'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))

color_by='external_SpOT.upreg'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "Inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_colourRevspatial_plot.png',sep=''))

color_by='external_Ron.downreg'

plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "Inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'colourRev_spatial_plot.png',sep=''))

plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))

color_by='external_Ron.upreg'

plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))

color_by='external_Neftel.MES.like'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))

color_by='external_Neftel.MES.like'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "Inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_colourRevSpatial_plot.png',sep=''))

color_by='external_Neftel.MES2'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))

color_by='external_Neftel.MES1'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

color_by='external_Neftel.AC.like'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

color_by='external_Neftel.OPC.like'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

color_by='external_Neftel.NPC.like'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

color_by='external_Neftel.NPC1'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

color_by='external_Neftel.NPC2'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

color_by='external_SpOT.downreg'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

color_by = "HM_HYPOXIA"
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))

plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "Inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_colourRevSpatial_plot.png',sep=''))

color_by='external_Ravi.Radial.Glia'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))

color_by='external_Ravi.Reactive.Immune'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))

color_by='external_Ravi.Regional.NPC'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

color_by='external_Ravi.Regional.OPC'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

color_by='external_Ravi.Reactive.Hypoxia'
plotSurface(spata_obj,color_by=color_by,
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))

# Plot spatial trajectories-------
genes <- c("EGFR", "MBP", "SNAP25", "BCL9")
hypox_gene_sets <- c("BC_P53HYPOXIA_PATHWAY",
                  "BP.GO_REGULATION_OF_CELLULAR_RESPONSE_TO_HYPOXIA",
                  "BP.GO_HYPOXIA_INDUCIBLE_FACTOR_1ALPHA_SIGNALING_PATHWAY",
                  "BP.GO_NEGATIVE_REGULATION_OF_CELLULAR_RESPONSE_TO_HYPOXIA",
                  "BP.GO_NEGATIVE_REGULATION_OF_HYPOXIA_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
                  "BP.GO_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_HYPOXIA",
                  "RCTM_REGULATION_OF_GENE_EXPRESSION_BY_HYPOXIA_INDUCIBLE_FACTOR",
                  "RCTM_CELLULAR_RESPONSE_TO_HYPOXIA",
                  "HM_HYPOXIA")
neftel <- c("external_Neftel.MES.like","external_Neftel.OPC.like",
            "external_Neftel.AC.like","external_Neftel.NPC.like")

gene_colors <- color_vector(clrp = "npg", names = neftel)

trajectory <- 
  ggpLayerTrajectories(
    object = spata_obj, 
    ids = "low_high_HYPOXIA",
    size = 1
  )

tissue_outline <- ggpLayerTissueOutline(object = spata_obj)

plist <- 
  imap(
    .x = gene_colors, 
    .f = function(color, neftel){
      
      plotSurface(spata_obj, color_by = neftel, display_image = F) + 
        scale_color_gradient(low = alpha("white", 0), high = color) + 
        tissue_outline + 
        trajectory
      
    })
library(patchwork)
wrap_plots(plist, ncol = 2)
ggplot2::ggsave('neftel_along_traj.png')


# surface
# hist_plot_surface <-
#   plotSurface(
#     object = spata_obj,
#     color_by = "histology",
#     pt_clrp = "npg"
#   )

seurat_plot_surface <- 
  plotSurface(
    object = spata_obj,
    color_by = "seurat_clusters",
    pt_clrp = "lo"
  )

# bar plots
# hist_barplot <- 
#   plotTrajectoryBarplot(
#     object = spata_obj,
#     id = "low_high_HYPOXIA",
#     grouping_variable = "histology",
#     clrp = "npg"
#   )

seurat_barplot <- 
  plotTrajectoryBarplot(
    object = spata_obj,
    id = "low_high_HYPOXIA",
    grouping_variable = "seurat_clusters",
    clrp = "lo"
  )

# plot results
# hist_plot_surface + trajectory
# ggplot2::ggsave('histplotsurface_traj.png')
# 
# hist_barplot
# ggplot2::ggsave('hist_barplot.png') 

seurat_plot_surface + trajectory
ggplot2::ggsave('seuratplot_traj.png')

seurat_barplot
ggplot2::ggsave('seurat_barplot.png')

# Spatially Variable Genes-----
library(SPARK)
# spatially variable genes
spata_obj <- runSparkx(object = spata_obj)

spark_df <- getSparkxGeneDf(object = spata_obj, threshold_pval = 0.01)

# extract genes with a sparkx pvalue of 0.01 or lower
sparkx_genes <- spark_df[["genes"]]

STS_sparkx <- 
  spatialTrajectoryScreening(
    object = spata_obj, 
    id = "low_high_HYPOXIA", # ID of the spatial trajectory
    variables = sparkx_genes # the variables/genes to scree 
    # model_subset = model_subset 
  )

plotOverview(
  object = STS_sparkx,
  label_vars = 4, # label top 4 variables/genes per model
  label_size = 3
)
ggplot2::ggsave(paste(sig_trajectory_subdir,'/sparkx_STS.png',sep=''))

sparkxSTS_df <- 
  getResultsDf(object = STS_sparkx) %>% 
  head(200)
write.csv(sparkxSTS_df, paste0(sig_trajectory_subdir,'/sparkxSTS_results.csv'))

traj_layer <- 
  ggpLayerTrajectories(
    object = spata_obj, 
    ids = "low_high_HYPOXIA", 
    size = 2.5
  )

for(mod in unique(sparkxSTS_df[['models']])){
  genes <- unlist(subset(sparkxSTS_df, models  == mod, 
                        select=c(variables)))
  if(length(genes) < 2){
    plotSurfaceComparison(
      object = spata_obj, 
      color_by = genes
    ) + 
      traj_layer + 
      labs(subtitle = mod)
    ggplot2::ggsave(paste0(sig_trajectory_subdir,'/sigs_',mod,'.png'))
  } else{
    plotSurfaceComparison(
      object = spata_obj, 
      color_by = genes[1:3]
    ) + 
      traj_layer + 
      labs(subtitle = mod)
    ggplot2::ggsave(paste0(sig_trajectory_subdir,'/sigs_',mod,'.png'))
  }
}


# All Signature Trends Along Trajectory-----
sigs_to_plot <- getGeneSets(spata_obj, of_class = "external")

STS_sigs <- 
  spatialTrajectoryScreening(
    object = spata_obj, 
    id = "low_high_HYPOXIA", # ID of the spatial trajectory
    variables = sigs_to_plot # the variables/genes to scree 
    # model_subset = model_subset 
  )
plotOverview(
  object = STS_sigs,
  label_vars = 4, # label top 4 variables/genes per model
  label_size = 3
)
ggplot2::ggsave(paste(sig_trajectory_subdir,'/allsigs_STS.png',sep=''))

res_df <- 
  getResultsDf(STS_sigs)
write.csv(res_df, paste0(sig_trajectory_subdir,'/externalSTS_results.csv'))

for(mod in unique(res_df[['models']])){
  sets <- unlist(subset(res_df, models  == mod, 
                        select=c(variables)))
  if(length(sets) < 2){
    plotSurfaceComparison(
      object = spata_obj, 
      color_by = sets
    ) + 
      traj_layer + 
      labs(subtitle = mod)
    ggplot2::ggsave(paste0(sig_trajectory_subdir,'/sigs_',mod,'.png'))
  } else{
    plotSurfaceComparison(
      object = spata_obj, 
      color_by = sets[1:3]
    ) + 
      traj_layer + 
      labs(subtitle = mod)
    ggplot2::ggsave(paste0(sig_trajectory_subdir,'/sigs_',mod,'.png'))
  }
}


# STS for Individual Genes in a Set-----
for (sig in sigs_to_plot){
sig=sigs_to_plot[16] #TEMP PLACEHOLDER FOR DEBUGGING 1 GENE SIGNATURE
#Get a list of genes matched in our sample with an external signature
  genes_detected <- getGenes(object=spata_obj,of_gene_sets = sig)
  if (length(genes_detected)>1){
  #Here, you need to pull a list of genes picked up in our dataset from signature "sig"
  #Determine the length, and have a case or elseif statement handle it
  #Plot trajectory fit derived from upregulated genes
  
    # spatial trajectory screening
    STS_genes <- 
      spatialTrajectoryScreening(
        object = spata_obj, 
        id = "low_high_HYPOXIA", # ID of the spatial trajectory
        variables = genes_detected # the variables/genes to scree 
        # model_subset = model_subset 
      )
    
    plotOverview(
      object = STS_genes,
      label_vars = 4, # label top 4 variables/genes per model
      label_size = 3
    )
    ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig,'_STS.png',sep=''))
    
    res_df <- getResultsDf(object = STS_genes)
    res_df <- subset(res_df, sts_score > 0.5 & p_value < 0.05)
    write.csv(res_df, paste0(sig_trajectory_subdir,'/',sig,'_STSresults.csv'))
    
    for(mod in unique(res_df[['models']])){
      genes <- unlist(subset(res_df, models  == mod, 
                            select=c(variables)))
      if(length(genes) < 3){
        plotSurfaceComparison(
          object = spata_obj, 
          color_by = genes
        ) + 
          traj_layer + 
          labs(subtitle = mod)
        ggplot2::ggsave(paste0(sig_trajectory_subdir,'/',sig,'_',mod,'.png'))
      } else{
        plotSurfaceComparison(
          object = spata_obj, 
          color_by = genes[1:3]
        ) + 
          traj_layer + 
          labs(subtitle = mod)
        ggplot2::ggsave(paste0(sig_trajectory_subdir,'/',sig,'_',mod,'.png'))
      }
    }
  

  
  #Plotting trajectory heatmap will show >10 genes (signature has around 50)
  # hm_colors <- viridis::inferno(n = 100)
  # 
  # xx=plotTrajectoryHeatmap(object = spata_obj,
  #                          trajectory_name = "low_high_HYPOXIA",
  #                          variables = genes_of_interest,
  #                          arrange_rows = "maxima",
  #                          colors = hm_colors,
  #                          show_rownames = TRUE,
  #                          split_columns = FALSE,
  #                          smooth_span = 0.5,
  #                          fontsize=6)
  # png(filename=paste(sig_trajectory_subdir,'/',sig,'_heatmap.png',sep=''),res=res,width = 1024, height = 768)
  # 
  # grid::grid.newpage()
  # grid::grid.draw(xx$gtable)
  # dev.off()
  } else if (length(genes_detected)==1){
    #Plot trajectory fit derived from upregulated genes
    plotTrajectoryFit(object = spata_obj,
                      trajectory_name = "low_high_HYPOXIA",
                      variable = genes_detected,
                      display_residuals = TRUE) +
      legendTop()
    ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig,'_trajectory_trends.png',sep=''))
    
    #Color trajectory by signature
    plotTrajectory(object = spata_obj,
                   trajectory_name = "low_high_HYPOXIA",
                   color_by = genes_detected,
                   smooth = TRUE,
                   pt_alpha = 0.5,
                   display_image = TRUE) +
      legendTop()
    ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig,'_trajectory_line.png',sep=''))
    
    #Plot specific genes constituting signature along trajectory
    genes_of_interest = dplyr::pull(external_geneset, strsplit(sig,split='_')[[1]][2])
    plotTrajectoryGenes(object = spata_obj,
                        trajectory_name = "low_high_HYPOXIA",
                        genes = genes_detected,
                        smooth_span = 0.2,
                        smooth_se = TRUE,
                        display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                        nrow = 2# align the sub plots in two rows
    )
    ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig,'_examples.png',sep=''))
    
    #Generate one line plot per signature along the trajectory
    plotTrajectoryGenes(spata_obj,trajectory_name ='low_high_HYPOXIA' ,
                           genes=genes_detected,
    ) +
      ggplot2::labs(title = sig)
    ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig,'_trajectory_line_plot.png',sep=''))
    
    #Plotting trajectory heatmap will show >10 genes (signature has aroung 50)
    hm_colors <- viridis::inferno(n = 100)
    #FOR COLOUR REVERSE, USE:  hm_colors <- viridis::inferno(n = 100,direction=-1)

    xx=plotTrajectoryHeatmap(object = spata_obj,
                             trajectory_name = "low_high_HYPOXIA",
                             variables = genes_detected,
                             arrange_rows = "maxima",
                             colors = hm_colors,
                             show_rownames = TRUE,
                             split_columns = FALSE,
                             smooth_span = 0.5,
                             fontsize=6)
    png(filename=paste(sig_trajectory_subdir,'/',sig,'_heatmap.png',sep=''),res=res,width = 1024, height = 768)

    grid::grid.newpage()
    grid::grid.draw(xx$gtable)
    dev.off()
  }
}



# Save SpataObj----
saveSpataObject(object = spata_obj)#Save SpataObj
#saveCorrespondingCDS(cds = sample_cds, object = spata_obj) 


 
# Neftel Plots----
neftel <- c("external_Neftel.AC.like","external_Neftel.MES.like",
            "external_Neftel.NPC.like","external_Neftel.OPC.like")

gene_colors <- color_vector(clrp = "npg", names = neftel)
tissue_outline <- ggpLayerTissueOutline(object = spata_obj)

plist <-
  imap(
    .x = gene_colors,
    .f = function(color, neftel){

      plotSurface(spata_obj, color_by = neftel, display_image = T,
                  smooth = TRUE, 
                  #smooth_span = 0.2, 
                  pt_alpha = 0.8,
                  pt_clrsp = "inferno",
                  display_title = F)

    })
library(patchwork)
wrap_plots(plist, ncol = 2)
ggplot2::ggsave('neftel_colourful.png')

plotSurfaceComparison(
  object = spata_obj, 
  color_by = neftel,
  smooth = TRUE, 
  #smooth_span = 0.2, 
  pt_alpha = 0.8,
  display_image = TRUE,
  pt_clrsp = "inferno"
  )
ggplot2::ggsave('neftel_normal.png')
