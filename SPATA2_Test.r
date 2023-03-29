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
#install.packages("viridisLite")
#install.packages("msigdbr")
# 
# 
# Declare constants ---------
 library(SPATA2)
 library(magrittr)
 library(tidyverse)
 library(monocle3)
 library(viridisLite)

#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__56_A1" #Need to reverse slashes when copied in Windows
directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/50D" #Need to reverse slashes when copied in Windows
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/16_1" #Need to reverse slashes when copied in Windows
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/45B"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/47"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/A1 - 15A"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/B1(missing currently) - 19"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/C1 - 33B"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/D1_-_35"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__60_B1"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__61_C1"
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__63_D1"

sample_name <- substr(directory_10X, 67, nchar(directory_10X))
#sample_name="50D"
subdir='SPATA2_figures'
sig_trajectory_subdir='signature_trajectories'
custom_fig_subdir='custom_figures'
directory_genesets="C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/all_SPATA2_external.csv"
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
#     directory_10X = directory_10X, # the directory from which to load the data
#     sample_name = sample_name
#   )
# 
# 
# Or alternatively, load premade spata obj --------
spata_obj <- loadSpataObject(spata_obj_dir, verbose = TRUE, update = TRUE)
# Adjust directory instructions for spata object ---------
spata_obj <-
  adjustDirectoryInstructions(
    object = spata_obj,
    to = "spata_object",
    directory_new = spata_obj_dir, # the directory under which to store the object
  )
# Load external gene sets -------
#Load in a .csv where each row is a gene, and each column is a gene set name

#semi-manually add gene sets to the spata_obj
#We only have two gene sets to add for now, but in the future, we may want to package this in a loop
#spata_obj<-addGeneSet(spata_obj,class_name="external",gs_name="SpOT.downreg",genes=dplyr::pull(external_geneset, "SpOT.downreg"))
#spata_obj<-addGeneSet(spata_obj,class_name="external",gs_name="SpOT.upreg",genes=dplyr::pull(external_geneset, "SpOT.upreg"))
#Add neftel
sigs_to_add=colnames(external_geneset)
for (sig in sigs_to_add)
{
  print(sig)
  spata_obj<-addGeneSet(spata_obj,class_name="external",gs_name=sig,genes=dplyr::pull(external_geneset, sig),overwrite=TRUE)
}

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
spata_obj <- createTrajectories(object = spata_obj)

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
###---- 

# Generate misc figures ----
#for colour spectra (pt_clrsp), use "inferno" for regular colour scale; "Inferno" for reverse colour scale
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

color_by="RPL24"
plotSurface(spata_obj,color_by="RPL24",
            smooth = TRUE,
            pt_size = 1.8,
            display_image = TRUE,
            pt_clrsp = "inferno",
            pt_alpha = 0.8) +
  legendTop() +
  ggplot2::labs(title = color_by)
ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))

# Plot spatial trajectories-------
# Color by cluster
plotTrajectory(object = spata_obj,
               trajectory_name = "low_high_HYPOXIA",
               color_by = "seurat_clusters",
               pt_alpha = 0.25, # reduce alpha to highlight the trajectory's course
               display_image = TRUE) +
  legendTop()
ggplot2::ggsave('seurat_low_high_HYPOXIA.png')

#Color by HM_HYPOXIA
plotTrajectory(object = spata_obj,
               trajectory_name = "low_high_HYPOXIA",
               color_by = "HM_HYPOXIA",
               smooth = TRUE,
               pt_alpha = 0.5,
               display_image = TRUE) +
  legendTop()
ggplot2::ggsave('HM_HYPOXIA_low_high_HYPOXIA.png')


hypox_gene_sets=c("BC_P53HYPOXIA_PATHWAY",
                  "BP.GO_REGULATION_OF_CELLULAR_RESPONSE_TO_HYPOXIA",
                  "BP.GO_HYPOXIA_INDUCIBLE_FACTOR_1ALPHA_SIGNALING_PATHWAY",
                  "BP.GO_NEGATIVE_REGULATION_OF_CELLULAR_RESPONSE_TO_HYPOXIA",
                  "BP.GO_NEGATIVE_REGULATION_OF_HYPOXIA_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
                  "BP.GO_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_HYPOXIA",
                  "RCTM_REGULATION_OF_GENE_EXPRESSION_BY_HYPOXIA_INDUCIBLE_FACTOR",
                  "RCTM_CELLULAR_RESPONSE_TO_HYPOXIA",
                  "HM_HYPOXIA")
#Plot trajectory of HM_HYPOXIA
ggp=plotTrajectoryGeneSets(object = spata_obj,
                       trajectory_name = "low_high_HYPOXIA",
                       gene_sets = hypox_gene_sets,
                       smooth_method = "loess",
                       display_facets=TRUE
)
ggp +                                                                # Change font size
  ggplot2::theme(strip.text.x = ggplot2::element_text(size = 3))
ggplot2::ggsave('Hypoxia_gene_sets_trajectory.png')

#Plot cluster breakdown along trajectory
plotTrajectoryFeaturesDiscrete(object = spata_obj,
                               trajectory_name = "low_high_HYPOXIA",
                               discrete_feature = "seurat_clusters",
                               display_trajectory_parts = FALSE)
ggplot2::ggsave('Seurat_clusters_trajectory.png')

#Plot specific genes along trajectory
genes_of_interest <- c("VEGFA", "CA9", "LDHA","MKI67")
plotTrajectoryGenes(object = spata_obj,
                    trajectory_name = "low_high_HYPOXIA",
                    genes = genes_of_interest,
                    smooth_span = 0.2,
                    smooth_se = TRUE,
                    display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                    nrow = 2 # align the sub plots in two rows
)
ggplot2::ggsave('HIF_examples.png')

# Plot specific genes along trajectory----
genes_of_interest <- c("ID3", "ENDOD1")
plotTrajectoryGenes(object = spata_obj,
                    trajectory_name = "low_high_HYPOXIA",
                    genes = genes_of_interest,
                    #smooth_span = 0.2,
                    smooth_se = TRUE,
                    display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                    nrow = 2, # align the sub plots in two rows
) + theme(strip.text.x = element_text(size = 30))
ggplot2::ggsave('EndoD1_ID3_trajectory.png')

# Plot general figures relating to external signatures-----
sigs_to_add=colnames(external_geneset)
for (sig in sigs_to_add)
{
#sig=sigs_to_add[16] #TEMP PLACEHOLDER FOR DEBUGGING 1 GENE SIGNATURE
sig_to_plot=paste('external_',sig,sep='')
#Get a list of genes matched in our sample with an external signature
genes_detected=getGenes(object=spata_obj,of_gene_sets = sig_to_plot)
  if (length(genes_detected)>1){
  #Here, you need to pull a list of genes picked up in our dataset from signature "sig"
  #Determine the length, and have a case or elseif statement handle it
  
  
  #Plot trajectory fit derived from upregulated genes
  plotTrajectoryFit(object = spata_obj,
                    trajectory_name = "low_high_HYPOXIA",
                    variable = sig_to_plot,
                    display_residuals = TRUE) +
    legendTop()
  ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_trajectory_trends.png',sep=''))
  
  #Color trajectory by signature
  plotTrajectory(object = spata_obj,
                 trajectory_name = "low_high_HYPOXIA",
                 color_by = sig_to_plot,
                 smooth = TRUE,
                 pt_alpha = 0.5,
                 display_image = TRUE) +
    legendTop()
  ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_trajectory_line.png',sep=''))
  
  #Plot specific genes constituting signature along trajectory
  genes_of_interest = dplyr::pull(external_geneset, strsplit(sig_to_plot,split='_')[[1]][2])
  plotTrajectoryGenes(object = spata_obj,
                      trajectory_name = "low_high_HYPOXIA",
                      genes = genes_of_interest,
                      smooth_span = 0.2,
                      smooth_se = TRUE,
                      display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                      nrow = 2# align the sub plots in two rows
  )
  ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_examples.png',sep=''))
  
  #Generate one line plot per signature along the trajectory
  plotTrajectoryGeneSets(spata_obj,trajectory_name ='low_high_HYPOXIA' ,
                         gene_sets=sig_to_plot,
  ) +
    ggplot2::labs(title = sig_to_plot)
  ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_trajectory_line_plot.png',sep=''))
  
  #Plotting trajectory heatmap will show >10 genes (signature has around 50)
  hm_colors <- viridis::inferno(n = 100)
  
  xx=plotTrajectoryHeatmap(object = spata_obj,
                           trajectory_name = "low_high_HYPOXIA",
                           variables = genes_of_interest,
                           arrange_rows = "maxima",
                           colors = hm_colors,
                           show_rownames = TRUE,
                           split_columns = FALSE,
                           smooth_span = 0.5,
                           fontsize=6)
  png(filename=paste(sig_trajectory_subdir,'/',sig_to_plot,'_heatmap.png',sep=''),res=res,width = 1024, height = 768)
  
  grid::grid.newpage()
  grid::grid.draw(xx$gtable)
  dev.off()
  } else if (length(genes_detected)==1){
    #Plot trajectory fit derived from upregulated genes
    plotTrajectoryFit(object = spata_obj,
                      trajectory_name = "low_high_HYPOXIA",
                      variable = genes_detected,
                      display_residuals = TRUE) +
      legendTop()
    ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_trajectory_trends.png',sep=''))
    
    #Color trajectory by signature
    plotTrajectory(object = spata_obj,
                   trajectory_name = "low_high_HYPOXIA",
                   color_by = genes_detected,
                   smooth = TRUE,
                   pt_alpha = 0.5,
                   display_image = TRUE) +
      legendTop()
    ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_trajectory_line.png',sep=''))
    
    #Plot specific genes constituting signature along trajectory
    genes_of_interest = dplyr::pull(external_geneset, strsplit(sig_to_plot,split='_')[[1]][2])
    plotTrajectoryGenes(object = spata_obj,
                        trajectory_name = "low_high_HYPOXIA",
                        genes = genes_detected,
                        smooth_span = 0.2,
                        smooth_se = TRUE,
                        display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                        nrow = 2# align the sub plots in two rows
    )
    ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_examples.png',sep=''))
    
    #Generate one line plot per signature along the trajectory
    plotTrajectoryGenes(spata_obj,trajectory_name ='low_high_HYPOXIA' ,
                           genes=genes_detected,
    ) +
      ggplot2::labs(title = sig_to_plot)
    ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_trajectory_line_plot.png',sep=''))
    
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
    png(filename=paste(sig_trajectory_subdir,'/',sig_to_plot,'_heatmap.png',sep=''),res=res,width = 1024, height = 768)

    grid::grid.newpage()
    grid::grid.draw(xx$gtable)
    dev.off()
  }
}

# # Compute spatial pseudotime embeddings ----------
# # compile a cell-data-set (as required by monocle)
# # note, calling this function will prompt you to interactively label
# sample_cds <- transformSpataToCDS(object = spata_obj)
# #Adjust directory instructions for a place to store the cds objects
# spata_obj <- adjustDirectoryInstructions(object = spata_obj,
#                                          to = "cell_data_set",
#                                          directory_new = paste(getwd(),"CDS_obj.RDS",sep="/"))
# #Create a standard spatial plot, with the same clrp as the following monocle plot
# plotSurface(spata_obj, color_by = "seurat_clusters", pt_clrp = "npg", pt_size = 1.4, display_image = TRUE) + 
#   legendBottom()
# ggplot2::labs(title = "Clusters")
# ggplot2::ggsave('Monocle_Clusters.png')
# 
# #Create a monocle plot showing some kind of trajectory
# plot_cells(sample_cds,
#              reduction_method = "UMAP",
#              color_cells_by = "seurat_clusters" # visualize a transferred variable
#   ) +
#   scale_color_add_on(aes = "color", variable = "discrete", clrp = "npg") +
#   ggplot2::labs(title = "Transcriptomic Velocities")
#   ggplot2::ggsave('Monocle_Transcriptomic_Velocities.png')
# 
# 
# 
# 
# 
# 
#   
# Replicate Neftel cellular states plots --------
#Manually state what the 4 cell states should be, from the external signature list
four_states <- c("external_Neftel.OPC.like", "external_Neftel.NPC.like", "external_Neftel.AC.like", "external_Neftel.MES.like")
plotSurfaceComparison(object = spata_obj, 
                      color_by = four_states, 
                      smooth = TRUE, 
                      smooth_span = 0.25,
                      display_image = TRUE,
                      pt_size = 1.5) + 
  legendNone()
ggplot2::ggsave('FourCellularStates_surface.png')
#Replicate PCA 4 quadrant plot
plotFourStates(object = spata_obj, 
               states = four_states,
               method_gs = "mean",
               # using color here to underline that the mathematical algorithm is valid
               color_by = "external_Neftel.NPC.like", 
               pt_alpha = 0.75,
               pt_size = 1) +
  legendNone()
ggplot2::ggsave('FourCellularStates_PCPlot_MES.png')
#Show that MES is correlated with Hypoxia
plotScatterplot(object = spata_obj, 
                variables = c("external_Neftel.MES.like", "external_Ron.upreg"),
                smooth = TRUE, 
                smooth_method = "lm",
                pt_clr = "black", 
                pt_size = 1.8,
                pt_alpha = 0.75)
ggplot2::ggsave('FourCellularStates_MESvsRonHypoxia.png')



# Save SpataObj----
saveSpataObject(object = spata_obj)#Save SpataObj
#saveCorrespondingCDS(cds = sample_cds, object = spata_obj) 

# Compare Surface Plots (for reverse colour plots e.g. down DAEGs)----
#for colour spectra (pt_clrsp), use "inferno" for regular colour scale; "Inferno" for reverse colour scale

# compare SINGLE GENE expression on the surface
# store example genes of interest as character vectors
  genes_b <- c("ID3", "ENDOD1")
  
  plotSurfaceComparison(object = spata_obj, 
                        color_by = genes_b,
                        smooth = TRUE, 
                        #smooth_span = 0.2, 
                        pt_alpha = 0.8,
                        display_image = TRUE,
                        pt_clrsp = "inferno") +
    ggplot2::labs(title = "ID3 and EndoD1 Spatial Maps")
  ggplot2::ggsave(paste(custom_fig_subdir,'/ID3_ENDOD1_colourNormalspatial_plot.png',sep=''))
  
  plotSurfaceComparison(object = spata_obj, 
                        color_by = genes_b,
                        smooth = TRUE, 
                        #smooth_span = 0.2, 
                        pt_alpha = 0.8,
                        display_image = TRUE,
                        pt_clrsp = "Inferno") +
    ggplot2::labs(title = "ID3 and EndoD1 Spatial Maps")
  ggplot2::ggsave(paste(custom_fig_subdir,'/ID3_ENDOD1_colourRevspatial_plot.png',sep=''))
  
# compare GENESET expression
  sigs_b <- c("external_SpOT.upreg","external_SpOT.downreg","external_Ron.downreg","external_Ron.upreg")
  
  plotSurfaceComparison(object = spata_obj, 
                        color_by = sigs_b,
                        smooth = TRUE,
                        pt_alpha = 0.5,
                        display_image = TRUE,
                        pt_clrsp = "inferno") +
    ggplot2::labs(title = "PIMO and DAEG Geneset Spatial Maps")
  ggplot2::ggsave(paste(custom_fig_subdir,'/DAEGs_PIMO_colourNormalspatial_plot.png',sep=''))
  
  plotSurfaceComparison(object = spata_obj, 
                        color_by = sigs_b,
                        smooth = TRUE,
                        pt_alpha = 0.5,
                        display_image = TRUE,
                        pt_clrsp = "Inferno") +
    ggplot2::labs(title = "PIMO and DAEG Geneset Spatial Maps")
  ggplot2::ggsave(paste(custom_fig_subdir,'/DAEGs_PIMO_colourRevspatial_plot.png',sep=''))
  

# open application to obtain a list of plots
#plots <- plotSurfaceInteractive(object = spata_obj)

  
# Transcript Counts Surface Plot----
  plots <- plotSurfaceInteractive(object = spata_obj)
  plots$transcriptCounts + 
       ggplot2::labs(title = "Transcript Counts")
  ggplot2::ggsave('transcript_counts.png')
  

  

# ---- Next two blocks are duplicates (pretty much) of "general figures relating----
# to external signatures". Can be used if you want to generate specific 
# trajectory-related plots, but trajectory-related plots for all external
# signatures are already generated in "general figures relating to external
# signatures" above.

  
  
#Plot trajectory trends identified on our signatures-----
  #Plot figures relating to external_Ron.upreg signature-----
  # sig_to_plot="external_Ron.upreg"
  # 
  # #Plot trajectory fit derived from upregulated genes
  # plotTrajectoryFit(object = spata_obj,
  #                   trajectory_name = "low_high_HYPOXIA",
  #                   variable = sig_to_plot,
  #                   pt_clrsp = "inferno",
  #                   display_residuals = TRUE) +
  #   legendTop()
  # ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_trajectory_trends.png',sep=''))
  # 
  # #Color trajectory by signature
  # plotTrajectory(object = spata_obj,
  #                trajectory_name = "low_high_HYPOXIA",
  #                color_by = sig_to_plot,
  #                pt_clrsp = "inferno",
  #                smooth = TRUE,
  #                pt_alpha = 0.5,
  #                display_image = TRUE) +
  #   legendTop()
  # ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_trajectory_line.png',sep=''))
  # 
  # #Plot specific genes constituting signature along trajectory
  # genes_of_interest = dplyr::pull(external_geneset, strsplit(sig_to_plot,split='_')[[1]][2])
  # plotTrajectoryGenes(object = spata_obj,
  #                     trajectory_name = "low_high_HYPOXIA",
  #                     genes = genes_of_interest,
  #                     smooth_span = 0.2,
  #                     smooth_se = TRUE,
  #                     display_facets = TRUE, # use facet_wrap() to split the plot in four parts
  #                     nrow = 2# align the sub plots in two rows
  # )
  # ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_examples.png',sep=''))
  # hm_colors <- viridis::inferno(n = 100)
  # #Plotting trajectory heatmap will show >10 genes (signature has around 50)
  # xx=plotTrajectoryHeatmap(object = spata_obj,
  #                          trajectory_name = "low_high_HYPOXIA",
  #                          variables = genes_of_interest,
  #                          arrange_rows = "maxima",
  #                          colors = hm_colors,
  #                          show_rownames = TRUE,
  #                          split_columns = FALSE,
  #                          smooth_span = 0.5)
  # png(filename=paste(sig_trajectory_subdir,'/',sig_to_plot,'_heatmap.png',sep=''),res=res,width = 1024, height = 768)
  # 
  # grid::grid.newpage()
  # grid::grid.draw(xx$gtable)
  # dev.off()
  # 
  #Plot figures relating to external_Ron.downreg signature-----
  sig_to_plot="external_Ron.downreg"
  
  #Color trajectory by signature
  plotTrajectory(object = spata_obj,
                 trajectory_name = "low_high_HYPOXIA",
                 color_by = sig_to_plot,
                 smooth = TRUE,
                 pt_alpha = 0.5,
                 pt_clrsp = "Inferno",
                 display_image = TRUE) +
    legendTop()
  ggplot2::ggsave(paste(sig_trajectory_subdir,'/',sig_to_plot,'_colourRev_trajectory_line.png',sep=''))
  
  hm_colors <- viridis::inferno(n = 100,direction=-1)
  
  #Plotting trajectory heatmap will show >10 genes (signature has around 50)
  #SpOT_allgenes=external_geneset$SpOT.ALLGENES
  Ron_down=external_geneset$Ron.downreg
  xx=plotTrajectoryHeatmap(object = spata_obj,
                           trajectory_name = "low_high_HYPOXIA",
                           variables = Ron_down,
                           arrange_rows = "minima",
                           colors = hm_colors,
                           show_rownames = TRUE,
                           split_columns = FALSE,
                           smooth_span = 0.5
  )
  png(filename=paste(sig_trajectory_subdir,'/',sig_to_plot,'_colourRev_heatmap.png',sep=''),res=res,width = 1024, height = 768)
  
  grid::grid.newpage()
  grid::grid.draw(xx$gtable)
  dev.off()
  

##Generate user defined trajectory lines figures ----
# #color_by='SPARCL1'
# color_by='external_SpOT.upreg'
# #color_by='external_Neftel.OPC.like'
# #color_by='external_Ron.downreg'
# plotTrajectoryGeneSets(spata_obj,trajectory_name ='low_high_HYPOXIA',
#                        gene_sets=color_by
# ) +
#   ggplot2::labs(title = color_by)
# ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_trajectory_line_plot.png',sep=''))
# #Plot surface and highlight the spots that were used in making the trajectory gene plot
# plotTrajectory(spata_obj,trajectory_name ='low_high_HYPOXIA' ,
#                color_by=color_by,
#                pt_alpha = 0.25, # reduce alpha to highlight the trajectory's course
#                display_image = TRUE
# ) +
#   ggplot2::labs(title = color_by)
# ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_trajectory_surface_plot.png',sep=''))

###---- 
  
# load premade spata_obj2-----
spata_obj2 <- loadSpataObject(spata_obj_dir2, verbose = TRUE, update = TRUE)
#Differential Expression Analysis----
#clustering
pythonObsClusters=paste(directory_10X,"/outs/obsClusters.csv",sep='')
leiden_genedf=read.csv(pythonObsClusters)
#Remove empty columns from external_geneset
empty_columns <- colSums(is.na(leiden_genedf) | leiden_genedf == "") == nrow(leiden_genedf)
leiden_genedf=leiden_genedf[, !empty_columns]

# Rename barcodes column
names(leiden_genedf)[names(leiden_genedf) == "X"] <- "barcodes"
leiden_genedf$leiden_gene <- as.character(leiden_genedf$leiden_gene)

#subset the spata object to include only spots included after scanpy preprocessing
scanpy_barcodes <- leiden_genedf$barcodes
spata_obj2 <- subsetByBarcodes_CountMtr(
  spata_obj,
  scanpy_barcodes
)

spata_obj2 <- 
  addFeatures(object = spata_obj2, 
              feature_names = c("leiden_gene"), 
              feature_df = leiden_genedf,
              overwrite = TRUE,
              key = "barcodes")

# spata_obj <-
#   addFeatures(object = spata_obj,
#               feature_names = c("leiden_gene"),
#               feature_df = leiden_genedf,
#               overwrite = TRUE,
#               key = "barcodes")

# monocle_clusters <- findMonocleClusters(object = spata_obj, 
#                                       preprocess_method = "PCA", 
#                                       reduction_method = "UMAP", 
#                                       cluster_method = "leiden", 
#                                       k = 5, 
#                                       num_iter = 5)
# # output
# examineClusterResults(monocle_clusters)
# 
# spata_obj <-
#   addFeatures(object = spata_obj,
#               feature_names = c("leiden_gene"),
#               feature_df = leiden_genedf,
#               overwrite = TRUE,
#               key = "barcodes")

plotSurface(object = spata_obj2, 
          color_by = "leiden_gene",
          display_image = TRUE,
          #pt_clrp = "jama", 
          pt_size = 1.8) +
  ggplot2::labs(title = "Scanpy leiden clusters", color = "Scanpy Leiden")
ggplot2::ggsave('scanpy_leiden.png')

spata_obj2 <- 
  runDeAnalysis(object = spata_obj2,
                across = "leiden_gene", # across which identity groups
                method_de = c("wilcox") # with which methods
  )

# get an overview about the de-analysis results stored in your object
#printDeaOverview(spata_obj2)

DeaDf <- getDeaResultsDf(object = spata_obj2, 
                across = "leiden_gene", 
                method_de = "wilcox",
                max_adj_pval = 0.025)
write.csv(DeaDf, paste(directory_10X,"\\",subdir,"\\DeaDf.csv",sep=""),row.names = FALSE)

#Gene Set Enrichment Analysis-----
leiden_gene_colors = c("0"='#1f77b4',
                       "1"='#ff7f0e',
                       "2"='#279e68',
                       "3"='#d62728',
                       "4"='#aa40fc',
                       "5"='#8c564b',
                       "6"= '#e377c2',
                       "7"='#b5bd61',
                       "8"= '#17becf',
                       "9"= '#aec7e8',
                       "10"='#ffbb78',
                       "11"='#98df8a'
)
keyword <- c("HM_")
keywords_cellstates <- c("Neftel","Verhaak","Pugh","external_Ravi","external_SpOT.upreg","external_SpOT.downreg","external_Ron.upreg","external_Ron.downreg")
geneset_dataframe <- getGeneSetDf(spata_obj2)

# geneset_list <- geneset_dataframe %>%
#   filter(str_detect(ont, "HM_")) %>%
#   select(ont)
# geneset_names <- unique(geneset_list$ont)

geneset_names <- geneset_dataframe %>%
  filter(str_detect(ont, paste(keywords_cellstates, collapse = "|"))) %>%
  pull(ont)
geneset_names <- unique(geneset_names)

#run GSEA
spata_obj2 <- 
  runGsea(
    object = spata_obj2,
    gene_set_names = geneset_names,
    across = "leiden_gene",
    methods_de = "wilcox",
    reduce = FALSE
  )

#get GSEA results
gseaDf <- getGseaDf(
  object = spata_obj2, 
  across = "leiden_gene",
  method_de = "wilcox"
)
write.csv(gseaDf, paste(directory_10X,"\\",subdir,"\\cellstates_gseaDf.csv",sep=""),row.names = FALSE)

#GSEA dot plot
p3 <- plotGseaDotPlot(
  object = spata_obj2,
  across = "leiden_gene",
  #across_subset = c("1", "2", "4", "5", "6"),
  color_by = "fdr",
  #n_gsets = 5,
  pt_alpha = 0.8,
  #transform_with = list("fdr" = c("log10")),
  by_group = FALSE, # merge in one plot
  size_by = "fdr",
  pt_clrsp = "Inferno" ,
)
p3
ggplot2::ggsave('GSEAdotplot_cellstates.png')

four_states <- c("external_Neftel.MES1", "external_Neftel.MES2", "external_Neftel.NPC1", "external_Neftel.NPC2","external_Neftel.AC.like","external_Neftel.OPC.like")
p1 <- plotSurfaceComparison(object = spata_obj, 
                      color_by = four_states, 
                      smooth = TRUE, 
                      smooth_span = 0.25,
                      pt_size = 1.9,
                      display_image = TRUE,
                      pt_alpha = 1
                      #method_gs =  'ssgsea',
                      ) + 
  legendNone()
p2 <- 
  plotSurface(object = spata_obj2,
              color_by = "leiden_gene",
              pt_size = 2.2,
              clrp_adjust = leiden_gene_colors,
              display_image = TRUE) +
  ggplot2::labs(color = "Clusters")

# combine with patchwork 
p1 + legendTop() +
  p2 + legendTop() 

p2 + legendRight() +
  p3 + legendRight() 
ggplot2::ggsave('cluster_GSEAdotplot_cellstates.png')
#Adjust directory instructions for spata object2 ---------
spata_obj2 <-
  adjustDirectoryInstructions(
    object = spata_obj2,
    to = "spata_object",
    directory_new = spata_obj_dir2, # the directory under which to store the object
  )

saveSpataObject(object = spata_obj2)#Save SpataObj

#Trying GSEA again with clusterprfiler----
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(msigdbr)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

leidenclusters <- unique(DeaDf$leiden_gene)

# this is old and inefficient but it does work

# df_master=NULL
# 
# for(cluster in leidenclusters){
#   dea_df <-
#       getDeaResultsDf(
#         object = spata_obj2,
#         across = "leiden_gene",
#         across_subset = cluster,
#         method_de = "wilcox",
#         max_adj_pval = 0.025
#       )
#   # sort (required for clusterProfiler)
#   dea_df <-dea_df[order(-dea_df$avg_logFC),]
#   # We will lose some genes here because not all IDs will be converted
#   ids<-bitr(dea_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
#   # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
#   dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
#   dea_df <- dea_df[dea_df$gene %in% dedup_ids$SYMBOL,]
#   # Create a new column in df2 with the corresponding ENTREZ IDs
#   dea_df$Entrez_ID = dedup_ids$ENTREZID
#   #append to master df
#   df_master = rbind(df_master, dea_df)
#   # Create a vector of the gene unuiverse
#   # kegg_gene_list <- df2$avg_logFC
#   # Name vector with ENTREZ ids
#   # names(kegg_gene_list) <- df2$Entrez_ID
# }

# this does the exact same thing
ids<-bitr(DeaDf$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
#dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
DeaDf <- DeaDf[DeaDf$gene %in% ids$SYMBOL,]
map = setNames(c(ids$ENTREZID),c(ids$SYMBOL))
DeaDf$Entrez_ID = map[DeaDf$gene]

geneList_master <- vector(mode = "list")
for(cluster in leidenclusters){
  geneList <- DeaDf %>%
    filter(leiden_gene==cluster) %>% 
    pull(Entrez_ID)
  geneList_master[[cluster]] <- geneList
}

m_t2g1 <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
                   dplyr::select(gs_name, entrez_gene)
mt2g2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
mt2g = rbind(m_t2g1, mt2g2)
mt2g3 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
cell_marker_data <- read.csv('C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/brain_cell_markers.csv')
cells <- cell_marker_data %>%
  dplyr::select(Cell.name, Cell.marker) %>%
  dplyr::mutate(Cell.marker = strsplit(Cell.marker, ', ')) %>%
  tidyr::unnest(cols = c(Cell.marker))
ids <- bitr(cells$Cell.marker, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
cells <- cells[cells$Cell.marker %in% dedup_ids$SYMBOL,]
cells = cells[!duplicated(cells[c("Cell.marker")]),]
cells$Entrez_ID = dedup_ids$ENTREZID
t2g <- cells %>%
  dplyr::select(Cell.name, Entrez_ID)%>%
  dplyr::mutate(Entrez_ID = strsplit(Entrez_ID, ', ')) %>%
  tidyr::unnest(cols = c(Entrez_ID))

ck <- compareCluster(geneCluster = geneList_master,
                     fun = enricher,
                     TERM2GENE=mt2g,
                     minGSSize = 10,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
ck <- setReadable(ck, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")

write.csv(ck@compareClusterResult, paste(directory_10X,"\\",subdir,"\\C5C2_enricher.csv",sep=""),row.names = FALSE)
dotplot(ck, font.size = 10,showCategory = 5,label_format = 80)
ggplot2::ggsave('C5C2_dotplot_enricher.png')

# get proper data in numeric matrix form for import to heatmap plots
gseaDf_pivot <- gseaDf %>%
  dplyr::select(leiden_gene, label, pval)
gseaDf_pivot <- gseaDf_pivot[, c("label","leiden_gene","pval")]
#gseaDf_pivot$log10pval <- apply(gseaDf_pivot[,"pval"], 2, function(x) { -1*log10(x) } )
gseaDf_pivot <- transform(gseaDf_pivot, log10pval = (-1*log10(pval)))
gseaDf_pivot <- gseaDf_pivot %>% dplyr::select(-pval)
gseaDf_pivot <- gseaDf_pivot %>% pivot_wider(
  names_from = leiden_gene,
  values_from = log10pval)
# gseaDf_pivot <- gseaDf_pivot %>% 
#   mutate(label = factor(label, label))
# mat <- gseaDf_pivot
row_names <- dplyr::pull(gseaDf_pivot, label)
#mat <- mat %>% dplyr::select(-label)
#mat <- as.matrix(mat)
gseaDf_pivot <- gseaDf_pivot %>% dplyr::select(-label)
rownames(gseaDf_pivot) <- c(row_names)
mat <- as.matrix(gseaDf_pivot)

# add -log10(pval) as a column in gseaDf
newgseaDf <- transform(gseaDf, log10pval = (-1*log10(pval)))
row_names <- dplyr::pull(newgseaDf, label)

c2c5 <- data.frame(ck@compareClusterResult)

c2c5_forhm=NULL

for(cluster in leidenclusters){
  clusterdf <- c2c5_df %>%
    filter(Cluster==cluster)
  clusterdf <- clusterdf[1:5, ]
  
  c2c5_forhm = rbind(c2c5_forhm,clusterdf)
}

pb<-ggplot(newgseaDf, aes(leiden_gene, label, fill = log10pval)) + 
  geom_tile(colour="gray20", linewidth=1.5, stat="identity") + 
  scale_fill_viridis_c(option="B") +
  #scale_y_continuous(breaks=1:23, labels=waiver())+
  xlab("") + 
  ylab("") +
  ggtitle("Enrichment - Transcriptional Cell State Signatures") +
  theme(
    plot.title = element_text(color="white",hjust=0,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="gray20"),
    panel.background = element_rect(fill="gray20"),
    panel.border = element_rect(fill=NA,color="gray20", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(), 
    axis.text = element_text(color="white", size=rel(1.5)),
    axis.text.y  = element_text(hjust=1),
    legend.text = element_text(color="white", size=rel(1.3)),
    legend.background = element_rect(fill="gray20"),
    legend.position = "bottom",
    legend.title=element_blank()
  )

c2c5_df <- transform(gseaDf, log10pval = (-1*log10(pvalue)))
row_names <- dplyr::pull(c2c5_df, label)

ggplot(newgseaDf, aes(leiden_gene, label, fill= log10pval)) + 
  geom_tile(colour="gray20", linewidth=1.5, stat="identity") +
  scale_fill_viridis_c(option="B",direction = 1,
                       na.value = "grey50",guide = "colourbar") +
  ggtitle("Enrichment - Transcriptional Cell State Signatures")

#add -log10(pval) 
c2c5_df <- transform(c2c5_df, log10pval = (-1*log10(pvalue)))
#to prep for plotting, take the top 10 genesets for each cluster
c2c5_forhm=NULL
for(cluster in leidenclusters){
  clusterdf <- c2c5_df %>%
    filter(Cluster==cluster)
  clusterdf <- clusterdf[1:10, ]
  
  c2c5_forhm = rbind(c2c5_forhm,clusterdf)
}
row_names <- unique(dplyr::pull(c2c5_forhm, ID))
c2c5_forhm <- c2c5_df %>%
  filter(ID %in% row_names)

#cut down the genesets names for plotting
c2c5_forhm <- as_tibble(c2c5_forhm)
nameTooLong <- function(ID) {
  case_when(
    nchar(ID) > 40 ~ as.character(gsub("_"," ", str_trunc(ID, 40, "right"))),
    .default = as.character(gsub("_"," ", ID))
  )
}
c2c5_forhm <- c2c5_forhm %>%
  mutate(ID = nameTooLong(ID))

hmc2c5 <- ggplot(c2c5_forhm, aes(Cluster, ID, fill= log10pval)) + 
  geom_tile(colour="gray20", linewidth=0.5, stat="identity") +
  scale_fill_viridis_c(option="B",direction = 1,
                       na.value = "white",guide = "colourbar") +
  ggtitle("Enrichment - C2 and C5 Genesets")
hmc2c5
#ggplot2::ggsave('C5C2_heatmap_enricher.png')

# heatmap(mat,
#         #dendrogram = "row",
#         xlab = "leiden_gene",
#         ylab = "label",
#         #col =  viridis::inferno(),
#         labRow = row_names,
#         keep.dendro = FALSE,
#         Rowv = NA,
#         Colv = NA,
#         #margins = c(60,100,40,20),
#         #grid_color = "white",
#         #grid_width = 0.00001,
#         #file = "heatmaply_gseacellstates.png",
#         #plot_method = "ggplot",
#         #titleX = FALSE,
#         #hide_colorbar = FALSE,
#         #branches_lwd = 0.1,
#         #label_names = c("Signature", "Cluster", "-log10(pval)"),
#         #fontsize_row = 5, fontsize_col = 5,
#         #labCol = colnames(mat),
#         #labRow = rownames(mat),
#         #heatmap_layers = theme(axis.line=element_blank())
# ) + scale_color_viridis(option="magma")

# heatmaply interactive heatmap. Can't save without installing a driver so kind of impractical. But looks cool! 
heatmaply(gseaDf_pivot,
          #dendrogram = "row",
          xlab = "leiden_gene",
          ylab = "label", 
          dendrogram = "none",
          labRow = row_names,
          #margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.00001,
          file = "heatmaply_gseacellstates.png",
          plot_method = "ggplot",
          #titleX = FALSE,
          #hide_colorbar = FALSE,
          #branches_lwd = 0.1,
          label_names = c("Signature", "Cluster", "-log10(pval)"),
          fontsize_row = 5, fontsize_col = 5,
          #labCol = colnames(mat),
          #labRow = rownames(mat),
          #heatmap_layers = theme(axis.line=element_blank())
)

r <- cor(mtcars)
## We use this function to calculate a matrix of p-values from correlation tests
## https://stackoverflow.com/a/13112337/4747043
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
p <- cor.test.p(mtcars)

heatmaply_cor(
  r,
  node_type = "scatter",
  point_size_mat = -log10(p), 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation")
)

#replicating ravi supplemental fig 4
# first they did Go enrichment and pathway enrichment in cluster profiler for each cluster
genelist_allgenes <-bitr(getGenes(spata_obj2), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
ego <- compareCluster(geneCluster = geneList_master,
                      fun = enrichGO,
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",
                      keyType = "ENTREZID",
                      readable = TRUE)
write.csv(ego, paste(subdir,"\\comparecluster_ego.csv"))
dotplot(ego, font.size = 10,showCategory = 4,label_format = 100)
ggplot2::ggsave('dotplot_enrichGO.png')

leidenclusters <- as.character(leidenclusters)
for(cluster in leidenclusters){
  ego <- enrichGO(gene = geneList_master[[cluster]],
                  #universe = genelist_allgenes$ENTREZID,
                  OrgDb = "org.Hs.eg.db",
                  ont = "ALL",
                  keyType = "ENTREZID",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable = TRUE)
  if(nrow(ego@result)>5){
    # bar plot
    mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory=20)
    #ggplot2::ggsave('cnetplot_6.png')
    ggplot2::ggsave(paste(subdir,'/','barplot_',cluster,'.png',sep=''))
  }
}


