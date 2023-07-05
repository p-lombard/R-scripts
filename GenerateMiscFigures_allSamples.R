# Declare constants ---------
library(SPATA2)
library(magrittr)
library(tidyverse)
library(monocle3)
library(viridisLite)

directories = c("C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__56_A1",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/50D",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/16_1", 
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/45B",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/47",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/A1 - 15A",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/B1(missing currently) - 19",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/C1 - 33B",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/D1_-_35",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__60_B1",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__61_C1",
                "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__63_D1"
)

directory_genesets="C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/GBmap_lv3.csv"
external_geneset=read.csv(directory_genesets)
#Remove empty columns from external_geneset
empty_columns <- colSums(is.na(external_geneset) | external_geneset == "") == nrow(external_geneset)
external_geneset=external_geneset[, !empty_columns]

for(dir in directories){
  directory_10X <- dir
  sample_name <- substr(directory_10X, 67, nchar(directory_10X))
  #sample_name="50D"
  subdir='SPATA2_figures'
  sig_trajectory_subdir='signature_trajectories'
  custom_fig_subdir='custom_figures'
  directory_genesets="C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/GBmap_lv3.csv"
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
  spata_obj <- loadSpataObject(spata_obj_dir2, verbose = TRUE, update = TRUE)
  # Adjust directory instructions for spata object ---------
  spata_obj <-
    adjustDirectoryInstructions(
      object = spata_obj,
      to = "spata_object",
      directory_new = spata_obj_dir2, # the directory under which to store the object
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

  # Generate misc figures ----
  #for colour spectra (pt_clrsp), use "inferno" for regular colour scale; "Inferno" for reverse colour scale
  # color_by='SOX4'
  # # color_by='external_SpOT.upreg'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))
  # 
  # color_by='external_SpOT.upreg'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "Inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_colourRevspatial_plot.png',sep=''))
  # 
  # color_by='external_Ron.downreg'
  # 
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "Inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'colourRev_spatial_plot.png',sep=''))
  # 
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))
  #   
  # color_by='external_Ron.upreg'
  # 
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))
  # 
  # color_by='external_Neftel.MES.like'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))
  # 
  # color_by='external_Neftel.MES.like'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "Inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_colourRevSpatial_plot.png',sep=''))
  # 
  # color_by='external_Neftel.MES2'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))
  # 
  # color_by='external_Neftel.MES1'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))
  # 
  # color_by='external_Neftel.AC.like'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))
  # 
  # color_by='external_Neftel.OPC.like'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))
  # 
  # color_by='external_Neftel.NPC1'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))
  # 
  # color_by='external_Neftel.NPC2'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))
  # 
  # color_by='external_SpOT.downreg'
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_Spatial_plot.png',sep=''))
  # 
  # color_by = "HM_HYPOXIA"
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))
  # 
  # plotSurface(spata_obj,color_by=color_by,
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "Inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_colourRevSpatial_plot.png',sep=''))
  # 
  # color_by="RPL24"
  # plotSurface(spata_obj,color_by="RPL24",
  #             smooth = TRUE,
  #             pt_size = 1.8,
  #             display_image = TRUE,
  #             pt_clrsp = "inferno",
  #             pt_alpha = 0.8) +
  #   legendTop() +
  #   ggplot2::labs(title = color_by)
  # ggplot2::ggsave(paste(custom_fig_subdir,'/',color_by,'_spatial_plot.png',sep=''))
  # 
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
  
  
}