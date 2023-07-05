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
# devtools::install_github(repo = "theMILOlab/SPATAData")
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
#install.packages("gridExtra")
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# 
# Set Constants -----
library(SPATA2)
library(magrittr)
library(tidyverse)
library(monocle3)
library(plotly)

# directories = c("C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__56_A1",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/50D",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/16_1",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/45B",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/47",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/A1 - 15A",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/B1(missing currently) - 19",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/C1 - 33B",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/D1_-_35",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__60_B1",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__61_C1",
#   "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__63_D1"
# )

directories = c("C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__56_A1",
  "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/16_1",
  "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/45B",
  "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/A1 - 15A",
  "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/B1(missing currently) - 19",
  "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/C1 - 33B",
  "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/D1_-_35",
  "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__60_B1",
  "C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__61_C1"
)


directory_genesets="C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/GBmap_lv3.csv"
external_geneset=read.csv(directory_genesets)
#Remove empty columns from external_geneset
empty_columns <- colSums(is.na(external_geneset) | external_geneset == "") == nrow(external_geneset)
external_geneset=external_geneset[, !empty_columns]

maps.names <- substr(directories, 67, nchar(directories))
maps <- vector("list", length(maps.names))
names(maps) <- maps.names

# seurat_objs.names <- c()
# for (dir in directories){
#   name = substr(dir, 67, nchar(dir))
#   seurat_objs.names <- append(seurat_objs.names,name)
# }
# seurat_objs = vector(mode='list', length=length(seurat_objs.names))
# names(seurat_objs) <- seurat_objs.names

# Main----
for(dir in directories){
  directory_10X <- dir
  sample_name <- substr(directory_10X, 67, nchar(directory_10X))
  #sample_name="50D"
  subdir='PrimIDHwt_superclusters'
  #sig_trajectory_subdir='signature_trajectories'
  #custom_fig_subdir='custom_figures'
  directory_genesets="C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/GBmap_lv3.csv"
  external_geneset=read.csv(directory_genesets)
  #Remove empty columns from external_geneset
  empty_columns <- colSums(is.na(external_geneset) | external_geneset == "") == nrow(external_geneset)
  external_geneset=external_geneset[, !empty_columns]
  
  res=100 #Resolution for figures
  
  supercluster_colours = c('#d62728',"#61D04F","#2297E6","#28E2E5","#CD0BBC","#F5C710")
  #supercluster_colours[supercluster_colours == "black"] <- "#E69F00"
  names(supercluster_colours) <- c("A","B","C","D","E","F") 
  
  # Create output directory ------------
  dir.create(file.path(directory_10X, subdir), showWarnings = FALSE) #create base figure directory
  setwd(file.path(directory_10X, subdir))
  #dir.create(sig_trajectory_subdir, showWarnings = FALSE) #create subdirectory for iterative sig plotting
  #dir.create(custom_fig_subdir, showWarnings = FALSE) #create subdirectory for iterative sig plotting
  
  # Create directory to read/write spata objects
  spata_obj_dir=paste(getwd(),"spata_obj.RDS",sep="/")
  spata_obj_dir2=paste(getwd(),"spata_obj2.RDS",sep="/")
  spata_obj_dir3=paste(getwd(),"spata_obj3.RDS",sep="/")
  
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
  # # Or alternatively, load premade spata obj --------
  # spata_obj <- loadSpataObject(paste(directory_10X,"SPATA2_figures/spata_obj.RDS",sep="/"), verbose = TRUE, update = TRUE)
  # # Adjust directory instructions for spata object ---------
  # spata_obj <-
  #   adjustDirectoryInstructions(
  #     object = spata_obj,
  #     to = "spata_object",
  #     directory_new = spata_obj_dir, # the directory under which to store the object
  #   )
  # # Load external gene sets -------
  # #Load in a .csv where each row is a gene, and each column is a gene set name
  # 
  # #semi-manually add gene sets to the spata_obj
  # #We only have two gene sets to add for now, but in the future, we may want to package this in a loop
  # #spata_obj<-addGeneSet(spata_obj,class_name="external",gs_name="SpOT.downreg",genes=dplyr::pull(external_geneset, "SpOT.downreg"))
  # #spata_obj<-addGeneSet(spata_obj,class_name="external",gs_name="SpOT.upreg",genes=dplyr::pull(external_geneset, "SpOT.upreg"))
  # #Add neftel
  # sigs_to_add=colnames(external_geneset)
  # for (sig in sigs_to_add)
  # {
  #   print(sig)
  #   spata_obj<-addGeneSet(spata_obj,class_name="external",gs_name=sig,genes=dplyr::pull(external_geneset, sig),overwrite=TRUE)
  # }
  # load premade spata_obj2-----
  spata_obj2 <- loadSpataObject(paste(directory_10X,"SPATA2_figures/spata_obj2.RDS",sep="/"), verbose = TRUE, update = TRUE)
  # # create spata_obj2----
  # #clustering
  # # pythonObsClusters=paste(directory_10X,"/outs/obsClusters.csv",sep='')
  # pythonObsClusters=paste(directory_10X,"/outs/obsClusters.csv",sep='')
  # leiden_genedf=read.csv(pythonObsClusters)
  # #Remove empty columns from external_geneset
  # empty_columns <- colSums(is.na(leiden_genedf) | leiden_genedf == "") == nrow(leiden_genedf)
  # leiden_genedf=leiden_genedf[, !empty_columns]
  # 
  # # Rename barcodes column
  # names(leiden_genedf)[names(leiden_genedf) == "X"] <- "barcodes"
  # leiden_genedf$leiden_gene <- as.character(leiden_genedf$leiden_gene)
  # 
  # #subset the spata object to include only spots included after scanpy preprocessing
  # scanpy_barcodes <- leiden_genedf$barcodes
  # spata_obj2 <- subsetByBarcodes_CountMtr(
  #   spata_obj,
  #   scanpy_barcodes
  # )
  # 
  # spata_obj2 <-
  #   addFeatures(object = spata_obj2,
  #               feature_names = c("leiden_gene"),
  #               feature_df = leiden_genedf,
  #               overwrite = TRUE,
  #               key = "barcodes")
  # 
  # plotSurface(object = spata_obj2,
  #             color_by = "leiden_gene",
  #             display_image = TRUE,
  #             #pt_clrp = "jama",
  #             pt_size = 1.8) +
  #   ggplot2::labs(title = "Scanpy leiden clusters (low)", color = "Scanpy Leiden")
  # ggplot2::ggsave('scanpy_leidenLow.png')
  # Get DEA Df----
  DeaDf <- getDeaResultsDf(object = spata_obj2,
                           across = "leiden_gene",
                           method_de = "wilcox",
                           max_adj_pval = 0.025)
  write.csv(DeaDf, paste(directory_10X,"\\",subdir,"\\DeaDf.csv",sep=""),row.names = FALSE)
  
  # plot leiden clusters-----
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
  p2 <-
    plotSurface(object = spata_obj2,
                color_by = "leiden_gene",
                pt_size = 1.8,
                clrp_adjust = leiden_gene_colors,
                display_image = TRUE) +
    ggplot2::labs(title = "Non-integrated Gene-Based Clustering",color = "Leiden Clusters")
  
  # Save spata_obj2 ---------
  spata_obj2 <-
    adjustDirectoryInstructions(
      object = spata_obj2,
      to = "spata_object",
      directory_new = spata_obj_dir2, # the directory under which to store the object
    )

  saveSpataObject(object = spata_obj2)#Save SpataObj

  # # create spata_obj3----
  # # load cluster supercluster assignments
  # ppath = "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\primaryIDHwt\\c2c5 gsea results"
  # # ppath = "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\whole_GBM_dataset_analysis\\primary wt gsea results"
  # load(paste0(ppath, "\\clusterAssignments_NES058_FDR001.RData"))
  # 
  # # turn this into a dataframe shape we can add to spata_obj2 (2 columns,
  # # one representing barcodes, one representing cluster label)
  # leiden_barcodes_df = getFeatureDf(spata_obj2)
  # leiden_barcodes_df <- leiden_barcodes_df[,c('barcodes','leiden_gene')]
  # leiden_barcodes_df$leiden_gene <- paste0(leiden_barcodes_df$leiden_gene," - ",sample_name)
  # leiden_barcodes_df <- leiden_barcodes_df %>%
  #   add_column(supercluster = rep("0",times=length(leiden_barcodes_df$leiden_gene)))
  # 
  # names(clusterAssignments) <- c("A","B","C","D","E","F")
  # 
  # for(i in 1:length(clusterAssignments)){
  #   leiden_barcodes_df$supercluster[leiden_barcodes_df$leiden_gene %in% clusterAssignments[[i]]] <- names(clusterAssignments)[i]
  # }
  # 
  # leiden_barcodes_df <- leiden_barcodes_df %>%
  #   filter(!supercluster == "0")
  # 
  # supercluster_barcodes <- leiden_barcodes_df$barcodes
  # spata_obj3 <- subsetByBarcodes_CountMtr(
  #   spata_obj2,
  #   supercluster_barcodes
  # )
  # 
  # # add as feature to spata_obj3
  # spata_obj3 <-
  #   addFeatures(object = spata_obj3,
  #               feature_names = c("supercluster"),
  #               feature_df = leiden_barcodes_df,
  #               overwrite = TRUE,
  #               key = "barcodes")
  # load premade spata_obj3-----
  spata_obj3 <- loadSpataObject(spata_obj_dir3, verbose = TRUE, update = TRUE)

  # plot superclusters----
  p4 <- plotSurface(object = spata_obj3, 
                    color_by = "supercluster",
                    display_image = TRUE,
                    #pt_clrp = "jama", 
                    pt_size = 1.8,
                    clrp_adjust = supercluster_colours) +
    ggplot2::labs(title = "Integrated Pathway-Based Clustering", color = "Consensus Clusters")
  
  # p2 <-
  #   plotSurface(object = spata_obj2,
  #               color_by = "leiden_gene",
  #               pt_size = 2.5,
  #               clrp_adjust = leiden_gene_colors,
  #               display_image = TRUE) +
  #   ggplot2::labs(title = "Non-integrated Gene-Based Clustering",color = "Leiden Clusters")
  
  # combine with patchwork 
  p2 + legendTop() +
    p4 + legendTop() 
  
  # ggplot2::ggsave('leidenAndSuperclusters_NES1FDRlow_map.png')
  ggplot2::ggsave('leidenAndSuperclusters_map.png')
  
  
  p4 <- plotSurface(object = spata_obj3, 
                    color_by = "supercluster",
                    display_image = TRUE,
                    #pt_clrp = "jama", 
                    pt_size = 1.8,
                    clrp_adjust = supercluster_colours) 
  p4
  
  # ggplot2::ggsave('superclusters_NES1FDRlow_map.png')
  ggplot2::ggsave('superclusters_map.png')
  maps[[sample_name]] <- p4 + theme(legend.position="none")
  
  # Load external gene sets -------
  #Load in a .csv where each row is a gene, and each column is a gene set name
  sigs_to_add=colnames(external_geneset)
  for (sig in sigs_to_add)
  {
    print(sig)
    spata_obj3<-addGeneSet(spata_obj3,class_name="gbmap",gs_name=sig,genes=dplyr::pull(external_geneset, sig),overwrite=TRUE)
  }
  
  # enrichment in superclusters----
  # spata_obj3 <-
  #      runDeAnalysis(object = spata_obj3,
  #                    across = "supercluster", # across which identity groups
  #                    method_de = c("wilcox") # with which methods
  #      )
    
  DeaDf <- getDeaResultsDf(object = spata_obj3,
                            across = "supercluster",
                            method_de = "wilcox",
                            max_adj_pval = 0.025)
  write.csv(DeaDf, paste(directory_10X,"\\",subdir,"\\supercluster_DeaDf.csv",sep=""),row.names = FALSE)
  
  #keywords_cellstates <- c("Neftel","Verhaak","Pugh","external_Ravi","external_SpOT.upreg","external_SpOT.downreg","external_Ron.upreg","external_Ron.downreg")
  keywords_gbmap <- c("gbmap")
  geneset_dataframe <- getGeneSetDf(spata_obj3)

  geneset_names <- geneset_dataframe %>%
    filter(str_detect(ont, paste(keywords_gbmap, collapse = "|"))) %>%
    pull(ont)
  geneset_names <- unique(geneset_names)

  #run GSEA
  # spata_obj3 <-
  #   runGsea(
  #     object = spata_obj3,
  #     gene_set_names = geneset_names,
  #     across = "supercluster",
  #     methods_de = "wilcox",
  #     reduce = FALSE
  #   )

  #get GSEA results
  gseaDf <- getGseaDf(
    object = spata_obj3,
    across = "supercluster",
    method_de = "wilcox"
  )
  #write.csv(gseaDf, paste(directory_10X,"\\",subdir,"\\superclustergbmap_gseaDf.csv",sep=""),row.names = FALSE)

  #GSEA dot plot
  p3 <- plotGseaDotPlot(
    object = spata_obj3,
    across = "supercluster",
    #across_subset = c("1", "2", "4", "5", "6"),
    color_by = "pval",
    #n_gsets = 5,
    pt_alpha = 0.8,
    #transform_with = list("fdr" = c("log10")),
    by_group = FALSE, # merge in one plot
    size_by = "fdr",
    pt_clrsp = "Inferno" ,
  )
  p3
  ggplot2::ggsave('superclusterGSEAdotplot_gbmap.png')
  
  #combine with patchwork
  p4 + legendRight() +
    p3 + legendRight()
  ggplot2::ggsave('supercluster_GSEAdotplot_gbmap.png')
  
  # Save spata_obj3 ---------
  spata_obj3 <-
    adjustDirectoryInstructions(
      object = spata_obj3,
      to = "spata_object",
      directory_new = spata_obj_dir3, # the directory under which to store the object
    )

  saveSpataObject(object = spata_obj3) #Save SpataObj


  # # convert spata_obj3 to Seurat object and append to list----
  # #obj_seurat <- transformSpataToSeurat(object = spata_obj3)
  # seurat_objs[[sample_name]] <- transformSpataToSeurat(object = spata_obj3, assay_name = "Spatial",)
}

ppath = "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\primaryIDHwt\\c2c5 gsea results"

library(gridExtra)
legend = supercluster_colours
pdf(paste0(ppath, "\\allSuperclusters.pdf")) # Open a new pdf file
grid.arrange(grobs = maps)
dev.off()

ggplot2::ggsave(paste0(ppath, "\\allSuperclusters.png"))


# "#DF536B" "#61D04F" "#2297E6" "#28E2E5" "#CD0BBC" "#F5C710"
# "#E69F00" "#56B4E9" "#009E73" "#F0E442" "#0072B2" "#D55E00" "#CC79A7"

#Seurat Neighbourhood Enrichment----
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)

# Load Visium H&E data
# transform each spata_obj3 to a Seurat object, then merge


# vst <- Load10X_Visium(data.dir = "path/to/visium/data", 
#                       assay = "Spatial", 
#                       img.format = "tif")

# Preprocessing steps
vst <- SCTransform(vst, assay = "Spatial")
vst <- RunPCA(vst, assay = "Spatial")
vst <- FindNeighbors(vst, dims = 1:30)
vst <- FindClusters(vst, resolution = 0.5)

# Compute spatial neighborhood enrichment
n_genes <- 20 # number of genes to use
top_genes <- head(VariableFeatures(vst, n = Inf)$gene, n_genes)
sne <- SpatialNeighborhoodEnrichment(vst, cluster = "RNA_snn_res.0.5", features = top_genes, method = "overlap")

# Visualize results
sne_df <- as.data.frame(sne$`RNA_snn_res.0.5`$enrichment)
sne_df$row_cluster <- FALSE
sne_df$col_cluster <- FALSE
sns <- sns::heatmap(sne_df, cmap = "coolwarm", xlab = "Leiden clusters", ylab = "Top genes", 
                    cex.axis = 0.8, cex.main = 0.9, margins = c(10, 10))

