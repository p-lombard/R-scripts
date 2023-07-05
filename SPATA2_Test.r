# https://themilolab.github.io/SPATA2/
# Install packages -----------------
install.packages("devtools")

if (!base::requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'EBImage'))

install.packages("Seurat")

devtools::install_github(repo = "kueckelj/confuns")
devtools::install_github(repo = "theMILOlab/SPATA2")

# if you want to use monocle3 related wrappers
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
# Default installation misses hdf5r which is required to read the filtered feature matrix
install.packages("hdf5r")
install.packages('pracma')
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
install.packages("viridis")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")
install.packages("viridisLite")
install.packages("msigdbr")


# Declare constants ---------
library(SPATA2)
library(magrittr)
library(tidyverse)
library(monocle3)
library(viridisLite)

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
#directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/V21004-FF__63_D1"

# MILO Lab samples:
directory_10X = "C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/242_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/243_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/248_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/251_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/255_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/256_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/259_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/260_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/262_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/269_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/275_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/296_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/304_T"
# directory_10X="C:/Users/mzaidi/Documents/Visium Datasets/SPATA2 Data/334_T"


#sample_name <- substr(directory_10X, 67, nchar(directory_10X))
sample_name <- substr(directory_10X, 55, nchar(directory_10X))
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
# Create spata object -----
#Note, this also runs various preprocessing, DR, and clustering techniques too.
spata_obj <-
  initiateSpataObject_10X(
    directory_10X = paste(directory_10X,"/outs",sep=""), # the directory from which to load the data
    sample_name = sample_name,
    directory_spata = spata_obj_dir
  )


# # Or alternatively, load premade spata obj --------
# spata_obj <- loadSpataObject(spata_obj_dir2, verbose = TRUE, update = TRUE)
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
for (sig in sigs_to_add)
{
  print(sig)
  genes=dplyr::pull(external_geneset, sig)
  genes = setdiff(genes, c('SGOL1-AS1', 'MNX1', 'LOC102546294', 'LOC102723885', 
                           'KIAA0226L', 'C2orf48', 'LMAN1L', 'CD207', 'MIR3145', 
                           'GAREML', 'KRT18P55', 'PI4KAP1', 'LOC148709', 
                           'LOC101928973', 'LINCR-0001', 'LIMD1-AS1', 'PRSS50', 
                           'TMSB15A', 'MIR4697HG', 'SMIM17', 'FAM150B', 'MIR1268A', 
                           'MIR6804', 'C19orf80', 'LOC101928509', 'NAMA', 'PNMAL2', 
                           'LOC101927292', 'LOC100128288', 'BCRP3', 'MIR6871', 
                           'ALOX15P1', 'LOC728392', 'MIR573', 'PFN1P2', 'C12orf79', 
                           'MAP6D1', 'PROX1-AS1', 'LINC01250', 'MIR548AR', 
                           'LOC100996251', 'C1orf210', 'LOC100506725', 'LOC100506497', 
                           'LOC100130476', 'HAND2-AS1', 'C4orf45', 'PRSS46', 
                           'LOC100129083', 'FAM35BP', 'LOC101927391', 'LINC01279', 
                           'FLJ46906', 'ASAH2', 'MIR378I', 'LOC101926935', 
                           'TBC1D26', 'SNORA47', 'MANEA-AS1', 'SEPT5-GP1BB', 
                           'LOC389705', 'CA6', 'LOC101929224', 'STX16-NPEPL1', 
                           'RPL34-AS1', 'AQP7P1', 'ZNF157', 'LINC00437', 
                           'LOC102723927', 'FAM105A', 'ZAR1L', 'HS1BP3-IT1', 
                           'SNORA76A', 'LINC00377', 'SNORA74B', 'LOC100132356', 
                           'LOC286297', 'LINC00202-1', 'LOC100288152','SFTPA2', 
                           'ZNF33BP1', 'LOC101927248', 'LOC339539', 'LOC101927043', 
                           'CHRM4', 'LOC100128573', 'CRHR1-IT1', 'MIR6757', 
                           '02-Mar', 'ZNF32-AS2', 'PCDHGC4', 'SNORA11', 'LOC729506', 
                           'C15orf26', 'FAM134B', 'TMEM180', 'LOC102724596', 'CBY3', 
                           'PCDHAC2', 'DBNDD2', 'PCBP2-OT1', 'LOC101928885', 'TMEM35', 
                           'CEBPA-AS1', 'RAPGEF4-AS1', 'BSN-AS2', 'LOC729732', 
                           'LOC101928766', 'MIR1276', 'HNRNPU-AS1', 'MIR4312', 
                           'EGFEM1P', 'MIR548AQ', 'LOC100507250', 'MIR548AI', 
                           'PDCD6IPP2', 'MST1L', 'LOC101928283', 'TOB1-AS1', 
                           'LOC339529', 'LOC728743', 'HERC2P7', 'SNORA25', 
                           'TEX15', 'MIR3662', 'LOC400927', 'ZNF37BP', 'MFSD4', 
                           'LOC100130705', 'SLC5A1', 'TEX19', 'MIR3125', 'MIR7851', 
                           'LOC101928847', 'ANP32A-IT1', 'FLJ37453', 'MRPL42P5', 
                           'RNU5E-1', 'LOC101927394', 'OSER1-AS1', 'LOC100507073', 
                           'ABCA17P', 'GLRA4', 'LOC401127', 'SCARNA2', 'LINC00282', 
                           'PCDHA9', 'LOC101929679', 'SPRY4-IT1', 'TMEFF1', 'LINC01126', 
                           'LOC440311', 'MIR5581', 'TYMSOS', 'MGC2889', 'TTC3P1', 
                           'PCDHGB2', 'LOC101929541', 'GS1-24F4.2','SSSCA1-AS1', 
                           '04-Sep', 'ZNF559-ZNF177', 'FAM27E2', 'LAMTOR5-AS1', 
                           'SNORA54', 'MIR641', 'LOC100272216', 'RPS6KA2-IT1', 
                           'FLJ26850', 'RNF224', 'CHST5', 'LOC100506127', 'RNU11', 
                           'LINC00597', 'LOC284023', 'LOC101593348', 'ACTG1P17', 
                           'RNU4ATAC', 'MIR4324', 'LOC100134368', 'C2orf66', 
                           'MST1P2', 'LOC101928063', 'MIR92B', 'MIR3685', 
                           'LOC100133286', 'MMP23A', 'TMEM133', 'FAM212B-AS1', 
                           'MIR548AN', 'LENEP', 'VN1R1', 'TAS2R50', 'PCDHB19P', 
                           'KIAA1715', 'DLEU2L', 'LOC100506834', 'PARK2', 
                           'RNU2-2P', 'SNORD59A', 'LOC100288181', 'FOXE3', 
                           'LINC00936', 'LOC101927164', 'SNORD4B', 'EGR1', 
                           'PPAP2B', 'MIR3613', 'SNORD25', 'MIR5010', 'LOC643072', 
                           'LOC101928324', 'ZNF197-AS1', 'ALG1L9P', 'PLA2G4B', 
                           'FAM35DP', 'CECR5-AS1', 'LOC100129138', 'LOC100507194', 
                           'GUCA1B', 'LOC100294362', 'HOXD4', 'AMZ2P1', 
                           'N4BP2L2-IT2', 'LINC00673', 'LOC728485', 'LOC403312', 
                           'LOC100128568', 'PTENP1', 'ASB16', 'SNORA74A', 'RPS6KA2', 
                           'CDIPT-AS1', 'TAS2R3', 'SRD5A3-AS1', 'LOC643770', 
                           'C1orf234', 'GBAS', 'LOC100129973','GBAT2', 'ADCK3', 
                           'LINC00641', 'NUP50-AS1', 'SNORD59B', 'LINC00618', 
                           'C10orf131', 'ALMS1-IT1', 'LINC00515', 'MIR4730', 
                           'EFTUD1P1', 'IPO5P1', 'CBX3P2', 'LOC101409256', 
                           'C20orf196', 'HIST1H1E', 'LOC102724009', 'RP9P', 
                           'DBIL5P2', 'CLK2P1', 'DPY19L1P1', 'LOC100289230', 
                           'LOC102723703', 'LOC101928279', 'SUV420H2', 'SNORD5', 
                           'LOC730183', 'LOC100506022', 'ZNF528-AS1', 'SNORD95', 
                           'BRWD1-AS1', 'ZNF705E', 'GBAP1', 'CFL1P1', 
                           'LOC105376671', 'LINC00173', 'ZNRF2P2', 'CBR3-AS1', 
                           'PCBP1-AS1', 'LOC283440', 'DPY19L2P3', 'RN7SK', 
                           'CELF6', 'LOC100132057', 'ZNF815P', 'ZNF833P', 'PNMAL1', 
                           'MIR5087', 'LOC101928438', 'LOC401320', 'C1orf106', 
                           'LOC101927768', 'LOC100128398', 'THCAT158', 
                           'LOC100505666', 'KLK14', 'C1QTNF3-AMACR', 
                           'LOC101928222', 'LINC01530', 'SNORA8', 'LOC100270804', 
                           'SLMO1', 'MIR1271', 'LOC100134868', 'LOC101927765', 
                           'AKR7L', 'NIPBL-AS1', 'LOC728730', 'NDUFA6-AS1', 
                           'ERI3-IT1', 'RRP7BP', 'LOC100506472', 'ZNF192P1', 
                           'PTGES2-AS1', 'SNORA5A', 'RALGAPA1P', 'LOC105373383',
                           'LOC101929577', 'LOC646762', 'SNORD35B', 'LINC00938', 
                           'LINC01604', 'GLTSCR2', 'C17orf96', 'LINC01578', 
                           'MFI2-AS1', 'LOC730202', 'LOC102724919', 'LOC101928053', 
                           'NUTM2B', 'CROCCP3', 'MIR600HG', 'SMA4', 'RPL13AP5', 
                           'LOC100507557', 'LINC00441', 'LOC101929231', 'PDCL3P4', 
                           'LOC644656', 'CROCCP2', 'LOC100289495', 'SNORA1', 
                           'LOC101928812', 'LOC101927811', 'SDHAP3', 'KIAA1467', 
                           'AHSA2', 'LOC100506125', 'LOC642846', 'LOC100288778', 
                           'FAM132B', 'MIR4712', 'GNG12-AS1', 'LOC100506124', 
                           'DFNB59', 'LOC100505912', 'C9orf41-AS1', 'SH3GL1P1', 
                           'HRAT5', 'LOC100506603', 'LOC283140', 'ASAP1-IT2', 
                           'DUSP5P1', 'SNORD87', 'C9orf173-AS1', 'MIR1914', 
                           'GUSBP9', 'LOC100507283', 'SNORD31', 'TMEM27', 'ZNF767P', 
                           'DKFZP434I0714', 'LCMT1-AS1', 'RAB40AL', 'LOC100506548', 
                           'MIR4720', 'SMG1P7', 'LINC00174', 'LOC389765', 
                           'SNORA71E', 'ZNF878', 'KIAA1407', 'LINC00869', 'ADCK4', 
                           'GCSHP3', 'CCDC101', 'GATSL2', 'KLHL11', 'FAM27E3', 
                           'SNORD22', 'SCARNA17', 'SNORD56', 'SNORD26', 
                           'AGPAT4-IT1', 'ZNF724P', 'NCK1-AS1','RAD51-AS1', 
                           'ZNF460', 'GTF2IP1', 'GUSBP2', 'ATP5D', 'TMEM198B',
                           'LINC00969', 'NSUN5P2', 'LOC100506730', 'SNORD42B', 
                           'LOC150776', 'LOC100507424', 'SNORA5C', 'LOC81691', 
                           'GLUD1P7', 'LINC01590', 'DDX12P', 'C5orf45', 'KLRAP1', 
                           'RPS11', 'SUZ12P1', 'SNORD97', 'RPL18A', 'LOC105373300', 
                           'LOC101926911', 'CECR5', 'RPS15AP10', 'WDFY3-AS2', 
                           'WASH2P', '08-Sep', 'MIR4517', 'LOC728752', 'KLF3-AS1', 
                           'SCARNA15', 'FAM195A', 'MLLT4', 'SEPW1', 'LOC101929066', 
                           'MGC16142', 'HERC2P9', 'SNORA32', 'SNHG20', 
                           'SLC25A5-AS1', 'TTC28-AS1', 'NPHP3-ACAD11', 
                           'LINC-PINT', 'ADAM1A', 'GNRHR2', '07-Sep', 'FLJ10038', 
                           'FAM195B', 'PMS2CL', 'MGC27345', 'LOC554206', 
                           'LOC100130691', 'C7orf55', 'LOC101928068', 
                           'LOC101929767', 'LOC100507291', 'KIAA0195', 'MIR3960', 
                           'SNORD96A', 'TEN1-CDK3', 'TPT1-AS1', 'DKFZP586I1420', 
                           'ILF3-AS1', 'GUSBP11', 'LOC101927415', 'ARHGAP31-AS1', 
                           'MGC57346', 'MIR7-1', 'PARGP1', 'ARMCX5-GPRASP2', 
                           'MTRNR2L2', 'C14orf80', 'LOC100289333', 'MIR25', 
                           'POLR2J4', 'ZBTB12', 'LOC100419583','LOC220729', 
                           'NGFRAP1', '06-Mar', 'MIR4292', 'GLUD1P3', 'FLJ37035', 
                           'PP7080', 'MIR3916', 'LOC101929715', 'MIR570', 'TMEM183B', 
                           'SEPT7P2', 'BMS1P5', 'PKI55', 'NPTN-IT1', 'C2orf44', 
                           'LOC100130987', 'MSTO2P', 'ZNF137P', 'SDHAP1', 'ZNRD1-AS1', 
                           'LOC90784', 'LOC101929140', 'PRUNE', 'ZNF625-ZNF20', 
                           'FAM134A', 'LOC654841', 'WASH1', 'FOXO3B', 'MRPS31P5', 
                           'LINC01000', 'SNORD4A', 'LOC103344931', 'JUN', 'ATP5A1', 
                           'NICN1', 'LRRC37BP1', 'FRG1HP', 'CRAMP1L', 'RPL32P3', 
                           'PMS2P3', 'KIAA0226', 'SELK', 'ATP5B', 'LOC100506639', 
                           'LOC645513', 'LOC283922', 'SUV420H1', 'C5orf42', 'KIAA0141', 
                           'LDOC1L', 'DSCR3', 'ATP5SL', 'ZCCHC11', 'LOC93622', 
                           'C3orf17','HIAT1', 'RPL7L1', 'C2orf47', 'CIRH1A', 'KIAA0196', 'KIAA1429', 'RRN3P3', 'KIAA1033', 'VIMP', 'RPL36AL', 'ZZZ3', 'LOC652276', 'C9orf69', 'TCEB1', 'RPS17', 'C14orf1', 'SHFM1', 'HN1L', 'C17orf59', 'ASUN', 'CCDC53', 'C20orf24', 'CRYBB2P1', 'LOC101929147', 'APITD1', 'FAM175A', 'BRE', 'MINA', 'LOC101927204', 'ATP8B5P', 'LOC100129461', 'SPG20', 'PDF', 'RPSAP9', 'PYCRL', 'C17orf89', 'MIR635', 'MESDC2', 'SEC14L1P1', 'ATP5EP2', 'C10orf12', 'SLC16A13', 'LOC389247', 'MIR7641-2', 'MESDC1', 'PDIA3P1', 'LINC00883', 'LINC01420', 'KIAA0368', 'FAM103A1', 'PPAP2A', 'E2F8', 'OBFC1', 'LHFP', 'PLEKHA8P1', 'DTX2P1-UPK3BP1-PMS2P11', 'ABCC2', 'SNORD101', 'DDX26B', 'RPS10', 'C16orf52', 'PVRL2', 'METTL21B', 'SLC7A5P2', 'TOB2P1', 'QTRTD1', 'FLJ23867', 'HLA-L', 'SNORD17', 'FLJ20021', 'SRPR', 'KIAA0754', 'LOC100288798', 'SNORD43', 'SNORA73B', 'FAM127C', 'LOC284454', 'RPS6KA4', 'GNAT2', 'LRRC48', 'RPL17', 'LOC100128361', 'SNORD50B', 'KIAA0020', 'FAM86DP', 'LINC01160', 'ZCCHC6', 'CTB-113P19.1',
                           'GRAMD3', 'FTO-IT1', 'SNORA73A', 'SNORD83A', 'ZAK', 'C8orf31', '06-Sep', 'KIAA1656', 'N6AMT2', 'LINC01057', 'TMEM2', 'PRR7-AS1', 'LOC100288911', 'LVCAT8', 'GALNT4', 'SUGT1P1', 'BRE-AS1', 'SSR4P1', 'LOC102724927', 'LOC101927040', 'FLJ13224', 'BBOX1-AS1', 'PPP3CB-AS1', 'ZMYM6NB', 'RALY-AS1', 'LOC101927157', 'DLG5-AS1', 'LOC100507487', 'ANKRD26P1', 'LINC00152', 'GABARAPL3', 'LBHD1', 'HIST2H2BA', 'WBSCR27', 'MIR616', 'LINC00933', 'DGCR11', 'C15orf52', 'MIR24-1', 'LOC101927056', 'PCDHGB5', 'RPS6KA3', 'MIR4785', 'LOC100507431', 'TNRC18P1', 'HIST1H2BK', 'RPS16P5', 'LOC100507156', 'POTEE', 'SLCO5A1', 'LOC644838', 'LOC100506990', 'MIR4263', 'LHFPL3-AS2', 'DENND5B-AS1', 'LINC00899', 'LOC101929579', 'EPGN', 'ABALON', 'HTR7P1', 'LOC101929295', 'PVRL3', 'MAGEB17', 'LOC100499489', 'KIRREL', 'MMP24-AS1', 'LOC101927322', 'MIR27B', 'LOC646626', 'LOC102724094', 'HLA-J', 'C9orf152', 'AOX2P', 'C14orf169', 'LOC100507195', 'LOC101928445', 'LOC388813', 'LOC101929441', 'RNF217-AS1',
                           'LOC101927588', 'LINC00312', 'CTD-2201I18.1', 'HIST1H2AI', 'MIR4767', 'LOC101928650', 'DEPDC1-AS1', 'LOC145783', 'MEIS3P1', 'THEGL', 'LOC100506100', 'LOC100130331', 'C10orf10', 'FAM65B', 'DCDC5', 'LOC645166', 'LOC101928858', 'RNU6-35P', 'HIST1H3D', 'LINC00563', 'SLC7A5P1', 'LOC285043', 'LOC105616981', 'MAGOH2P', 'MGC72080', 'VWA7', 'PLA2G12B', 'KCNQ5-IT1', 'CMAHP', 'WBSCR28', 'MIR548F3', 'GCM1', 'GTSF1L', 'TPTE2P6', 'KLKP1', 'SNX29P2', 'LINC01489', 'LOC101928674', 'HYPM', 'LOC100652768', 'LOC101928782', 'MIR4685', 'LOC100049716', 'MIR4731', 'TMEM191A', 'MFI2', 'FSCN2', 'LOC101928489', 'LINC01599', 'LOC101928794', 'MIR4450', 'LOC100131496', 'LOC100506368', 'C1orf204', 'SLC9A2', 'SLC6A3', 'FAM71E2', 'SERHL', 'LOC101928034', 'MIR6819', 'HBE1', 'LMNTD1', 'LINC01305', 'POTEF', 'MIR23A', 'GUSBP5', 'LOC730102', 'DIAPH2-AS1', 'LOC101929555', 'LINC00856', 'BLID', 'KCCAT333', 'RGAG1', 'LINC00536', 'LOC641746', 'FAM109B', 'C8orf4', 'SYNJ2-IT1', 'AFAP1-AS1', 'GRAMD2', 'HMCN2',
                           'MIR1206', 'XG', 'TMEM190', 'CCDC178', 'CRYBA1', 'LOC105375734', 'SBK2', 'PKP3', 'LOC100507144', 'LOC145845', 'SLC9A4', 'XIRP2-AS1', 'LMO7DN-IT1', 'LOC102724084', '10-Mar', 'MUC22', 'LOC101929762', 'ESR2', 'LOC101927934', 'TNRC6C-AS1', 'LOC100507002', 'SVILP1', 'STON1-GTF2A1L', 'ENTHD1', 'KIAA1804', 'DMBX1', 'LOC730668', 'LOC100240735', 'LOC100287728', 'LINC00398', 'C14orf105', 'TICAM2', 'DPPA4', 'TRABD2B', 'C10orf67', 'LOC102723824', 'MBL1P', 'LOC101928436', 'PNLIPRP3', 'AKAP2', 'CSRP3', 'APCDD1L-AS1', 'LINC01392', 'C12orf56', 'ANKRD20A5P', 'LOC101928880', 'FAM83C', 'OGFRP1', 'FAM155A-IT1', 'LOC442497', 'LOC101929596', 'LOC100287225', 'IL21-AS1', 'NBPF18P', 'LRRC37A6P', 'FLJ43879', 'EMBP1', 'ATP6V0D2', 'LOC101927740', 'LOC101927911', 'LOC100507351', 'ACTL8', 'NALCN-AS1', 'KU-MEL-3', 'CPXCR1', 'CTD-2297D10.2', 'LOC100506860', 'LOC284080', 'LINC00552', 'DPRXP4', 'LINC00942', 'LOC101926955', 'DUSP1', 'ZNF702P', 'MGAM2', 'FNDC8', 'METTL24', 'ZFPM2-AS1', 'LINC01031',
                           'MUC17', 'DPCR1', 'LOC101927630', 'TRIML2', 'DSC1', 'IL17B', 'KANK4', 'ANKRD30B', 'DLG1-AS1', 'RNF144A-AS1', 'LINC00578', 'LOC389641', 'RNU6-26P', 'LOC100128233', 'LOC100507391', 'MIR2116', 'LOC101926889', 'FER1L6-AS2', 'HCG26', 'PTRF', 'VTRNA2-1', 'RNASE7', 'SLC1A7', 'ATP6V0A4', 'COL6A6', 'LOC101929294', 'CCDC63', 'EPYC', 'AIM1L', 'MB21D1', 'PSORS1C3', 'LOC643339', 'ANKRD18B', 'RPS6KA6', 'CCDC64B', 'TMEM26-AS1', 'LOC102724933', 'AGAP11', 'POM121L10P', 'C2orf78', 'MUC2', 'MUC13', 'FAM46C', 'CHAT', 'LHX5', 'HYMAI', 'CATIP-AS1', 'LINC00521', 'MIR4664', 'LINC00592', 'SAA2', 'PAX8-AS1', 'C17orf102', 'ABHD11-AS1', 'GJB3', 'LINC01594', 'C15orf32', 'MIR2682', 'LOC101929470', 'ANXA2P1', 'LOC101928841', 'PKD1L2', 'WFDC10B', 'LINC01589', 'LOC200772', 'FER1L4', 'ANXA2P2', 'LIN28B', 'LINC00944', 'LINC00961', 'MT1L', 'LOC149684', 'LOC100505478', 'CCL1', 'GVINP1', 'ADAMTS16', 'LOC541472', 'C1orf140', 'LOC102546229', 'LOC101927476', 'KRT32', 'PHKA2-AS1', 'FLJ31356', 'LOC344887',
                           'ST8SIA2', 'MIR6730', 'FYB', 'PMCH', 'LOC100132111', 'LOC729987', 'HCP5', 'LOC101928569', 'CCDC185', 'LOC101928891', 'ENAM', 'LOC400548', 'ANKRD63', 'LOC101927780', 'TEX41', 'LOC102467080', 'DGAT2L6', 'ANXA2P3', 'TARID', 'LOC100507639', 'AADACP1', 'FAM26F', 'KCNU1', 'LOC401585', 'LOC101927653', 'LINC01501', 'BTG4', 'ILDR1', 'SELM', 'CCDC70', 'MMP20', 'GSN-AS1', 'GCNT3', 'LOC440028', 'LINC01581', 'HOXB9', 'NLGN4Y-AS1', 'IL24', 'MMP8', 'CFHR1', 'SPERT', 'CTAGE11P', 'PGAM1P5', 'ALS2CR11', 'MIR221', 'LINC01048', 'LINC00341', 'LINC01393', 'LINC01522', 'FAM180A', 'AIM1', 'HNF4A', 'RAET1G', 'LOC643201', 'MIR6753', 'LOC101927166', 'KRT9', 'RPLP0P2', 'PKP1', 'HTATSF1P2', 'RAD21L1', 'AIPL1', 'C4orf22', 'LOC100507334', 'LOC101927229', 'FLJ41941', 'MROH9', 'GJD3', 'LOC101929723', 'KPRP', 'C15orf56', 'LINC00626', 'KRT39', 'LOC101929427', 'LOC101927560', 'MEOX1', 'GRPR', 'FAM196B', 'CPA1', 'CCDC60', 'LOC283731', 'RPS10P7', 'AGR2', 'PSG5', 'SBK3', 'LINC01151', 'ROPN1', 'SPATA31E1',
                           'LINC00605', 'UGT3A2', 'LOC729739', 'IFNB1', 'LOC100126784', 'DCDC2C', 'CCKAR', 'HEPHL1', 'MTL5', 'LOC389332', 'LOC115110', 'GATA6-AS1', 'DSG3', 'LOC101928994', 'LOC100129434', 'C7orf69', 'LOC100131289', 'LOC100287036', 'SPANXC', 'LINC00243', 'LINC00311', 'FER1L6', 'MIR196A1', 'LOC389033', 'SPANXB1', 'LOC101928358', 'GAS6-AS2', 'LOC285626', 'OR7E14P', 'SIGLEC15', 'GBA3', 'KRT79', 'GP1BA', 'AICDA', 'ST6GALNAC2', 'LINC01040', 'SLC2A2', 'CR1L', 'MIR4645', 'LOC101927523', 'SCGB3A2', 'PRKCDBP', 'IL36RN', 'IGSF23', 'LRRN4', 'TCN1', 'IFNE', 'KRT80', 'LOC90246', 'MYL1', 'LOC100507477', 'LOC101928461', 'HOXB5', 'C11orf86', 'C10orf54', 'HCRTR2', 'LINC01111', 'LOC101929199', 'FAM163A', 'CHRNB3', 'GPR20', 'LINC00954', 'C7orf65', 'LOC643711', 'MRGPRF-AS1', 'LOC105747689', 'SLC5A8', 'LINC00475', 'CFLAR-AS1', 'C15orf54', 'CGB8', 'ZDHHC8P1', 'NTSR1', 'KRT14', 'HOXB8', 'GPC6-AS1', 'RUNX1-IT1', 'EPN3', 'LOC100506188', 'TRIM43B', 'CPN2', 'C1orf110', 'SPTA1', 'BLK', 'PSG9', 'TMEM40',
                           'PLCE1-AS1', 'HYAL4', 'BPI', 'CCL24', 'LOC101929710', 'ZBED2', 'GATA6', 'KRTAP2-3', 'LAD1', 'TRAPPC3L', 'MIR222', 'MUC6', 'CRLF2', 'HAND1', 'GPR87', 'HTR3C', 'MC4R', 'GSTA1', 'ROR1-AS1', 'TMEM155', 'LINC01444', 'C16orf47', 'RAET1L', 'PRKG1-AS1', 'LINC01050', 'LRTM1', 'KCCAT198', 'SLC6A20', 'LUCAT1', 'LOC100996351', 'GATA3-AS1', 'EXOC3L4', 'KRT13', 'TYRP1', 'LOC101060542', 'FBXL21', 'SELP', 'MYPN', 'LGALS17A', 'PAEP', 'ANKRD1', 'LOC284344', 'CECR1', 'MUC5AC', 'GBP1P1', 'LOC101926940', 'LOC643733', 'OR2S2', 'CCL11', 'KRT4', 'C11orf44', 'MAGEC2', 'PADI3', 'EGOT', 'KIRREL3-AS3', 'BNC1', 'SLCO4C1', 'KLHL38', 'TMPRSS11E', 'MFSD7', 'C5orf46', 'MYH16', 'SLC9A7P1', 'FOXA2', 'MAGEC1', 'C20orf166-AS1', 'FAM222A-AS1', 'EML2-AS1', 'IL36B', 'LOC102724434', 'LOC100507600', 'PLEKHS1', 'LOC101929106', 'ELF5', 'MIR146A', 'POM121L9P', 'GIP', 'CST2', 'PRSS22', 'SERPIND1', 'SPINK1', 'LOC101927851', 'IRS4', 'LINC01204', 'HOXB-AS3', 'RASSF6', 'LINC01272', 'MIR31HG', 'CPZ', 'ESRP1', 'HP',
                           'CCDC67', 'LINC00525', 'FOXI1', 'SLFN12L', 'CD5L', 'WTAPP1', 'LOC100505622', 'ABCB5', 'MYOCD', 'CTSE', 'SIM1', 'PLEKHN1', 'TMEM215', 'RPSAP52', 'SDPR', 'STYK1', 'ZP4', 'MIR614', 'IL31RA', 'PAPL', 'TGM4', 'CBLC', 'THRA1/BTR', 'TMEM239', 'LECT1', 'SDR42E1', 'LOC101927230', 'SH3PXD2A-AS1', 'TMPRSS2', 'ADH1C', 'LINC01179', 'LINC01605', 'CSF2', 'ANO1-AS2', 'PSG4', 'IL1RL1', 'NLRP10', 'CCL7', 'PI3', 'CSF3', 'LINC00520', 'AP1M2', 'LINC00452', 'LOC100507065', 'FLJ22447', 'SAA1', 'ROS1', 'LOC101927356', 'FGF5', 'LCN2', 'SPINK5', 'ACTBL2', 'LOC284798', 'LOC101927482', 'NUDT16P1', 'SELE', 'LRRC15', 'GPR101', 'AREG', 'LOC101927884', 'ADRA1B', 'FAM26E', 'IGF2BP1', 'LINC01611', 'LINC00968', 'MMP10', 'C20orf141', 'TNIP3', 'FAM65C', 'HRAT17', 'LHCGR', 'TCF21', 'ITIH3', 'KRT15', 'KRT19', 'LOC152225', 'WFDC21P', 'SERPINB7', 'MMP13', 'TRHDE-AS1', 'ARHGEF34P', 'LOC100505817', 'FRMD7', 'EREG', 'IDO1', 'KCCAT211', 'PCSK9', 'XDH', 'SLC28A3', 'LOC100128317', 'FGFBP1', 'LINC00052', 'C4orf32',
                           'HSD17B2', 'ZCCHC5', 'LINC00704', 'MMP3', 'LINC00857', 'MUC16', 'MMP1', 'LINC00261',
                           '09-Sep', 'SEPTIN3','03-Sep','C19orf22', 'CD97', 'CENTD1', 'CENTD3', 'FLJ11286', 'FLJ21963', 'GPR56', 'LRRC16', 'M6PRBP1', 'MGC2752', 'NOS2A', '11-Sep', 'ZNF228', 'POLR1G',
                           'C1orf38', 'CCDC109B', 'DRAM', 'ECGF1', 'FER1L3', 'FLJ20273', 'FLJ22662', 'FPRL2', 'GLT25D1', 'LYPLA3', 'PGCP', 'PSCD4', 'PSCDBP', 'SQRDL',
                           'AGXT2L1', 'ATP5F1', 'ATP5L', 'EDG1', 'FLJ22655', 'KIAA1598', 'LOC201229', 'ORC4L', 'PSCD1', 'SEPP1', 'ZNF323',
                           'AOF2', 'BAI3', 'C20orf42', 'C6orf134', 'CAMSAP1L1', 'FAM125B', 'FAM77C', 'GPR23', 'HN1', 'KIAA1166', 'KLRK1', 'LOC55565', 'LPHN3', 'MYB', 'MYST2', 'PCDH11X', 'PCDH11Y', 'PHF16', 'PHLPP', 'RP11-35N6.1', 'SNX26', 'TMEM118', 'TMSL8', 'TNRC4', 'WDR68', 'ZNF643',
                           'ATF3', 'ERO1L', 'RPL21','LOC150568'
                           ))
  spata_obj<-addGeneSet(spata_obj,class_name="external",gs_name=sig,genes=genes[genes != ""],overwrite=TRUE)
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
###---- 

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
  

## Generate user defined trajectory lines figures ----
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
# pythonObsClusters=paste(directory_10X,"/outs/obsClusters.csv",sep='')
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
  ggplot2::labs(title = "Scanpy leiden clusters (low)", color = "Scanpy Leiden")
ggplot2::ggsave('scanpy_leidenLow.png')

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
write.csv(gseaDf, paste(directory_10X,"\\",subdir,"\\cellstates_gseaDfLow.csv",sep=""),row.names = FALSE)

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
  ggplot2::labs(title = "Non-integrated Gene-Based Clustering",color = "Leiden Clusters")

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

#Trying GSEA again with clusterprofiler----
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

#replicating ravi supplemental fig 4 cnet plots
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



# Plot Garofano Superclusters----
# load cluster supercluster assignments
ppath = "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\whole_GBM_dataset_analysis\\c2c5 gsea results"
load(paste0(ppath, "\\clusterAssignments_NESpos_FDRlow.RData"))

# turn this into a dataframe shape we can add to spata_obj2 (2 columns, 
# one representing barcodes, one representing cluster label)
leiden_barcodes_df = getFeatureDf(spata_obj2)
leiden_barcodes_df <- leiden_barcodes_df[,c('barcodes','leiden_gene')]
leiden_barcodes_df$leiden_gene <- paste0(leiden_barcodes_df$leiden_gene," - ",sample_name)
leiden_barcodes_df <- leiden_barcodes_df %>%
  add_column(supercluster = rep("0",times=length(leiden_barcodes_df$leiden_gene)))

names(clusterAssignments) <- c("A","B","C","D","E","F")

for(i in 1:length(clusterAssignments)){
  leiden_barcodes_df$supercluster[leiden_barcodes_df$leiden_gene %in% clusterAssignments[[i]]] <- names(clusterAssignments)[i]
}

# add as feature to spata_obj2
spata_obj2 <- 
  addFeatures(object = spata_obj2, 
              feature_names = c("supercluster"), 
              feature_df = leiden_barcodes_df,
              overwrite = TRUE,
              key = "barcodes")

p4 <- plotSurface(object = spata_obj2, 
            color_by = "supercluster",
            display_image = TRUE,
            #pt_clrp = "jama", 
            pt_size = 1.8) +
  ggplot2::labs(title = "Integrated Pathway-Based Clustering", color = "Consensus Clusters")

# combine with patchwork 
p2 + legendTop() +
  p4 + legendTop() 

ggplot2::ggsave('superclusters_map.png')
  

