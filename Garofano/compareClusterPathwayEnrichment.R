#https://github.com/lucgar/GBMstates/blob/master/functions.R
# In this script we want to compare gene set enrichment results of all our clusters across all samples
# The code is based on workflow outlined in
# Garofano, L., Migliozzi, S., Oh, Y.T. et al. "Pathway-based classification of glioblastoma 
# uncovers a mitochondrial subtype with therapeutic vulnerabilities". Nat Cancer 2, 141–156 (2021). 
# https://doi.org/10.1038/s43018-020-00159-4
#
# Garofano github: https://github.com/lucgar/GBMstates
# Relevant R scripts from the github are saved at:
# 1. C:/Users/mzaidi/Documents/R Scripts/Garofano/Garofano_Figure1.Rmd -- an R markdown notebook showing how to create the plots
#    for figure 1 in their paper
# 2. C:/Users/mzaidi/Documents/R Scripts/Garofano/functions.R -- R file containing the functions they used to perform analysis.

# # Install Packages----
library(devtools)
# install_github("miccec/yaGST")
# BiocManager::install("ConsensusClusterPlus")
# Set Constants----
source("C:\\Users\\mzaidi\\Documents\\R Scripts\\Garofano\\functions.R")

# they use this library for parallelization. Not compatible with Windows 
#so we could try to find an alternative. Don't really need it for the parts of their workflow we are replicating though. 
#library(doMC)

# this is a library from another one of their lab github pages (https://github.com/miccec/yaGST)
library(yaGST)
library(msigdbr)
# import a gmt file containing all the pathways you want to interrogate
c2 <- gmt2GO("C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Gene signatures\\c2.all.v2023.1.Hs.symbols.gmt")
c5 <- gmt2GO("C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Gene signatures\\c5.all.v2023.1.Hs.symbols.gmt")
geneSet <- c(c2,c5)

# load NES_allClusters, pValue_allClusters, FDR_allClusters. Each dataframe consists 
# of rows = pathways, columns = cluster+sample ID, values = respective measurements.

# myfiles = c("C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\whole_GBM_dataset_analysis\\c2c5 gsea results\\NES_primarywt.csv",
#            "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\whole_GBM_dataset_analysis\\c2c5 gsea results\\pval_primarywt.csv",
#            "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\whole_GBM_dataset_analysis\\c2c5 gsea results\\fdr_primarywt.csv")

# myfiles = c("C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\primaryIDHwt\\c2c5 gsea results\\NES_all.csv",
#             "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\primaryIDHwt\\c2c5 gsea results\\pval_all.csv",
#             "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\Sheila_GBM\\GBM_datasets\\primaryIDHwt\\c2c5 gsea results\\fdr_all.csv")

# myfiles = c("C:\\Users\\mzaidi\\Documents\\Visium Datasets\\SPATA2 Data\\c2c5 gsea results\\NES_all.csv",
#             "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\SPATA2 Data\\c2c5 gsea results\\pval_all.csv",
#             "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\SPATA2 Data\\c2c5 gsea results\\fdr_all.csv")

myfiles = c("C:\\Users\\mzaidi\\Documents\\Visium Datasets\\SPATA2 Data\\c2c5 MILO only gsea results\\june12\\NES_all.csv",
            "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\SPATA2 Data\\c2c5 MILO only gsea results\\june12\\pval_all.csv",
            "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\SPATA2 Data\\c2c5 MILO only gsea results\\june12\\fdr_all.csv")


ppath = "C:\\Users\\mzaidi\\Documents\\Visium Datasets\\SPATA2 Data\\c2c5 MILO only gsea results\\june12"

# STEP 3:----
# probably have to rename columns after converting to matrices
mydata <- lapply(myfiles, read.csv)
names(mydata) <- c("NES_all","pvalue_all","FDR_all")

rowNames <- substr(mydata[[1]]$Term,4,length(mydata[[1]]$Term))
colNames <- colnames(mydata[[1]])[2:ncol(mydata[[1]])]

for(file in 1:length(mydata)){
  mydata[[file]] <- mydata[[file]] %>% dplyr::select(-Term)
  rownames(mydata[[file]]) <- rowNames
  colnames(mydata[[file]]) <- colNames
}

NES_allClusters <- as.matrix(mydata[[1]])
pValue_allClusters <- as.matrix(mydata[[2]])
FDR_allClusters <- as.matrix(mydata[[3]])

pValue_allClusters <- pValue_allClusters[rowSums(is.na(NES_allClusters)) == 0, ]
FDR_allClusters <- FDR_allClusters[rowSums(is.na(NES_allClusters)) == 0, ]
NES_allClusters <- NES_allClusters[rowSums(is.na(NES_allClusters)) == 0, ]

ffile <- paste0(ppath, "/aMwwGSTs_allClusters.RData")
save(NES_allClusters, pValue_allClusters, FDR_allClusters, file = ffile)

NES <- NES_allClusters
# decide what parameters we want to use for this one
#NES[NES_allClusters < 1 & NES_allClusters > -1] <- 0
NES[NES_allClusters < 0.58] <- 0
NES[FDR_allClusters > 0.01] <- 0
#NES[pValue_allClusters > 0.01] <- 0
NES <- NES[which(rowSums(NES > 0) != 0), ]
NES <- NES[, which(colSums(NES > 0) != 0)]

#write.csv(NES, paste0(ppath, "/binarization_matrix_NESall.csv"), row.names=TRUE)
#write.csv(NES, paste0(ppath, "/binarization_matrix_NESpos_FDRlow.csv"), row.names=TRUE)
write.csv(NES, paste0(ppath, "/binarization_matrix_NES058_FDR001.csv"), row.names=TRUE)


jaccard <- matrix(0, nrow = ncol(NES), ncol = ncol(NES))
rownames(jaccard) <- colnames(jaccard) <- colnames(NES)
for(i in 1:ncol(NES))
  for(j in 1:ncol(NES)){
    a <- rownames(NES[NES[, i] != 0, i, drop = F])
    b <- rownames(NES[NES[, j] != 0, j, drop = F])
    jaccard[i, j] <- length(intersect(a, b))/length(union(a, b))
    #print(i)
  }
jaccard <- 1 - jaccard
#write.csv(NES, paste0(ppath, "/jaccard_NESall.csv"), row.names=TRUE)
#write.csv(NES, paste0(ppath, "/jaccard_NESpos_FDRlow.csv"), row.names=TRUE)
write.csv(NES, paste0(ppath, "/jaccard_NES058_FDR001.csv"), row.names=TRUE)


# Clustering ----
calinsky <- function(hhc, dist = NULL, gMax = round(1 + 3.3 * log(length(hhc$order), 10))) {
  msg <- ""
  if (is.null(dist)) {
    require(clue)
    dist <- sqrt(as.cl_ultrametric(hhc))
    # message(msg <- "The distance matrix not is provided, using the cophenetic matrix")
  } else if (attr(dist, "method") != "euclidean") {
    require(clue)
    dist <- sqrt(as.cl_ultrametric(hhc))
    # message(msg <- "The distance matrix is not euclidean, using the cophenetic matrix")
  }
  dist <- as.matrix(dist)^2
  A <- -dist/2
  A_bar <- apply(A, 1, mean)
  totalSum <- sum(diag(A) - 2 * A_bar + mean(A))
  n <- length(hhc$order)
  ans <- rep(0, gMax)
  for (g in 2:gMax) {
    cclust <- cutree(hhc, k = g)
    withinSum <- 0
    for (k in 1:g) {
      if (sum(cclust == k) == 1) 
        next
      A <- as.matrix(-dist/2)[cclust == k, cclust == k]
      A_bar <- apply(A, 1, mean)
      withinSum <- withinSum + sum(diag(A) - 2 * A_bar + 
                                     mean(A))
    }
    betweenSum <- totalSum - withinSum
    betweenSum <- betweenSum/(g - 1)
    withinSum <- withinSum/(n - g)
    ans[g] <- betweenSum/withinSum
  }
  class(ans) <- "calinsky"
  attr(ans, "message") <- msg
  return(ans)
}

sHc <- hclust(as.dist(jaccard), method = "average")
aCalinsky <- calinsky(sHc, gMax = 10)
nCls <- which.max(aCalinsky)

aConsensus <- ConsensusClusterPlus(d = as.dist(jaccard), maxK = ifelse(nCls == 2, nCls + 1, nCls), reps = 10000, pItem = 0.7,
                                   innerLinkage = "ward.D2", finalLinkage = "ward.D2", distance = "euclidean")
aConsensus <- aConsensus[[nCls]]

obj <- aConsensus
G <- nCls

consHc <- obj$consensusTree
consClust <- as.factor(cutree(consHc, k = G))
names(consClust) <- colnames(jaccard)
# clusterAssignments contains the lists of what sample-clusters belong to which hierarchical groups
clusterAssignments <- vector(mode = "list")
for(i in 1:G){
  clusterAssignments[[i]] <- names(consClust[consClust==i])
}
# save for use in future analysis (or other workflows)
#ffile <- paste0(ppath, "\\clusterAssignments_NESpos_FDRlow.RData")
ffile <- paste0(ppath, "\\clusterAssignments_NES058_FDR001.RData")
save(clusterAssignments, file = ffile)

library(tibble)
#to export cluster assignments as a df
max.len <- max(sapply(clusterAssignments,function(x) length(x)))

for(i in 1:G){
  clusterAssignments[[i]] = c(clusterAssignments[[i]], rep("null", max.len - length(clusterAssignments[[i]])))
}
# names(clusterAssignments) <- c("A", "B","C","D","E","F","G","H")
names(clusterAssignments) <- letters[1:G]

clusterAssignments <- as_tibble(clusterAssignments)
write.csv(clusterAssignments, paste0(ppath, "/clusterAssignments_df_NES058_FDR001.csv"))

thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
             "#bd18ea", #magenta
             "#2ef4ca", #aqua
             "#f4cced", #pink,
             "#f4cc03", #lightorange
             "#05188a", #navy,
             "#e5a25a", #light brown
             "#06f106", #bright green
             "#85848f", #med gray
             "#000000", #black
             "#076f25", #dark green
             "#93cd7f",#lime green
             "#4d0776", #dark purple
             "#ffffff" #white
)

#levels(consClust) <- palette("Set 3")[1:G]
levels(consClust) <- thisPal[1:G]
col <- colorRampPalette(c("lemonchiffon2", "#053061"))(51)
Colv  <- as.dendrogram(consHc)

colNames <- colnames(jaccard)
consMatrix <- obj$consensusMatrix
colnames(consMatrix) <- colnames(jaccard)
consMatrix <- consMatrix[, consHc$order]
RowSideColors <- as.character(consClust[colnames(consMatrix)])

newColors = c(thisPal[consClust])
colori=2
colorList=list(newColors,colori,unique(newColors) )

#library(heatmap3)
llabRow <- colnames(consMatrix)
ffile <- paste0(ppath, "/heatmap_NES058_FDR001.pdf")
pdf(file = ffile, width = 12, height = 8)
heatmap3(t(consMatrix), col = col, scale = "none",
         labCol = NA, labRow = llabRow, cexRow = 0.5,
         Colv = Colv, Rowv = NA,
         RowSideColors = RowSideColors)
dev.off()

#ffile <- paste0(ppath, "\\aConsensus_NESpos_FDRlow.RData")
ffile <- paste0(ppath, "\\aConsensus_NES058_FDR001.RData")
save(aConsensus, consHc, consClust, file = ffile)


#levels(consClust)[levels(consClust) == "black"] <- "cyan"

#library(heatmap3)
#ffile <- paste0(ppath, "/1a_NESpos_FDRlow.pdf")
ffile <- paste0(ppath, "/1a_NES058_FDR001.pdf")
pdf(file = ffile, width = 12, height = 8)
heatmap3(t(consMatrix),
         col = colorRampPalette(c("lemonchiffon2", "#053061"))(51),
         scale = "none",
         Colv = rev(as.dendrogram(consHc)),
         Rowv = NA,
         labCol = NA,
         labRow = llabRow,#colnames(consMatrix),
         cexRow = 0.7,
         RowSideColors = RowSideColors,
         RowSideLabs = NA)
legend("bottomleft",legend=unique(names(clusterAssignments)),fill=unique(RowSideColors),horiz=FALSE )
dev.off()


# Fig 1b----
# load("C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/primaryIDHwt/c2c5 gsea results/aMwwGSTs_allClusters.RData", verbose = T)
# load("C:/Users/mzaidi/Documents/Visium Datasets/Sheila_GBM/GBM_datasets/primaryIDHwt/c2c5 gsea results/aConsensus_NES058_FDR001.RData", verbose = T)

names <- letters[1:G]

levels(consClust) <- names

tmp <- as.character(consClust)
names(tmp) <- names(consClust)
consClust <- tmp

consClust <- consClust[consHc$order][length(consClust):1]
#consClust <- c(consClust[consClust == "A"], consClust[consClust == "B"])
#consClust <- c(consClust[consClust == "A"], consClust[consClust == "B"], 
               # consClust[consClust == "C"], consClust[consClust == "D"], 
               # consClust[consClust == "E"], consClust[consClust == "F"],
               # consClust[consClust == "G"], consClust[consClust == "H"],
               # consClust[consClust == "I"], consClust[consClust == "J"])
# consClust <- c(consClust[consClust == "A"], consClust[consClust == "B"], 
#                consClust[consClust == "C"])

consClust <- c(consClust[consClust == "a"], consClust[consClust == "b"],
               consClust[consClust == "c"], consClust[consClust == "d"],
               consClust[consClust == "e"], consClust[consClust == "f"],
               consClust[consClust == "g"])

#consClust <- c(consClust[consClust == "A"], consClust[consClust == "B"])

NES_allClusters <- NES_allClusters[, names(consClust)]

aDEA <- DEAgroups(ddata = NES_allClusters, groups = consClust, method = "MWW")
# edited these filtering params because nothing was passing through (got rid of qValue filter)
aDEA <- lapply(aDEA, function(x){
  x <- x[x$logFC > 0.3, ]
  # x <- x[x$qValue < 0.01, ]
  x <- x[order(x$pValue), ]
  x <- x[order(x$logFC, decreasing = T),]
  return(x)
})
aDEA <- lapply(aDEA, rownames)

for (cl in names(aDEA)){
  aPath <- as.character(aDEA[[cl]])
  toPlot <- NES_allClusters[aPath, ]
  classCol <- thisPal[1:G]
  names(classCol) <- names
  ColSideColors <- classCol[consClust]
  RowSideColors[aPath %in% aDEA[[cl]]] <- classCol[cl]
  toPlot <- (toPlot - rowMeans(toPlot))/apply(toPlot, 1, sd)
  toPlot[toPlot <= quantile(toPlot, 0.15)] <- quantile(toPlot, 0.15)
  toPlot[toPlot >= quantile(toPlot, 0.85)] <- quantile(toPlot, 0.85)
  
  ffile <- paste0(ppath, "/1b_", cl,".pdf")
  pdf(file = ffile,width = 12, height = 8)
  heatmap3(toPlot[nrow(toPlot):1, ], showRowDendro = F, showColDendro = F,
           Rowv = NA,
           Colv = NA,
           ColSideColors = ColSideColors, ColSideLabs = NA,
           RowSideColors = RowSideColors[nrow(toPlot):1], RowSideLabs = NA,
           col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75),
           #labCol = NA,
           labCol = colnames(toPlot),
           #labRow = NA,
           labRow = rownames(toPlot),
           cexCol = 0.5,#3,
           cexRow = 0.6,
           scale = 'none', useRaster = F)
  legend("bottomleft",legend=cl,fill=classCol[cl],horiz=FALSE )
  dev.off()
  
}

allPath <- as.character(unlist(aDEA))
dup <- sort(unique(allPath[duplicated(allPath)]))
allPath <- allPath[!allPath %in% dup]
aDEA <- lapply(aDEA, function(x) {x[x %in% allPath]})

toPlot <- NES_allClusters[allPath, ]

# classCol <- palette("Set 3")[1:G]
classCol <- thisPal[1:G]
names(classCol) <- names

ColSideColors <- classCol[consClust]
RowSideColors <- rep(c("#8DD3C7"), length(allPath))
for(n in 1:length(names)){
  RowSideColors[allPath %in% aDEA[[n]]] <- classCol[n]
}

toPlot <- (toPlot - rowMeans(toPlot))/apply(toPlot, 1, sd)
toPlot[toPlot <= quantile(toPlot, 0.15)] <- quantile(toPlot, 0.15)
toPlot[toPlot >= quantile(toPlot, 0.85)] <- quantile(toPlot, 0.85)

#library(heatmap3)
#library(gplots)
#ffile <- paste0(ppath, "/1b_NESpos_FDRlow.pdf")
ffile <- paste0(ppath, "/1b_NES058_FDR001.pdf")
pdf(file = ffile,width = 12, height = 8)
heatmap3(toPlot[nrow(toPlot):1, ], showRowDendro = F, showColDendro = F,
         Rowv = NA,
         Colv = NA,
         ColSideColors = ColSideColors, ColSideLabs = NA,
         RowSideColors = RowSideColors[nrow(toPlot):1], RowSideLabs = NA,
         col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75),
         #labCol = NA,
         labCol = colnames(toPlot),
         #labRow = NA,
         labRow = rownames(toPlot),
         cexCol = 0.5,#3,
         cexRow = 0.6,
         scale = 'none', useRaster = F)
legend("bottomleft",legend=unique(names(clusterAssignments)),fill=unique(ColSideColors),horiz=FALSE )
dev.off()


library (tibble)
#aDEA_df <- tibble("A"=aDEA[[1]],"B"=aDEA[[2]],"C"=aDEA[[3]])
aDEA_df = list(A=aDEA[[1]],B=aDEA[[2]],C=aDEA[[3]], D=aDEA[[4]], E=aDEA[[5]],F=aDEA[[6]],G=aDEA[[7]],H=aDEA[[8]])
max.len = max(length(aDEA[[1]]), length(aDEA[[2]]),length(aDEA[[3]]),
              length(aDEA[[4]]),length(aDEA[[5]]),length(aDEA[[7]]),length(aDEA[[6]]),length(aDEA[[8]]))
aDEA_df$A = c(aDEA_df$A, rep("0", max.len - length(aDEA_df$A)))
aDEA_df$B = c(aDEA_df$B, rep("0", max.len - length(aDEA_df$B)))
aDEA_df$C = c(aDEA_df$C, rep("0", max.len - length(aDEA_df$C)))
aDEA_df$D = c(aDEA_df$D, rep("0", max.len - length(aDEA_df$D)))
aDEA_df$E = c(aDEA_df$E, rep("0", max.len - length(aDEA_df$E)))
aDEA_df$F = c(aDEA_df$F, rep("0", max.len - length(aDEA_df$F)))
aDEA_df$G = c(aDEA_df$G, rep("0", max.len - length(aDEA_df$G)))
aDEA_df$H = c(aDEA_df$H, rep("0", max.len - length(aDEA_df$H)))


attributes(aDEA_df) = list(names = names(aDEA_df),
                      row.names=1:max(length(aDEA[[1]]), length(aDEA[[2]]),length(aDEA[[3]]),length(aDEA[[4]]),length(aDEA[[5]]),length(aDEA[[6]]),length(aDEA[[7]]),length(aDEA[[8]])), class='data.frame')


#write.csv(aDEA_df, paste0(ppath, "/aDEAdf_NESpos_FDRlow.csv"))
write.csv(aDEA_df, paste0(ppath, "/aDEAdf_NES058_FDR001.csv"))


#ffile <- paste0(ppath, "/1b_NESpos.png")
#ffile <- paste0(ppath, "/1b_NESall.png")

#ggplot2::ggsave(ffile)

#keywords----
# **************CLEAN THIS UP**************
keywords_list <- vector(mode = "list")
# keywords_count <- vector(mode = "list")

for(col in colnames(aDEA_df)){
  aDEA_df[[col]] <- chartr("_", " ", aDEA_df[[col]])
  t <- as.character(unlist(aDEA_df[[col]]))
  for(i in 1:length(t)){
    keywords_list[[col]] <- append(keywords_list[[col]],strsplit(t[i], " +")[[1]])
  }
  #keywords_list[[col]] <- sort(keywords_list[[col]])
  #keywords_count[[col]] <- sort(unique(keywords_list[[col]][duplicated(keywords_list[[col]])]))
}
max.len = max(length(keywords_list$A), length(keywords_list$B),length(keywords_list$C),length(keywords_list$D),length(keywords_list$E),length(keywords_list$F))
keywords_list$A = c(keywords_list$A, rep("0", max.len - length(keywords_list$A)))
keywords_list$B = c(keywords_list$B, rep("0", max.len - length(keywords_list$B)))
keywords_list$C = c(keywords_list$C, rep("0", max.len - length(keywords_list$C)))
keywords_list$D = c(keywords_list$D, rep("0", max.len - length(keywords_list$D)))
keywords_list$E = c(keywords_list$E, rep("0", max.len - length(keywords_list$E)))
keywords_list$F = c(keywords_list$F, rep("0", max.len - length(keywords_list$F)))
keywords_list$G = c(keywords_list$G, rep("0", max.len - length(keywords_list$G)))
keywords_list$H = c(keywords_list$H, rep("0", max.len - length(keywords_list$H)))



library(dplyr)
keywords_list <- as_tibble(keywords_list)

# for(i in 1:length(colnames(keywords_list))){
#   keywords_count %>% add_count(colnames(keywords_list)[i],sort = TRUE)
# }
keywordsA <- keywords_list %>% count(A,sort = TRUE)
keywordsB <- keywords_list %>% count(B,sort = TRUE)
keywordsC <- keywords_list %>% count(C,sort = TRUE)
keywordsD <- keywords_list %>% count(D,sort = TRUE)
keywordsE <- keywords_list %>% count(E,sort = TRUE)
keywordsF <- keywords_list %>% count(F,sort = TRUE)
keywordsG <- keywords_list %>% count(G,sort = TRUE)
keywordsH <- keywords_list %>% count(H,sort = TRUE)


#keywords <- c(keywordsA,keywordsB,keywordsC,keywordsD,keywordsE)

keywordsA <- rename(keywordsA, keyword = A) %>%
  add_column(group = rep("A",times=length(keywordsA[[1]])))
keywordsB <- rename(keywordsB, keyword = B) %>%
  add_column(group = rep("B",times=length(keywordsB[[1]])))
keywordsC <- rename(keywordsC, keyword = C) %>%
  add_column(group = rep("C",times=length(keywordsC[[1]])))
keywordsD <- rename(keywordsD, keyword = D) %>%
  add_column(group = rep("D",times=length(keywordsD[[1]])))
keywordsE <- rename(keywordsE, keyword = E) %>%
  add_column(group = rep("E",times=length(keywordsE[[1]])))
keywordsF <- rename(keywordsF, keyword = F) %>%
  add_column(group = rep("F",times=length(keywordsF[[1]])))
keywordsG <- rename(keywordsF, keyword = G) %>%
  add_column(group = rep("G",times=length(keywordsG[[1]])))
keywordsH <- rename(keywordsH, keyword = H) %>%
  add_column(group = rep("H",times=length(keywordsH[[1]])))

toRemove <- c("0","GOBP","HP","OF","GOMF","GOCC","AND","THE","IN","TO","BY",
              "VIA","I","1","1ST","2ND","A","II","2","3","6","D","8","KEGG",
              "REACTOME","WP", "PID")
keywords <- rbind(keywordsA,keywordsB,keywordsC,keywordsD,keywordsE,keywordsF,keywordsG,keywordsH,make.row.names = TRUE)
keywords <- keywords %>%
  filter(!keyword %in% toRemove) %>%
  # filter(length(keyword) > 1) %>%
  filter(n > 2)
#write.csv(keywords, paste0(ppath, "/keywords_NESpos_FDRlow.csv"))
write.csv(keywords, paste0(ppath, "/keywords_NES058_FDR001.csv"))


temp <- aDEA_df$B[grepl("PID", aDEA_df$B, fixed = TRUE)]
#temp <- aDEA_df$C[grepl("BETA", aDEA_df$C, fixed = TRUE)]
temp




