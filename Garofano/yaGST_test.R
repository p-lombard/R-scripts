library(devtools)
install_github("miccec/yaGST")

nr <- 100; nc <- 1000
# generate a data-matrix with nr samples, and nc features
exprData <- matrix(rpois(nc * nr, 100), nrow = nr, ncol = nc)
colnames(exprData) <- paste0("feat", 1:nc)
rownames(exprData) <- paste0("sam", 1:nr)

# increase the first 3 samples (minority set) of 10\% of the original intensity
# of the first 30 features (later the gene-set)
exprData[1, 1:30] <- exprData[1, 1:30]* runif(30, min = 1, max = 1.10)
exprData[2, 1:30] <- exprData[1, 1:30]* runif(30, min = 1, max = 1.10)
exprData[3, 1:30] <- exprData[1, 1:30]* runif(30, min = 1, max = 1.10)
samples_of_interest <- rownames(exprData)[1:3] # minority set

# running in parallel
library(doParallel)
# adjust the number of CPUs as needed
cl <- makePSOCKcluster(3)
clusterApply(cl, floor(runif(length(cl), max = 10000000)), set.seed)
registerDoParallel(cl)
ans_eeMWW <- yaGST::eeMWW(exprData, samples_of_interest)
stopCluster(cl)

# set the gene-set and run the enrichment analysis
geneSet <- colnames(exprData)[1:30]
(tmp <- mwwGST(ans_eeMWW, geneSet))
plot(tmp, rankedList = ans_eeMWW)


means <- rowMeans(H3K27)
sds <- apply(H3K27, 1, sd)

for (col in ncol(H3K27)){
  currentSample <- (H3K27[, col] - means)/sds
  rankedList <- sort(currentSample, decreasing = T)
  print(rankedList)
}

currsamp = (H3K27[, 1] - means)/sds
rl = sort(currsamp, decreasing = T)

geneSet <- yaGST::gmt2GO("C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/c5.all.v2023.1.Hs.symbols.gmt")
aMwwGST <- lapply(geneSet, function(x) yaGST::mwwGST(rl, geneSet = x, minLenGeneSet = 15, alternative = "two.sided", verbose = F))
aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0]
tmp_NES <- sapply(aMwwGST, function(x) x$log.pu)
tmp_pValue <- sapply(aMwwGST, function(x) x$p.value)
ans <- list(tmp_NES = tmp_NES, tmp_pValue = tmp_pValue)

originalGeneSetCount <- sapply(aMwwGST, function(x) x$originalGeneSetCount)
actualGeneSetCount <- sapply(aMwwGST, function(x) x$actualGeneSetCount)
NES <- sapply(aMwwGST, function(x) x$nes)
odd_NES <- sapply(aMwwGST, function(x) x$pu)
logit2_NES <- sapply(aMwwGST, function(x) x$log.pu)
pValue <- sapply(aMwwGST, function(x) x$p.value)
qValue <- p.adjust(pValue, method = "BH")

aMwwGST <- cbind(originalGeneSetCount, actualGeneSetCount, NES, odd_NES, logit2_NES, pValue, qValue)
aMwwGST <- as.data.frame(aMwwGST)


  