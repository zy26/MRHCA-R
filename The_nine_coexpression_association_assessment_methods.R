load("Ecoli_RNAseq_data.RData")

#Pearson Correlation
cor_c <- cor(t(Ecoli_RNASeq_top1000_randomized))
#Spearman Correlation
cor_c <- cor(t(Ecoli_RNASeq_top1000_randomized), method = "spearman")
#Kendall Correlation
cor_c <- cor(t(Ecoli_RNASeq_top1000_randomized), method = "kendall")

#Mutual Information
library(entropy)
cor_c <- matrix(0, nrow(Ecoli_RNASeq_top1000_randomized), nrow(Ecoli_RNASeq_top1000_randomized))
rownames(cor_c) <- rownames(Ecoli_RNASeq_top1000_randomized)
colnames(cor_c) <- rownames(Ecoli_RNASeq_top1000_randomized)
for (i in 1:nrow(Ecoli_RNASeq_top1000_randomized)) {
    for (j in 1:nrow(Ecoli_RNASeq_top1000_randomized)) {
        if (i < j) {
            disc <- discretize2d(Ecoli_RNASeq_top1000_randomized[i,], Ecoli_RNASeq_top1000_randomized[j,], numBins1 = 5, numBins2 = 5)
            cor_c[i, j] <- cor_c[j, i] <- mi.empirical(disc)
        }
    }
}

#Hoeffding's measure of dependence
library(Hmisc)
cor_c <- hoeffd(t(Ecoli_RNASeq_top1000_randomized))

library(MASS)
library(class
library(cluster)
library(impute)
library(Hmisc)
library(WGCNA)

#Polynomial-Regression based dependence
polyReg_dist <- adjacency.polyReg(t(Ecoli_RNASeq_top1000_randomized), degree = 3, symmetrizationMethod = "mean")
#Spline-Regression based dependence
splineReg_dist <- adjacency.splineReg(t(Ecoli_RNASeq_top1000_randomized), symmetrizationMethod = "mean")

#Distance Covariance
library(energy)
cor_c <- matrix(0, nrow(Ecoli_RNASeq_top1000_randomized), nrow(Ecoli_RNASeq_top1000_randomized))
rownames(cor_c) <- rownames(Ecoli_RNASeq_top1000_randomized)
colnames(cor_c) <- rownames(Ecoli_RNASeq_top1000_randomized)
for (i in 1:nrow(Ecoli_RNASeq_top1000_randomized)) {
    for (j in 1:nrow(Ecoli_RNASeq_top1000_randomized)) {
        if (i < j) {
            cor_c[j, i] <- cor_c[i, j] <- dcov(Ecoli_RNASeq_top1000_randomized[i,], Ecoli_RNASeq_top1000_randomized[j,])
        }
    }
}

#Weighted Rank Correlation
WRC <- function(x, y) {
    RR <- rank(x)
    QQ <- rank(y)
    N <- length(x)
    cc <- 1 - 6 * sum((RR - QQ) ^ 2 * (N - RR + 1 + N - QQ + 1)) / (N ^ 4 + N ^ 3 - N ^ 2 - N)
    return(cc)
}
cor_c <- matrix(0, nrow(Ecoli_RNASeq_top1000_randomized), nrow(Ecoli_RNASeq_top1000_randomized))
rownames(cor_c) <- rownames(Ecoli_RNASeq_top1000_randomized)
colnames(cor_c) <- rownames(Ecoli_RNASeq_top1000_randomized)
for (i in 1:nrow(Ecoli_RNASeq_top1000_randomized)) {
    for (j in 1:nrow(Ecoli_RNASeq_top1000_randomized)) {
        if (i < j) {
            cor_c[j, i] <- cor_c[i, j] <- WRC(Ecoli_RNASeq_top1000_randomized[i,], Ecoli_RNASeq_top1000_randomized[j,])
        }
    }
}
