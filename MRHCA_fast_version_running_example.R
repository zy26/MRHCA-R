load("Ecoli_RNAseq_data.RData")
data <- Ecoli_RNASeq_top1000
data_randomized <- Ecoli_RNASeq_top1000_randomized

#Step 1. Compute Rank Index matrix with certain K for certain association method, e.g. Pearson Correlation. K controls computational consumption O(N^2*K) and memory consumption O(N*K)
K <- min(500, nrow(data))
Rank_index <- c()
for (i in 1:nrow(data)) {
    cor_cc <- cor(t(data), data[i,])
    cor_cc[i] <- 0
    Rank_index <- rbind(Rank_index, order(-abs(cor_cc))[1:K])
}
rownames(Rank_index) <- rownames(data)
#Or use the following function for compute rank_index for Pearson Correlation
Rank_index <- Compute_rank_index(data, K = K)

#Step 2. Simulated Empirical Null Distribution of MR, i.e. MR_E
#Parameters:
#Rounds: Number of runs in generation of MR_E
#K
#Inputs:
#Gene expression data, e.g. Ecoli_RNASeq_top1000_randomized
MR_E <- Empirical_null_distribution_MR_F(Rank_index, K = K, Rounds = 1000)

#Step 3. Identify Hub Genes and Estimate size of Co-expression Modules
#Parameters:
#step_size0=50, Step size in computation of growth rate
#p_sig=0.05,  p value cut-off
#Hit_score_cutoff=100,  Hits Score Cutoff
#K
#Inputs:
#Rank_index matrix, e.g. Rank_index
#Empirical Null Distribution of MR, e.g. MR_E
MRHCA_output <- MRHCA_fast_version(Rank_index, MR_E, K = K)

MRHCA_output[[1]]
#Hub and Module size information. First column: 1 for hub, 0 for non-hub; Second column, non-zero value for estimated Module Size (<K); Third column, lower bound of estimated Module Size (>=K)
MRHCA_output[[2]]
#MR matrix

#Step 4. Identify Hub Genes and Estimate size of Co-expression Modules for Small hubs
#Recommended after Step 3 if the stepsize in computing MR_growth_rate is larger than 20
MRHCA_small_hub <- MRHCA_small_hub_identificaiton(MRHCA_output, data = data_randomized, MR_E)
MRHCA_small_hub #Hub and Module size information of small hubs. First column: 1 for small hub; Second column?? Estimated Module Size. 
