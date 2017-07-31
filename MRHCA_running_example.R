load("Ecoli_RNAseq_data.RData")

#Step 1. Compute C matrix: e.g. Pearson Correlation
cor_C0 <- cor(t(Ecoli_RNASeq_top1000))
cor_C <- cor(t(Ecoli_RNASeq_top1000_randomized))

#Step 2. Simulated Empirical Null Distribution of MR, i.e. MR_E
#Parameters:
#Rounds: Number of runs in generation of MR_E
#Inputs:
#Gene expression data, e.g. Ecoli_RNASeq_top1000_randomized
MR_E <- Empirical_null_distribution_MR(Ecoli_RNASeq_top1000_randomized, Rounds = 1000)

#Step 3. Identify Hub Genes and Estimate size of Co-expression Modules
#Parameters:
#step_size0=50, Step size in computation of growth rate
#p_sig=0.05,  p value cut-off
#Hit_score_cutoff=100,  Hits Score Cutoff
#Inputs:
#C matrix, e.g. cor_C
#Empirical Null Distribution of MR, e.g. MR_E
MRHCA_output <- MRHCA(cor_C0, MR_E, step_size0 = 50, p_sig = 0.05, Hit_score_cutoff = 100)
MRHCA_output <- MRHCA(cor_C, MR_E, step_size0 = 50, p_sig = 0.05, Hit_score_cutoff = 100)

MRHCA_output[[1]] #Hub and Module size information. First column: 1 for hub, 0 for non-hub; Second column£º Estimated Module Size.
MRHCA_output[[2]] #MR matrix

#Step 4. Identify Hub Genes and Estimate size of Co-expression Modules for Small hubs
#Recommended after Step 3 if the stepsize in computing MR_growth_rate is larger than 20
MRHCA_small_hub <- MRHCA_small_hub_identificaiton(MRHCA_output, data = Ecoli_RNASeq_top1000_randomized, MR_E)

MRHCA_small_hub #Hub and Module size information of small hubs. First column: 1 for small hub; Second column£º Estimated Module Size. 
