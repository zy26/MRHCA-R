MRHCA <- function(cor_M, MR_E, step_size0 = 50, p_sig = 0.05, Hit_score_cutoff = 100) {
    cor_c <- cor_M
    diag(cor_c) <- 0
    cor_c <- abs(cor_c)
    cor_r_rank1 <- apply(cor_c, 1, rank)
    cor_r_rank2 <- t(apply(cor_c, 2, rank))
    cor_r_rank1 <- nrow(cor_r_rank1) - cor_r_rank1 + 1
    cor_r_rank2 <- nrow(cor_r_rank2) - cor_r_rank2 + 1
    MR_M <- cor_r_rank1 * cor_r_rank2
    ddd <- MR_E
    null_growth_rate <- c()
    step_size <- step_size0 #parameter
    print("Generating Empirical Null Growth Rate")
    for (i in 1:nrow(ddd)) {
        x0 <- sqrt(sort(ddd[i,]))
        null_growth_rate <- rbind(null_growth_rate, calculate_growth_rate(x0, step = step_size))
    }
    print("Hub Identification")
    hub_info_all <- c()
    sig_N <- ncol(null_growth_rate) #parameter
    p_sig <- p_sig0 #parameter
    hub_sig_cut <- Hit_score_cutoff #parameter
    for (i in 1:nrow(MR_M)) {
        x0 <- sqrt(sort(MR_M[i,]))
        tg_growth_rate <- calculate_growth_rate(x0, step = 50)
        is_sig <- growth_rate_significance(tg_growth_rate, null_growth_rate, sig_N) < p_sig
        S_all <- growth_rate_scoring(is_sig)
        cluster_bound <- 0
        if_hub <- 0
        if ((max(S_all) > hub_sig_cut) & (sum(tg_growth_rate[1:50] > 0.7) == 0)) {
            cluster_bound <- max(which(S_all == max(S_all)))
            if_hub <- 1
        }
        hub_info <- c(if_hub, cluster_bound)
        hub_info_all <- rbind(hub_info_all, hub_info)
    }
    rownames(hub_info_all) <- rownames(MR_M)
    return(list(hub_info_all, MR_M))
}



Empirical_null_distribution_MR <- function(data, Rounds = 1000) {
    N <- nrow(data)
    aaa <- runif(N * Rounds, 0, 1)
    bbb1 <- rep(0, length(aaa))
    bbb2 <- rep(0, length(aaa))
    for (i in 1:length(aaa)) {
        bbb1[i] <- rbinom(1, N - 1, aaa[i])
        bbb2[i] <- rbinom(1, N - 1, aaa[i])
    }
    ccc <- (bbb1 + 1) * (bbb2 + 1)
    ddd <- matrix(ccc, Rounds, N)
    return(ddd)
}



growth_rate_scoring <- function(is_sig, up_score = 1, down_score = -2, delete_top = 10) {
    S <- 0
    S_all <- c()
    for (i in 1:length(is_sig)) {
        if (i < delete_top) {
            S <- S + 0
            S_all <- c(S_all, S)
        }
        else {
            if (is_sig[i] == 1) {
                S <- S + up_score
                S_all <- c(S_all, S)
            }
            if (is_sig[i] == 0) {
                S <- S + down_score
                S_all <- c(S_all, S)
            }
        }
    }
    return(S_all)
}


growth_rate_significance <- function(tg_growth_rate, null_growth_rate, sig_N) {
    growth_rate_sig <- rep(0, sig_N)
    for (i in 1:sig_N) {
        growth_rate_sig[i] <- sum(tg_growth_rate[i] >= null_growth_rate[, i])
    }
    growth_rate_sig <- growth_rate_sig / nrow(null_growth_rate)
    return(growth_rate_sig)
}


calculate_growth_rate <- function(x0, step = 20) {
    gr_all <- rep(0, length(x0))
    for (i in 1:length(x0)) {
        sid <- max(1, i - step)
        eid <- min(length(x0), i + step)
        tg_id <- c(sid:eid)
        l <- eid - sid + 1
        gr_all[i] <- (x0[eid] - x0[sid]) / l
    }
    return(gr_all)
}


calculate_MR <- function(tg_id, Rank_index) {
    N <- nrow(Rank_index)
    K <- ncol(Rank_index)
    MR_result <- rep(0, N)
    for (i in 1:N) {
        R1 <- K + 1
        R2 <- K + 1
        if (sum(tg_id == Rank_index[i, ]) > 0) {
            R1 <- which(Rank_index[i,] == tg_id)
        }
        if (sum(i == Rank_index[tg_id, ]) > 0) {
            R2 <- which(Rank_index[tg_id,] == i)
        }
        MR_result[i] <- R1 * R2
        #if(i%%1000==1)
        #{
        #print(i)
        #}
    }
    MR_result[tg_id] <- K ^ 2
    return(MR_result)
}


MRHCA_small_hub_identificaiton <- function(MRHCA_output, data, MR_E) {
    hubs_small_result <- MRHCA_output
    data_c <- data
    MR_M <- hubs_small_result[[2]]
    MR_M0 <- MR_M
    for (i in 1:nrow(MR_M0)) {
        MR_M0[i,] <- sort(MR_M0[i,])
    }
    ############
    MR_E0 <- MR_E
    for (i in 1:nrow(MR_E0)) {
        MR_E0[i,] <- sort(MR_E0[i,])
    }
    #############
    ccc <- hubs_small_result[[1]]
    tg_cut_id <- 100
    tg_hubs <- c()
    if (sum(ccc[, 1] == 1) > 0) {
        tg_cut_id <- max(min(ccc[which(ccc[, 1] == 1), 2]), 100)
        tg_hubs <- names(which(hubs_small_result[[1]][, 1] == 1))
    }
    print("Cut Smallest Module Size!!!\n")
    print(tg_cut_id)
    MR_M_c <- sqrt(MR_M0[, 1:tg_cut_id])
    MR_E_c <- sqrt(MR_E0[, 1:tg_cut_id])
    p_cut <- 0.01
    stat_MR <- MR_M_c * 0
    for (i in 1:ncol(stat_MR)) {
        stat_MR[, i] <- MR_M_c[, i] < quantile(MR_E_c[, i], p_cut)
    }
    cumsum_MR <- stat_MR * 0
    for (i in 1:nrow(stat_MR)) {
        cumsum_MR[i,] <- cumsum(stat_MR[i,])
    }

    for (i in 1:nrow(cumsum_MR)) {
        if (length(which(cumsum_MR[i, ] != 0)) == 0) {
            fff <- rep(1, length(cumsum_MR[i,]))
        }
        else {
            tg_ids <- min(which(cumsum_MR[i,] != 0))
            fff <- c(rep(1, tg_ids - 1), 1:(length(cumsum_MR[i,]) - tg_ids + 1))
        }
        cumsum_MR[i,] <- cumsum_MR[i,] / fff
    }

    #intersect(names(which(apply(cumsum_MR>0.8,1,sum)>5)),tg_hubs)
    tg_id_ccc <- setdiff(names(which(apply(cumsum_MR > 0.8, 1, sum) > 2)), tg_hubs)
    small_hub_info <- c()
    small_hubs <- c()
    small_hub_info_selected <- c()
    if (length(tg_id_ccc) > 0) {
        for (i in 1:length(tg_id_ccc)) {
            x0 <- sqrt(sort(MR_M[tg_id_ccc[i],]))
            tg_growth_rate <- calculate_growth_rate(x0, step = 10)
            ccc <- cumsum(tg_growth_rate < 0.5) - cumsum(tg_growth_rate > 1) * 100 - cumsum(tg_growth_rate > 0.85) * 5
            ddd <- ccc / c(1:length(ccc))
            if ((sum(ccc[1:10] < 0) == 0) & (max(ccc) > 0.8)) {
                if (sum(ddd < 0) > 0) {
                    small_hub_info <- rbind(small_hub_info, c(1, min(which(ddd < 0))))
                }
                else {
                    small_hub_info <- rbind(small_hub_info, c(1, ncol(MR_M)))
                }
                small_hubs <- c(small_hubs, tg_id_ccc[i])
            }
        }
        rownames(small_hub_info) <- small_hubs
        small_hub_info_selected <- small_hub_info[which(small_hub_info[, 2] <= 300),]
        if (length(which(small_hub_info[, 2] <= 300)) == 1) {
            tg_id <- rownames(small_hub_info)[which(small_hub_info[, 2] <= 300)]
            small_hub_info_selected <- t(as.matrix(small_hub_info_selected))
            rownames(small_hub_info_selected) <- tg_id
        }
        print("small hub selected: ")
        print(length(which(small_hub_info[, 2] <= 300)))
        print("small hub removed: ")
        print(length(which(small_hub_info[, 2] > 300)))
        #hist(small_hub_info[,2],breaks=100)
    }
    return(small_hub_info_selected)
}


########################
#Fast version functions
########################


Empirical_null_distribution_MR_F <- function(data, K = 1000, Rounds = 1000) {
    N <- nrow(data)
    aaa <- runif(N * Rounds, 0, 1)
    bbb1 <- rep(0, length(aaa))
    bbb2 <- rep(0, length(aaa))
    for (i in 1:length(aaa)) {
        bbb1[i] <- rbinom(1, N - 1, aaa[i])
        bbb2[i] <- rbinom(1, N - 1, aaa[i])
    }
    bbb1[which(bbb1 > (M))] <- M
    bbb2[which(bbb2 > (M))] <- M
    ccc <- (bbb1 + 1) * (bbb2 + 1)
    ddd <- matrix(ccc, Rounds, N)
    return(ddd)
}


Compute_rank_index <- function(data, K) {
    Rank_index <- c()
    for (i in 1:nrow(data)) {
        cor_cc <- cor(t(data), data[i,])
        cor_cc[i] <- 0
        Rank_index <- rbind(Rank_index, order(-abs(cor_cc))[1:K])
    }
    rownames(Rank_index) <- rownames(data)
    return(Rank_index)
}

MR_E <- Empirical_null_distribution_MR_F(Rank_index, K, Rounds = 1000)

MRHCA_fast_version <- function(Rank_index, MR_E, K = 500, step_size0 = 50, p_sig0 = 0.05, Hit_score_cutoff = 100, TN_p = 1) {
    sig_N0 <- floor(K * 0.75)
    hub_sig_cut0 <- Hit_score_cutoff
    MR_id <- 1:nrow(Rank_index)
    if (TN < 1) {
        MR_id <- names(sort(table(Rank_index), decreasing = T))
        MR_id <- as.numeric(MR_id)
    }
    MR_EM <- c()
    TN <- floor(nrow(Rank_index) * TN_p)
    rnames <- c()
    print("Compute MR")
    for (i in 1:TN) {
        tg_id <- MR_id[i]
        a1 <- calculate_MR(tg_id, Rank_index)
        rnames <- c(rnames, tg_id)
        MR_EM <- rbind(MR_EM, a1)
    }
    rownames(MR_EM) <- rnames
    colnames(MR_EM) <- c(1:length(a1))
    #print("Compute Empirical Null Distribution")
    ####calculate null distribution of MR
    ddd2 <- MR_E
    ####calculate null distribution of MR growth rate
    null_growth_rate <- c()
    step_size <- step_size0 #parameter
    for (i in 1:nrow(ddd2)) {
        x0 <- sqrt(sort(ddd2[i,]))
        null_growth_rate <- rbind(null_growth_rate, calculate_growth_rate(x0, step = step_size))
    }
    print("Hub Identification")
    hub_info_all <- c()
    sig_N <- sig_N0 #parameter
    p_sig <- p_sig0 #parameter
    hub_sig_cut <- hub_sig_cut0 #parameter
    for (i in 1:nrow(MR_EM)) {
        x0 <- sqrt(sort(MR_EM[i,]))
        #x0<-sqrt(sort(calculate_MR(i,Rank_index)))
        tg_growth_rate <- calculate_growth_rate(x0, step = 50)
        is_sig <- growth_rate_significance(tg_growth_rate, null_growth_rate, sig_N) < p_sig
        S_all <- growth_rate_scoring(is_sig)
        cluster_bound <- 0
        if_hub <- 0
        if ((max(S_all) > hub_sig_cut) & (sum(tg_growth_rate[1:50] > 0.7) == 0)) {
            cluster_bound <- max(which(S_all == max(S_all)))
            if_hub <- 1
        }
        hub_info <- c(if_hub, cluster_bound, 0)
        if (cluster_bound == sig_N) {
            hub_info[3] <- max(which(tg_growth_rate != 0))
        }
        hub_info_all <- rbind(hub_info_all, hub_info)
    }
    hub_info_all[which(hub_info_all[, 2] == sig_N), 2] <- 0
    rownames(hub_info_all) <- rownames(MR_EM)
    return(list(hub_info_all, MR_EM))
}
