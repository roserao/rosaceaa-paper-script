compute_fdr <- function(truth, test, bool = FALSE) {
  if (!bool) {
    truth <- (truth != "Neutral")
    test <- (test != "Neutral")
  }
  
  if (length(truth) != length(test)) {
    stop("length different")
  }
  
  l <- length(truth)
  TP <- sum(truth & test, na.rm = TRUE)/l
  FN <- sum(truth & (!test), na.rm = TRUE)/l
  FP <- sum((!truth) & test, na.rm = TRUE)/l
  TN <- sum((!truth) & (!test), na.rm = TRUE)/l
  
  if (FP + TP > 0) {
    FDR <-  FP / (FP + TP)
  } else {
    FDR <- NA
  }
  
  if (TP + FN > 0) {
    POWER <- TP / (TP + FN)
  } else {
    POWER <- NA
  }
  
  return(c(TP, FN, FP, TN, FDR, POWER))
}

compute_rankfdr <- function(true_label, test_stats, test_value, resolution) {
  df_rank <- data.frame(
    true_label = true_label, 
    test_stats = test_stats,
    test_value = abs(test_value))
  df_rank$true_label <- (df_rank$true_label != "Neutral")
  df_rank <- df_rank %>% dplyr::arrange(test_stats, desc(test_value))
  
  rank_fdr <- data.frame(rank = (1:resolution)/resolution,
                         TP = -1, FN = -1, FP = -1, TN = -1,
                         FDR = -1, Power = -1, Sensitivity = -1)
  cutoff <- floor(nrow(df_rank)/resolution * (1:resolution))
  for (i in 1:resolution) {
    rank_fdr[i, 2:7] <- compute_fdr(truth = df_rank$true_label[1:cutoff[i]], 
                                    test = rep(TRUE, cutoff[i]), bool = TRUE)
    rank_fdr[i, 8] <- sum(df_rank$true_label[1:cutoff[i]])/sum(df_rank$true_label)
  }
  return(tibble(rank_fdr))
}


compute_fdrsense <- function(true_label, test_stats) {
  
  df <- data.frame(
    true_label = true_label, 
    test_stats = test_stats)
  df$true_label <- (df$true_label != "Neutral")
  
  fdr_seq <- c(seq(0, 0.1, by = 0.001), seq(0.11, 0.5, by = 0.01))
  
  fdr_sense <- data.frame(fdr = fdr_seq,
                          TP = -1, FN = -1, FP = -1, TN = -1,
                          FDR = -1, Power = -1)
  for (i in 1:nrow(fdr_sense)) {
    fdr <- fdr_sense[i, 1]
    fdr_sense[i, 2:7] <- compute_fdr(df$true_label, (df$test_stats <= fdr),
                                     bool = TRUE)
  }
  
  return(tibble(fdr_sense))
}
