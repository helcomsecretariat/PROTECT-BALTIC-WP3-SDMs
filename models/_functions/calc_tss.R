calc_tss <- function(obs, pred, na.rm = TRUE) {

  # Checks
  if (!all(obs %in% c(0, 1, NA)))  stop("Observed column must be binary 0/1 (or NA).")
  if (!all(pred %in% c(0, 1, NA))) stop("Predicted column must be binary 0/1 (or NA).")
  
  # Handle NAs
  if (na.rm) {
    keep <- !is.na(obs) & !is.na(pred)
    obs  <- obs[keep]
    pred <- pred[keep]
  } else if (anyNA(obs) || anyNA(pred)) {
    stop("NAs present. Set na.rm = TRUE to drop NA pairs.")
  }
  
  # Confusion matrix components
  TP <- sum(obs == 1 & pred == 1)
  TN <- sum(obs == 0 & pred == 0)
  FP <- sum(obs == 0 & pred == 1)
  FN <- sum(obs == 1 & pred == 0)
  
  # Rates
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_  # sensitivity / TPR
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_  # specificity / TNR
  
  # TSS = sensitivity + specificity - 1
  tss <- if (is.na(sens) || is.na(spec)) NA_real_ else sens + spec - 1
  
  data.frame(
    TSS = tss,
    sensitivity = sens,
    specificity = spec,
    TP = TP, TN = TN, FP = FP, FN = FN,
    n = length(obs)
  )
}