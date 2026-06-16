# x = binary presence absence data 0s and 1s only
# train_perc = proportion of data for training, e.g. 0.8
# n = number of times to repeat the generation of cross-validation set

cv_partition <- function(x, train_perc, n){
  
  if(any(is.na(x))){stop("NAs in x!")}
  ux <- unique(x)
  if (!all(ux %in% c(0, 1))) stop("x must contain only 0 and 1.")
  
  a <- which(x==0)
  p <- which(x==1)
  
  prev <- sum(x)/length(x)
  
  adj_p_train <- round(train_perc*length(p))
  adj_a_train <- round(train_perc*length(a))
  
  for(i in 1:n){
    p_train <- sample(p, adj_p_train, replace = FALSE)
    a_train <- sample(a, adj_a_train, replace = FALSE)
    train <- c(p_train, a_train)
    
    cv <- rep(FALSE, length(x))
    cv[train] <- TRUE
    
    df <- data.frame("col" = cv)
    colnames(df) <- paste0("_allData_RUN", i)
    
    if(i == 1){dff <- df}
    if(i > 1){dff <- cbind(dff, df)}
  }
  
  return(dff)
  
}