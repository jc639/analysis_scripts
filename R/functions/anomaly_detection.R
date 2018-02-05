# function to determine p(x) in multivariate gaussian distribution
# parameterised by N~mu;Sigma

normalParams <- function(x=NULL){
  mu <- matrix(apply(x, 2, FUN = mean))
  sigma <- cov(x)
  
  return(list(mu=mu, sigma=sigma))
}

gaussianConfInt <- function(x = NULL, mu=NULL, Sigma=NULL, alpha=NULL){
  if(any(c(is.null(x), is.null(mu), is.null(sigma), is.null(alpha)))){
    stop("pass all parameters: one or more missing
         x = named list for each condition
         mu = column vector of means (length n) for each type
         sigma = covariance matrix
         alpha = value between 0 and 1")
  }
  
  truth_df <- as.data.frame(matrix(nrow=length(names(x)), ncol=2))
  colnames(truth_df) <- c("condition", "signif")
  truth_df$condition <- names(x)
  
  for(name in names(x)){
    values <- matrix(x[[name]])
    outside_confidence <- as.vector(t(values - mu) %*% solve(Sigma) %*% (values - mu)) > qchisq(1-alpha, df = dim(values)[1])
    truth_df$signif[truth_df$condition==name] <- outside_confidence
  }
  
  return(truth_df)
}