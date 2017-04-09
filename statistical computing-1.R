###---------------------------Code of a statistical computing assignment--------------------###

#--------------------- Problem 1 -------------------------
# calculate gradient and Hessian matrix of posterior log-likelihood
grad.hessian <- function(x, y, beta) {
  a <- exp(beta[1] + beta[2] * x)
  pai <- a / (1 + a)
  ## gradient
  grad <- c(- beta[1] + sum(y - pai), - beta[2] + sum((y - pai) * x))
  ## Hessian matrix
  H <- matrix(0, 2, 2)
  H[1,1] <- - 1 - sum(a / (1 + a)^2)
  H[1,2] <- - sum((a * x)/(1 + a)^2)
  H[2,1] <- H[1,2]
  H[2,2] <- - 1 - sum((a * x * x)/(1 + a)^2)
  ## check if hessian matrix is singular
  if(rcond(H) < .Machine$double.eps)         
  {                                              
    stop("Hessian matrix is singular!")
  }
  ## posterior log-likelihood
  l.star <- - log(2 * pi) - (beta[1]^2 + beta[2]^2)/2 + sum(y * log(pai) + (1 - y) * log(1 - pai))
  
  return(list(pai = pai, grad = grad, hessian = H, l.star = l.star))
}

# calculate the Laplace approximation of marginal log-likelihood
Laplace.approx <- function(x, y, tol = 1e-4) {
  ## mode of posterior distribution by Newton-Raphson
  beta.current <- rep(0, 2)
  gh <- grad.hessian(x, y, beta.current)
  beta.new <- beta.current - solve(gh$hessian) %*% (gh$grad)
  iter <- 0
  while(abs(beta.new - beta.current)[1] > tol | abs(beta.new - beta.current)[2] > tol) {
    iter <- iter + 1
    beta.current <- beta.new
    gh <- grad.hessian(x, y, beta.current)
    beta.new <- beta.current - solve(gh$hessian) %*% (gh$grad)
  }
  
  beta <- beta.new
  mode.gh <- grad.hessian(x, y, beta)
  pai <- mode.gh$pai
  H <- mode.gh$hessian
  
  ## calculate marginal log-likelihood
  logmarglik <- - (beta[1]^2 + beta[2]^2)/2 + sum(y*log(pai) + (1 - y) * log(1 - pai)) - log(det(H)) / 2
  
  return(list(iter = iter, beta = beta, hessian = H, logmarglik = logmarglik))
}



#----------------------- Problem 2 --------------------------
MH <- function(x, y, k) {
  laplace <- Laplace.approx(x, y)
  beta.current <- laplace$beta
  beta.current <- as.numeric(beta.current)   ## turn beta into a row vector
  sigma <- - solve(laplace$hessian)
  stepNumber <- 0
  MHbeta <- NULL
  
  while(stepNumber <= k) {
    stepNumber = stepNumber + 1  
    mu <- beta.current
    library(MASS)
    beta.candi <- mvrnorm(n = 1, mu = mu, Sigma = sigma)
    
    if(grad.hessian(x, y, beta.candi)$l.star >= grad.hessian(x, y, beta.current)$l.star) {
      beta.current <- beta.candi
    }
    else{
      u <- runif(1)
      if(log(u) <= (grad.hessian(x, y, beta.candi)$l.star - grad.hessian(x, y, beta.current)$l.star)) {
        beta.current <- beta.candi
      }
    }
    MHbeta <- rbind(MHbeta, beta.current)
  }
  
  beta.est <- apply(MHbeta, 2, mean)
  
  return(beta.est = beta.est)
  
}



#------------------------- Problem 3 -----------------------------
bayesLogistic <- function(apredictor, response, data, NumberOfIterations) {
  x <- data[, apredictor]
  y <- data[, response]
  
  laplace <- Laplace.approx(x, y)
  logmarglik <- laplace$logmarglik
  beta.mle <- as.numeric(coef(glm(y ~ x, family = binomial(link = "logit"))))
  beta0mle <- beta.mle[1]
  beta1mle <- beta.mle[2]
  
  beta.est <- MH(x, y, NumberOfIterations)
  beta0bayes <- beta.est[1]
  beta1bayes <- beta.est[2]
  
  return(list(apredictor = apredictor, logmarglik = logmarglik, beta0bayes = beta0bayes, 
              beta1bayes = beta1bayes, beta0mle = beta0mle, beta1mle = beta1mle))
}

# PARALLEL VERSION
# datafile = the name of the file with the data
# NumberOfIterations = number of iterations of the Metropolis-Hastings algorithm
# clusterSize = number of separate processes; each process performs one or more
# univariate regressions
main <- function(datafile, NumberOfIterations, clusterSize) {
  # read the data
  data <-  read.table(datafile, header = F)
  
  # the sample size is 148 (number of rows)
  # the explanatory variables are the first 60 columns for '534binarydata.txt'
  # the last column is the binary response
  response <- ncol(data)
  lastPredictor <- ncol(data) - 1
  
  # initialize a cluster for parallel computing
  cluster <- makeCluster(clusterSize, type = "SOCK")
  
  clusterExport(cluster, list("grad.hessian", "Laplace.approx", "MH", "bayesLogistic",
                              "mvrnorm"))
  
  # run the MC3 algorithm from several times
  results <- clusterApply(cluster, 1:lastPredictor, bayesLogistic,
                          response, data, NumberOfIterations)
  
  # print out the results
  for(i in 1:lastPredictor) {
    cat('Regression of Y on explanatory variable ', results[[i]]$apredictor,
        ' has log marginal likelihood ', results[[i]]$logmarglik,
        ' with beta0 = ', results[[i]]$beta0bayes, ' (',results[[i]]$beta0mle,')',
        ' and beta1 = ', results[[i]]$beta1bayes, ' (',results[[i]]$beta1mle,')',
        '\n')   
  }
  
  # destroy the cluster
  stopCluster(cluster)
}

# THE PACKAGE 'SNOW' IS NEEDED FOR PARALLEL COMPUTING
require(snow)
require(MASS)

# this is where the program starts
main('534binarydata.txt',10000,10)


