# It takes as input:
#   - Y: a response vector
#   - phi: a design matrix
#   - lam: a tuning parameter controlling the degree of sparsity in the solution
#   - etaini: an initial value for the optimization procedure
#   - admmAbsTol: absolute tolerance for the ADMM algorithm
#   - admmRelTol: relative tolerance for the ADMM algorithm
#   - Gamma: a parameter controlling the weight assigned to non-zero coefficients in the penalty term
# The output is a list containing:
#   - one.step.x: the solution obtained after one step of the CMAL algorithm
#   - x: the final solution obtained after the CMAL algorithm
#   - t2-t1: the time taken for the first step of the algorithm (multi_signal_lasso)
#   - t2-t1+t4-t3: the total time taken for the algorithm

# Some Remarks:
#   multi_signal_lasso is called to perform the first step of the CMAL algorithm. Its inputs are Y, phi, lam, etaini, admmAbsTol and admmRelTol. The output is x, which is stored in the variable x.
#   The weights for the penalty term are computed using apply and cbind functions. The resulting vector wei contains the minimum values of the abs(x) and abs(x-1) columns, where x is the solution obtained from the first step of the CMAL algorithm.
#   If all elements of wei are zero, then the solution is already binary and no further processing is needed. The function returns a list with one.step.x and x equal to the same binary solution.
#   If there are non-zero elements in wei, then the indices ind and indi are computed based on the values of wei.
#   Next, the design matrix phi2 is computed by weighting the columns of phi corresponding to non-zero elements in wei.
#   Depending on whether ind contains a single or multiple indices, the response vector Y2 is computed using matrix operations.


# The function constrained_MDSL is called to perform the second step of the CMAL algorithm. Its inputs are:
#   - Y2: the modified response vector after removing the effect of the zero-weighted design matrix columns
#   - phi2: the modified design matrix
#   - lam: the tuning parameter controlling the degree of sparsity in the solution
#   - rho: a parameter for the ADMM algorith
#   - CMALweights: the weights for the penalty term
#   - etaini: an initial value for the optimization procedure
#   - admmAbsTol: absolute tolerance for the ADMM algorithm
#   - admmRelTol: relative tolerance for the ADMM algorithm
# Finally, the resulting solution x2 is combined with the first step solution one.step.x to obtain the final solution. The output is a list containing one.step.x, x, and two time measurements t2-t1 and t2-t1+t4-t3.



CMAL <- function(Y,phi,lam, # non-refit
                 etaini=NULL,rho=10,
                 admmAbsTol=0.1,admmRelTol=0.1,
                 Gamma_tuning=1){

  t1 <- proc.time() #timer
  
  one.step.x <- as.numeric(multi_signal_lasso(Y,phi,lam, 
                                     etaini=etaini,
                                     admmAbsTol=1,admmRelTol=1))
  t2 <- proc.time()
  
  #print(paste("time：",t[3][[1]],'second'))
  wei <- as.numeric(apply(cbind(abs(one.step.x),abs(one.step.x-1)),1,min)) 
  one.step.x  
  
  # if one step over, refit and return
  iiind <- 0
  if(sum(wei)==0){#print("completed compressing");
    print("ALL 01")
    iiind <- 1
    return(list(one.step.x,
                one.step.x,
                t2-t1,
                t2-t1,
                 iiind))}  

  if(sum(wei)!=0){
    print("【NOT】———————————— 01")

  ind <- which(wei==0)  
  indi <- which(wei!=0) #   ind for Fixed coef; indi for non-fixed coef     
  
  phi2 <- t( t(phi[,indi]) * wei[indi]^Gamma_tuning )  # Fix the coefficients has already been estimated as 0 or 1
  
  if(length(indi)==1){
    # phi2 <- matrix( phi2, ncol = 1 )
    # r <- Y- as.matrix(phi[,ind]) %*% matrix(x[ind],ncol=1) # residuals
    return(list(one.step.x,
                one.step.x,
                t2-t1,
                t2-t1,
                iiind))
    }  
  
  if(length(ind)==1){
    Y2 <- Y- matrix(phi[,ind],ncol=1) %*% matrix(one.step.x[ind],ncol=1) 
  }
  
  if(length(ind)!=1){
    Y2 <- Y- as.matrix(phi[,ind]) %*% matrix(one.step.x[ind],ncol=1)  
  }
  
  t3 <- proc.time()
  x2 <- constrained_MDSL(as.numeric(Y2),phi2,
                         lam,
                         rho = rho, 
                         CMALweights = 1 / (wei[indi]^Gamma_tuning ) ,        
                         etaini =  etaini[indi] , # 
                         admmAbsTol=1e-2,
                         admmRelTol=1e-2)
  t4 <- proc.time()
  final <- c()
  final[indi] <- x2   *   wei[indi]^Gamma_tuning   
  final[ind] <- one.step.x[ind]
  # final[which(final!=1 & final !=0)] <- one.step.x[which(final!=1 & final !=0)]
  
  final
  
 

  return(list(as.numeric(one.step.x),
              as.numeric(final),
              t2-t1,
              t4-t3,
              iiind))}
}









CMAL_cv <- function(Y,
                    phi,
                    lam.list=seq(0.01,10,length = 50),
                    etaini,rho,x,ct="cv",c2=1 # ct for crateria 
){
  require(sBIC)
  n = 5

  require(glmnet)
  group <- c(1:n)  
 
  
  mse <- c()
  for(lam in 1:length(lam.list)){
    mse.tem <- c()
    for(g in 1:n){
      trainY <- Y[-which(group==g)]
      testY  <- Y[which(group==g)]
      trainPhi <- phi[-which(group==g),]
      testPhi <- phi[which(group==g),] 
      xhat <- constrained_MDSL(trainY,trainPhi,
                               lam=lam.list[lam],
                               rho=rho,
                               etaini=etaini)
      #mse.tem[g] <- mean((x-xhat)^2)
      if(ct=="cv"){mse.tem[g] <- mean((testY-testPhi%*%xhat)^2)}
      if(ct=="bic"){
        # modified BIC cerition
        sbic = log(sum((testY-testPhi%*%xhat)^2)) + c2* log(length(testY))*(length(xhat)-sum(xhat==0|xhat ==1))
        mse.tem[g] <- sbic}
      #print(mse.tem)
    }
    mse[lam] <- mean(na.omit(mse.tem))
  }
  (lambda_final <- min(which(mse==min(mse))))
  #print(lam.list[lambda_final])
  xhat_final <- CMAL(Y,phi,
                     lam=lam.list[lambda_final],
                     etaini=etaini)
  return(xhat_final)
}




