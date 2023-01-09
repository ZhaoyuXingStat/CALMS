 
CMAL <- function(Y,phi,lam, # non-refit
                 etaini,
                 admmAbsTol=1e-2,admmRelTol=1e-2,
                 Gamma=1){
  
  t1 <- proc.time() #timer
  
  x <- as.numeric(multi_signal_lasso(Y,phi,lam, 
                                     etaini=etaini,
                                     admmAbsTol=1e-2,admmRelTol=1e-2))
  t2 <- proc.time()
  
  wei <- as.numeric(apply(cbind(abs(x),abs(x-1)),1,min))  
  one.step.x <- x
  
  if(sum(wei)==0){
    print("ALL 01")
    return(list(one.step.x,
                one.step.x,
                t2-t1,
                t2-t1))}  

  if(sum(wei)!=0){
  ind <- which(wei==0)   
  indi <- which(wei!=0)     
  
  phi2 <- t( t(phi[,indi]) * wei[indi]^Gamma )   
  
  if(length(indi)==1){
 
    return(list(one.step.x,
                one.step.x,
                t2-t1,
                t2-t1))
    }  
  
  if(length(ind)==1){
    Y2 <- Y- matrix(phi[,ind],ncol=1) %*% matrix(x[ind],ncol=1) 
  }
  
  if(length(ind)!=1){
    Y2 <- Y- as.matrix(phi[,ind]) %*% matrix(x[ind],ncol=1)  
  }
  
  t3 <- proc.time()
  x2 <- constrained_MDSL(Y2,phi2,
                         lam,
                         rho=10,  
                         CMALweights = 1 / (wei[indi]^Gamma ) ,        
                         etaini =  etaini[indi] , 
                         admmAbsTol=1e-2,
                         admmRelTol=1e-2)
  t4 <- proc.time()
  final <- c()
  final[indi] <- x2   
  final[ind] <- one.step.x[ind]
  # final[which(final!=1 & final !=0)] <- one.step.x[which(final!=1 & final !=0)]
  
  x <- final
  
 

  return(list(as.numeric(one.step.x),
              as.numeric(x),
              t2-t1,
              t2-t1+t4-t3))}
}









CMAL_cv <- function(Y,
                    phi,
                    lam.list=seq(0.01,10,length = 50),
                    etaini,rho,x,ct="cv",c2=1 # ct for crateria 
){
  require(sBIC)
  n = nrow(phi)

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
      
      if(ct=="cv"){mse.tem[g] <- mean((testY-testPhi%*%xhat)^2)}
      if(ct=="bic"){
        
        sbic = log(sum((testY-testPhi%*%xhat)^2)) + c2* log(length(testY))*(length(xhat)-sum(xhat==0|xhat ==1))
        mse.tem[g] <- sbic}
      
    }
    mse[lam] <- mean(na.omit(mse.tem))
  }
  (lambda_final <- min(which(mse==min(mse))))
  
  xhat_final <- CMAL(Y,phi,
                     lam=lam.list[lambda_final],
                     etaini=etaini)
  return(xhat_final)
}



