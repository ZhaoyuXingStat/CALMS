 
ALMS <- function(Y,phi,lam, 
                etaini,fixed.adaptive=F,
                admmAbsTol=1e-2,admmRelTol=1e-2){
  t1 <- proc.time() 
  x <- as.numeric(multi_signal_lasso(Y,phi,lam, 
                                   etaini=etaini,
                                   admmAbsTol=1e-2,admmRelTol=1e-2))
  t2 <- proc.time() 
  wei <- apply(cbind(abs(x),abs(x-1)),1,min)
  one.step.x <- x
  
  if(sum(wei)==0){ 
    return(list(one.step.x,
                one.step.x,
                t2-t1,
                t2-t1))}  
  
  
  ind <- which(wei==0)   
  indi <- which(wei!=0) 
  phi2 <- t(t(phi[,indi]) * wei[indi])  
  if(length(indi)==1){
    phi2 <- matrix(phi2   , ncol = length(indi))}  
  
  if(length(ind)==1){
    Y2 <- Y- matrix(phi[,ind],ncol=1) %*% matrix(x[ind],ncol=1) 
  }
  
  if(length(ind)!=1){
    Y2 <- Y- as.matrix(phi[,ind]) %*% matrix(x[ind],ncol=1)  
  }
  
  
  
  t3 <- proc.time()
  x2 <- as.numeric(multi_signal_lasso(Y2,phi2,lam,
                                      etaini=as.numeric(qr.solve(Y2,phi2)),
                                      admmAbsTol=1e-2,admmRelTol=1e-2))
  t4 <- proc.time()
  final <- c()
  final[indi] <- x2 * wei[indi]
  final[ind] <- one.step.x[ind]
  
  
  
  return(list(as.numeric(one.step.x),
              as.numeric(final),
              t2-t1,
              t2-t1+t4-t3))
}


MAL_cv <- function(Y,
                   phi,
                   lam.list=seq(0.01,10,length = 50),
                   etaini 
){
  n = nrow(phi)
  
  require(glmnet)
  group <- c(1:2, sample(c(1,2),n-2,replace = T))  
  
  pre <- 0
  for(i in 1:2){
    pre = pre + 0 %in% colSums(phi[which(group==i),])
  }
  
  while(pre > 0 ){
    group <- c(1:2, sample(c(1,2),n-2,replace = T))  
    pre <- 0
    for(i in 1:2){
      pre = pre + 0 %in% colSums(phi[which(group==i),])
    }
  }
  
  mse <- c()
  for(lam in 1:length(lam.list)){
    mse.tem <- c()
    for(g in 1:2){
      trainY <- Y[-which(group==g)]
      testY  <- Y[which(group==g)]
      trainPhi <- phi[-which(group==g),]
      testPhi <- phi[which(group==g),] 
      xhat <- multi_signal_lasso(trainY,trainPhi,lam=lam.list[lam],etaini=etaini) 
      mse.tem[g] <- mean((testY-testPhi%*%xhat)^2) 
    }
    mse[lam] <- mean(mse.tem)
  }
  (lambda_final <- min(which(mse==min(mse))))   
  xhat_final <- MAL(Y,phi,lam=lam.list[lambda_final],etaini=etaini)
  return(xhat_final)
}



