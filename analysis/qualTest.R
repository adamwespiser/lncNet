entropy <- function(...){
  x = c(...)
  x = x / sum(x)
  y = log2(x) * x
  y[is.na(y)]=0
  sum(y) * -1
}

entropyVec <- function(x,logFn=log2){
  x = x / sum(x)
  y = logFn(x) * x
  y[is.na(y)]=0
  sum(y) * -1
}

JS <- function(x,y,logFn=log2){
  x = x / sum(x)
  y = y / sum(y)
  a = (x + y)/2
  entropyVec(a ,logFn)- ((entropyVec(x,logFn)+entropyVec(y,logFn))/2) -> z
  #  entropyVec(a )- ((entropyVec(x)+entropyVec(y))/2) -> z
  z}

JSsp <- function(e1,e2,logFn=log2){
  1 - sqrt(JS(e1,e2,logFn))
}

calcTissSpec <- function(x,cols,logFn=log2){
  x <- sapply(x,as.numeric)
  profile <- diag(length(cols))
  specEn <- sapply(1:length(cols),function(y)JSsp(x,profile[y,],logFn))
  # cols[which(specEn == max(specEn))]
  max(specEn)
  #specEn
}

calcTissSpecVector <- function(x,logFn=log2){
  x <- sapply(x,as.numeric)
  profile <- diag(length(x))
  specEn <- sapply(seq_along(x),function(y)JSsp(x,profile[y,],logFn))
  # cols[which(specEn == max(specEn))]
  max(specEn)
  #specEn
}


# over individuals
testTissueConsistency <- function(ncols=5,distro=c("unif","pois")){
  nrows=500
  X <- NULL 
  if(identical(pmatch(distro,"uniform"), as.integer(1))){
    X=matrix(runif(nrows * ncols),ncol=ncols)
  }
  else if(identical(pmatch(distro,"poisson"), as.integer(1))){
    X=matrix(rpois((nrows * ncols),lambda=0.5) + runif((nrows * ncols),min=0,max=0.1),ncol=ncols)
  }else {
    stop("data distribution not defined")
  }
  maxExprVec <- apply(X,1,function(x)which(x==max(x)))
  majorityCellType <- as.numeric(names(table(maxExprVec))[which(table(maxExprVec) == max(table(maxExprVec)))])
  JS.compare <- rep(0,ncols)
  JS.compare[majorityCellType] <- 1
  tissueConsistency <- apply(X, 1,function(x)JSsp(x,JS.compare,log10))
  mean(tissueConsistency)
}


variationOfTissConsistency <- function(){
  reps=1:10
  range=4:50
  results <- data.frame(trial=seq_len(length(reps) * length(range) ))
  results$rep <- NA
  results$range <- NA
  results$meanConsistency <- NA
  results$localConsistencyPois <- NA
  for(i in reps){
    for(j in range){
      row = (i - 1) * length(range) + (j - min(range) + 1 )
      print(paste(i,j,row))
      localConsistency <- testTissueConsistency(ncols=j,distro="unif")
      localConsistencyPois <- testTissueConsistency(ncols=j,distro="pois")
      results$rep[row] = i
      results$range[row] = j
      results$meanConsistency[row] = mean(localConsistency)  
      results$localConsistencyPois[row] = mean(localConsistencyPois)  
      
    }
  }
  
  stats <- data.frame(
  mean=tapply(results$meanConsistency,results$range,mean),
  sd=tapply(results$meanConsistency,results$range,sd),
  meanP=tapply(results$localConsistencyPois,results$range,mean),
  sdP=tapply(results$localConsistencyPois,results$range,sd),
  range=as.numeric(names(tapply(results$meanConsistency,results$range,sd))))
  
  
  ggplot(stats, aes(x=range,y=mean))+geom_line()+
    xlab("number of entries")+ylab("mean tissue consistency")+
    ggtitle("tissue consistency mean by vector length\nuniform samplin")
  
  
  ggplot(stats, aes(x=range,y=meanP))+geom_line()+
    xlab("number of entries")+ylab("mean tissue consistency")+
    ggtitle("tissue consistency mean by vector length\nPois + uniform noise sampling")
  
}





