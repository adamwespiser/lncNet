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
require(network)
# make toy random network
x <- 100
ndyads <- x * (x - 1)
density <- x / ndyads
nw.mat <- matrix(0, nrow = x, ncol = x)
dimnames(nw.mat) <- list(1:x, 1:x)
nw.mat[row(nw.mat) != col(nw.mat)] <- runif(ndyads) < density
nw.mat[row(nw.mat) != col(nw.mat)] <- runif(ndyads)
nw.mat
rnd <- network(mcl(nw.mat))
rnd

category = LETTERS[rbinom(x, 4, .5)]
ggnet(rnd, label = TRUE, color = "white", segment.color = "grey10", node.group = category)


x = 100
cutoff <- 3.763007
mat <- matrix(rexp(n=x*x,rate=0.5),nrow=x)
dimnames(mat) <- list(1:x, 1:x)
rnd <- network(mat)
mclNet <- network(mcl(applyThresh(X=mat,cutoff)))
category = LETTERS[rbinom(x, 4, .5)]

ggnet(rnd, label = TRUE, color = "white", segment.color = "grey10", node.group = category)
ggnet( network(mcl(applyThresh(X=mat,5))), label = TRUE, color = "white", segment.color = "grey10", node.group = category)

ggplot(melt(mat), aes(x=Var1,y=Var2,fill=value))+geom_tile()
ggplot(melt(applyThresh(X=mat,5)), aes(x=Var1,y=Var2,fill=value))+geom_tile()
ggplot(melt(mcl(applyThresh(X=mat,5))), aes(x=Var1,y=Var2,fill=value))+geom_tile()


n = 50
block1 <- matrix(runif(n=n*n),nrow=n)
block2 <-matrix( runif(n=n*n)*0.1,nrow=n)
block3<- matrix(runif(n=n*n)*0.1,nrow=n)
block4<- matrix(runif(n=n*n),nrow=n)

blocks <- cbind(rbind(block2,block1),rbind(block4,block3))
blocks <- blocks + t(blocks)
dimnames(blocks) <- list(1:(2*n), 1:(2*n))
blockClus <- mcl(applyThresh(blocks,0))

ggplot(melt(blocks), aes(x=Var1,y=Var2,fill=value))+geom_tile()
ggplot(melt(blockClus), aes(x=Var1,y=Var2,fill=value))+geom_tile()
ggnet( network(blockClus), label = TRUE, 
       color = "white", segment.color = "grey10",
       node.group =  factor(apply(blockClus,2,which.max),labels=c('a','b')))
ggnet( network(blocks), label = TRUE, 
       color = "white", segment.color = "grey10",
       node.group =  factor(apply(blockClus,2,which.max),labels=c('a','b')))

#MCL 
# square rows, normalize
# set matrix as square of self
# calculate chnage


getLabels <- function(X){
  apply(X,1,which(X > 0))
}

t = matrix(runif(25),nrow=5)

applyThresh <- function(X,cutoff){
  X[which(X < cutoff)] <- 0
  minVal <- min(X[X>0])
  X <- X - minVal
  X <- X/max(X)
  X[X < 0] <- 0
  X
}


mcl <- function(X){
  
  for(i in 1:50){
    X1 <- X %*% X
    
    X <- apply(X1,2,function(x)x^2/sum(x^2))

  change <- sum(abs(X - X1))
    print(change)
  if(change < 2e-4){
    return(X)
  }
  }
  X
}



degreeDistro <- function(X){
  dis <- unlist(apply(X,1,function(x)sum(x > 0)))
  data.frame(mean=mean(dis),median=median(dis),singletons=sum(dis == 1))
}

getCutoff <- function(X){
  for(i in seq(min(X),max(X),length.out=20)){
    X.local <- X
    X.local[X.local < i] <- 0
    print(i)
    print(degreeDistro(X.local))  
  }
  
}






