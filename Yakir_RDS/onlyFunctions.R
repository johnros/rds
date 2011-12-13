#----------------- Utility functions---------------#

# Old versions #
#load('/home/johnros/Projects/Yakir/oldFunctions.Rdata')

####################################################




import <- function (file) {
  raw1<- scan(file, what=character())
  return(as.numeric(strsplit(raw1, split=",")[[1]]))
}


sum.ranks<- function(data){
  data.table<- table(data)[-1]
  as.numeric(rownames(data.table))*c(data.table)
}



comparison <- function (estimation) {
  non.na<- !is.na(estimation)
  print(signif(cbind(sapply(estimation[non.na], function(x) return(x$Nj)))))
  print(sapply(estimation[non.na], function(x) sum(x$Nj)))
  print(sapply(estimation[non.na], function(x) x$theta))
  print(sapply(estimation[non.na], function(x) x[[4]]$value))
}


# Maps the line to the range: minimum+range
inv.qnorm.theta<- function(qnorm.theta, const, range=20, minimum=-10){
  normalized.theta<- pnorm(qnorm.theta / const)
  normalized.theta2<- (normalized.theta*range) + minimum
  return(normalized.theta2)
}
## Test:
# inv.qnorm.theta(10, 100)
# inv.qnorm.theta(10, 10)
#inv.qnorm.theta(20000,10)
# inv.qnorm.theta(25,200)
# inv.qnorm.theta(25,200, range=2)
# inv.qnorm.theta(25,200, range=1)
# curve(inv.qnorm.theta(x,const=0.6, range=2, minimum=1), -10,10)

qnorm.theta<- function(theta, const, range=20, minimum=-10){
  normalized.theta<- (theta-minimum)/range
  const*qnorm(normalized.theta)
}
## Test:
# qnorm.theta(25, range=30, minimum=-1, const=0.1)
# qnorm.theta(theta=-0.5, const=10, range=2, minimum=-1)
# inv.qnorm.theta(qnorm.theta(5,100),100)
# curve(qnorm.theta(x,const=0.6, range=1, minimum=0), -2,2)
# curve(qnorm.theta(x,const=30, range=1, minimum=-1), -2,2)






#---------------------- CPP likelihod version-------------#
 
 # Optimization stage:

## TODO: Fix estimation for a single degree!
dyn.load("/home/johnros/Projects/Yakir/likelihoodVer9.32.so")
dyn.load("/home/johnros/Projects/Yakir/likelihoodVer9.64.so")



likelihood.cpp.wrap.general <- function (func, data, init, initial.Nj, const, arc, maxit=10000, theta.range, theta.minimum,...) {
  # Look for degrees in the data, so their estimates are nony vanishing
  N.j<- rep(0, max(data)) 
  uniques<- unique(data)
  uniques<- uniques[uniques!=0]
  s.uniques<- sort(uniques)
  param.size<- length(uniques)
  data.table<-table(data)[-1]
  
  likelihood.wrap1<- function(par){
	  final.result<- -Inf
	  ccc<- exp(par[1])
	  theta<- inv.qnorm.theta(par[2], const = const, range=theta.range, minimum=theta.minimum)
	  N.j<- rep(0, max(data))
	  N.j[s.uniques]<- exp(tail(par,-2))
	  if(any( N.j[s.uniques] < data.table )) return(final.result)	  
	  
	  if(is.numeric(ccc) && is.numeric(theta) && !is.infinite(theta) && !is.infinite(ccc) ) {
		  result<- .C(func, 
				  sample=as.integer(data), 
				  c=as.numeric(ccc), 
				  theta=as.numeric(theta), 
				  Nj=as.numeric(N.j), 
				  constant=as.numeric(const),
				  n=as.integer(length(data)), 
				  N=as.integer(length(N.j)),
				  arc=arc,
				  result=as.double(0))
		  final.result<- result$result
	  }	  
	  return(final.result)
  }
  
  cap<- max(data)
  cap.arcs<-  sum.ranks(data)
  max.cap.arcs<- max(cap.arcs)
  log.c<- -log(max(sum.ranks(data)*sum(data)) * max(data)) -1
  
  if(missing(init)) { 
	  init<- list( 
			  theta0.2= c(log.c=log.c, qnorm.theta=qnorm.theta(0.2,const, range=theta.range, minimum=theta.minimum), log(data.table)+1 ),
        theta0.5= c(log.c=log.c, qnorm.theta=qnorm.theta(0.5,const, range=theta.range, minimum=theta.minimum), log(data.table)+1 ),
			  theta1= c(log.c=log.c, qnorm.theta=qnorm.theta(1,const,range=theta.range, minimum=theta.minimum), log(data.table)+1 ),
        theta2= c(log.c=log.c, qnorm.theta=qnorm.theta(2,const, range=theta.range, minimum=theta.minimum), log(data.table)+1 ),
			  theta3= c(log.c=log.c, qnorm.theta=qnorm.theta(3,const, range=theta.range, minimum=theta.minimum), log(data.table)+1 ),
			  theta4= c(log.c=log.c, qnorm.theta=qnorm.theta(4,const, range=theta.range, minimum=theta.minimum), log(data.table)+1 ),
			  theta10= c(log.c=log.c, qnorm.theta=qnorm.theta(10,const, range=theta.range, minimum=theta.minimum), log(data.table)+1 )			  			  
	  )}
  
  #range.ind<- sapply(init, function(x) x['qnorm.theta']>theta.minimum && x['qnorm.theta.']<theta.minimum+range)
  likelihood.optim<-lapply(init, function(x) {
			  try(optim(par=x, fn=likelihood.wrap1, control=list(fnscale=-1, maxit=maxit),...)) }  )  
  
  prepare.result<- function(x){    
    result<- NA
    
    if(is.list(x)){
      N.j[s.uniques]<- exp(tail(x$par,-2))
      result<- list( 
        c=exp(x$par[[1]]), 
        theta=inv.qnorm.theta(x$par[[2]],const, range=theta.range, minimum=theta.minimum), 
        Nj=N.j, 
        x  )
    }    
    return(result)
  }   
  return(lapply(likelihood.optim, function(x) try(prepare.result(x))))
}
    
#dyn.load("/specific/a/home/cc/students/math/rosenbla/Projects/Yakir/likelihoodVer8.so")
dyn.load("../likelihoodVer8.so")


# ## Test:
#  (tempData<- sample(c(t(data1[1,])), size=1000))
#  tempData<- data1[1,]
#  (temp21<-likelihood.cpp.wrap.general(func="likelihood", data=tempData, const=100, arc=TRUE, theta.minimum=-1, theta.range=4, maxit=10000))
#  comparison(temp21)




likelihood.cpp2.wrap.general<- function (func, data, init, initial.Nj, const, arc, maxit=50000,...) {
  # Look for degrees in the data, so their estimates are nony vanishing
  N.j<- rep(0, max(data)) 
  uniques<- unique(data)
  uniques<- uniques[uniques!=0]
  s.uniques<- sort(uniques)
  param.size<- length(uniques)
  data.table<-table(data)[-1]
  
  likelihood.wrap1<- function(par){
    final.result<- -Inf
	  ccc<- exp(par[1])
	  theta<- inv.qnorm.theta(par[2], const = const)
	  N.j<- rep(0, max(data))
	  N.j[s.uniques]<- exp(tail(par,-2))
	  if(any( N.j[s.uniques] < data.table )) return(final.result)	  
	  
	  if(is.numeric(ccc) && is.numeric(theta) && !is.infinite(theta) && !is.infinite(ccc) ) {
		  result<- .C(func, 
				  sample=as.integer(data), 
				  c=as.numeric(ccc), 
				  theta=as.numeric(theta), 
				  Nj=as.numeric(N.j), 
				  constant=as.numeric(const),
				  n=as.integer(length(data)), 
				  N=as.integer(length(N.j)),
				  arc=arc,
				  result=as.double(0))
		  final.result<- result$result
	  }	  
	  #print(c(as.integer(data),log.c=as.numeric(par[1]),qnorm.theta=as.numeric(par[2]),Nj=as.numeric(N.j),constant=as.numeric(const)))
	  
	  return(final.result)
  }
  cap<- max(data)
  cap.arcs<-  sum.ranks(data)
  max.cap.arcs<- max(cap.arcs)
  log.c<- -log(max(sum.ranks(data)*sum(data)) * max(data)) -1
  
  if(missing(init)) { 
	  init<- list( 
			  #one= c(log.c=log.c, qnorm.theta=qnorm.theta(1,const), log(rep(10*cap, param.size))+1),
			  #two= c(log.c=log.c, qnorm.theta=qnorm.theta(1,const), log(rep(initial.Nj,param.size))+1),
			  #three= c(log.c=log.c, qnorm.theta=qnorm.theta(1,const), log(runif(min=cap, max=3*cap ,n=param.size))+1) ,
			  #four= c(log.c=log.c, qnorm.theta=qnorm.theta(1,const), log(cap.arcs)+1),
			  #five= c(log.c=log.c, qnorm.theta=qnorm.theta(1,const), log(rep(max.cap.arcs, param.size))+1 ),
			  six0.2= c(log.c=log.c, qnorm.theta=qnorm.theta(0.2,const), log(data.table)+1 ),
			  #six0.5= c(log.c=log.c, qnorm.theta=qnorm.theta(0.5,const), log(data.table)+1 ),
			  #six0.9= c(log.c=log.c, qnorm.theta=qnorm.theta(0.9,const), log(data.table)+1 ),
			  six1= c(log.c=log.c, qnorm.theta=qnorm.theta(1,const), log(data.table)+1 ),
			  #six1.1= c(log.c=log.c, qnorm.theta=qnorm.theta(1.1,const), log(data.table)+1 ),
			  #six1.5= c(log.c=log.c, qnorm.theta=qnorm.theta(1.5,const), log(data.table)+1 ),
			  six2= c(log.c=log.c, qnorm.theta=qnorm.theta(2,const), log(data.table)+1 )			  
	  )}
  
  likelihood.optim<-lapply(init, function(x) {
			  try(optim(par=x, fn=likelihood.wrap1, control=list(fnscale=-1, maxit=maxit),...)) }  )  
  
  prepare.result<- function(x){    
    N.j[s.uniques]<- exp(tail(x$par,-2))
    result<- list( 
      c=exp(x$par[[1]]), 
      theta=inv.qnorm.theta(x$par[[2]],const), 
      Nj=N.j, 
	  x  )
    return(result)
  }  
  return(lapply(likelihood.optim, function(x) try(prepare.result(x))))
}





likelihood.cpp2.wrap.fix.theta<- function (func, data, Sij, init, initial.Nj, const, arc, maxit=1000, theta, ...) {
  # Look for degrees in the data, so their estimates are nony vanishing
  N.j<- rep(0, max(data)) 
  uniques<- unique(data)
  uniques<- uniques[uniques!=0]
  s.uniques<- sort(uniques)
  param.size<- length(uniques)
  data.table<-table(data)[-1]
  
  S<- compute.S(Sij)
  
  likelihood.wrap1<- function(par){
    final.result<- -Inf
	  ccc<- exp(par[1])
	  N.j<- rep(0, max(data))
	  N.j[s.uniques]<- exp(tail(par,-1))
	  if(any( N.j[s.uniques] < data.table )) return(final.result)	  
	  
	  if(is.numeric(ccc) && is.numeric(theta) && !is.infinite(theta) && !is.infinite(ccc) ) {
		  result<- .C(func, 
				  sample=as.integer(data), 
				  Sij=as.integer(as.matrix(Sij)),
				  S=as.integer(S),				  
				  c=as.numeric(ccc), 
				  theta=as.numeric(theta), 
				  Nj=as.numeric(N.j), 
				  constant=as.numeric(const),
				  observed_degrees=as.integer(rownames(Sij)),
				  n=as.integer(length(data)), 
				  N=as.integer(length(N.j)),
				  N_observed=as.integer(nrow(Sij)),				  
				  arc=arc,
				  result=as.double(0))		  
		  final.result<- result$result
	  }	  	  
	  return(final.result)
  }
  cap<- max(data)
  cap.arcs<-  sum.ranks(data)
  max.cap.arcs<- max(cap.arcs)
  log.c<- -log(max(sum.ranks(data)*sum(data)) * max(data)) -1
  
  if(missing(init)) { 
	  init<- list( 
			  six= c(log.c=log.c, log(data.table)+1 )

	  )}
  
   likelihood.optim<-lapply(init, function(x) {
			  try(optim(par=x, fn=likelihood.wrap1, control=list(fnscale=-1, maxit=maxit),...)) 
		  }  ) 
  
  prepare.result<- function(x){    
    N.j[s.uniques]<- exp(tail(x$par,-1))
    result<- list( 
      c=exp(x$par[[1]]), 
      theta=theta,
      Nj=N.j, 
	  x  )
    return(result)
  }  
  return(lapply(likelihood.optim, function(x) try(prepare.result(x))))
}






















######################################################



exctract.best <- function (y) {
  # Extractor of resutls from Condor in **wierd** format.
  location.index<- which(!sapply(y, is.null)) # this is needed due to the stupid output of condor.R
  max.ind<-which.max(sapply(y[[location.index]], function(x) x[[4]]$value))
  return(y[[location.index]][max.ind])
}


exctract.best2 <- function (y) {
  # Extractor of resutls from Condor in **reasonbale** format.
  max.ind<-which.max(sapply(y, function(x) x[[4]]$value))
  return(y[[max.ind]])
}


exctract.best3 <- function (y) {
   # this is needed due to the stupid output of condor.R
    location.index<- which(!sapply(y, is.null))
    #browser()
    max.ind<-which.max( as.numeric(sapply(y[[location.index]], function(x){
        if(length(x)==4) return(x[[4]]$value)
        else return(NA)}) ))                              
    return(y[[location.index]][max.ind])
}
## Test:    
#lapply(results, exctract.best3)  

save.output <- function(estimation) {
  max.nodes<- max(sapply(estimation, function(x) length(x[[1]]$Nj)))
  Nj.matrix<- sapply(estimation, function(x) c(x[[1]]$Nj, rep(0, max.nodes-length(x[[1]]$Nj))))
  write.csv(Nj.matrix,file=paste(substitute(result),".nj.estimates.csv",sep=""))
}



save.output2 <- function(result) {
  #capture.output(sapply(estimation, comparison), file=paste(substitute(result),".node.txt",sep=""))
	max.nodes<- max(sapply(estimation, function(x) length(x[[1]]$Nj)))
  Nj.matrix<- sapply(estimation, function(x) c(x[[1]]$Nj, rep(0, max.nodes-length(x[[1]]$Nj))))
  write.csv(Nj.matrix,file=paste(substitute(result),".nj.estimates.csv",sep=""))
  theta.vector<- unlist(sapply(estimation, function(x) c(x[[1]]$theta))  )
  write(theta.vector,file=paste(substitute(result),".theta.vactor.txt",sep=""), ncolumns=1)
}




comparison <- function (estimation) {
  print(signif(cbind(sapply(estimation, function(x) return(x$Nj)))))
  print(sapply(estimation, function(x) sum(x$Nj)))
  print(sapply(estimation, function(x) x$theta))
  print(sapply(estimation, function(x) x[[4]]$value))
}
