resample = TRUE
method = "logistic"
plotToFile = TRUE
makeStandard = TRUE
library("sn")

head1 = c("Pleiotropic","Null","First_Only", "PleioMissing")
head2 = c("(Score-test)","(Multi-inverse-logistic)", "(Single-inverse-logistic)", "(Single-standard)")
#"(Multi-inverse-probit)","(Multi-inverse-cloglog)","(Multi-inverse-cauchit)",


run = TRUE
if(!run){
plotToFile = FALSE

 modMult=100
 sampleIndex = 2
 randDir = FALSE
 phenoMissing = 0.00
 sampSize = 50
 numRep = 1
 vall = 1:2
 p = 0.1
 beta =0.05  #fraction of variance
 vall = 1:2

#beta =0.02  #fraction of variance
#sampSize = 1000  
#p = 0.5  # allele freq
#numRep = 3 
#vall = 1:3
#SBP DBP BMI GLU INS HDL LDL TG CRP

}



degF = 1


library(MASS)
#library(sn)

getA<-function(sigma, rho2){
 rho = sqrt(rho2)
	A = diag(sigma)
	A[1,2] = rho*sigma[1]*sigma[2]
	A[2,1] = A[1,2]
A
}
###FUNCTIONS#

##rho2 is squared correlation coffeciient on error
sampJoint<-function(x,betaSa,ch){
 	
	
	num = length(x)
	z1 = matrix(rnorm(num*length(betaSa),0,1),ncol =length(betaSa), nrow = num)
	#print(paste(dim(ch),dim(z1)))
	

 res = z1 %*% ch
       # betas = matrix(0,nrow = dim(res)[1], ncol = dim(res)[2])
        for(i in 1:length(beta)){
#	  betas[,i] = x*betaSa[i]
	  	res[,i] = res[,i]  + x*beta[i]
        }
#
	res
	#x1 = mu1 + sigma1*z1
	#x2 = mu2 + sigma2*(z1*rho + z2*(sqrt(1-rho^2)))
	#cbind(x1,x2)
}

standardize<-function(LRR,mult){
colmn=apply(LRR,2,mean,na.rm=TRUE)
colsd=apply(LRR,2,sd,na.rm=TRUE)

sLRR=LRR
for (i in 1:dim(LRR)[2]) {
  sLRR[,i]=(LRR[,i]-colmn[i])/colsd[i]
  }
if(sampleIndex>0){
 sLRR[,sampleIndex] = mult*sLRR[,sampleIndex]
}
#risk = sLRR %*% beta
#cat = rep(NA,length(risk))
#for(i in 1:length(risk)){
sLRR  
}


#}





sampleGeno<-function(p,n){
 x =sample((0:(length(p)-1)),n,prob = p,replace=TRUE)
}



estSigma <- function(x){
  d = ncol(x)
  v = matrix(ncol=d,nrow=d)
  for( i in 1:d ){
    ptrI = !is.na(x[,i])
    v[i,i] = var(x[ptrI,i])
    if( i < d ){
      for( j in (i+1):d ){
        ptrJ = !is.na(x[,j])
        ptr = as.logical(ptrI*ptrJ)
        v[i,j] = cov( x[ptr,i], x[ptr,j] )
        v[j,i] = v[i,j]
      }
    }
  }
  return(v)
}

#fillMissing <- function( x, m, lambda ){
#  xx = x
#  ptr.ind = which( apply(is.na(x),1,sum)>0 )
#  for( ii in 1:length(ptr.ind) ){
#    print(ii)
#    i = ptr.ind[ii]
#    ptr.miss = which(is.na(x[i,]))
#    fill =  m[ptr.miss] - solve(lambda[ptr.miss,ptr.miss]) %*% lambda[ptr.miss,-ptr.miss] %*% t(as.matrix(x[i,-ptr.miss] - m[-ptr.miss]))
#    xx[i,ptr.miss] = fill
#  }
#  return(xx)
#}

fillMissing <- function( x, m, lambda ){
  xx = x
  ptr.ind = which( apply(is.na(x),1,sum)>0 )
  for( ii in 1:length(ptr.ind) ){
    i = ptr.ind[ii]
    ptr.miss = which(is.na(x[i,]))
    fill =  m[ptr.miss] - solve(lambda[ptr.miss,ptr.miss]) %*% lambda[ptr.miss,-ptr.miss] %*% as.matrix(x[i,-ptr.miss] - m[-ptr.miss])
    xx[i,ptr.miss] = fill
  }
  return(xx)
}

#IN
#y - values 0,1,2 (also works for binary outcomes)
#x - covariates/phenotypes

#OUT
#z statistic for each covariate/phenotype
#p-value for each covariate/phenotype
#p-value for overall model fit

ordTest <- function(  y,x){
  cats = max(y) + 1
  x = x - apply(x,2,mean,na.rm=T)
  d = ncol(x)
  gamma = as.vector(table(y)/length(y))
  phi = log(gamma)
  U = matrix(ncol=1,nrow=d)
  num1 = list(length=d)
  for( j in 1:d ){
    U[j] = sum(y*x[,j])
    num1[[j]]=0
    for( i in 1:cats ){
      num1[[j]] = num1[[j]] + x[,j] * (i-1) * exp(phi[i])
    }
    U[j,1] = U[j,1] - sum(num1[[j]])
  }

  I = matrix(ncol=d,nrow=d)
  num2 = matrix(ncol=d,nrow=d)
  for( j in 1:d ){
    for( k in j:d ){
      num = 0
      for( i in 1:cats ){
        num = num + (i-1)*(i-1)*x[,j]*x[,k]*exp(phi[i])
      }
      I[j,k] = sum(num1[[j]]*num1[[k]])
      I[j,k] = I[j,k] - sum(num)
      I[k,j] = I[j,k]
    }
  }
  #ret = list()
  #ret[[1]] = U/sqrt(diag(-I))
  #ret[[2]] = 2*(1-pnorm(abs(U/sqrt(diag(-I)))))
#  ret[[3]] = 1 - pchisq(t(U) %*% solve(-I) %*% U, d)
  c(2*(1-pnorm(abs(U/sqrt(diag(-I))))),
	1 - pchisq(t(U) %*% solve(-I) %*% U, d))

  #return(ret)
#  return( 1 - pchisq(t(U) %*% solve(-I) %*% U, d) )
}

getSig<-function(x,pheno, method){
#print(paste(dim(pheno),dim(pheno_)))
x1 =as.factor(x)
numPhen = dim(pheno)[2]
numPhen1 = numPhen+1
pv = rep(1,numPhen1)
beta = rep(0,numPhen1)
if(length(levels(x1))==2){
	res = glm(x1 ~ pheno,family="binomial")
        coeff = summary(res)$coeff[,4]
        coeff1 = summary(res)$coeff[,1]
        llhFull <- logLik(res)
        llhNull <- logLik(update(res, ~1))
        gstat <- -2*(llhNull-llhFull)
        fitP <- pchisq(gstat[1],attr(llhFull,"df")[1]-attr(llhNull,"df")[1],lower.tail=F)
	  for(j in 1:numPhen){
	    pv[j+1] = coeff[j+1]
	    beta[j+1] = coeff1[j+1]
	}
}else if(length(levels(x1))>2){
	res =  polr(x1 ~ pheno,Hess=T,method = method)

	llhFull <- logLik(res) 
	llhNull <- logLik(update(res, ~1))
        gstat <- -2*(llhNull-llhFull)
        fitP <- pchisq(gstat[1],attr(llhFull,"df")[1]-attr(llhNull,"df")[1],lower.tail=F)
	tvals =  summary(res)$coefficients[,3];   
  
  for(j in 1:numPhen){
           pv[j+1] <- 2*(pt(abs(tvals[j]),df.residual(res),lower.tail=F))
	   beta[j+1]<-res$coeff[j]
        }
}

  pv[1] = fitP

rbind(pv,beta)
}

getSigStd<-function(x,pheno){
numPhen = dim(pheno)[2]
pv = rep(1,numPhen+1)
for(i in 1:numPhen){
	ind = which(c(1:numPhen)!=i)
#	print(ind)
#        phen1 = glm(pheno[,i]~pheno[,ind])$res
#	pv[i+1] = summary(glm(phen1~x,family="gaussian"))$coeff[2,4]
   #     pv[i+1] = max(1e-100,summary(glm(pheno[,i]~x+pheno[,ind],family="gaussian"))$coeff[2,4])
	    pv[i+1] = max(1e-100,summary(glm(pheno[,i]~x,family="gaussian"))$coeff[2,4])
	
}
  pv[1] = 1-(1 -min(pv))^degF
pv
}



getSig1<-function(x,pheno,method){
  numPhen = dim(pheno)[2]
  pv = rep(1,numPhen+1)
  for( k in 1:numPhen){
    pv[k+1] = getSig(x,as.matrix(pheno[,k]),method)[1,2]
  }
  pv[1] =  1-(1 -min(pv))^degF
  pv
}


logQQ<-function(pv, conf=NA, yEqx=TRUE, cols=NA, pointcols=NA, title=expression(paste("-log"[10],"QQ")), ymax=NA)
  {
    pv=pv[!is.na(pv)]
   pv[pv<1e-50] = 1e-50   
    if(sum(pv==0)>0)
    {
      print("Warning: Ommitting pvalue=0")
      pv=pv[pv!=0]
    }
    if(length(pv)==0){
      print("warning all pvals were zero")
	}
    pv=sort(pv)
    np=length(pv)

    if(is.na(pointcols))
      {
        pointcols=rep(1,np)
      }
    if(is.na(sum(cols)))
      {

        if(yEqx)
          cols=1:(length(conf)+1)+1
        else
          cols=1:length(conf)+2
      }

    expec=rep(NA, np)
    for(s in 1:np)

      expec[s]=s/(np+1)

    plotpv=rep(FALSE, np)
    plotpv[1:min(1000,length(pv))]=TRUE
    plotpv[floor(10^(runif(n=4000,min=0,max=log10(length(pv)))))]=TRUE

    expec=expec[plotpv]
    pv=pv[plotpv]
    if(is.na(ymax))
      ymax=-log10(min(pv))              
    plot(cbind(-log10(expec), -log10(pv)),  xlim=c(0,log10(np+1)), ylim=c(0,ymax), xlab=expression(paste("-log"[10]," (Expected p-value)")), ylab=expression(paste("-log"[10]," (Observed p-value)")), main=title, pch=20, cex=0.5)
 if(yEqx)
      abline(0,1, lw=2, col=cols[1])
    
    if(!is.na(sum(conf)))
      {
        Obs=rep(NA, np)
        for(c in 1:length(conf))
          {
            for(s in 1:np)
              {
                Obs[s]=qbeta(conf[c], s, np-s+1)
              }
            Obs=Obs[plotpv]
            #expec=expec[plotpv]
            lines(-log10(expec), -log10(Obs), col=cols[c+as.numeric(yEqx)], lw=2)
          }
      }
  }


plotPheno<-function(pheno_to_plot,x, i){
  zero = which(x==0)
  one = which(x==1)
  two = which(x==2)
  plot(pheno_to_plot[,i], pheno_to_plot[,i+1], type='p', col="white")

      title = paste("Beta = ",beta,"Effect") #,paste(t(uinv)[,1],collapse="_"))
      lines(pheno_to_plot[zero,i], pheno_to_plot[zero,i+1], type='p', col=1,main=title,xlab=nmes[1],ylab=nmes[2])
      lines(pheno_to_plot[one,i], pheno_to_plot[one,i+1], type='p', col=2,main = title,xlab=nmes[1],ylab=nmes[2])
      lines(pheno_to_plot[two,i], pheno_to_plot[two,i+1], type='p', col=3,main=title,xlab=nmes[1],ylab=nmes[2])
}

applyMissing<-function(vec,rate){
 vec[sample(c(TRUE,FALSE),length(vec), repl = TRUE,c(rate,1-rate))] = NA
 vec
}

getNaInd<-function(pheno){
apply(is.na(pheno),1,sum)>0 
}
getMissingMatrix<-function(pheno){
phenom = pheno
if(phenoMissing>0){ 
phenom = apply(pheno,2,applyMissing, phenoMissing)
 inds = apply(is.na(phenom),1,sum)<dim(pheno)[2]
 v = estSigma(phenom[inds,])
 m = apply(phenom[inds,],2,mean,na.rm=TRUE)
 phenom[inds,] = fillMissing(phenom[inds,],m,solve(v))
 phenom
}
phenom
}

runf<-function(x, beta0,beta1,sampind, corrM,corrM1,uinv,iter,nmes,pdist){
mu_null = x*0
#print(max(x))
betaNull = rep(0,dim(corrM)[2])
betaFirst = rep(0,length(betaNull))
betaFirst[sampind] = beta1
betaN = rep(0,length(betaNull))
betaN[1] = beta0

##same effect on bothres = read.table("metabolic-residuals.txt",stringsAsFactors=F)

pheno_joint = t(uinv %*% t(sampJoint(x,betaN,corrM))) #+means
##no effect on either
pheno_null = sampJoint(x,betaNull,corrM)

#effect only for first one
pheno_firstOnly = sampJoint(x, betaFirst,corrM1)


#conditional effect
#pheno_cond = sampJoint(x,betaN,corrM)

if(iter==1){
pheno_to_plot = pheno_firstOnly
  dimnames(pheno_to_plot)[[2]] = nmes
  print(nmes)
  if(plotToFile){
   pdf(file="example.pdf")
}
  len1 = dim(pheno_firstOnly)[2]-1
  for(i in 1:len1){
	plotPheno(pheno_firstOnly,x,i)
  }
  if(plotToFile){
  dev.off();
}
}

pheno_missing = getMissingMatrix(pheno_joint)
inds2 = apply(is.na(pheno_missing),1,sum)<dim(pheno_missing)[2]

xrand = x
if(resample){
 xrand = sampleGeno(pdist,length(x))
 pheno_null = pheno_firstOnly
}

#print(naind)
#x = x_[!naind]
#pheno = as.matrix(pheno_[!naind,])

method = c("logistic","probit","cloglog","cauchit")

results = rbind(
ordTest(x,pheno_joint),
ordTest(xrand,pheno_null),
ordTest(x,pheno_firstOnly),
ordTest(x[inds2],pheno_missing[inds2,]),
getSig(x,pheno_joint,method[1])[1,],
getSig(xrand,pheno_null,method[1])[1,],
getSig(x,pheno_firstOnly,method[1])[1,],
getSig(x[inds2],pheno_missing[inds2,],method[1])[1,],
#getSig(x,pheno_joint,method[2])[1,],
#getSig(xrand,pheno_null,method[2])[1,],
#getSig(x,pheno_firstOnly,method[2])[1,],
#getSig(x[inds2],pheno_missing[inds2,],method[2])[1,],
#getSig(x,pheno_joint,method[3])[1,],
#getSig(xrand,pheno_null,method[3])[1,],
#getSig(x,pheno_firstOnly,method[3])[1,],
#getSig(x[inds2],pheno_missing[inds2,],method[3])[1,],
#getSig(x,pheno_joint,method[4])[1,],
#getSig(xrand,pheno_null,method[4])[1,],
#iggetSig(x,pheno_firstOnly,method[4])[1,],
#getSig(x[inds2],pheno_missing[inds2,],method[4])[1,],
getSig1(x,pheno_joint,method[1]),
getSig1(xrand,pheno_null,method[1]),
getSig1(x,pheno_firstOnly,method[1]),
getSig1(x[inds2],pheno_missing[inds2,],method[1]),
getSigStd(x,pheno_joint),
getSigStd(xrand,pheno_null),
getSigStd(x,pheno_firstOnly),
getSigStd(x[inds2],pheno_missing[inds2,])
)
if(results[7,1]<1e-10){
  print(paste("low pv ",results[7,1],sampind))
  print(corrM1)
  print(betaFirst)
}
results
}

empirical1<-function(x,pnull){
max(1e-5,(1+length(which(pnull<x)))/length(pnull))
}

trans<-function(x){
log10(x)
}

empirical<-function(mat, null_ind){
#print(dim(mat))
pnull = mat[null_ind,]
pnull = pnull[order(pnull)]

## fit = cp.to.dp(sn.mle(y=trans(pnull),plot.it = FALSE)$cp)
##  mat1= apply(trans(mat),c(1,2), psn,dp=fit)

mat1 = apply(mat, c(1,2),empirical1,pnull)
mat1
}

empiricalall<-function(pvals){
len2 = length(head2)
len1 = length(head1)
for(k in 1:length(pvals)){
	for(i in 1:len2){
		start = len1*(i-1)+1
		end = len1*i
  		pvals[[k]][start:end,] = empirical(pvals[[k]][start:end,],2)
	}	
}
pvals
}

getDistr<-function(pdist,numRep,res1,v,beta,sampSize,nmes, randDir){
leng = dim(res1)[2]
numCat = length(head1)*length(head2)
numPhen1 = length(v)+1
pvals<-vector("list",numPhen1)
for(k in 1:numPhen1){
 pvals[[k]] <- matrix(0,nrow =numCat,ncol=numRep)
}
u = gs(v)
uinv = t(u)
res2 = t(u %*% t(res1))
covar1 = cov(res1,use = "complete.obs")
#means = apply(res1,2,mean)
corrM1 = chol(covar1)
covar = cov(res2,use = "complete.obs")
corrM=  chol(covar)
rho2 = covar[1,2]/(covar[1,1]*covar[2,2])
beta0 = sqrt(covar[1,1])*beta

sampind = sampleIndex
                if(sampleIndex < 0){
			sampind = sample(leng,1);
                }
beta1 = sqrt(covar1[sampind,sampind])*beta
for(i in 1:numRep){

	if(randDir){
		v = runif(leng)-0.5
		u = gs(v)
		uinv = t(u)
		res2 = t(u %*% t(res1))
		
		covar = cov(res2,use = "complete.obs")
		corrM=  chol(covar)
		rho2 = covar[1,2]/(covar[1,1]*covar[2,2])
		beta0 = sqrt(covar[1,1])*beta
                if(sampleIndex < 0){
			sampind = sample(leng,1);
                }
               	beta1 = sqrt(covar1[sampind,sampind])*beta
        }
	x = sampleGeno(pdist,sampSize)
 
	res = runf(x,beta0,beta1,sampind, corrM,corrM1,uinv,i,nmes,pdist)
	print(paste(dim(res),numPhen1, length(pvals)))
	for(k in 1:numPhen1){
	   pvals[[k]][,i] = res[,k]  
	}	
}
pvals
}

#### PLOTTING RESULTS
plotting<-function(pvals,mod){
header1 = matrix(nrow = length(head1), ncol = length(head2))
for(i in 1:length(head1)){
  for (j in 1: length(head2)){
	header1[i,j] = paste(head1[i], head2[j])
	}
}
header = as.vector((header1))
numPhen1 =length(pvals)
for(j in 1:length(head2)){
 pdf(file=paste(mod,paste(head2[j],"pdf",sep="."),sep=""))
 start = (j-1)*length(head1) +1
 end = j*length(head1)
 for(i in start:end){
       for(k in 1:numPhen1){
           logQQ(pvals[[k]][i,], title = paste(header[i],"pheno",(k-1)))
	}	
	
 }
 dev.off()
}
}


gs<-function(v){
len = length(v)
basis=diag(rep(1,len))
svd(cbind(v,basis))$u
}
 
#covar = getA(c(0.5,0.5),rho2)


if(run){
index = vall
v = vall
#degF = max(1,length(v)-3)  ## approx correct for using 9 metabolic traits but does not affect empirical distribution anyway which is based on ordering
for(i in 1:length(vall)){
if(i%%2==0) v[i] = 1
else v[i] = -1
}
print("index")
print(index)
#index = which(vall!=0)
#v = vall[index]

res = read.table("metabolic-residuals.txt",stringsAsFactors=F)
dims = dim(res)
nmes1 = res[1,(2:dims[2])]
res = res[2:(dims[1]),(2:dims[2])]
 for(i in 1:(dim(res)[2])){
res[,i] = as.numeric(res[,i])
}
res1 = res[,index]
nmes = nmes1[index]
sv = sqrt(svd(res1)$d)
degF = length(which(sv/sv[1] >0.12))
#degF = sum(sv$d/sv$d[1])
#print(paste("degF",degF,sum(sv/sv[1])))

q = 1-p
if(makeStandard){
res1 = standardize(res1,modMult)
}
#getDistr<-function(p,numRep,res1,v,beta0,sampSize,nmes, randDir){
pdist = c(q^2, 2*p*q,p^2)
pvals = getDistr(pdist,numRep,res1,v,beta, sampSize,nmes,TRUE)
plotting(pvals,"")
save(pvals,file="objects.r_save")
plotting(empiricalall(pvals),"adjusted_")
}
