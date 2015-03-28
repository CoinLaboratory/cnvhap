library(MASS)
library(abind)
library(epitools)
geno = NULL
genotmp = NULL
weights = NULL

###NOTE: TO DEBUG THE SCRIPT, YOU NEED TO MAKE A DIRECTORY CONTAINING FILES
### pheno.txt, limit.txt, exclude.txt,  Samples (unzip from zip file), Name (unzip from zip file) ,
#snp (which is also from zip)
debug = TRUE


##USER DEFINED VARS
#first columns from zip file to be associated
step = 1
maf_thresh = 0.001
inverseRegress = FALSE
useContinuous = FALSE
## whether to use normal linear regression for genotypes

###NOTE if we choose inverseRegress, we should be careful with resid for genotypes as outcome (pref covariates)
###
varianceIsCov= NULL ##is covariate, or residual, if null, it is ignored
multiPhen = TRUE
multiGen = FALSE
if(multiPhen | inverseRegress) multiGen = FALSE
fillMissingPhens = FALSE
scoreTest = FALSE
exactTest = FALSE
exactMethod = "wald"
#  method = c("midp", "fisher", "wald", "small"),
##ARE GENOTYPE COVARS TREATED AS RESIDUALS, OR COVARIATES IN REGRESSION
genoCovIsRes = TRUE

#if(scoreTest) fillMissingPhens = TRUE
corThresh = 0.0 ##used to only choose correlated phenotypes
expandData = FALSE
if(multiGen) expandData = FALSE
if(exactTest) expandData = TRUE


rescale = 10  ##how much to rescale in regression (for weights)
if(!expandData) rescale = 1

max_indiv = 1000000



metares = "metaresInvVarianceFixed" 
#metares="metaresFisher"

extralevelp=c(0.0,0.0)
#ids = c("GType")#,"Log","Freq")
#ids = c("Log") #,"Freq")
#ids= c("GType")
ids = c("geno")
#ids = c("Allele[12] - Top")
#ids = c("countAll" ,"state.0","state.2")
#ids = c("Freq") #,"Allele","GType"

baselevel=c(0,0,0,0)  ##this is used for maf calc.  For count all base (ie. normal) is 2
baselevel[ids=="countAll"] = 2
singleOutput=TRUE
zipOutput=FALSE

##now define plate and variance columns for the Samples entry (if they exist)
plate_id = "Plate"
variance_id ="variance"
##define the variance threshold to apply (if the variance_id column is present)
varthresh =0.25


###FUNCTIONS 
expandTodo<-function(tod){
tod1 =  paste(tod[3],"",sep="")
vec = as.numeric(strsplit(strsplit(tod1,"_")[[1]][2],"\\.\\.")[[1]])
len = vec[2] - vec[1] +1
res = matrix("",ncol=3,nrow = len)
for(i in vec[1]:vec[2]) res[i-vec[1]+1,] = c(tod[1:2],paste("factor",i,sep="_"))
res
}
expandTodoMat<-function(tod){
 res =matrix("",ncol=3,nrow = 0)
 for(i in 1:(dim(tod)[1])) res = rbind(res,expandTodo(tod[i,]))
  res
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


fillMissing <- function( x){
  v = estSigma(x)
 m = apply(x,2,mean,na.rm=TRUE)
  lambda = solve(v)
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



convAvgToUnc<-function(vec,n,sc=1){
res = rep(0,n)
fl = floor(vec)
if(fl==vec) res[vec+1] = sc
 if(fl!=vec){
  ceil = fl+1
  val = round(sc*(ceil-vec))
  res[fl+1] = val
  res[ceil+1] = sc-val
 }
 res
}

convAvgToUncM<-function(vec,sc,max){
#print(vec)
 matr = apply(as.matrix(vec),1,convAvgToUnc,max)
 as.vector(t(apply(as.matrix(vec),1,convAvgToUnc,max,sc)))
}

convAvgToUncMat<-function(mat,sc){
 max = max(mat,na.rm=T)+1
 if(max%%1>0) max = floor(max)+1
 apply(mat,2,convAvgToUncM,sc,max)
}

getIncl<-function(phenv, na_cov) !((is.na(phenv) | na_cov))

getFamily<-function(phenv,incl){
 if(length(levels(as.factor(phenv[incl])))>2) "gaussian" else "binomial" 
}

expandCoeff<-function(coeff,ind){
 diff =  max(ind) - dim(coeff)[1]
 if(diff>0) coeff = rbind(coeff,matrix(NA, nrow = diff,ncol = dim(coeff)[2]))
 coeff
}
exp10<-function(x) 10^x

extractSig<-function(res,ind = 2,totweight = 1){
# print("extractSig")
# print(summary(res)$coeff)
 #print(ind)
 vec  = summary(res)$coeff[ind,]
 if(length(vec)==3) vec = c(vec, 2*pt(abs(vec[3]),1,lower.tail=F)) 
 vec[c(1,4)]
}

extractSigP<-function(res,ind,totweight = 1,hasCov=FALSE,family="binomial"){
 #print("extractSigP")
coeff = expandCoeff(summary(res)$coeff,ind)
 #print(coeff)
 #print(family)
 vec  = coeff[ind,,drop=F]
 ll1 = logLik(res)
 if(dim(vec)[2]==3) vec = cbind(vec, 2*pt(abs(vec[ind,3]),1,lower.tail=F))
 if(hasCov)ll2 = logLik(update(res, ~covar,family=family)) else ll2 = logLik(update(res, ~1,family=family))
 pv = pchisq((2*(ll1 - ll2))/(totweight),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=F)
 #print(vec[,4])
 #print(pv)
 cbind(c(vec[,1],NA),c(vec[,4],pv))
}


extractSig1<-function(res,ind=2,totweight=1,hasCov=FALSE,family="binomial"){
 #print("extractSig1")
 ll1 = logLik(res)
 beta = summary(res)$coeff[ind,1]
nullfam = is.null(family)
 if(hasCov & !nullfam) ll2 = logLik(update(res, ~covar,family=family)) else if(!nullfam) ll2 = logLik(update(res, ~1,family=family)) else if(hasCov) ll2 = logLik(update(res, ~covar)) else ll2 = logLik(update(res, ~1))
 df = attr(ll1,"df")[1]-attr(ll2,"df")[1]
# print(paste("df",df))
 pv = pchisq((2*(ll1 - ll2))/(totweight),df,lower.tail=F)
 c(beta,pv) 
}

extractSig1P<-function(res,ind=2,totweight=1,hasCov=FALSE,family="binomial"){
# print("extractSig1P")
 ll1 = logLik(res)
 coeff = summary(res)$coeff
 beta = coeff[ind,1]
 ve = abs(as.numeric(coeff[ind,3]))
 pvs = 2*pt(ve,1,lower.tail=F)
 nullfam = is.null(family)
 if(hasCov & !nullfam) ll2 = logLik(update(res, ~covar,family=family)) else if(!nullfam) ll2 = logLik(update(res, ~1,family=family)) else if(hasCov) ll2 = logLik(update(res, ~covar)) else ll2 = logLik(update(res, ~1))
 pv = pchisq((2*(ll1 - ll2))/(totweight),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=F)
# print(paste("beta",beta))
# print(coeff)
 cbind(c(beta,NA),c(pvs,pv))
}

seFromPBeta<-function(v) abs(v[1])/qnorm(v[2]/2,lower.tail=F)
pFromSEBeta<-function(v) pnorm(abs(v[1]),sd = v[2],lower.tail=F)*2

findInd<-function(vec,header,ind) grep(vec[ind],header)[1]
centralise<-function(vec)  (vec-mean(vec,na.rm=T))/sd(vec,na.rm=T)



expand<-function(matr, repl){
 res = matr
 for(i in 2:repl) res = rbind(res,matr)
 res
}

expand2<-function(matr, repl){
 res = matr
 for(i in 2:repl)res = abind(res,matr,along = 3)
 res
}

nullDim<-function(matr) length(dim(matr))

expandDat<-function(data,repl){
 inds = sapply(data,nullDim)
 data[inds==2] = lapply(data[inds==2],expand,repl)
 data[inds==3] = lapply(data[inds==3],expand2,repl)
 data
}

getGlmRevBinom<-function(vars,phenvec,covar,weights, family) glm(vars ~ phenvec + covar,family=family,weights=weights)

getGlmRevOrd<-function( vars, phenvec, covar,weights) {

 polr(vars ~ phenvec + covar,Hess=T,method="logistic",weights = weights)
}
getGlmRevBinomNoCov<-function( vars, phenvec, weights, family){
  glm(vars ~ phenvec,family=family,weights=weights)
}
getGlmRevOrdNoCov<-function( vars, phenvec, weights){
 polr(vars ~ phenvec,Hess=T,method="logistic",weights = weights)
}
getGlm<-function(phenvec,vars,covar, fam1,weights,offset){
 glm(phenvec~vars+covar,family = fam1,weights = weights,offset = offset)
}

getGlmNoCovar<-function(phenvec,vars, fam1,weights,offset){
  glm(phenvec~vars,weights = weights,offset = offset, family = fam1)
}


calcOffset<-function(phenv,family,pheno_cov){
 offsets1 = rep(NA,dim(pheno_cov)[2])
 incl = phenv[2,]==1
 dimres = dim(pheno_cov)[2]
 if(dimres>0) offsets1[incl] =  glm(phenv[1,incl]~pheno_cov[incl,],family=family)$lin else offsets1[incl] = 0
 offsets1
}

maxCatForOrd = 5

##note subtracting offset only works for continous outcomes
getPvRev<-function(phend, gvar, family, inclu,covar,totweight){ 
 vars = gvar[1,]
 weights = as.numeric(gvar[2,])
include = phend[2,]>0 & inclu & weights > 0 & gvar[3,]>0
 var1 = as.factor(vars[include])
 phenvec = phend[1,include] - phend[4,include]
    emptyCov = dim(covar)[2]==0
binom =  length(levels(var1))<=2
cont =  useContinuous | length(levels(var1))>maxCatForOrd
family = "binomial"
ord = !binom & !cont
if(cont) var1 = vars[include]
if(cont) family = "gaussian"
if(ord)family=NULL
 todo = TRUE
# print(paste("running pv rev",binom,cont,family,ord, totweight))

 if(length(which(include))==0 | var(var1,na.rm=T)==0) todo=FALSE
 if(todo & !ord & !emptyCov) glm = getGlmRevBinom(var1,phenvec,covar[include,],weights[include], family) else if(todo & !emptyCov) glm = getGlmRevOrd(var1,phenvec,covar[include,],weights[include]) else if(todo & !ord) glm = getGlmRevBinomNoCov(var1,phenvec,weights[include],family) else if(todo) glm = getGlmRevOrdNoCov(var1,phenvec,weights[include])
 if(!todo) rep(NA,2) else extractSig1(glm,ind=1,totweight=totweight,hasCov=!emptyCov,family=family)
} #####FINISH




getPvLrr<-function(phend, gvar, family,inclu,covar,totweight){
  phenvec = phend[1,]
  vars = gvar[1,]
  weights = gvar[2,]
  offset = phend[4,]
  include = phend[2,] & inclu & weights > 0 & gvar[3,]>0 
   emptyCov = dim(covar)[2]==0
   todo = TRUE
   if(length(which(include))==0 | var(vars[include],na.rm=T)==0) todo=FALSE
 if(todo & emptyCov) glm = getGlmNoCovar(phenvec[include] ,as.numeric(vars[include]),family,weights[include],offset[include]) else if(todo) glm = getGlm(phenvec[include] , as.numeric(vars[include]), covar[include,],family,weights[include],offset[include])
  if(!todo) rep(NA,2) else extractSig(glm,ind=2,totweight=totweight)
}

tabulateWithWeights<-function(phen,vars,weights,totweight){
 tab = table(list(vars,phen))
 nme = dimnames(tab)
 nme1 = as.numeric(nme[[1]])
 nme2 = as.numeric(nme[[2]]) 
 tab1 = matrix(0,nrow = length(nme1), ncol  = length(nme2))
 for(i in 1:length(nme1)) for(j in 1:length(nme2)) tab1[i,j] = sum(weights[vars==nme1[i] & phen==nme2[j]])
 tab2 = tab1/totweight
 tab2
}



extractSigExact<-function(t){
# print(t) note, not sure this is getting things in right order for good effect size estimate
   or1 = fisher.test(t)
  or = oddsratio(t,method=exactMethod)
#  ind = grep(exactMethod,dimnames(or$p.value))[1]
#  print(or)
  c(log(or$measure[2,1]),or1$p.value)
}

getPvTable<-function(phend, gvar, family,inclu,covar,totweight){
  phenvec = phend[3,]
  vars = as.factor(gvar[1,])
  weights = gvar[2,]
  offset = phend[4,]
  include = phend[2,] & inclu & weights > 0 & gvar[3,]>0
   emptyCov = dim(covar)[2]==0
   todo = TRUE
   if(length(which(include))==0 | var(vars[include],na.rm=T)==0 | family!="binomial") todo=FALSE
 if(!todo) rep(NA,2) else extractSigExact( tabulateWithWeights(phenvec[include],vars[include],weights[include],totweight))
}




getPvLrrMult<-function(vars,phen,families,include,covar,totweight, functi){
   res = matrix(NA, nrow = length(families),ncol = 2)
   fams = levels(as.factor(families))
   for(k in 1:length(fams))   if(length(which(families==fams[k]))>0) res[families==fams[k],]  =  t(apply(phen[families==fams[k],,,drop=F],1,functi,vars, fams[k],include,covar,totweight))
   res
}


getPvLrrAllGen<-function(vars,phen,families,include,covar,totweight,functiAll,functi){
  t(apply(vars,1,functiAll,phen,families,include,covar,totweight,functi))
}

getPvLrrMultiGen<-function(gvar,phend, family,inclu,covar,totweight,functiAll,functi){
t(apply(phend,1,getPvLrrMultiGen1,gvar,family,inclu,covar,totweight,functiAll,functi))
}

getPvLrrMultiGen1<-function(phend, gvar,family,inclu,covar,totweight,functiAll,functi){
  phenvec = phend[1,]
  vars = gvar[,1,]
  weights1 = gvar[1,2,]  
  nainds =  apply(as.matrix(gvar[,3,]),2,min)
  offset = phend[4,]
  include = phend[2,] & inclu & nainds >0 & weights1 > 0
  # & nainds # & weights > 0
  emptyCov = dim(covar)[2]==0
  ind = 1:(dim(vars)[1]) + 1
   todo = TRUE
  if(length(which(include))==0 ) todo=FALSE
  if(todo & emptyCov) glm = getGlmNoCovar(phenvec[include],t(vars[,include]),family,weights[include],offset[include]) else if(todo) glm = getGlm(phenvec[include] , t(vars[,include]), covar[include,],family,weights[include],offset[include])
   if(!todo) rep(NA,2) else extractSigP(glm,ind=ind,totweight=totweight,family=family,hasCov=!emptyCov)
}


getPvRevPleio<-function(gvar,phen,families,include,covar,totweight,functi){
 phenvec = phen[,1,]-phen[,4,]
 vars = gvar[1,]
 weights = gvar[2,]
 nainds =  gvar[3,]
 phenincl = apply(as.matrix(phen[,2,]),2,min)
 include =  phenincl>0 & inclu & weights>0 & nainds >0
  var1 = as.factor(vars[include])
 binom = length(levels(var1))<=2
 cont = useContinuous | length(levels(var1))>maxCatForOrd

 family = "binomial"
 ord = !binom & !cont
 if(cont) family = "gaussian"
if(ord)family=NULL
if(cont) var1 = vars[include]
 emptyCov = dim(covar)[2]==0
 ind = 1:(dim(phenvec)[1])
 if(!ord) ind = ind+1

 todo = TRUE
 if(length(which(include))==0 | var(var1,na.rm=T)==0) todo=FALSE
 if(todo & !ord & emptyCov) glm = getGlmRevBinomNoCov(var1,t(phenvec[,include]),weights[include],family) else if(todo & emptyCov) glm = getGlmRevOrdNoCov(var1,t(phenvec[,include]),weights[include]) else if(todo & !ord) getGlmRevBinom(var1,t(phenvec[,include]),covar[include],weights[include],family) else if(todo)  glm = getGlmRevOrd(var1,t(phenvec[,include]),covar[include],weights[include])
    if(!todo) rep(NA,2) else extractSig1P(glm,ind=ind,totweight=totweight,family=family)
}

getIMatrix<-function(d,phi,cats,num1,x){
 I = matrix(ncol=d,nrow=d)
 for( j in 1:d ){
    for( k in j:d ){
      num = 0
      for( i in 1:cats ) num = num + (i-1)*(i-1)*x[,j]*x[,k]*exp(phi[i])
      I[j,k] = sum(num1[[j]]*num1[[k]])
      I[j,k] = I[j,k] - sum(num)
      I[k,j] = I[j,k]
    }
 }
 I
}
getUNum1Matrix<-function(d,y,x,phi,num1,cats){
  U = matrix(ncol=1,nrow=d)
  num1 = list(length=d)
  for( j in 1:d ){U[j] = sum(y*x[,j]);    num1[[j]]=0; for( i in 1:cats ) num1[[j]] = num1[[j]] + x[,j] * (i-1) * exp(phi[i]); U[j,1] = U[j,1] - sum(num1[[j]])}
  list(U,num1)
}



ordTest2 <- function(  y,x){
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
  zstats  = U/sqrt(diag(-I))
  pvs = 2*(1-pnorm(abs(U/sqrt(diag(-I)))))
  overallp = 	1 - pchisq(t(U) %*% solve(-I) %*% U, d)
  cbind(c(NA,zstats),c(overallp,pvs))
}
ordTes1 <- function(gvar,phend,families, inclu,covar,totweight,functi){
  include =  apply(as.matrix(phend[,2,]),2,min) & inclu & as.numeric(gvar[2,])>0 &  gvar[3,]
   y = as.numeric(gvar[1,include])
  x = t(as.matrix(phend[,1,include] - phend[,4,include]))
  ordTest(y,x)
}
ordTest <- function(gvar,phend,families, inclu,covar,totweight,functi){
  include =  apply(as.matrix(phend[,2,]),2,min) & inclu & as.numeric(gvar[2,])>0 &  gvar[3,]
  y = as.numeric(gvar[1,include])
  x = t(as.matrix(phend[,1,include] - phend[,4,include]))
  print(is.numeric(x[,1]))
  cats = max(y) + 1
  x = x - apply(x,2,mean,na.rm=T)
  d = ncol(x)
  gamma = as.vector(table(y)/length(y))
  phi = log(gamma)
  Unum1 = getUNum1Matrix(d,y,x,phi,num1,cats)
  U = Unum1[[1]]
  num1 = Unum1[[2]]
  I = getIMatrix(d,phi,cats,num1,x)
  zstats  = U/sqrt(diag(-I))
  pvs = 2*(1-pnorm(abs(U/sqrt(diag(-I)))))
  overallp = 	1 - pchisq(t(U) %*% solve(-I) %*% U, d)
  cbind(c(NA,zstats),c(overallp,pvs))
}




nonNumeric<-function(vec) sum(is.na(as.numeric(vec)))==length(vec)
nonNumeric1<-function(mat) apply(mat,2,nonNumeric)
makeNumeric<-function(vec) as.numeric(as.factor(vec))
getHist<-function(vec, spl){ 
le = vec < spl[1]
gt = vec > spl[length(spl)]
#print(spl);print(max(vec));print(min(vec))
c(length(which(le)),hist(as.numeric(vec[!le & !gt]),br=spl,plot=F)$count,length(which(gt)))
}
findmax<-function(vec) max(0,max(as.numeric(vec),na.rm=TRUE))
iscase<-function(vec) !is.na(vec) & vec>mean(vec,na.rm=TRUE)
minOneMinus<-function(x) min(x,1-x)
getMaf<-function(vec,base,noallele){
if(length(levels(as.factor(vec[!is.na(vec)])))>5) var(vec,na.rm=T) else minOneMinus(sum(abs(vec[!is.na(vec)] -base))/(length(which(!is.na(vec)))*noallele))
}
getDim<-function(arraynames) c(length(arraynames[[1]]), length(arraynames[[2]]), length(arraynames[[3]]))
getArray<-function(arraynames) array(NA,dim= getDim(arraynames), dimnames = arraynames)


getOnes<-function(vec,f) f==vec[1]
convertFactorToMatrix<-function(f,nme){
   res= apply(as.matrix(levels(as.factor(f))),1,getOnes,f)
   dimnames(res) = list(names(f),paste(nme,levels(as.factor(f)),sep="_"))
   res
}
makeFactor<-function(mat1,nme){
  append = dim(mat1)[2]>0
  fact = factor(apply(mat1,1,paste,collapse=".")) 
  lev = levels(fact) 
  matrix = convertFactorToMatrix(fact,paste(nme,sep="."))
  nmes = dimnames(matrix)[[2]]
  txt = dimnames(matrix)[[2]]
  index = setdiff(seq(length(txt)),grep("_NA",txt))
  res = matrix[,index,drop=F]
  dimnames(res) = list(NULL,nmes[index])
  if(append)  res = cbind(res, rep(TRUE,dim(res)[1]))
  if(append)   dimnames(res) = list(NULL,c(nmes[index],"all")) 
  res 
}
 #if(dim(mat1)[2]>0){
 #   res = cbind(res, rep(TRUE,dim(res)[1]))
 #   dimnames(res) = list(NULL,c(nmes[index],"all")) 
 # }

fisherPvalue<-function(p){ pchisq(-2*sum(log(p)),df = 2*length(p),lower.tail=FALSE)}

metaresFisher<-function(matr){
  inds = !is.na(matr[,2] )
 fisherp = fisherPvalue(matr[inds,2])
 c(mean(matr[inds,1]),fisherp)
}

metaresInvVarianceFixed<-function(matr){
 inds = !is.na(matr[,2])
 metres = metagen(matr[inds,1,drop=F],apply(matr[inds,,drop=F],1,seFromPBeta))
 pv =  pFromSEBeta(c(metres$TE.fixed,metres$seTE.fixed))
##c(metres$Q,pchisq(metres$Q,1,lower.tail=F)) 
c(metres$TE.fixed,pv)
}
#  matr = mat[1:(dim[1]-1),]

makeQuantile<-function(vec)  qqnorm(vec,plot=F)$x
makeTopTail<-function(vec,perc1){
 perc = as.numeric(perc1)/100
 vec1 = rep(NA,length(vec))
 quants = ecdf(sort(vec))(vec)
 vec1[quants<=perc] = 0
 vec1[quants>=1-perc] = 1
 vec1
}

makeThresh<-function(vec,thresh){
 vec1 = rep(NA,length(vec))
 vec1[vec<=thresh[1]] = 0
 vec1[vec>=thresh[2]] = 1
 vec1
}

convertMatrixToFactor<-function(mat){
 nme = as.matrix(data.frame(strsplit(dimnames(mat)[[2]],"_")))
 res = rep("",dim(mat)[1])
 for(i in 1:length(res)) res[i] = paste(nme[2,mat[i,]==1],collase="")
 res
}

joinResVec<-function(vec,form) paste(sprintf(form,vec),collapse=";")

joinRes<-function(mat,ncol,form){
 if(is.null(dim(mat))) mat = matrix(mat,ncol=ncol)
 res1 =  apply(mat,2,joinResVec,form)
 res1
}


fixNonNumeric<-function(pheno){
 nonnumvec = nonNumeric1(pheno)
 pheno1 = matrix(0,nrow=dim(pheno)[1],ncol=dim(pheno)[2])
 pheno1[,!nonnumvec] = apply(as.matrix(pheno[,!nonnumvec]),2,as.numeric)
 pheno1[,nonnumvec] = apply(as.matrix(pheno[,nonnumvec]),2,makeNumeric)
 dimnames(pheno1) = dimnames(pheno)
 pheno = pheno1
}
##different transformations of data (currently not used)
proc<-function(vec) (vec - 0.5)*vec*(1-vec)
proc1<-function(vec) vec
cnt<-function(str) if(length(grep("N",str))>0) NA else nchar(gsub("A","",str))
cntVec<-function(str) apply(as.matrix(str,nrow = length(str)),1,cnt)

numLev<-function(vec) length(levels(as.factor(vec)))


##for hwe
vecm = c(2,1,0)
allFreq<-function(y) (vecm[1:length(y)] %*% y) / (2*(sum(y)))
expDist<-function(y,p)   c(p^2, 2*p*(1-p), (1-p)^2)[1:length(y)]*sum(y)
hwestat<-function(obs,exp) sum((obs-exp)^2/exp) 
hwep<-function(y) if(length(which(y>0))<=1) 1.0 else 1-pchisq(hwestat(y,expDist(y,allFreq(y))),1) 
hwepall<-function(y) min(hwep(y[1:min(3,length(y))]), hwep(y[3:min(5,length(y))]),na.rm=TRUE) 
mafCount<-function(y) max(minOneMinus(allFreq(y[1:min(3,length(y))])), minOneMinus(allFreq(y[3:min(5,length(y))])),na.rm=TRUE)

##for extracting ids
getin<-function(stend,vec) vec>=stend[1] & vec<stend[2]   #NEW
getin1<-function(genvec,spl) apply(cbind(spl[1:(length(spl)-1)],spl[2:length(spl)]),1,getin, genvec) #NEW
getin2<-function(genvec,spl,cc) apply(cbind(spl[1:(length(spl)-1)],spl[2:length(spl)]),1,getin, genvec) & cc
extractIds<-function(vec,nme) {paste(nme[vec],collapse=",")[[1]]}
extractIdsVec<-function(vec,spl,caseInd,nme) apply(getin2(as.numeric(vec),spl,caseInd),2,extractIds,nme)
extractIdsMat<-function(mat,spl,caseInd) apply(mat,1,extractIdsVec,spl,caseInd,dimnames(mat)[[3]])


 mod<-function(vec){ vec[vec==""] = 0; vec[vec=="-"] = 1000 - sum(as.numeric(vec[vec!="-"]),na.rm=T); (as.numeric(vec)%*% 0:2)/1000}

mod1 <-function(vec){
vec[vec==""] = 0
 vec[vec=="-"] = 1000 - sum(as.numeric(vec[vec!="-"]),na.rm=T);
#vec[vec<200]  =0
vec1 = floor(rescale*(as.numeric(vec)/1000))
ind = vec1==max(vec1)
vec1[ind] = vec1[ind]+ rescale - sum(vec1)
vec1
}

fixGenoCols<-function(geno, gCols){
 geno1 = matrix(0,nrow=dim(geno)[1],ncol=dim(geno)[2])
 geno1[,gCols] = apply(geno[,gCols,drop=F],2,cntVec)
 geno1[,-gCols] = apply(geno[,-gCols,drop=F],2,as.numeric)
 geno1
}

mergeAlleleCols<-function(snp,cols){
len = 1+ dim(snp)[2] - length(cols) 
snp1 = matrix(0,nrow = dim(snp)[1],ncol =len) 
snp1[,1] = cntVec(apply(snp[,cols],1,paste,collapse=""))
matr = snp[,-cols,drop=F]
nme1 = rep("",len)
nme1[1] = "geno"
if(len>1) nme1[2:len] = dimnames(snp)[[2]][-cols]
if(len>1) for(i in 2:len) snp1[,i] = matr[,i-1]
dimnames(snp1) = list(dimnames(snp)[[1]],nme1)
#print(snp1[1:5,])
snp1
}

#getuniq<-function(samples,i) samples[match(unique(samples[,i]),samples[,i]),]



#######

#ONLY TO USE TO DEBUG SCRIPT.  MUST HAVE EXTRACTED Name and Samples and snp file
if(debug){
 max_indiv = 100000
 names  = read.table("Name",fill=T,as.is=T,header=FALSE,sep="\t")
 samples = read.table("Samples",header=F,as.is=T)
print(samples[1:10,])
 dimnames(samples)[[2]] = names[3,1:dim(samples)[2]]
 dataUnc = TRUE
 #if(!expandData) dataUnc = FALSE
 if(!dataUnc) snp = read.table("snp",as.is = T,sep="\t")
 if(dataUnc){
  snp = read.table("snp",as.is = T,sep=",")
  if(expandData) weights = as.vector(t(apply(snp,1,mod1)))
  snp = as.matrix(apply(snp,1,mod))
if(inverseRegress & !expandData) snp = apply(snp,c(1,2),round)
 }
 print("SNP!!!!!")
 print(head(snp))
 dimnames(snp) = list(1:length(samples[,1]),names[1,1:dim(snp)[2]])
# dimnames(snp)[[2]] =names[1,1:dim(snp)[2]]
geno_header = dimnames(snp)[[2]]
geno = snp
genotmp = snp[,1,drop=F]
if(max_indiv<dim(geno)[1]) geno = geno[1:max_indiv,,drop=F]
if(expandData & !dataUnc) weights = convAvgToUncMat(geno,rescale) 
 workingDir = "."
# print(weights)
}

#
#
#END OF DEBUGGING SECTION
#JRI will define maxg, workingDir, samples and geno_header 







###INITIALISATION SECTION
#JRI will provide a  `workingDir' variable
#note any variable with`datanme' is replaced with name of directory!!!!###


##READ PHENO FILES

if(dim(samples)[1]>max_indiv) samples = samples[1:max_indiv,,drop=FALSE]

phenoFile = paste(workingDir,"pheno.txt",sep="/")
limitFile = paste(workingDir,"limit.txt",sep="/")
excludeFile = paste(workingDir,"exclude.txt",sep="/")


###NOW PREPARING PHENO FILE
pheno = as.matrix(read.table(phenoFile,header=T,row.names=1,as.is=T))

##following is to convert character columns to strings



pheno = fixNonNumeric(pheno)

##APPLY exclude.txt if it exists
if(length( grep("exclude.txt",dir(workingDir)))>0) pheno[match(read.table(excludeFile,as.is=T,header=T)[,1],dimnames(pheno)[[1]]),] = NA
	


##apply limit.txt to extract covars and phenotypes we are interested in
phenNames = dimnames(pheno)[[2]]
if(length(grep("limit.txt",dir(workingDir)))>0) todo_ = read.table(limitFile,as.is=T,fill=T)
#if(dim(todo_)[2]>=3 & length(grep("..", todo_[,3])) >0){
# ind = grep("..", todo_[,3])
# todo_ = rbind(expandTodoMat(as.matrix(todo_[ind,,drop=F])),as.matrix(todo_[-ind,,drop=F]))
#}
todo = todo_[grep("^pheno",todo_[,1]),2]
excl = todo_[grep("^exclpheno",todo_[,1]),2]
covar = todo_[grep("^covar",todo_[,1]),2]
stratify=todo_[grep("^strat",todo_[,1]),2]
resid=todo_[grep("^resid",todo_[,1]),2]
genocov=todo_[grep("^genocov",todo_[,1]),2]

if(length(stratify)>0 & metares=="metaresInvVarianceFixed") library(meta)

index = c()
index_cov = c()
index_resid = c()
index_strat = c()
exclindex = c()
for ( ik in todo) index =c(index,which(phenNames==ik))

#for ( ik in todo) index =c(index,grep(ik,phenNames))
for ( ik in excl) exclindex =c(exclindex,grep(ik,phenNames))
for ( ik in covar) index_cov = c(index_cov,which(phenNames==ik))
for ( ik in resid) index_resid = c(index_resid,which(phenNames==ik))
for ( ik in stratify) index_strat = c(index_strat,which(phenNames==ik))

#index = unique(index)
exclindex = unique(exclindex)
print(paste("indices",exclindex))
print(paste("indices",index))

if(length(exclindex)>0)index = index[-match(exclindex,index)]

index_cov = unique(index_cov)
index_resid = unique(index_resid)
index_strat = unique(index_strat)

print(paste("indices",exclindex))
print(paste("indices",index))


if(dim(todo_)[2]>2)phenoTrans = todo_[grep("^pheno",todo_[,1]),3] else phenoTrans = NULL
if(length(genocov)>0) genoTrans = apply(todo_[grep("^genocov",todo_[,1]),,drop=F],1,findInd,geno_header,3) else genoTrans = NULL

if(length(genocov)==0) genotmp=NULL

if(length(index)==1) multiPhen = FALSE
if(length(index_strat)>0 | length(index_cov) > 0) exactTest = FALSE
if(exactTest){
 multiPhen = FALSE
 multiGen = FALSE
 inverseRegress = FALSE
}
funct = getPvLrr
pvFunct = getPvLrrMult
pvFunctMult = getPvLrrAllGen


if(inverseRegress) funct = getPvRev


if(multiPhen) pvFunct = getPvRevPleio
if(scoreTest & multiPhen) pvFunct = ordTest
if(multiGen) pvFunctMult = getPvLrrMultiGen
if(exactTest) funct = getPvTable

##define the regression function to use (see function section)

##if no pheno index was specified, we consider all columns in the phenotype file
if(length(index)==0) index = 1:(dim(pheno)[[2]])



# toStrat = as.matrix(rep(TRUE,dim(pheno)[[1]]))
# stratNames = c("")

 toStrat = pheno[,index_strat,drop=F]
 stratNames = dimnames(pheno)[[2]][index_strat]



pheno_cov = pheno[,index_cov,drop=F]
pheno_resid =pheno[,index_resid,drop=F]
pheno = pheno[,index,drop=F]
dimnames(pheno)[[2]] = phenNames[index]




if(multiPhen){
 if(fillMissingPhens) pheno = fillMissing(pheno)
 corel = (cor(pheno,use="complete.obs"))^2
 index1 = corel[1,]>corThresh
 pheno = pheno[,index1,drop=F]
 index = index[index1]
 dimnames(pheno)[[2]] = phenNames[index]
}
 



##find the variance and plate (if it exists)
var_id = which(dimnames(samples)[[2]]==variance_id)
plate_ind = which(dimnames(samples)[[2]]==plate_id)
plate = matrix(0, ncol=1, nrow = dim(samples)[1])
if(length(plate_ind)>0){
  plate = makeFactor(as.matrix(samples[,plate_ind]),c("plate"))
}


##apply variance filter
variance = samples[,var_id]
if(length(var_id)>0) varna =  variance>varthresh else varna = rep(FALSE,dim(samples)[1]) 
indiv = samples[,1]
###now match up the correct entries in pheno file, with the entries from the Samples file
mat = match(samples[,1],dimnames(pheno)[[1]])

pheno1 = pheno[mat,,drop=F]
dimnames(pheno1)[[2]] = dimnames(pheno)[[2]]
pheno_cov1 = pheno_cov[mat,,drop=F]
pheno_resid1 = pheno_resid[mat,,drop=F]
stratMatrix1 = makeFactor(toStrat[mat,,drop=F],stratNames)


if(!is.null(phenoTrans)){
  print(phenoTrans)
 qqfam =  grep("^quantile",phenoTrans)
 if(length(qqfam)>0) pheno1[,qqfam] = apply(pheno1[,qqfam,drop=F],2,makeQuantile)
 toptail = grep("^toptail",phenoTrans)
 if(length(toptail)>0) pheno1[,toptail] = apply(pheno1[,toptail,drop=F],2,makeTopTail,strsplit(phenoTrans[toptail[1]],"_")[[1]][2])
  thresh = grep("^thresh",phenoTrans)
  if(length(thresh)>0) pheno1[,thresh] = apply(pheno1[,thresh,drop=F],2,makeThresh,strsplit(phenoTrans[thresh[1]],"_")[[1]][1:2])
 exps = grep("^exp",phenoTrans)
  print(exps)
  if(length(exps)>0) pheno1[,exps] = apply(pheno1[,exps,drop=F],2,exp10)
save.image()
 #facts = grep("^factor",phenoTrans)
 #pheno2 = pheno1
 #factInds = as.numeric(as.matrix(data.frame(strsplit(phenoTrans[facts],"_")))[2,])
 #for(ik in 1:length(factInds)) pheno1[,facts[ik]] =  makeFactor(pheno2[,facts[ik],drop=F],dimnames(pheno2)[[2]][facts[ik]])[,factInds[ik]]

dimnames(pheno1)[[2]] = dimnames(pheno)[[2]]
}


#print(stratMatrix1[1:10,])


##append the variance column and plate column
if(!is.null(varianceIsCov)){
if(varianceIsCov & length(var_id)>0) pheno_cov1 = cbind(pheno_cov1,as.numeric(variance))
if(varianceIsCov & length(plate_ind)>0) pheno_cov1 = cbind(pheno_cov1,plate)
if(!varianceIsCov & length(var_id)>0) pheno_resid1 = cbind(pheno_resid1,as.numeric(variance))
if(!varianceIsCov & length(plate_ind)>0) pheno_resid1 = cbind(pheno_resid1,plate)
}

##GETGENOCOV #DO NOT REMOVE LINE, USED BY JRI
#IN THIS SECTION, JRI PROVIDES genotmp matrix
#ITERATES OVER ALL genocov
if(!is.null(genotmp)){
  if(genoCovIsRes) pheno_resid1 = cbind(pheno_resid1,fixNonNumeric(genotmp)) else pheno_cov1 = cbind(pheno_cov1,fixNonNumeric(genotmp))
}
##ENDGETGENOCOV ##DO NOT REMOVE LINE, USED BY JRI



if(dim(pheno_cov1)[2]==0) phenoCovna = rep(0,dim(samples)[1]) else phenoCovna = is.na(apply(pheno_cov1,1,sum))
if(dim(pheno_resid1)[2]==0) phenoResna = rep(0,dim(samples)[1]) else phenoResna = is.na(apply(pheno_resid1,1,sum))




families = rep(0, dim(pheno1)[2])


##calculate which rows to include based on NAs in pheno, and in cov matrices
 inclM = as.matrix(apply(pheno1,2,getIncl, phenoCovna | phenoResna | varna))
 

#figure out whether to use binomial or gaussian
for(i in 1:(dim(pheno1)[2])) families[i] = getFamily(pheno1[,i],inclM[,i])


#if(inverseRegress) families = rep("gaussian",dim(pheno1)[2])

binom = families=="binomial"

##calculate indices of cases or controls
caseInd = apply(pheno1,2,iscase)





###FOLLOWING LINE IS TO MAKE 1 2 phenotypes 0 1
for(i in 1:(dim(pheno1)[2])) if(binom[i]) pheno1[,i] = pheno1[,i]-min(pheno1[,i],na.rm=TRUE)

if(dim(pheno_cov1)[2]>0)for(i in 1:(dim(pheno_cov1)[2])) pheno_cov1[,i] = centralise(pheno_cov1[,i])

# if(getFamily(pheno_cov1[,i],inclM[,i])=="binomial") pheno_cov1[,i]= pheno_cov1[,i]-min(pheno_cov1[,i],na.rm=TRUE)
#print(todo[binom])
#print(pheno1[1:10,])

#nagained = length(which(is.na(pheno1[,1]))) - length(which(is.na(pheno[,1])))
#print(paste("extra NA through mismatch of ids ",nagained,"of", dim(pheno1)[1]))

#calculate which indicies of geno file to consider
inds = c()
for(i in 1:length(ids)) inds = c(inds,grep(ids[i],geno_header))
ids = geno_header[inds]


geno_header = geno_header[inds]

alleleCols = grep("Allele[12]",geno_header)
gCols = grep("GType",geno_header)
if(length(alleleCols)>0){
  geno_header = c("geno",geno_header[!alleleCols])
}


if(!is.null(geno)) {
 geno = geno[,inds,drop=F]
 if(length(alleleCols)>0) geno = mergeAleleCols(geno,alleleCols)
 if(length(gCols)>0) geno = fixGenoCols(geno, gCols)

 maxg = apply(geno,2,max,na.rm=T)
 ming = apply(geno,2,min,na.rm=T)
 print(maxg)
 print(ming)
 
}

nonint = maxg%%1 > 0
nonint1 = ming%%1 > 0
maxg[nonint] = floor(maxg[nonint])+1
ming[nonint1] = floor(ming[nonint1])

print(maxg)
print(ming)

spl = seq(min(ming)+1,(max(maxg)),step) - 0.5


#make an array for the results



typenames = c("pheno","incl","caseInd","offsets")
phenN = phenNames
if(!is.null(phenoTrans)) phenN = apply(rbind(phenN,phenoTrans),2,paste,collapse=".")
arraynames = list(phenN, typenames,indiv )
if(debug)save.image()
###NOTE FOLLOWING LINE COULD ALSO BE PREMERGE IF WE WANTED TO QUANTILE
###NORMALISE EACH DATASET SEPARATELY
phendata = getArray(arraynames)
phendata[,1,] = t(pheno1)
phendata[,2,] = t(inclM)
phendata[,3,] = t(caseInd)



genoweights =getArray(list(geno_header,c("geno","weights","include"),samples[,1]))
#print(dim(geno))
#print(dim(genoweights))

if(!is.null(geno)) genoweights[,1,] = t(geno)

genoweights[,2,] =  t( matrix(1,nrow = dim(samples)[1],ncol = length(geno_header)))




##make the data object to store everything important
datanme = list("pheno"=phendata,"pheno_cov" = pheno_cov1, "pheno_res" = pheno_resid1, "geno" = genoweights, "strat"=stratMatrix1, "alleleCols" = alleleCols,"gCols" = gCols, "ind" = inds)

datanme2 = datanme
if(length(geno_header)==1) multiGen = FALSE
if(expandData){
 len = length(indiv)
 datanme2 =expandDat(datanme,length(spl)+1)
 # if(debug)save.image()
 if(!is.null(weights))datanme2$geno[,2,] = weights 

 for(i in 1:(length(spl)+1)) datanme2$geno[,1,((i-1)*len +1):(i*len)] = rep((i-1),length(geno_header))
 datanme2$geno[,3,] = matrix(TRUE,ncol = length(geno_header), nrow = length(indiv)*(length(spl)+1))
}




##FOLLOWING SECTION IS CALLED AFTER MERGING DIFFERENT DATASETS
##POSTMERGE DO NOT REMOVE COMMENT: JRI USES THIS TO KNOW WHEN TO RUN POSTMERGE
##calculate the offset


if(length(index)>0) for(i in 1:(length(index))) datanme$pheno[i,4,] = calcOffset(datanme$pheno[i,,],families[i], datanme$pheno_res)
if(debug) save.image()
#print(offsets)
 if(expandData){
   if(length(index)>0) for(i in 1:(length(index))) datanme2$pheno[i,4,] = calcOffset(datanme2$pheno[i,,],families[i], datanme2$pheno_res)
 }
if(!expandData) datanme2 = datanme


fams = levels(as.factor(families))
dimcounts = c(dim(datanme$strat)[2],length(index),length(geno_header),2)
if(length(index_strat)>0) dimcounts[1] = dimcounts[1]+1
if(multiPhen) dimcounts[2] = dimcounts[2]+1
if(multiGen) dimcounts[3] = dimcounts[3]+1
res = array(NA,dim=dimcounts)
#res2 = array(NA,dim=dimcounts)


phenNames1 = dimnames(pheno1)[[2]]
phenSubInds = 1:length(phenNames1)

#REPEAT   DO NOT REMOVE COMMENT JRI USES THIS TO KNOW WHERE TO FIND CODE TO RUN EVERY ITERATION
###TO DO EVERY LOOP
##JRI will provide a matrix geno 
##also rsid
##THIS PART MUST DEFINE res, countsControl and countsCase,maf,hwe_control, hwe_case  (to print out)
geno1 = datanme$geno




geno1[,3,] = !apply(as.matrix(geno1[,1,]),c(1,2),is.na)
geno2 = geno1

if(expandData) geno2 = datanme2$geno




##for OLR



maf = rep(NA,dim(geno1)[1])
for(k in 1:length(geno_header)) maf[k] = getMaf(as.numeric(geno1[k,1,]),baselevel[k],2)



#CONDITIONONMAF  ###DONT REMOVE THIS LINE, and must calculate maf before this line


if(multiPhen) out_inds = c(phenSubInds,dimcounts[2]) else out_inds = phenSubInds
if(debug)  save.image()
#print(geno2[1,1,1:10])
for(j in 1:(dim(datanme$strat)[2])){
    inclu = datanme2$strat[,j]
    res[j,out_inds,,] =  pvFunctMult(geno2,datanme2$pheno[phenSubInds,,,drop=F],families,inclu,datanme2$pheno_cov,rescale,pvFunct,funct)
    res[j,-out_inds,,] = NA
}
 



##meta-analysis
#print(length(index_strat))
if(length(index_strat)>0) res[dim(datanme$strat)[2]+1,,,] =apply(apply(res[1:(dim(datanme$strat)[2]-1),,,,drop=F],c(2,3),metares),c(3,1),t)


#EXTRALEVEL_1  ###DONT REMOVE THIS LINE
dimcounts = c(dim(datanme$strat)[2],length(families),length(geno_header),(length(spl)+1))

countsCase = array(NA,dim=dimcounts)
countsControl = array(NA,dim=dimcounts)
hwe_control = array(NA,dim=dimcounts[1:3])
hwe_case = array(NA,dim=dimcounts[1:3])
maf_control = array(NA,dim=dimcounts[1:3])
maf_case = array(NA,dim=dimcounts[1:3])





#print(dimcounts[1])


for(j in 1:(dimcounts[1])){
 stra = datanme$strat[,j]
 for(k in 1:length(families)){
  #  print(paste(j,k))
    inclu = datanme$pheno[k,2,] & stra
    incluCont = inclu &  !datanme$pheno[k,3,]
    incluCase = inclu &  datanme$pheno[k,3,]
    if(length(which(incluCont))>0)  countsControl[j,k,,] = t(apply(geno1[,1,incluCont,drop=F],1,getHist,spl))
    if(length(which(incluCont))>0)  countsCase[j,k,,] = t(apply(geno1[,1,incluCase,drop=F],1,getHist,spl))
    hwe_control[j,k,] =  apply(countsControl[j,k,,,drop=F],3,hwepall)
    hwe_case[j,k,] =  apply(countsCase[j,k,,,drop=F],3,hwepall)
    maf_control[j,k,] =  apply(countsControl[j,k,,,drop=F],3,mafCount)
    maf_case[j,k,] =  apply(countsCase[j,k,,,drop=F],3,mafCount)
 }
}

#EXTRALEVEL_2    ###DONT REMOVE THIS LINE


#uncomment following section if you want to print out the ids of cases and controls (subject to pvalue threshold
#specified above
#idsCase = array(NA,dim=dimcounts)
#idsControl = array(NA,dim=dimcounts)
#for(j in (dimcounts[1])){
# stra = datanme$strat[,j]
# for(k in 1:length(families)){ 
#    inclu = datanme$pheno[k,2,] & stra
#    incluCont = inclu &  !datanme$pheno[k,3,]
#    incluCase = inclu &  datanme$pheno[k,3,]
#        if(length(which(incluCase))>0) idsCase[j,k,,] = t(extractIdsMat(geno1[,1,,drop=F],spl,incluCase))
#        if(length(which(incluCont))>0) 	idsControl[j,k,,] = t(extractIdsMat(geno1[,1,,drop=F],spl,incluCont))
# }
#}

if(debug) save.image()
if(debug) print(res[dim(res)[1],,,])

