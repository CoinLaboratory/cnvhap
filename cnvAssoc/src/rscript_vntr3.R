run=TRUE
#doabs=TRUE
maf_thresh =0.01
preCalc = TRUE
thresh = 0.06
overlapThresh = 0
splitAdultGender = FALSE
getPlots=TRUE
cohort_inds = c(1,2,4)
plotToFile = TRUE
secondLevel=TRUE
lambda = 0
vntrname = c("VNTRA","VNTRB")
codefirst = FALSE  ##if true then it searches through all genos for a given model first, before trying diff model
offset_noInfo=NA
startpos = 0
do_enet=0
 todo =   c("child", "adult", "adultbC")
if(!run){
 allowNested = FALSE
 numRep=1
 phen_ind =3
 train_ind = 1
 stepind = "1"
 do_enet=0
 outdir<-paste("results",numRep,phen_ind,stepind,todo[train_ind],sep="_")
 dir.create(outdir)
 gen_inds = "1"
 code_inds = "1" 
 qt=FALSE
 central = TRUE
}

restrictToControls=getPlots & phen_ind>2
if(qt){
  phencolsmax = c(1,2,3)
}else{
 phencolsmax = c(1,2,3)
}
allcodes_full = cbind(c(0,1,2),c(0,0,1),c(0,1,1),c(0,1,0)) 
allcodes = as.matrix(allcodes_full[,as.numeric(strsplit(code_inds,',')[[1]])])
steps = as.numeric(strsplit(stepind,',')[[1]])
codesstring = apply(allcodes,2,paste,collapse="_")
codesstring1 = as.numeric(as.matrix(as.data.frame(strsplit(code_inds,",")))[,1])

geno_inds = as.numeric(strsplit(gen_inds,',')[[1]])
#allcodes = as.matrix(c(0,1,2))
library("sn")
library("mnormt")
library("elasticnet")

dir.create(outdir)
restrictCNV<-FALSE
watson_recoding = c(0,0,1)
saveTable = FALSE

options(contrasts=c("contr.treatment","contr.treatment"))

#### the input files

vntrb_file = c("vntrb_children.txt","vntrb_adult.txt")
vntra_file = c("vntra_children.txt","vntra_adult.txt")
if(qt){
 phen_file = c("pheno_child_qt.txt","pheno_adult_qt.txt","pheno_desir_qt.txt")
}else{
 phen_file = c("Child_pheno.txt","Adult1_pheno.txt","Adult1_pheno.txt")
}
cnv_file = c("cnv12003PP1_child","cnv12003PP1_adult")
cnv_file1 = c("cnv_child_all","cnv_adult_only","cnv_adult")
watson_file = c("watson_child","watson_adult")
#cnv_adult restricted to probes on diab dataset as well
#"pheno_adult.txt")


tokeep = 1
dotable=F

#################functions

# plots an x y1 y2 using left and right axes for the different y's
# rescales y2 to fit in the same space as y1
TwoVarPlot <- function(xvar, y1var, y2var, labels, noincs = 5, marks=c(1,2), legpos, leglabs, title="",cols = 1:2)
{
 # rescale to fit on same axis
 scaledy2var <- (y2var - min(y2var,na.rm=T)) / (max(y2var,na.rm=T) - min(y2var, na.rm=T))
 scaledy2var <- (scaledy2var * (max(y1var,na.rm=T) - min(y1var,na.rm=T))) + min(y1var,na.rm=T)

 # plot it up and add the points
 plot(xvar, y1var, xlab=labels[1], ylab="", axes=F, pch=marks[1],
 main=title,col=cols[1])
 points(xvar, scaledy2var, pch=marks[2],col=cols[2])

 # make up some labels and positions
 y1labs <- round(seq(min(y1var,na.rm=T), max(y1var,na.rm=T), length=noincs),2)

 # convert these to the y2 axis scaling
 y2labs <- (y1labs - min(y1var,na.rm=T)) / (max(y1var,na.rm=T) - min(y1var,na.rm=T))
 y2labs <- (y2labs * (max(y2var,na.rm=T) - min(y2var,na.rm=T))) + min(y2var,na.rm=T)
 y2labs <- round(y2labs, 2)

 axis(1)
 axis(2, at=y1labs, labels=y1labs)
 axis(4, at=y1labs, labels=y2labs)
 mtext(labels[3], side=4, line=2)
 mtext(labels[2], side=2, line=2)
 box()

 #legend(legpos[1], legpos[2], legend=leglabs, pch=marks, bty="n")
}



plotRes<-function(results, ind=1, main = NULL){
 plot(dimnames(results)[[1]],-log10(results[,ind]),main=main)
}

plotRes1<-function(results, results1,ind=1,be=0,fi=Inf){
 vec1=as.numeric(dimnames(results$r)[[1]])
 vec2=as.numeric(dimnames(results1$r)[[1]])
 m = match(vec2,vec1)
 inds = !is.na(m) 
 vec11 = vec1[m[inds]]
 vec22 = vec2[inds]
 st = which(vec11>=be)[1]
 ve = which(vec11<=fi)
 end =ve[length(ve)]
 ind1 = m[inds][st:end]
 ind2 = which(inds)[st:end]
 TwoVarPlot(dimnames(results$r)[[1]][ind1],-log10(results$r[ind1,ind]),-log10(results1$r[ind2,ind]),labels=c("length of vntr",results$name, results1$name),title=dimnames(results$r)[[2]][ind])
}


getInRange<-function(genvec, min, max){
 genvec >min & genvec <= max 
}

getOdds<-function(t){
 print(t[,,1])
 print(t[,,2])
 t1n = rbind(t[,,1],apply(t[,,1],2,sum))
 t1n = cbind(t1n,apply(t1n,1,sum))
 t2n = rbind(t[,,2],apply(t[,,2],2,sum))
 t2n = cbind(t2n,apply(t2n,1,sum))
 res = t2n/(t2n+t1n)
 dn = dimnames(t)
 dimnames(res) = list(c(dn[[1]],"sum"), c(dn[[2]],"sum"))
 res
}

getSig<-function(t){
 dimt = dim(t)
 if(length(dimt)>2){
   t2 = t[,,1]+t[,,2]
   ratio = t2/sum(t2)
   t1_exp = ratio*sum(t[,,1])
   t2_exp = ratio*sum(t[,,2])
   df = length(as.vector(t2))-1
   v = c(as.vector((t1_exp-t[,,1])^2/t1_exp),as.vector((t2_exp-t[,,2])^2/t2_exp))
   su = sum(v,na.rm = T)
   pchisq(su,df = df,lower.tail=T)
  }else{
    # print(dimt)
   	if(dimt[1]>=2 & dimt[2]>=2){
			prop.trend.test(t[,1],(t[,1]+t[,2]))$p.value
		#fisher.test(t)$p.value
   	}else{
    	      1.0
   	}
 }
}

innerProduct<-function(vec1,vec2){
  sum(vec1 > vec2[1] & vec1 <= vec2[2]) 
}

getVectorGLInner<-function(mat1,vec2){
 apply(mat1,1,innerProduct,vec2=as.numeric(vec2))
}

getVectorGL1<-function(mat1,mat2){
  apply(mat2,1,getVectorGLInner,mat1 = mat1)
}

getVectorGL22<-function(vlist1,j, pre = preCalc){
 vlist = as.numeric(vlist1)
 ind = vlist[3]
 step_ind = length(steps)
  print(ind)
 ind1 = which(rowhead[[step_ind]][[ind]]==vlist[1])
 ind2 = which(rowhead[[step_ind]][[ind]]==vlist[2])
 mat = allallvars[[step_ind]][[j]][[ind]]

if(pre){
   top1 = getVariable(mat[[1]],ind1,ind2)
   top2 = getVariable(mat[[2]],ind1,ind2)
   res =  getMatrix(top1,top2,allcodes[,vlist[4]])[,1]
 }else{
   geno = recode(vars[[1]][,ind1]+vars[[2]][,ind1],allcodes[,vlist[4]],0:2)
   geno1 = recode(vars[[1]][,ind2]+vars[[2]][,ind2],allcodes[,vlist[4]],0:2)
   res = geno-geno1
 }

##old 
## mat = allallvars[[step_ind]][[j]][[ind]][[vlist[4]]]
## res = mat[,ind2] - mat[,ind1]
 res
}

getVectorGL22All<-function(vlist,j){
 as.matrix(apply(vlist,2,getVectorGL22,j))
}

getVectorGL3<-function(datas,vlist){
   lev = levels(as.factor(vlist[3,]))
   li = matrix(NA,dim(datas$pheno)[1],dim(vlist)[2])
   for(i in as.numeric(lev)){
     inds = vlist[3,]==i
     mat = vlist[,inds]
     if(length(which(inds))==1) mat = as.matrix(mat)
     
     li[,inds] = getVectorGL1(datas[[i]][,2:3],t(mat))
    
   }
    dimnames(li) = list(NULL, apply(as.matrix(vlist[1:3,]), 2,paste,collapse="_"))

	res = recodeMats(li,allcodes)
    res1 = getY1Inner1(res,0,dimnames(li)[[2]])
   res1
   li
}

odds<-function(vec){
   c(vec,vec[2]/sum(vec))
}
overlap<-function(v,v1){
  vmin = min(v[1:2])
  vmax = max(v[1:2])
 # if(vmin==Inf & vmax==Inf){
 #   res =TRUE
#  }else{
# v1 = exclList[[i]]
  v1min = min(v1[1:2])
  v1max = max(v1[1:2])
  minmax = min(v1max,vmax)
  maxmin = max(v1min, vmin)
  leng = min(c(v1max - v1min,vmax - vmin))
  overlap = minmax - maxmin
  res = FALSE
 # print("overl")
 # print(c(v,v1))
  if(is.na(overlap)){
    res = TRUE
  }else if(overlap>overlapThresh){
     res = TRUE
    if(leng==overlap & allowNested){
      res = FALSE
      if(v1max==vmax | v1min==vmin){
        res  = TRUE
      }
    }
  }
 # }
  res
}
##geno_ind -> rowhead end index -> code_index

getRowHead<-function(step){
 rowhead = list();
 data = alldata[[1]]
 for(j in geno_inds){
 #resallj = list()
   dats = data[[j]][,2:3]
 level = as.numeric(levels(as.factor(c(dats[,1],dats[,2]))))
 rowhead[[j]] = c(level[ seq(1,length(level),step)],Inf)

  
 }	
 rowhead
}

getAllAllVars<-function(rowhead){
 res = list();
 for(i in 1:length(alldata)){
   res[[i]] = getAllVars(rowhead, alldata[[i]],geno_inds)
 }
 names(res) = names(alldata)
 res
}


getAllAllVars1<-function(rowhead){
 res = list();
 for(i in 1:length(alldata)){
   res[[i]] = getAllVars1(rowhead, alldata[[i]],geno_inds)
 }
 names(res) = names(alldata)
 res
}


getAllVars<-function(rowhead, data, geno_ind){
resall = list()
for(j in geno_ind){
 #resallj = list()
  
  rowheadj = rowhead[[j]]
 	
   len = length(rowheadj)
   i = len
 #for(i in 1:length(rowheadj)){
   matrix = cbind(rowheadj,rowheadj[i])  
   swap_ind = matrix[,1]> matrix[,2]
   matrix[swap_ind,] = matrix[swap_ind,c(2,1)] 
   vars1 = getVectorGL1(data[[j]][,2:3],matrix)
    dimnames(vars1) = list(NULL, rowheadj)

   vars = recodeMats(vars1,allcodes)   
  # resallj[[1+ len-i]] = vars   
 #} 
 resall[[j]] = vars  
#resallj
}
names(resall) = geno_ind
resall
}



getGt<-function(matj, num){
num>matj[1]
}

getAllVars1<-function(rowhead, data, geno_ind){
resall = list()
for(j in geno_ind){
  rowheadj = rowhead[[j]]	
   len = length(rowheadj)
   i = len
   matrix = cbind(rowheadj,rowheadj[i])  
   swap_ind = matrix[,1]> matrix[,2]
   matrix[swap_ind,] = matrix[swap_ind,c(2,1)]
   vars = list(apply(matrix,1,getGt,data[[j]][,2]),apply(matrix,1,getGt,data[[j]][,3]))  
   dimnames(vars[[1]]) = list(NULL, rowheadj)
   dimnames(vars[[2]]) = list(NULL, rowheadj)
   resall[[j]] = vars
}
names(resall) = geno_ind
resall
}


###data is global variable
getResults<-function(phenvec, len, colhead,  rowheadj,rowindex, geno_ind = 1,
  covar=NULL, tosubtract = rep(0,length(phenvec)), inds = 1:dim(allcodes)[2]){ 
 #  geno = data[[geno_ind]]
     offset = getOffsetL(phenvec,covar,len) 
   nme = paste(colhead,geno_ind,sep="_")
   nmes = names(data)
  
   vars = allvars[[geno_ind]][inds]
   results = getSigListMat(tosubtract,phenvec, vars,len, offset,rownames = rowheadj, rowindex = rowindex, inds = inds)
   toreturn = list("results"=results,"name"=paste(data$name,names(data)[geno_ind],sep="_"),bound=colhead)
   toreturn
}

readData<-function(phe_file, vn_file_list,phen_n="Obesity",watson_recoding = NULL){
	phen = read.table(phe_file,header=T,sep="\t")
        if(restrictToControls){
           phen = phen[phen[,2]==1,]
        }
        print(phen[1:3,])
       print(phen[,2])
      phen = phen[,phencolsmax]
	ind = c(1,which(names(phen)==phen_n))
	pheno = phen[,ind]
	nocol = dim(pheno)[2]
      resu = list()
      for(i in 1:length(vn_file_list)){
        vn_file = vn_file_list[[i]]
        nmei = names(vn_file_list)[i]
        if(nmei=="watson"){
		geno = procWatson(read.table(vn_file,header=TRUE,as.is=TRUE))
            if(!is.null(watson_recoding)) {
			geno[,2] =recode(geno[,2],watson_recoding)
            }
	  }else{
	    geno = read.table(vn_file,as.is=TRUE)
	    geno[,2] = round(as.numeric(geno[,2]))
	    geno[,3] = round(as.numeric(geno[,3]))
	    names(geno) = c("PATIENT","allele1","allele2")
        }
        mat = match(phen[,1],geno[,1])
        resu[[i]] = geno[mat,]
        names(resu)[i] = nmei
      }
      len = length(resu)
      resu[[len+1]] = phen
      names(resu)[len+1] = "pheno"
       resu[[len+2]] = phen_n
      names(resu)[len+2] = "name"
    
     
        if(length(resu$pheno)>2){
         resu$covar = as.matrix(resu$pheno[3:length(resu$pheno)] )
         dimnames(resu$covar)[[2]] = paste("covar",1:(length(resu$pheno)-2),sep="_")
          resu$pheno = resu$pheno[1:2]
        }
       
     
	resu
}

getOverlapExclusion<-function(excl,rowheadj, bnd,geno_ind){
 # print(paste("excl",geno_ind))
 # print(excl)
   v = rep(FALSE, length(rowheadj))
  if(dim(excl)[2]>0){
  excl1 = as.matrix(excl[,excl[3,]==geno_ind])
  if(dim(excl1)[2]>0){  
 if(bnd==Inf){
  mat = cbind(rowheadj,rowheadj)
}else{
  mat = cbind(rowheadj,bnd)
}
  ind = mat[,1]> mat[,2]
  mat[ind,] = mat[ind,2:1]
 
  for(i in 1: length(v)){
     v[i] = max( apply(excl1,2,overlap,mat[i,]))
  }
  }
  }
 !v
}
 
##data is global variable
train<-function(phenvec,len,covar=NULL,excl = data.frame(row.names = c("start","end","geno_ind","code_ind"))
,geno_ind = 1, inds = 1:dim(allcodes)[2]){
      
      rowheadj =  rowhead[[geno_ind]]
      tosubtract = rep(0,length(phenvec))
      if(allowNested){
        rowheadj1 = rep(TRUE, length(rowheadj)) 
      }else{
        rowheadj1 = getOverlapExclusion(excl,rowheadj,Inf,geno_ind)
      }

#getOverlapExclusion(excl,rowheadj,Inf, geno_ind)
    if(length(which(rowheadj1))>0){
       #   print(paste("hh",colhead[1]))

      res = getResults(phenvec,len, Inf, rowheadj, rowheadj1,geno_ind = geno_ind, covar = covar, tosubtract = tosubtract,inds)
        
      if(secondLevel){
            colhead =  as.numeric(findMax(res$res)[1,])   
            mat = allvars[[geno_ind]][[colhead[2]]]
		tosubtract1 = mat[,which(dimnames(mat)[[2]]==paste(colhead[1],colhead[2],sep="_"))]
             rowheadj1 = getOverlapExclusion(excl,rowheadj,colhead[1],geno_ind)
         #  if(length(which(rowheadj1))==0){
         #         rowheadj1 = c(getOverlapExclusion(excl,rowheadj[1:(length(rowheadj)-1)],Inf,geno_ind),FALSE)
         #         colhead[1] = Inf
         #  }
           if(length(which(rowheadj1))>0){
          #   print(paste("h",colhead[1]))
            # print(rowheadj1)
      	 res = getResults(phenvec, len,colhead[1],rowheadj,rowheadj1, covar = covar,geno_ind = geno_ind,tosubtract = tosubtract1, inds =colhead[2] )
            }else{
                res = NULL
          
            }
      }
      }else{
      res = NULL
    } 
	res
}

trainMultiple<-function(len,phenvec, covar = NULL,thresh = 0.05,geno_ind = c(1),targetNum=100){
 resultsall = list()

 excl = data.frame(row.names = c("start","end","geno_ind","code_ind","step_ind"))
 k=1 
  first = 1:length(codesstring1)
  second = geno_ind
  for(iks in second){ 
   pv = 0.0
   while(pv<thresh & length(excl)<targetNum){
    print(excl)
     if(length(excl)>0){
                          
		 covar1 =   cbind(getVectorGL22All(as.matrix(excl),train_ind),covar)
      }else{
 		covar1 = covar
      }
      minp =1.0
      max = rep(thresh+1,4)
      for(ijs in first){
        ij = ijs
 #     for(iks in second){
         ik = iks

         max21 = rep(0,2)
         max22 = rep(0,2)
         for(il in 1:length(steps)){
         print(paste("ij ik il", ij,ik,il))
          vars = allvars[[il]][[ik]]
          rowh =  rowhead[[il]][[ik]]
          if(il==1){
           rowindexL = rep(TRUE, length(rowh))
           rowindexU = rep(TRUE, length(rowh))
          }else{
	   rowindexL = rep(FALSE, length(rowh))
           rowindexU = rep(FALSE, length(rowh))
        #   print(max21)
        #   print(max22)
           rowindexL[(which(rowh == max21[1])[1]) :(which(rowh == max21[2])[1]) ] = TRUE
	   rowindexU[(which(rowh == max22[1])[1]) :(which(rowh == max22[2])[1]) ] = TRUE
          }
	  offset = getOffsetL(phenvec,covar1,len)
          res1 = getSigListMatMatFor(phenvec,vars, len, offset= offset,inds = ij,
			rownames = rowh, rowindexL = rowindexL,
	      rowindexU = rowindexU,sigfun=getSigListInner,excl=excl,geno_ind = ik)[[1]]
          max1 = as.numeric(findMax(res1,1)[1,])
          if(max1[2]<max1[1]) max1[1:2] = max1[2:1]
          inds2 = c(which(rowh == max1[1]),which(rowh == max1[2]))
          max21 = c(rowh[max(1,inds2[1]-1)],rowh[min(length(rowh),inds2[1]+1)])
          max22 = c(rowh[max(1,inds2[2]-1)],rowh[min(length(rowh),inds2[2]+1)])
	 }
         if(!is.null(res1)){    
          max1 = as.numeric(findMax(res1,1)[1,])
          if(max1[2]<max1[1]) max1[1:2] = max1[2:1]
          max1 = c(max1[1:2], ij,max1[3])
         if(!is.na(max1[4])){ 
         if(max1[4]<minp & max1[4]>=0){
           minp = max1[4]
           max = max1
           res = res1
           ind1 = ik
           ind2 = max1[3]
          }
         }
#        }
      }
     } 
     pv = max[4]
      print(paste("pv is ",pv))
     if(pv<thresh || dim(excl)[2]==0 ){
      v = c(as.numeric(max[1:2]),ind1,ind2,length(steps))
      v[1:2] = v[order(v[1:2])]
      print(v)
      print(names(data))
      excl[[k]] = v
      resultsall[[k]] = res
      nme = paste(v[1],v[2],v[3],v[4],sep="_")   
      mat = length(which(names(resultsall)==nme))
      if(mat >0){
          pv =1
       }    
      names(resultsall)[k] = nme
       k = k+1 
     }
   
   }
 #  }  
 }
 list("res" = resultsall,"nme"=names(resultsall))
}

merg<-function(datachilda,datachildb,nme){
	m = match(datachilda$g[,1], datachildb$g[,1])
	data = cbind(datachilda$g[,],datachildb$g[m,2:3])
	list("geno" = datachilda$g,"geno1" = datachildb$g[m,2:3],pheno = datachilda$p,name=nme)
}

join<-function(datac, dataa){
 res = list()
 for(i in 1:(length(datac))){
    print(i)
    nmei = names(datac)[i]
    print(nmei)
    if(nmei=="name" ){
       res[[i]] = "both"
    }else if(nmei=="rsid"){
 	 res[[i]] = datac$rsid
    }else{
     nmechild = names(datac[[i]])
     nmeadult = names(dataa[[i]])
     if(is.null(nmechild)){
       nmechild = dimnames(datac[[i]])[[2]]
        nmeadult = dimnames(dataa[[i]])[[2]]

     }
     if(length(nmechild)==length(nmeadult) && length(which(nmechild==nmeadult))){
	  mat = rbind(datac[[i]],dataa[[i]])
     } else{
        adult_ind = match(names(datac[[i]]),names(dataa[[i]]))
        na_ind = is.na(adult_ind)
        adult_ind[na_ind] = 1
        adult1 = dataa[[i]][adult_ind]
        adult1[na_ind] = NA
        names(adult1) = nmechild
        mat = rbind(datac[[i]],adult1)
     }
     res[[i]] = mat
   }
 } 
 names(res) = names(datac)
 res
}


findMax<-function(res,m=1){
	vec = as.vector(res)
	dimt = dim(res)
	order = order(vec)[1:m]
	res1 = matrix(nrow=m,ncol=3)
	nme = dimnames(res)
	for(i in 1:m){
	  row = order[i]%%dimt[1]
        if(row==0) row = dimt[1]
	  col = 1+(order[i] - row)/dimt[1]
	  res1[i,] = (c(nme[[1]][row], nme[[2]][col], res[row,col]))
	}
	#swap = res1[,1]>res1[,2]
	#res1[swap,c(1,2)] = res1[swap,c(2,1)]
      res1
}

as.factor.mat<-function(mat){
apply(mat,2,as.factor)
}

#datachild = mergeIn(datachild,cnvs[[1]], rsids = default_rsids)
calcEigs<-function(LRR){
colmn=apply(LRR,2,mean,na.rm=TRUE)
colsd=apply(LRR,2,sd,na.rm=TRUE)
sLRR=LRR

for (i in 1:dim(LRR)[2]) {
  sLRR[,i]=(LRR[,i]-colmn[i])/colsd[i]
  }
svd = svd(sLRR)
svd$u[,1]
}

mergeIn<-function(data, cnv,rsids=NULL){
	names(cnv)[[1]][1] = "PATIENT"
	mat = match(as.vector(data$pheno$PATIENT),as.vector(cnv$PATIENT))
	data$pheno = data$pheno[,1:2]
	data$cnv = cnv[mat,]
      leng = dim(data$pheno)[1]
      hasNull = FALSE
	if(is.null(rsids)){
        	rsids = getAllAssocs(data,1e-5)
        }
	data$rsid=rsids
      #  pvals = rep(0,length(rsids))
       
        for(i in 1:length(rsids)){
	        ind = which(names(data$cnv)==rsids[i])
              if(length(ind)==0){
			data$pheno[[2+i]] = rep(NA, leng)
			hasNull = TRUE
              }else{
			data$pheno[[2+i]] = data$cnv[,ind]
			
              }
          
	#	names(data$pheno)[2+i]="LRR"
#	        pvals[i] = assoc(
	}
       len = 2+length(rsids)
     names(data$pheno)[3:len] = rsids
    if(hasNull){
       offs = rep(NA,leng)
	}else{
      offs =  getOffset(data$pheno[3:len],data$pheno[[2]])
      }
      data$pheno[[len+1]] = offs
      names(data$pheno)[len+1] = "combined"
      
	data
}

transToPCA<-function(data, tospl, pca_name){
              eigsall = list()
            for(kk in 1:(dim(tospl)[1])){
                   ve = tospl[kk,1]:tospl[kk,2]
                   ve1 = ve +2
                   mat = data$pheno[,ve1]
                   inds = apply(is.na(mat),1,sum)==0
                   eigsa = rep(NA,length(inds))
                  
                  
                   eigsa[inds] =  calcEigs(as.matrix(mat[inds,]))
                   eigsall[[kk]] = eigsa
                   
               ##  summary(  glm(data$pheno[inds,2] ~ eigsa[inds], family="binomial"))
            }
            data$pheno[3:(2+length(eigsall))] = eigsall
          
          
        names(data$pheno)[3:(2+length(pca_nme))] = pca_nme
  data
}

getAllAssocs<-function(data,thresh=1e-5){
  dimt = dim(data$cnv)
  offset = rep(0,dimt[1])
  results = matrix(0,dimt[2],3)
  k=1
  minv =0
  rsid="none"
  len = 0
  while(minv<thresh && len==0 ){
  	vec  = assocAll(data,offset)
#phen =data$cnv$cnv12003PP1 )
	minv = min(vec[which(names(vec)!="cnv12003p2")])
  	rsid =  names(which(vec==minv))[1]
	
	len = length(which(results[,1]==rsid))
	if(len==0){
    	
        if(length(vec)==1){
		 ind = 2
             rsid = names(data$cnv)[2]
	  }else{
	      ind = which(names(data$cnv)==rsid)
	  }
        results[k,] = c(rsid,ind,minv)
        cnvvec = data$cnv[,as.numeric(results[1:k,2])]
	  if(k==1) cnvvec = as.matrix(cnvvec)
 	  offset = getOffset(cnvvec,data$pheno[[2]])
     	  k = k+1
	}
  }
  print(results[1:(k-1),])
  results[1:(k-1),1]
}




getLM<-function(phenvec,matr,len = length(levels(as.factor(phenvec)))){

	 ind = !is.na(phenvec) & apply(apply(matr,2,is.na),1,sum)==0
       
	  if(len>=3){
		family="gaussian"
  	  }else{
		family="binomial"
  	  }
	lm = glm(phenvec[ind] ~ matr[ind,],family=family)
	list("lm" =lm,"ind"=ind)
}


getOffsetL<-function(phenvec,list,len = length(levels(as.factor(phenvec)))){
	
 if(!is.null(list) && dim(list)[2]>0 ){
      offset = rep(offset_noInfo,length(phenvec))
	lms = getLM(phenvec,list,len)
	offset[lms$ind]=lms$lm$lin
  }else{
 offset = rep(0,length(phenvec))

} 
  offset
}

getOffset<-function(cnvvec,phenvec,offset=NULL){
   ind = !is.na(phenvec) & apply(apply(cnvvec,2,is.na),1,sum)==0
   offset1 = rep(0,length(ind))   
    if(qt){
    fam = "gaussian"
  }else{
     fam = "binomial"
   }

lm = glm(phenvec[ind]~as.matrix(cnvvec[ind,]),family=fam,offset = offset)
   offset1[ind] = lm$lin
   offset1
}

assoc<-function(cnvvec, phenvec,offset=NULL, fam="binomial"){
  ind = !is.na(phenvec) & !is.na(cnvvec)
 
  lm = glm(as.numeric(phenvec[ind])~as.numeric(cnvvec[ind]),family=fam,offset = offset[ind])
  summ =  summary(lm)

# coeff =   summary(glm(as.numeric(cnvvec[ind])~as.numeric(phenvec[ind])))$coeff
  coeff = summ$coeff 
dimt = dim(coeff)
      		if(dimt[1]>=2 & dimt[2]>=2){ 
	p = coeff[2,4]
}else{
 p =1.0
}
p
}


assocAll<-function(data,offset = NULL, phen = data$pheno[[2]]){
 len = length(levels(as.factor(phen)))
 if(len<=2){
   family = "binomial"
 }else{
    family = "gaussian"
 }

mat = data$cnv[,2:(dim(data$cnv)[2])]
if(dim(data$cnv)[2]==2){
 mat = as.matrix(mat)
}
 vec = apply(mat,2,assoc,phen,offset,family)
vec
}


procWatson<-function(wts){
 #for(j in 1:length(watson)){
    #wts = watson[[j]]
    vec = as.factor(wts[,2])
	l = levels(vec)
	res = rep(0,length(vec))
	v = c(0,1,2,NA)
	for(i in 1:length(l)){
	   res[vec==l[i]] = v[i]
	}
    wts[,2] =res
  ##  watson[[j]] = wts    
 #}
 #watson
 wts
}


getSigListInner<-function(var,pheno,family,offset){
  #print(paste(length(var), length(pheno), length(offset)))
  ind = !(is.na(pheno) | is.na(var))
  
 #print(tosubtract)
  x = var[ind]
  maf = (length(which(x==1)) + 2*length(which(x==2)))/(2*length(x))
    #-tosubtract[ind]
  ##if(doabs) x = abs(x)
  pv = 1.0
  if(maf>maf_thresh & (1-maf)>maf_thresh){
   lm = glm(pheno[ind] ~ x,family=family,offset=offset[ind])
   coeff = summary(lm)$coeff[,4]
   pv = min(coeff[2:length(coeff)])
  }
  pv
}

getAvgInner<-function(var,pheno,family,offset, tosubtract=NULL){
  #print(paste(length(var), length(pheno), length(offset)))
   ma  = max(var - tosubtract,na.rm=T)
  mi  = min(var - tosubtract,na.rm=T)

  ind = !(is.na(pheno) | is.na(var))
 #print(tosubtract)
  cases = pheno[ind]==1
  controls = pheno[ind] == 0
 # print(length(which(cases)))
#print(length(which(controls)))

  x = var[ind]
   #-tosubtract[ind]
  ##if(doabs) x = abs(x)
if(plotdiff){
   res = mean(x[controls]) - mean(x[cases])
}else{
  if(ma>0){
    res = mean(x[controls])
  }else if(mi<0){
   res = mean(x[cases])
  }else{
    res = 0
  }
}
res
}


  #results = getSigListMat(tosubtract,phenvec, vars,len, offset,rownames = rowheadj, rowindex = rowindex, inds = inds)
  
getVariable<-function(vars, ind1,ind2){
  vars1 = vars[,ind1]
  vars2 = vars[,ind2]
  #abs(vars2-vars1)
(vars1 & !vars2) | (vars2 & !vars1)
}


getMatrix<-function(top1,top2,codes){
###NOTE NEED TO CHANGE THIS BACK IF WE ADD IN CODES
apply(as.matrix(top1+top2),2,recode,codes,0:2)
#as.matrix(top1+top2)
#top1+top2
}


#tosubtract
getSigListMat<-function(index, pheno, vars, len = length(levels(as.factor(pheno))),offset = rep(0, length(pheno)),
rownames = dimnames(vars[[1]])[[2]],rowindex = rep(TRUE, length(rownames)), inds = 1:length(codestring1),
    sigfun=getSigListInner, pre = preCalc){
  if(len>=3){
	family="gaussian"
  } else{
	family="binomial"
  }

  resu  = matrix(NA, length(which(rowindex)),length(inds))
  if(pre){
   top1 = getVariable(vars[[1]],rowindex,index)
   top2 = getVariable(vars[[2]],rowindex,index)
   for(i in length(inds)){
     matr = getMatrix(top1,top2,allcodes[,inds[i]])
     resu[,i] = #getSigListInner(matr,pheno,family,offset) 
        apply(matr,2,getSigListInner,pheno,family,offset)   
   }
  }else{
    top1 = as.matrix(vars[[1]][,rowindex]+vars[[2]][,rowindex])
    for(i in length(inds)){
      geno = apply(top1,2,recode,allcodes[,inds[i]],0:2)
      geno1 = recode(vars[[1]][,index]+vars[[2]][,index],allcodes[,inds[i]],0:2)
      matr = geno-geno1
      resu[,i] = apply(matr,2,getSigListInner,pheno,family,offset)
    }
  }

 ###  dimnames(resu) = list(rownames[rowindex],codesstring1[inds[i]])

####OLD   
##  for(i in 1:length(vars)){
     
##     resu[,i] = apply(as.matrix(vars[[i]][,rowindex]),2,sigfun,pheno,family,offset, tosubtract)   
##     resu[is.na(resu[,i]),i] = 1 
 ## }
   dimnames(resu) = list(rownames[rowindex],codesstring1[inds])
 
  resu  
}

getSigListMatMat<-function( pheno, vars,len = length(levels(as.factor(pheno))),offset = rep(0,length(pheno)),
  inds = 1:length(codestring1),rownames = dimnames(vars[[1]])[[2]], rowindexL = rep(TRUE, length(rownames)),
      rowindexU = rowindexL, sigfun=getSigListInner){
  print(length(pheno))
  resu  = list()
   rowindexL1 = which(rowindexL)

  for(i in 1:length(inds)){
     ##old
    matri = 
     apply(as.matrix(rowindexL1),1,getSigListMat,
     ## as.matrix(vars[[i]][,rowindexL]),2,getSigListMat, 
             pheno,vars,len,offset,rownames=rownames,
                   rowindex = rowindexU,inds = inds[i], sigfun = sigfun )  
     for(k in 1: length(rowindexL1)){
            matri[k,is.na(matri[k,])] = 1.0
      } 
      resu[[i]] =    matri
    dimnames(resu[[i]]) = list(rownames[rowindexU], rownames[rowindexL])
    
  }
   names(resu) =codesstring1[inds]
   resu  
}

getSigListMatMatFor<-function( pheno, vars,len = length(levels(as.factor(pheno))),offset = rep(0,length(pheno)),
  inds = 1:length(codestring1),rownames = dimnames(vars[[1]])[[2]], rowindexL = rep(TRUE, length(rownames)),
      rowindexU = rowindexL, sigfun=getSigListInner, excl=NULL, geno_ind = 1){
  print(length(pheno))
  resu  = list()
  rowIndexL1 = which(rowindexL)
  rownamesin = rownames[rowindexL]
  dimy = length(rownamesin)
  for(i in 1:length(inds)){   
     #matr = as.matrix(vars[[i]][,rowindexL])
     #dimm = dim(matr)
    
     resmatr =  matrix(NA, length(rownames),dimy) 
     for(j in 1:dimy){
     #  print(j)
       rowindex = rowindexU
       if(!getPlots){
       #   print(rownames)
       #   print(rowindexL)
       #   print(rowindexU)    
       #   print(rownamesin[j])
          rowindex[rownames <= rownamesin[j]] = FALSE
       #   print(rowindex)
       }
       if(!is.null(excl) & length(excl)>0){
     		rowindex  = rowindex & getOverlapExclusion(excl,rownames,rownames[j],geno_ind)
       }
  ##should change here to make sure we don't calculate both i,j and j,i
      # print(rowIndexL1[j])
      # print(length(v
      if(length(which(rowindex))>0){
       resmatr[rowindex,j] = getSigListMat(rowIndexL1[j], pheno,vars,len,offset,rownames=rownames,
                   rowindex = rowindex,inds = inds[i], sigfun = sigfun )   
       }
     }
     dimnames(resmatr) = list(rownames, rownames[rowindexL])
     resu[[i]] = resmatr    
  }
   names(resu) =codesstring1[inds]
   resu  
}



seq1<-function(vec,l){
rgb = seq(vec[1],vec[2],length.out=l)
}

rgbfun<-function(rgb){
paste('#',paste(sprintf("%X",rgb),collapse=""),sep="")
}

plotContour<-function(m,main,range,type = 1,cont=T,numlev = 16){
print(paste("range", range))
#numlev1 = 10   #(floor(numlev/10))*10
#if(!cont & plotdiff)numlev1 = 50
#range[2] = (round(range[2]*numlev1))/numlev1
#print(paste("range", range))
#range = range*1.01
levels = seq(range[1],range[2],(range[2] - range[1])/numlev)
#seq(range[1], range[2], length.out = numlev)
#levels = pretty(range,numlev)
print(levels)
#start = 2/6 #(green)
start = 4/6 #(blue)
end = 1/6 #(yellow)
levels1 = seq(range[1],range[2],(range[2] - range[1])/10)


cols = colmap(numlev)
numext = 10
 #paste('#',paste(sprintf("%X",rgb),collapse=""),sep="")
cols1 = 
  rbind(as.numeric(paste('0x',(c(substr(cols[length(cols)],2,3), substr(cols[length(cols)],4,5),substr(cols[length(cols)],6,7))),sep="")),
    c(255,255,255))
vec1 = apply(as.matrix(apply(cols1,2,seq1,numext)),1,rgbfun)


leng = length(levels)
cols = c(cols[1:leng-1],vec1)
differ = levels[leng] - levels[leng-1]
levels = c(levels[1:leng-1],seq(levels[leng],levels[leng] + (numext-1)*differ, differ))
for(ik in 1:length(m)){
mk = m[[ik]]
nme = dimnames(mk)
x = as.numeric(nme[[1]])
y = as.numeric(nme[[2]])

x[length(x)] = x[length(x)-1]+1
y[length(y)] = y[length(y)-1]+1

#image(x, x, -log10(mk))
if(type==1){
filled.contour(x, y, -log10(mk),col = cols,
      plot.title = nme,levels = levels,main = main,title = main)
}else if(type==2){
  image(x,y,-log10(mk), col = cols)
 # x1 = rep(x,length(y))
 # y1 = 
 # text(x,y)
  if(cont) contour(x,y, -log10(mk),add=TRUE,levels = levels1)
}else{
 heatmap(-log10(mk),Rowv=NA, Colv=NA, sym=T,col = cols,main=main)
}
}
}

plotContours<-function(m,cohort_inds,vntr_type_ind, mode_ind,exp=F,type = 1,numlev = 16){
inds1 = cohort_inds
minpv = c(Inf,-Inf)
for(kj in inds1){
matkj = m[[kj]]
if(min(matkj[[1]])<minpv[1]){
  minpv[1] = min(matkj[[1]])
 }
 if(max(matkj[[1]]) > minpv[2]){
   minpv[2] = max(matkj[[1]])
 }
}
range = -log10(minpv[2:1])
if(range[1]>0) range[1] = 0
if(!exp) range[2] = rangemax
#range = 1.05*range
if(exp){
 nmeout = paste("cnts",plotdiff,sep="_")
}else{
 nmeout = "log10sig"
}
nmeout1 =  paste(names(allallvars[[1]])[cohort_inds],collapse=".")
nmeout = paste(nmeout,"_",nmeout1,"_",vntrname[vntr_type_ind],"_",codesstring1[mode_ind],".pdf",sep="")
if(plotToFile) pdf(file=paste(outdir,nmeout,sep="/"))
for(ik in inds1){
  if(!plotToFile & ik>1) dev.new()
   plotContour(m[[ik]], names(allallvars[[1]])[ik],range, type  = type,cont= !exp,numlev = numlev)
}
if(plotToFile) dev.off()
}

getContour<-function(allallvars_, rowhead_, vntr_type_ind, cohort_inds, 
        mode_inds, ranges = NULL, sigfun = getSigListInner, exp=F){
m = list()
j = vntr_type_ind
inds = mode_inds
inds1 = cohort_inds
rownames = rowhead_[[j]]
rowindex = rep(TRUE, length(rownames))
rowindexL = rowindex
rowindexU = rowindex
if(!is.null(ranges)){
  lowrange = ranges[1,]
  upperrange = ranges[2,]
  rowindexL = (rownames <= lowrange[2] & rownames >=lowrange[1]) 
   rowindexU = (rownames <= upperrange[2] & rownames >=upperrange[1]) 
}
for(kj in inds1){

 matkj = getSigListMatMat(alldata[[kj]]$pheno[,phen_ind], allallvars_[[kj]][[j]], inds = inds,rownames = rownames, 
              rowindexL = rowindexL, rowindexU = rowindexU,sigfun = sigfun)
 if(exp){
   matkj[[1]] = 10^(-1*abs(matkj[[1]]))
 }
 
 m[[kj]] =matkj
}
m
}


getTable<-function(phenvec,vars){
  results = matrix(0,dim(vars)[2],7)
 nmes = dimnames(vars)[[2]]
 dimnames(results) = list(
 nmes,
c("nme","cntrl0","cntrl1", "cntrl2", "case0", "case1", "case2")
 )
 len = length(levels(as.factor(phenvec)))
 ind1 = !is.na(phenvec)
 for(i in 1:(dim(vars)[2])){
    ind = ind1 & !is.na(vars[,i])
    vari = vars[ind,i]
    lev = as.numeric(levels(as.factor(vari[ind]))) 
    if(length(lev)<=3){
     t  =  table(list(phenvec[ind],vari))
    # print(i)
     dimt = dim(t)
     results[i,1] = nmes[i]
     ext = min(dimt[2],3)
     results[i,2:(1+ext)] = t[1,1:ext]
     results[i,5:(4+ext)] = t[2,1:ext]
    }
 }
 results
}

##alldata is global variable
getSigTableAll<-function(phenvec,vars, startpos = startpos){
 results = matrix(0,dim(vars)[2],6)
 #tosubtract = rep(0,length(phenvec))
 nmes = dimnames(vars)[[2]]
 dimnames(results) = list(
 nmes,
c("nme","pv","pv_model","var_expl_dev","var_explMZ","N")
 )
 len = length(levels(as.factor(phenvec)))
 prev_var = 0
prev_var1 = 0
 if(len>=3){
	family="gaussian"
  } else{
	family="binomial"
  }
for(i in 1:(dim(vars)[2])){
 
 if(i==1){
   covar1 = NULL
   
 }else if(i==2){
   covar1 = as.matrix(vars[,1:(i-1)])
 }else{
  covar1 = vars[,1:(i-1)]
 }
# print(len)
 offset = getOffsetL(phenvec, covar1,len)
#print(covar1[1:10,])
#print(offset[1:10])
 pv = getSigListInner(vars[,i],phenvec,family,offset)
 lms = getLM(phenvec,as.matrix(vars[,1:i]),len)
 hyper.lm2 = getLM(phenvec,as.matrix(vars[,i]),len)$lm

 ind = lms$ind
 hyper.lm =  lms$lm 
avg =  mean(phenvec[ind])
a = hyper.lm$lin
if(len<=2){
 a = exp(a)/(1+exp(a))
}
 x = (phenvec[ind] - a) %*% (phenvec[ind] - a)
 x1 = (phenvec[ind] - avg) %*% (phenvec[ind] - avg)
 varexpl1 =  (x1-x)/x1
 #print(hyper.lm$family)
if(i<=startpos || startpos==0 ){
 df = i
 hyper.lm1 = update(hyper.lm,~1,family=hyper.lm$family)
}else{
  df = i-startpos
  hyper.lm1 =getLM(phenvec[ind],as.matrix(vars[ind,1:startpos]),len)$lm
}
 llhFull <- logLik(hyper.lm) 
 llhNull <- logLik(hyper.lm1)
 gstat <- -2*(llhNull-llhFull)
 fitP <- pchisq(gstat[1],df,lower.tail=F)
 if(length(fitP)==0){
 fitP = 1.0
} 
 varexpl = (hyper.lm$null.deviance - hyper.lm$dev)/hyper.lm$null.dev
 Var = var(hyper.lm2$lin)
 
varexpl1 = Var/(Var + (pi^2)/3)
  varexpl = varexpl - prev_var
 # varexpl1 = varexpl1 - prev_var1
prev_var = varexpl+prev_var
prev_var1 = varexpl1+prev_var1
if(len>2){
 varexpl = cor(phenvec,vars[,i],use="pairwise.complete.obs")^2
}
# print(varexpl)
  vec = c(nmes[i],c(pv,fitP,varexpl,varexpl1,length(which(ind))))
 # print(vec)
 # print(fitP)
 results[i,] = vec
}
 lms = getLM(phenvec, vars,len)
 coeff = summary(lms$lm)$coeff
 coeff1 = matrix(NA,dim(results)[1],3)
 dimnames(coeff1) = list(dimnames(results)[[1]], c("beta","se","pv_t"))
 if(length(dimnames(vars)[[2]])+1 == dim(coeff)[1]){
    coeff1[1:length(dimnames(vars)[[2]]),] = coeff[2:dim(coeff)[1],c(1,2,4)]
 }else{
 inds2 = match(dimnames(vars)[[2]],as.matrix(data.frame(strsplit(dimnames(coeff)[[1]],']')))[2,])
 inds3 = which(!is.na(inds2))
 inds4 = which(is.na(inds2))
 #inds2 =match(as.matrix(data.frame(strsplit(dimnames(coeff)[[1]],']')))[2,],dimnames(vars)[[2]])
 

coeff1[inds3,]= coeff[inds2[inds3],c(1,2,4)]
le = length(inds4)
coeff1[inds4,1] = rep(0,le)
coeff1[inds4,2] = rep(Inf,le)
coeff1[inds4,3] = rep(1,le)
}
list("results"=cbind(results,coeff1),"lm" = lms)
}


getOrder<-function(nump,phenvec){
 orderall = list()
 order1 = 1:length(phenvec)
 for(i in 1:nump){
   if(i==1){
     phenvec1 = phenvec
     covar1 = covar
     order = order1
   }else{
     order = sample(order1, length(order1))
     phenvec1 = phenvec[order]
     covar1 = covar[order,]
   }
   orderall[[i]] = order
 }
 orderall
}

##data is global variable
getPermutedPv<-function(orderall,geno_ind = c(1), thresh = thresh,covar=NULL){
 tables = list()
 nump = length(orderall)
 result = matrix(0,nump,3)
 phenvec = data$pheno[,phen_ind]
 len = length(levels(as.factor(phenvec)))
 thresh1 = thresh
 target = 100
 resultsall = list()
 for(i in 1:nump){
 
  order = orderall[[i]]
   print(paste("HERE",i))
   print(order[1:10])
  phenvec1 = phenvec[order]
  #print(phenvec1)
  covar1 = covar[order,]
   if(do_enet>0){
	   resu = trainENet(phenvec1)
         resu$nme = names(resu$nme)
         inds = getInds(resu$nme)
	   resu$nme = resu$nme[inds[[1]]]
   }else{
        resu = trainMultiple(len,phenvec1,geno_ind = geno_ind,covar = covar1,thresh =thresh1,targetNum = target)
      
   }
   nmes = resu$nme
 
   matr = as.matrix(as.data.frame(strsplit(nmes,'_')))
 
   if(length(matr)>0){
    vars = getVectorGL22All(matr,train_ind)
    dimnames(vars) = list(NULL,resu$nme)
    covarvar  = cbind(covar1,vars)
   }else{
     covarvar = covar1
   }
   if(is.null(covar1)){
	dimcc=0
   }else{
      dimcc = dim(covar1)[2]
   }
      tt = getSigTableAll(phenvec1,covarvar, startpos = startpos)
    if(i==1){
     nme = tt$results[,7]
     target = dim(matr)[2]
     thresh1 = 2.0
     print(paste("matrix",target))
     print(matr)
    }
   t = tt$results
   tables[[i]] = t
   dimt = dim(t)
   leng = dimt[1]
   resultsall[[i]] = tt$results
   result[i,] = c(t[1+dimcc,2],t[leng,3],min(t[(dimcc+1):leng,2],na.rm=TRUE))
 }
res=list("simul" = result,"nme" = nme, resultsall=resultsall)
res
}

getPermPv<-function(perm){
 minl = 3
 resultsall = perm$resultsall
 matr = resultsall[[1]]
 nump = length(resultsall)
 len = length(resultsall)
 result = matrix(NA,len,dim(matr)[1])
 result1 = matrix(NA,len,dim(matr)[1])
 dimnames(result) = list(NULL,dimnames(resultsall[[1]])[[1]])
 dimnames(result1) = list(NULL,dimnames(resultsall[[1]])[[1]])
 for(i in 1:(dim(matr)[1])){
     for(j in 1:len){
       resj = resultsall[[j]]
       if(i<=(dim(resj)[1])){
	       result[j,i] = resj[i,3] 
               result1[j,i] = resj[i,2] 
       }
     }
 }
   resall = list(result1,result)
   pv = matrix(NA,dim(result)[2],length(resall))
  # print(result)
   if(nump>=minl){
   for(k in 1:(dim(pv)[1])){
    for(j in 1:(dim(pv)[2])){
      print(paste("resvec",j,k))
      resvec = as.numeric(resall[[j]][,k])
      resvec[resvec==0] = 2e-16
      resvec = resvec[!is.na(resvec)]
     
      print(resvec)
      len = length(resvec)
     
      if(len>=minl & var(resvec)>0.001){
        fit = cp.to.dp(sn.mle(y=log10(resvec[2:len]),plot.it = FALSE)$cp)
        pv[k,j] =  psn(as.numeric(log10(as.numeric(resvec[1]))),dp=fit)
      }
     }
   }
  }
  dimnames(pv) = list(dimnames(matr)[[1]],dimnames(matr)[[2]][2:3])
  perm$simulpv = resall[[1]]
  perm$simulpvmodel = resall[[2]]
  perm$pv = pv
  perm
}

getInds<-function(nmeall){
list(ind1 = grep("covar", nmeall,invert=T,value=F),
indcov = grep("covar", nmeall,invert=F,value=F))
}

#alldata is global var
getAllRes<-function(permnme, phenosToIncl = NULL){
    nmeall = names(permnme)
    inds_all = getInds(nmeall) 
    ind1 = inds_all[[1]]
    nme = strsplit(nmeall[ind1],'_') 
    indcov = inds_all[[2]]
    nmecov  = strsplit(nmeall[indcov],'_')  
    beta = permnme[ind1]
    betacov = permnme[indcov]

    if(length(nmecov)>0){
      cov_inds = as.numeric(as.matrix(as.data.frame(nmecov))[2,])
    }
    all_res = list()
      allvars = list()

    for(i in 1:length(alldata)){
      print(names(alldata)[i])
      allresj = list()
      allvarsj = list()

      datasi = alldata[[i]]
      vars = getVectorGL22All(as.matrix(as.data.frame(nme)),i)
      
       vars3 = vars
   ##getVectorGL3(datasi,as.matrix(as.data.frame(nme))) 
      varssum = vars%*%as.numeric(beta)
      dimnames(vars) = list(NULL,nmeall[ind1])   
      
      if(!is.null(cov_inds)){
         cov = as.matrix(allcovar[[i]][,cov_inds])
         dimnames(cov) = list(NULL,nmeall[indcov])
         covsum1 = as.matrix(cov[,1:(startpos-1)])%*%as.numeric(betacov[1:(startpos-1)])
         covsum2 = as.matrix(cov[,(startpos):length(betacov)])%*%as.numeric(betacov[(startpos):length(betacov)])
         dimnames(covsum2) = list(NULL,"covsum2")
         dimnames(covsum1) = list(NULL,"covsum1")
         vars = cbind(allcovar[[i]],vars)
         if(!is.null(phenosToIncl)){
            vars = cbind(vars,datasi$pheno[,phenosToIncl])
         }
         vars3 = cbind(allcovar1[[i]],vars3)
         dimnames(varssum) = list(NULL,"weighted")
         varssum = cbind(covsum1,covsum2,varssum)
      }
      phen_namesi = names(datasi$phen)[2:(dim(datasi$phen)[2])]
      for(j in 1:length(phen_namesi)){
           phenvec = datasi$pheno[,(j+1)]  
          tabular =  getTable(phenvec, vars3)
	    tt = getSigTableAll(phenvec,vars, startpos = startpos)
          tt1 = getSigTableAll(phenvec,varssum, startpos = 1)
          tabular1 = getTable(phenvec,varssum)
     #     tt2 = getSigTableAll1(phenvec,allcovar[[i]],vars3)
           allvarsj[[j]] = cbind(phenvec,vars)

          allresj[[j]] =cbind( rbind(tt1$res,tt$res), rbind(tabular1, tabular))
#,tt2$res)
      }
      names(allresj) = phen_namesi
       names(allvarsj) = phen_namesi

      all_res[[i]] = allresj   
      allvars[[i]] = allvarsj
    } 
    names(all_res) = names(alldata)
    names(allvars) = names(alldata)

    list(all_res = all_res,allvars = allvars)
}
writePermRes<-function(res,file){
 append=F
 dir.create(file)
 file1 = paste(file,"perms.txt",sep="/")
file2 = paste(file,"simulpv.txt",sep="/")
file3 = paste(file,"simulpvmodel.txt",sep="/")

 #results = perm$res
 write.table(t(as.matrix(res$nme)),file=file1,append=F,row.names=FALSE,quote=F,sep="\t")
 write("perm pvals",file=file1,append=T)
 write.table(res$pv,file=file1,append=F,quote=F,sep="\t",row.names=TRUE,col.names=TRUE) 
# write("simul_results",file=file2,append=F)
 write.table(res$simulpv,file=file2,append=F,row.names=FALSE,quote=F,sep="\t") 
  write.table(res$simulpvmodel,file=file3,append=F,row.names=FALSE,quote=F,sep="\t") 
}
writeRes<-function(results,file){
 append=F
 dir.create(file)
 len = length(results)
 for(i in 1:len){
    namei = paste(file,names(results)[i],sep="/")
    dir.create(namei)
    resj = results[[i]]
    nme = names(resj)
    for(j in 1:length(nme)){
       fileout = paste(namei,nme[j],sep="/")
       write.table(resj[[j]],file=fileout,append=F,row.names=FALSE,quote=F,sep="\t")
    }
 }
 
}

getY1Inner1<-function(allv,var_thresh,nme_base){
  dimt =  dim(allv[[1]])
 matr = matrix(0,dimt[1],dimt[2]*length(allv))
 names = rep("",dimt[2]*length(allv))
 for(i in 1:dimt[2]){
    st = length(allv)*(i-1) 
 #   print(paste("i",i,st))
    for(j in 1:length(allv)){
      matr[,(st+j)] = allv[[j]][,i]
      names[(st+j)] = paste(nme_base[i],j,sep="_")
    }
    
 }
 
 dimnames(matr) = list(NULL, names)
 matr
}


getY1Inner<-function(allv,var_thresh){
  dimt =  dim(allv[[1]])
 matr = matrix(0,dimt[1],dimt[2]*length(allv))
 names = rep("",dimt[2]*length(allv))
 for(i in 1:length(allv)){
    st = dimt[2]*(i-1) 
    print(paste("i",i,st))
    matr[,(st+1):(st+dimt[2])] = allv[[i]]
    names[(st+1):(st+dimt[2])] = dimnames(allv[[i]])[[2]]
 }
 var = apply(matr,2,var,na.rm=T)
 ma  = matr[,var>var_thresh]
 dimnames(ma) = list(NULL, names[var>var_thresh])
 ma
}
getY1<-function(allvars,j,var_thresh){
 allv = allvars[[j]]
 ma = getY1Inner(allv,var_thresh)
 dimnames(ma) = list(NULL, paste(dimnames(ma)[[2]],j,sep="_"))
 ma
}

 
prepareMat<-function(covar,var_thresh){
 m1 = getY1(allvars,1,var_thresh)
 m2 =  getY1(allvars,2,var_thresh)
 mat = cbind(m1,m2)
 matname = as.matrix(as.data.frame(strsplit(dimnames(mat)[[2]],'_')))
 matname = rbind(matname,"Inf")
 matname = matname[c(1,4,3,2),]
 dimnames(mat)[[2]] = apply(matname,2,paste,collapse="_")
 allmat = cbind(covar,mat)
 allmat
}
 


cntNonZero<-function(vec, thresh,inds){
sum(abs(vec[inds])>thresh)
}

trainENet<-function(phenvec1){
 ind = apply(apply(allmat,c(1,2),is.na),1,sum)==0 & !is.na(phenvec1)
 x=allmat[ind,]
 y = phenvec1[ind]
 print(dim(x))
 print(length(y))
 #en = enet(x, y,lambda = 0.5)
s = 0
max = 1
#while(s==0){
pdf(file=paste(outdir,"cv.pdf",sep="/"))
 en = cv.enet(x, y,K=10,lambda = lambda,s=seq(0,max,length=100),mode="fraction",trace=TRUE,max.steps=80,plot.it=TRUE)
# en = cv.enet(x, y,K=10,lambda = 0.5,s=1:50,mode="step")
 vec3 = which(en$cv == min(en$cv))
 s = en$s[vec3[length(vec3)]]
# max = max/50
#}
 en1 = enet(x, y,lambda = lambda)
 plot(en1)
 dev.off()
 beta = en1$beta
 noncovind = grep("covar", dimnames(beta)[[2]],invert=T,value=F)
 cntvec = apply(beta,1,cntNonZero,0.001, noncovind)
 vec2 = abs(en1$L1/max(en1$L1)-s)
 ind = which(vec2==min(vec2))[1]
 ind = max(ind, which(cntvec>0)[1])
 print(paste("ind is",ind))  
 ind1 =  beta[ind,]!=0
 print(which(ind1))
 v1 =  beta[ind,ind1]
 print(dimnames(beta)[[2]])
names(v1) = (dimnames(beta)[[2]])[ind1]

 list(nme=v1)
}

restrict<-function(cnvs, nmes){
  cnvs1 = list()
	for(i in 1:length(cnvs)){
		cnvs1[[i]] = cnvs[[i]][match(nmes,names(cnvs[[i]]))]
            names(cnvs1[[i]]) = nmes
	}	
names(cnvs1) = names(cnvs)
cnvs1
}

recode<-function(vec,code,lev){
 res = rep(NA,length(vec))
 for(i in 1:length(lev)){
   res[vec==lev[i]] = code[i]
 }
 res
}
########################
#

recodeMat<-function(mat,codes){
apply(mat,2,recode,codes,c(0,1,2))
}
recodeMats<-function(mat, allcodes){
  res = list()
  len = dim(allcodes)[2]
  for(i in 1:len){
     res[[i]] =  recodeMat(mat, allcodes[,i])
     dimnames(res[[i]]) = list(NULL,paste(dimnames(mat)[[2]],i,sep="_"))
  }
  #dimnames(res) = codesstring1
  res
}

centralise<-function(covar){
mean = mean(covar,na.rm=T)
naind = is.na(covar)
covar[!naind] = covar[!naind] - mean
covar[naind] = 0
covar
}

centraliseAll<-function(covar){
for(i in 1:dim(covar)[2]){
covar[,i] = centralise(covar[,i])
}
covar
}
prop<-function(vec){
ind = !is.na(vec)
sum(vec[ind])/length(vec[ind])
}
fishert<-function(vec,ph){
tab = table(vec, ph)
dimt = dim(tab)
pv = 1.0
if(dimt[1]>=2 & dimt[2]>=2){
 pv = fisher.test(tab)$p.v
}
pv
}
plotDistr<-function(){
type = "l"
dim = dim(allallvars[[1]][[1]][[1]])
nme = c("VNTRA","VNTRB")
for(j in 1:length(allallvars[[1]])){##j  vntra/b
 for(i in 1:2){ ##i child, adult, etc 
  
    pheno = alldata[[i]]$pheno[,2]
   case = alldata[[i]]$pheno[,2]==1
    control = alldata[[i]]$pheno[,2]==0
    x =  rowhead[[j]]  
    x[length(x)] = x[(length(x)-1)]+1
    mat = allallvars[[i]][[j]][[1]]
     
      ind = !is.na(pheno) & !is.na(mat[,1])
    fisher = apply(mat[ind,],2,fishert,pheno[ind])
   #plot(x,-log10(fisher))
    allv = apply(mat,2,prop)  
    main = paste(names(alldata)[i], nme[j],sep="_")
   
      matcase = mat[case,]
      matcontrol = mat[control,]
    if(FALSE){
     
      matcasek = apply(matcase,2,prop) - allv
      matcontrolk = apply(matcontrol,2,prop) - allv
      k1 = i
       if(i==1){
    
         TwoVarPlot(x,matcasek,-log10(fisher), lables = c("length VNTRA (bp)",ylab="cumulative fraction","-log10p"),title=main)
        # plot(x,matcasek,type=type,col=k1,lty=1,pty=1,xlab="length VNTRA (bp)",ylab="cumulative fraction",main=main)
       }else{
	   lines(x,matcasek,type=type,col=k1,lty=1,pty=1,xlab="length VNTRA (bp)",ylab="cumulative fraction",main=main)
       }
       lines(x,matcontrolk,type=type,col=k1,lty=2,pty=2,xlab="length VNTRA (bp)",ylab="cumulative fraction",main=main)
      
     }else{
        pdf(file = paste(outdir,"/",main,".pdf",sep=""))
    for(k in 0:2){
      k1 = k+1
      matcasek = apply(matcase==k,2,prop)
      matcontrolk = apply(matcontrol==k,2,prop) 
      if(k==0){
        TwoVarPlot(x,matcasek,-log10(fisher), labels = c("length VNTRA (bp)",ylab="cumulative fraction","-log10p"),title=main)
       
        #plot(x,matcasek,type=type,col=k1,lty=1,pty=1,xlab="length VNTRA (bp)",ylab="cumulative fraction",main=main)
      }else{
   	   lines(x,matcasek,type=type,col=k1,lty=1,pty=1,xlab="length VNTRA (bp)",ylab="cumulative fraction",main=main)
      }
       lines(x,matcontrolk,type=type,col=k1,lty=2,pty=2,xlab="length VNTRA (bp)",ylab="cumulative fraction",main=main)
     }
     
     dev.off() 
    }   
  }
}

}

getCovars<-function(data,codes, central){
#res1 =getOffset(as.matrix(data$watson[,2]),data$pheno[,2])
#res2 =getOffset(as.factor.mat(as.matrix(data$watson[,2])),data$pheno[,2])
#as.matrix(res2)
len = length(data$watson[,2])
dimt = dim(codes)
res = matrix(NA,len,dimt[2])
for(i in 1:dimt[2]){
  res[,i] = recode(as.numeric(data$watson[,2]),codes[,i],c(0,1,2))
}
le = dim(data$cov)[2]
codes1 = rbind("covar",(1+le):(dimt[2]+le))
dimnames(res) =  list(NULL,apply(codes1,2,paste,collapse="_"))

#apply(codes,2,paste,collapse="_"))
if(central){
 data$cov = centraliseAll(data$cov)
#resu1 = centraliseAll(resu1)
}
resu1 = cbind(data$cov,res)

resu1
}

splitBySex<-function(data, covar){
v = (as.numeric(levels(as.factor(covar[,1]))))
ind1 = abs(covar[,1]-min(v))<0.01
ind2 = abs(covar[,1]-max(v))<0.01
data1 = list()
data2 = list()
nme = c("vntra","vntrb","watson","pheno","covar","cnv")
nme1 = names(data)
for(i in 1:length(nme)){
  
  ind = grep(nme[i],nme1)
 
    data1[[i]] = data[[ind]][ind1,]
   data2[[i]] = data[[ind]][ind2,]


}
names(data1) = nme
names(data2) = nme
list(data1=data1,data2=data2,covar1 = covar[ind1,],covar2 = covar[ind2,])
}

if(run){
  
  cnvs = list(read.table(cnv_file1[1],header=T,as.is=TRUE),read.table(cnv_file1[2],header=T,as.is=TRUE),
                                                     read.table(cnv_file1[3],header=T,as.is=TRUE) )
  cnvs[[4]] = rbind(cnvs[[1]],cnvs[[2]])
  names(cnvs[[4]]) = names(cnvs[[1]])
  tokeep = 1
  datachild = readData(phen_file[1], 
       list("vntra" = vntra_file[1],"vntrb" = vntrb_file[1], "watson" = watson_file[1]),"child")    
  dataadult = readData(phen_file[2], 
       list("vntra" = vntra_file[2],"vntrb" = vntrb_file[2], "watson" = watson_file[2]),"adult")                                  
  dataadult_comb =readData(phen_file[3], 
       list("vntra" = vntra_file[2],"vntrb" = vntrb_file[2], "watson" = watson_file[2]),"adultComb")                                  
  codes = cbind(c(0,0,1)) 
  codes1 = cbind(c(0,1,2))
  covar_child = getCovars(datachild,codes, central)
  covar_adult = getCovars(dataadult,codes, central)
  
  covar_adultC = getCovars(dataadult_comb,codes, central)
  
  covar_child1 = getCovars(datachild,codes1, FALSE)
  covar_adult1 = getCovars(dataadult,codes1, FALSE)
  covar_adultC1 = getCovars(dataadult_comb,codes1, FALSE)

  #covar_both = getCovars(databoth,codes)
  startpos  = dim(covar_child)[2] 
   v = dimnames(cnvs[[1]])[[2]]
  default_rsids = v[grep("cnv12",v)] 
 tospl = rbind(c(1,3),c(3,length(default_rsids)),c(1,length(default_rsids)))
pca_nme = c("PCA VNTRA", "PCA VNTRB", "PCA BOTH")
  datachild = mergeIn(datachild,cnvs[[1]], rsids = default_rsids)
  dataadult = mergeIn(dataadult,cnvs[[2]],rsids = datachild$rsid)
# dataadult_comb = 
   
  if(phen_ind==2){
 databoth = join(datachild,mergeIn(dataadult_comb,cnvs[[3]],rsids = datachild$rsid))   
  } else{
 databoth = join(datachild,dataadult) 
  }
 
  dataadult_comb = mergeIn(dataadult_comb,cnvs[[3]]) 
  covar_both = getCovars(databoth,codes, central)
  covar_both1 = getCovars(databoth,codes1,FALSE)

 if(qt){
   covar_both = cbind(databoth$covar,covar_both)
 }

datachild = transToPCA(datachild ,tospl, pca_nme)
dataadult = transToPCA(dataadult, tospl, pca_nme)
databoth = transToPCA(databoth, tospl, pca_nme)

  alldata = list("child" = datachild,"adult"=dataadult,"adultbC" = dataadult_comb,"both" = databoth)
  allcovar= list(covar_child,covar_adult, covar_adultC, covar_both)
    allcovar1= list(covar_child1,covar_adult1, covar_adultC1, covar_both1)

  if(splitAdultGender){
     dataadultsplit = splitBySex(dataadult,covar_adult)

    alldata = list("child" = datachild,"adult1"=dataadultsplit$data1,"adult2" = dataadultsplit$data2,"both" = databoth)
    allcovar= list(covar_child,dataadultsplit$covar1, dataadultsplit$covar2, covar_both) 
 }
  
  
  

  rowhead = list()
  allallvars = list()
  allvars = list()

for(j in 1:length(steps)){
  rowhead[[j]] = getRowHead(steps[j])
  allallvars[[j]] = getAllAllVars1(rowhead[[j]]);
  allvars[[j]] = allallvars[[j]][[train_ind]] 
}
 
 

 
  
 

dir.create(outdir)
if(!getPlots){
  data = alldata[[train_ind]]
  covar = allcovar[[train_ind]]
  orderall = getOrder(numRep,data$pheno[,phen_ind])
  perm = getPermutedPv(orderall, geno_ind = geno_inds, thresh = thresh,covar = allcovar[[train_ind]])
  perm1 = getPermPv(perm) 
  perm = perm1


#res1 = getAllRes(permnme)
  res = getAllRes(perm$nme) 
obj = list(perm=perm,res =res)
  save(obj,file=paste(outdir,"robj.save",sep="/"))
 writeRes(res$all_res,outdir)
 writePermRes(perm,outdir)

 outdir1 = paste(outdir,"vars",sep="/")
dir.create(outdir1)
 writeRes(res$allvars,outdir1)
}else{
 colmap = topo.colors
rangemax = 3.5

for(mode_ind in codesstring1){
 for(vni in geno_inds){
    plotdiff=FALSE
  conts = getContour(allallvars[[1]],rowhead[[1]],vntr_type_ind = vni,cohort_inds =cohort_inds,
                         mode_inds = mode_ind, sigfun = getSigListInner)
  plotContours(conts,cohort_inds, vni,mode_ind,type =1,numlev =40)



  #conts = getContour(allallvars[[1]],rowhead[[1]],
  #      vntr_type_ind = vni,cohort_inds = cohort_inds, mode_inds = mode_ind,ranges = NULL,sigfun = getAvgInner,exp=T)# getSigListInner) #
  #plotContours(conts,cohort_inds, vni,mode_ind, exp=T,type = 1,numlev = 40)

  #plotdiff=TRUE
  #conts = getContour(allallvars[[1]],rowhead[[1]],
   #    vntr_type_ind = vni,cohort_inds = c(1,2,4), mode_inds = mode_ind,ranges = NULL,sigfun = getAvgInner,exp=T)
 # plotContours(conts, cohort_inds, vni,mode_ind, exp=T,type = 1,numlev = 40)# getSigListInner) #
 }
}
}




getPermNme<-function(l1){
l = c("covar_1","covar_2",l1)
#st,end,ind
#l = c("covar_1","covar_2",paste("671","741",ind, codesstring1,sep="_"))
beta = as.vector(as.factor(rnorm(length(l))))
names(beta) = l
beta
}

getTotalVarExpl<-function(res){
inds = which(dimnames(res$all_res$child$O)[[2]]=="var_explMZ")
vec1 = res$all_res$child$O[,inds]
vec2 = res$all_res$adult$O[,inds]
s1 = sum(as.numeric(vec1[5:length(vec1)]))
s2 = sum(as.numeric(vec2[5:length(vec1)]))
print(paste("child",s1))
print(paste("adult",s2))
}

getNewAssoc<-function(dats){
 resu = list()
for(i in 1:length(dats$allvars)){
mat = dats$allvars[[i]][[1]]
mat1 = as.data.frame(mat)
names(mat1) = gsub("^","x",names(mat1))
 form = paste("xphenvec~",
    paste(names(mat1[2:length(mat1)]),collapse = "+"),
 # "+",paste(names(mat1[4:length(mat1)]),collapse = "*"),
   sep=" ")
for(k in 4:(length(mat1)-1)){
  for(j in (k+1): length(mat1)){
     print(paste(j,k))
    form = paste(form, paste(names(mat1[c(k,j)]),collapse = "*"),sep="+")
  }
}
summ = summary(glm(form,data = mat1, family="binomial"))
resu[[i]] = summ$coeff
}
resu

}

getStEnd<-function(vec){
mat = as.matrix(as.data.frame(strsplit(vec,"_")))
v1 = as.numeric(mat[1:2,])
v1 = v1[order(v1)]
res = matrix(0,4,(length(v1)-1))
for(i in 1:(length(v1)-1)){
  res[1,i] = v1[i]
  res[2,i] = v1[i+1]
}
res[3,] = rep(mat[3,1],(length(v1)-1))
res[4,] = rep(mat[4,1],(length(v1)-1))
apply(res, 2,paste,collapse = "_")
}

getMultiTable<-function(dats,vec,flip=FALSE){
resu = list()
for(i in 1:length(dats$allvars)){
mat = dats$allvars[[i]][[1]]
ve = c(which(dimnames(mat)[[2]]==vec[1]),which(dimnames(mat)[[2]]==vec[2]))
if(flip) ve = ve[2:1]
mat = mat[,c(ve,1)]
mat = mat[apply(is.na(mat),1,sum)==0,]
tab = table(as.data.frame(mat))
resu[[i]] = apply(tab,1,oddsr,"")
}
resu
}

getMultiTable1<-function(dats,vec,flip=FALSE){
resu = list()
for(i in 1:length(dats$allvars)){
mat = dats$allvars[[i]][[1]]
ve = c(which(dimnames(mat)[[2]]==vec[1]),which(dimnames(mat)[[2]]==vec[2]))
if(flip) ve = ve[2:1]
mat = mat[,c(1,ve)]
mat = mat[apply(is.na(mat),1,sum)==0,]
tab = table(as.data.frame(mat))

case = tab[1,,]
control = tab[2,,]
p = case/(case+control)
odds = p/(1-p)
resu[[i]] = odds
}
resu
}

if(FALSE){
phen_ind = 2
adultnme = getAllRes(getPermNme(c("619_671_1_1","589_619_1_1","995_1033_2_1")), 3:4)
childnme = getAllRes(getPermNme(c("590_640_1_1","944_1022_2_1","1112_1127_2_1",
    "1073_1084_2_1","1099_1103_2_1","1075_1094_2_1")))


resad = getMultiTable(adultnme,  c("589_619_1_1","995_1033_2_1"),flip=TRUE)
reschild = getMultiTable1(childnme1,  c("590_640_1_1","976_1041.9_2_1"),flip=TRUE)
res1  = getNewAssoc(childnme1)
res2 = getNewAssoc(adultnme)



}

}
