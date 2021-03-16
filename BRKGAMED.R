BRKGAMED<-function(dataset,p=50,NG=50,pe=0.2,pm=0.7,pc=0.85,NGSM=0.35*NG,k=3)
{
#dataset = file with separate columns with semicolons 
#p = size of population
#NG = Number of generations
#pe = size of elite population
#pm = size of mutant population
#pc = crossover probability
#NGSM = generations without improvement
#k = number of clusters  
solution<-function(x,D,k)
{
 Med<-x[1:k]
 fobj<-Fobj(Med,D,k)
 G<-rep(0,nrow(D))
 for(i in 1:k) {G[fobj$C[[i]]]<-i}
 return(list(Me=Med,Fobj=fobj$Fobj,Gk=G))
}     
###############################################################
###########Cálculo Função Medoid###############################
###############################################################
Fobj<-function(medoids,D,nclust)
{ n<-nrow(D)
  Y<-1:n
  X<-setdiff(Y,medoids)
  Y[medoids]<-1:nclust
  Fm<-rep(0,nclust)
  Gm<-replicate(nclust,list(NULL))
  for(i in X) {Y[i]<-which.min(D[i,medoids])}
  for(j in 1:nclust) 
    {Gj<-which(Y==j)
     Fm[j]<-sum(D[medoids[j],Gj])
     Gm[[j]]<-Gj
    }
  return(list(Me=medoids,Fmed=Fm,Fobj=sum(Fm)/n,C=Gm))
}    

Update_medoids<-function(M,D,k)
{
S<-Fobj(M,D,k)
n<-nrow(D)
Me<-rep(0,k)
Fobj<-0
Fmed<-rep(0,k)
for(i in 1:k) 
 {g<-S$C[[i]]
  if (length(g)>1)
   {z=cluster::pam(D[g,g],1,diss=TRUE)
    Me[i]<-as.numeric(z$medoids)
    Fmed[i]<-z$objective[2]*z$clusinfo[1,1]
   } else {Fmed[i]<-0;Me[i]<-g}  
  } 
fobj<-sum(Fmed)/n
return(c(Me,Fmed,fobj))
}
  
###########################################################################
Generation<-function(p,k,D,n)
{
Me<-t(replicate(p,sample(n,k)))
#XO<-t(apply(Me,1,function(xm) Fobj(xm,D,k)))
#XO<-t(parSapply(cl=clust,XO,function(M) Update_medoids(M$Me,D,k)))
XO<-t(parApply(cl=clust,Me,1,function(M) Update_medoids(M,D,k)))
#XO<-t(sapply(XO,function(M) Update_medoids(M$Me,D,k)))
return(XO)
}  
###########################################################################

crossover<-function(c1,c2,D,k)
{
nm<-matrix(0,ncol=k,nrow=2*k^2)
li<-1
ls<-k
for(i in 1:k)
  {nm[li:ls,]<-cbind(c2,matrix(rep(c1[-i],k),ncol=k-1,byrow=TRUE))
   li<-ls+1
   ls<-ls+k
  }
for(i in 1:k)
  {nm[li:ls,]<-cbind(c1,matrix(rep(c2[-i],k),ncol=k-1,byrow=TRUE))  
   li<-ls+1
   ls<-ls+k
 }
cf<-which(apply(nm,1,function(x) any(duplicated(x)==TRUE))==FALSE)

if (length(cf)==0) 
  {nm<-rbind(c1,c2);pc=1}
else {nm<-nm[cf,];if(is.matrix(nm)==FALSE){nm<-t(as.matrix(nm));pc=1}}
px<-which(runif(nrow(nm))<=pc)
if (length(px)>0) {nm<-nm[px,]}
if (is.matrix(nm)==FALSE) {nm<-t(as.matrix(nm))}
Mec<-nm
Mec<-t(apply(Mec,1,function(x) Update_medoids(x,D,k)))
Mec<-Mec[which.min(Mec[,2*k+1]),]
return(Mec)  
}  



####################Main Program###########################################
B<-read.table(dataset,sep=";",dec=",") 
D<-as.matrix(dist(matrix(scale(B),ncol=ncol(B))))
n<-nrow(B)
tempo<-proc.time()
library(parallel)
nucleos<-detectCores(logical=F)
clust<-makeCluster(nucleos)
clusterEvalQ(clust, library(cluster))
FGLOBAL<-Inf
iter<-0
PE<-round(pe*p) #Size of elite set
PM<-round(pm*p) #Size of mutation set
PC<-p-PE-PM     #Size of crossover set 
fbest<-Inf
sm<-0
nb<-1
s<-Generation(p,k,D,n)
while((iter<NG) & (sm<NGSM))
  {iter<-iter+1  
   sm<-sm+1
   inx<-order(s[,2*k+1])
   fmin<-s[inx[1],2*k+1]
   if (fmin<fbest)
      {sbest<-s[inx[1],]
       fbest<-fmin
       cat("Generation ",iter," FOBJ = ",fbest,"\n")
       sm<-0
       ibest<-iter
      }
    CE<-inx[1:PE] #Elite solutions
    CNE<-inx[(PE+1):p] #NonElite Solutions
    sn<-s[CE,]
    smutants<-Generation(PM,k,D,n) ### Mutation
    sn<-rbind(sn,smutants) 
    ce<-cbind(sample(CE,PC,replace = TRUE),sample(CNE,PC,replace=TRUE))
  
    sc<-t(parApply(cl=clust,ce,1,function(y) crossover(s[y[1],1:k],s[y[2],1:k],D,k)))
    sn<-rbind(sn,sc)
    s<-sn
}
stopCluster(clust)
tempo<-(proc.time()-tempo)[3]
cat("Tempo em segundos ",tempo,"\n")
s=solution(sbest,D,k)

#Fobj = Value of objective function
#clustering = vector with the allocation of objects to the clusters
#Med = Medoids
#ibest = Generation in which the best solution was obtained
#cpu = BRKGA runtime

return(list(Fobj=s$Fobj,clustering=s$Gk,Med=sbest[1:k],ibest=ibest,cputime=tempo))
  
}




