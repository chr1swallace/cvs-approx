## qR.rb -r -y 1-500 -a tesla ./approx-generate-data.R
library(mvtnorm) 
library(magrittr) 

fdir <- "/home/ja628/causal_sim/results/IL2RA/BF"

##' .. content for \description{} (no empty lines) .. 
##' ##' .. content for \details{} .. 
##' @title 
##' @param pred degree of predictability for the phenotype. between 0 (completely random) and 1 (completely determined) 
##' @return ##' @author Chris Wallace 

simone <- function(pred) {    
    x <- replicate(nblocks,rmvnorm(n=nsamp,sigma=S),simplify=FALSE)    
    xpred <- lapply(x,function(xi) xi[,1]) %>% do.call("cbind",.) %>% rowSums(.)    
    x %<>% do.call("cbind",.)    
    y <- ifelse(runif(nsamp) > pred,                
                sample(c(0,1),nsamp,replace=TRUE),         
         ifelse(xpred>median(xpred),1,0))   
    
    as.data.frame(cbind(y=y,x=x)) 
} 

nsamp <- 4000
##' .. content for \description{} (no empty lines) ..
##'
##' @title simulate data for a pair of diseases
##' @param pred proportion of cases influenced by genetics
##' @param nsamp number of samples
##' @param p1 proportion of samples that are cases
##' @param nshare number of causal variants shared between traits
##' @return data.frame of case/control (y) and genotypes (x)
##' @author Chris Wallace
simpair <- function(pred,nsamp,p1,nshare=2) {    
    x <- replicate(nblocks,
                   rmvnorm(n=nsamp,sigma=S),
                   simplify=FALSE)    
    xpred1 <- x[[1]][,1] + x[[2]][,1] #lapply(x,function(xi) xi[,1]) %>% do.call("cbind",.) %>% rowSums(.)    
    xpred2 <- if(nshare==2) {
                  x[[1]][,1] + x[[2]][,1]
              } else if(nshare==1) {
                  x[[1]][,1] + x[[3]][,1]
              } else {
                  x[[4]][,1] + x[[3]][,1]
              }
    x %<>% do.call("cbind",.)    
    y1 <- ifelse(runif(nsamp) > (pred),
                sample(c(0,1),nsamp,replace=TRUE),         
         ifelse(xpred1>median(xpred1),1,0))    
    y2 <- ifelse(runif(nsamp) > (pred),
                sample(c(0,1),nsamp,replace=TRUE),         
         ifelse(xpred2>median(xpred2),1,0))    
    cases <- which(y1==1 & y2==1)
    y <- rep(0,nsamp)
    ## randomly reorder
    y[y1==1] <- 1
    y[y2==1] <- 2
    y[y1==1 & y2==1] <- sample(c(1,2),length(cases),prob=c(p1,1-p1),replace=TRUE)
    as.data.frame(cbind(y=y,x=x))
}

library(mlogitBMA)
library(BMA)
library(annotSnpStats)
library(GUESSFM)

library(Rcpp)
library(RcppArmadillo)

source("/scratch/wallace/IL2RA/myglib.R")

modlist.fn <- function(mT1,s) {
    modT1 <- vector("list",mT1)
    for(i in 1:mT1) {
        mci <- combn(letters[1:s],i,simplify=TRUE)  # all combinations of i msnps
        nci <- dim(mci)[2]
        modT1[[i]] <- matrix(0,nrow=nci,ncol=s)
        for(j in 1:nci)  modT1[[i]][j,match(mci[,j],myLetters)] <- 1 
    }
    modsT1 <- rep(0,s)
    for(i in 1:mT1) modsT1 <- rbind(modsT1,modT1[[i]])
    return(modsT1)
}
                                        #s=length(msnps) # number of snps considered in models for the 2 traits

format.mod.fn <- function(k,out) {
    #' called by T1mods.fn
    ind <- which(out[k,]==1)
    if(length(ind)>0) {
        mod <- paste(names(out[k,][ind]),sep="",collapse="%")
    } else {mod <- ""}
    return(mod)
}


mT1=2
nperblock <- 2 
nblocks <- 5 ## only blocks 1-2 have anything linked to causal
nsamp <- 2000 
S <- matrix(0.3,nperblock,nperblock) 
diag(S) <- 1 
s <- nperblock*nblocks
myLetters <- letters[1:26]
T1mods <- modlist.fn(mT1=mT1,s=s)
T2mods <- T1mods
nT1 <- dim(T1mods)[1]
nT2 <- dim(T2mods)[1]
T1modsrep <- matrix(rep(t(T1mods),nT2),ncol=ncol(T1mods),byrow=TRUE)
T2modsrep <- matrix(rep(T2mods,each=nT1),ncol=ncol(T2mods),byrow=FALSE)
T1T2mods <- cbind(T1modsrep,T2modsrep)
                                        # add column of 1's to the models for the trait2*effect variable
T1T2mods1 <- cbind(T1T2mods,1)
T1T2mods <- T1T2mods1
idx <- seq(1,((nblocks-1)*nperblock+1),by=nperblock)
nocausal <- rowSums(T1modsrep[, idx])==0 & rowSums(T2modsrep[,idx])==0
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param pred
##' @param nsamp total sample size
##' @param p1
##' @param nshare number of shared CVs
##' @return 
##' @author Chris Wallace
check.fn <- function(pred,nsamp,p1,nshare=2) {
    
    mult <- simpair(pred,nsamp,p1,nshare=nshare)
    y <- c("CONTROL","trait1","trait2")[ mult$y+1 ]
    
    mult1 <- data.frame(Y=y,mult[,-1])
    m1 <- mlogit2logit(Y ~ 1|. -Y,mult1,choices=c("CONTROL","trait1","trait2"),base.choice=1)
    mod1 <- glib.1(x=m1$data[,(4+s):(dim(m1$data)[2])],
                   y=m1$data$Y.star,
                   error="binomial", link="logit",models=T1T2mods)
    LOGABF12 <- mod1$bf$twologB10[,2]*0.5
    
    hi1 <- mult[mult$y !=1,]
    hi1$y <- hi1$y/2
    lo1 <- mult[mult$y !=2,]
    
    colnames(T1mods) <- names(mult[,-1])
    j.mat <- matrix(1:dim(T1mods)[1],ncol=1) 
    t1mod <- apply(j.mat,1,format.mod.fn,as.matrix(T1mods))
    t1size <- apply(T1mods,1,sum)
    T1mod <- data.frame(mod=t1mod,size=t1size,stringsAsFactors=FALSE)
    
    mod <- glib.1(x=hi1[,-1],y=hi1[,1],error="binomial", link="logit",models=T1mods)
    LOGABF2 <- mod$bf$twologB10[,2]*0.5
    mod <- glib.1(x=lo1[,-1],y=lo1[,1],error="binomial", link="logit",models=T1mods)
    LOGABF1 <- mod$bf$twologB10[,2]*0.5
    bf1bf2 <- cbind(bf12=LOGABF12, bf1=rep(LOGABF1,nT2), bf2=rep(LOGABF2,each=nT1))
    bf1bf2 <- cbind(bf1bf2,bfsum=bf1bf2[,"bf1"] + bf1bf2[,"bf2"])
    
    ## min p
    idx <- which(lo1[,1]==0)
    pv1 <- min(sapply(2:ncol(lo1), function(i) t.test(lo1[idx,i],lo1[-idx,i])$p.value))
    idx <- which(hi1[,1]==0)
    pv2 <- min(sapply(2:ncol(hi1), function(i) t.test(hi1[idx,i],hi1[-idx,i])$p.value))
    res <- cbind(as.data.frame(bf1bf2),pred=pred,nsamp=nsamp,p1=p1,minp1=pv1,
                 nshare=nshare,
          minp2=pv2,n1=sum(lo1[,1]),n2=sum(hi1[,1]),nocausal=nocausal,m1=rowSums(T1modsrep),m2=rowSums(T2modsrep))
}

resbak <- NULL

library(magrittr)
library(parallel)
options(mc.cores=14)
library(data.table)

res <- mclapply(seq(1000,3000,by=150), function(n) {
    message(n)
    p <- sample(seq(0.1,0.5,by=0.1),1)
    p1 <- sample(c(0.2,0.4,0.5,0.6,0.8),1)
    ns <- sample(0:2,1)
    tmp <- check.fn(pred=p,nsamp=n,p1=p1,nshare=ns)
    tmp$bf12 <- tmp$bf12-tmp$bf12[1]
    tmp$n0 <- n-tmp$n1-tmp$n2
    tmp$bf12.off <- tmp$bf12[1]
    tmp$bfdiff <- tmp$bf12 - tmp$bf12[1] - tmp$bfsum
    as.data.table(tmp)
    ## cat(".")
    ## tmp[,.(bfdiff=mean(bfdiff),
    ##        bfdiffsd=sd(bfdiff)),by=c("n0","n1","n2","nshare","m1","m2")]
})  %>% do.call("rbind",.)
res[,nsamp:=NULL]
head(res)


f <- tempfile(tmpdir = "/mrc-bsu/scratch/cew54/abf", fileext = ".rds")
saveRDS(res,f)
