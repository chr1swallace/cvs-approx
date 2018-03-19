
library(magrittr)
library(data.table)
length(files <- list.files("/mrc-bsu/scratch/cew54/abf",full=TRUE))
data <- lapply(files, function(f) {
    x <- readRDS(f)
    x$f <- basename(f)
    x
})  %>% rbindlist()
head(data)
data[,minp:=pmax(minp1,minp2)]

## how many?
nsim <- nrow(unique(data,by=c("f","n0","n1","n2")))
message("number of completed simulations: ",nsim)

## our predicted intercept
data[,off12:=0.5*(m1 + m2) * log(n2+n1+n0)]
data[,off1:=0.5*(m1 ) * log(n1+n0)]
data[,off2:=0.5*(m2) * log(n2+n0)]
data[,int:=-off12 + off1 + off2]
data[,bfpred:=bfsum+int]
## with(data,plot(bfsum,bf12)); abline(0,1,col="red")
## with(data,plot(bfpred,bf12)); abline(0,1,col="red")
data[,minp:=pmax(minp,1e-40)]
data <- data[order(minp),]

head(data[,cor(bfsum,bf12),by=c("minp")])
tail(data[,cor(bfsum,bf12),by=c("minp")])
head(data[,cor(bfpred,bf12),by=c("minp")])
tail(data[,cor(bfpred,bf12),by=c("minp")])

f <- function(m) {
    cf <- coefficients(m)
    r2 <- summary(m)$r.squared
    list(int=cf[1],s=cf[2],r2=r2)
}
res.c <- data[,f(lm(bf12 ~ bfsum)),by=c("n0","p1","f","pred","minp","n1","n2")]
res.p <- data[,f(lm(bf12 ~ bfpred)),by=c("n0","p1","f","pred","minp","n1","n2")]

library(ggplot2)
library(cowplot)
library(randomFunctions)
nicecow <- function(p) randomFunctions::nicecow(p) + theme(legend.position="bottom")


if(FALSE) {
    tmp <- data[bf12!=0,.(bf12=bf12,err0=bf12-bfsum, err1=bf12-bfpred, err2=bf12-(bfpred-log(pi)/2))]
    tmp[,bagg:=round(bf12,1)]
    length(unique(tmp$bagg))
    dim(tmp)
    tmp[,n:=.N,by=bagg]
    s <- sample(1:nrow(tmp),10000,prob=1/tmp$n)
    summary(tmp[s,]$bf12)
    summary(tmp$bf12)

    dim(tmp)
    m <- melt(tmp[s,],c("bf12","bagg","n"))
    head(m)
    ## ggplot(m,aes(x=bf12,y=value)) + geom_hex() + facet_wrap(~variable)

    p <- ggplot(m,aes(x=bf12,y=value)) + geom_hex() + geom_density2d(col="red",alpha=0.5)  + facet_wrap(~variable) +
      ggtitle("Error in approximation") + labs(x="Multnomial ABF",y="Multinom ABF - approx") +
      geom_hline(yintercept = 0) + scale_y_continuous(limits=c(-1.6,0.8)) 
    nicecow(p)
}

    res.c[,approx:="sum"]
    res.p[,approx:="sumimproved"]
if(FALSE) {
    res <- rbind(res.c,res.p)

p1 <- ggplot(res,aes(x=-log10(minp),y=s)) + geom_point(size=.3,alpha=0.5,col="grey") + geom_smooth(aes(col=factor(p1)),se=FALSE) +
  geom_hline(yintercept = 1) + ggtitle("Slope") + facet_grid(p1 ~ approx) 
p2 <- ggplot(res,aes(x=-log10(minp),y=r2)) + geom_point(size=0.3,alpha=0.5,col="grey") + geom_smooth(aes(col=factor(p1)),se=FALSE) +
  geom_hline(yintercept = 1) + ggtitle("Rsq")+ facet_grid(p1 ~ approx)
plots <- list(p1,p2)  %>%  lapply(., function(p)
    randomFunctions::nicecow(p + #xlim(0,30) +
                             theme(legend.position="none", axis.title.y=element_blank())))
p <- plot_grid(plotlist=plots)
p


p11 <- ggplot(res.c,aes(x=-log10(minp),y=s)) + geom_point(size=3,alpha=0.5,col="grey") + geom_smooth(aes(col=factor(p1)),se=FALSE) +
  geom_hline(yintercept = 1) + ggtitle("Slope") + facet_grid(p1 ~ .) 
p12 <- ggplot(res.p,aes(x=-log10(minp),y=s)) + geom_point(size=3,alpha=0.5,col="grey") + geom_smooth(aes(col=factor(p1)),se=FALSE) +
  geom_hline(yintercept = 1) + ggtitle("Slope") + facet_grid(p1 ~ .) 
p21 <- ggplot(res.c,aes(x=-log10(minp),y=r2)) + geom_point(size=0.3,alpha=0.5,col="grey") + geom_smooth(aes(col=factor(p1)),se=FALSE) +
  geom_hline(yintercept = 1) + ggtitle("Rsq")+ facet_grid(p1 ~ .)
p22 <- ggplot(res.p,aes(x=-log10(minp),y=r2)) + geom_point(size=0.3,alpha=0.5,col="grey") + geom_smooth(aes(col=factor(p1)),se=FALSE) +
  geom_hline(yintercept = 1) + ggtitle("Rsq")+ facet_grid(p1 ~ .)
plots <- list(p11,p12,p21,p22)  %>%  lapply(., function(p)
    randomFunctions::nicecow(p + #xlim(0,30) +
                             theme(legend.position="none", axis.title.y=element_blank())))
p <- plot_grid(plotlist=plots)
p
}

res.p[,minp:=pmax(minp,1e-20)]
p12 <- ggplot(res.p,aes(x=-log10(minp),y=s)) + geom_point(size=1,alpha=0.5,col="grey") + geom_smooth(aes(col=factor(p1)),se=FALSE) +
  geom_hline(yintercept = 1) + ggtitle("Slope") + facet_grid(p1 ~ .) 
p22 <- ggplot(res.p,aes(x=-log10(minp),y=r2)) + geom_point(size=1,alpha=0.5,col="grey") + geom_smooth(aes(col=factor(p1)),se=FALSE) +
  geom_hline(yintercept = 1) + ggtitle("Rsq")+ facet_grid(p1 ~ .)
plots <- list(p12,p22)  %>%  lapply(., function(p)
    randomFunctions::nicecow(p + #xlim(0,30) +
                             theme(legend.position="none", axis.title.y=element_blank())))
p <- plot_grid(plotlist=plots)
p


ggsave(file="/scratch/wallace/IL2RA/approx-bf.png",plot=p,height=8,width=8)

