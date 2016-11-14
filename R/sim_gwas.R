library(snpStats)
library(mvtnorm)
library(optparse)

## HARDCODED PARAMS
DATA.DIR<-'/home/ob219/scratch/bs_sim/support' # where sigma + support data is
ncp<-80 # non centrality parameter for computing mean (sqrt(ncp)) of MVN under alternative
n.perms<-20 # number of causal variants to simulate
types<-c('null','no.pir','pir','test') # feature spike ins to simulate

sim_alt<-function(s,gt,sigma,type=c('null','test','control','pir','no.pir'),ncp=80,n=1){
    ## this is under the null where there is no association
    if(type=='null'){
        ncp=0
    }
    ## else we can simulate diff scenarios
    else if(type=='no.pir'){
        caus.var<-sample(which(!b$support[['pir']]),1)
    }else{
        caus.var<-sample(which(b$support[[type]]),1)
    }
    idx.pir.snps<-which(s$pir)
    if(ncp>0){
        ## ok so compute r with the causal variant for all other SNPs we intend to simulate
        to.sim<-setdiff(idx.pir.snps,caus.var)
        r<-as.vector(ld(x=b$gt[,caus.var],y=b$gt[,idx.pir.snps],stats="R"))
        mvn.mean<-ceiling(sqrt(ncp))
        cv.pos<-b$support[caus.var,]$start
    }else{
        mvn.mean<-0
        r<-rep(0,nrow(sigma))
        cv.pos<-0
    }
    rd<-as.vector(rmvnorm(n,mean=r*mvn.mean,sigma=sigma,method="chol"))
    return(list(cv.pos=cv.pos,rd=rd,type=type))
}


##MAIN

## OPTION PROCESSING

option_list = list(
    make_option(c("-f","--file"),type="character",default = NULL, help = "Input RData object see - create_filtered_sigma.R"),
    make_option(c("-o","--out_dir"),type="character",default = NULL, help = "Directory to store result data objects"),
    make_option(c("-i","--it_no"),type="numeric",default = NULL, help = "iteration number for perms")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$file)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n",call.=FALSE)
}



b<-get(load(opt$file))

p<-lapply(types,function(t){
    message(t)
    mat<-matrix(NA,nrow=nrow(b$sigma),ncol=n.perms)
    cv.vec<-numeric(length = n.perms)
    for(i in 1:n.perms){
        tmp<-sim_alt(b$support,b$gt,b$sigma,t,ncp=ncp)
        mat[,i]<-tmp$rd
        cv.vec[i]<-tmp$cv.pos
    }
    list(perms=cbind(chr='1',b$support[b$support$pir,c('start','end')],mat),caus.var.pos<-cv.vec)
    
})
names(p)<-types
save(p,file=file.path(opt$out_dir,paste(opt$it_no,basename(opt$file),sep='_')))
message(paste("Written",file.path(opt$out_dir,paste(opt$it_no,basename(opt$file),sep='_'))))


## TEST CODE FOR DEBUGGING sim_alt function not run except for debugging purposes
if(1==0){
pl<-do.call('rbind',lapply(1:20,function(i){
    message(i)
    sz.obj<-sim_alt(b$support,b$gt,b$sigma,t,ncp=ncp)
    lp<--log10(pnorm(abs(sz.obj$rd), lower.tail = FALSE))
    cvp<-b$support[b$support$pir,]$start %in% sz.obj$cv.pos
    cv.type<-paste('cv',t,sep='.')
    sol<-factor(character(length=sum(b$support$pir)),levels=c('test','control',cv.type))
    sol[b$support[b$support$pir,]$test]<-'test'
    sol[b$support[b$support$pir,]$control]<-'control'
    sol[cvp]<-cv.type
    out<-data.frame(pos=b$support[b$support$pir,]$start,sim.z=sz.obj$rd,lp=lp,it=i,sol=sol)
}))
ggplot(pl,aes(x=pos,y=lp,color=sol)) + geom_point() + facet_wrap(~it) + theme_bw() + geom_hline(yintercept=8)
}