library(data.table)
library(parallel)

### code to prepare recipes for generating simulations under the null

DATA.DIR<-'/home/ob219/scratch/bs_sim/sims/'

create_perm_support_df<-function(d){
    files<-list.files(path=d,pattern='*.RData')

    ## organise files

    all.sims<-mclapply(files,function(x){
        itno<-sub("([0-9]+)\\_.*","\\1",x)
        s<-sub("[0-9]+\\_1:([^\\-]+)\\-.*","\\1",x)
        e<-sub("[0-9]+\\_1:[^\\-]+\\-([^\\.]+).*","\\1",x)
        r<-sub("[0-9]+\\_1:([^\\-]+\\-[^\\.]+).*","\\1",x)
        data.table(filename=x,itno=itno,start=s,end=e,r=r)
    },mc.cores = 8)
    all.sims<-rbindlist(all.sims)
    ## first off check and annotate all.sims
    chk<-mclapply(all.sims$filename,function(f){
        message(f)
        t<-get(load(file.path(DATA.DIR,f)))
        ## check perm counts for each type
        data.frame(lapply(t,function(x){
            if(!is.data.frame(x$perms))
                return(0)
            return(ncol(x$perms)-3)
        }))
    },mc.cores=8)

    chk.res<-cbind(all.sims,do.call("rbind",chk))
    chk.res<-chk.res[order(chk.res$r,chk.res$itno),]
    return(chk.res)
}

## robustly sum logs
logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}

## compute variance shrinkage for case control study
Var.data.cc <- function(f, N, s) {
    1 / (2 * N * f * (1 - f) * s * (1 - s))
}

## compute approx bayes factors and resultant posterior probabilities
## based on the assumption of one causal variant in a region
approx.bf.p <- function(z,f,type, N, s,pi_i,suffix=NULL) {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    sBF <- logsum(lABF + log(pi_i))
    exp(lABF + log(pi_i))/(exp(sBF) + 1)
}


## MAIN
ps.file<-'/home/ob219/scratch/bs_sim/support/perm_support.RData'
if(!file.exists(ps.file)){
    chk.res<-create_perm_support_df(DATA.DIR)
    save(chk.res,file=ps.file)
}else{
    chk.res<-get(load(ps.file))
}


library(optparse)

##MAIN

## OPTION PROCESSING

option_list = list(
    make_option(c("-b","--blockno"),type="character",default = NULL, help = "Number of blocks to contain a causal variant"),
    make_option(c("-o","--out_dir"),type="character",default = NULL, help = "Directory to store result data objects"),
    make_option(c("-i","--it_no"),type="numeric",default = NULL, help = "iteration number for perms")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$blockno)){
    print_help(opt_parser)
    stop("blockn is a required parameter",call.=FALSE)
}

if(is.null(opt$out_dir)){
    print_help(opt_parser)
    stop("please supply an output dir",call.=FALSE)
}

if(is.null(opt$it_no)){
    print_help(opt_parser)
    stop("please supply an iteration number (or block of perms)",call.=FALSE)
}

opt<-list(it_no=1,out_dir='/home/ob219/scratch/bs_sim/tmp',blockno=50)

## 8 regions don't contain any test but that is not important for
## initial analysis

## FIRST ANALYSIS - DO WE GET AN ACCURATE NULL UNDER THE ALT ?
## simulate different numbers of blocks having p = 0.5 for cv being in PIR or in NON-PIR
## Try 1K permutations

# do 200 at a time
# number of blocks
m<-opt$blockno
it_no<-opt$it_no
no.perms.it<-200
N.ok<-14361+43923
case.ratio<-14361/N.ok
support.dir<-'/home/ob219/scratch/bs_sim/sigma_pir/'

bf<-subset(chk.res,itno==it_no)
bf<-bf[order(as.numeric(bf$start)),]
## ok we need to select m blocks 200 times
bar<-do.call('cbind',lapply(1:no.perms.it,function(i){1:nrow(bf) %in% sample(1:nrow(bf),m)}))
recipe<-cbind(bf,bar)
pnames<-names(recipe)[grep("^V[0-9]+",names(recipe))]
all.perms<-lapply(recipe$filename,function(f){
    message(f)
    t<-get(load(file.path(DATA.DIR,f)))
    perm.vec<-as.logical(recipe[1,pnames,with=FALSE])
    sel<-character(length=no.perms.it)
    sel[!perm.vec]<-'null'
    sel[perm.vec]<-sample(c('pir','no.pir'),sum(perm.vec),replace=TRUE,prob=rep(0.5,2))
    ##create matrix of Z scores based on this
    bp<-do.call('cbind',lapply(1:no.perms.it,function(i){
        sel[i]
        t[[sel[i]]]$perms[i+3]
    }))
    ## conver to Posterior Probabilities for this we need to load support file which has the MAF
    ## first load in the support object
    support.file<-file.path(support.dir,sub("[^\\_]+\\_(.*)","\\1",f))
    supp<-get(load(support.file))     
    bp<-apply(bp,2,approx.bf.p,f=supp$support[supp$support$pir,]$MAF,type='CC',N=N.ok,s=case.ratio,pi_i=1e-4)
    if(sum(supp$support$pir)==1)
	bp<-matrix(bp,ncol=no.perms.it)
    bp<-data.table(bp)
    setnames(bp,paste0('V',1:no.perms.it))
    bp<-cbind(t$null$perms[,1:3],bp)
})
perms<-rbindlist(all.perms)
of<-paste0(paste(it_no,m,sep='-'),'.RData')
save(perms,file=file.path(opt$out_dir,of))
message(paste("Written",file.path(opt$out_dir,of)))
