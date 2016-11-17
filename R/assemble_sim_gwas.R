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


## MAIN
ps.file<-'/home/ob219/scratch/bs_sim/support/perm_support.RData'
if(!file.exists('/home/ob219/scratch/bs_sim/support')){
    chk.res<-create_perm_support_df(DATA.DIR)
    save(chk.res,ps.file)
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

## 8 regions don't contain any test but that is not important for
## initial analysis

## FIRST ANALYSIS - DO WE GET AN ACCURATE NULL UNDER THE ALT ?
## simulate different numbers of blocks having p = 0.5 for cv being in PIR or in NON-PIR
## Try 1K permutations

# do 200 at a time
# number of blocks
m<-opt$it_no
it_no<-opt$itno
no.perms.it<-200

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
    bp<-cbind(t$null$perms[,1:3],bp)
})
perms<-rbindlist(all.perms)
of<-paste0(paste(it_no,m,sep='-'),'.RData')
save(perms,file=file.path(opt$out_dir,of))
message(paste("Written",file.path(opt$out_dir,of)))






