library(data.table)
library(parallel)

### code to prepare recipes for generating simulations under the null

DATA.DIR<-'/home/ob219/scratch/bs_sim/sims/'

files<-list.files(path=DATA.DIR,pattern='*.RData')

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



## FIRST ANALYSIS - DO WE GET AN ACCURATE NULL UNDER THE ALT ?
## simulate different numbers of blocks having p = 0.5 for cv being in PIR or in NON-PIR
## Try 1K permutations

# do 200 at a time
# number of blocks
#m<-50
#it_no<-1
#pos.blocks<-



