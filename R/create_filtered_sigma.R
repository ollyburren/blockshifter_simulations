library(GenomicRanges)
library(snpStats)
library(corpcor)
library(mvtnorm)
library(optparse)
## compute sigmas and annotate support files for simulations

DATA.DIR<-'/home/ob219/scratch/bs_sim/support'


## creates a cutdown contact file that we can use to filter and annotate snpStats
## objects when we are computing and storing sigma (covariance matrices)
make_simple_contact_file<-function(){
    library(data.table)
    test.set<-c('Total_CD4_Activated','Total_CD4_NonActivated')
    control.set<-c('Megakaryocytes','Erythroblasts')
    contacts<-fread(file.path(DATA.DIR,'merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab'))
    chic.thresh<-5
    test<-contacts$Total_CD4_Activated>chic.thresh | contacts$Total_CD4_NonActivated>chic.thresh
    control<-contacts$Megakaryocytes>chic.thresh | contacts$Erythroblasts>chic.thresh
    contacts<-cbind(contacts[,1:13,with=FALSE],test,control)
    ## remove contacts which are not found in these tissues
    contacts.f<-contacts[rowSums(cbind(test,control))!=0,]
    ## only care about chr1
    contacts.f<-subset(contacts.f,oeChr=='1' & baitChr=='1')
    contacts.gr<-with(contacts.f,GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,end=oeEnd),test=test,control=control))
    save(contacts.gr,file=file.path(DATA.DIR,'pirs.RData'))
}

attempt.pos.def<-function(mat,diag.val=1.0001){
    print(paste("diag.val",diag.val))
    if(!is(mat,"Matrix"))
        stop("mat is not a Matrix!")
    if(diag.val >= 1.1){
        print("Matrix is not positive definite. Finding closest approximation..")
        diag(mat)<-1
        return(as(make.positive.definite(mat),"Matrix"))
    }
    diag(mat)<-diag.val
    if(is.positive.definite(mat,,method="chol")==FALSE){
        new.diag<-signif(1+((diag.val-trunc(diag.val))*10))
        mat<-attempt.pos.def(mat,new.diag)
    }else{
        return(mat)
    }
}

mvs.sigma.r2<-function(r2){
    diag(r2)<-1
    if(!is.positive.definite(r2,,method="chol")){
        #this recurses through various values of diag if we exceed 1 then
        #we compute the closest matrix that is positive definite.
        r2<-attempt.pos.def(r2)
    }
    r2
}

mvs.perm<-function(sigma,n=1000){
    if(!is.matrix(sigma))
        stop("sigma parameter is not a matrix")		
    if(!is.positive.definite(sigma,,method="chol"))
        stop("sigma is not positive definite")
    ## in original paper method="chol" was not defined so I assume used eigen default
    ## this is slower than the choleski decomp ! Perhaps we should contact the author ?
    rd<-rmvnorm(n,mean=rep(0,ncol(sigma)),sigma=sigma,method="chol")
    t(rd)
}

option_list = list(
    make_option(c("-f","--file"),type="character",default = NULL, help = "Input RData object see - prepare_eur_1kg_files.R"),
    make_option(c("-o","--out_dir"),type="character",default = NULL, help = "Directory to store result data objects")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$file)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n",call.=FALSE)
}

## this is generated from using make_simple_contact_file above
if(!file.exists(file.path(DATA.DIR,'pirs.RData')))
    make_simple_contact_file()

contact.file<-file.path(DATA.DIR,'pirs.RData')
c.gr<-get(load(contact.file))
## next load in LD block file

b<-get(load(opt$file))
## first generate sigma for just SNPs that overlap either test or control
support<-b$support
support.gr<-with(support,GRanges(seqnames=Rle('1'),ranges=IRanges(start=start,end=end)))
for(cn in c('test','control')){
    tmp<-logical(length=nrow(support))
    cgr<-contacts.gr[mcols(contacts.gr)[[cn]],]
    ol<-as.matrix(findOverlaps(support.gr,cgr))
    tmp[ol[,1]]<-TRUE
    b$support[[cn]]<-tmp
}
b$support$pir<-rowSums(b$support[,c('test','control')])!=0
## next compute and store sigma cov matrix for just SNPs that overlap pirs
gt<-b$gt[,which(b$support$pir==TRUE)]
r2<-forceSymmetric(ld(x=gt,depth=ncol(gt)-1,stats="R.squared",symmetric=TRUE))
r2[which(is.na(r2))]<-0
b$sigma<-as.matrix(mvs.sigma.r2(r2))
save(b,file=file.path(opt$out_dir,basename(opt$file)))
message(paste("Written",file.path(opt$out_dir,basename(opt$file))))
    


