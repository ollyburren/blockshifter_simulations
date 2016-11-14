## code to simulate  a GWAS under the NULL using code from VSEAMS

library(snpStats)
library(data.table)
library(GenomicRanges)
## use chr1 EUR data 
DATA.DIR<-'/Users/oliver/DATA/1KG/'
## loads into snps object
(load(file.path(DATA.DIR,'chr1.RData')))
## clean to remove snps below 1% poor call rate and violating HWE
snp.sum<-col.summary(snps$genotypes)
maf<-split(snp.sum$MAF,row.names(snp.sum))
keep<-with(snp.sum,which(MAF>0.01 & z.HWE^2<25 & Call.rate>0.95))
snp.sum<-snp.sum[keep,]
gt<-snps$genotypes[,keep]
sup<-snps$support[keep,]
sup$MAF<-unlist(maf[sup$snp.names])
##next load in 0.1cM LD blocks.
ld<-fread('/Users/oliver/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37.bed')
setnames(ld,c('chr','start','end','region'))
ld.1<-subset(ld,chr=='1')
ld.gr<-with(ld.1,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start+1,end=end),region=region))
snp.gr<-with(sup,GRanges(seqnames=Rle('1'),ranges=IRanges(start=start,end=end),snp.name=snp.names))
m<-as.matrix(findOverlaps(ld.gr,snp.gr))
m<-cbind(m,ld.gr[m[,1],]$region)
by.block<-lapply(split(m[,2],m[,3]),function(x){
    x<-as.numeric(x)
    list(support=sup[x,],gt=gt[,x])
})

## make blockname a variable

out.dir<-file.path(DATA.DIR,'chr1_EUR_LD')
by.block<-lapply(names(by.block),function(n){
    b<-by.block[[n]]
    b$block<-n
    save(b,file=file.path(out.dir,paste0(n,'.RData')))
    b
})



