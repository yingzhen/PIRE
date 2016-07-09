#download and install SNPRelate
#source("http://bioconductor.org/biocLite.R")
#biocLite("SNPRelate")

library(gdsfmt)
library(SNPRelate)
setwd("/Users/Ying/Desktop/analysis/")
setwd("/Users/Ying/Desktop/analysis/greenbul/RAD/all-rerun/population-3/")
#format conversion from VCF
vcf.fn = "batch_4.vcf"
snpgdsVCF2GDS(vcf.fn, "greenbul.gds", method="biallelic.only")
snpgdsSummary("greenbul.gds")

#data analysis
genofile <- snpgdsOpen("greenbul.gds")
pca <- snpgdsPCA(genofile, autosome.only = F, maf=0.02)
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))

popmap = read.table("../greenbul-good-popmap-habitat-rm2TibatiAllLope")
sample.id = as.character(popmap$V1)
pop_code = as.character(popmap$V2)
habitat = as.character(popmap$V3)

#color by habitat
palette("default")
tab1 <- data.frame(sample.id = pca$sample.id, hab = factor(habitat)[match(pca$sample.id, sample.id)], EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
#tab1 = subset(tab1, )
pdf("PCA-byhabitat-maf02.pdf", width=6, height=5)
#pdf("PCA-alltags-byhabitat.pdf", width=6, height=5)
plot(tab1$EV2, tab1$EV1, col=as.integer(tab1$hab), xlab=paste("PC2 (", pc[2],"%)"), ylab=paste("PC1 (", pc[1],"%)"))
legend("bottomleft", legend=levels(tab1$hab), pch="o", col=1:nlevels(tab1$hab))
dev.off()
pdf("PCA-byhabitat-maf02-PC3.pdf", width=6, height=5)
lbls = paste("PC", 3:6, "\n", format(pc.percent[3:6], digits=2),"%", sep="")
pairs(pca$eigenvect[,3:6], col=tab1$hab, labels=lbls)
dev.off()

#color by population
pdf("PCA-bypop-maf02.pdf", width=6, height=5)
palette(c("royalblue1", "red1", "gold","darkblue", "darkorange", "gray58", "purple", "greenyellow", "violetred1", "black", "turquoise1", "forestgreen", "salmon4", "darkorchid4", "seagreen1"))
#palette(rainbow(18))
tab <- data.frame(sample.id = pca$sample.id, pop = factor(pop_code)[match(pca$sample.id, sample.id)], EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab=paste("PC2 (", pc[2],"%)"), ylab=paste("PC1 (", pc[1],"%)"))
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop), cex=0.6)
dev.off()

pdf("PCA-bypop-maf02-PC3.pdf", width=6, height=5)
lbls = paste("PC", 3:6, "\n", format(pc.percent[3:6], digits=2),"%", sep="")
pairs(pca$eigenvect[,3:6], col=tab$pop, labels=lbls)
dev.off()

#snploading of each pc
PCARV <- snpgdsPCA(genofile, eigen.cnt=8, autosome.only = F)
SnpLoad <- snpgdsPCASNPLoading(PCARV, genofile)
names(SnpLoad)
#dim(SnpLoad$snploading)
pdf("PCA-alltags-rmOutlierLope-loading.pdf", width=6, height=5)
plot(SnpLoad$snploading[1,], type="h", ylab="PC 1")
plot(SnpLoad$snploading[2,], type="h", ylab="PC 2")
dev.off()

which(SnpLoad$snploading[2,] < -0.07)

diss <- snpgdsDiss(genofile, autosome.only=F, maf=0.02)
hc <- snpgdsHCluster(diss)
rv <- snpgdsCutTree(hc, label.H=TRUE, label.Z=TRUE)
pdf("PCA-alltags-rmOutlierLope-tree.pdf", width=6, height=5)
snpgdsDrawTree(rv, edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"), leaflab="perpendicular")
dev.off()

#To calculate the SNP correlations between eigenvactors and SNP genotypes:
CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)
savepar <- par(mfrow=c(3,1), mai=c(0.3, 0.55, 0.1, 0.25))
for (i in 1:3)
{
  plot(abs(CORR$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste("PC", i), pch="+")
}
