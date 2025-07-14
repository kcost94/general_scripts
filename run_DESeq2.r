#!/usr/bin/env R

library('DESeq2')

## Directory with all counts files

directory<-"./"

## link to all files in directory -- call them sample files

sampleFiles<- grep(".counts",list.files(directory),value=TRUE)
sampleFiles

sampleNames <- sampleFiles
sampleNames

## setup groups in order of file names

sampleCondition<-c('control','control','control','fat','fat','fat')

sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, condition=sampleCondition)

treatments = c("control","fat")

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels = treatments)

dds = DESeq(ddsHTSeq, test = "LRT", reduced = ~1)
res <- results(dds)
f <-res[which(res$padj<0.05),]
write.table(file='diff-exp-fdr-0.05.txt', f, quote=F, col.names=T, sep='\t')
f <-res[which(res$padj<0.1),]
write.table(file='diff-exp-fdr-0.1.txt', f, quote=F, col.names=T, sep='\t')
f <-res[which(res$pvalue<0.05),]
write.table(file='diff-exp-pval.05.txt', f, quote=F, col.names=T, sep='\t')
f <-res[which(res$padj<=1),]
write.table(file='diff-exp-fdr.txt', f, quote=F, col.names=T, sep='\t')

head(res)
plotMA(dds,ylim=c(-2,2),main="DESeq2")
dev.copy(png,"deseq2_MAplot.png")
dev.off()


rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

pdf('deseq2_pca.pdf')
plotPCA(rld, intgroup=c("condition"))
dev.off()


norm.counts <-as.data.frame(counts(dds, normalized=T))
write.table(norm.counts, file='norm-all-data.txt', quote=F, col.names=T, sep='\t')

