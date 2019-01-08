##############################################################################
## PROJECT: DATA ANALYSIS: NGS (TASK SOPHIA GENETICS)
## AUTHOR: PAURUSH PRAVEEN
## DATE: 10 JANUARY 2015
## CONTACT: praveen@bit.uni-bonn.de
## URL: www.paurushpraveen.co.in
##############################################################################
## Libraries

library(Rsamtools)
library(ShortRead)
library(GenomicAlignments)

## Input files
DIR = "/home/praveen/Paurush/SG_tasks"
fastqFiles=list.files(DIR, pattern="*.gz")
bamFile = list.files(DIR, pattern="*.bam")
OutDir = file.path(DIR, "Results")

## Read fastq files 


seq=list()
for(i in 1:length(fastqFiles)){
	seq[[i]] = readFastq(file.path(DIR, fastqFiles[i]))
}

## Read bam files 
bam = scanBam(file.path(DIR, bamFile[1]))
countBam(file.path(DIR, bamFile[1]))
aln = readGAlignments(file.path(DIR, bamFile[1]), use.names=TRUE)



## Phred Quality Assesment

assesPhredQuality=function(seqs, outdir){
	for(i in 1:length(seqs)){
		readM=as(FastqQuality(quality(quality(seqs[[i]]))), "matrix")
		pdf(paste("outdir", "phreqQualityBoxPlot", i, ".pdf", sep=""))
		boxplot(as.data.frame((readM)), outline = FALSE, main=paste("Per Cycle Read Quality fastq:",i, sep=""), xlab="Cycle", ylab="Phred Quality", las=2)
		lines(x=c(1:ncol(readM)), y=rep(20,ncol(readM)), col="red")
		dev.off()
	}
}

assesPhredQuality(seqs, outdir=OutDir)



##

qas = lapply(seq_along(file.path(DIR, fastqFiles)), function(i, file.path(DIR, fastqFiles)) qa(readFastq(file.path(DIR, fastqFiles)), namesfile.path(DIR, fastqFiles)), file.path(DIR, fastqFiles))

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/fastqQuality.R")



which <- RangesList(seq1=IRanges(1000, 2000), seq2=IRanges(c(100, 1000), c(1000, 2000)))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)
bamFile <-system.file("extdata", "ex1.bam", package="Rsamtools")
bam2 <- scanBam(bamFile, param=param)
bam <- scanBam("bamExample.bam")
countBam("bamExample.bam")
what <- c("rname", "strand", "pos", "seq")
param<-ScanBamParam(what=what)
bam2<-scanBam("bamExample.bam", param=param)
samp <- BamSampler("bamExample.bam", yieldSize=1000)
library(GenomicAlignments)
gal2 <- readGAlignments(bamFile, use.names=TRUE)
gal2
head(names(gal2))





fastq <- list.files("data", "*.fastq$")
fastq <- paste("data/", fastq, sep="")
names(fastq) <- paste("flowcell6_lane", 1:length(fastq), sep="_") # Values in name slot will be used for sample labeling in plot!
## Compute quality stats
fqlist <- seeFastq(fastq=fastq, batchsize=100000, klength=8) # Returns summary stats in a relatively small list object that can be saved to disk with 'save()' and reloaded with 'load()' for later plotting. The argument 'klength' specifies the k-mer length and 'batchsize' the number of reads to random sample from each fastq file.  
## Plot quality stats 
seeFastqPlot(fqlist)
seeFastqPlot(fqlist[4:1], arrange=c(1,2,3,4,6,7)) # Example for plotting specific rows and columns.
## Output plot to PDF
pdf("fastqReport.pdf", height=18, width=4*length(fastq))
seeFastqPlot(fqlist)
dev.off()







library(ShortRead)
fls <- list.files(fastqDir, "fastq$", full=TRUE)
names(fls) <- sub(".fastq", "", basename(fls))
## use FastqSampler if fastq files are large
qas <- lapply(seq_along(fls), function(i, fls) qa(readFastq(fls[i]), names(fls)[i]), fls)
qa <- do.call(rbind, qas)
save(qa, file=file.path(outputDir, "qa.rda")
browseURL(report(qa))

library(TxDb.Scerevisiae.UCSC.sacCer2.sgdGene)
txdb <- Scerevisiae_UCSC_sacCer2_ensGene_TxDb
gnModel <- exonsBy(txdb, "gene")

bamFls <- list.files(bamDir, "bam$", full=TRUE)
names(bamFls) <- sub("\\..*", "", basename(bamFls))
counter <- function(fl, gnModel)
{
aln <- readGappedAlignments(fl)
strand(aln) <- "*" # for strand-blind sample prep protocol
hits <- countOverlaps(aln, gnModel)
counts <- countOverlaps(gnModel, aln[hits==1])
names(counts) <- names(gnModel)
counts
}
counts <- sapply(bamFls, counter, gnModel)
save(counts, file=file.path(outputDir, "counts.rda"))

library(edgeR)

## identify treatment groups
grp <- factor(<...>)

## create data structures
dge <- DGEList(counts, group=grp)
dge <- calcNormFactors(dge)

## filter uniformative genes
m <- 1e6 * t(t(dge$counts) / dge$samples$lib.size)
ridx <- rowSums(m > 1) >= 2
dge <- dge[ridx,]

## comparison between groups
design <- model.matrix( ~ grp )
dge <- estimateCommonDisp(dge, design)
fit <- glmFit(dge, design, dispersion=dge$common.dispersion)
lrTest <- glmLRT(dge, fit, coef=2)
tt <- topTags(lrTest, Inf)
save(tt, file=file.path(dataDir, "tt.rda"))


myvcf=read.csv('/home/praveen/Downloads/GalaxyFreeBayes.tab', head=TRUE, sep="\t")
library(ggplot2)
p=ggplot(myvcf, aes(x=factor(POS), y=QUAL, fill=CHROM))+geom_bar(stat="identity", position="dodge")+ylab("Quality score")+xlab("Position in respective chromosome")+theme(panel.background=element_blank())

