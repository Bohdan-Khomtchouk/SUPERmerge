#!/usr/bin/env/Rscript


suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
option_list = list(
  make_option(c("-r", "--rmdup"), action="store", default="no", type='character',
              help="Remove Duplicates"),
  make_option(c("-d", "--depth"), action="store", default="20", type='integer',
              help="Depth filter"),
  make_option(c("-i", "--interval"), action="store", default="500", type='integer',
              help="Interval specification"),
  make_option(c("-f", "--gtf"), action="store", default="NA", type='character',
              help="GTF file location"),
  make_option(c("-c", "--control"), action="store", default="NA", type='character',
              help="Control bam file location"),
  make_option(c("-p", "--ip_bam"), action="store", default="NA", type='character',
              help="IP bam file location"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Should the program print extra stuff out? [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))
if (opt$v) {
    cat("Please cite SUPERmerge script at : \n")
    cat("Input Parameters\n================================================================")
    cat("\nr:\n")
    cat(opt$r)
    cat("\n\nd:\n")
    cat(opt$d)
    cat("\n\ni:\n")
    cat(opt$interval)
    cat("\n\nf:\n")
    cat(opt$gtf)
    cat("\n\nc:\n")
    cat(opt$control)
    cat("\n\np:\n")
    cat(opt$ip_bam)
    cat("\n\n")
}

if (is.na(opt$interval)) {
	cat("The interval number you have specified is not a number\n")
	quit()
}
if (is.na(opt$depth)) {
	cat("The depth number you have specified is not a number\n")
	quit()
}
# Check if the ip and gtf files exist
if (!file.exists(opt$ip_bam)) {
	cat("The IP bam file does not exist\n")
	quit()
}
if (!file.exists(opt$gtf)) {
	cat("The gtf file you have specified does not exist\n")
	quit()
}
# check if control file specified and if so does it exist
if (opt$control != "NA") {
	if (!file.exists(opt$control)) {
		cat("The control bam file you have specified does not exist\n")
		quit()
	}
}
cat("Parameters pass initial sanity check\n")
cat("Checking if remove duplicates is necessary\n")

if (opt$rmdup == "yes") {
	require(dupRadar)
	bamDuprm <- markDuplicates(dupremover="bamutil", bam=opt$ip_bam, path="/hihg/smoke/applications/bamUtil/bamUtil/bin", rminput=FALSE, out=gsub("\\.bam$","_duprm.bam --rmDups", opt$ip_bam))
	if (opt$control != "NA") {
		bamDuprmC <- markDuplicates(dupremover="bamutil", bam=opt$control, path="/hihg/smoke/applications/bamUtil/bamUtil/bin", rminput=FALSE, out=gsub("\\.bam$","_duprm.bam --rmDups", opt$ip_bam))
	}
}

cat("Past removal duplicates now to pileup function\n")

require(Rsamtools)

# index bam file and use the pileup command to calculate depth
indexBam(opt$ip_bam)
sparam <- PileupParam(max_depth=250, min_nucleotide_depth=opt$depth)
res <- pileup(opt$ip_bam)
# read in the gtf file for annotation
require(rtracklayer)
grd<-import.gff(opt$gtf)
# start counting the ranges and reduce it to unique pseudo-peaks
require(GenomicRanges)
ir <- GRanges(seqnames=Rle(res$seqnames), ranges=IRanges(start=res$pos, width=1))
gr<-GRanges(seqnames=Rle(res$seqnames), ranges=IRanges(start=res$pos-opt$depth, width=(opt$depth * 2))
gr_reduce<-reduce(gr)
# find the nearest peak for each of the positions that meet criteria (this will eventually give us frequency)
n<-nearest(ir, gr_reduce)
dfn<-as.data.frame(table(n))
# subset the gtf file by gene_biotype of gene.  Make this a requirement.
txs<-grd[grd$type=="gene"]
# annotate the reduced pseudo peaks to the newly formed gtf GRanges object
n2<-nearest(gr_reduce, txs)
# create final
final<-txs[n2]

output_df<-data.frame(df$seqnames, df$start, df$end, dfn$Freq, (dfn$Freq / df$width), (dfn$Freq / (df$width - opt$interval)), final$seqnames, final$start, final$end, final$gene_id, final$gene_type)
colnames(output_df)<-c("Chr", "Start", "Stop", "Bases", "% bases meeting criteria", "% interval - extensions", "GTF chr", "GTF start", "GTF stop", "GTF gene_id", "GTF gene_type")
write.csv(file="supermerge_results.csv", output_df)
cat ("END OF RUN\n")


