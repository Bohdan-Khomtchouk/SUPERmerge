#!/usr/bin/env/Rscript


suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
option_list = list(
  make_option(c("-r", "--rmdup"), action="store", default="no", type='character',
              help="Remove Duplicates"),
  make_option(c("-d", "--depth"), action="store", default="20", type='integer',
              help="Depth filter"),
  make_option(c("-i", "--interval"), action="store", default="500", type='integer',
              help="Interval specification"),
  make_option(c("-p", "--threads"), action="store", default="1", type='integer',
              help="# of Threads to Run"),
  make_option(c("-f", "--gtf"), action="store", default="NA", type='character',
              help="GTF file location"),
  make_option(c("-c", "--ip_bam"), action="store", default="NA", type='character',
              help="IP bam file location"),
  make_option(c("-s", "--stranded"), action="store", default="NA", type='character',
              help="Enforce strandedness"),
  make_option(c("-l", "--filelist"), action="store", default="NA", type='character',
              help="List of bam files"),
  make_option(c("-o", "--output"), action="store", default="NA", type='character',
              help="Output prefix"),
  make_option(c("-t", "--read_type"), action="store", default="single", type='character',
              help="Read Type"),
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
    cat("\n\ns:\n")
    cat(opt$stranded)
    cat("\n\np:\n")
    cat(opt$ip_bam)
    cat("\n\nl:\n")
    cat(opt$filelist)
    cat("\n\ns:\n")
    cat(opt$read_type)
    cat("\n\no:\n")
    cat(opt$output)
    cat("\n\n")
}

if (is.na(opt$threads)) {
	opt$threads = 1
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
	if (!file.exists(opt$filelist)) {
		cat("Invalid Input\n")
		quit()
	}
	else {
		files2=scan(opt$filelist, character(), quote="")
		cat("Input files : \n")
		cat(files2)
		cat("\n\n")
	}
}
if (!file.exists(opt$gtf)) {
	cat("The gtf file you have specified does not exist\n")
	quit()
}
cat("Parameters pass initial sanity check\n")

require(Rsamtools)
require(rtracklayer)
require(GenomicRanges)
require(encodeChIPqc)
require(dplyr)
require(snow)

# first thing to do is read in the annotation
grd<-import.gff(opt$gtf)

# Main function ...
supermerge<-function(bamname, interval, depth, dup, out, strand, cflag) {
	
        cat("Checking if remove duplicates is necessary\n")

        if (dup == "yes") {
                bamDuprm <- markDuplicates(dupremover="bamutil", bam=bamloc, path="/hihg/smoke/applications/bamUtil/bamUtil/bin", rminput=FALSE, out=gsub("\\.bam$","_duprm.bam --rmDups", bamname))
                bamname <- gsub("\\.bam$","_duprm.bam", bamname)
        }

        # Check if index exists, and if it does not then create one with Rsamtools
        if (!file.exists(paste0(bamname, ".bai"))) {
                indexBam(bamname)
        }

        cat("Running pileup function\n")

        # create pileup for interval building
        sparam <- PileupParam(max_depth=250, min_nucleotide_depth=depth, distinguish_strands=FALSE)
        res <- pileup(bamname, pileupParam=sparam)
        # start counting the ranges and reduce it to unique pseudo-peaks
        cat("Past pileup now annotating\n")
        # possibility exists where chromosomes do not exist in gtf file so remove these from the object
        chrs<-as.data.frame(runValue(grd@seqnames))
        names(chrs)[1] <- "chr"
        names(res)[1]<-"chr"
        res_final<-semi_join(res, chrs, by="chr")
        res<-res_final
        names(res)[1]<-"seqnames"

        ir <- GRanges(seqnames=Rle(res$seqnames), ranges=IRanges(start=res$pos, width=1))
        gr<-GRanges(seqnames=Rle(res$seqnames), ranges=IRanges(start=res$pos-interval, width=(interval * 2)))
        gr_reduce<-reduce(gr)

        grdf<-as.data.frame(gr_reduce)
        gr_reduce2<-GRanges(seqnames=Rle(grdf$seqnames), ranges=IRanges(start=grdf$start+interval, width=(grdf$width - (2 * interval))))

        # find the nearest peak for each of the positions that meet criteria (this will eventually give us frequency)
        if (strand=="yes") {
                n<-nearest(ir, gr_reduce, ignore.strand=FALSE)
        } else {
                n<-nearest(ir, gr_reduce)
        }
        df<-as.data.frame(gr_reduce)
        dfn<-as.data.frame(table(n))
        # subset the gtf file by gene_biotype of gene.  Make this a requirement.
        txs<-grd[grd$type=="gene"]
        # annotate the reduced pseudo peaks to the newly formed gtf GRanges object
        if (strand=="yes") {
                n2<-nearest(gr_reduce, txs, ignore.strand=FALSE)
        } else {
                n2<-nearest(gr_reduce, txs)
        }
        # create final
        n_d2<-as.data.frame(n2)
        n_d2[is.na(n_d2)] <- 0
        n2<-as.integer(n_d2[,1])
        final2<-txs[n2]
        final<-as.data.frame(final2)

	bamView<-BamViews(bamname)
        bamRanges(bamView) <- gr_reduce
        aln<-countBam(bamView)
        a<-aln@listData
        adf<-as.data.frame(a)
        names(adf)[6]<-"records"

        output_df<-data.frame(df$seqnames, df$start+interval, df$end-interval, dfn$Freq, (dfn$Freq / (df$width - (2 * interval) + 2)), final$seqnames, final$start, final$end, final$gene_id, final$gene_type, adf$records)
        colnames(output_df)<-c("Chr", "Start", "Stop", "Bases", "% bases meeting criteria", "GTF chr", "GTF start", "GTF stop", "GTF gene_id", "GTF gene_type", "Reads in Expanded Interval")
	return(output_df)
}

counted<-function(bamname, interval, depth, dup, out, strand, cflag) {

	aln<-countBam(bamname)
	output_info<-aln
	return(output_info)
}

pfrip<-0
cfrip<-0
cat("Run the main bam file through the function\n")
cat(paste("Creating a cluster with ",opt$threads," threads\n"))
cluster<-makeCluster(opt$threads)

inte<-opt$interval
dep<-opt$depth
rmd<-opt$rmdup
outp<-opt$output
stra<-opt$stranded


library(parallel)
ip_info<-mclapply(files2, function(x) counted(x, inte, dep, rmd, outp, stra), mc.cores=opt$threads)
ip_df<-mclapply(files2, function(x) supermerge(x, inte, dep, rmd, outp, stra), mc.cores=opt$threads)
names(ip_df)<-files2

# Now we have a list of the returned results from the main part of the program.  Write output here
# First calculate complete number of peaks
# Create the proper statistics in data frames
usq<-0
usq2<-0
usq3<-0
for(i in 1:length(ip_df)) {
	# Get the total number of reads in the file
	usq[i]<-as.integer(ip_info[[i]][,6])
	# Get the sum of the  number of reads
	b<-as.character(ip_info[[i]][,5])
	gr<-GRanges(seqnames=ip_df[[b]]$Chr, IRanges(start=ip_df[[b]]$Start, end=ip_df[[b]]$Stop))
	usq2[i]<-sum(width(gr))
	usq3[i]<-length(gr)
}
st<-data.frame(TotalReads=usq, PeakReads=usq2, FileName=files2, PeakNumber=usq3)

# Draw the graphical output here
library(ggplot2)
bar1<-data.frame(st$PeakReads / st$TotalReads, st$FileName)
names(bar1)<-c("Perc", "FileName")
p<-ggplot(data=bar1, aes(x=FileName, y=Perc))+
geom_bar(stat="identity", fill="steelblue")+
theme_minimal()
p
bar2<-data.frame(st$PeakNumber, st$FileName)
names(bar2)<-c("PeakNumber", "FileName")
p2<-ggplot(data=bar2, aes(x=FileName, y=PeakNumber))+
geom_bar(stat="identity", fill="steelblue")+
theme_minimal()
p2
sbar1<-data.frame(FileName=st$FileName, ReadType=c("peak"), Reads=st$PeakReads)
sbar2<-data.frame(FileName=st$FileName, ReadType=c("other"), Reads=st$TotalReads - st$PeakReads)
sbar<-rbind(sbar1, sbar2)
p3<-ggplot(sbar, aes(fill=ReadType, y=Reads, x=FileName))+
geom_bar(position="stack", stat="identity")
p3
p4<-ggplot(sbar, aes(fill=ReadType, y=Reads, x=FileName))+
geom_bar(position="fill", stat="identity")
p4
dev.off()

cat ("END OF RUN\n")
