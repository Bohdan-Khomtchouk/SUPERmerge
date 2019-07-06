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
  make_option(c("-s", "--stranded"), action="store", default="NA", type='character',
              help="Enforce strandedness"),
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
    cat("\n\nc:\n")
    cat(opt$control)
    cat("\n\np:\n")
    cat(opt$ip_bam)
    cat("\n\nt:\n")
    cat(opt$read_type)
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

require(Rsamtools)
require(rtracklayer)
require(GenomicRanges)
require(encodeChIPqc)
require(dplyr)

# first thing to do is read in the annotation
grd<-import.gff(opt$gtf)

# Main function ...
supermerge<-function(bamname, cflag) {
	
	cat("Checking if remove duplicates is necessary\n")

	if (opt$rmdup == "yes") {
		bamDuprm <- markDuplicates(dupremover="bamutil", bam=bamloc, path="/hihg/smoke/applications/bamUtil/bamUtil/bin", rminput=FALSE, out=gsub("\\.bam$","_duprm.bam --rmDups", bamname))
		bamname <- gsub("\\.bam$","_duprm.bam", bamname)
	}

	# Check if index exists, and if it does not then create one with Rsamtools	
	if (!file.exists(paste0(bamname, ".bai"))) {
		indexBam(bamname)
	}

	cat("Running pileup function\n")
	# create pileup for interval building
	sparam <- PileupParam(max_depth=250, min_nucleotide_depth=opt$depth, distinguish_strands=FALSE)
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
	gr<-GRanges(seqnames=Rle(res$seqnames), ranges=IRanges(start=res$pos-opt$interval, width=(opt$interval * 2)))
	gr_reduce<-reduce(gr)

	grdf<-as.data.frame(gr_reduce)
	gr_reduce2<-GRanges(seqnames=Rle(grdf$seqnames), ranges=IRanges(start=grdf$start+opt$interval, width=(grdf$width - (2 * opt$interval))))

	# find the nearest peak for each of the positions that meet criteria (this will eventually give us frequency)
	if (opt$stranded=="yes") { 
		n<-nearest(ir, gr_reduce, ignore.strand=FALSE)
	} else { 
		n<-nearest(ir, gr_reduce)
	}
	df<-as.data.frame(gr_reduce)
	dfn<-as.data.frame(table(n))
	# subset the gtf file by gene_biotype of gene.  Make this a requirement.
	txs<-grd[grd$type=="gene"]
	# annotate the reduced pseudo peaks to the newly formed gtf GRanges object
	if (opt$stranded=="yes") { 
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


	#Sanity check 
	#final2<-txs[n2]
	#final<-as.data.frame(final2)
	#rows_to_output<-nrow(final)
	#df<-head(df, rows_to_output)
	#dfn<-head(dfn, rows_to_output)

	#calculate the # of reads overlapping the interval
	bamView<-BamViews(bamname)
	bamRanges(bamView) <- gr_reduce
	aln<-countBam(bamView)
	a<-aln@listData
	adf<-as.data.frame(a)
	names(adf)[6]<-"records"

	output_df<-data.frame(df$seqnames, df$start+opt$interval, df$end-opt$interval, dfn$Freq, (dfn$Freq / (df$width - (2 * opt$interval) + 2)), final$seqnames, final$start, final$end, final$gene_id, final$gene_type, adf$records)
	colnames(output_df)<-c("Chr", "Start", "Stop", "Bases", "% bases meeting criteria", "GTF chr", "GTF start", "GTF stop", "GTF gene_id", "GTF gene_type", "Reads in Expanded Interval")
	if (cflag == "yes") {
		write.csv(file=paste0(opt$output, "_c_supermerge_results.csv"), output_df)
		if (opt$read_type == "paired") {
			cfrip<<-frip(bamname, gr_reduce2, singleEnd=F)
		} else {
			cfrip<<-frip(bamname, gr_reduce2, singleEnd=T)
		}
	} else {
		write.csv(file=paste0(opt$output, "_supermerge_results.csv"), output_df)
		if (opt$read_type == "paired") {
			pfrip<<-frip(bamname, gr_reduce2, singleEnd=F)
		} else {
			pfrip<<-frip(bamname, gr_reduce2, singleEnd=T)
		}
	}
	return(output_df)
}

pfrip<-0
cfrip<-0
cat("Run the main bam file through the function\n")
ip_df<-supermerge(opt$ip_bam, "no")
if (opt$control != "NA") {
	cat("Now run the comparison bam file through the function\n")
	control_df<-supermerge(opt$control, "yes")
}

cat("pfrip is \n")
cat(pfrip)
cat("\ncfrip is \n")
cat(cfrip)

# Draw the graphical output here

cat ("END OF RUN\n")
