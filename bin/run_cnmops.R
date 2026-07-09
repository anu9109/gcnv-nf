#!/usr/local/bin Rscript

library(cn.mops)

# parse arguments
args = commandArgs(trailingOnly = TRUE)
samples_bams_file = args[1]
nparallel = args[2]
#outpath = args[3]
outprefix = args[3]

# samples_bams_file = "/data1/papaemme/isabl/home/amallra/repos/germline_svs/data/nbl.bams.txt"
# nparallel = 7
# outpath = "/data1/papaemme/isabl/home/amallra/repos/germline_svs/analysis/cnmops/"
# outprefix = "01_NBL"

samples_bams = read.delim(samples_bams_file, sep = "\t", header = FALSE, col.names = c("samples", "bams"))
BAMFiles <- samples_bams$bams

# get read counts
message("1. Getting read counts from BAM files:", Sys.time())
bamDataRanges <- getReadCountsFromBAM(BAMFiles, 
                                      sampleNames = samples_bams$samples, 
                                      refSeqName = c(as.character(1:22), "X", "Y"),
                                      WL = 1000, 
                                      parallel = nparallel)

# run cn.mops 
message("2. Running cn.mops..", Sys.time())
res <- cn.mops(bamDataRanges)
result <- calcIntegerCopyNumbers(res)

# outputs
message("3. Compiling results..", Sys.time())
segm <- as.data.frame(segmentation(result))
CNVs <- as.data.frame(cnvs(result))
CNVRegions <- as.data.frame(cnvr(result))

# write to file
#message("4. Writing to path:", outpath, Sys.time())
#write.table(segm,file=paste0(outpath, "/", outprefix, "/", outprefix, "_cohort_segmentation.tsv"), sep = "\t", quote = FALSE)
#write.table(CNVs,file=paste0(outpath, "/", outprefix, "/", outprefix, "_cohort_cnvs.tsv"), sep = "\t", quote = FALSE)
#write.table(CNVRegions,file=paste0(outpath, "/", outprefix, "/", outprefix, "_cohort_cnvr.tsv"), sep = "\t", quote = FALSE)

message("4. Writing output files.", Sys.time())
write.table(segm,file=paste0(outprefix, "_cohort_segmentation.tsv"), sep = "\t", quote = FALSE)
write.table(CNVs,file=paste0(outprefix, "_cohort_cnvs.tsv"), sep = "\t", quote = FALSE)
write.table(CNVRegions,file=paste0(outprefix, "_cohort_cnvr.tsv"), sep = "\t", quote = FALSE)

message("Done!", Sys.time())
