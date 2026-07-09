 #!/usr/local/bin Rscript

# library(readr)
library(tidyverse)

options(scipen=999)
vcf_header = read_lines("/data1/papaemme/isabl/home/amallra/repos/germline_svs/scripts/cnmops.vcf.header.txt")

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
#out_dir <- args[2]
sample_cnvs <- args[2]

# function to create VCF
#create_cnmops_vcf <- function(sample_name, out_dir) {
create_cnmops_vcf <- function(sample_name, sample_cnvs) {
  
  cnmops_sample_name = paste0(sample_name, "_cnmops")
  
  # sample_cnmops = read.delim(paste0(out_dir, "/", sample_name, "/", sample_name, "_cohort_cnvs.tsv"), sep = "\t", header = TRUE)
  sample_cnmops = read.delim(sample_cnvs, sep = "\t", header = TRUE)
  sample_cnmops = sample_cnmops %>% 
    dplyr::filter(sampleName == sample_name) %>%
    dplyr::rename(CHROM = seqnames, POS = start) 
  
  sample_cnmops = sample_cnmops %>%
    mutate(ID = paste0("CNV_", CHROM, "_", POS, "_", end)) %>%
    mutate(REF = "N") %>% 
    mutate(ALT = case_when(CN %in% c("CN0", "CN1") ~ "<DEL>", 
                           CN %in% c("CN3", "CN4", "CN5", "CN6", "CN7", "CN8") ~ "<DUP>")) %>% 
    mutate(QUAL = mean) %>%
    mutate(FILTER = "PASS") %>% 
    mutate(CN != "CN2") %>% 
    mutate(INFO = case_when(CN %in% c("CN0", "CN1") ~ paste0("END=", as.numeric(end), ";", "SVLEN=", as.numeric(width), ";", "SVTYPE=", gsub("[<>]", "", ALT), ";", "STRANDS=", "+-"), 
                            CN %in% c("CN3", "CN4", "CN5", "CN6", "CN7", "CN8") ~ paste0("END=", as.numeric(end), ";", "SVLEN=", as.numeric(width), ";", "SVTYPE=", gsub("[<>]", "", ALT), ";", "STRANDS=", "-+"))) %>% 
    mutate(FORMAT = "CN:Median:Mean") %>%
    mutate(!!cnmops_sample_name := case_when(CN %in% c("CN0") ~ paste0(gsub("CN", "", CN),":", median, ":", mean), 
                                             CN %in% c("CN1") ~ paste0(sub("CN", "", CN), ":", median, ":", mean), 
                                             CN %in% c("CN3", "CN4", "CN5", "CN6", "CN7", "CN8") ~ paste0(gsub("CN", "", CN), ":", median, ":", mean))) %>%
    dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, eval(cnmops_sample_name)) %>%
    dplyr::rename(`#CHROM` = CHROM)
  
  #write_lines(vcf_header, 
  #            file = file.path(paste0(out_dir, "/", sample_name, "/", sample_name, ".vcf")), 
  #            append = FALSE)
  #write_tsv(sample_cnmops, 
  #          file = file.path(paste0(out_dir, "/", sample_name, "/", sample_name, ".vcf")), 
  #          col_names = TRUE, 
  #          append = TRUE)

  write_lines(vcf_header, 
              file = file.path(paste0(sample_name, ".vcf")), 
              append = FALSE)
  write_tsv(sample_cnmops, 
            file = file.path(paste0(sample_name, ".vcf")), 
            col_names = TRUE, 
            append = TRUE)
}

# call function
#create_cnmops_vcf(sample_name, out_dir)
create_cnmops_vcf(sample_name, sample_cnvs)

# individual samples
# create_cnmops_vcf(sample_name="I-H-133671-N1-1-D1-1", out_dir="/data1/papaemme/isabl/home/amallra/repos/germline_svs/analysis/cnmops_per_sample")
# create_cnmops_vcf(sample_name="IID_H201573_N01_01_WG01", out_dir="/data1/papaemme/isabl/home/amallra/repos/germline_svs/analysis/cnmops_per_sample")
# create_cnmops_vcf(sample_name="IID_H210166_N01_01_WG01", out_dir="/data1/papaemme/isabl/home/amallra/repos/germline_svs/analysis/cnmops_per_sample")
# create_cnmops_vcf(sample_name="IID_H158916_N01_01_WG01", out_dir="/data1/papaemme/isabl/home/amallra/repos/germline_svs/analysis/cnmops_per_sample")
# create_cnmops_vcf(sample_name="IID_H208287_N01_01_WG01", out_dir="/data1/papaemme/isabl/home/amallra/repos/germline_svs/analysis/cnmops_per_sample")
# create_cnmops_vcf(sample_name="IID_H210868_N01_01_WG01", out_dir="/data1/papaemme/isabl/home/amallra/repos/germline_svs/analysis/cnmops_per_sample")
       
  
