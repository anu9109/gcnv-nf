#!/usr/bin/env Rscript
library(ggplot2)
library(data.table) 
library(tidyverse)

# plotting function
plot_genome_cov <- function(depth_file, plot_chr, plot_start, plot_end, plot_label, locus_label, plot_savepath) {
  
  # read depth data
  depth_data = fread(depth_file, header = FALSE)
  setnames(depth_data, c("chromosome", "start", "stop", "depth"))
  
  # filter to chr of interest
  depth_chr <- depth_data[chromosome == plot_chr]
  
  # create 10,000 bp window groups
  depth_chr[, window := floor(start / 100000)]
  
  # calculate average depth per 10,000 bp window
  avg_depth <- depth_chr[, .(mean_start = mean(start), mean_depth = mean(depth)), by = window]
  
  # plot
  plot_pad = (plot_end - plot_start) 
  viz = ggplot() +
    geom_point(data = depth_chr, aes(x = start, y = as.numeric(depth)), alpha = 0.4, color = "darkblue", size = 0.5) +
    geom_rect(aes(xmin = plot_start, xmax = plot_end, ymin = 90, ymax = 100), fill = "darkred", alpha = 0.5) + 
    annotate("text", x = plot_start + (plot_pad*0.4), y = 95, label = locus_label) + 
    labs(
      title = plot_label,
      x = "Genomic Position (bp)",
      y = "Read Depth"
    ) +
    scale_x_continuous(
      limits = c(plot_start - 2*plot_pad, plot_end + 2*plot_pad),
      labels = scales::label_comma()  # avoids scientific notation
    ) + 
    scale_y_continuous(limits = c(-5, max(avg_depth$mean_depth, 110)+5)) + 
    theme_minimal()
  
  # save plot to file
  pdf(file = plot_savepath, width = 10, height = 5)
  print(viz)
  dev.off()
}

# define arguments
args         = commandArgs(trailingOnly = TRUE)
sample_id    = args[1]
depth_file   = args[2]
priority_tsv = args[3]

# read priority TSV file
tsv = fread(priority_tsv, sep = "\t", header = TRUE)

ranking_labels = c("1" = "Benign", "2" = "Likely Benign", "3" = "VOUS",
                   "4" = "Likely Pathogenic", "5" = "Pathogenic")

# create plots
for (i in seq_len(nrow(tsv))) {
    row           = tsv[i, ]
    chr           = as.character(row[["SV chrom"]])
    start         = as.integer(row[["SV start"]])
    end           = as.integer(row[["SV end"]])
    sv_type       = as.character(row[["SV type"]])
    info          = as.character(row[["INFO"]])
    ranking       = as.integer(row[["AnnotSV ranking"]])

    supp_match    = regmatches(info, regexpr("SUPP=[0-9]+", info))
    supp          = if (length(supp_match) > 0) as.integer(sub("SUPP=", "", supp_match)) else 0

    ranking_label = ifelse(!is.na(ranking_labels[as.character(ranking)]),
                           ranking_labels[as.character(ranking)], "Unknown")

    plot_label    = sprintf("%s | %s:%d_%d %s | SUPPORTED BY %d CALLERS | RANKING: %s",
                            sample_id, chr, start, end, sv_type, supp, ranking_label)
    plot_savepath = sprintf("%s_%s_%d_%d_%s.pdf", sample_id, chr, start, end, sv_type)

    plot_genome_cov(
        depth_file    = depth_file,
        plot_chr      = chr,
        plot_start    = start,
        plot_end      = end,
        plot_label    = plot_label,
        locus_label   = sv_type,
        plot_savepath = plot_savepath
    )
}
