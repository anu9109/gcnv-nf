#!/usr/bin/env Rscript
library(ggplot2)
library(data.table) 
library(tidyverse)

# ── raw depth plot ────────────────────────────────────────────────────────────
plot_genome_cov = function(depth_file, plot_chr, plot_start, plot_end, plot_label, locus_label, plot_savepath) {
  
  # read depth data
  depth_data = read.table(gzfile(depth_file), header = FALSE, sep = "\t")
  setDT(depth_data)

  # depth_data = fread(depth_file, header = FALSE)
  setnames(depth_data, c("chromosome", "start", "stop", "depth"))
  
  # filter to chr of interest
  depth_chr = depth_data[chromosome == plot_chr]
  
  # create 10,000 bp window groups
  depth_chr[, window := floor(start / 100000)]
  
  # calculate average depth per 10,000 bp window
  avg_depth = depth_chr[, .(mean_start = mean(start), mean_depth = mean(depth)), by = window]
  
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
    scale_y_continuous(limits = c(-5, max(avg_depth$mean_depth, 150)+5)) + 
    theme_minimal()
  
  # save plot to file
  pdf(file = plot_savepath, width = 10, height = 5)
  print(viz)
  dev.off()
}


# ── GC + mappability normalised depth plot ────────────────────────────────────
plot_genome_cov_normalized = function(depth_file, gc_file, map_file,
                                       plot_chr, plot_start, plot_end,
                                       plot_label, locus_label, plot_savepath) {

  # ─ 1. read depth data ─────────────────────────────────────────────
  depth_data = read.table(gzfile(depth_file), header = FALSE, sep = "\t")
  setDT(depth_data)
  setnames(depth_data, c("chromosome", "start", "stop", "depth"))
  depth_chr = depth_data[chromosome == plot_chr]

  # ─ 2. read GC content (bedtools nuc output: chr, start, end, pct_at, pct_gc, ...) ─
  gc_data = fread(gc_file, skip = 1, header = FALSE, sep = "\t",
                   select = c(1L, 2L, 5L))
  setnames(gc_data, c("chromosome", "start", "gc"))
  gc_data[, chromosome := as.character(chromosome)]
  gc_chr = gc_data[chromosome == plot_chr]

  # ─ 3. read mappability BED (high-mappability intervals, gzipped) ──────────
  map_data = fread(cmd = paste("zcat", shQuote(map_file)),
                    header = FALSE, sep = "\t", select = c(1L, 2L, 3L))
  setnames(map_data, c("chromosome", "map_start", "map_end"))
  map_data[, chromosome := as.character(chromosome)]
  map_chr = map_data[chromosome == plot_chr]

  # ─ 4. mappability filter via range overlap (foverlaps) ────────────────
  setkey(map_chr, chromosome, map_start, map_end)
  depth_chr[, stop2 := stop]
  setkey(depth_chr, chromosome, start, stop2)
  depth_mappable = foverlaps(
    depth_chr, map_chr,
    by.x = c("chromosome", "start", "stop2"),
    by.y = c("chromosome", "map_start", "map_end"),
    type = "any", nomatch = NA
  )[!is.na(map_start), .(chromosome, start, stop, depth)]

  # ─ 5. join per-bin GC content by chromosome + start position ─────────
  merged = merge(depth_mappable, gc_chr[, .(start, gc)], by = "start")

  if (nrow(merged) < 10) {
    message(sprintf("  Skipping normalised plot for %s:%d-%d — insufficient bins after mappability filter",
                    plot_chr, plot_start, plot_end))
    return(invisible(NULL))
  }

  # ─ 6. LOESS GC normalisation (chromosome-wide) ───────────────────────
  loess_fit = loess(depth ~ gc, data = merged, span = 0.3, na.action = na.omit, family="symmetric")
  merged[, pred := pmax(predict(loess_fit, newdata = merged), 1e-6)]
  median_depth = median(merged$depth, na.rm = TRUE)
  merged[, depth_norm := (depth / pred) * median_depth]
  merged[depth_norm < 0, depth_norm := 0]

  # ─ 7. plot ──────────────────────────────────────────────────────────
  plot_pad = plot_end - plot_start
  avg_norm = merged[, .(mean_depth_norm = mean(depth_norm, na.rm = TRUE)),
                     by = .(window = floor(start / 100000))]

  viz = ggplot() +
    geom_point(data = merged, aes(x = start, y = depth_norm),
               alpha = 0.4, color = "darkgreen", size = 0.5) +
    geom_rect(aes(xmin = plot_start, xmax = plot_end, ymin = 90, ymax = 100),
              fill = "darkred", alpha = 0.5) +
    annotate("text", x = plot_start + (plot_pad * 0.4), y = 95, label = locus_label) +
    labs(
      title = paste0(plot_label, " [GC & Mappability Normalised]"),
      x     = "Genomic Position (bp)",
      y     = "Normalised Read Depth"
    ) +
    scale_x_continuous(
      limits = c(plot_start - 2*plot_pad, plot_end + 2*plot_pad),
      labels = scales::label_comma()
    ) +
    scale_y_continuous(limits = c(-5, max(avg_norm$mean_depth_norm, 150) + 5)) +
    theme_minimal()

  pdf(file = plot_savepath, width = 10, height = 5)
  print(viz)
  dev.off()
}

# define arguments
args         = commandArgs(trailingOnly = TRUE)
sample_id    = args[1]
depth_file   = args[2]
priority_tsv = args[3]
gc_file      = args[4]
map_file     = args[5]

run_normalized = !is.na(gc_file) && !is.na(map_file) && file.exists(gc_file) && file.exists(map_file)

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

    supp_vec_match = regmatches(info, regexpr("SUPP_VEC=[^;]+", info))
    supp_vec       = if (length(supp_vec_match) > 0) sub("SUPP_VEC=", "", supp_vec_match) else ""
    caller_names   = c("cn.mops", "GATK gCNV")
    callers        = paste(caller_names[strsplit(supp_vec, "")[[1]] == "1"], collapse = ", ")
    if (callers == "") callers = "unknown"

    ranking_label = ifelse(!is.na(ranking_labels[as.character(ranking)]),
                           ranking_labels[as.character(ranking)], "Unknown")

    plot_label    = sprintf("%s | %s:%d_%d %s | %s | RANKING: %s",
                            sample_id, chr, start, end, sv_type, callers, ranking_label)
    plot_savepath = sprintf("%s_%s_%d_%d_%s.pdf", sample_id, chr, start, end, sv_type)

    message(sprintf("Plotting event %d/%d: %s:%d-%d %s (ranking: %s, SUPP: %d)", i, nrow(tsv), chr, start, end, sv_type, ranking_label, supp))

    plot_genome_cov(
        depth_file    = depth_file,
        plot_chr      = chr,
        plot_start    = start,
        plot_end      = end,
        plot_label    = plot_label,
        locus_label   = sv_type,
        plot_savepath = plot_savepath
    )

    if (run_normalized) {
        norm_savepath = sprintf("%s_%s_%d_%d_%s_normalized.pdf", sample_id, chr, start, end, sv_type)
        message(sprintf("  -> GC+mappability normalised plot: %s", norm_savepath))
        plot_genome_cov_normalized(
            depth_file    = depth_file,
            gc_file       = gc_file,
            map_file      = map_file,
            plot_chr      = chr,
            plot_start    = start,
            plot_end      = end,
            plot_label    = plot_label,
            locus_label   = sv_type,
            plot_savepath = norm_savepath
        )
    }
}
