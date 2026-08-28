#!/usr/bin/env Rscript
library(ggplot2)
library(data.table) 
library(tidyverse)

build_key_from_input = function(build_label) {
  b = gsub("[^a-z0-9]", "", tolower(trimws(as.character(build_label))))
  if (b %in% c("grch37", "hg19")) return("grch37")
  if (b %in% c("grch38", "hg38")) return("grch38")
  stop(sprintf("Unsupported genome build: %s. Supported: GRCh37/hg19 or GRCh38/hg38", build_label))
}

load_genome_reference = function(build_label, data_dir) {
  key = build_key_from_input(build_label)
  lengths_path = file.path(data_dir, sprintf("%s_chrom_lengths.tsv", key))
  centromeres_path = file.path(data_dir, sprintf("%s_centromeres.tsv", key))

  if (!file.exists(lengths_path)) {
    stop(sprintf("Chromosome lengths file not found: %s", lengths_path))
  }
  if (!file.exists(centromeres_path)) {
    stop(sprintf("Centromere file not found: %s", centromeres_path))
  }

  lengths_dt = fread(lengths_path, sep = "\t", header = TRUE)
  cent_dt = fread(centromeres_path, sep = "\t", header = TRUE)

  req_lengths = c("chrom", "length")
  req_cent = c("chrom", "start", "end")
  if (!all(req_lengths %in% names(lengths_dt))) {
    stop(sprintf("Invalid lengths file schema in %s. Expected columns: chrom,length", lengths_path))
  }
  if (!all(req_cent %in% names(cent_dt))) {
    stop(sprintf("Invalid centromere file schema in %s. Expected columns: chrom,start,end", centromeres_path))
  }

  lengths_dt[, chrom := toupper(gsub("^chr", "", as.character(chrom), ignore.case = TRUE))]
  cent_dt[, chrom := toupper(gsub("^chr", "", as.character(chrom), ignore.case = TRUE))]

  chr_lengths = setNames(as.numeric(lengths_dt$length), lengths_dt$chrom)
  chr_centromeres = setNames(
    lapply(seq_len(nrow(cent_dt)), function(i) c(as.numeric(cent_dt$start[i]), as.numeric(cent_dt$end[i]))),
    cent_dt$chrom
  )

  list(key = key, chr_lengths = chr_lengths, chr_centromeres = chr_centromeres)
}

normalize_chr = function(chr_label) {
  chr = gsub("^chr", "", as.character(chr_label), ignore.case = TRUE)
  toupper(chr)
}

default_data_dir = function() {
  cmd_args = commandArgs(trailingOnly = FALSE)
  script_arg = cmd_args[grep("^--file=", cmd_args)]
  if (length(script_arg) > 0) {
    script_path = sub("^--file=", "", script_arg[1])
    return(normalizePath(file.path(dirname(script_path), "..", "data"), mustWork = FALSE))
  }
  file.path(getwd(), "data")
}

plot_chromosome_ideogram = function(plot_chr, plot_start, plot_end, plot_label, locus_label, plot_savepath) {
  chr_key = normalize_chr(plot_chr)
  chr_len = CHR_LENGTHS[[chr_key]]
  if (is.na(chr_len) || is.null(chr_len)) {
    message(sprintf("  Skipping ideogram for %s: no chromosome length available", plot_chr))
    return(invisible(NULL))
  }

  cen = CHR_CENTROMERES[[chr_key]]
  if (is.null(cen) || length(cen) != 2) {
    cen = c(chr_len * 0.45, chr_len * 0.55)
  }
  cen_start = max(1, min(as.numeric(cen[1]), chr_len))
  cen_end = max(cen_start + 1, min(as.numeric(cen[2]), chr_len))
  cen_mid = (cen_start + cen_end) / 2

  event_start = max(1, min(plot_start, plot_end))
  event_end   = max(1, max(plot_start, plot_end))
  event_start = min(event_start, chr_len)
  event_end   = min(event_end, chr_len)

  ideogram_df = data.frame(
    xmin = c(0, cen_end),
    xmax = c(cen_start, chr_len),
    arm = c("p", "q")
  )

  centromere_top = data.frame(
    x = c(cen_start, cen_mid, cen_end),
    y = c(0.35, 0, 0.35)
  )
  centromere_bottom = data.frame(
    x = c(cen_start, cen_mid, cen_end),
    y = c(-0.35, 0, -0.35)
  )

  viz = ggplot() +
    geom_rect(
      data = ideogram_df,
      aes(xmin = xmin, xmax = xmax, ymin = -0.35, ymax = 0.35, fill = arm),
      color = "black", linewidth = 0.25
    ) +
    geom_polygon(
      data = centromere_top,
      aes(x = x, y = y),
      fill = "grey50", color = "black", linewidth = 0.2
    ) +
    geom_polygon(
      data = centromere_bottom,
      aes(x = x, y = y),
      fill = "grey50", color = "black", linewidth = 0.2
    ) +
    geom_rect(aes(xmin = event_start, xmax = event_end, ymin = -0.28, ymax = 0.28),
              fill = "firebrick2", alpha = 0.8, color = "darkred", linewidth = 0.2) +
    geom_vline(xintercept = event_start, color = "darkred", linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = event_end, color = "darkred", linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = chr_len * 0.2, y = 0.6, label = "p-arm", size = 3.2, color = "grey20") +
    annotate("text", x = chr_len * 0.8, y = 0.6, label = "q-arm", size = 3.2, color = "grey20") +
    annotate("text", x = cen_mid, y = -0.6, label = "centromere", size = 3.0, color = "grey30") +
    annotate("text", x = chr_len * 0.5, y = 0.82, label = paste0("Chr ", chr_key), size = 4) +
    labs(
      title = paste0("Ideogram: ", plot_label),
      x = "Chromosome position (bp)",
      y = ""
    ) +
    scale_x_continuous(limits = c(0, chr_len), labels = scales::label_comma()) +
    scale_fill_manual(values = c("p" = "grey86", "q" = "grey78"), guide = "none") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )

  png(filename = plot_savepath, width = 10, height = 2.5, units = "in", res = 150)
  print(viz)
  dev.off()
}

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
  y_upper = max(avg_depth$mean_depth, 150)

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
    scale_y_continuous(limits = c(-5, y_upper + 5)) + 
    theme_minimal()
  
  png(filename = plot_savepath, width = 10, height = 5, units = "in", res = 150)
  print(viz)
  dev.off()

  invisible(y_upper + 5)
}


# ── GC + mappability normalised depth plot ────────────────────────────────────
plot_genome_cov_normalized = function(depth_file, gc_file, map_file,
                                       plot_chr, plot_start, plot_end,
                                       plot_label, locus_label, plot_savepath,
                                       y_upper_limit = NULL) {

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
    scale_y_continuous(
      limits = c(
        -5,
        ifelse(is.null(y_upper_limit), max(avg_norm$mean_depth_norm, 150) + 5, y_upper_limit)
      )
    ) +
    theme_minimal()

  png(filename = plot_savepath, width = 10, height = 5, units = "in", res = 150)
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
genome_build = ifelse(length(args) >= 6 && nzchar(args[6]), args[6], "GRCh37")
data_dir = default_data_dir()

genome_ref = load_genome_reference(genome_build, data_dir)
CHR_LENGTHS = genome_ref$chr_lengths
CHR_CENTROMERES = genome_ref$chr_centromeres
message(sprintf("Using genome build %s from %s", toupper(genome_ref$key), data_dir))

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
    plot_savepath = sprintf("%s_%s_%d_%d_%s.png", sample_id, chr, start, end, sv_type)
    ideogram_path = sprintf("%s_%s_%d_%d_%s_ideogram.png", sample_id, chr, start, end, sv_type)

    message(sprintf("Plotting event %d/%d: %s:%d-%d %s (ranking: %s, SUPP: %d)", i, nrow(tsv), chr, start, end, sv_type, ranking_label, supp))

    plot_chromosome_ideogram(
        plot_chr      = chr,
        plot_start    = start,
        plot_end      = end,
        plot_label    = plot_label,
        locus_label   = sv_type,
        plot_savepath = ideogram_path
    )

    depth_y_limit = plot_genome_cov(
        depth_file    = depth_file,
        plot_chr      = chr,
        plot_start    = start,
        plot_end      = end,
        plot_label    = plot_label,
        locus_label   = sv_type,
        plot_savepath = plot_savepath
    )

    if (run_normalized) {
        norm_savepath = sprintf("%s_%s_%d_%d_%s_normalized.png", sample_id, chr, start, end, sv_type)
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
            plot_savepath = norm_savepath,
            y_upper_limit = depth_y_limit
        )
    }
}
