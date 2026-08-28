#!/usr/bin/env python3
"""Generate HTML CNV analysis report for gcnv-nf."""

import sys
import re
import gzip
import base64
from pathlib import Path
from datetime import date


# ── helpers ───────────────────────────────────────────────────────────────────

def count_vcf_calls(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    n = 0
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.strip() and not line.startswith('#'):
                n += 1
    return n


def parse_vcf_supp_vec(path):
    """Return list of SUPP_VEC values from a VCF (handles gzip)."""
    opener = gzip.open if str(path).endswith('.gz') else open
    vecs = []
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            m = re.search(r'SUPP_VEC=([^;]+)', line)
            if m:
                vecs.append(m.group(1))
    return vecs


def info_field(info_str, key):
    m = re.search(rf'{key}=([^;]+)', info_str)
    return m.group(1) if m else ''


def parse_tsv(path):
    with open(path) as fh:
        lines = [l.rstrip('\n') for l in fh if l.strip()]
    if len(lines) < 2:
        return [], []
    headers = lines[0].split('\t')
    return headers, [dict(zip(headers, l.split('\t'))) for l in lines[1:]]


def parse_tsv_raw(path):
    """Return (headers_list, rows_list_of_lists) for index-based access."""
    with open(path) as fh:
        lines = [l.rstrip('\n') for l in fh if l.strip()]
    if not lines:
        return [], []
    return lines[0].split('\t'), [l.split('\t') for l in lines[1:]]


def build_gene_lookup(annotated_tsv_path):
    """
    Read the full annotated TSV and, for every 'split' row, record:
      - the gene name (col 17, 1-indexed)
      - whether that gene has a pathogenic classification (col 53, 1-indexed)
    Returns dict keyed by (chrom, start, end, svtype) ->
      {'genes': set, 'pathogenic_genes': set}
    """
    _, rows = parse_tsv_raw(annotated_tsv_path)
    lookup = {}
    for row in rows:
        if len(row) < 17:
            continue
        ann_mode = row[15].strip() if len(row) > 15 else ''   # col 16
        if ann_mode != 'split':
            continue
        chrom  = row[1].strip()  if len(row) > 1  else ''    # col 2
        start  = row[2].strip()  if len(row) > 2  else ''    # col 3
        end    = row[3].strip()  if len(row) > 3  else ''    # col 4
        svtype = row[5].strip()  if len(row) > 5  else ''    # col 6
        gene   = row[16].strip() if len(row) > 16 else ''    # col 17
        pathcol = row[52].strip() if len(row) > 52 else ''   # col 53
        key = (chrom, start, end, svtype)
        if key not in lookup:
            lookup[key] = {'genes': [], 'pathogenic_genes': []}
        if gene and gene not in lookup[key]['genes']:
            lookup[key]['genes'].append(gene)
        if gene and re.search(r'pathogenic', pathcol, re.IGNORECASE) \
                and gene not in lookup[key]['pathogenic_genes']:
            lookup[key]['pathogenic_genes'].append(gene)
    return lookup


SIZE_BINS = [
    ('<1 kb',      0,          1_000),
    ('1–10 kb',    1_000,      10_000),
    ('10–100 kb',  10_000,     100_000),
    ('100kb–1Mb',  100_000,    1_000_000),
    ('>1 Mb',      1_000_000,  float('inf')),
]


def parse_merged_vcf_sv_info(path):
    """Return list of {svtype, size, supp_vec} from a merged VCF."""
    opener = gzip.open if str(path).endswith('.gz') else open
    svs = []
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            info     = fields[7] if len(fields) > 7 else ''
            svtype   = info_field(info, 'SVTYPE') or 'OTHER'
            supp_vec = info_field(info, 'SUPP_VEC')
            try:
                size = abs(int(info_field(info, 'SVLEN')))
            except (ValueError, TypeError):
                try:
                    end_val = int(info_field(info, 'END'))
                    size = abs(end_val - int(fields[1]))
                except (ValueError, TypeError):
                    size = 0
            svs.append({'svtype': svtype, 'size': size, 'supp_vec': supp_vec})
    return svs


def make_sv_type_table(svs):
    """Return HTML table of SV type counts broken down by caller."""
    from collections import defaultdict
    counts = defaultdict(lambda: {'cnmops_only': 0, 'gatk_only': 0, 'both': 0, 'total': 0})
    for sv in svs:
        t = sv['svtype']
        v = sv['supp_vec']
        counts[t]['total'] += 1
        if v == '10':   counts[t]['cnmops_only'] += 1
        elif v == '01': counts[t]['gatk_only'] += 1
        elif v == '11': counts[t]['both'] += 1
    if not counts:
        return '<p class="muted">No SV data available.</p>'
    rows = ''
    totals = {'cnmops_only': 0, 'gatk_only': 0, 'both': 0, 'total': 0}
    for svtype in sorted(counts):
        c = counts[svtype]
        rows += (f'<tr><td>{svtype}</td><td>{c["cnmops_only"]}</td>'
                 f'<td>{c["gatk_only"]}</td><td>{c["both"]}</td>'
                 f'<td><strong>{c["total"]}</strong></td></tr>')
        for k in totals: totals[k] += c[k]
    rows += (f'<tr class="total"><td><strong>Total</strong></td>'
             f'<td>{totals["cnmops_only"]}</td><td>{totals["gatk_only"]}</td>'
             f'<td>{totals["both"]}</td><td><strong>{totals["total"]}</strong></td></tr>')
    return (f'<table><thead><tr><th>SV Type</th><th>cn.mops only</th>'
            f'<th>GATK only</th><th>Both (consensus)</th><th>Total</th>'
            f'</tr></thead><tbody>{rows}</tbody></table>')


def make_size_dist_chart(svs):
    """Return an HTML CSS bar chart of SV size distribution by type."""
    from collections import defaultdict
    counts = defaultdict(lambda: {'DEL': 0, 'DUP': 0, 'OTHER': 0})
    for sv in svs:
        t = sv['svtype'] if sv['svtype'] in ('DEL', 'DUP') else 'OTHER'
        for label, lo, hi in SIZE_BINS:
            if lo <= sv['size'] < hi:
                counts[label][t] += 1
                break
    max_total = max(
        (counts[label]['DEL'] + counts[label]['DUP'] + counts[label]['OTHER']
         for label, *_ in SIZE_BINS),
        default=1
    ) or 1
    MAX_H = 150
    SV_COLORS = {'DEL': '#e53e3e', 'DUP': '#3182ce', 'OTHER': '#a0aec0'}
    groups = ''
    for label, lo, hi in SIZE_BINS:
        c = counts[label]
        total = c['DEL'] + c['DUP'] + c['OTHER']
        bars = ''
        for t, color in SV_COLORS.items():
            if c[t] > 0:
                h = max(4, int((c[t] / max_total) * MAX_H))
                bars += (f'<div class="bar" style="height:{h}px;background:{color}" '
                         f'title="{t}: {c[t]}"></div>')
        groups += (f'<div class="bar-group">'
                   f'<div class="bars">{bars}</div>'
                   f'<div class="bar-count">{total}</div>'
                   f'<div class="bar-label">{label}</div>'
                   f'</div>')
    legend = ''.join(
        f'<div class="legend-item">'
        f'<div class="legend-dot" style="background:{c}"></div>{t}'
        f'</div>'
        for t, c in SV_COLORS.items()
    )
    return (f'<div class="chart">{groups}</div>'
            f'<div class="chart-legend">{legend}</div>')


def parse_image_kind(path):
    name = Path(path).name
    if name.endswith('_ideogram.png'):
        return 'ideogram'
    if name.endswith('_normalized.png'):
        return 'normalized'
    return 'depth'


def image_src_png(img_path):
    """Return data URL for a PNG path, or the original path on read failure."""
    path = resolve_image_path(img_path)
    try:
        with open(path, 'rb') as fh:
            b64 = base64.b64encode(fh.read()).decode('ascii')
        return f'data:image/png;base64,{b64}'
    except OSError:
        return img_path


def resolve_image_path(img_path):
    """Resolve image path as-is, then by basename in current working directory."""
    p = Path(str(img_path))
    if p.exists():
        return str(p)
    alt = Path(p.name)
    if alt.exists():
        return str(alt)
    return str(p)


def parse_event_image_name(img_path):
    """Parse event PNG naming into sample/chrom/start/end/svtype/kind.

    Supports names like:
      sample_chr_start_end_svtype.png
      sample_chr_start_end_svtype_ideogram.png
      sample_chr_start_end_svtype_normalized.png
    """
    name = Path(img_path).name
    if not name.lower().endswith('.png'):
        return None

    stem = name[:-4]
    kind = 'depth'
    if stem.endswith('_ideogram'):
        stem = stem[:-len('_ideogram')]
        kind = 'ideogram'
    elif stem.endswith('_normalized'):
        stem = stem[:-len('_normalized')]
        kind = 'normalized'

    m = re.match(r'^(?P<sample_chrom>.+)_(?P<start>\d+)_(?P<end>\d+)_(?P<svtype>.+)$', stem)
    if not m:
        return None

    sample_chrom = m.group('sample_chrom')
    if '_' not in sample_chrom:
        return None
    sample, chrom = sample_chrom.rsplit('_', 1)

    return {
        'sample': sample,
        'chrom': chrom,
        'start': m.group('start'),
        'end': m.group('end'),
        'svtype': m.group('svtype'),
        'kind': kind,
    }


def embed_png_html(img_path, title=None):
    stem = Path(img_path).stem
    anchor_id = f'plot-{stem}'
    title_html = f'<p class="plot-title">{title or stem}</p>'
    src = image_src_png(img_path)
    return (
        f'<div class="plot" id="{anchor_id}">'
        f'{title_html}'
        f'<img src="{src}" alt="{stem}" style="width:100%;max-width:760px;display:block;margin:0 auto;">'
        f'</div>'
    )


RANKING_LABEL = {
    '1': 'Benign', '2': 'Likely Benign', '3': 'VOUS',
    '4': 'Likely Pathogenic', '5': 'Pathogenic'
}
CALLER_NAME = ['cn.mops', 'GATK gCNV']

CSS = """
  body{font-family:-apple-system,Arial,sans-serif;max-width:1100px;margin:0 auto;padding:24px;color:#2d3748}
  h1{color:#1a365d;border-bottom:3px solid #3182ce;padding-bottom:8px}
  h2{color:#2b6cb0;margin-top:40px}
  .meta{color:#718096;margin-bottom:24px}
  .boxes{display:flex;gap:16px;flex-wrap:wrap;margin:16px 0}
  .box{background:#ebf8ff;border:1px solid #90cdf4;border-radius:8px;padding:16px 28px;text-align:center;min-width:140px}
  .box .num{font-size:2.2em;font-weight:700;color:#2b6cb0}
  .box .lbl{font-size:.85em;color:#4a5568;margin-top:4px}
  table{border-collapse:collapse;width:100%;margin:12px 0;font-size:.95em}
  th{background:#2b6cb0;color:#fff;padding:10px 14px;text-align:left}
  td{padding:8px 14px;border-bottom:1px solid #e2e8f0}
  tr:nth-child(even) td{background:#f7fafc}
  tr.total td{background:#ebf8ff;border-top:2px solid #90cdf4}
  .rank-4{color:#dd6b20;font-weight:600}
  .rank-5{color:#c53030;font-weight:600}
  .muted{color:#a0aec0;font-style:italic}
  .genes{font-size:.85em;color:#2d3748}
  .gene-tag{display:inline-block;background:#e9d8fd;color:#553c9a;
            border-radius:3px;padding:1px 6px;margin:1px;font-size:.8em}
  .path-tag{display:inline-block;background:#fed7d7;color:#9b2335;
            border-radius:3px;padding:1px 6px;margin:1px;font-size:.8em;font-weight:600}
  a.plot-link{color:#2b6cb0;text-decoration:none;border-bottom:1px dashed #90cdf4}
  a.plot-link:hover{color:#1a365d;border-bottom-style:solid}
  .chart{display:flex;align-items:flex-end;gap:8px;padding:10px 0;
         border-bottom:2px solid #e2e8f0;margin-bottom:4px}
  .bar-group{display:flex;flex-direction:column;align-items:center;gap:4px;min-width:90px}
  .bars{display:flex;align-items:flex-end;gap:3px;height:160px}
  .bar{min-width:22px;min-height:2px;border-radius:3px 3px 0 0;cursor:default}
  .bar:hover{opacity:.75}
  .bar-label{font-size:.8em;color:#4a5568;text-align:center}
  .bar-count{font-size:.8em;color:#718096}
  .chart-legend{display:flex;gap:16px;margin:6px 0 20px}
  .legend-item{display:flex;align-items:center;gap:6px;font-size:.85em}
  .legend-dot{width:12px;height:12px;border-radius:2px}
  .plot{margin:20px 0;border:1px solid #e2e8f0;border-radius:6px;overflow:hidden}
  .plot-title{background:#f7fafc;margin:0;padding:8px 14px;font-size:.85em;
              color:#4a5568;border-bottom:1px solid #e2e8f0}
  .event-plot-block{margin:26px 0;padding:0 0 10px 0;border:1px solid #e2e8f0;border-radius:8px;overflow:hidden}
  .event-plot-block .event-header{background:#f7fafc;padding:10px 14px;border-bottom:1px solid #e2e8f0;font-weight:600;color:#2d3748}
  .image-stack{display:flex;flex-direction:column;gap:16px;padding:18px 18px 8px}
  .image-stack img{display:block;width:100%;max-width:760px;margin:0 auto;border:1px solid #e2e8f0;border-radius:8px;background:#fff}
  footer{margin-top:50px;font-size:.8em;color:#a0aec0;
         border-top:1px solid #e2e8f0;padding-top:12px}
"""


# ── report builder ────────────────────────────────────────────────────────────

def build_report(sample_id, cnmops_vcf, gatk_vcf, merged_vcf, annotated_tsv, priority_tsv, image_files):

    n_cnmops = count_vcf_calls(cnmops_vcf)
    n_gatk   = count_vcf_calls(gatk_vcf)
    n_merged = count_vcf_calls(merged_vcf)

    merged_vecs     = parse_vcf_supp_vec(merged_vcf)
    n_merged_cnmops = sum(1 for v in merged_vecs if len(v) > 0 and v[0] == '1')
    n_merged_gatk   = sum(1 for v in merged_vecs if len(v) > 1 and v[1] == '1')
    n_consensus     = sum(1 for v in merged_vecs if v == '11')

    gene_lookup = build_gene_lookup(annotated_tsv)

    svs        = parse_merged_vcf_sv_info(merged_vcf)
    qc_sv_type = make_sv_type_table(svs)
    qc_size    = make_size_dist_chart(svs)

    _, priority_rows  = parse_tsv(priority_tsv)
    n_priority        = len(priority_rows)
    n_priority_cnmops = sum(
        1 for r in priority_rows
        if info_field(r.get('INFO', ''), 'SUPP_VEC')[0:1] == '1'
    )
    n_priority_gatk = sum(
        1 for r in priority_rows
        if info_field(r.get('INFO', ''), 'SUPP_VEC')[1:2] == '1'
    )

    # summary boxes
    boxes = ''.join(
        f'<div class="box"><div class="num">{n}</div>'
        f'<div class="lbl">{lbl}</div></div>'
        for n, lbl in [
            (n_cnmops,    'cn.mops calls'),
            (n_gatk,      'GATK gCNV calls'),
            (n_merged,    'Merged calls'),
            (n_consensus, 'Consensus calls'),
            (n_priority,  'Priority calls'),
        ]
    )

    # caller breakdown table
    breakdown = f"""
        <tr><td>cn.mops</td><td>{n_cnmops}</td>
            <td>{n_merged_cnmops}</td><td>{n_priority_cnmops}</td></tr>
        <tr><td>GATK gCNV</td><td>{n_gatk}</td>
            <td>{n_merged_gatk}</td><td>{n_priority_gatk}</td></tr>
        <tr class="total">
            <td><strong>Merged total</strong></td><td>—</td>
            <td><strong>{n_merged}</strong></td>
            <td><strong>{n_priority}</strong></td>
        </tr>
    """

    # priority events table
    if priority_rows:
        event_rows = ""
        for r in priority_rows:
            chrom   = r.get('SV chrom', '')
            start   = r.get('SV start', '')
            end     = r.get('SV end', '')
            sv_type = r.get('SV type', '')
            info    = r.get('INFO', '')
            ranking = r.get('AnnotSV ranking', '')
            try:
                size = f"{int(end) - int(start):,} bp"
            except (ValueError, TypeError):
                size = '—'
            supp_vec = info_field(info, 'SUPP_VEC')
            callers  = ', '.join(
                CALLER_NAME[i] for i, c in enumerate(supp_vec) if c == '1'
            )
            rnk_lbl = RANKING_LABEL.get(ranking, 'Unknown')
            rnk_cls = f'rank-{ranking}' if ranking in ('4', '5') else ''

            # gene info from annotated TSV split rows
            key = (chrom.strip(), start.strip(), end.strip(), sv_type.strip())
            gene_info = gene_lookup.get(key, {'genes': [], 'pathogenic_genes': []})
            genes_html = ', '.join(gene_info['genes']) or '—'
            path_genes_html = ', '.join(gene_info['pathogenic_genes']) or '—'

            # anchor link to coverage plot
            plot_anchor = f'event-{sample_id}_{chrom}_{start}_{end}_{sv_type}'
            locus_html = f'<a href="#{plot_anchor}">{chrom}:{start}–{end}</a>'

            event_rows += (
                f'<tr>'
                f'<td>{locus_html}</td>'
                f'<td>{size}</td>'
                f'<td>{sv_type}</td>'
                f'<td>{callers}</td>'
                f'<td class="{rnk_cls}">{ranking} – {rnk_lbl}</td>'
                f'<td class="genes">{genes_html}</td>'
                f'<td class="genes">{path_genes_html}</td>'
                f'</tr>'
            )
        events_html = f"""
        <table>
          <thead><tr>
            <th>Locus</th><th>Size</th><th>Type</th>
            <th>Caller(s)</th><th>AnnotSV Ranking</th>
            <th>Genes</th><th>Pathogenic Genes</th>
          </tr></thead>
          <tbody>{event_rows}</tbody>
        </table>"""
    else:
        events_html = '<p class="muted">No priority events found.</p>'

    # coverage plots grouped by event (ideogram -> depth -> normalized depth)
    event_images = {}
    unmatched_images = []
    for img_path in sorted(image_files):
        parsed = parse_event_image_name(img_path)
        if not parsed:
            unmatched_images.append(img_path)
            continue
        key = (parsed['chrom'], parsed['start'], parsed['end'], parsed['svtype'])
        kind = parsed['kind']
        event_images.setdefault(key, {'chrom': parsed['chrom'], 'start': int(parsed['start']), 'end': int(parsed['end']), 'svtype': parsed['svtype'], 'images': {}})
        event_images[key]['images'][kind] = img_path

    if event_images:
        blocks = []
        for key in sorted(event_images, key=lambda k: (event_images[k]['chrom'], event_images[k]['start'], event_images[k]['end'], event_images[k]['svtype'])):
            event = event_images[key]
            chrom = event['chrom']
            start = event['start']
            end = event['end']
            svtype = event['svtype']
            anchor_id = f'event-{sample_id}_{chrom}_{start}_{end}_{svtype}'
            image_stack = []
            for kind, label in [('ideogram', 'Ideogram'), ('depth', 'Depth Plot'), ('normalized', 'Normalized Depth Plot')]:
                img_path = event['images'].get(kind)
                if img_path:
                    img_src = image_src_png(img_path)
                    image_stack.append(f'<img src="{img_src}" alt="{label}">')
            if image_stack:
                blocks.append(f'<div class="event-plot-block" id="{anchor_id}"><div class="event-header">{chrom}:{start}–{end} {svtype}</div><div class="image-stack">{"".join(image_stack)}</div></div>')
        if unmatched_images:
            fallback_imgs = ''.join(
                f'<img src="{image_src_png(p)}" alt="{Path(p).name}">'
                for p in unmatched_images
            )
            blocks.append(
                '<div class="event-plot-block">'
                '<div class="event-header">Additional Plots (unmatched naming)</div>'
                f'<div class="image-stack">{fallback_imgs}</div>'
                '</div>'
            )
        plots_html = '\n'.join(blocks) if blocks else '<p class="muted">No coverage plots generated.</p>'
    else:
        if unmatched_images:
            plots_html = ''.join(
                f'<img src="{image_src_png(p)}" alt="{Path(p).name}">'
                for p in unmatched_images
            )
        else:
            plots_html = '<p class="muted">No coverage plots generated.</p>'

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>CNV Report – {sample_id}</title>
<style>{CSS}</style>
</head>
<body>
<h1>CNV Analysis Report</h1>
<p class="meta">
  <strong>Sample:</strong> {sample_id} &nbsp;|&nbsp;
  <strong>Date:</strong> {date.today()}
</p>

<h2>Summary</h2>
<div class="boxes">{boxes}</div>

<h2>Quality Control</h2>
<h3>SV Type Breakdown</h3>
{qc_sv_type}
<h3>Size Distribution</h3>
{qc_size}

<h2>Calls per Caller</h2>
<table>
  <thead><tr>
    <th>Caller</th><th>Total Calls</th>
    <th>Calls in Merged Set</th><th>Priority Calls</th>
  </tr></thead>
  <tbody>{breakdown}</tbody>
</table>

<h2>Priority Events</h2>
{events_html}

<h2>Coverage Plots</h2>
{plots_html}

<footer>Generated by gcnv-nf &mdash; {date.today()}</footer>
</body>
</html>"""


# ── entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    if len(sys.argv) < 7:
        print(
            'Usage: create_report.py <sample_id> <cnmops_vcf> <gatk_vcf> '
            '<merged_vcf> <annotated_tsv> <priority_tsv> [image_files ...]',
            file=sys.stderr
        )
        sys.exit(1)

    sample_id     = sys.argv[1]
    cnmops_vcf    = sys.argv[2]
    gatk_vcf      = sys.argv[3]
    merged_vcf    = sys.argv[4]
    annotated_tsv = sys.argv[5]
    priority_tsv  = sys.argv[6]
    raw_image_args = sys.argv[7:]
    image_files = []
    for token in raw_image_args:
        # Handle list-like argument rendering: [a.png, b.png]
        for part in str(token).split(','):
            cleaned = part.strip().strip('[]').strip('"\'')
            if cleaned:
                image_files.append(cleaned)

    html = build_report(
        sample_id, cnmops_vcf, gatk_vcf, merged_vcf, annotated_tsv, priority_tsv, image_files
    )

    out = f'{sample_id}.cnv_report.html'
    with open(out, 'w') as fh:
        fh.write(html)
    print(f'Report written: {out}')
