#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Interactive chromosome visualization with tandem repeat density and telomere status.
Generates standalone HTML with CSS/JS — no Plotly or matplotlib dependencies.
"""

import csv
import json
import logging
import os
import sys

logger = logging.getLogger(__name__)

BIN_SIZE = 100_000  # 100 kb bins for density


def load_fai(fai_path):
    """Load chromosome names and lengths from .fai index."""
    chroms = []
    with open(fai_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chroms.append({"name": parts[0], "length": int(parts[1])})
    return chroms


def load_telomeres(tsv_path):
    """Load telomere check results."""
    telo = {}
    with open(tsv_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('chr\t'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 11:
                telo[parts[0]] = {
                    "left_status": parts[5],
                    "left_density": float(parts[4]),
                    "right_status": parts[9],
                    "right_density": float(parts[8]),
                    "t2t": parts[10] == "T2T",
                }
    return telo


def load_its(bed_path):
    """Load interstitial telomere sites."""
    its = {}
    if not os.path.exists(bed_path):
        return its
    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                its.setdefault(chrom, []).append({
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "count": int(parts[3]),
                })
    return its


def bin_repeats(sat_path, chroms, bin_size=BIN_SIZE):
    """Bin tandem repeats into density tracks per chromosome.

    Returns dict: {chr_name: {bin_index: {"micro": bp, "complex": bp, "total": bp}}}
    """
    # Initialize bins
    bins = {}
    for c in chroms:
        n_bins = (c["length"] // bin_size) + 1
        bins[c["name"]] = [{"micro": 0, "complex": 0, "total": 0} for _ in range(n_bins)]

    csv.field_size_limit(sys.maxsize)

    with open(sat_path, 'r') as f:
        reader = csv.DictReader(
            (line for line in f if not line.startswith('#')),
            delimiter='\t'
        )
        for row in reader:
            chrom = row.get("trf_head", "")
            if chrom not in bins:
                continue
            try:
                start = int(row.get("trf_l_ind", 0))
                end = int(row.get("trf_r_ind", 0))
                period = int(row.get("trf_period", 0))
            except (ValueError, TypeError):
                continue

            arr_len = end - start
            category = "micro" if period < 10 else "complex"

            start_bin = start // bin_size
            end_bin = end // bin_size

            for b in range(start_bin, min(end_bin + 1, len(bins[chrom]))):
                b_start = b * bin_size
                b_end = (b + 1) * bin_size
                overlap = min(end, b_end) - max(start, b_start)
                if overlap > 0:
                    bins[chrom][b][category] += overlap
                    bins[chrom][b]["total"] += overlap

    return bins


def generate_chromosome_html(chroms, bins, telomeres, its, assembly_name="", output_path=None):
    """Generate the full chromosome visualization HTML."""

    # Sort chromosomes intelligently
    def chr_sort_key(c):
        name = c["name"].replace("chr", "").replace("Chr", "")
        if name.isdigit():
            return (0, int(name), "")
        return (1, 0, name)

    chroms_sorted = sorted(chroms, key=chr_sort_key)
    max_len = max(c["length"] for c in chroms_sorted) if chroms_sorted else 1

    # Build chromosome data for JS
    chrom_data = []
    for c in chroms_sorted:
        name = c["name"]
        telo = telomeres.get(name, {})
        chrom_its = its.get(name, [])

        # Compress bins to density values (0-1)
        chr_bins = bins.get(name, [])
        density = []
        for b in chr_bins:
            d = min(b["total"] / BIN_SIZE, 1.0) if BIN_SIZE > 0 else 0
            micro_d = min(b["micro"] / BIN_SIZE, 1.0) if BIN_SIZE > 0 else 0
            complex_d = min(b["complex"] / BIN_SIZE, 1.0) if BIN_SIZE > 0 else 0
            density.append([round(micro_d, 3), round(complex_d, 3)])

        chrom_data.append({
            "name": name,
            "length": c["length"],
            "pct": round(c["length"] / max_len * 100, 2),
            "telo_left": telo.get("left_status", "UNKNOWN"),
            "telo_right": telo.get("right_status", "UNKNOWN"),
            "t2t": telo.get("t2t", False),
            "its": chrom_its,
            "density": density,
        })

    data_json = json.dumps(chrom_data)
    bin_size_kb = BIN_SIZE // 1000

    html = _build_html(data_json, bin_size_kb, assembly_name, max_len)

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        logger.info(f"Chromosome visualization: {output_path}")

    return html


def _build_html(data_json, bin_size_kb, assembly_name, max_len):
    return f'''<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Chromosomes — {assembly_name}</title>
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;500&family=Playfair+Display:wght@700;900&family=Source+Sans+3:wght@300;400;600&display=swap" rel="stylesheet">
<style>
  :root, [data-theme="light"] {{
    --bg: #f8f6f1; --card: #fff; --border: #e2ddd5;
    --text: #1a1a1a; --text2: #555; --muted: #888;
    --cyan: #0891b2; --emerald: #059669; --rose: #be123c;
    --violet: #7c3aed; --amber: #b45309;
    --shadow: rgba(0,0,0,.08);
    --chr-bg: #e8e5df; --chr-border: #d4d0c8;
    --micro-color: #0891b2; --complex-color: #7c3aed;
    --telo-ok: #059669; --telo-partial: #b45309; --telo-absent: #be123c;
    --its-color: #db2777;
  }}
  [data-theme="dark"] {{
    --bg: #0a0e17; --card: #111827; --border: #1e293b;
    --text: #e2e8f0; --text2: #94a3b8; --muted: #64748b;
    --cyan: #22d3ee; --emerald: #34d399; --rose: #fb7185;
    --violet: #a78bfa; --amber: #fbbf24;
    --shadow: rgba(0,0,0,.3);
    --chr-bg: #1e293b; --chr-border: #334155;
    --micro-color: #22d3ee; --complex-color: #a78bfa;
    --telo-ok: #34d399; --telo-partial: #fbbf24; --telo-absent: #fb7185;
    --its-color: #f472b6;
  }}
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{ background:var(--bg); font-family:'Source Sans 3',sans-serif; color:var(--text); padding:40px 24px; }}

  .container {{ max-width:1200px; margin:0 auto; }}

  .header {{
    text-align:center; margin-bottom:48px;
    opacity:0; animation:fadeUp .6s ease-out forwards;
  }}
  .header-label {{
    font-family:'JetBrains Mono',monospace; font-size:10px; font-weight:500;
    letter-spacing:3px; text-transform:uppercase; color:var(--cyan); margin-bottom:12px;
  }}
  .header h1 {{
    font-family:'Playfair Display',serif; font-size:clamp(28px,4vw,42px); font-weight:900;
    background:linear-gradient(135deg,var(--text) 0%,var(--cyan) 50%,var(--emerald) 100%);
    -webkit-background-clip:text; -webkit-text-fill-color:transparent; background-clip:text;
  }}

  .toggle {{
    position:fixed; top:20px; right:20px; z-index:100;
    padding:5px 12px; border-radius:20px;
    background:var(--card); border:1px solid var(--border); cursor:pointer;
    font-family:'JetBrains Mono',monospace; font-size:10px; color:var(--muted);
    box-shadow:0 2px 8px var(--shadow); user-select:none;
    display:flex; align-items:center; gap:6px;
  }}
  .toggle:hover {{ transform:translateY(-1px); }}
  .toggle-icon {{ font-size:14px; }}

  .legend {{
    display:flex; gap:20px; justify-content:center; flex-wrap:wrap;
    margin-bottom:32px; opacity:0; animation:fadeUp .5s ease-out .1s forwards;
  }}
  .legend-item {{
    display:flex; align-items:center; gap:6px;
    font-family:'JetBrains Mono',monospace; font-size:11px; color:var(--text2);
  }}
  .legend-dot {{ width:10px; height:10px; border-radius:2px; }}

  .chr-list {{
    display:flex; flex-direction:column; gap:6px;
  }}

  .chr-row {{
    display:grid; grid-template-columns:80px 1fr 60px; align-items:center; gap:12px;
    opacity:0; transform:translateY(12px);
    animation:fadeUp .4s ease-out forwards;
    padding:4px 0;
  }}

  .chr-label {{
    font-family:'JetBrains Mono',monospace; font-size:12px; font-weight:500;
    color:var(--text); text-align:right; white-space:nowrap;
  }}

  .chr-bar-wrap {{
    position:relative; height:22px; cursor:pointer;
  }}

  .chr-bar {{
    position:absolute; top:0; left:0; height:100%;
    background:var(--chr-bg); border:1px solid var(--chr-border);
    border-radius:11px; overflow:hidden;
    transition:box-shadow .2s;
  }}
  .chr-bar:hover {{
    box-shadow:0 2px 12px var(--shadow);
  }}

  .chr-density {{
    position:absolute; top:0; left:0; height:100%;
  }}

  .density-bin {{
    position:absolute; top:0; height:100%;
  }}

  .telo-cap {{
    position:absolute; top:-1px; width:8px; height:24px;
    border-radius:4px; z-index:2;
    transition:transform .2s;
  }}
  .telo-cap:hover {{ transform:scaleY(1.3); }}
  .telo-cap.left {{ left:-1px; border-top-right-radius:0; border-bottom-right-radius:0; }}
  .telo-cap.right {{ right:-1px; border-top-left-radius:0; border-bottom-left-radius:0; }}
  .telo-cap.PRESENT {{ background:var(--telo-ok); }}
  .telo-cap.PARTIAL {{ background:var(--telo-partial); }}
  .telo-cap.ABSENT {{ background:var(--telo-absent); }}
  .telo-cap.UNKNOWN {{ background:var(--muted); }}

  .its-mark {{
    position:absolute; top:0; height:100%;
    background:var(--its-color); opacity:0.7;
    min-width:2px; z-index:1;
  }}

  .chr-size {{
    font-family:'JetBrains Mono',monospace; font-size:10px; color:var(--muted);
    white-space:nowrap;
  }}

  .tooltip {{
    position:fixed; z-index:1000; pointer-events:none;
    background:var(--card); border:1px solid var(--border);
    border-radius:8px; padding:10px 14px; box-shadow:0 4px 16px var(--shadow);
    font-size:12px; color:var(--text2); max-width:260px;
    opacity:0; transition:opacity .15s;
    font-family:'Source Sans 3',sans-serif;
  }}
  .tooltip.visible {{ opacity:1; }}
  .tooltip .tt-title {{
    font-family:'JetBrains Mono',monospace; font-size:11px;
    font-weight:500; color:var(--text); margin-bottom:4px;
  }}
  .tooltip .tt-row {{ display:flex; justify-content:space-between; gap:12px; }}
  .tooltip .tt-val {{ font-family:'JetBrains Mono',monospace; font-size:11px; }}

  .footer {{
    text-align:center; margin-top:40px; padding-top:24px;
    border-top:1px solid var(--border);
    font-size:11px; color:var(--muted);
    opacity:0; animation:fadeUp .5s ease-out 1s forwards;
  }}

  @keyframes fadeUp {{
    from {{ opacity:0; transform:translateY(12px); }}
    to {{ opacity:1; transform:translateY(0); }}
  }}
</style>
</head>
<body>

<div class="toggle" onclick="toggleTheme()">
  <span class="toggle-icon" id="tIcon">&#9790;</span>
  <span id="tLabel">dark</span>
</div>

<div class="tooltip" id="tooltip"></div>

<div class="container">
  <div class="header">
    <div class="header-label">Chromosome Map</div>
    <h1>{assembly_name}</h1>
  </div>

  <div class="legend">
    <div class="legend-item"><div class="legend-dot" style="background:var(--micro-color)"></div>Microsatellites (&lt;10 bp)</div>
    <div class="legend-item"><div class="legend-dot" style="background:var(--complex-color)"></div>Complex repeats (&ge;10 bp)</div>
    <div class="legend-item"><div class="legend-dot" style="background:var(--telo-ok)"></div>Telomere present</div>
    <div class="legend-item"><div class="legend-dot" style="background:var(--telo-partial)"></div>Telomere partial</div>
    <div class="legend-item"><div class="legend-dot" style="background:var(--telo-absent)"></div>Telomere absent</div>
    <div class="legend-item"><div class="legend-dot" style="background:var(--its-color)"></div>ITS</div>
  </div>

  <div class="chr-list" id="chrList"></div>

  <div class="footer">Bin size: {bin_size_kb} kb &middot; satellome</div>
</div>

<script>
var DATA = {data_json};
var BIN_KB = {bin_size_kb};
var MAX_LEN = {max_len};

function render() {{
  var list = document.getElementById('chrList');
  list.innerHTML = '';
  DATA.forEach(function(c, idx) {{
    var row = document.createElement('div');
    row.className = 'chr-row';
    row.style.animationDelay = (0.15 + idx * 0.03) + 's';

    var sizeMb = (c.length / 1e6).toFixed(1);

    row.innerHTML =
      '<div class="chr-label">' + c.name + '</div>' +
      '<div class="chr-bar-wrap">' +
        '<div class="chr-bar" style="width:' + c.pct + '%">' +
          '<div class="chr-density" id="density-' + idx + '"></div>' +
          '<div class="telo-cap left ' + c.telo_left + '" title="Left: ' + c.telo_left + '"></div>' +
          '<div class="telo-cap right ' + c.telo_right + '" title="Right: ' + c.telo_right + '"></div>' +
        '</div>' +
      '</div>' +
      '<div class="chr-size">' + sizeMb + ' Mb</div>';

    list.appendChild(row);

    // Render density bins
    var densityEl = document.getElementById('density-' + idx);
    var binCount = c.density.length;
    if (binCount === 0) return;

    for (var i = 0; i < binCount; i++) {{
      var micro = c.density[i][0];
      var complex = c.density[i][1];
      if (micro + complex < 0.01) continue;

      var leftPct = (i / binCount * 100);
      var widthPct = (1 / binCount * 100);

      if (complex > 0.01) {{
        var bin = document.createElement('div');
        bin.className = 'density-bin';
        bin.style.left = leftPct + '%';
        bin.style.width = Math.max(widthPct, 0.3) + '%';
        bin.style.background = 'var(--complex-color)';
        bin.style.opacity = Math.min(0.2 + complex * 0.8, 1);
        densityEl.appendChild(bin);
      }}
      if (micro > 0.01) {{
        var bin2 = document.createElement('div');
        bin2.className = 'density-bin';
        bin2.style.left = leftPct + '%';
        bin2.style.width = Math.max(widthPct, 0.3) + '%';
        bin2.style.background = 'var(--micro-color)';
        bin2.style.opacity = Math.min(0.2 + micro * 0.8, 1);
        bin2.style.height = '50%';
        bin2.style.top = '50%';
        densityEl.appendChild(bin2);
      }}
    }}

    // Render ITS marks
    var bar = row.querySelector('.chr-bar');
    c.its.forEach(function(site) {{
      var mark = document.createElement('div');
      mark.className = 'its-mark';
      var pos = site.start / c.length * 100;
      var w = Math.max((site.end - site.start) / c.length * 100, 0.15);
      mark.style.left = pos + '%';
      mark.style.width = w + '%';
      mark.title = 'ITS: ' + site.count + ' repeats';
      bar.appendChild(mark);
    }});

    // Tooltip on hover
    bar.addEventListener('mousemove', function(e) {{
      var rect = bar.getBoundingClientRect();
      var x = (e.clientX - rect.left) / rect.width;
      var pos = Math.floor(x * c.length);
      var binIdx = Math.floor(x * binCount);
      if (binIdx >= binCount) binIdx = binCount - 1;
      var d = c.density[binIdx] || [0,0];

      var tt = document.getElementById('tooltip');
      tt.innerHTML =
        '<div class="tt-title">' + c.name + '</div>' +
        '<div class="tt-row"><span>Position</span><span class="tt-val">' + (pos/1e6).toFixed(2) + ' Mb</span></div>' +
        '<div class="tt-row"><span>Microsatellites</span><span class="tt-val">' + (d[0]*100).toFixed(1) + '%</span></div>' +
        '<div class="tt-row"><span>Complex</span><span class="tt-val">' + (d[1]*100).toFixed(1) + '%</span></div>' +
        '<div class="tt-row"><span>Left telomere</span><span class="tt-val">' + c.telo_left + '</span></div>' +
        '<div class="tt-row"><span>Right telomere</span><span class="tt-val">' + c.telo_right + '</span></div>';
      tt.style.left = (e.clientX + 14) + 'px';
      tt.style.top = (e.clientY - 10) + 'px';
      tt.classList.add('visible');
    }});
    bar.addEventListener('mouseleave', function() {{
      document.getElementById('tooltip').classList.remove('visible');
    }});
  }});
}}

function toggleTheme() {{
  var h = document.documentElement;
  var ic = document.getElementById('tIcon');
  var lb = document.getElementById('tLabel');
  if (h.getAttribute('data-theme') === 'dark') {{
    h.setAttribute('data-theme', 'light'); ic.innerHTML='\\u263E'; lb.textContent='dark';
  }} else {{
    h.setAttribute('data-theme', 'dark'); ic.innerHTML='\\u263C'; lb.textContent='light';
  }}
}}

render();
</script>
</body>
</html>'''
'''

Minimal usage:
    from satellome.core_functions.tools.chromosome_viz import create_chromosome_visualization
    create_chromosome_visualization(fai_path, sat_path, telomere_path, output_path)
'''


def create_chromosome_visualization(fai_path, sat_path, output_path,
                                     telomere_path=None, its_path=None,
                                     assembly_name=""):
    """
    Create interactive chromosome visualization HTML.

    Args:
        fai_path: Path to .fai index file
        sat_path: Path to .sat file (1kb recommended)
        output_path: Path to output HTML file
        telomere_path: Path to telomeres.tsv (optional)
        its_path: Path to telomeres.its.bed (optional)
        assembly_name: Assembly name for title
    """
    logger.info("Loading chromosome data...")
    chroms = load_fai(fai_path)
    if not chroms:
        logger.warning("No chromosomes found in .fai file")
        return

    logger.info(f"Loaded {len(chroms)} chromosomes")

    telomeres = {}
    if telomere_path and os.path.exists(telomere_path):
        telomeres = load_telomeres(telomere_path)
        logger.info(f"Loaded telomere data for {len(telomeres)} chromosomes")

    its = {}
    if its_path and os.path.exists(its_path):
        its = load_its(its_path)
        logger.info(f"Loaded {sum(len(v) for v in its.values())} ITS sites")

    logger.info("Binning repeat density...")
    bins = bin_repeats(sat_path, chroms)

    logger.info("Generating visualization...")
    generate_chromosome_html(chroms, bins, telomeres, its, assembly_name, output_path)
