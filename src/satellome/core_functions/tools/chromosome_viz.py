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

BIN_SIZE = 500_000  # 500 kb bins for density


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
    # Categories by period length
    CATS = ["micro", "mini", "satellite", "macro"]
    bins = {}
    for c in chroms:
        n_bins = (c["length"] // bin_size) + 1
        bins[c["name"]] = [{cat: 0 for cat in CATS} for _ in range(n_bins)]

    csv.field_size_limit(sys.maxsize)

    def _sat_lines(fh):
        """Yield lines, stripping '#' from header but skipping comment lines."""
        for line in fh:
            if line.startswith('#project') or line.startswith('#Project'):
                yield line[1:]  # strip '#' from header row
            elif line.startswith('#'):
                continue  # skip comments
            else:
                yield line

    with open(sat_path, 'r') as f:
        reader = csv.DictReader(_sat_lines(f), delimiter='\t')
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

            if period < 6:
                category = "micro"
            elif period < 10:
                category = "mini"
            elif period <= 100:
                category = "satellite"
            else:
                category = "macro"

            start_bin = start // bin_size
            end_bin = end // bin_size

            for b in range(start_bin, min(end_bin + 1, len(bins[chrom]))):
                b_start = b * bin_size
                b_end = (b + 1) * bin_size
                overlap = min(end, b_end) - max(start, b_start)
                if overlap > 0:
                    bins[chrom][b][category] += overlap

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

        # Compress bins to density values (0-1) per category
        chr_bins = bins.get(name, [])
        density = []
        for b in chr_bins:
            density.append([
                round(min(b.get("micro", 0) / BIN_SIZE, 1.0), 3),
                round(min(b.get("mini", 0) / BIN_SIZE, 1.0), 3),
                round(min(b.get("satellite", 0) / BIN_SIZE, 1.0), 3),
                round(min(b.get("macro", 0) / BIN_SIZE, 1.0), 3),
            ])

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

    # Count totals for summary
    total_bp = sum(c["length"] for c in chroms_sorted)
    t2t_count = sum(1 for c in chrom_data if c["t2t"])
    its_count = sum(len(c["its"]) for c in chrom_data)

    html = _build_html(data_json, bin_size_kb, assembly_name, max_len,
                       len(chroms_sorted), total_bp, t2t_count, its_count)

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        logger.info(f"Chromosome visualization: {output_path}")

    return html


def _build_html(data_json, bin_size_kb, assembly_name, max_len,
                n_chroms, total_bp, t2t_count, its_count):
    total_mb = round(total_bp / 1e6, 1)

    return f'''<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Chromosomes — {assembly_name}</title>
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;500&family=Playfair+Display:wght@700;900&family=Source+Sans+3:wght@300;400;600&display=swap" rel="stylesheet">
<style>
  :root, [data-theme="light"] {{
    --bg:#f8f6f1; --card:#fff; --border:#e2ddd5;
    --text:#1a1a1a; --text2:#555; --muted:#888;
    --shadow:rgba(0,0,0,.08);
    --chr-bg:#e8e5df; --chr-border:#d4d0c8;
    --c-micro:#0891b2; --c-mini:#0d9488; --c-sat:#7c3aed; --c-macro:#be123c;
    --c-telo-ok:#059669; --c-telo-part:#b45309; --c-telo-abs:#dc2626;
    --c-its:#db2777; --c-accent:#0891b2;
    --panel-bg:rgba(255,255,255,.92);
  }}
  [data-theme="dark"] {{
    --bg:#0a0e17; --card:#111827; --border:#1e293b;
    --text:#e2e8f0; --text2:#94a3b8; --muted:#64748b;
    --shadow:rgba(0,0,0,.3);
    --chr-bg:#1e293b; --chr-border:#334155;
    --c-micro:#22d3ee; --c-mini:#2dd4bf; --c-sat:#a78bfa; --c-macro:#fb7185;
    --c-telo-ok:#34d399; --c-telo-part:#fbbf24; --c-telo-abs:#f87171;
    --c-its:#f472b6; --c-accent:#22d3ee;
    --panel-bg:rgba(17,24,39,.92);
  }}
  *{{margin:0;padding:0;box-sizing:border-box;}}
  body{{background:var(--bg);font-family:'Source Sans 3',sans-serif;color:var(--text);}}
  .layout{{display:flex;min-height:100vh;}}

  /* === SIDE PANEL === */
  .panel{{
    width:220px; flex-shrink:0; position:sticky; top:0; height:100vh;
    background:var(--panel-bg); backdrop-filter:blur(12px); -webkit-backdrop-filter:blur(12px);
    border-right:1px solid var(--border); padding:24px 16px;
    display:flex; flex-direction:column; gap:20px; overflow-y:auto;
    z-index:10;
  }}
  .panel-title{{
    font-family:'Playfair Display',serif; font-size:18px; font-weight:700;
    color:var(--text); line-height:1.2;
  }}
  .panel-sub{{font-size:11px;color:var(--muted);margin-top:2px;font-family:'JetBrains Mono',monospace;}}
  .panel-section{{display:flex;flex-direction:column;gap:4px;}}
  .panel-section-title{{
    font-family:'JetBrains Mono',monospace; font-size:9px; font-weight:500;
    letter-spacing:1.5px; text-transform:uppercase; color:var(--muted); margin-bottom:4px;
  }}
  .panel-stat{{display:flex;justify-content:space-between;font-size:12px;color:var(--text2);}}
  .panel-stat .val{{font-family:'JetBrains Mono',monospace;font-size:11px;color:var(--text);}}

  /* Layer toggles */
  .layer-toggle{{
    display:flex;align-items:center;gap:8px;cursor:pointer;
    padding:5px 8px;border-radius:6px;transition:background .15s;
    font-size:12px;color:var(--text2);user-select:none;
  }}
  .layer-toggle:hover{{background:var(--border);}}
  .layer-toggle .dot{{width:10px;height:10px;border-radius:2px;flex-shrink:0;transition:opacity .2s;}}
  .layer-toggle.off .dot{{opacity:.2;}}
  .layer-toggle.off{{opacity:.5;}}
  .layer-toggle .lbl{{flex:1;}}

  .theme-btn{{
    margin-top:auto; padding:6px 10px; border-radius:6px;
    background:var(--card); border:1px solid var(--border); cursor:pointer;
    font-family:'JetBrains Mono',monospace; font-size:10px; color:var(--muted);
    text-align:center; transition:all .2s;
  }}
  .theme-btn:hover{{color:var(--c-accent);border-color:var(--c-accent);}}

  /* === MAIN CONTENT === */
  .main{{flex:1;padding:40px 32px;max-width:1100px;}}

  .header{{text-align:center;margin-bottom:40px;opacity:0;animation:fadeUp .6s ease-out forwards;}}
  .header-label{{
    font-family:'JetBrains Mono',monospace;font-size:10px;font-weight:500;
    letter-spacing:3px;text-transform:uppercase;color:var(--c-accent);margin-bottom:10px;
  }}
  .header h1{{
    font-family:'Playfair Display',serif;font-size:clamp(24px,3.5vw,36px);font-weight:900;
    background:linear-gradient(135deg,var(--text) 0%,var(--c-accent) 100%);
    -webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text;
  }}

  .chr-list{{display:flex;flex-direction:column;gap:5px;}}
  .chr-row{{
    display:grid;grid-template-columns:70px 1fr 55px;align-items:center;gap:10px;
    opacity:0;transform:translateY(10px);animation:fadeUp .35s ease-out forwards;padding:3px 0;
  }}
  .chr-label{{font-family:'JetBrains Mono',monospace;font-size:11px;font-weight:500;color:var(--text);text-align:right;white-space:nowrap;}}
  .chr-bar-wrap{{position:relative;height:20px;cursor:crosshair;}}
  .chr-bar{{
    position:absolute;top:0;left:0;height:100%;
    background:var(--chr-bg);border:1px solid var(--chr-border);
    border-radius:10px;overflow:hidden;transition:box-shadow .2s;
  }}
  .chr-bar:hover{{box-shadow:0 2px 12px var(--shadow);}}
  .chr-density{{position:absolute;top:0;left:0;width:100%;height:100%;}}
  .density-bin{{position:absolute;height:100%;}}

  .telo-cap{{
    position:absolute;top:-1px;width:7px;height:22px;border-radius:3px;z-index:2;
    transition:transform .15s;
  }}
  .telo-cap:hover{{transform:scaleY(1.25);}}
  .telo-cap.left{{left:-1px;border-top-right-radius:0;border-bottom-right-radius:0;}}
  .telo-cap.right{{right:-1px;border-top-left-radius:0;border-bottom-left-radius:0;}}
  .telo-cap.PRESENT{{background:var(--c-telo-ok);}}
  .telo-cap.PARTIAL{{background:var(--c-telo-part);}}
  .telo-cap.ABSENT{{background:var(--c-telo-abs);}}
  .telo-cap.UNKNOWN{{background:var(--muted);}}

  .its-mark{{position:absolute;top:0;height:100%;background:var(--c-its);opacity:.7;min-width:2px;z-index:1;}}

  .chr-size{{font-family:'JetBrains Mono',monospace;font-size:9px;color:var(--muted);white-space:nowrap;}}

  /* Tooltip */
  .tooltip{{
    position:fixed;z-index:1000;pointer-events:none;
    background:var(--card);border:1px solid var(--border);
    border-radius:8px;padding:10px 14px;box-shadow:0 4px 16px var(--shadow);
    font-size:12px;color:var(--text2);max-width:280px;
    opacity:0;transition:opacity .12s;
  }}
  .tooltip.visible{{opacity:1;}}
  .tooltip .tt-title{{font-family:'JetBrains Mono',monospace;font-size:11px;font-weight:500;color:var(--text);margin-bottom:6px;}}
  .tooltip .tt-row{{display:flex;justify-content:space-between;gap:12px;line-height:1.6;}}
  .tooltip .tt-val{{font-family:'JetBrains Mono',monospace;font-size:11px;}}

  .footer{{text-align:center;margin-top:36px;padding-top:20px;border-top:1px solid var(--border);font-size:10px;color:var(--muted);opacity:0;animation:fadeUp .5s ease-out 1s forwards;}}

  @keyframes fadeUp{{from{{opacity:0;transform:translateY(10px);}}to{{opacity:1;transform:translateY(0);}}}}
  @media(max-width:768px){{.panel{{display:none;}}.main{{padding:20px 12px;}}}}
</style>
</head>
<body>

<div class="tooltip" id="tooltip"></div>

<div class="layout">
  <!-- SIDE PANEL -->
  <div class="panel">
    <div>
      <div class="panel-title">{assembly_name}</div>
      <div class="panel-sub">Chromosome Map</div>
    </div>

    <div class="panel-section">
      <div class="panel-section-title">Assembly</div>
      <div class="panel-stat"><span>Chromosomes</span><span class="val">{n_chroms}</span></div>
      <div class="panel-stat"><span>Genome size</span><span class="val">{total_mb} Mb</span></div>
      <div class="panel-stat"><span>T2T chromosomes</span><span class="val">{t2t_count}/{n_chroms}</span></div>
      <div class="panel-stat"><span>ITS sites</span><span class="val">{its_count}</span></div>
    </div>

    <div class="panel-section">
      <div class="panel-section-title">Layers</div>
      <div class="layer-toggle" data-layer="micro" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-micro)"></div><span class="lbl">Micro (&lt;6 bp)</span>
      </div>
      <div class="layer-toggle" data-layer="mini" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-mini)"></div><span class="lbl">Mini (6-9 bp)</span>
      </div>
      <div class="layer-toggle" data-layer="satellite" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-sat)"></div><span class="lbl">Satellite (10-100 bp)</span>
      </div>
      <div class="layer-toggle" data-layer="macro" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-macro)"></div><span class="lbl">Macro (&gt;100 bp)</span>
      </div>
      <div class="layer-toggle" data-layer="its" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-its)"></div><span class="lbl">ITS</span>
      </div>
      <div class="layer-toggle" data-layer="telomere" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-telo-ok)"></div><span class="lbl">Telomeres</span>
      </div>
    </div>

    <div class="theme-btn" onclick="toggleTheme()" id="themeBtn">&#9790; Dark mode</div>
  </div>

  <!-- MAIN -->
  <div class="main">
    <div class="header">
      <div class="header-label">Tandem Repeat Landscape</div>
      <h1>{assembly_name}</h1>
    </div>
    <div class="chr-list" id="chrList"></div>
    <div class="footer">Bin size: {bin_size_kb} kb &middot; satellome</div>
  </div>
</div>

<script>
var DATA={data_json};
var BIN_KB={bin_size_kb};
var LAYERS={{micro:true,mini:true,satellite:true,macro:true,its:true,telomere:true}};
var COLORS=['--c-micro','--c-mini','--c-sat','--c-macro'];
var CAT_NAMES=['Micro','Mini','Satellite','Macro'];

function toggleLayer(el){{
  var layer=el.getAttribute('data-layer');
  LAYERS[layer]=!LAYERS[layer];
  el.classList.toggle('off',!LAYERS[layer]);
  render();
}}

function render(){{
  var list=document.getElementById('chrList');
  list.innerHTML='';
  DATA.forEach(function(c,idx){{
    var row=document.createElement('div');
    row.className='chr-row';
    row.style.animationDelay=(0.08+idx*0.025)+'s';
    var sizeMb=(c.length/1e6).toFixed(1);

    row.innerHTML=
      '<div class="chr-label">'+c.name+'</div>'+
      '<div class="chr-bar-wrap">'+
        '<div class="chr-bar" style="width:'+c.pct+'%">'+
          '<div class="chr-density" id="d-'+idx+'"></div>'+
          (LAYERS.telomere?'<div class="telo-cap left '+c.telo_left+'"></div><div class="telo-cap right '+c.telo_right+'"></div>':'')+
        '</div>'+
      '</div>'+
      '<div class="chr-size">'+sizeMb+' Mb</div>';

    list.appendChild(row);

    var de=document.getElementById('d-'+idx);
    var bc=c.density.length;
    if(bc===0)return;

    for(var i=0;i<bc;i++){{
      var lp=i/bc*100;
      var wp=1/bc*100;

      // Stack categories: each gets a vertical slice
      for(var cat=0;cat<4;cat++){{
        var key=['micro','mini','satellite','macro'][cat];
        if(!LAYERS[key])continue;
        var val=c.density[i][cat];
        if(val<0.005)continue;
        var bin=document.createElement('div');
        bin.className='density-bin layer-'+key;
        bin.style.left=lp+'%';
        bin.style.width=Math.max(wp,0.5)+'%';
        bin.style.background='var('+COLORS[cat]+')';
        bin.style.opacity=Math.min(0.3+val*0.7,1);
        // Stack vertically: each category gets a quarter
        bin.style.top=(cat*25)+'%';
        bin.style.height='25%';
        de.appendChild(bin);
      }}
    }}

    // ITS
    if(LAYERS.its){{
      var bar=row.querySelector('.chr-bar');
      c.its.forEach(function(s){{
        var m=document.createElement('div');
        m.className='its-mark';
        m.style.left=(s.start/c.length*100)+'%';
        m.style.width=Math.max((s.end-s.start)/c.length*100,0.15)+'%';
        bar.appendChild(m);
      }});
    }}

    // Tooltip
    var bar=row.querySelector('.chr-bar');
    bar.addEventListener('mousemove',function(e){{
      var rect=bar.getBoundingClientRect();
      var x=(e.clientX-rect.left)/rect.width;
      var pos=Math.floor(x*c.length);
      var bi=Math.min(Math.floor(x*bc),bc-1);
      var d=c.density[bi]||[0,0,0,0];
      var tt=document.getElementById('tooltip');
      tt.innerHTML=
        '<div class="tt-title">'+c.name+' : '+(pos/1e6).toFixed(2)+' Mb</div>'+
        '<div class="tt-row"><span>Micro (&lt;6 bp)</span><span class="tt-val">'+(d[0]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>Mini (6-9 bp)</span><span class="tt-val">'+(d[1]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>Satellite (10-100)</span><span class="tt-val">'+(d[2]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>Macro (&gt;100 bp)</span><span class="tt-val">'+(d[3]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>Telomere L / R</span><span class="tt-val">'+c.telo_left+' / '+c.telo_right+'</span></div>';
      tt.style.left=(e.clientX+14)+'px';
      tt.style.top=(e.clientY-10)+'px';
      tt.classList.add('visible');
    }});
    bar.addEventListener('mouseleave',function(){{
      document.getElementById('tooltip').classList.remove('visible');
    }});
  }});
}}

function toggleTheme(){{
  var h=document.documentElement,b=document.getElementById('themeBtn');
  if(h.getAttribute('data-theme')==='dark'){{
    h.setAttribute('data-theme','light');b.innerHTML='\\u263E Dark mode';
  }}else{{
    h.setAttribute('data-theme','dark');b.innerHTML='\\u263C Light mode';
  }}
}}

render();
</script>
</body>
</html>'''


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
