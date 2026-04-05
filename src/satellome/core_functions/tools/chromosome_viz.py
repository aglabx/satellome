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

BIN_SIZE = 500_000       # 500 kb bins for genome/chromosome view
FINE_BIN_SIZE = 50_000   # 50 kb bins for region view (x10 zoom)


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


def load_families(families_path):
    """Load sat_families.tsv and return list of family records + top family names."""
    families = []  # [{chr, start, end, length, family, consensus}, ...]
    family_counts = {}
    if not os.path.exists(families_path):
        return families, []

    with open(families_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('chr\t'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            fam = parts[5]
            families.append({
                "chr": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "length": int(parts[3]),
                "family": fam,
                "consensus": parts[6] if len(parts) > 6 else "",
            })
            family_counts[fam] = family_counts.get(fam, 0) + 1

    # Top families by count (for color assignment)
    top = sorted(family_counts.keys(), key=lambda k: -family_counts[k])
    return families, top


def bin_families(families, chroms, top_families, bin_size=BIN_SIZE):
    """Bin family assignments into density tracks per chromosome.

    Returns dict: {chr_name: [{family_name: bp, ...}, ...]}
    """
    # Use top N families, rest goes to "other"
    MAX_FAMILIES = 12
    tracked = set(top_families[:MAX_FAMILIES])

    bins = {}
    for c in chroms:
        n_bins = (c["length"] // bin_size) + 1
        bins[c["name"]] = [{} for _ in range(n_bins)]

    for f in families:
        chrom = f["chr"]
        if chrom not in bins:
            continue
        start, end = f["start"], f["end"]
        fam = f["family"] if f["family"] in tracked else "other"

        start_bin = start // bin_size
        end_bin = end // bin_size
        for b in range(start_bin, min(end_bin + 1, len(bins[chrom]))):
            b_start = b * bin_size
            b_end = (b + 1) * bin_size
            overlap = min(end, b_end) - max(start, b_start)
            if overlap > 0:
                bins[chrom][b][fam] = bins[chrom][b].get(fam, 0) + overlap

    return bins, list(tracked)


def _sat_lines(fh):
    """Yield lines from SAT file, fixing '#project' header."""
    for line in fh:
        if line.startswith('#project') or line.startswith('#Project'):
            yield line[1:]
        elif line.startswith('#'):
            continue
        else:
            yield line


# Categories by array length (matching satellome classification approach)
CATEGORIES = ["lt1kb", "1-10kb", "10-100kb", "gt100kb"]
CAT_DISPLAY = ["< 1 kb", "1-10 kb", "10-100 kb", "> 100 kb"]


def _classify_by_array_length(array_length):
    """Classify repeat by array length."""
    if array_length < 1000:
        return "lt1kb"
    elif array_length < 10000:
        return "1-10kb"
    elif array_length < 100000:
        return "10-100kb"
    else:
        return "gt100kb"


def _add_to_bins(bins_list, start, end, category, bin_size):
    """Add overlap to bins."""
    start_bin = start // bin_size
    end_bin = end // bin_size
    for b in range(start_bin, min(end_bin + 1, len(bins_list))):
        b_start = b * bin_size
        b_end = (b + 1) * bin_size
        overlap = min(end, b_end) - max(start, b_start)
        if overlap > 0:
            bins_list[b][category] += overlap


def bin_repeats(sat_path, chroms):
    """Bin tandem repeats at two resolutions: coarse (500kb) and fine (50kb).

    Returns:
        (coarse_bins, fine_bins) — both dicts {chr_name: [{cat: bp}, ...]}
    """
    coarse = {}
    fine = {}
    chroms_set = set()
    for c in chroms:
        n_coarse = (c["length"] // BIN_SIZE) + 1
        n_fine = (c["length"] // FINE_BIN_SIZE) + 1
        coarse[c["name"]] = [{cat: 0 for cat in CATEGORIES} for _ in range(n_coarse)]
        fine[c["name"]] = [{cat: 0 for cat in CATEGORIES} for _ in range(n_fine)]
        chroms_set.add(c["name"])

    csv.field_size_limit(sys.maxsize)

    with open(sat_path, 'r') as f:
        reader = csv.DictReader(_sat_lines(f), delimiter='\t')
        for row in reader:
            chrom = row.get("trf_head", "")
            if chrom not in chroms_set:
                continue
            try:
                start = int(row.get("trf_l_ind", 0))
                end = int(row.get("trf_r_ind", 0))
            except (ValueError, TypeError):
                continue

            category = _classify_by_array_length(end - start)
            _add_to_bins(coarse[chrom], start, end, category, BIN_SIZE)
            _add_to_bins(fine[chrom], start, end, category, FINE_BIN_SIZE)

    return coarse, fine


def _load_monomers_compact(monomers_tsv, min_array_length=10000):
    """Load monomer lengths for large arrays. Returns dict for JSON embedding.

    Returns: {array_id: {c: chr, s: start, e: end, p: period, m: [lengths]}}
    """
    csv.field_size_limit(sys.maxsize)
    arrays = {}

    logger.info(f"Loading monomers from {monomers_tsv}...")
    with open(monomers_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            aid = row['array_id']
            parts = aid.split('_')
            if len(parts) < 5:
                continue
            arr_len = int(parts[3])
            if arr_len < min_array_length:
                continue

            if aid not in arrays:
                arrays[aid] = {
                    "c": parts[0],
                    "s": int(parts[1]),
                    "e": int(parts[2]),
                    "p": int(parts[4]),
                    "m": [],
                }

            if row['type'] == 'monomer':
                arrays[aid]["m"].append(int(row.get('length', 0)))

    logger.info(f"Loaded {len(arrays)} arrays with monomers")
    return arrays


def generate_chromosome_html(chroms, bins, fine_bins, telomeres, its, assembly_name="",
                             output_path=None, monomers_data=None,
                             family_bins=None, tracked_families=None):
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
                round(min(b.get("lt1kb", 0) / BIN_SIZE, 1.0), 3),
                round(min(b.get("1-10kb", 0) / BIN_SIZE, 1.0), 3),
                round(min(b.get("10-100kb", 0) / BIN_SIZE, 1.0), 3),
                round(min(b.get("gt100kb", 0) / BIN_SIZE, 1.0), 3),
            ])

        # Fine density for region zoom
        fine_chr_bins = fine_bins.get(name, [])
        fine_density = []
        for b in fine_chr_bins:
            fine_density.append([
                round(min(b.get("lt1kb", 0) / FINE_BIN_SIZE, 1.0), 3),
                round(min(b.get("1-10kb", 0) / FINE_BIN_SIZE, 1.0), 3),
                round(min(b.get("10-100kb", 0) / FINE_BIN_SIZE, 1.0), 3),
                round(min(b.get("gt100kb", 0) / FINE_BIN_SIZE, 1.0), 3),
            ])

        # Family density
        fam_density = []
        if family_bins and name in family_bins:
            fam_chr = family_bins[name]
            tf = tracked_families or []
            for b in fam_chr:
                fam_bin = {}
                for fam_name in tf:
                    val = b.get(fam_name, 0)
                    if val > 0:
                        fam_bin[fam_name] = round(min(val / BIN_SIZE, 1.0), 3)
                other = b.get("other", 0)
                if other > 0:
                    fam_bin["other"] = round(min(other / BIN_SIZE, 1.0), 3)
                fam_density.append(fam_bin)

        chrom_data.append({
            "name": name,
            "length": c["length"],
            "pct": round(c["length"] / max_len * 100, 2),
            "telo_left": telo.get("left_status", "UNKNOWN"),
            "telo_right": telo.get("right_status", "UNKNOWN"),
            "t2t": telo.get("t2t", False),
            "its": chrom_its,
            "density": density,
            "fine": fine_density,
            "fam": fam_density,
        })

    data_json = json.dumps(chrom_data, separators=(',', ':'))
    bin_size_kb = BIN_SIZE // 1000

    # Count totals for summary
    total_bp = sum(c["length"] for c in chroms_sorted)
    t2t_count = sum(1 for c in chrom_data if c["t2t"])
    its_count = sum(len(c["its"]) for c in chrom_data)

    fine_bin_kb = FINE_BIN_SIZE // 1000
    monomers_json = json.dumps(monomers_data or {}, separators=(',', ':'))
    families_list_json = json.dumps(tracked_families or [], separators=(',', ':'))
    html = _build_html(data_json, monomers_json, families_list_json,
                       bin_size_kb, fine_bin_kb, assembly_name, max_len,
                       len(chroms_sorted), total_bp, t2t_count, its_count)

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        logger.info(f"Chromosome visualization: {output_path}")

    return html


def _build_html(data_json, monomers_json, families_list_json,
                bin_size_kb, fine_bin_kb, assembly_name, max_len,
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

  /* Nav levels */
  .nav-level{{
    display:flex;align-items:center;gap:8px;cursor:pointer;
    padding:5px 8px;border-radius:6px;transition:all .15s;
    font-size:12px;color:var(--text2);user-select:none;
  }}
  .nav-level:hover:not(.disabled){{background:var(--border);}}
  .nav-level.active{{color:var(--c-accent);font-weight:500;}}
  .nav-level.active .nav-icon{{color:var(--c-accent);}}
  .nav-level.disabled{{opacity:.35;cursor:default;}}
  .nav-icon{{font-size:10px;width:12px;text-align:center;color:var(--muted);}}

  /* View toggles */
  .view-toggle{{
    display:flex;align-items:center;gap:6px;cursor:pointer;
    padding:5px 8px;border-radius:6px;transition:all .15s;
    font-size:12px;color:var(--text2);user-select:none;
  }}
  .view-toggle:hover{{background:var(--border);}}
  .view-toggle.active{{color:var(--c-accent);font-weight:500;
    background:rgba(8,145,178,.08);
  }}
  [data-theme="dark"] .view-toggle.active{{background:rgba(34,211,238,.08);}}
  .view-toggle .soon{{
    font-family:'JetBrains Mono',monospace;font-size:8px;
    letter-spacing:.5px;text-transform:uppercase;
    color:var(--muted);background:var(--border);
    padding:1px 5px;border-radius:3px;margin-left:auto;
  }}

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
    border-radius:10px;transition:box-shadow .2s;
  }}
  .chr-bar:hover{{box-shadow:0 2px 12px var(--shadow);}}
  .chr-density{{position:absolute;top:0;left:0;width:100%;height:100%;overflow:hidden;border-radius:10px;}}
  .density-bin{{position:absolute;height:100%;}}

  .telo-cap{{
    position:absolute;top:0;width:7px;height:100%;border-radius:3px;z-index:2;
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

  /* === CHROMOSOME PICKER === */
  .chr-picker{{
    display:flex;flex-wrap:wrap;gap:4px;margin-bottom:16px;
    justify-content:center;
  }}
  .chr-pick{{
    display:flex;align-items:center;gap:3px;
    padding:4px 8px;border-radius:5px;cursor:pointer;
    font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--text2);
    border:1px solid transparent;transition:all .15s;user-select:none;
  }}
  .chr-pick:hover{{background:var(--border);}}
  .chr-pick.active{{
    background:var(--card);border-color:var(--c-accent);color:var(--text);
    font-weight:500;box-shadow:0 1px 4px var(--shadow);
  }}
  .chr-pick-dot{{width:5px;height:5px;border-radius:50%;flex-shrink:0;}}
  .chr-pick-name{{}}

  /* === CHROMOSOME DETAIL VIEW === */
  .chr-detail{{display:none;}}
  .chr-detail.active{{display:block;}}

  .chr-detail-header{{
    display:flex;align-items:center;gap:12px;margin-bottom:20px;
  }}
  .back-btn{{
    padding:6px 12px;border-radius:6px;cursor:pointer;
    background:var(--card);border:1px solid var(--border);
    font-family:'JetBrains Mono',monospace;font-size:11px;color:var(--text2);
    transition:all .15s;user-select:none;
  }}
  .back-btn:hover{{color:var(--c-accent);border-color:var(--c-accent);}}

  .chr-detail-title{{
    font-family:'Playfair Display',serif;font-size:24px;font-weight:700;color:var(--text);
  }}
  .chr-detail-sub{{font-family:'JetBrains Mono',monospace;font-size:11px;color:var(--muted);margin-left:8px;}}

  .chr-grid{{
    display:flex;flex-direction:column;gap:0px;
    background:var(--card);border:1px solid var(--border);border-radius:8px;
    padding:6px;overflow:hidden;
  }}
  .chr-grid-row{{
    display:flex;position:relative;border-radius:1px;overflow:hidden;
    background:var(--chr-bg);cursor:crosshair;
  }}
  .chr-grid-cell{{position:absolute;top:0;height:100%;}}

  .chr-grid-labels{{
    display:flex;justify-content:space-between;
    font-family:'JetBrains Mono',monospace;font-size:9px;color:var(--muted);
    margin-top:6px;padding:0 8px;
  }}

  /* === REPEAT VIEW === */
  .repeat-card{{
    background:var(--card);border:1px solid var(--border);border-radius:8px;
    margin-bottom:10px;overflow:hidden;
  }}
  .repeat-header{{
    padding:10px 14px;display:flex;justify-content:space-between;align-items:center;
    border-bottom:1px solid var(--border);
  }}
  .repeat-title{{
    font-family:'JetBrains Mono',monospace;font-size:11px;font-weight:500;color:var(--text);
  }}
  .repeat-info{{
    font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--muted);
  }}
  .repeat-monomers{{
    display:flex;height:24px;padding:4px;gap:1px;align-items:stretch;
  }}
  .mono-block{{
    height:100%;border-radius:1px;min-width:1px;
    transition:opacity .1s;cursor:default;
  }}
  .mono-block:hover{{opacity:0.7;}}

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
      <div class="layer-toggle" data-layer="lt1kb" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-micro)"></div><span class="lbl">&lt; 1 kb</span>
      </div>
      <div class="layer-toggle" data-layer="1-10kb" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-mini)"></div><span class="lbl">1 &ndash; 10 kb</span>
      </div>
      <div class="layer-toggle" data-layer="10-100kb" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-sat)"></div><span class="lbl">10 &ndash; 100 kb</span>
      </div>
      <div class="layer-toggle" data-layer="gt100kb" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-macro)"></div><span class="lbl">&gt; 100 kb</span>
      </div>
      <div class="layer-toggle" data-layer="its" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-its)"></div><span class="lbl">ITS</span>
      </div>
      <div class="layer-toggle" data-layer="telomere" onclick="toggleLayer(this)">
        <div class="dot" style="background:var(--c-telo-ok)"></div><span class="lbl">Telomeres</span>
      </div>
    </div>

    <div class="panel-section">
      <div class="panel-section-title">Navigation</div>
      <div class="nav-level active" data-level="genome" onclick="setLevel(this)">
        <span class="nav-icon">&#9673;</span><span class="lbl">Genome</span>
      </div>
      <div class="nav-level" data-level="chromosome" onclick="setLevel(this)">
        <span class="nav-icon">&#9656;</span><span class="lbl">Chromosome</span>
      </div>
      <div class="nav-level" data-level="region" onclick="setLevel(this)">
        <span class="nav-icon">&#9656;</span><span class="lbl">Region</span>
      </div>
      <div class="nav-level" data-level="repeat" onclick="setLevel(this)">
        <span class="nav-icon">&#9656;</span><span class="lbl">Repeat</span>
      </div>
      <div class="nav-level disabled" data-level="monomer">
        <span class="nav-icon">&#9656;</span><span class="lbl">Monomer</span>
      </div>
    </div>

    <div class="panel-section">
      <div class="panel-section-title">View</div>
      <div class="view-toggle active" data-view="size" onclick="setView(this)">
        <span class="lbl">By array size</span>
      </div>
      <div class="view-toggle" data-view="family" onclick="setView(this)">
        <span class="lbl">By superfamily</span>
      </div>
      <div class="view-toggle" data-view="class" onclick="setView(this)">
        <span class="lbl">By family</span><span class="soon">soon</span>
      </div>
      <div class="view-toggle" data-view="monomer" onclick="setView(this)">
        <span class="lbl">By monomer</span><span class="soon">soon</span>
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
    <div class="chr-detail" id="chrDetail"></div>
    <div class="footer">Bin size: {bin_size_kb} kb &middot; satellome</div>
  </div>
</div>

<script>
var DATA={data_json};
var BIN_KB={bin_size_kb};
var FAMILIES={families_list_json};
var CURRENT_VIEW='size';
var LAYERS={{'lt1kb':true,'1-10kb':true,'10-100kb':true,'gt100kb':true,its:true,telomere:true}};
var COLORS=['--c-micro','--c-mini','--c-sat','--c-macro'];
var CAT_KEYS=['lt1kb','1-10kb','10-100kb','gt100kb'];
var CAT_NAMES=['< 1 kb','1-10 kb','10-100 kb','> 100 kb'];

function toggleLayer(el){{
  var layer=el.getAttribute('data-layer');
  LAYERS[layer]=!LAYERS[layer];
  el.classList.toggle('off',!LAYERS[layer]);
  if(CURRENT_CHR!==null){{
    showChromosome(CURRENT_CHR);
  }}else{{
    render();
  }}
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
      '<div class="chr-label" style="cursor:pointer;text-decoration:underline;text-decoration-color:var(--border);text-underline-offset:3px" onclick="showChromosome('+idx+')">'+c.name+'</div>'+
      '<div class="chr-bar-wrap">'+
        '<div class="chr-bar" style="width:'+c.pct+'%">'+
          '<div class="chr-density" id="d-'+idx+'"></div>'+
          (LAYERS.telomere&&c.pct>1?'<div class="telo-cap left '+c.telo_left+'"></div><div class="telo-cap right '+c.telo_right+'"></div>':'')+
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

      if(CURRENT_VIEW==='family' && c.fam && c.fam[i]){{
        // Family view: color by family
        var famBin=c.fam[i];
        var famKeys=Object.keys(famBin);
        var sliceH=famKeys.length>0?100/famKeys.length:100;
        famKeys.forEach(function(fk,fi){{
          var val=famBin[fk];
          if(val<0.005)return;
          var famIdx=FAMILIES.indexOf(fk);
          var color=famIdx>=0?FAM_COLORS[famIdx%FAM_COLORS.length]:'#888';
          var bin=document.createElement('div');
          bin.className='density-bin';
          bin.style.left=lp+'%';
          bin.style.width=Math.max(wp,0.5)+'%';
          bin.style.background=color;
          bin.style.opacity=Math.min(0.3+val*0.7,1);
          bin.style.top=(fi*sliceH)+'%';
          bin.style.height=sliceH+'%';
          de.appendChild(bin);
        }});
      }} else if(CURRENT_VIEW==='size'){{
        // Size view: color by array length category
        for(var cat=0;cat<4;cat++){{
          var key=CAT_KEYS[cat];
          if(!LAYERS[key])continue;
          var val=c.density[i][cat];
          if(val<0.005)continue;
          var bin=document.createElement('div');
          bin.className='density-bin layer-'+key;
          bin.style.left=lp+'%';
          bin.style.width=Math.max(wp,0.5)+'%';
          bin.style.background='var('+COLORS[cat]+')';
          bin.style.opacity=Math.min(0.3+val*0.7,1);
          bin.style.top=(cat*25)+'%';
          bin.style.height='25%';
          de.appendChild(bin);
        }}
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
        '<div class="tt-row"><span>&lt; 1 kb</span><span class="tt-val">'+(d[0]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>1-10 kb</span><span class="tt-val">'+(d[1]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>10-100 kb</span><span class="tt-val">'+(d[2]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>&gt; 100 kb</span><span class="tt-val">'+(d[3]*100).toFixed(1)+'%</span></div>'+
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

var CURRENT_CHR=null;

function showChromosome(idx){{
  CURRENT_CHR=idx;
  var c=DATA[idx];
  document.getElementById('chrList').style.display='none';
  document.querySelector('.header').style.display='none';
  document.querySelector('.footer').style.display='none';
  var detail=document.getElementById('chrDetail');
  detail.classList.add('active');

  var bc=c.density.length;
  var bpPerBin=BIN_KB*1000;

  // Calculate grid to fill viewport as a square
  // Available height: viewport minus header/padding (~100px for header bar)
  var availH=window.innerHeight-140;
  var availW=detail.offsetWidth||800;

  // We want a roughly square grid: cols/rows ≈ availW/availH
  // Each cell: cellW = availW/cols, cellH = availH/rows
  // We want cellW ≈ cellH (square cells)
  // cols * rows = bc, cols/rows = availW/availH
  // cols = sqrt(bc * availW/availH)
  var aspect=availW/availH;
  var cols=Math.round(Math.sqrt(bc*aspect));
  if(cols<5)cols=5;
  if(cols>bc)cols=bc;
  var numRows=Math.ceil(bc/cols);
  var bpPerRow=cols*bpPerBin;

  // Cell height to fill available vertical space
  var rowH=Math.max(Math.floor(availH/numRows),2);

  // Telomere status colors
  function teloColor(status){{
    if(status==='PRESENT')return 'var(--c-telo-ok)';
    if(status==='PARTIAL')return 'var(--c-telo-part)';
    return 'var(--c-telo-abs)';
  }}

  // Chromosome picker bar
  var picker='<div class="chr-picker">';
  DATA.forEach(function(ch,i){{
    var cls=(i===idx)?'chr-pick active':'chr-pick';
    var tc=ch.t2t?'var(--c-telo-ok)':'var(--muted)';
    picker+='<div class="'+cls+'" onclick="showChromosome('+i+')" title="'+ch.name+' ('+( ch.length/1e6).toFixed(0)+' Mb)">'+
      '<span class="chr-pick-dot" style="background:'+tc+'"></span>'+
      '<span class="chr-pick-name">'+ch.name.replace("chr","")+'</span></div>';
  }});
  picker+='</div>';

  var html=picker+
    '<div class="chr-detail-header">'+
    '<div class="back-btn" onclick="showGenome()">&#8592; Genome</div>'+
    '<div class="chr-detail-title">'+c.name+'</div>'+
    '<div class="chr-detail-sub">'+( c.length/1e6).toFixed(1)+' Mb'+
    ' &middot; <span style="color:'+teloColor(c.telo_left)+'">&#9646;</span> '+c.telo_left+
    ' / '+c.telo_right+' <span style="color:'+teloColor(c.telo_right)+'">&#9646;</span>'+
    (c.t2t?' &middot; <span style="color:var(--c-telo-ok)">T2T</span>':'')+
    ' &middot; '+c.its.length+' ITS'+
    '</div></div>';

  html+='<div class="chr-grid">';

  // First row: telomere indicator at start
  for(var row=0;row<numRows;row++){{
    var startBin=row*cols;
    var endBin=Math.min(startBin+cols,bc);
    var isFirst=(row===0);
    var isLast=(row===numRows-1);

    // Last row: clip width to actual bins
    var rowWidthPct=(endBin-startBin)/cols*100;
    html+='<div class="chr-grid-row" data-row="'+row+'" style="height:'+rowH+'px;width:'+rowWidthPct+'%">';

    // Telomere caps on first/last row
    if(isFirst && LAYERS.telomere){{
      html+='<div class="chr-grid-cell" style="left:0;width:2px;background:'+teloColor(c.telo_left)+';z-index:2;border-radius:1px 0 0 1px"></div>';
    }}
    if(isLast && LAYERS.telomere){{
      html+='<div class="chr-grid-cell" style="right:0;left:auto;width:2px;background:'+teloColor(c.telo_right)+';z-index:2;border-radius:0 1px 1px 0"></div>';
    }}

    var actualCols=endBin-startBin;
    for(var i=startBin;i<endBin;i++){{
      var d=c.density[i];
      var lp=((i-startBin)/actualCols*100);
      var wp=(1/actualCols*100);

      // Render all active categories stacked
      for(var cat=0;cat<4;cat++){{
        var key=CAT_KEYS[cat];
        if(!LAYERS[key])continue;
        var val=d[cat];
        if(val<0.005)continue;
        html+='<div class="chr-grid-cell" style="'+
          'left:'+lp+'%;width:'+Math.max(wp,0.2)+'%;'+
          'top:'+(cat*25)+'%;height:25%;'+
          'background:var('+COLORS[cat]+');'+
          'opacity:'+Math.min(0.25+val*0.75,1)+
          '"></div>';
      }}
    }}

    // ITS marks
    if(LAYERS.its){{
      c.its.forEach(function(s){{
        var sBin=Math.floor(s.start/bpPerBin);
        var eBin=Math.ceil(s.end/bpPerBin);
        if(sBin>=startBin && sBin<endBin){{
          var lp=((sBin-startBin)/actualCols*100);
          var wp=Math.max(((eBin-sBin)/actualCols*100),0.3);
          html+='<div class="chr-grid-cell" style="left:'+lp+'%;width:'+wp+'%;top:0;height:100%;background:var(--c-its);opacity:0.8;z-index:1"></div>';
        }}
      }});
    }}

    html+='</div>';
  }}
  html+='</div>';

  // Row labels on left side
  html+='<div class="chr-grid-labels">'+
    '<span>0 Mb</span>'+
    '<span>'+(bpPerRow/1e6).toFixed(1)+' Mb / row &middot; '+cols+' bins &times; '+numRows+' rows</span>'+
    '<span>'+(c.length/1e6).toFixed(1)+' Mb</span>'+
    '</div>';

  detail.innerHTML=html;

  // Add tooltip to grid rows
  detail.querySelectorAll('.chr-grid-row').forEach(function(rowEl){{
    rowEl.addEventListener('mousemove',function(e){{
      var rect=rowEl.getBoundingClientRect();
      var x=(e.clientX-rect.left)/rect.width;
      var rowIdx=parseInt(rowEl.getAttribute('data-row'));
      var rStart=rowIdx*cols;
      var rEnd=Math.min(rStart+cols,bc);
      var rActual=rEnd-rStart;
      var binIdx=rStart+Math.floor(x*rActual);
      if(binIdx>=bc)binIdx=bc-1;
      var pos=binIdx*bpPerBin;
      var d=c.density[binIdx]||[0,0,0,0];
      var tt=document.getElementById('tooltip');
      tt.innerHTML=
        '<div class="tt-title">'+c.name+' : '+(pos/1e6).toFixed(2)+' Mb</div>'+
        '<div class="tt-row"><span>&lt; 1 kb</span><span class="tt-val">'+(d[0]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>1-10 kb</span><span class="tt-val">'+(d[1]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>10-100 kb</span><span class="tt-val">'+(d[2]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>&gt; 100 kb</span><span class="tt-val">'+(d[3]*100).toFixed(1)+'%</span></div>';
      tt.style.left=(e.clientX+14)+'px';
      tt.style.top=(e.clientY-10)+'px';
      tt.classList.add('visible');
    }});
    rowEl.addEventListener('mouseleave',function(){{
      document.getElementById('tooltip').classList.remove('visible');
    }});
    // Click to zoom into region
    rowEl.addEventListener('click',function(e){{
      var rect=rowEl.getBoundingClientRect();
      var x=(e.clientX-rect.left)/rect.width;
      var rowIdx=parseInt(rowEl.getAttribute('data-row'));
      var rowStart=rowIdx*cols;
      var rowEnd=Math.min(rowStart+cols,bc);
      var rowActual=rowEnd-rowStart;
      var binIdx=rowStart+Math.floor(x*rowActual);
      if(binIdx>=bc)binIdx=bc-1;
      showRegion(idx,binIdx);
    }});
  }});

  // Update nav
  setNavLevel('chromosome');
}}

var CURRENT_REGION=null;
var FINE_BIN_KB={fine_bin_kb};

function showRegion(chrIdx, coarseBinIdx){{
  CURRENT_REGION=coarseBinIdx;
  var c=DATA[chrIdx];
  var bpPerCoarse=BIN_KB*1000;
  var bpPerFine=FINE_BIN_KB*1000;
  var regionStart=coarseBinIdx*bpPerCoarse;
  var regionEnd=Math.min(regionStart+bpPerCoarse,c.length);

  // Fine bins for this region
  var fineStart=Math.floor(regionStart/bpPerFine);
  var fineEnd=Math.ceil(regionEnd/bpPerFine);
  var fineBins=c.fine.slice(fineStart,fineEnd);
  var fb=fineBins.length;
  if(fb===0)return;

  document.getElementById('chrList').style.display='none';
  document.querySelector('.header').style.display='none';
  document.querySelector('.footer').style.display='none';
  var detail=document.getElementById('chrDetail');
  detail.classList.add('active');

  var availH=window.innerHeight-140;
  var availW=detail.offsetWidth||800;
  var aspect=availW/availH;
  var cols=Math.round(Math.sqrt(fb*aspect));
  if(cols<3)cols=3;
  if(cols>fb)cols=fb;
  var numRows=Math.ceil(fb/cols);
  var rowH=Math.max(Math.floor(availH/numRows),3);

  // Chromosome picker
  var picker='<div class="chr-picker">';
  DATA.forEach(function(ch,i){{
    var cls=(i===chrIdx)?'chr-pick active':'chr-pick';
    var tc=ch.t2t?'var(--c-telo-ok)':'var(--muted)';
    picker+='<div class="'+cls+'" onclick="showChromosome('+i+')" title="'+ch.name+'">'+
      '<span class="chr-pick-dot" style="background:'+tc+'"></span>'+
      '<span class="chr-pick-name">'+ch.name.replace("chr","")+'</span></div>';
  }});
  picker+='</div>';

  var html=picker+
    '<div class="chr-detail-header">'+
    '<div class="back-btn" onclick="showChromosome('+chrIdx+')">&#8592; '+c.name+'</div>'+
    '<div class="chr-detail-title">'+(regionStart/1e6).toFixed(1)+' &ndash; '+(regionEnd/1e6).toFixed(1)+' Mb</div>'+
    '<div class="chr-detail-sub">'+c.name+' &middot; '+(( regionEnd-regionStart)/1e3).toFixed(0)+' kb region &middot; '+FINE_BIN_KB+' kb bins</div>'+
    '</div>';

  html+='<div class="chr-grid">';
  for(var row=0;row<numRows;row++){{
    var s=row*cols;
    var e=Math.min(s+cols,fb);
    var rowWidthPct=(e-s)/cols*100;
    html+='<div class="chr-grid-row" data-row="'+row+'" style="height:'+rowH+'px;width:'+rowWidthPct+'%">';
    var actualCols=e-s;
    for(var i=s;i<e;i++){{
      var d=fineBins[i];
      var lp=((i-s)/actualCols*100);
      var wp=(1/actualCols*100);
      for(var cat=0;cat<4;cat++){{
        var key=CAT_KEYS[cat];
        if(!LAYERS[key])continue;
        var val=d[cat];
        if(val<0.005)continue;
        html+='<div class="chr-grid-cell" style="'+
          'left:'+lp+'%;width:'+Math.max(wp,0.2)+'%;'+
          'top:'+(cat*25)+'%;height:25%;'+
          'background:var('+COLORS[cat]+');'+
          'opacity:'+Math.min(0.25+val*0.75,1)+
          '"></div>';
      }}
    }}
    html+='</div>';
  }}
  html+='</div>';

  html+='<div class="chr-grid-labels">'+
    '<span>'+(regionStart/1e6).toFixed(2)+' Mb</span>'+
    '<span>'+FINE_BIN_KB+' kb / bin</span>'+
    '<span>'+(regionEnd/1e6).toFixed(2)+' Mb</span>'+
    '</div>';

  detail.innerHTML=html;

  // Tooltip
  detail.querySelectorAll('.chr-grid-row').forEach(function(rowEl){{
    rowEl.addEventListener('mousemove',function(ev){{
      var rect=rowEl.getBoundingClientRect();
      var x=(ev.clientX-rect.left)/rect.width;
      var rowIdx=parseInt(rowEl.getAttribute('data-row'));
      var binIdx=rowIdx*cols+Math.floor(x*cols);
      if(binIdx>=fb)binIdx=fb-1;
      var pos=regionStart+(fineStart+binIdx)*bpPerFine-regionStart;
      var absPos=fineStart*bpPerFine+binIdx*bpPerFine;
      var d=fineBins[binIdx]||[0,0,0,0];
      var tt=document.getElementById('tooltip');
      tt.innerHTML=
        '<div class="tt-title">'+c.name+' : '+(absPos/1e6).toFixed(3)+' Mb</div>'+
        '<div class="tt-row"><span>&lt; 1 kb</span><span class="tt-val">'+(d[0]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>1-10 kb</span><span class="tt-val">'+(d[1]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>10-100 kb</span><span class="tt-val">'+(d[2]*100).toFixed(1)+'%</span></div>'+
        '<div class="tt-row"><span>&gt; 100 kb</span><span class="tt-val">'+(d[3]*100).toFixed(1)+'%</span></div>';
      tt.style.left=(ev.clientX+14)+'px';
      tt.style.top=(ev.clientY-10)+'px';
      tt.classList.add('visible');
    }});
    rowEl.addEventListener('mouseleave',function(){{
      document.getElementById('tooltip').classList.remove('visible');
    }});
    // Click to show repeats in this fine bin
    rowEl.addEventListener('click',function(ev){{
      var rect=rowEl.getBoundingClientRect();
      var x=(ev.clientX-rect.left)/rect.width;
      var rowIdx=parseInt(rowEl.getAttribute('data-row'));
      var rActual=Math.min(rowIdx*cols+cols,fb)-rowIdx*cols;
      var bi=rowIdx*cols+Math.floor(x*rActual);
      if(bi>=fb)bi=fb-1;
      var bStart=(fineStart+bi)*bpPerFine;
      var bEnd=bStart+bpPerFine;
      showRepeat(chrIdx,bStart,bEnd);
    }});
  }});

  setNavLevel('region');
}}

var MONOMERS={monomers_json};
var CURRENT_REPEAT=null;

function showRepeat(chrIdx, binStart, binEnd){{
    CURRENT_REPEAT={{chr:chrIdx,s:binStart,e:binEnd}};
    var c=DATA[chrIdx];
    document.getElementById('chrList').style.display='none';
    document.querySelector('.header').style.display='none';
    document.querySelector('.footer').style.display='none';
    var detail=document.getElementById('chrDetail');
    detail.classList.add('active');

    // Find arrays in this region from monomers data
    var regionArrays=[];
    for(var aid in MONOMERS){{
      var a=MONOMERS[aid];
      if(a.c===c.name && a.s<binEnd && a.e>binStart){{
        regionArrays.push({{id:aid,chr:a.c,start:a.s,end:a.e,period:a.p,monomers:a.m}});
      }}
    }}
    regionArrays.sort(function(a,b){{return a.start-b.start;}});

    var html='<div class="chr-detail-header">'+
      '<div class="back-btn" onclick="showRegion('+chrIdx+','+Math.floor(binStart/(BIN_KB*1000))+')">&#8592; Region</div>'+
      '<div class="chr-detail-title">Repeats</div>'+
      '<div class="chr-detail-sub">'+c.name+' : '+(binStart/1e6).toFixed(3)+' &ndash; '+(binEnd/1e6).toFixed(3)+' Mb &middot; '+regionArrays.length+' arrays</div>'+
      '</div>';

    if(regionArrays.length===0){{
      html+='<div style="color:var(--muted);font-size:13px;padding:20px">No arrays &ge; 10 kb in this region</div>';
    }}

    regionArrays.forEach(function(arr){{
      var arrLen=arr.end-arr.start;
      var nMono=arr.monomers.length;
      html+='<div class="repeat-card">'+
        '<div class="repeat-header">'+
          '<span class="repeat-title">'+arr.chr+':'+arr.start.toLocaleString()+'&ndash;'+arr.end.toLocaleString()+'</span>'+
          '<span class="repeat-info">'+(arrLen/1000).toFixed(1)+' kb &middot; period '+arr.period+' bp &middot; '+nMono+' monomers</span>'+
        '</div>'+
        '<div class="repeat-monomers">';

      // Render monomers as blocks
      if(nMono>0){{
        var maxMono=Math.max.apply(null,arr.monomers);
        arr.monomers.forEach(function(mlen,mi){{
          var w=Math.max(mlen/arrLen*100*nMono,0.3);
          // Color by deviation from period
          var dev=Math.abs(mlen-arr.period)/arr.period;
          var color=dev<0.1?'var(--c-telo-ok)':dev<0.3?'var(--c-mini)':'var(--c-macro)';
          html+='<div class="mono-block" style="width:'+w+'%;background:'+color+'" title="M'+mi+': '+mlen+' bp"></div>';
        }});
      }}

      html+='</div></div>';
    }});

    detail.innerHTML=html;
    setNavLevel('repeat');
}}

function showGenome(){{
  CURRENT_CHR=null;
  document.getElementById('chrList').style.display='';
  document.querySelector('.header').style.display='';
  document.querySelector('.footer').style.display='';
  var detail=document.getElementById('chrDetail');
  detail.classList.remove('active');
  detail.innerHTML='';
  setNavLevel('genome');
}}

function setNavLevel(level){{
  document.querySelectorAll('.nav-level').forEach(function(n){{
    n.classList.remove('active');
    if(n.getAttribute('data-level')===level)n.classList.add('active');
  }});
}}

function setLevel(el){{
  if(el.classList.contains('disabled'))return;
  var level=el.getAttribute('data-level');
  if(level==='genome')showGenome();
  else if(level==='chromosome')showChromosome(CURRENT_CHR!==null?CURRENT_CHR:0);
  else if(level==='region'){{
    var ci=CURRENT_CHR!==null?CURRENT_CHR:0;
    var midBin=Math.floor(DATA[ci].density.length/2);
    showRegion(ci,CURRENT_REGION!==null?CURRENT_REGION:midBin);
  }}
  else if(level==='repeat'&&CURRENT_REPEAT){{
    showRepeat(CURRENT_REPEAT.chr,CURRENT_REPEAT.s,CURRENT_REPEAT.e);
  }}
}}

// 12 distinct colors for families
var FAM_COLORS=['#0891b2','#7c3aed','#be123c','#059669','#b45309','#2563eb',
  '#db2777','#0d9488','#dc2626','#4338ca','#65a30d','#c026d3'];

function setView(el){{
  if(el.querySelector('.soon'))return;
  document.querySelectorAll('.view-toggle').forEach(function(v){{v.classList.remove('active');}});
  el.classList.add('active');
  CURRENT_VIEW=el.getAttribute('data-view');
  if(CURRENT_CHR!==null){{
    showChromosome(CURRENT_CHR);
  }}else{{
    render();
  }}
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
                                     monomers_path=None, families_path=None,
                                     assembly_name=""):
    """
    Create interactive chromosome visualization HTML.

    Args:
        fai_path: Path to .fai index file
        sat_path: Path to .sat file
        output_path: Path to output HTML file
        telomere_path: Path to telomeres.tsv (optional)
        its_path: Path to telomeres.its.bed (optional)
        monomers_path: Path to monomers.tsv (optional, for repeat detail view)
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

    logger.info("Binning repeat density (coarse + fine)...")
    bins, fine_bins = bin_repeats(sat_path, chroms)

    # Load monomers for repeat detail view (inline in HTML)
    monomers_data = {}
    if monomers_path and os.path.exists(monomers_path):
        monomers_data = _load_monomers_compact(monomers_path)

    # Load family assignments for "By family" view
    family_bins = {}
    tracked_families = []
    if families_path and os.path.exists(families_path):
        logger.info("Loading family assignments...")
        families, top_fams = load_families(families_path)
        family_bins, tracked_families = bin_families(families, chroms, top_fams)
        logger.info(f"Loaded {len(families)} family assignments, {len(tracked_families)} tracked families")

    logger.info("Generating visualization...")
    generate_chromosome_html(chroms, bins, fine_bins, telomeres, its, assembly_name,
                             output_path, monomers_data, family_bins, tracked_families)
