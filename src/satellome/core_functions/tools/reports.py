#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.10.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
HTML report generation with embedded image and interactive chart support.

Provides utilities for creating self-contained HTML reports from visualization
images and Plotly interactive charts. Converts image files (PNG, SVG) to
Base64-encoded data URIs and embeds them directly in HTML, eliminating
external file dependencies. Plotly HTML charts are embedded as iframes.

Functions:
    image_to_data_uri: Convert PNG image to Base64 data URI
    svg_to_data_uri: Convert SVG file to Base64 data URI
    create_html_report: Generate HTML report with embedded images from folder
"""

import base64
import logging
import os

logger = logging.getLogger(__name__)

def image_to_data_uri(image_path):
    """Convert PNG image file to Base64-encoded data URI."""
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        return f"data:image/png;base64,{encoded_string}"

def svg_to_data_uri(svg_path):
    """Convert SVG vector graphics file to Base64-encoded data URI."""
    with open(svg_path, "rb") as svg_file:
        encoded_string = base64.b64encode(svg_file.read()).decode('utf-8')
        return f"data:image/svg+xml;base64,{encoded_string}"


# All visualizations are now interactive Plotly HTML charts.
# Section definitions: (section_id, title, description, file_patterns)
# Patterns are matched against filenames (without directory)
REPORT_SECTIONS = []


# Chart sections: (section_id, title, description, filename_keywords)
# Charts are grouped by matching keywords in their filenames
CHART_SECTIONS = [
    (
        "karyotypes",
        "Chromosome Karyotypes",
        "Genome-wide distribution of tandem repeat families across chromosomes. "
        "Enhanced views amplify small arrays for visibility.",
        [".raw.", ".nosing.", ".sing."],
        [".gaps.", ".repeats."],
    ),
    (
        "gaps",
        "Gap Analysis",
        "Assembly gaps and their proximity to tandem repeat arrays. "
        "Gap-adjacent repeats may indicate unresolved satellite regions.",
        [".gaps.", ".repeats.with.gaps.", ".repeats.nogaps."],
        [],
    ),
    (
        "scatter",
        "Repeat Characteristics",
        "Scatter plots showing relationships between GC content, period length, "
        "and percent match across all detected tandem repeat arrays.",
        [".3D.", ".2D."],
        [],
    ),
    (
        "flow",
        "Classification Flow",
        "Sankey diagram showing how tandem repeats are classified into families "
        "and the relative abundance of each category.",
        ["trs_flow", "sankey"],
        [],
    ),
]


def _classify_charts(image_folder):
    """Classify HTML chart files into sections based on filename patterns.

    Returns:
        list of (section_id, title, desc, chart_files) and list of uncategorized charts
    """
    if not os.path.isdir(image_folder):
        return [], []

    html_files = sorted([
        f for f in os.listdir(image_folder)
        if f.endswith('.html')
    ])

    classified = {}  # section_id -> list of (fpath, display_name)
    used = set()

    for section_id, title, desc, include_keywords, exclude_keywords in CHART_SECTIONS:
        section_files = []
        for fname in html_files:
            if fname in used:
                continue
            fname_lower = fname.lower()
            matches = any(kw in fname_lower for kw in include_keywords)
            excluded = any(kw in fname_lower for kw in exclude_keywords) if exclude_keywords else False
            if matches and not excluded:
                fpath = os.path.join(image_folder, fname)
                display_name = os.path.splitext(fname)[0].replace('.', ' ').replace('_', ' ')
                section_files.append((fpath, display_name, fname))
                used.add(fname)
        if section_files:
            classified[section_id] = (section_id, title, desc, section_files)

    # Uncategorized
    uncategorized = []
    for fname in html_files:
        if fname not in used:
            fpath = os.path.join(image_folder, fname)
            display_name = os.path.splitext(fname)[0].replace('.', ' ').replace('_', ' ')
            uncategorized.append((fpath, display_name, fname))

    sections = [classified[sid] for sid, _, _, _, _ in CHART_SECTIONS if sid in classified]
    return sections, uncategorized


def _build_chart_html(fpath, display_name, anim_delay):
    """Build HTML for a single chart card with lazy-loaded iframe."""
    with open(fpath, 'r', encoding='utf-8', errors='replace') as f:
        chart_content = f.read()
    encoded = base64.b64encode(chart_content.encode('utf-8')).decode('utf-8')
    return f'''
            <div class="chart-card" style="animation-delay: {anim_delay:.2f}s;">
              <div class="chart-label">{display_name}</div>
              <iframe srcdoc="" data-src="{encoded}" class="chart-iframe lazy-iframe" loading="lazy"></iframe>
            </div>
'''


def _load_results_yaml(yaml_path):
    """Load results.yaml and return repeat stats dict."""
    import yaml
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    repeats = {}
    try:
        dataset = data['work_files']['ref_assembly_name_for_trf']
        repeats = data['work_files']['repeats'][dataset]['trevis']
    except (KeyError, TypeError):
        pass
    return repeats


def _build_summary_html(assembly_name, genome_size, repeats):
    """Build the summary statistics section HTML."""

    # Repeat categories with display names and descriptions
    categories = [
        ("pmicro", "Perfect Microsatellites", "Period < 6 bp, 100% match"),
        ("micro", "Microsatellites", "Period < 6 bp"),
        ("tSSR", "True SSR", "Simple sequence repeats"),
        ("fSSR", "Fuzzy SSR", "Imperfect simple repeats"),
        ("compex", "Complex Repeats", "Period > 9 bp"),
        ("1kb", "Arrays > 1 kb", ""),
        ("10kb", "Arrays > 10 kb", ""),
        ("100kb", "Arrays > 100 kb", ""),
        ("1000kb", "Arrays > 1 Mb", ""),
    ]

    genome_size_display = f"{genome_size:,} bp" if genome_size else "N/A"
    genome_mb = f"{genome_size / 1e6:.1f} Mb" if genome_size else ""

    rows_html = ""
    for key, display_name, desc in categories:
        if key not in repeats:
            continue
        r = repeats[key]
        n = r.get('n', 0)
        pgenome = r.get('pgenome', 0)
        bar_width = min(pgenome * 10, 100)  # scale for visual bar
        rows_html += f'''
        <tr>
          <td class="stat-name">{display_name}<span class="stat-desc">{desc}</span></td>
          <td class="stat-count">{n:,}</td>
          <td class="stat-pct">
            <div class="pct-bar-wrap">
              <div class="pct-bar" style="width: {bar_width}%"></div>
            </div>
            <span>{pgenome}%</span>
          </td>
        </tr>'''

    if not rows_html:
        return ""

    return f'''
    <section id="summary" class="section" style="animation-delay: 0.1s;">
      <h2 class="section-title">Assembly Summary</h2>
      <div class="summary-grid">
        <div class="summary-card">
          <div class="summary-label">Assembly</div>
          <div class="summary-value">{assembly_name}</div>
        </div>
        <div class="summary-card">
          <div class="summary-label">Genome Size</div>
          <div class="summary-value">{genome_mb}</div>
          <div class="summary-sub">{genome_size_display}</div>
        </div>
      </div>
      <div class="stats-table-wrap">
        <table class="stats-table">
          <thead>
            <tr>
              <th>Category</th>
              <th>Count</th>
              <th>% Genome</th>
            </tr>
          </thead>
          <tbody>
            {rows_html}
          </tbody>
        </table>
      </div>
    </section>
'''


def _generate_report_html(sections, uncategorized, taxon_name=None,
                          assembly_name=None, genome_size=0, results_yaml=None):
    """Generate the full HTML report string."""

    display_name = assembly_name or taxon_name or "Satellome Report"
    title = f"Satellome Report — {display_name}"

    sections_html = ""
    section_nav = ""
    anim_delay = 0.1

    # Build summary table from results.yaml
    summary_html = ""
    if results_yaml and os.path.exists(results_yaml):
        repeats = _load_results_yaml(results_yaml)
        if repeats:
            summary_html = _build_summary_html(
                assembly_name or display_name, genome_size, repeats
            )
            section_nav += '<a href="#summary" class="nav-link">Summary</a>\n'
            anim_delay = 0.3

    for section_id, section_title, section_desc, items in sections:
        section_nav += f'<a href="#{section_id}" class="nav-link">{section_title}</a>\n'

        charts_html = ""
        for fpath, display_name, fname in items:
            if fname.endswith('.html'):
                charts_html += _build_chart_html(fpath, display_name, anim_delay)
            elif fname.endswith('.svg'):
                data_uri = svg_to_data_uri(fpath)
                charts_html += f'''
            <div class="image-card" style="animation-delay: {anim_delay:.2f}s;">
              <div class="image-label">{display_name}</div>
              <img src="{data_uri}" alt="{display_name}" loading="lazy">
            </div>
'''
            elif fname.endswith('.png'):
                data_uri = image_to_data_uri(fpath)
                charts_html += f'''
            <div class="image-card" style="animation-delay: {anim_delay:.2f}s;">
              <div class="image-label">{display_name}</div>
              <img src="{data_uri}" alt="{display_name}" loading="lazy">
            </div>
'''
            anim_delay += 0.05

        sections_html += f'''
    <section id="{section_id}" class="section" style="animation-delay: {anim_delay:.2f}s;">
      <h2 class="section-title">{section_title}</h2>
      <p class="section-desc">{section_desc}</p>
      <div class="charts-grid">
        {charts_html}
      </div>
    </section>
'''
        anim_delay += 0.1

    # Uncategorized charts
    if uncategorized:
        section_nav += '<a href="#other" class="nav-link">Other</a>\n'
        other_html = ""
        for fpath, display_name, fname in uncategorized:
            if fname.endswith('.html'):
                other_html += _build_chart_html(fpath, display_name, anim_delay)
            elif fname.endswith(('.png', '.svg')):
                data_uri = svg_to_data_uri(fpath) if fname.endswith('.svg') else image_to_data_uri(fpath)
                other_html += f'''
            <div class="image-card" style="animation-delay: {anim_delay:.2f}s;">
              <div class="image-label">{display_name}</div>
              <img src="{data_uri}" alt="{display_name}" loading="lazy">
            </div>
'''
            anim_delay += 0.05

        sections_html += f'''
    <section id="other" class="section" style="animation-delay: {anim_delay:.2f}s;">
      <h2 class="section-title">Additional Visualizations</h2>
      <div class="charts-grid">
        {other_html}
      </div>
    </section>
'''

    return f'''<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;500;700&family=Playfair+Display:wght@400;700;900&family=Source+Sans+3:wght@300;400;600&display=swap" rel="stylesheet">
<style>
  /* === LIGHT THEME (default) === */
  :root, [data-theme="light"] {{
    --bg-deep: #f8f6f1;
    --bg-card: #ffffff;
    --bg-card-hover: #fafaf8;
    --border-subtle: #e2ddd5;
    --text-primary: #1a1a1a;
    --text-secondary: #555555;
    --text-muted: #888888;
    --accent-cyan: #0891b2;
    --accent-emerald: #059669;
    --overlay-soft: rgba(0,0,0,0.02);
    --overlay-medium: rgba(0,0,0,0.03);
    --overlay-strong: rgba(0,0,0,0.04);
    --shadow-card: rgba(0,0,0,0.08);
    --grain-opacity: 0.015;
    --title-start: #1a1a1a;
    --smart-bg: rgba(8, 145, 178, 0.05);
    --smart-border: rgba(8, 145, 178, 0.15);
    --nav-bg: rgba(255,255,255,0.85);
    --nav-border: rgba(0,0,0,0.06);
  }}

  /* === DARK THEME === */
  [data-theme="dark"] {{
    --bg-deep: #0a0e17;
    --bg-card: #111827;
    --bg-card-hover: #1a2332;
    --border-subtle: #1e293b;
    --text-primary: #e2e8f0;
    --text-secondary: #94a3b8;
    --text-muted: #64748b;
    --accent-cyan: #22d3ee;
    --accent-emerald: #34d399;
    --overlay-soft: rgba(255,255,255,0.02);
    --overlay-medium: rgba(255,255,255,0.03);
    --overlay-strong: rgba(255,255,255,0.04);
    --shadow-card: rgba(0,0,0,0.3);
    --grain-opacity: 0.03;
    --title-start: #e2e8f0;
    --smart-bg: rgba(34, 211, 238, 0.04);
    --smart-border: rgba(34, 211, 238, 0.12);
    --nav-bg: rgba(17,24,39,0.9);
    --nav-border: rgba(255,255,255,0.06);
  }}

  * {{ margin: 0; padding: 0; box-sizing: border-box; }}

  body {{
    background: var(--bg-deep);
    color: var(--text-primary);
    font-family: 'Source Sans 3', sans-serif;
    min-height: 100vh;
    overflow-x: hidden;
  }}

  /* === GRAIN OVERLAY === */
  .grain {{
    position: fixed;
    top: -50%; left: -50%;
    width: 200%; height: 200%;
    pointer-events: none;
    z-index: 999;
    opacity: var(--grain-opacity);
    background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 256 256' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='noise'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.9' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23noise)'/%3E%3C/svg%3E");
  }}

  /* === THEME TOGGLE === */
  .theme-toggle {{
    position: fixed;
    top: 24px;
    right: 24px;
    z-index: 1001;
    display: flex;
    align-items: center;
    gap: 10px;
    padding: 6px 14px;
    border-radius: 40px;
    background: var(--bg-card);
    border: 1px solid var(--border-subtle);
    cursor: pointer;
    transition: all 0.3s ease;
    box-shadow: 0 2px 12px var(--shadow-card);
    user-select: none;
  }}

  .theme-toggle:hover {{
    transform: translateY(-1px);
    box-shadow: 0 4px 16px var(--shadow-card);
  }}

  .theme-toggle .toggle-icon {{
    font-size: 16px;
    line-height: 1;
    transition: transform 0.3s ease;
  }}

  .theme-toggle:hover .toggle-icon {{ transform: rotate(20deg); }}

  .theme-toggle .toggle-label {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 11px;
    color: var(--text-muted);
    letter-spacing: 0.5px;
  }}

  /* === NAV SIDEBAR === */
  .nav {{
    position: fixed;
    top: 50%;
    left: 20px;
    transform: translateY(-50%);
    z-index: 1000;
    display: flex;
    flex-direction: column;
    gap: 6px;
    background: var(--nav-bg);
    backdrop-filter: blur(12px);
    -webkit-backdrop-filter: blur(12px);
    border: 1px solid var(--nav-border);
    border-radius: 10px;
    padding: 10px 12px;
    animation: fadeUp 0.6s ease-out 0.3s both;
  }}

  .nav-link {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 10px;
    letter-spacing: 0.5px;
    color: var(--text-muted);
    text-decoration: none;
    padding: 4px 8px;
    border-radius: 4px;
    transition: all 0.2s ease;
    white-space: nowrap;
  }}

  .nav-link:hover {{
    color: var(--accent-cyan);
    background: var(--overlay-soft);
  }}

  /* === CONTAINER === */
  .container {{
    position: relative;
    z-index: 1;
    max-width: 1100px;
    margin: 0 auto;
    padding: 60px 40px 100px;
  }}

  /* === HEADER === */
  .header {{
    text-align: center;
    margin-bottom: 70px;
    animation: fadeUp 0.8s ease-out;
  }}

  .header-label {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 11px;
    font-weight: 500;
    letter-spacing: 4px;
    text-transform: uppercase;
    color: var(--accent-cyan);
    margin-bottom: 16px;
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 12px;
  }}

  .header-label::before,
  .header-label::after {{
    content: '';
    width: 40px;
    height: 1px;
    background: linear-gradient(90deg, transparent, var(--accent-cyan));
  }}

  .header-label::after {{
    background: linear-gradient(90deg, var(--accent-cyan), transparent);
  }}

  .header h1 {{
    font-family: 'Playfair Display', serif;
    font-size: clamp(32px, 4vw, 48px);
    font-weight: 900;
    letter-spacing: -1px;
    line-height: 1.1;
    background: linear-gradient(135deg, var(--title-start) 0%, var(--accent-cyan) 50%, var(--accent-emerald) 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    margin-bottom: 16px;
  }}

  .header p {{
    font-size: 17px;
    color: var(--text-secondary);
    max-width: 600px;
    margin: 0 auto;
    line-height: 1.6;
    font-weight: 300;
  }}

  /* === SECTIONS === */
  .section {{
    margin-bottom: 64px;
    opacity: 0;
    transform: translateY(24px);
    animation: fadeUp 0.6s ease-out forwards;
  }}

  .section-title {{
    font-family: 'Playfair Display', serif;
    font-size: 28px;
    font-weight: 700;
    margin-bottom: 8px;
    color: var(--text-primary);
  }}

  .section-desc {{
    font-size: 14px;
    color: var(--text-secondary);
    font-weight: 300;
    line-height: 1.6;
    margin-bottom: 28px;
    max-width: 700px;
  }}

  /* === IMAGE GRID === */
  .image-grid {{
    display: flex;
    flex-direction: column;
    gap: 24px;
  }}

  .image-card {{
    background: var(--bg-card);
    border: 1px solid var(--border-subtle);
    border-radius: 12px;
    overflow: hidden;
    transition: all 0.35s ease;
    position: relative;
    opacity: 0;
    transform: translateY(20px);
    animation: fadeUp 0.5s ease-out forwards;
  }}

  .image-card::before {{
    content: '';
    position: absolute;
    top: 0; left: 0; right: 0;
    height: 2px;
    background: linear-gradient(90deg, transparent, var(--accent-cyan), transparent);
    opacity: 0;
    transition: opacity 0.35s ease;
  }}

  .image-card:hover {{
    border-color: var(--overlay-medium);
    box-shadow: 0 8px 32px var(--shadow-card);
    transform: translateY(-2px);
  }}

  .image-card:hover::before {{ opacity: 1; }}

  .image-label {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 11px;
    font-weight: 500;
    letter-spacing: 1px;
    text-transform: uppercase;
    color: var(--text-muted);
    padding: 16px 20px 8px;
  }}

  .image-card img {{
    width: 100%;
    display: block;
    padding: 4px 16px 16px;
    background: var(--bg-card);
  }}

  /* === CHARTS GRID === */
  .charts-grid {{
    display: flex;
    flex-direction: column;
    gap: 24px;
  }}

  .chart-card {{
    background: var(--bg-card);
    border: 1px solid var(--border-subtle);
    border-radius: 12px;
    overflow: hidden;
    transition: all 0.35s ease;
    opacity: 0;
    transform: translateY(20px);
    animation: fadeUp 0.5s ease-out forwards;
  }}

  .chart-card:hover {{
    border-color: var(--overlay-medium);
    box-shadow: 0 8px 32px var(--shadow-card);
  }}

  .chart-label {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 11px;
    font-weight: 500;
    letter-spacing: 1px;
    text-transform: uppercase;
    color: var(--text-muted);
    padding: 16px 20px 8px;
  }}

  .chart-iframe {{
    width: 100%;
    height: 600px;
    border: none;
    display: block;
  }}

  /* === SUMMARY === */
  .summary-grid {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 16px;
    margin-bottom: 28px;
  }}

  .summary-card {{
    background: var(--bg-card);
    border: 1px solid var(--border-subtle);
    border-radius: 10px;
    padding: 20px 24px;
  }}

  .summary-label {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 10px;
    font-weight: 500;
    letter-spacing: 1.5px;
    text-transform: uppercase;
    color: var(--text-muted);
    margin-bottom: 6px;
  }}

  .summary-value {{
    font-family: 'Playfair Display', serif;
    font-size: 22px;
    font-weight: 700;
    color: var(--text-primary);
  }}

  .summary-sub {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 11px;
    color: var(--text-muted);
    margin-top: 4px;
  }}

  .stats-table-wrap {{
    background: var(--bg-card);
    border: 1px solid var(--border-subtle);
    border-radius: 10px;
    overflow: hidden;
  }}

  .stats-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 14px;
  }}

  .stats-table thead th {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 10px;
    font-weight: 500;
    letter-spacing: 1px;
    text-transform: uppercase;
    color: var(--text-muted);
    text-align: left;
    padding: 14px 20px 10px;
    border-bottom: 1px solid var(--border-subtle);
  }}

  .stats-table tbody tr {{
    transition: background 0.15s;
  }}

  .stats-table tbody tr:hover {{
    background: var(--overlay-soft);
  }}

  .stats-table td {{
    padding: 12px 20px;
    border-bottom: 1px solid var(--border-subtle);
    color: var(--text-secondary);
  }}

  .stats-table tbody tr:last-child td {{
    border-bottom: none;
  }}

  .stat-name {{
    color: var(--text-primary);
    font-weight: 400;
  }}

  .stat-desc {{
    display: block;
    font-size: 11px;
    color: var(--text-muted);
    font-weight: 300;
    margin-top: 2px;
  }}

  .stat-count {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 13px;
    text-align: right;
    white-space: nowrap;
  }}

  .stat-pct {{
    display: flex;
    align-items: center;
    gap: 10px;
    min-width: 160px;
  }}

  .stat-pct span {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 12px;
    white-space: nowrap;
    min-width: 45px;
    text-align: right;
  }}

  .pct-bar-wrap {{
    flex: 1;
    height: 6px;
    background: var(--overlay-soft);
    border-radius: 3px;
    overflow: hidden;
  }}

  .pct-bar {{
    height: 100%;
    background: linear-gradient(90deg, var(--accent-cyan), var(--accent-emerald));
    border-radius: 3px;
    transition: width 0.4s ease;
  }}

  @media (max-width: 860px) {{
    .summary-grid {{ grid-template-columns: 1fr; }}
    .stat-pct {{ min-width: 100px; }}
  }}

  /* === FOOTER === */
  .footer {{
    text-align: center;
    margin-top: 60px;
    padding-top: 32px;
    border-top: 1px solid var(--border-subtle);
    opacity: 0;
    animation: fadeUp 0.6s ease-out 1.5s forwards;
  }}

  .footer p {{
    font-size: 12px;
    color: var(--text-muted);
    font-weight: 300;
  }}

  .footer a {{
    color: var(--accent-cyan);
    text-decoration: none;
  }}

  /* === ANIMATIONS === */
  @keyframes fadeUp {{
    from {{ opacity: 0; transform: translateY(24px); }}
    to {{ opacity: 1; transform: translateY(0); }}
  }}

  /* === RESPONSIVE === */
  @media (max-width: 860px) {{
    .container {{ padding: 40px 16px 60px; }}
    .nav {{ display: none; }}
    .image-card img {{ padding: 4px 8px 8px; }}
  }}

  /* === BACK TO TOP === */
  .back-top {{
    position: fixed;
    bottom: 24px;
    right: 24px;
    z-index: 1000;
    width: 40px;
    height: 40px;
    border-radius: 50%;
    background: var(--bg-card);
    border: 1px solid var(--border-subtle);
    color: var(--text-muted);
    font-size: 18px;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
    box-shadow: 0 2px 12px var(--shadow-card);
    opacity: 0;
    transition: all 0.3s ease;
    pointer-events: none;
  }}

  .back-top.visible {{
    opacity: 1;
    pointer-events: auto;
  }}

  .back-top:hover {{
    transform: translateY(-2px);
    color: var(--accent-cyan);
  }}
</style>
</head>
<body>

<div class="grain"></div>

<div class="theme-toggle" onclick="toggleTheme()" title="Switch theme">
  <span class="toggle-icon" id="themeIcon">&#9790;</span>
  <span class="toggle-label" id="themeLabel">dark</span>
</div>

<nav class="nav">
  {section_nav}
</nav>

<button class="back-top" id="backTop" onclick="window.scrollTo({{top:0,behavior:'smooth'}})" title="Back to top">&#8593;</button>

<div class="container">

  <header class="header">
    <div class="header-label">Analysis Report</div>
    <h1>{title}</h1>
    <p>Tandem repeat detection, classification, and genome-wide distribution analysis.</p>
  </header>

  {summary_html}

  {sections_html}

  <footer class="footer">
    <p>Generated by <a href="https://github.com/ad3002/satellome">Satellome</a></p>
  </footer>

</div>

<script>
// Theme toggle
function toggleTheme() {{
  var html = document.documentElement;
  var icon = document.getElementById('themeIcon');
  var label = document.getElementById('themeLabel');
  var isDark = html.getAttribute('data-theme') === 'dark';

  if (isDark) {{
    html.setAttribute('data-theme', 'light');
    icon.innerHTML = '\\u263E';
    label.textContent = 'dark';
    localStorage.setItem('satellome-theme', 'light');
  }} else {{
    html.setAttribute('data-theme', 'dark');
    icon.innerHTML = '\\u263C';
    label.textContent = 'light';
    localStorage.setItem('satellome-theme', 'dark');
  }}
}}

// Restore saved theme
(function() {{
  var saved = localStorage.getItem('satellome-theme');
  if (saved === 'dark') {{
    document.documentElement.setAttribute('data-theme', 'dark');
    document.getElementById('themeIcon').innerHTML = '\\u263C';
    document.getElementById('themeLabel').textContent = 'light';
  }}
}})();

// Back to top button visibility
window.addEventListener('scroll', function() {{
  var btn = document.getElementById('backTop');
  if (window.scrollY > 400) {{
    btn.classList.add('visible');
  }} else {{
    btn.classList.remove('visible');
  }}
}});

// Lazy-load iframes for interactive charts
(function() {{
  var observer = new IntersectionObserver(function(entries) {{
    entries.forEach(function(entry) {{
      if (entry.isIntersecting) {{
        var iframe = entry.target;
        var encoded = iframe.getAttribute('data-src');
        if (encoded) {{
          var decoded = atob(encoded);
          iframe.srcdoc = decoded;
          iframe.removeAttribute('data-src');
        }}
        observer.unobserve(iframe);
      }}
    }});
  }}, {{ rootMargin: '200px' }});

  document.querySelectorAll('.lazy-iframe').forEach(function(iframe) {{
    observer.observe(iframe);
  }});
}})();
</script>

</body>
</html>'''


def create_html_report(image_folder, report_file, taxon_name=None,
                       assembly_name=None, genome_size=0, results_yaml=None):
    """
    Generate self-contained HTML report with embedded images and charts.

    Args:
        image_folder (str): Path to folder containing image and chart files
        report_file (str): Path to output HTML file to create
        taxon_name (str, optional): Species/taxon name for report title
        assembly_name (str, optional): Input assembly filename
        genome_size (int, optional): Total genome size in bp
        results_yaml (str, optional): Path to results.yaml with classification stats
    """
    # Classify charts into sections
    sections, uncategorized = _classify_charts(image_folder)

    total_items = sum(len(items) for _, _, _, items in sections) + len(uncategorized)

    # Generate and write report (even if no charts, we may have summary)
    html = _generate_report_html(
        sections, uncategorized, taxon_name,
        assembly_name=assembly_name, genome_size=genome_size,
        results_yaml=results_yaml,
    )

    os.makedirs(os.path.dirname(report_file), exist_ok=True)
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(html)

    logger.info("HTML report created successfully!")
    logger.info(f"File: {report_file}")
    logger.info(f"Embedded: {total_items} visualizations")
