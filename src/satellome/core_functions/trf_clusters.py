#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 11.01.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import os
import logging
from collections import Counter

import plotly.express as px
import plotly.graph_objects as go
from intervaltree import IntervalTree
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server environments
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Configure logging
logger = logging.getLogger(__name__)

from satellome.core_functions.trf_drawing import (get_gaps_annotation, read_trf_file,
                                  scaffold_length_sort_length)
from satellome.core_functions.trf_embedings import get_disances
from satellome.constants import (
    CANVAS_WIDTH_DEFAULT, CANVAS_HEIGHT_DEFAULT, CANVAS_HEIGHT_MIN, CANVAS_HEIGHT_MAX,
    CHROMOSOME_HEIGHT, VERTICAL_SPACER, BASE_HEIGHT,
    MARGIN_TOP, MARGIN_BOTTOM, MARGIN_LEFT, MARGIN_RIGHT,
    ENHANCE_DEFAULT, ENHANCE_LARGE, GAP_CUTOFF_DEFAULT, GAP_SEARCH_WINDOW,
    SAMPLE_SIZE_FOR_CLUSTERING, START_CUTOFF_MAX,
    TR_CUTOFF_LARGE, MIN_SCAFFOLD_LENGTH_FILTER,
    SEPARATOR_LINE, TR_SIZE_RANGES, GAP_SIZE_RANGES,
    RECURSION_LIMIT_DEFAULT
)

import sys
sys.setrecursionlimit(RECURSION_LIMIT_DEFAULT)


def safe_write_figure(fig, output_file, width=None, height=None, engine="kaleido"):
    """
    Write plotly figure to HTML file.

    Static PNG export is handled separately by matplotlib (via _create_matplotlib_karyotype).
    This function only creates interactive HTML visualizations.

    Args:
        fig: Plotly figure object
        output_file: Output file path (extension will be changed to .html)
        width: Figure width (optional, ignored - kept for compatibility)
        height: Figure height (optional, ignored - kept for compatibility)
        engine: Image export engine (ignored - kept for compatibility)

    Returns:
        str: Path to the created HTML file
    """
    # Change extension to .html
    html_file = os.path.splitext(output_file)[0] + ".html"
    fig.write_html(html_file)
    logger.debug(f"Exported interactive plot to {html_file}")
    return html_file


def save_matplotlib_figure(fig, output_file, dpi=150):
    """
    Save matplotlib figure to file (PNG or SVG).

    Args:
        fig: Matplotlib figure object
        output_file: Output file path
        dpi: DPI for raster formats (default: 150)

    Returns:
        str: Path to the created file
    """
    try:
        ext = os.path.splitext(output_file)[1].lower()
        if ext in ['.png', '.jpg', '.jpeg']:
            fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        elif ext == '.svg':
            fig.savefig(output_file, format='svg', bbox_inches='tight')
        elif ext == '.pdf':
            fig.savefig(output_file, format='pdf', bbox_inches='tight')
        else:
            # Default to PNG
            output_file = os.path.splitext(output_file)[0] + '.png'
            fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')

        plt.close(fig)
        logger.info(f"Exported matplotlib plot to {output_file}")
        return output_file
    except Exception as e:
        logger.error(f"Failed to save matplotlib figure: {e}")
        plt.close(fig)
        return None


class Graph:

    # init function to declare class variables
    def __init__(self, V):
        self.V = []
        self.id2node = {}
        self.node2id = {}
        for i, id1 in enumerate(V):
            self.id2node[id1] = i
            self.node2id[i] = id1
            self.V.append(i)
        self.adj = [[] for i in V]

    def DFSUtil(self, temp, v, visited):

        # Mark the current vertex as visited
        visited[v] = True

        # Store the vertex to list
        temp.append(v)

        # Repeat for all vertices adjacent
        # to this vertex v
        for i, dist in self.adj[v]:
            if visited[i] == False:
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp

    # method to add an undirected edge
    def addEdge(self, v, w, dist):
        id1 = self.id2node[v]
        id2 = self.id2node[w]
        self.adj[id1].append((id2, dist))
        self.adj[id2].append((id1, dist))

    def remove_edges_by_distances(self, cutoff):
        for id1 in self.V:
            new_adj = []
            for id2, dist in self.adj[id1]:
                if dist < cutoff:
                    new_adj.append((id2, dist))
            self.adj[id1] = new_adj

    # Method to retrieve connected components
    # in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in self.V:
            visited.append(False)
        for v in self.V:
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc


def name_clusters(distances, tr2vector, df_trs, level=1):

    all_distances = list(distances.values())
    all_distances = list(set(map(int, all_distances)))
    all_distances.sort(reverse=True)
    start_cutoff = min(int(all_distances[0]), START_CUTOFF_MAX)

    G = Graph(list(tr2vector.keys()))
    for (id1, id2) in distances:
        G.addEdge(id1, id2, distances[(id1, id2)])

    for i in tqdm(range(start_cutoff, level - 1, -1), desc="Naming clusters"):
        G.remove_edges_by_distances(i)
        comps = G.connectedComponents()
        items = []
        singl = []
        for c in comps:
            ids = [(G.node2id[id1], df_trs[G.node2id[id1]].get("period")) for id1 in c]
            if len(ids) > 3:  # Why 3? It should be 1
                items.append(ids)
            else:
                singl += ids
        items.sort(key=lambda x: len(x))
        for class_name, d in enumerate(items):
            median_monomer = [x[1] for x in d]
            median_monomer.sort()
            median_monomer = median_monomer[int(len(median_monomer) / 2)]

            name = f"{class_name}_{median_monomer}"
            for id1, period in d:
                df_trs[id1]["family_name"] = name
                df_trs[id1]["locus_name"] = name

        for id1, period in singl:
            df_trs[id1]["family_name"] = "SING"

    return df_trs, tr2vector, distances, all_distances


def _draw_sankey(output_file_name, title_text, labels, source, target, value):
    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=1),
                    label=labels,
                    color="blue",
                ),
                link=dict(
                    source=source,
                    target=target,
                    value=value,
                ),
            )
        ]
    )
    fig.update_layout(
        title_text=title_text,
        font_size=10,
        height=2000,
        width=2000,
    )
    safe_write_figure(fig, output_file_name)


def _deprecated_draw_sankey(
    output_file_name,
    title_text,
    df_trs,
    tr2vector,
    distances,
    all_distances,
    skip_singletons=True,
):

    G = Graph(list(tr2vector.keys()))
    for (id1, id2) in distances:
        G.addEdge(id1, id2, distances[(id1, id2)])

    steps = []

    start_cutoff = int(all_distances[0])

    name2trs = {}
    last_n_comp = 0
    id2names = {}

    for i in tqdm(range(start_cutoff, 0, -1), desc="Cut distances"):
        G.remove_edges_by_distances(i)
        comps = G.connectedComponents()
        if len(comps) == last_n_comp:
            continue
        last_n_comp = len(comps)

        items = []
        singl = []

        id2INstep = {}
        name2size = {}
        name2ids = {}
        name2id = {}

        for c in comps:
            ids = [(G.node2id[id1], df_trs[G.node2id[id1]].get("period")) for id1 in c]
            if len(ids) > 3:
                items.append(ids)
            else:
                singl += ids

        items.sort(key=lambda x: len(x))
        for class_name, d in enumerate(items):
            median_monomer = [x[1] for x in d]
            median_monomer.sort()
            median_monomer = median_monomer[int(len(median_monomer) / 2)]
            name = f"{i}_{class_name}_{median_monomer}"
            for id1, period in d:
                id2INstep[id1] = name
                name2id[name] = id1
            name2size[name] = len(d)
            name2ids[name] = d
            name2trs[name] = d

            for id_, _ in d:
                id2names.setdefault(id_, [])
                id2names[id_].append(name)

        if not skip_singletons and singl:
            name = f"{i}_SING"
            for id1, period in singl:
                id2INstep[id1] = name
                name2id[name] = id1

                id2names.setdefault(id1, [])
                id2names[id1].append(name)

            name2size[name] = len(singl)
            name2ids[name] = singl

        steps.append((id2INstep, name2size, name2ids, name2id))

    labels = []
    source = []
    target = []
    value = []
    name2monomers = {}
    name2lid = {}
    lid = 0
    prev_id2INstep, prev_name2size, name2ids, prev_name2id = steps[0]
    for name in prev_name2size:
        labels.append(name)
        name2lid[name] = lid
        lid += 1

        name2monomers[name] = name2ids[name]

    for id2INstep, name2size, name2ids, name2id in steps[1:]:

        for name in name2size:
            labels.append(name)
            name2lid[name] = lid
            lid += 1

            start = name2lid[prev_id2INstep[name2id[name]]]
            end = name2lid[name]

            source.append(start)
            target.append(end)
            value.append(name2size[name])

            name2monomers[name] = name2ids[name]

        prev_id2INstep = id2INstep

    _draw_sankey(output_file_name, title_text, labels, source, target, value)

    return name2monomers, name2lid, name2trs, id2names


def _deprecated_draw_spheres(output_file_name_prefix, title_text, df_trs):

    fig = px.scatter_3d(
        df_trs,
        x="gc",
        y="period",
        z="pmatch",
        color="family_name",
        size="log_length",
    )
    fig.update_layout(
        title={
            "text": title_text,
            "y": 0.99,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        }
    )
    fig.update_layout(width=CANVAS_WIDTH_DEFAULT, height=CANVAS_HEIGHT_DEFAULT)
    output_file_name = output_file_name_prefix + ".3D.svg"
    safe_write_figure(fig, output_file_name)

    fig = px.scatter_3d(
        df_trs[df_trs["family_name"] != "SING"],
        x="gc",
        y="period",
        z="pmatch",
        color="family_name",
        size="log_length",
    )
    fig.update_layout(
        title={
            "text": title_text + " No Singletons",
            "y": 0.99,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        }
    )
    fig.update_layout(width=CANVAS_WIDTH_DEFAULT, height=CANVAS_HEIGHT_DEFAULT)
    output_file_name = output_file_name_prefix + ".3D.nosingl.svg"
    safe_write_figure(fig, output_file_name)

    fig = px.scatter(df_trs, x="gc", y="period", color="family_name", size="log_length")
    output_file_name = output_file_name_prefix + ".2D.gc_period.svg"
    safe_write_figure(fig, output_file_name)

    fig = px.scatter(df_trs, x="gc", y="period", color="family_name", size="log_length")
    output_file_name = output_file_name_prefix + ".2D.gc_period.svg"
    safe_write_figure(fig, output_file_name)

    fig = px.scatter(df_trs, x="gc", y="pmatch", color="family_name", size="log_length")
    output_file_name = output_file_name_prefix + ".2D.gc_pmatch.svg"
    safe_write_figure(fig, output_file_name)

    fig = px.scatter(
        df_trs, x="pmatch", y="period", color="family_name", size="log_length"
    )
    output_file_name = output_file_name_prefix + ".2D.period_period.svg"
    safe_write_figure(fig, output_file_name)


def _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=False):

    if use_chrm:
        scaffold_items = scaffold_for_plot["chrm"]
        yaxis_title = "Chromosome name"
    else:
        scaffold_items = scaffold_for_plot["scaffold"]
        yaxis_title = "Scaffold name"
    
    # Apply intelligent sorting for better chromosome organization
    scaffold_items = _sort_chromosomes_intelligent(scaffold_items)
    
    # Reorder scaffold_for_plot data to match the new chromosome order
    if use_chrm:
        # Create mapping from old to new order
        old_chrm_list = scaffold_for_plot["chrm"]
        old_end_list = scaffold_for_plot["end"]
        
        # Create dictionary for quick lookup
        chrm_to_end = dict(zip(old_chrm_list, old_end_list))
        
        # Reorder end values according to new chromosome order
        scaffold_end_values = [chrm_to_end[chrm] for chrm in scaffold_items]
    else:
        # Create mapping from old to new order
        old_scaffold_list = scaffold_for_plot["scaffold"]
        old_end_list = scaffold_for_plot["end"]
        
        # Create dictionary for quick lookup
        scaffold_to_end = dict(zip(old_scaffold_list, old_end_list))
        
        # Reorder end values according to new scaffold order
        scaffold_end_values = [scaffold_to_end[scaffold] for scaffold in scaffold_items]
    
    # Calculate dynamic height based on number of scaffolds/chromosomes
    num_items = len(scaffold_items)
    
    # New sizing logic: minimum 50px per chromosome + 20px spacer
    # Use constants for chromosome dimensions
    chromosome_height = CHROMOSOME_HEIGHT
    vertical_spacer = VERTICAL_SPACER
    base_height = BASE_HEIGHT
    
    # Calculate total height needed
    dynamic_height = base_height + (num_items * (chromosome_height + vertical_spacer))
    
    # Set reasonable bounds
    dynamic_height = max(CANVAS_HEIGHT_MIN, min(CANVAS_HEIGHT_MAX, dynamic_height))
    
    # Calculate dynamic margins based on scaffold names length if available
    if len(scaffold_items) > 0:
        max_name_length = max(len(str(name)) for name in scaffold_items)
        # Base margin of 120px, then 10px per character, with maximum of 400px
        left_margin = max(120, min(400, max_name_length * 10))
    else:
        left_margin = 150
    
    # Calculate appropriate width - make it wider for better readability
    canvas_width = 1400  # default width, can be adjusted if needed
    
    logger.info(f"Drawing {num_items} {yaxis_title.lower()}s:")
    logger.info(f"  Canvas size: {canvas_width}x{dynamic_height}px")
    logger.info(f"  Left margin: {left_margin}px")
    logger.info(f"  Height per item: {chromosome_height + vertical_spacer}px")
    
    # Calculate dynamic font size based on number of items
    if num_items <= 20:
        font_size = 15
    elif num_items <= START_CUTOFF_MAX:
        font_size = 12
    else:
        font_size = 10

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=scaffold_end_values,
            y=scaffold_items,
            orientation="h",
            name="Scaffold",
            marker_color="#f3f4f7",
        )
    )
    fig.update_layout(barmode="overlay")
    fig.update_layout(
        title={
            "text": title_text,
            "y": 0.99,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        xaxis_title="bp",
        yaxis_title=yaxis_title,
    )

    fig.update_layout(
        xaxis=dict(
            automargin=True,
            showline=True,
            showgrid=False,
            showticklabels=True,
            linecolor="rgb(204, 204, 204)",
            linewidth=1,
            ticks="outside",
            rangemode="nonnegative",
            tickfont=dict(
                family="Arial",
                size=font_size,
                color="rgb(82, 82, 82)",
            ),
        ),
        # Turn off everything on y axis
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showline=False,
            showticklabels=True,
            ticklabelstep=1,
            tickwidth=15,
            tickfont=dict(
                family="Arial",
                size=font_size,
                color="rgb(82, 82, 82)",
            ),
        ),
        width=canvas_width,
        height=dynamic_height,
        margin=dict(
            autoexpand=True,
            l=left_margin,
            r=MARGIN_RIGHT // 2,  # Use half of standard right margin
            t=MARGIN_TOP + 10,    # Slightly more than standard
            b=MARGIN_BOTTOM,
        ),
        showlegend=True,
        plot_bgcolor="white",
    )

    fig.update_layout(legend=dict(font=dict(family="Arial", size=font_size, color="black")))

    fig.update_xaxes(range=[0, max(scaffold_end_values) + 1000])

    return fig, canvas_width, dynamic_height


def _create_matplotlib_karyotype(scaffold_for_plot, title_text, output_file, use_chrm=False, traces=None):
    """
    Create karyotype visualization using matplotlib for PNG/SVG export without external dependencies.

    Args:
        scaffold_for_plot: Scaffold data for plotting
        title_text: Plot title
        output_file: Output file path (will change extension to .png)
        use_chrm: Whether to use chromosome names
        traces: List of trace configurations (dict with plotly trace parameters)

    Returns:
        None
    """
    # Extract data
    if use_chrm:
        scaffold_items = scaffold_for_plot["chrm"]
        yaxis_title = "Chromosome name"
    else:
        scaffold_items = scaffold_for_plot["scaffold"]
        yaxis_title = "Scaffold name"

    # Apply same sorting as plotly version
    scaffold_items = _sort_chromosomes_intelligent(scaffold_items)

    # Reorder end values
    if use_chrm:
        chrm_to_end = dict(zip(scaffold_for_plot["chrm"], scaffold_for_plot["end"]))
        scaffold_end_values = [chrm_to_end[chrm] for chrm in scaffold_items]
    else:
        scaffold_to_end = dict(zip(scaffold_for_plot["scaffold"], scaffold_for_plot["end"]))
        scaffold_end_values = [scaffold_to_end[scaffold] for scaffold in scaffold_items]

    # Calculate figure dimensions
    num_items = len(scaffold_items)
    fig_height = max(6, min(30, 2 + num_items * 0.3))  # Height in inches
    fig_width = 14  # Width in inches

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Create y-axis positions
    y_positions = list(range(len(scaffold_items)))

    # Draw background bars (scaffolds)
    ax.barh(y_positions, scaffold_end_values, height=0.8,
            color='#f3f4f7', edgecolor='none', zorder=1)

    # Add traces if provided
    has_labeled_artists = False
    if traces:
        for trace_config in traces:
            # Extract matplotlib-compatible parameters from plotly trace config
            x_values = trace_config.get('x', [])
            base_values = trace_config.get('base', [])
            y_items = trace_config.get('y', [])
            color = trace_config.get('marker_color', 'black')
            label = trace_config.get('name', '')

            # Map y items to positions
            y_item_to_pos = {item: pos for pos, item in enumerate(scaffold_items)}
            trace_y_positions = [y_item_to_pos.get(item, -1) for item in y_items]

            # Filter out invalid positions
            valid_indices = [i for i, pos in enumerate(trace_y_positions) if pos >= 0]
            if valid_indices:
                valid_x = [x_values[i] for i in valid_indices]
                valid_base = [base_values[i] for i in valid_indices]
                valid_y_pos = [trace_y_positions[i] for i in valid_indices]

                ax.barh(valid_y_pos, valid_x, left=valid_base, height=0.8,
                       color=color, label=label, zorder=2)
                if label:  # Track if we added a labeled artist
                    has_labeled_artists = True

    # Formatting
    ax.set_yticks(y_positions)
    ax.set_yticklabels(scaffold_items)
    ax.set_xlabel('bp', fontsize=12)
    ax.set_ylabel(yaxis_title, fontsize=12)
    ax.set_title(title_text, fontsize=14, pad=20)

    # Style axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.grid(False)
    ax.set_axisbelow(True)

    # Set x-axis to start at 0
    ax.set_xlim(left=0)

    # Add legend only if we actually added labeled artists
    if has_labeled_artists:
        ax.legend(loc='best', frameon=False)

    # Save as PNG
    png_output = os.path.splitext(output_file)[0] + '.png'
    save_matplotlib_figure(fig, png_output, dpi=150)


def _create_and_save_bar_chart(scaffold_for_plot, title_suffix, output_file, use_chrm=False, traces=None):
    """
    Helper function to create a bar chart with given traces and save it.
    Saves both interactive HTML (plotly) and static PNG (matplotlib) versions.

    Args:
        scaffold_for_plot: Scaffold data for plotting
        title_suffix: Suffix to add to the title
        output_file: Output file path
        use_chrm: Whether to use chromosome names
        traces: List of trace configurations (dict with trace parameters)

    Returns:
        None
    """
    # Create plotly version for interactive HTML
    fig, canvas_width, canvas_height = _draw_chromosomes(scaffold_for_plot, title_suffix, use_chrm=use_chrm)

    if traces:
        for trace_config in traces:
            fig.add_trace(go.Bar(**trace_config))

    safe_write_figure(fig, output_file, width=canvas_width, height=canvas_height)

    # Create matplotlib version for static PNG/SVG export
    try:
        _create_matplotlib_karyotype(scaffold_for_plot, title_suffix, output_file, use_chrm, traces)
    except Exception as e:
        logger.warning(f"Failed to create matplotlib version: {e}")
        logger.warning("Continuing with plotly HTML only")


def _prepare_tr_traces(df_trs, names_filter=None, enhance=None):
    """
    Helper function to prepare TR traces for matplotlib rendering.

    Args:
        df_trs: List of dicts with tandem repeat data
        names_filter: Optional function to filter family names (e.g., lambda x: x != "SING")
        enhance: Optional minimum size to enhance small repeats

    Returns:
        List of trace configurations (dicts) for matplotlib
    """
    # Get unique family names from list of dicts
    names = set(record.get("family_name") for record in df_trs if record.get("family_name"))

    traces = []
    for name in names:
        if names_filter and not names_filter(name):
            continue

        # Filter items by family_name
        items = [record for record in df_trs if record.get("family_name") == name]

        # Extract data
        starts = [item.get("start") for item in items]
        lengths = [max(item.get("length", 0), enhance) for item in items] if enhance else [item.get("length", 0) for item in items]
        chrms = [item.get("chrm") for item in items]

        traces.append({
            'base': starts,
            'x': lengths,
            'y': chrms,
            'name': name,
            'marker_color': None,  # Let matplotlib choose colors
        })

    return traces


def _add_tr_families_by_name(fig, df_trs, names_filter=None, enhance=None):
    """
    Helper function to add tandem repeat families as traces to a figure.

    Args:
        fig: Plotly figure object
        df_trs: DataFrame with tandem repeat data
        names_filter: Optional function to filter family names (e.g., lambda x: x != "SING")
        enhance: Optional minimum size to enhance small repeats

    Returns:
        None (modifies fig in place)
    """
    # Get unique family names from list of dicts
    names = set(record.get("family_name") for record in df_trs if record.get("family_name"))

    for name in names:
        if names_filter and not names_filter(name):
            continue

        # Filter items by family_name
        items = [record for record in df_trs if record.get("family_name") == name]

        # Apply enhancement if specified
        if enhance:
            for item in items:
                item["length"] = max(item.get("length", 0), enhance)

        # Extract data for plotly
        starts = [item.get("start") for item in items]
        lengths = [item.get("length") for item in items]
        chrms = [item.get("chrm") for item in items]

        fig.add_trace(
            go.Bar(
                base=starts,
                x=lengths,
                y=chrms,
                orientation="h",
                name=name,
            )
        )


def _create_tr_visualization(scaffold_for_plot, title_text, df_trs, output_suffix,
                            use_chrm=False, names_filter=None, enhance=None):
    """
    Helper function to create and save a TR visualization with specified parameters.
    Saves both interactive HTML (plotly) and static PNG (matplotlib) versions.

    Args:
        scaffold_for_plot: Scaffold data for plotting
        title_text: Base title text
        df_trs: DataFrame with tandem repeat data
        output_suffix: Suffix for output file
        use_chrm: Whether to use chromosome names
        names_filter: Optional function to filter family names
        enhance: Optional minimum size to enhance small repeats

    Returns:
        Output file path
    """
    # Create plotly version for interactive HTML
    fig, canvas_width, canvas_height = _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=use_chrm)
    _add_tr_families_by_name(fig, df_trs, names_filter=names_filter, enhance=enhance)
    safe_write_figure(fig, output_suffix, width=canvas_width, height=canvas_height)

    # Create matplotlib version for static PNG export
    try:
        # Prepare traces for matplotlib
        traces = _prepare_tr_traces(df_trs, names_filter, enhance)
        _create_matplotlib_karyotype(scaffold_for_plot, title_text, output_suffix, use_chrm, traces)
    except Exception as e:
        logger.warning(f"Failed to create matplotlib version: {e}")
        logger.warning("Continuing with plotly HTML only")

    return output_suffix


def _draw_repeats_with_gaps(scaffold_for_plot, title_text, output_file, repeats_with_gap, size, use_chrm):
    """
    Helper function to draw repeats with gaps visualization.
    Saves both interactive HTML (plotly) and static PNG (matplotlib) versions.

    Args:
        scaffold_for_plot: Scaffold data for plotting
        title_text: Title text for the plot
        output_file: Output file path
        repeats_with_gap: List of repeats with gaps
        size: Enhancement size
        use_chrm: Whether to use chromosome names
    """
    # Create plotly version for interactive HTML
    fig, canvas_width, canvas_height = _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=use_chrm)

    # Convert repeats_with_gap (list of tuples/lists) to list of dicts
    repeats_with_gap_records = [
        {
            "chrm": row[0],
            "start": row[1],
            "end": row[2],
            "family_name": row[3],
            "gap_type": row[4],
            "length": row[5]
        }
        for row in repeats_with_gap
    ]

    # Filter by gap type and size
    aN_records = [r for r in repeats_with_gap_records if r["gap_type"] == "aN" and r["length"] < size]
    Na_records = [r for r in repeats_with_gap_records if r["gap_type"] == "Na" and r["length"] < size]
    aNa_records = [r for r in repeats_with_gap_records if r["gap_type"] == "aNa" and r["length"] < size]

    # Prepare traces for both plotly and matplotlib
    traces = []

    # Add aN traces
    if len(aN_records) > 0:
        aN_starts = [r["start"] for r in aN_records]
        aN_chrms = [r["chrm"] for r in aN_records]
        fig.add_trace(go.Bar(
            base=aN_starts, x=[size] * len(aN_records), y=aN_chrms,
            orientation="h", name="Tandem Repeat_with gap aN", marker_color="#FF00FF"
        ))
        fig.add_trace(go.Bar(
            base=[s + size for s in aN_starts], x=[size] * len(aN_records), y=aN_chrms,
            orientation="h", name="gaps aN", marker_color="#663399"
        ))
        traces.extend([
            {'base': aN_starts, 'x': [size] * len(aN_records), 'y': aN_chrms, 'name': "TR gap aN", 'marker_color': "#FF00FF"},
            {'base': [s + size for s in aN_starts], 'x': [size] * len(aN_records), 'y': aN_chrms, 'name': "gaps aN", 'marker_color': "#663399"}
        ])

    # Add Na traces
    if len(Na_records) > 0:
        Na_starts = [r["start"] for r in Na_records]
        Na_chrms = [r["chrm"] for r in Na_records]
        fig.add_trace(go.Bar(
            base=Na_starts, x=[size] * len(Na_records), y=Na_chrms,
            orientation="h", name="Tandem Repeat_with gap Na", marker_color="#00CED1"
        ))
        fig.add_trace(go.Bar(
            base=[s + size for s in Na_starts], x=[size] * len(Na_records), y=Na_chrms,
            orientation="h", name="gaps Na", marker_color="#00BFFF"
        ))
        traces.extend([
            {'base': Na_starts, 'x': [size] * len(Na_records), 'y': Na_chrms, 'name': "TR gap Na", 'marker_color': "#00CED1"},
            {'base': [s + size for s in Na_starts], 'x': [size] * len(Na_records), 'y': Na_chrms, 'name': "gaps Na", 'marker_color': "#00BFFF"}
        ])

    # Add aNa traces
    if len(aNa_records) > 0:
        third_size = size * 2 / 3
        aNa_starts = [r["start"] for r in aNa_records]
        aNa_chrms = [r["chrm"] for r in aNa_records]
        fig.add_trace(go.Bar(
            base=aNa_starts, x=[third_size] * len(aNa_records), y=aNa_chrms,
            orientation="h", name="Tandem Repeat_with gap aNa", marker_color="#00FF7F"
        ))
        fig.add_trace(go.Bar(
            base=[s + third_size for s in aNa_starts], x=[third_size] * len(aNa_records), y=aNa_chrms,
            orientation="h", name="gaps aNa", marker_color="#228B22"
        ))
        fig.add_trace(go.Bar(
            base=[s + third_size * 2 for s in aNa_starts], x=[third_size] * len(aNa_records), y=aNa_chrms,
            orientation="h", marker_color="#00FF7F"
        ))
        traces.extend([
            {'base': aNa_starts, 'x': [third_size] * len(aNa_records), 'y': aNa_chrms, 'name': "TR gap aNa", 'marker_color': "#00FF7F"},
            {'base': [s + third_size for s in aNa_starts], 'x': [third_size] * len(aNa_records), 'y': aNa_chrms, 'name': "gaps aNa", 'marker_color': "#228B22"},
            {'base': [s + third_size * 2 for s in aNa_starts], 'x': [third_size] * len(aNa_records), 'y': aNa_chrms, 'name': "", 'marker_color': "#00FF7F"}
        ])

    safe_write_figure(fig, output_file, width=canvas_width, height=canvas_height)

    # Create matplotlib version for static PNG export
    try:
        _create_matplotlib_karyotype(scaffold_for_plot, title_text, output_file, use_chrm, traces)
    except Exception as e:
        logger.warning(f"Failed to create matplotlib version: {e}")
        logger.warning("Continuing with plotly HTML only")


def _draw_repeats_without_gaps(scaffold_for_plot, title_text, output_file, repeats_without_gaps, size, use_chrm):
    """
    Helper function to draw repeats without gaps visualization.
    Saves both interactive HTML (plotly) and static PNG (matplotlib) versions.

    Args:
        scaffold_for_plot: Scaffold data for plotting
        title_text: Title text for the plot
        output_file: Output file path
        repeats_without_gaps: List of repeats without gaps
        size: Enhancement size
        use_chrm: Whether to use chromosome names
    """
    # Create plotly version for interactive HTML
    fig, canvas_width, canvas_height = _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=use_chrm)

    # Prepare traces for both plotly and matplotlib
    traces = []
    names = set([x["family_name"] for x in repeats_without_gaps])
    for name in names:
        items = [x for x in repeats_without_gaps if x["family_name"] == name]
        starts = [x.get("start") for x in items]
        lengths = [max(x.get("length", 0), size) for x in items]
        chrms = [x.get("chrm") for x in items]

        fig.add_trace(
            go.Bar(
                base=starts,
                x=lengths,
                y=chrms,
                orientation="h",
                name=name,
            )
        )

        traces.append({
            'base': starts,
            'x': lengths,
            'y': chrms,
            'name': name,
            'marker_color': None,  # Let matplotlib choose colors
        })

    safe_write_figure(fig, output_file, width=canvas_width, height=canvas_height)

    # Create matplotlib version for static PNG export
    try:
        _create_matplotlib_karyotype(scaffold_for_plot, title_text, output_file, use_chrm, traces)
    except Exception as e:
        logger.warning(f"Failed to create matplotlib version: {e}")
        logger.warning("Continuing with plotly HTML only")


def draw_karyotypes(
    output_file_name_prefix,
    title_text,
    df_trs,
    scaffold_for_plot,
    gaps_df,
    repeats_with_gap,
    repeats_without_gaps,
    use_chrm=False,
    enhance=ENHANCE_LARGE,
    gap_cutoff=GAP_CUTOFF_DEFAULT,
):
    """
    Draw various karyotype visualizations with tandem repeats and gaps.

    This function creates 10 different visualizations showing gaps and tandem repeats
    in various configurations (raw, enhanced, with/without singletons, etc.).
    """

    # Filter df_trs by scaffolds in scaffold_for_plot
    if use_chrm:
        allowed_scaffolds = set(scaffold_for_plot["chrm"])
        _df_trs = [record for record in df_trs if record.get("chrm") in allowed_scaffolds]
    else:
        allowed_scaffolds = set(scaffold_for_plot["scaffold"])
        _df_trs = [record for record in df_trs if record.get("chrm") in allowed_scaffolds]

    ### 1. Raw gaps
    _create_and_save_bar_chart(
        scaffold_for_plot,
        title_text + "(raw gaps)",
        output_file_name_prefix + ".gaps.svg",
        use_chrm=use_chrm,
        traces=[{
            "base": gaps_df["start"],
            "x": gaps_df["length"],
            "y": gaps_df["scaffold"],
            "orientation": "h",
            "name": "gaps",
            "marker_color": "rgba(0, 0, 0)",
        }]
    )

    ### 2. Enhanced gaps
    # Filter gaps by cutoff and enhance length
    filtered_indices = [i for i, length in enumerate(gaps_df["length"]) if length > gap_cutoff]
    _gaps_df = {
        "scaffold": [gaps_df["scaffold"][i] for i in filtered_indices],
        "start": [gaps_df["start"][i] for i in filtered_indices],
        "end": [gaps_df["end"][i] for i in filtered_indices],
        "length": [max(gaps_df["length"][i], enhance) for i in filtered_indices]
    }

    _create_and_save_bar_chart(
        scaffold_for_plot,
        title_text + "(enlarged gaps)",
        output_file_name_prefix + f".gaps.{gap_cutoff}bp.enhanced.svg",
        use_chrm=use_chrm,
        traces=[{
            "base": _gaps_df["start"],
            "x": _gaps_df["length"],
            "y": _gaps_df["scaffold"],
            "orientation": "h",
            "name": "gaps",
            "marker_color": "rgba(0, 0, 0)",
        }]
    )

    ### 3. Enhanced repeats_with_gap
    _draw_repeats_with_gaps(
        scaffold_for_plot,
        title_text + "(enlarged TRs with gaps)",
        output_file_name_prefix + ".repeats.with.gaps.enhanced.svg",
        repeats_with_gap,
        enhance,
        use_chrm
    )

    ### 4. Enhanced TRs without gaps
    _draw_repeats_without_gaps(
        scaffold_for_plot,
        title_text + " (TRs without gaps)",
        output_file_name_prefix + ".repeats.nogaps.enhanced.svg",
        repeats_without_gaps,
        enhance,
        use_chrm
    )

    ### 5-10. Various TR visualizations
    visualizations = [
        # (title_suffix, output_suffix, names_filter, enhance_value)
        (" (all)", ".raw.svg", None, None),
        (" (no singletons)", ".nosing.svg", lambda x: x != "SING", None),
        (" (only singletons)", ".sing.svg", lambda x: x == "SING", None),
        (" (enlarged, all)", ".raw.enhanced.svg", None, enhance),
        (" (enlarged, no singletons)", ".nosing.enchanced.svg", lambda x: x != "SING", enhance),
        (" (enlarged, only singletons)", ".sing.enchanced.svg", lambda x: x == "SING", enhance),
    ]

    for title_suffix, output_suffix, names_filter, enhance_value in visualizations:
        _create_tr_visualization(
            scaffold_for_plot,
            title_text + title_suffix,
            _df_trs,
            output_file_name_prefix + output_suffix,
            use_chrm=use_chrm,
            names_filter=names_filter,
            enhance=enhance_value
        )


def _sort_chromosomes_intelligent(scaffold_items):
    """
    Intelligent sorting of chromosomes/scaffolds:
    1. If diploid pattern (chr1_pat, chr1_mat) detected - group by chromosome number
    2. If simple pattern (chr1, chr2, chrX, chrZ) detected - sort by number
    3. Otherwise - keep original order (size-based)
    """
    import re
    
    # Check for diploid pattern: chr1_pat, chr1_mat, chr16_pat, chr16_mat, etc.
    diploid_pattern = re.compile(r'^chr(\d+|[XYZW])_([pm]at)$', re.IGNORECASE)
    simple_pattern = re.compile(r'^chr(\d+|[XYZW])$', re.IGNORECASE)
    
    diploid_matches = 0
    simple_matches = 0
    
    for item in scaffold_items:
        item_str = str(item)
        if diploid_pattern.match(item_str):
            diploid_matches += 1
        elif simple_pattern.match(item_str):
            simple_matches += 1
    
    # If majority are diploid pattern
    if diploid_matches > len(scaffold_items) * 0.6:
        logger.info("  📋 Detected diploid chromosome pattern - grouping maternal and paternal")
        return _sort_diploid_chromosomes(scaffold_items)
    
    # If majority are simple chr pattern
    elif simple_matches > len(scaffold_items) * 0.6:
        logger.info("  📋 Detected simple chromosome pattern - sorting by number")
        return _sort_simple_chromosomes(scaffold_items)
    
    else:
        logger.info("  📋 Using size-based chromosome order")
        return scaffold_items

def _sort_diploid_chromosomes(scaffold_items):
    """Sort diploid chromosomes: chr1_mat, chr1_pat, chr2_mat, chr2_pat, etc."""
    import re
    
    diploid_pattern = re.compile(r'^chr(\d+|[XYZW])_([pm]at)$', re.IGNORECASE)
    
    def get_sort_key(item):
        match = diploid_pattern.match(str(item))
        if match:
            chr_num, parent = match.groups()
            # Convert chromosome number to integer for proper sorting
            if chr_num.isdigit():
                chr_sort = int(chr_num)
            else:
                # Sex chromosomes come after autosomes
                sex_order = {'X': 1000, 'Y': 1001, 'Z': 1002, 'W': 1003}
                chr_sort = sex_order.get(chr_num.upper(), 2000)
            
            # Maternal first, then paternal
            parent_sort = 0 if parent.lower() == 'mat' else 1
            return (chr_sort, parent_sort)
        else:
            # Non-matching items go to the end
            return (TR_CUTOFF_LARGE, 0)
    
    return sorted(scaffold_items, key=get_sort_key)

def _sort_simple_chromosomes(scaffold_items):
    """Sort simple chromosomes: chr1, chr2, ..., chrX, chrY, etc."""
    import re
    
    simple_pattern = re.compile(r'^chr(\d+|[XYZW])$', re.IGNORECASE)
    
    def get_sort_key(item):
        match = simple_pattern.match(str(item))
        if match:
            chr_num = match.group(1)
            if chr_num.isdigit():
                return int(chr_num)
            else:
                # Sex chromosomes come after autosomes
                sex_order = {'X': 1000, 'Y': 1001, 'Z': 1002, 'W': 1003}
                return sex_order.get(chr_num.upper(), 2000)
        else:
            # Non-matching items go to the end
            return TR_CUTOFF_LARGE
    
    return sorted(scaffold_items, key=get_sort_key)

def draw_all(
    trf_file,
    fasta_file,
    chm2name,
    output_folder,
    taxon,
    genome_size,
    lenght_cutoff=10000000,
    enhance=ENHANCE_DEFAULT,
    gap_cutoff=1000,
    force_rerun=False,
):

    logger.info("Loading chromosomes...")
    scaffold_df = scaffold_length_sort_length(fasta_file, lenght_cutoff=lenght_cutoff)

    logger.info("Loading trs...")
    df_trs = read_trf_file(trf_file)
    # Filter by period > 5
    df_trs = [record for record in df_trs if float(record.get("period", 0)) > 5]
    logger.info(f"Quantity of TRs: {len(df_trs)}")

    if len(df_trs) > SAMPLE_SIZE_FOR_CLUSTERING:
        logger.warning("Too many TRs")
        logger.info("Filtering them...")
        # Sort by length descending and take top 2000
        df_trs = sorted(df_trs, key=lambda x: float(x.get("length", 0)), reverse=True)[:2000]
        logger.info(f"Updated quantity of TRs: 2000")

    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    ### TODO: save distances and tr2vector

    # Deprecated drawing functions removed - these visualizations are no longer generated
    # The data structures are preserved for potential future use
    name2monomers = {}
    name2lid = {}
    name2ids = {}
    id2names = {}

    # Draw gaps visualization

    # Extract project name from TRF file (e.g., GCF_000005845.2_ASM584v2_genomic.1kb.trf -> GCF_000005845.2_ASM584v2_genomic)
    trf_basename = os.path.basename(trf_file)
    # Remove .trf extension and any suffix like .1kb, .3kb, .10kb
    project_name = trf_basename.replace('.sat', '')
    for suffix in ['.1kb', '.3kb', '.10kb', '.micro', '.complex', '.pmicro', '.tssr']:
        project_name = project_name.replace(suffix, '')

    # Use BED file as the only storage format (no .pkl cache)
    bed_output_file = os.path.join(os.path.dirname(output_folder), f"{project_name}.gaps.bed")

    # Check if BED file exists and load from it, or compute and create it
    if os.path.isfile(bed_output_file) and os.path.getsize(bed_output_file) > 0 and not force_rerun:
        logger.info(f"Loading gaps data from existing BED file: {bed_output_file}")
        gaps_data = []
        try:
            with open(bed_output_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        chrm, start, end, name, length = parts[:5]
                        gaps_data.append((chrm, int(start), int(end), int(length)))
            logger.info(f"✓ Loaded {len(gaps_data)} gaps from BED file")
        except Exception as e:
            logger.warning(f"Failed to load gaps from BED file: {e}")
            logger.info("Computing gaps annotation...")
            gaps_data = get_gaps_annotation(fasta_file, genome_size, lenght_cutoff=lenght_cutoff)
    else:
        if force_rerun and os.path.isfile(bed_output_file):
            logger.info("Force rerun: Computing gaps annotation...")
        else:
            logger.info("Computing gaps annotation (this may take a while)...")
        gaps_data = get_gaps_annotation(fasta_file, genome_size, lenght_cutoff=lenght_cutoff)

    # Export/update gaps to BED format in output root directory (not in images/)
    # Only write if we computed new data or force_rerun
    if not os.path.isfile(bed_output_file) or force_rerun or len(gaps_data) == 0:
        logger.info(f"Exporting gaps to BED format: {bed_output_file}")
        try:
            with open(bed_output_file, 'w') as f:
                # BED format header
                f.write("# Gaps annotation from Satellome\n")
                f.write(f"# Project: {project_name}\n")
                f.write(f"# Taxon: {taxon}\n")
                f.write(f"# Total gaps: {len(gaps_data)}\n")
                f.write(f"# Format: chr\\tstart\\tend\\tname\\tscore\\tstrand\n")

                # Write gaps in BED6 format
                for chrm, start, end, length in gaps_data:
                    # BED format: chr, start, end, name, score, strand
                    # score = length of gap (for filtering)
                    # strand = . (not applicable for gaps)
                    f.write(f"{chrm}\t{start}\t{end}\tgap\t{length}\t.\n")

            logger.info(f"✓ Exported {len(gaps_data)} gaps to BED format")
        except Exception as e:
            logger.warning(f"Failed to export gaps to BED format: {e}")

    gaps_lengths = Counter([x[-1] for x in gaps_data])

    if gaps_lengths:
        logger.info("\n📊 Gaps Distribution Summary:")
        logger.info(SEPARATOR_LINE)
        total_gaps = sum(gaps_lengths.values())
        total_gap_length = sum(length * count for length, count in gaps_lengths.items())
        genome_coverage = (total_gap_length / genome_size) * 100
        
        logger.info(f"Total gaps found: {total_gaps:,}")
        logger.info(f"Total gap length: {total_gap_length:,} bp ({genome_coverage:.2f}% of genome)")
        logger.info(f"Average gap size: {total_gap_length // total_gaps:,} bp")
                
        logger.info("\nGap size distribution:")
        size_ranges = [
            (1, 100, "Very small (1-100 bp)"),
            (101, 1000, "Small (101-1,000 bp)"),
            (1001, TR_CUTOFF_LARGE, "Medium (1-10 kb)"),
            (10001, 100000, "Large (10-100 kb)"),
            (100001, float('inf'), "Very large (>100 kb)")
        ]
        
        for min_size, max_size, label in size_ranges:
            count = sum(gaps_lengths[size] for size in gaps_lengths 
                       if min_size <= size <= max_size)
            if count > 0:
                length_in_range = sum(size * gaps_lengths[size] for size in gaps_lengths 
                                    if min_size <= size <= max_size)
                percentage = (count / total_gaps) * 100
                coverage_percentage = (length_in_range / genome_size) * 100
                logger.info(f"  {label:<25}: {count:>6,} gaps ({percentage:>5.1f}%) - {length_in_range:>10,} bp ({coverage_percentage:>5.2f}% genome)")
        
        logger.info(SEPARATOR_LINE)
    else:
        logger.info("No gaps found in the genome!")

    # Convert gaps_data (list of tuples) to dict format
    gaps_df = {
        "scaffold": [row[0] for row in gaps_data],
        "start": [row[1] for row in gaps_data],
        "end": [row[2] for row in gaps_data],
        "length": [row[3] for row in gaps_data]
    }

    # Add chrm field to df_trs records (list of dicts)
    for record in df_trs:
        record["chrm"] = record.get("scaffold")

    chrms = set([x[0] for x in gaps_data])
    if chrms:
        logger.info(chrms)

    chrm2gapIT = {}
    for chrm in chrms:
        chrm2gapIT[chrm] = IntervalTree()

    for chrm, start, end, length in gaps_data:
        chrm2gapIT[chrm].addi(start, end, "gap")

    repeats_with_gap = []
    repeats_without_gaps = []
    for d in df_trs:
        chrm = d.get("chrm")
        start = d.get("start")
        end = d.get("end")
        family_name = d.get("family_name")

        if chrm not in chrm2gapIT:
            continue
        found_gaps = chrm2gapIT[chrm][start - GAP_SEARCH_WINDOW : end + GAP_SEARCH_WINDOW]
        if not found_gaps:
            repeats_without_gaps.append(d)
            continue
        found_gaps = list(found_gaps)
        start_gap = min([x.begin for x in found_gaps])
        end_gap = max([x.end for x in found_gaps])

        if start < start_gap and end < end_gap:
            gap_type = "aN"
        elif start < start_gap and end > end_gap:
            gap_type = "aNa"
        elif start_gap < start:
            gap_type = "Na"
        else:
            logger.debug("Unknown")
            logger.debug(f"{family_name}, {start}, {end}, {start_gap}, {end_gap}")
        repeats_with_gap.append(
            [
                chrm,
                min(start, start_gap),
                max(end, end_gap),
                family_name,
                gap_type,
                abs(min(start, start_gap) - max(end, end_gap)),
            ]
        )

    ### TODO: save gaps

    # Draw karyotypes    
    output_file_name_prefix = os.path.join(output_folder, f"{taxon}.karyo")
    draw_karyotypes(
        output_file_name_prefix,
        taxon,
        df_trs,
        scaffold_df,
        gaps_df,
        repeats_with_gap,
        repeats_without_gaps,
        use_chrm=chm2name,
        enhance=enhance,
        gap_cutoff=gap_cutoff,
    )
