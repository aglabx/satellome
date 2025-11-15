#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.10.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
HTML report generation with embedded image support.

Provides utilities for creating self-contained HTML reports from visualization
images. Converts image files (PNG, SVG) to Base64-encoded data URIs and embeds
them directly in HTML, eliminating external file dependencies.

Functions:
    image_to_data_uri: Convert PNG image to Base64 data URI
    svg_to_data_uri: Convert SVG file to Base64 data URI
    create_html_report: Generate HTML report with embedded images from folder

Key Features:
    - Self-contained HTML output (no external image files needed)
    - Base64 encoding for portable reports
    - SVG and PNG image support
    - Automatic folder scanning for visualization files
    - Responsive HTML template with 60% width images

Typical Use Case:
    1. Generate visualization plots (SVG/PNG) in output folder
    2. Call create_html_report() to bundle all images into single HTML
    3. Share HTML file without worrying about missing image dependencies

Example:
    >>> # Generate report from visualization folder
    >>> from satellome.core_functions.tools.reports import create_html_report
    >>> create_html_report("output/plots/", "report.html")
    INFO:satellome.core_functions.tools.reports:HTML file with embedded image created successfully!
    INFO:satellome.core_functions.tools.reports:File: report.html
    >>>
    >>> # Convert single image to data URI
    >>> from satellome.core_functions.tools.reports import svg_to_data_uri
    >>> uri = svg_to_data_uri("plot.svg")
    >>> print(uri[:50])  # First 50 chars
    data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iODAw...

See Also:
    satellome.steps.trf_draw: Generates SVG visualizations for reports
"""

import base64
import logging
import os

logger = logging.getLogger(__name__)

def image_to_data_uri(image_path):
    """
    Convert PNG image file to Base64-encoded data URI.

    Reads binary image data and encodes as Base64 string with PNG MIME type.
    Result can be embedded directly in HTML <img> src attribute.

    Args:
        image_path (str): Path to PNG image file to convert

    Returns:
        str: Data URI string in format "data:image/png;base64,{encoded_data}"

    Example:
        >>> uri = image_to_data_uri("plot.png")
        >>> print(uri[:30])
        data:image/png;base64,iVBORw...
        >>>
        >>> # Use in HTML
        >>> html = f'<img src="{uri}" alt="Plot">'

    Note:
        - Assumes input is PNG format (hardcoded MIME type)
        - For SVG files, use svg_to_data_uri() instead
        - Result string can be large for high-resolution images
        - Binary file read in 'rb' mode for correct Base64 encoding
    """
    with open(image_path, "rb") as image_file:
        # Convert binary data to Base64 encoded string
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        return f"data:image/png;base64,{encoded_string}"

def svg_to_data_uri(svg_path):
    """
    Convert SVG vector graphics file to Base64-encoded data URI.

    Reads binary SVG data and encodes as Base64 string with SVG MIME type.
    Result can be embedded directly in HTML <img> src attribute while
    preserving vector graphics quality.

    Args:
        svg_path (str): Path to SVG file to convert

    Returns:
        str: Data URI string in format "data:image/svg+xml;base64,{encoded_data}"

    Example:
        >>> uri = svg_to_data_uri("chromosome_plot.svg")
        >>> print(uri[:35])
        data:image/svg+xml;base64,PHN2Zy...
        >>>
        >>> # Use in HTML
        >>> html = f'<img src="{uri}" alt="Chromosome plot" width="800">'

    Note:
        - Uses image/svg+xml MIME type (not image/svg)
        - SVG remains scalable even after Base64 encoding
        - Generally smaller file size than raster PNG for plots
        - Binary file read in 'rb' mode for correct Base64 encoding
        - Primary format for Satellome visualization outputs
    """
    with open(svg_path, "rb") as svg_file:
        # Convert binary data to Base64 encoded string
        encoded_string = base64.b64encode(svg_file.read()).decode('utf-8')
        return f"data:image/svg+xml;base64,{encoded_string}"

def create_html_report(image_folder, report_file):
    """
    Generate self-contained HTML report with embedded images from folder.

    Scans folder for SVG visualization files, converts each to Base64 data URI,
    and embeds all images in single HTML file. Creates portable report that
    can be shared without external file dependencies.

    Args:
        image_folder (str): Path to folder containing SVG/PNG image files
        report_file (str): Path to output HTML file to create

    Example:
        >>> # Generate report from output folder
        >>> create_html_report("output/chromosomes/", "analysis_report.html")
        INFO:satellome.core_functions.tools.reports:HTML file with embedded image created successfully!
        INFO:satellome.core_functions.tools.reports:File: analysis_report.html
        >>>
        >>> # Creates HTML with all SVG files embedded
        >>> # Can open report.html in browser without image folder

    Processing Steps:
        1. Scan image_folder for files ending with '.svg'
        2. Convert each SVG to Base64 data URI
        3. Generate <img> tags with 60% width styling
        4. Wrap images in minimal HTML5 template
        5. Write complete HTML to report_file
        6. Log success message with file path

    Generated HTML Structure:
        - HTML5 DOCTYPE with UTF-8 charset
        - Responsive viewport meta tag
        - All images at 60% width (stacked vertically)
        - No external dependencies (CSS, JS, images)

    Note:
        - Currently only scans for .svg files (line 30 filter)
        - SVG conversion via svg_to_data_uri(), PNG via image_to_data_uri()
        - Images appear in filesystem order (not sorted)
        - No error handling for missing folder or read failures
        - Logs completion via module logger (INFO level)
        - Report file overwritten if already exists (mode 'w')
    """

    image_files = [os.path.join(image_folder, f) for f in os.listdir(image_folder) if f.endswith('.svg')]

    images_content = ""

    for image in image_files:
        if image.endswith('.svg'):
            data_uri = svg_to_data_uri(image)
        else:
            data_uri = image_to_data_uri(image)
        im = f'<img src="{data_uri}" alt="{image}" width="60%">'
        images_content += im
    
    html_content = f"""<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Embedded Image</title>
    </head>
    <body>
        {images_content}
    </body>
    </html>
    """

    with open(report_file, 'w') as file:
        file.write(html_content)

    logger.info("HTML file with embedded image created successfully!")
    logger.info(f"File: {report_file}")
