#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.10.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import base64
import os

def image_to_data_uri(image_path):
    """Converts an image to a Base64-encoded data URI."""
    with open(image_path, "rb") as image_file:
        # Convert binary data to Base64 encoded string
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        return f"data:image/png;base64,{encoded_string}"

def svg_to_data_uri(svg_path):
    """Converts an SVG file to a Base64-encoded data URI."""
    with open(svg_path, "rb") as svg_file:
        # Convert binary data to Base64 encoded string
        encoded_string = base64.b64encode(svg_file.read()).decode('utf-8')
        return f"data:image/svg+xml;base64,{encoded_string}"

def create_html_report(image_folder, report_file):

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

    print("HTML file with embedded image created successfully!")
    print(f"File: {report_file}")
