#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com


import logging
import urllib.request
import urllib.error
from xml.etree import ElementTree
import time

logger = logging.getLogger(__name__)

def get_taxon_name(taxid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={taxid}"

    attempts = 0

    while attempts < 7:
        attempts += 1
        try:
            # Use urllib instead of requests - no external dependencies needed
            with urllib.request.urlopen(url, timeout=30) as response:
                content = response.read()

            tree = ElementTree.fromstring(content)

            for item in tree.findall(".//Item[@Name='ScientificName']"):
                return item.text

            # If no scientific name found, log warning
            logger.warning(f"No ScientificName found for taxid: {taxid}")
            return None

        except urllib.error.HTTPError as e:
            # HTTP errors (404, 500, etc.)
            logger.warning(f"HTTP error fetching data from NCBI: {e.code} {e.reason}")
            logger.info(f"Attempt {attempts} of 7")
            time.sleep(5)
        except urllib.error.URLError as e:
            # Network issues, invalid URL, timeout
            logger.warning(f"Network error fetching data from NCBI: {e.reason}")
            logger.info(f"Attempt {attempts} of 7")
            time.sleep(5)
        except ElementTree.ParseError:
            # Handles issues when parsing the XML
            logger.error("Error parsing the XML response from NCBI.")
            time.sleep(5)
        except Exception as e:
            # General catch-all for any other exceptions
            logger.error(f"An unexpected error occurred: {e}")
            logger.info(f"Attempt {attempts} of 7")
            time.sleep(5)

    return None  # Return None if any errors occurred
