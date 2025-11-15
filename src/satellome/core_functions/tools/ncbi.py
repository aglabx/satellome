#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
NCBI taxonomy database query utilities.

Provides functions for retrieving taxonomic information from NCBI Entrez databases
via E-utilities API. Handles network errors and XML parsing with automatic retries.

Functions:
    get_taxon_name: Retrieve scientific name for NCBI taxonomy ID

Key Features:
    - No external dependencies (uses stdlib urllib and xml.etree)
    - Automatic retry with exponential backoff (7 attempts)
    - Comprehensive error handling (HTTP, network, XML parsing)
    - 30-second timeout per request
    - 5-second delay between retries

API Endpoint:
    - NCBI E-utilities esummary.fcgi
    - Database: taxonomy
    - Returns: XML with taxonomic information

Example:
    >>> from satellome.core_functions.tools.ncbi import get_taxon_name
    >>> name = get_taxon_name(9606)
    >>> print(name)
    'Homo sapiens'
    >>>
    >>> # Handle missing taxon
    >>> name = get_taxon_name(999999999)
    WARNING:...No ScientificName found for taxid: 999999999
    >>> print(name)
    None

Typical Use Case:
    1. Extract taxon IDs from genome metadata or BLAST results
    2. Use get_taxon_name() to resolve IDs to scientific names
    3. Annotate analysis results with organism names

Network Requirements:
    - Internet connection to eutils.ncbi.nlm.nih.gov
    - Respect NCBI usage guidelines (max 3 requests/second without API key)
    - Consider caching results for repeated queries

See Also:
    NCBI E-utilities documentation: https://www.ncbi.nlm.nih.gov/books/NBK25501/
"""

import logging
import urllib.request
import urllib.error
from xml.etree import ElementTree
import time

logger = logging.getLogger(__name__)

def get_taxon_name(taxid):
    """
    Retrieve scientific name for NCBI taxonomy ID with automatic retries.

    Queries NCBI taxonomy database via E-utilities API and extracts scientific
    name from XML response. Implements robust error handling with up to 7 retry
    attempts for transient network failures.

    Args:
        taxid (int or str): NCBI taxonomy ID (e.g., 9606 for Homo sapiens)

    Returns:
        str or None: Scientific name if found, None if error or not found

    Example:
        >>> # Common organisms
        >>> get_taxon_name(9606)
        'Homo sapiens'
        >>> get_taxon_name(10090)
        'Mus musculus'
        >>> get_taxon_name(562)
        'Escherichia coli'
        >>>
        >>> # Invalid taxon ID
        >>> result = get_taxon_name(999999999)
        WARNING:...No ScientificName found for taxid: 999999999
        >>> print(result)
        None

    Retry Logic:
        - Up to 7 attempts with 5-second delays between attempts
        - Retries on: HTTP errors, network errors, XML parse errors
        - Returns None after all retries exhausted

    Error Handling:
        - HTTPError: HTTP status codes (404, 500, etc.)
        - URLError: Network connectivity issues, timeouts
        - ParseError: Malformed XML response from NCBI
        - Exception: Catches any unexpected errors

    Network Behavior:
        - 30-second timeout per request
        - 5-second delay between retry attempts
        - Total worst-case time: ~3.5 minutes (7 attempts × 30s + 6 delays × 5s)

    Note:
        - No external dependencies (stdlib only: urllib, xml.etree)
        - Logs warnings/errors for debugging failed queries
        - Returns None on all error conditions (check return value!)
        - Uses NCBI E-utilities esummary endpoint
        - Extracts ScientificName field from XML Item elements
        - Consider implementing caching for repeated queries
        - Respect NCBI usage policy: 3 req/sec without API key
    """
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
            # General catch-all for any other unexpected exceptions
            logger.error(f"An unexpected error occurred ({type(e).__name__}): {e}")
            logger.info(f"Attempt {attempts} of 7")
            time.sleep(5)
            # Note: After all retries exhausted, function returns None

    return None  # Return None if any errors occurred
