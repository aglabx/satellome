#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com


import requests
from xml.etree import ElementTree

def get_taxon_name(taxid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={taxid}"
    
    try:
        response = requests.get(url, timeout=30)
        
        response.raise_for_status()
        tree = ElementTree.fromstring(response.content)
        
        for item in tree.findall(".//Item[@Name='ScientificName']"):
            return item.text
        
        if taxon_name_element is None:
            raise ValueError(f"Invalid taxid: {taxid} or the response format has changed.")
        
        return taxon_name_element.text
    
    except requests.RequestException as e:
        # Handles any kind of request issues (e.g. network issues, invalid URL, etc.)
        print(f"Error fetching data from NCBI: {e}")
    except ElementTree.ParseError:
        # Handles issues when parsing the XML
        print("Error parsing the XML response from NCBI.")
    except Exception as e:
        # General catch-all for any other exceptions
        print(f"An unexpected error occurred: {e}")

    return None  # Return None if any errors occurred
