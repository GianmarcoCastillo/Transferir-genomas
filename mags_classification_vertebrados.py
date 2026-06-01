#!/usr/bin/env python3
import csv
import time
import sys
import re
import requests
import xml.etree.ElementTree as ET
from typing import Optional, Tuple
from Bio import Entrez

Entrez.email = "lothor1106@gmail.com"

INPUT_FILE = "mags_filtrados_animales.tsv"
OUTPUT_ALL = "todos_con_clasificacion.tsv"
OUTPUT_INVERTEBRADOS = "solo_invertebrados.tsv"
LOG_FILE = "clasificacion.log"

INVERT_KEYWORDS = [
    "insect", "nematode", "mollusk", "crustacean", "worm", "arthropod",
    "annelid", "cnidarian", "echinoderm", "sponge", "invertebrate",
    "bee", "ant", "beetle", "butterfly", "moth", "mosquito", "fly",
    "spider", "tick", "mite", "crab", "shrimp", "lobster", "snail",
    "slug", "clam", "oyster", "octopus", "squid", "coral", "jellyfish",
    "starfish", "urchin", "drosophila", "caenorhabditis"
]

VERT_KEYWORDS = [
    "human", "mouse", "rat", "cattle", "pig", "chicken", "fish",
    "zebrafish", "vertebrate", "mammal", "bird", "reptile", "amphibian",
    "homo sapiens", "mus musculus", "rattus", "bos taurus"
]

def log_message(accession: str, status: str, msg: str):
    with open(LOG_FILE, "a") as f:
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"{timestamp},{accession},{status},{msg}\n")

def normalize(text: str) -> str:
    if not text:
        return ""
    return re.sub(r'[^\w\s]', '', text.lower().strip())

def classify_host(host_str: str) -> Tuple[str, str]:
    if not host_str or host_str == "unknown":
        return "unknown", "missing host information"
    norm = normalize(host_str)
    for kw in VERT_KEYWORDS:
        if kw in norm:
            return "vertebrate", f"keyword '{kw}' in host"
    for kw in INVERT_KEYWORDS:
        if kw in norm:
            return "invertebrate", f"keyword '{kw}' in host"
    if re.search(r'(idae|inae|oidea|iformes)$', norm):
        return "invertebrate", f"suffix 'idae/inae/oidea' in host"
    return "unknown", "no clear match"

def query_ncbi(biosample_id: str) -> Tuple[Optional[str], Optional[str]]:
    try:
        handle = Entrez.esearch(db="biosample", term=biosample_id, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            return None, f"ID {biosample_id} not found in NCBI"
        bs_id = record["IdList"][0]
        handle = Entrez.esummary(db="biosample", id=bs_id)
        summary = Entrez.read(handle)
        handle.close()
        if not summary or len(summary) == 0:
            return None, "Cannot read BioSample summary from NCBI"
        if hasattr(summary[0], "Attributes"):
            for attr in summary[0].Attributes:
                if attr.get("name") == "host" or attr.get("harmonized_name") == "host":
                    host_value = attr.get("value", "unknown")
                    if host_value:
                        return host_value, None
        return None, "NCBI BioSample record has no 'host' field"
    except Exception as e:
        return None, f"NCBI exception: {str(e) or 'unknown error'}"

def query_ena(biosample_id: str) -> Tuple[Optional[str], Optional[str]]:
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{biosample_id}"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        root = ET.fromstring(response.content)
        ns = {'sra': 'http://www.ebi.ac.uk/ena/browser/api/xml'}
        for sample_attr in root.findall('.//SAMPLE_ATTRIBUTE', ns):
            tag_elem = sample_attr.find('TAG', ns)
            if tag_elem is not None and tag_elem.text and tag_elem.text.lower() == 'host':
                value_elem = sample_attr.find('VALUE', ns)
                if value_elem is not None and value_elem.text:
                    return value_elem.text.strip(), None
        for sample_attr in root.findall('.//SAMPLE_ATTRIBUTE'):
            tag_elem = sample_attr.find('TAG')
            if tag_elem is not None and tag_elem.text and tag_elem.text.lower() == 'host':
                value_elem = sample_attr.find('VALUE')
                if value_elem is not None and value_elem.text:
                    return value_elem.text.strip(), None
        return None, "ENA record has no 'host' attribute"
    except requests.exceptions.RequestException as e:
        return None, f"ENA connection error: {str(e)}"
    except ET.ParseError as e:
        return None, f"ENA XML parse error: {str(e)}"
    except Exception as e:
        return None, f"ENA exception: {str(e)}"

def get_host_from_biosample(biosample_id: str) -> Tuple[Optional[str], Optional[str]]:
    host, err = query_ncbi(biosample_id)
    if host:
        return host, None
    host, err2 = query_ena(biosample_id)
    if host:
        return host, None
    if err2:
        return None, f"NCBI: {err or 'no data'} | ENA: {err2}"
    else:
        return None, err or "No host from NCBI or ENA"

def main():
    print("=== Metagenome classifier (NCBI + ENA) ===")
    print(f"Reading {INPUT_FILE}...")
    with open(LOG_FILE, "w") as f:
        f.write("timestamp,biosample,status,message\n")
    try:
        with open(INPUT_FILE, "r", encoding="utf-8") as infile:
            reader = csv.reader(infile, delimiter='\t')
            header = next(reader)
            required = ["assembly_accession", "bioproject", "biosample"]
            for col in required:
                if col not in header:
                    print(f"ERROR: Column '{col}' not found in TSV")
                    sys.exit(1)
            idx_bs = header.index("biosample")
            new_header = header + ["host_ncbi", "categoria", "razon"]
            todos = [new_header]
            invertebrados = [new_header]
            for i, row in enumerate(reader, start=1):
                if len(row) < len(header):
                    log_message("", "WARNING", f"Row {i}: {len(row)} columns")
                    continue
                biosample = row[idx_bs].strip()
                print(f"Processing {i} : {biosample}")
                host, error_msg = None, None
                if biosample and biosample.startswith(("SAMN", "SAME", "SAMD")):
                    host, error_msg = get_host_from_biosample(biosample)
                    time.sleep(0.34)
                else:
                    error_msg = "Invalid ID (not starting with SAMN/SAME/SAMD)"
                categoria = "unknown"
                razon = ""
                if host:
                    categoria, razon = classify_host(host)
                else:
                    razon = error_msg or "no host data"
                nueva_fila = row + [host or "desconocido", categoria, razon]
                todos.append(nueva_fila)
                if categoria == "invertebrate":
                    invertebrados.append(nueva_fila)
                log_message(biosample, categoria, razon)
    except FileNotFoundError:
        print(f"\nERROR: File '{INPUT_FILE}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}")
        log_message("SYSTEM", "CRITICAL", str(e))
        sys.exit(1)
    with open(OUTPUT_ALL, "w", encoding="utf-8", newline='') as f:
        csv.writer(f, delimiter='\t').writerows(todos)
    with open(OUTPUT_INVERTEBRADOS, "w", encoding="utf-8", newline='') as f:
        csv.writer(f, delimiter='\t').writerows(invertebrados)
    total = len(todos) - 1
    num_inv = len(invertebrados) - 1
    print("\n" + "="*50)
    print(f"Total processed: {total}")
    print(f"Invertebrates found: {num_inv}")
    print(f"Output files: {OUTPUT_ALL}, {OUTPUT_INVERTEBRADOS}, {LOG_FILE}")
    print("="*50)

if __name__ == "__main__":
    main()
