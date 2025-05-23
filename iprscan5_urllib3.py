#!/usr/bin/env python3
# $Id: iprscan5_urllib3.py 2106 2012-05-01 17:00:40Z hpm $
# ======================================================================
# InterProScan 5 (REST) Python 3 client using urllib3 and 
# xmltramp2 (https://pypi.python.org/pypi/xmltramp2/).
#
# Tested with:
#  Python 3.4.3
#
import os
import sys
import time
import argparse
import requests
from Bio import SeqIO

# See:
# http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_rest
# http://www.ebi.ac.uk/Tools/webservices/tutorials/python
# ======================================================================
# Base URL for service
BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"

def read_sequence(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    return record.id, str(record.seq)

def submit_job(sequence, email, goterms=True, pathways=False):
    data = {
        "email": email,
        "sequence": sequence,
        "goterms": "true" if goterms else "false",
        "pathways": "true" if pathways else "false",
        "outputFormat": "xml"
    }
    response = requests.post(f"{BASE_URL}/run", data=data)
    response.raise_for_status()
    return response.text.strip()

def poll_job(job_id, interval=10, max_tries=100):
    for attempt in range(max_tries):
        status = requests.get(f"{BASE_URL}/status/{job_id}").text.strip()
        if status in ["FINISHED", "FAILURE", "ERROR"]:
            return status
        time.sleep(interval)
    return "TIMEOUT"

def download_result(job_id, output_file):
    url = f"{BASE_URL}/result/{job_id}/xml"
    response = requests.get(url)
    response.raise_for_status()
    with open(output_file, "w") as out:
        out.write(response.text)

def main():
    parser = argparse.ArgumentParser(description="Submit a protein FASTA to InterProScan REST API and download XML result.")
    parser.add_argument("--sequence", required=True, help="Path to protein FASTA file (single sequence).")
    parser.add_argument("--email", required=True, help="Your email (required by InterProScan API).")
    args = parser.parse_args()

    fasta_file = os.path.abspath(args.sequence)
    out_dir = os.path.dirname(fasta_file)

    try:
        print(f"Reading sequence from {fasta_file} ...")
        protein_id, sequence = read_sequence(fasta_file)
        output_file = os.path.join(out_dir, f"{protein_id}.xml")

        if os.path.exists(output_file):
            print(f"Output already exists: {output_file}. Skipping.")
            return

        print(f"Submitting job for protein ID: {protein_id}")
        job_id = submit_job(sequence, args.email)
        print(f"Job submitted. Job ID: {job_id}")

        print("Polling job status...")
        status = poll_job(job_id)
        if status != "FINISHED":
            print(f"Job failed or timed out: {status}")
            return

        print(f"Downloading result to: {output_file}")
        download_result(job_id, output_file)
        print("Done.")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
