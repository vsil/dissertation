"""
Downloads fasta files in .json file outputed by extract_accession_numbers.py from Genbank
"""

from Bio import Entrez, SeqIO
import json
import os


Entrez.email = 'vicente.silvestre@tecnico.ulisboa.pt'

def download_fasta_from_genbank(accession_number, output_filename):
    try:
        # Fetch the GenBank record
        handle = Entrez.efetch(db='nucleotide', id=accession_number, rettype='gb', retmode='text')
        record = SeqIO.read(handle, 'genbank')
        handle.close()

        # Write the sequence to a FASTA file
        SeqIO.write(record, output_filename, 'fasta')
        print(f"FASTA file '{output_filename}' downloaded successfully.")

    except Exception as e:
        print(f"An error occurred: {e}")


# path to JSON file outputed by extract_accession_numbers.py
json_file_path = 'accession_numbers.json'

# Load the JSON file into a Python dictionary
with open(json_file_path, 'r') as json_file:
    data_dict = json.load(json_file)


for key in data_dict:
    for j in range(0,len(data_dict[key])):
        path = os.path.join('genomes', key)
        
        # Ensure the directory exists
        os.makedirs(path, exist_ok=True)

        path = os.path.join(path, data_dict[key][j]+'.fasta')
        print(path)
        download_fasta_from_genbank(data_dict[key][j], path)
        



