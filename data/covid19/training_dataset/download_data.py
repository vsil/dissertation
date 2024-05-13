from Bio import Entrez, SeqIO
import os

def download_genomes_by_locus(locus_query, output_path):
    """
    Downloads FASTA genomes from GenBank nucleotide database based on the LOCUS field.

    Args:
    - locus_query (str): Query string for the LOCUS field.
    - output_file (str): File path to save the downloaded genomes in FASTA format.
    """
    Entrez.email = "vicente.silvestre@tecnico.ulisboa.pt"  # Enter your email address here
    search_term = f'"{locus_query}"[LOCUS]'
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10)
    record = Entrez.read(handle)
    handle.close()


    if "IdList" in record and len(record["IdList"]) > 0:
        gb_ids = ",".join(record["IdList"])
        handle = Entrez.efetch(db="nucleotide", id=gb_ids, rettype="fasta", retmode="text")
        
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        output_file = "{}.fasta".format(locus_query)
        output_file = os.path.join(output_path, output_file)

        with open(output_file, "w") as f:
            f.write(handle.read())
        handle.close()
        print("Genomes successfully downloaded and saved to", output_file)
    else:
        print("No genomes found matching the query: "+locus_query)


def read_data_from_file(file_path):
    """
    Read data from a .txt file and store it in a dictionary.

    Args:
    - file_path (str): Path to the .txt file.

    Returns:
    - data_dict (dict): Dictionary containing data read from the file.
    """
    data_dict = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        key = None
        for line in lines:
            if line.startswith('>'):
                key = line.strip('>\n')
                data_dict[key] = []
            else:
                data_dict[key] = line.strip().split(',')
    return data_dict



# Example usage
file_path = "data/source/test4.txt"  # Path to your .txt file
data_dict = read_data_from_file(file_path)
#print(data_dict.keys())
#print(data_dict['Adenoviridae'][0:10])


for key, value in data_dict.items():
    for id in value:
        output_path = "data/Test4/{}".format(key)

        file_name = id+'.fasta'
        file_path = os.path.join(output_path, file_name)
        if(os.path.exists(file_path)):
            print("already downloaded {}, next.".format(id))
            continue

        download_genomes_by_locus(id, output_path)



