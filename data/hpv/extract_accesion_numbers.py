"""
Saves the accession numbers and genotypes to a .json file from original data file
"""


import csv
import pickle
from collections import Counter
import json

# Initialize empty lists for each column
column1 = []
column2 = []

# Specify the CSV file path
csv_file_path = 'data.csv'

j=0
data_dict = {}

# Open and read the CSV file
with open(csv_file_path) as csvfile:
    
    reader = csv.reader(csvfile)
    
    # Iterate through each row in the CSV
    for row in reader:
        if(j<3):
            j+=1
            continue
        if(row[2] not in data_dict.keys()):
            data_dict[row[2]] = []
        data_dict[row[2]].append(row[1])

"""         if row[8]!='-':
            data_dict[row[0]] = row[8]   #ampicillin    -> stores genome ID and ampicilin class in data_dict
            
 """


json.dump(data_dict, open("accession_numbers.json", 'w' ) )