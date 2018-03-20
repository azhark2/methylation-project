
import pandas as pd
import pickle
import os
import csv
#combine all mutation files into one file

files = [f for f in os.listdir('/data/khandekara2/mutation_data/processed_data') if os.path.isfile(f)]
with open('all_mutations_2.bed', 'w') as csvout:
    for f in files:
        if f.endswith('.noDuplicates.bed'):
            with open(f, 'r') as csvin:
                reader = csv.reader(csvin, delimiter='\t')
                writer = csv.writer(csvout, delimiter='\t')
                for row in reader:
                    writer.writerow(row)
