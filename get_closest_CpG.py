import os
import csv
test = 'test.bed'
with open('all_single_2.sorted', 'w') as out:
    writer = csv.writer(out, delimiter='\t')
    for file in os.listdir('/data/khandekara2/imputation'):
        if file.endswith('.single.sorted'):
            with open(file, 'r') as f:
                reader = f.read().splitlines()
                for i, row in enumerate(reader): #first 2 rows are a special case
                    ratio = -1
                    distance = -1
                    new = []
                    if i <= 1:
                        ratio = reader[i+2].split('\t')[4]
                        distance = int(reader[i+2].split('\t')[1]) - int(reader[i].split('\t')[1])
                    elif i >= len(reader) - 2: #last 2 rows are a special case
                        ratio = reader[i-2].split('\t')[4]
                        distance = int(reader[i].split('\t')[1]) - int(reader[i - 2].split('\t')[1])
                    else:
                        prev_prev_distance = abs(int(reader[i].split('\t')[1]) - int(reader[i-2].split('\t')[1]))
                        prev_distance = abs(int(reader[i].split('\t')[1]) - int(reader[i-1].split('\t')[1]))
                        next_distance = abs(int(reader[i+1].split('\t')[1]) - int(reader[i].split('\t')[1]))
                        next_next_distance = abs(int(reader[i+2].split('\t')[1]) - int(reader[i].split('\t')[1]))
                        distance = min(prev_distance, next_distance)
                        if distance != 1 and distance == prev_distance:
                            ratio = reader[i-1].split('\t')[4]
                        elif distance != 1 and distance == next_distance:
                            ratio = reader[i+1].split('\t')[4]
                        else: #we have the neighboring C's ratio(distance = 1)
                            if prev_distance == 1:
                                distance = min(next_distance, next_next_distance, prev_prev_distance)
                            if next_distance == 1:
                                distance = min(next_next_distance, prev_prev_distance, prev_distance)
                            if distance == next_distance:
                                ratio = reader[i+1].split('\t')[4]
                            if distance == prev_distance:
                                ratio = reader[i-1].split('\t')[4]
                            if distance == next_next_distance:
                                ratio = reader[i+2].split('\t')[4]
                            if distance == prev_prev_distance:
                                ratio = reader[i-2].split('\t')[4]

                    row = row.split('\t')
                    new.append(ratio)
                    new.append(distance)
                    output_row = row + new
                    writer.writerow(output_row)
