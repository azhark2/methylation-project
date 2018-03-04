import csv
import sys

#remove duplicates in bed files

def remove(*args):
        duplicates = set([])
        count = 0
        if len(args) > 0:
                for arg in args:
                        with open(arg, 'r') as f, open(arg + '.noDuplicates', 'w') as csvout:
                                reader = csv.reader(f, delimiter='\t')
                                writer = csv.writer(csvout, delimiter='\t')
                                for row in reader:
                                        info = (row[2], row[3], row[4])
                                        if info not in duplicates:
                                                duplicates.add(info)
                                                writer.writerow(row)
                                        else:
                                                count += 1

                                print("There were %d duplicates found" % count)


if len(sys.argv) > 1:
        file = sys.argv[1]
        remove(file)


