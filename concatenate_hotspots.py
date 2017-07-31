import os
import pandas as pd

frames = []
files = [f for f in os.listdir('.') if os.path.isfile(f)]
for f in files:
    if os.stat(f).st_size != 0:
        df = pd.read_csv(f, sep='\t')
        frames.append(df)

result = pd.concat(frames)

result.to_csv( 'all_hotspots.tsv', sep='\t')




