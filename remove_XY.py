import os

for file in os.listdir():
    if file.endswith(''):
        with open(file, 'r') as csvin, open(file + '.noXY', 'w') as csvout:
