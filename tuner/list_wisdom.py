import pandas as pd
import numpy as np
import glob
import json
import argparse

parser = argparse.ArgumentParser(
    description='Read wisdom files and print them in a table.')
parser.add_argument('-u', '--unique_values', help='Give the unique values in a column')
parser.add_argument('-s', '--select', help='Select a subset of the data based on a dictionary', type=json.loads)

args = parser.parse_args()

files = glob.glob('../wisdom/*.wisdom')
files.sort()

df = pd.DataFrame(
        columns=['device', 'precision', 'problem_size', 'kernel', 'config'])

for f in files:
    name = f.split('/')[-1].split('.')[0]
    prec = name.split('_')[-1]
    kernel = '_'.join(name.split('_')[:-1])

    with open (f, 'r') as f:
        for l in f.readlines():
            if 'device_name' in l:
                data = json.loads(l)
                env = data['environment']

                serie = pd.Series(
                        {'device': env['device_name'],
                         'precision': prec,
                         'problem_size': data['problem_size'],
                         'kernel': kernel,
                         'config': data['config']})

                df = pd.concat([df, serie.to_frame().T], ignore_index=True)

if args.unique_values is not None:
    print(df[args.unique_values].unique())
    exit()

df.sort_values(by=['kernel', 'device'], inplace=True)

if args.select is not None:
    query_string = ' & '.join([f'{k} == "{v}"' for k, v in args.select.items()])
    df = df.query(query_string).drop(args.select.keys(), axis=1)

pd.options.display.max_columns = None
pd.options.display.max_rows = None
pd.options.display.width = None
print(df)
