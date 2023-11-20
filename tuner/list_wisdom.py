import pandas as pd
import numpy as np
import glob
import json

files = glob.glob('../wisdom/*float.wisdom')
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

kernels = np.unique(df['kernel'])
for kernel in kernels:
    print('--------------------')
    print(df[df['kernel'] == kernel])
