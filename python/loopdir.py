import glob
import os

filelist = glob.glob('[u,v,w].0*00')
filelist += glob.glob('th.0*00')

for i in filelist:
    ntime = int( i.split('.')[-1] )
    if (ntime % 7200 != 0):
        print('Deleting: {0}'.format(i))
        #os.remove(i)

