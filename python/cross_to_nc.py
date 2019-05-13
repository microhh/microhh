import microhh_tools as mht     # available in microhh/python directory
import argparse

#Parse command line and namelist options
parser = argparse.ArgumentParser(description='Convert microHH binary cross-sections to netCDF4 files.')
parser.add_argument('-m', '--mode', help = 'mode of the cross section', required=True, choices = ['xy', 'xz', 'yz', 'surf'])
parser.add_argument('-f', '--filename', help='ini file name')
parser.add_argument('-v', '--vars', nargs='*', help='variable names')
parser.add_argument('-x', '--index', nargs='*', help='indices')
parser.add_argument('-i', '--iotimeprec', help='iotimeprec', type=int)
parser.add_argument('-e', '--endian', help='endianess', choices = ['little', 'big'],default='little')
parser.add_argument('-p', '--precision', help='precision', choices = ['single', 'double'],default='single')

args = parser.parse_args()

mode = args.mode
indexes = args.index

nl = mht.Read_namelist(args.filename)
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']
starttime  = nl['time'] ['starttime']
endtime    = nl['time'] ['endtime']
sampletime = nl['cross']['sampletime']

variables = args.vars if args.vars is not None else nl['cross']['crosslist']

if args.iotimeprec is None:
    try:
        iotimeprec = nl['time']['iotimeprec']
    except KeyError:
        iotimeprec = 0.
        pass
else:
    iotimeprec = args.iotimeprec
endian = args.endian
precision = args.precision
#End option parsing


# calculate the number of iterations
niter = int((endtime-starttime) / sampletime + 1)

grid = mht.Read_grid(itot, jtot, ktot, endian = endian, precision = precision)

# Loop over the different variables
for variable in variables:
    try:
        filename = "{0}.{1}.nc".format(variable,mode)
        if mode is not 'surf':
            if indexes == None:
                indexes_local = mht.get_cross_indices(variable, mode)
            else:
                indexes_local = indexes
        
        dim = {'time' : range(niter), 'z' : range(ktot), 'y' : range(jtot), 'x': range(itot)}
        if mode == 'xy':
            dim.update({'z' : indexes_local})
            n = itot * jtot
        elif mode == 'surf':
            dim.pop('z')
            n = itot * jtot
            indexes_local = [-1]
        elif mode == 'xz':
            dim.update({'y' : indexes_local})
            n = itot * ktot
        elif mode == 'yz':
            dim.update({'x' : indexes_local})
            n = ktot * jtot
            
        ncfile = mht.Create_ncfile(grid, filename, variable, dim)
        for t in range(niter):
            print(t)
            for k in range(len(indexes_local)):
                index = indexes_local[k]
                otime = int((starttime + t*sampletime) / 10**iotimeprec)
                f_in  = "{0:}.{1}.{2:05d}.{3:07d}".format(variable, mode, index, otime)
        
                try:
                    fin = mht.Read_binary(f_in, endian=endian, precision=precision)
                except Exception as ex:
                    print (ex)
                    raise Exception('Stopping: cannot find file {}'.format(f_in))
        
                print("Processing %8s, time=%7i, index=%4i"%(variable, otime, index))

                ncfile.dimvar['t'][t] = otime * 10**iotimeprec

                if mode == 'xy':
                    ncfile.var[t,k,:,:] = fin.read(n)
                if mode == 'surf':
                    ncfile.var[t,:,:]   = fin.read(n)
                if mode == 'xz':
                    ncfile.var[t,:,k,:] = fin.read(n)
                if mode == 'yz':
                    ncfile.var[t,:,:,k] = fin.read(n)
        
                fin.close()
        ncfile.close() 
    except Exception as ex:
        print(ex)
        print("Failed to create %s.nc"%variable)
     
