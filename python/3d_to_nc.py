import microhh_tools as mht     # available in microhh/python directory
import argparse

#Parse command line and namelist options
parser = argparse.ArgumentParser(description='Convert microHH 3D binary to netCDF4 files.')
parser.add_argument('-f', '--filename', help='ini file name')
parser.add_argument('-v', '--vars', nargs='*', help='variable names')
parser.add_argument('-i', '--iotimeprec', help='iotimeprec', type=int)
parser.add_argument('-e', '--endian', help='endianess', choices = ['little', 'big'],default='little')
parser.add_argument('-p', '--precision', help='precision', choices = ['single', 'double'],default='single')

args = parser.parse_args()


nl = mht.Read_namelist(args.filename)
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']
starttime  = nl['grid'] ['starttime']
endtime    = nl['grid'] ['endtime']
sampletime = nl['dump']['sampletime']

variables = args.vars if args.vars is not None else nl['dump']['dumplist']

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
        filename = "{0}.nc".format(variable)
        dim = {'time' : range(niter), 'x' : range(itot), 'y' : range(jtot), 'z': range(ktot)}
        n = itot * jtot * ktot
        ncfile = mht.Create_ncfile(grid, filename, variable, dim)

        # Loop through the files and read 3d field
        for t in range(niter):
            otime = int((starttime + t*sampletime) / 10**iotimeprec)
            f_in  = "{0:}.{1:07d}".format(variable, otime)

            try:
                fin = mht.Read_binary(f_in, endian=endian, precision=precision)
            except:
                print('Stopping: cannot find file {}'.format(f_in))
                ncfile.sync()
                stop = True
                break

            print("Processing %8s, time=%7i"%(variable, otime))

            ncfile.dimvar['t'] = otime * 10**iotimeprec
            ncfile.var[t,:,:,:] = fin.read(n)

            fin.close()
            ncfile.close() 
    except:
        print("Failed to create %s.nc"%variable)
        ncfile.close()
        