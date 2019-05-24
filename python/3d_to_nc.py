import microhh_tools as mht     # available in microhh/python directory
import argparse

#Parse command line and namelist options
parser = argparse.ArgumentParser(description='Convert microHH 3D binary to netCDF4 files.')
parser.add_argument('-f', '--filename', help='ini file name')
parser.add_argument('-v', '--vars', nargs='*', help='variable names')
parser.add_argument('-e', '--endian', help='endianess', choices = ['little', 'big'],default='little')
parser.add_argument('-p', '--precision', help='precision', choices = ['single', 'double'],default='double')
parser.add_argument('-t0',    '--starttime', help='first time step to be parsed')
parser.add_argument('-t1',    '--endtime', help='last time step to be parsed')
parser.add_argument('-tstep', '--sampletime', help='time interval to be parsed')
args = parser.parse_args()

nl = mht.Read_namelist(args.filename)
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']
starttime  = args.starttime  if args.starttime  is not None else nl['time']['starttime']
endtime    = args.endtime    if args.endtime    is not None else nl['time']['endtime']
sampletime = args.sampletime if args.sampletime is not None else nl['dump']['sampletime']

try:
    iotimeprec = nl['time']['iotimeprec']
except KeyError:
    iotimeprec = 0.

variables = args.vars if args.vars is not None else nl['dump']['dumplist']

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
            otime = round((starttime + t*sampletime) / 10**iotimeprec)
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
        print("Failed to create %s"%filename)
        ncfile.close()
