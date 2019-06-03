import os
import microhh_tools as mht     # available in microhh/python directory
import argparse

#Parse command line and namelist options
cross_modes = ['xy', 'xz', 'yz']
parser = argparse.ArgumentParser(description='Convert MicroHH binary cross-sections to netCDF4 files.')
parser.add_argument('-m', '--modes', nargs='*', help = 'mode of the cross section', choices = cross_modes)
parser.add_argument('-f', '--filename', help='ini file name')
parser.add_argument('-v', '--vars', nargs='*', help='variable names')
parser.add_argument('-x', '--index', nargs='*', help='indices')
parser.add_argument('-t0', '--starttime', help='first time step to be parsed')
parser.add_argument('-t1', '--endtime', help='last time step to be parsed')
parser.add_argument('-tstep', '--sampletime', help='time interval to be parsed')
parser.add_argument('-p', '--precision', help='precision', choices = ['single', 'double'])

args = parser.parse_args()

modes = args.modes
indexes = args.index

nl = mht.Read_namelist(args.filename)
itot = nl['grid']['itot']
jtot = nl['grid']['jtot']
ktot = nl['grid']['ktot']
starttime  = args.starttime  if args.starttime  is not None else nl['time']['starttime']
endtime    = args.endtime    if args.endtime    is not None else nl['time']['endtime']
sampletime = args.sampletime if args.sampletime is not None else nl['cross']['sampletime']

if args.modes is None:
    modes = list(nl['cross'].keys() & cross_modes)
else:
    modes = args.modes
try:
    iotimeprec = nl['time']['iotimeprec']
except KeyError:
    iotimeprec = 0.

variables = args.vars if args.vars is not None else nl['cross']['crosslist']
precision = args.precision
# End option parsing

# Calculate the number of iterations
niter = int((endtime-starttime) / sampletime + 1)

grid = mht.Read_grid(itot, jtot, ktot)

# Loop over the different variables and crosssections
for mode in modes:
    for variable in variables:
        try:
            otime = round(starttime / 10**iotimeprec)
            if os.path.isfile("{0}.xy.{1:07d}".format(variable, otime)):
                if mode is not 'xy':
                    continue
                at_surface = True
            else:
                at_surface = False

            filename = "{0}.{1}.nc".format(variable,mode)
            if not at_surface:
                if indexes == None:
                    indexes_local = mht.get_cross_indices(variable, mode)
                else:
                    indexes_local = indexes

            dim = {'time' : range(niter), 'z' : range(ktot), 'y' : range(jtot), 'x': range(itot)}
            if at_surface:
                dim.pop('z')
                n = itot * jtot
                indexes_local = [-1]
            elif mode == 'xy':
                dim.update({'z' : indexes_local})
                n = itot * jtot
            elif mode == 'xz':
                dim.update({'y' : indexes_local})
                n = itot * ktot
            elif mode == 'yz':
                dim.update({'x' : indexes_local})
                n = ktot * jtot

            ncfile = mht.Create_ncfile(grid, filename, variable, dim, precision)
            for t in range(niter):
                for k in range(len(indexes_local)):
                    index = indexes_local[k]
                    otime = round((starttime + t*sampletime) / 10**iotimeprec)
                    if at_surface:
                        f_in  = "{0}.{1}.{2:07d}".format(variable, mode, otime)
                    else:
                        f_in  = "{0:}.{1}.{2:05d}.{3:07d}".format(variable, mode, index, otime)
                    try:
                        fin = mht.Read_binary(grid, f_in)
                    except Exception as ex:
                        print (ex)
                        raise Exception('Stopping: cannot find file {}'.format(f_in))

                    print("Processing %8s, time=%7i, index=%4i"%(variable, otime, index))

                    ncfile.dimvar['time'][t] = otime * 10**iotimeprec

                    if at_surface:
                        ncfile.var[t,:,:]   = fin.read(n)
                    elif mode == 'xy':
                        ncfile.var[t,k,:,:] = fin.read(n)
                    elif mode == 'xz':
                        ncfile.var[t,:,k,:] = fin.read(n)
                    elif mode == 'yz':
                        ncfile.var[t,:,:,k] = fin.read(n)

                    fin.close()
            ncfile.close()
        except Exception as ex:
            print(ex)
            print("Failed to create %s"%filename)
