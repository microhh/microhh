import numpy as np
import struct  as st
import pylab as pl
import os
import shutil
import subprocess

# Custom functions to do fast tri-linear interpolation:
from trilin import *

# Execute command line task
def execute(task):
    subprocess.call(task, shell=True, executable='/bin/bash')

# Replace value in namelist
def replace(filein, searchstring, value):
    value = str(value)
    execute('sed -i -e "s/\(' +searchstring+ r'\).*/\1=' +value+ '/g" ' + filein)

# Read variable from .ini file
def read_ini(filein, variable):
    value = None
    with open(filein) as f:
        for line in f:
           if(line.split('=')[0]==variable):
               value = line.split('=')[1]
    if(value == None):
        print('no variable %s'%variable)
    return value

def key_nearest(array, value):
    return np.abs(array-value).argmin()

# Empty class to wrap things in
class Empty:
    def __init__(self):
        pass

# Class to read single field
class Read_field:
    def __init__(self, path, name, time, nx, ny, nz):
        path = '%s%s.%07i'%(path,name,time)
        print(path)
        self.data = np.zeros((nz+2, ny+2, nx+2))

        fin   = open(path, "rb")
        for k in range(nz):
            raw = fin.read(nx*ny*8)
            tmp = np.array(st.unpack('{0}{1}d'.format('<', nx*ny), raw))
            self.data[k+1,1:-1,1:-1] = tmp.reshape((ny, nx))[:,:]

        # Set cyclic and top/bottom boundaries: 
        if(name in ['u','v','w']):
            # No slip bottom boundary
            self.data[0,:,:] = -self.data[1,:,:]   
            # No gradient top boundary
            self.data[-1,:,:] = self.data[-2,:,:]   
    
        elif(name in ['th']):
            # Constant gradient to surface, would be better to interpolate using Ts?
            self.data[0,:,:] = self.data[1,:,:] - (self.data[2,:,:] - self.data[1,:,:])
            # Constant gradient at top
            self.data[-1,:,:] = self.data[-2,:,:] + (self.data[-2] - self.data[-3])  

        self.data[:,0,:]  = self.data[:,-2,:]
        self.data[:,-1,:] = self.data[:,1,:]
        self.data[:,:,0]  = self.data[:,:,-2]
        self.data[:,:,-1] = self.data[:,:,1]

# Class to write field to restart file
def Write_field(data, path):
    tmp  = data.reshape(data.size)
    tmp2 = st.pack('{0}{1}d'.format('<', tmp.size), *tmp) 
    fout = open(path, "wb")
    fout.write(tmp2*8)
    fout.close()  

# Class to read grid.0000000 file
class Read_grid:
    def __init__(self, filename, nx, ny, nz):
        n   = nx*ny*nz
        fin = open(filename, "rb")
        raw = fin.read(nx*8)
        self.x  = np.array(st.unpack('{0}{1}d'.format('<', nx), raw))
        raw = fin.read(nx*8)
        self.xh = np.array(st.unpack('{0}{1}d'.format('<', nx), raw))
        raw = fin.read(ny*8)
        self.y  = np.array(st.unpack('{0}{1}d'.format('<', ny), raw))
        raw = fin.read(ny*8)
        self.yh = np.array(st.unpack('{0}{1}d'.format('<', ny), raw))
        raw = fin.read(nz*8)
        self.z  = np.array(st.unpack('{0}{1}d'.format('<', nz), raw))
        raw = fin.read(nz*8)
        self.zh = np.array(st.unpack('{0}{1}d'.format('<', nz), raw))
        fin.close()

        # Add ghost cells to dim_in
        for var in ['x','xh','y','yh','z','zh']:
            tmp1 = getattr(self, var)
            tmp2 = np.zeros(tmp1.size+2)
            tmp2[1:-1] = tmp1
            tmp2[0 ] = tmp2[1]  - (tmp2[2]-tmp2[1])
            tmp2[-1] = tmp2[-2] + (tmp2[-2]-tmp2[-3]) 
            setattr(self, var+'_p', tmp2) 

# Class to read input soundings and write it back   
class Sounding:
    def __init__(self, filename):   
        self.filename = filename
        snd = np.loadtxt(filename, skiprows=1)
        self.z  = snd[:,0]
        self.th = snd[:,1]
        self.u  = snd[:,2]
        self.ug = snd[:,3]
        self.v  = snd[:,4]
        self.vg = snd[:,5]
  
    def write_back(self):
        proffile = open(self.filename, 'w')
        proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s}\n'.format('z','th','u','ug','v','vg'))
        for k in range(self.z.size):
            proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} \n'.format(self.z[k], self.th[k], self.u[k], self.ug[k], self.v[k], self.vg[k]))
        proffile.close()
    
# -----------------
# Main part of script
# -----------------
if(__name__ == "__main__"): 
    pl.close('all')

    # Settings of INPUT files
    dir_in  = os.getcwd()
    time_in = 32400
    nx_in   = 256
    ny_in   = 256
    nz_in   = 304
 
    # Settings of OUTPUT files
    dir_out = 'restart/'

    # Read settings from restart .ini file
    nx_out = int(read_ini(dir_out+'gabls4s3.ini', 'itot'))
    ny_out = int(read_ini(dir_out+'gabls4s3.ini', 'jtot'))
    nz_out = int(read_ini(dir_out+'gabls4s3.ini', 'ktot'))

    # Re-run init mode to create grid file and FFTW plan
    os.chdir(dir_out)
    execute('./microhh init gabls4s3')
    os.chdir(dir_in)

    # Read input and output grid:
    grid_in  = Read_grid('grid.0000000',         nx_in,  ny_in,  nz_in )  
    grid_out = Read_grid(dir_out+'grid.0000000', nx_out, ny_out, nz_out)  

    # Read input fields:
    fields_in     = Empty()
    fields_in.u   = Read_field('', 'u',  time_in, nx_in, ny_in, nz_in)
    fields_in.v   = Read_field('', 'v',  time_in, nx_in, ny_in, nz_in)
    fields_in.w   = Read_field('', 'w',  time_in, nx_in, ny_in, nz_in)
    fields_in.th  = Read_field('', 'th', time_in, nx_in, ny_in, nz_in)

    # Create new empty fields:
    fields_out    = Empty()
    fields_out.u  = np.zeros((nz_out, ny_out, nx_out))
    fields_out.v  = np.zeros((nz_out, ny_out, nx_out))
    fields_out.w  = np.zeros((nz_out, ny_out, nx_out))
    fields_out.th = np.zeros((nz_out, ny_out, nx_out))
   
    # Interpolate to new grid
    xi  = np.zeros(nx_out, dtype=np.int)
    xhi = np.zeros(nx_out, dtype=np.int)
    yi  = np.zeros(ny_out, dtype=np.int)
    yhi = np.zeros(ny_out, dtype=np.int)
    zi  = np.zeros(nz_out, dtype=np.int)
    zhi = np.zeros(nz_out, dtype=np.int)
    
    xo  = np.zeros(nx_out)
    xho = np.zeros(nx_out)
    yo  = np.zeros(ny_out)
    yho = np.zeros(ny_out)
    zo  = np.zeros(nz_out)
    zho = np.zeros(nz_out)
   
    # Function to find left most index (ii) and its relative offset (offset) from that 
    def index(array, value):
        ii = np.abs(array - value).argmin()
        if(array[ii] > value):
            ii -= 1
        offset = (value - array[ii]) / (array[ii+1] - array[ii])
        return ii, offset
    
    # Calculate left most index in original array, and offset from that point:
    for i in range(nx_out):
        xi [i], xo [i] = index(grid_in.x_p,  grid_out.x[i] ) 
        xhi[i], xho[i] = index(grid_in.xh_p, grid_out.xh[i]) 
    for i in range(ny_out):
        yi [i], yo [i] = index(grid_in.y_p,  grid_out.y[i] ) 
        yhi[i], yho[i] = index(grid_in.yh_p, grid_out.yh[i]) 
    for i in range(nz_out):
        zi [i], zo [i] = index(grid_in.z_p,  grid_out.z[i] ) 
        zhi[i], zho[i] = index(grid_in.zh_p, grid_out.zh[i]) 
    
    # Interpolate fields. Change to trilin_python if Cython isn't installed
    trilin_cython(fields_in.th.data, fields_out.th, xi,  yi,  zi,  xo,  yo,  zo )
    trilin_cython(fields_in.u.data,  fields_out.u,  xhi, yi,  zi,  xho, yo,  zo )
    trilin_cython(fields_in.v.data,  fields_out.v,  xi,  yhi, zi,  xo,  yho, zo )
    trilin_cython(fields_in.w.data,  fields_out.w,  xi,  yi,  zhi, xo,  yo,  zho)

    # Write interpolated field to binary restart file
    Write_field(fields_out.th, '%s%s.%07i'%(dir_out, 'th', time_in))
    Write_field(fields_out.u,  '%s%s.%07i'%(dir_out, 'u',  time_in))
    Write_field(fields_out.v,  '%s%s.%07i'%(dir_out, 'v',  time_in))
    Write_field(fields_out.w,  '%s%s.%07i'%(dir_out, 'w',  time_in))

    # Overwrite input sounding as the buffer acts on the initial profiles:
    snd = Sounding(dir_out+'gabls4s3.prof') 
 
    for k in range(snd.z.size):
        snd.th[k] = fields_out.th[k,:,:].mean() 
        snd.u[k]  = fields_out.u [k,:,:].mean() 
        snd.v[k]  = fields_out.v [k,:,:].mean() 
  
    snd.write_back() 
