import numpy
import struct
import netCDF4

from pylab import *

def plotstats(name, line):
  stats = netCDF4.Dataset("bomex.{}.0000000.nc".format(name),"r")
  t   = stats.variables["t"][:]
  z   = stats.variables["z"][:]
  zh  = stats.variables["zh"][:]

  areat  = stats.variables["area"][:,:]
  
  st  = stats.variables["s"][:,:]
  qtt = stats.variables["qt"][:,:]*1000.
  bt  = stats.variables["b"][:,:]
  ut  = stats.variables["u"][:,:]
  vt  = stats.variables["v"][:,:]
  qlt = stats.variables["ql"][:,:]*1000.
  cft = stats.variables["cfrac"][:,:]
  
  sfluxt = stats.variables["sflux"][:,:]
  bfluxt = stats.variables["bflux"][:,:]
  ufluxt = stats.variables["uflux"][:,:]
  vfluxt = stats.variables["vflux"][:,:]
  Ufluxt = ufluxt + vfluxt
  
  u2t  = stats.variables["u2"][:,:]
  v2t  = stats.variables["v2"][:,:]
  w2t  = stats.variables["w2"][:,:]
  tket = 0.5*(u2t + v2t + 0.5*(w2t[:,0:-1]+w2t[:,1::]))
  
  end   = t.size
  start = t.size - 12
  
  area = numpy.mean(areat[start:end,:], 0)

  s  = numpy.mean(st [start:end,:], 0)
  qt = numpy.mean(qtt[start:end,:], 0)
  b  = numpy.mean(bt [start:end,:], 0)
  u  = numpy.mean(ut [start:end,:], 0)
  v  = numpy.mean(vt [start:end,:], 0)
  ql = numpy.mean(qlt[start:end,:], 0)
  cf = numpy.mean(cft[start:end,:], 0)
  
  sflux = numpy.mean(sfluxt[start:end,:], 0)
  bflux = numpy.mean(bfluxt[start:end,:], 0)
  uflux = numpy.mean(ufluxt[start:end,:], 0)
  vflux = numpy.mean(vfluxt[start:end,:], 0)
  Uflux = numpy.mean(Ufluxt[start:end,:], 0)
  
  u2  = numpy.mean(u2t [start:end,:], 0)
  v2  = numpy.mean(v2t [start:end,:], 0)
  w2  = numpy.mean(w2t [start:end,:], 0)
  tke = numpy.mean(tket[start:end,:], 0)
  
  # enable LaTeX plotting
  rc('font',**{'family':'serif','serif':['Palatino']})
  rc('text', usetex=True)

  f = 1

  figure(f)
  for n in range(start,end):
    plot(areat[n,:], z, color='#eeeeee')
  plot(area, z, line, label=name)
  xlabel(r'area coverage [-]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)

  f += 1
  figure(f)
  for n in range(start,end):
    plot(st[n,:], z, color='#eeeeee')
  plot(s, z, line, label=name)
  plot(st[0,:], z, 'k:')
  xlabel(r'$\theta$ [K]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
 
  f += 1
  figure(f)
  for n in range(start,end):
    plot(qtt[n,:], z, color='#eeeeee')
  plot(qt, z, line, label=name)
  plot(qtt[0,:], z, 'k:')
  xlabel(r'q$_t$ [g~kg$^{-1}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)

  f += 1
  figure(f)
  for n in range(start,end):
    plot(bt[n,:], z, color='#eeeeee')
  plot(b, z, line, label=name)
  plot(bt[0,:], z, 'k:')
  xlabel(r'b [m~s$^{-2}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)

  f += 1
  figure(f)
  for n in range(start,end):
    plot(ut[n,:], z, color='#eeeeee')
  plot(u, z, line, label=name)
  plot(ut[0,:], z, 'k:')
  xlabel(r'u [m~s$^{-1}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
  
  f += 1
  figure(f)
  for n in range(start,end):
    plot(vt[n,:], z, color='#eeeeee')
  plot(v, z, line, label=name)
  plot(vt[0,:], z, 'k:')
  xlabel(r'v [m~s$^{-1}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
  
  f += 1
  figure(f)
  for n in range(start,end):
    plot(qlt[n,:], z, color='#eeeeee')
  plot(ql, z, line, label=name)
  xlabel(r'q$_l$ [g~kg$^{-1}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
  
  f += 1
  figure(f)
  for n in range(start,end):
    plot(cft[n,:], z, color='#eeeeee')
  plot(cf, z, line, label=name)
  xlabel(r'cloud fraction [-]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
  
  f += 1
  figure(f)
  for n in range(start,end):
    plot(sfluxt[n,:], zh, color='#eeeeee')
  plot(sflux, zh, label=name)
  xlabel(r'w`$\theta_l$` [K~m~s$^{-1}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
  
  f += 1
  figure(f)
  for n in range(start,end):
    plot(bfluxt[n,:], zh, color='#eeeeee')
  plot(bflux, zh, label=name)
  xlabel(r'w`b` [m$^2$~s$^{-3}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
  
  f += 1
  figure(f)
  for n in range(start,end):
    plot(ufluxt[n,:], zh, color='#eeeeee')
    #plot(vfluxt[n,:], zh, color='#eeeeee')
  plot(uflux, zh, label=name)
  #plot(vflux, zh)
  #plot(Uflux, zh)
  xlabel(r'u`w` [m$^2$~s$^{-2}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)

  f += 1
  figure(f)
  for n in range(start,end):
    plot(u2t[n,:], z, color='#eeeeee')
  plot(u2, z, label=name)
  xlabel(r'u`$^2$ [m$^2$~s$^{-2}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
 
  f += 1
  figure(f)
  for n in range(start,end):
    plot(v2t[n,:], z, color='#eeeeee')
  plot(v2, z, label=name)
  xlabel(r'v`$^2$ [m$^2$~s$^{-2}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
 
  f += 1
  figure(f)
  for n in range(start,end):
    plot(w2t[n,:], zh, color='#eeeeee')
  plot(w2, zh, label=name)
  xlabel(r'w`$^2$ [m$^2$~s$^{-2}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)
  
  f += 1
  figure(f)
  for n in range(start,end):
    plot(tket[n,:], z, color='#eeeeee')
  plot(tke, z, line, label=name)
  xlabel(r'TKE [m$^2$~s$^{-2}$]')
  ylabel(r'z [m]')
  legend(loc=0, frameon=False)

ioff()
plotstats("default", 'k-')
plotstats("wplus", 'b-')
plotstats("ql", 'r-')
plotstats("qlcore", 'c-')
ion()
show()

