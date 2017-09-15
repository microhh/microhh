import numpy as np
import matplotlib.pylab as pl
import pandas as pd
import glob
import datetime

from microhh_tools import Read_statistics

pl.close('all')

MDT_to_UTC = datetime.timedelta(hours=6)   # MDT to UTC conversion

def components_to_direction(u, v):
    """ Convert u and v components into direction and absolute velocity """
    U = (u**2. + v**2)**0.5
    D = 270 - np.rad2deg(np.arctan2(v, u))
    return D,U

def direction_to_components(D, U):
    """ Convert direction and absolute velocity into u and v components """
    u = -U * np.sin(np.deg2rad(D)) 
    v = -U * np.cos(np.deg2rad(D)) 
    return u,v

def esat(T):
    return 0.611e3 * np.exp(17.2694 * (T - 273.16) / (T - 35.86))

def qsat(T,p):
    return 0.622 * esat(T) / p

class Read_balloon:
    """ Class to read the tethersonde data into a Pandas DataFrame """
    def __init__(self, file_name):
        print('Reading {}'.format(file_name))

        elev_site = 1310    

        cols = ['z','p','T','RH','U','dir_n','th','thv','ri']
        self.df = pd.read_csv(file_name, delimiter='\t', names=cols, skiprows=1)
        
        self.df['T']    += 273.15
        self.df['p']    *= 100.
        self.df['z_agl'] = self.df['z'] - elev_site
        self.df['exner'] = (self.df['p'] / 1e5)**(287.05 / 1004.)
        self.df['q']     = qsat(self.df['T'], self.df['p']) * self.df['RH'] / 100.

        self.df['u'], self.df['v'] = direction_to_components(self.df['dir_n'], self.df['U'])

        date = '2001'+file_name.split('/')[-1].split('.')[0]
        self.datetime = datetime.datetime.strptime(date, '%Y%m%d%H%M') + MDT_to_UTC


if (__name__ == "__main__"):
    pl.close('all')

    # Read the reference (tethersonde) observations   
    files = glob.glob('reference_data/0924*')
    balloons = []
    for f in files:
        balloons.append( Read_balloon(f) )

    # MicroHH output
    mhh = Read_statistics('must.default.0000000.nc')
    
    # Define colors for each balloon time
    cc = pl.cm.viridis(np.linspace(0,1,len(files)))

    
    pl.figure(figsize=[10,7])

    pl.subplot(221)
    for i,b in enumerate(balloons):
        # Find nearest time in MicroHH data
        time = (b.datetime - balloons[0].datetime).total_seconds()
        t = np.abs(mhh.t - time).argmin()

        pl.plot(b.df['th'], b.df['z_agl'], color=cc[i], label=b.datetime)
        pl.plot(mhh['thl'][t,:], mhh['z'], color=cc[i], dashes=[2,2])

    pl.legend(loc='best', frameon=False, fontsize=8)
    pl.ylim(0,400)
    pl.ylabel('z (m)')
    pl.xlabel('theta (K)')

    pl.subplot(222)
    for i,b in enumerate(balloons):
        # Find nearest time in MicroHH data
        time = (b.datetime - balloons[0].datetime).total_seconds()
        t = np.abs(mhh.t - time).argmin()

        pl.plot(b.df['q']*1000, b.df['z_agl'], color=cc[i], label=b.datetime)
        pl.plot(mhh['qt'][t,:]*1000, mhh['z'], color=cc[i], dashes=[2,2])

    pl.ylim(0,400)
    pl.ylabel('z (m)')
    pl.xlabel('q (g kg-1)')

    pl.subplot(223)
    for i,b in enumerate(balloons):
        # Find nearest time in MicroHH data
        time = (b.datetime - balloons[0].datetime).total_seconds()
        t = np.abs(mhh.t - time).argmin()

        pl.plot(b.df['u'], b.df['z_agl'], color=cc[i], label=b.datetime)
        pl.plot(mhh['u'][t,:], mhh['z'], color=cc[i], dashes=[2,2])

    pl.legend(loc='best', frameon=False, fontsize=8)
    pl.ylim(0,400)
    pl.ylabel('z (m)')
    pl.xlabel('theta (K)')

    pl.subplot(224)
    for i,b in enumerate(balloons):
        # Find nearest time in MicroHH data
        time = (b.datetime - balloons[0].datetime).total_seconds()
        t = np.abs(mhh.t - time).argmin()

        pl.plot(b.df['v'], b.df['z_agl'], color=cc[i], label=b.datetime)
        pl.plot(mhh['v'][t,:], mhh['z'], color=cc[i], dashes=[2,2])

    pl.legend(loc='best', frameon=False, fontsize=8)
    pl.ylim(0,400)
    pl.ylabel('z (m)')
    pl.xlabel('theta (K)')
