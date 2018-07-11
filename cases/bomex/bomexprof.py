import numpy as np

# Get number of vertical levels and size from .ini file
with open('bomex.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z     = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl   = np.zeros(np.size(z))
qt    = np.zeros(np.size(z))
u     = np.zeros(np.size(z))
ug    = np.zeros(np.size(z))
wls   = np.zeros(np.size(z))
thlls = np.zeros(np.size(z))
qtls  = np.zeros(np.size(z))

for k in range(kmax):
    # temperature
    if(z[k] <= 520.):
        thl[k] = 298.7
    elif(z[k] <= 1480):
        thl[k] = 298.7 + (z[k]-520.)*(302.4-298.7)/(1480.-520.)
    elif(z[k] <= 2000):
        thl[k] = 302.4 + (z[k]-1480.)*(308.2-302.4)/(2000.-1480.)
    else:
        thl[k] = 308.2 + (z[k]-2000.)*(311.85-308.2)/(3000.-2000.)

    # specific humidity
    if(z[k] <= 520.):
        qt[k] = 1e-3*(17.0 + z[k]*(16.3-17.0)/520.)
    elif(z[k] <= 1480):
        qt[k] = 1.e-3*(16.3 + (z[k]-520.)*(10.7-16.3)/(1480.-520.))
    elif(z[k] <= 2000):
        qt[k] = 1.e-3*(10.7 + (z[k]-1480.)*(4.2-10.7)/(2000.-1480.))
    else:
        qt[k] = 1.e-3*(4.2 + (z[k]-2000.)*(3.-4.2)/(3000.-2000.))


    # u-wind component
    if(z[k] <= 700.):
        u[k] = -8.75
    else:
        u[k] = -8.75 + (z[k]-700.)*(-4.61+8.75)/(3000.-700.)

    # ug-wind component
    ug[k] = -10. + 1.8e-3*z[k]

    # large scale vertical velocity
    if(z[k] <= 1500):
        wls[k] = z[k]*(-0.65)/1500.
    elif(z[k] <= 2100):
        wls[k] = -0.65 + (z[k]-1500)*(0.65)/(2100.-1500.)

    # large scale temperature tendency
    if(z[k] <= 1500):
        thlls[k] = (-2.)
    else:
        thlls[k] = (-2.) + (z[k]-1500)*(2.)/(3000.-1500.)

    # large scale moisture tendency
    if(z[k] <= 300):
        qtls[k] = -1.2
    elif(z[k] <= 500):
        qtls[k] = -1.2 + (z[k]-300)*(1.2)/(500.-300)

# normalize profiles to SI
#qtls /= 1000.  # from g/kg to kg/kg
wls  /= 100.   # from cm/s to m/s
thlls  /= 86400. # from K/d to K/s
qtls *= 1.e-8

# write the data to a file
proffile = open('bomex.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s} {6:^20s} {7:^20s}\n'.format('z','thl','qt','u','ug','wls','thlls','qtls'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E} {6:1.14E} {7:1.14E}\n'.format(z[k], thl[k], qt[k], u[k], ug[k], wls[k], thlls[k], qtls[k]))
proffile.close()
