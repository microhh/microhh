from netCDF4 import Dataset
import pylab as plt
import numpy as np
import sys

slist = 'hno3,co,hcho,rooh,h2o2,rh,no2,no,o3,ro2,ho2'
specs = slist.split(',')
profs = []
error = []
flux = []
with Dataset("tm5_ifs.default.0000000.nc") as x:
    zax = x.variables['z'][:]
    zaxh = x.variables['zh'][:]
    z = x.groups['default']
    for spec in specs:
        profs.append(z.variables[spec][:])
        spece = spec+'_2'
        error.append(z.variables[spece][:])
        if spec != specs[-1]:
            spece = spec+'_flux'
            flux.append(z.variables[spece][:])
    z = x.groups['thermo']
    thl = z.variables['thl'][:]
    thl_2 = z.variables['thl_2'][:]
profs = np.array(profs)

with Dataset("tm5_ifs.chemistry.0000000.nc") as x:
    zax = x.variables['z'][:]
    zaxh = x.variables['zh'][:]
    z = x.groups['default']
    print(z)
    rfa = z.variables['chem_budget'][:]
    ntime = rfa.shape[0]
nz = zax.size
nreact = 33
rfa = rfa.reshape(ntime,nz,nreact)
#  Vdot[10] = -A[6]-A[7]-A[10]+A[19]+A[25]-A[27];
rf = rfa[10,3,:]
tend = -rf[6]-rf[7]-rf[10]+rf[19]+rf[25]-rf[27]
print(-rf[6]/tend,-rf[7]/tend,-rf[10]/tend,+rf[19]/tend,+rf[25]/tend,-rf[27]/tend)
print(-rf[6],-rf[7],-rf[10],+rf[19],+rf[25],-rf[27])
print(-rf[6]-rf[7]-rf[10]+rf[19]+rf[25]-rf[27])


# OH production and loss:
poh = rf[1]+rf[7]+0.8*rf[16]+2*rf[18]+rf[22]+2*rf[23]
loh = rf[0]+rf[2]+rf[5]+rf[4]+rf[8]+rf[9]+0.6*rf[13]+rf[14]+rf[15]+rf[17]
print(poh,loh,poh/loh)
pho2 = rf[0] + rf[4] + rf[5] + rf[10] + rf[14] + rf[15] + 2*rf[20] + rf[22]
lho2 = rf[1] + rf[2] + 2*rf[3] + rf[7] + rf[11] +rf[12]
print(pho2,lho2,pho2/lho2)
print(rf[0]/pho2 , rf[4]/pho2 , rf[5]/pho2 , rf[10]/pho2 , rf[14]/pho2 , rf[15]/pho2 , 2*rf[20]/pho2 , rf[22]/pho2)
print(-rf[1]/lho2 ,-rf[2]/lho2 ,-2*rf[3]/lho2 ,-rf[7]/lho2, -rf[11]/lho2 ,-rf[12]/lho2)
#  Vdot[7] = A[0]-A[1]-A[2]-2*A[3]+A[4]+A[5]-A[7]+A[10]-A[11]-A[12] +A[14]+A[15]+2*A[20]+A[22];

