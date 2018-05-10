dthetadz = 0.003
N = (10./300.*dthetadz)**.5
H = 1e3
f = 1e-3
shear = 1e-2

L_d = N*H/f
print("Rossby radius L_d = {} km".format(1e-3*L_d))

T_e = N / (f*shear)
print("Eady time scale T_e = {} s".format(T_e))

sigma_e = 0.31 * shear * f / N
print("Eady growth rate sigma_e = {} s-1, or {} d-1".format(sigma_e, sigma_e * 3600.))

### Stone 1971
u0 = shear*H
Ri = (H*N/u0)**2
print("(STONE71) Richardson number Ri = {}".format(Ri))

delta = f*H/u0
print("(STONE71) Aspect ratio to Rossby = {}".format(delta))

