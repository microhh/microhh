N = 1.
H = 100.
f = 1e-4
shear = 1e-4

L_d = N*H/f
print("Rossby radius L_d = {:.3E} m".format(L_d))

T_e = N / (f*shear)
print("Eady time scale T_e = {:.3E} s".format(T_e))

sigma_e = 0.31 * shear * f / N
print("Eady growth rate sigma_e = {:.3E} s-1".format(sigma_e))

