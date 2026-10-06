import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

D1 = np.linspace(0.5, 300, 256) * 1e-6              # Smooth curve
D2 = np.array([1, 6, 15, 39, 70.5, 261.5]) * 1e-6   # LES sizes

rho_p = 1500   # Density particles [kg m-3]
rho_a = 1.225  # Reference density air [kg m-3]
nu = 1e-5      # Kinematic viscosity air [m2 s-1]
g = 9.81       # Gravitational acceleration [m s-2]

tau_p1 = D1**2 * rho_p / (18 * nu * rho_a)
w_s1 = -tau_p1 * g

tau_p2 = D2**2 * rho_p / (18 * nu * rho_a)
w_s2 = -tau_p2 * g

plt.figure(figsize=(8,4))

ax = plt.subplot(121)
plt.plot(D1*1e6, -w_s1)
plt.scatter(D2*1e6, -w_s2, facecolor='none', edgecolor='k', label='LES sizes')
plt.xlabel('D ($\mu$m)')
plt.ylabel('w (m s$^{-1}$)')
plt.grid()
plt.legend()

ax = plt.subplot(122)
plt.plot(D1*1e6, -w_s1)
plt.scatter(D2*1e6, -w_s2, facecolor='none', edgecolor='k')
plt.xlabel('D ($\mu$m)')
plt.ylabel('w (m s$^{-1}$)')
plt.grid()
ax.set_yscale('log')
ax.set_xscale('log')

plt.tight_layout()
plt.savefig('sedimentation_velocity.png')
