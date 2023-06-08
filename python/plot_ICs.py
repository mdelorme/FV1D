import numpy as np
import matplotlib.pyplot as plt
import h5py

f = h5py.File('run.h5', 'r')

x   = np.array(f['x'])
rho = np.array(f['ite_0/rho'])
u   = np.array(f['ite_0/vel']) / rho

Ek = 0.5 * rho * u*u

gamma0 = 5.0/3.0
prs = (np.array(f['ite_0/prs']) - Ek) * (gamma0-1.0)

fig, ax = plt.subplots(3, 1, figsize=(7, 10))

ax[0].plot(x, rho)
ax[1].plot(x, u)
ax[2].plot(x, prs)

ax[2].set_xlabel('x')
ax[0].set_ylabel('Density')
ax[1].set_ylabel('Velocity')
ax[2].set_ylabel('Pressure')

f.close()

plt.show()