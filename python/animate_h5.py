import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import h5py

if '--interval' in sys.argv:
  i = sys.argv.index('--interval')
  interval = int(sys.argv[i+1])
else:
  interval = 100

no_anim = '--no-anim' in sys.argv
loop = '--loop' in sys.argv

f = h5py.File('run.h5', 'r')

T = []

N = len(f.keys())-1

x = np.array(f['x'])

first = True
rho_lim = [0]*2
vel_lim = [0]*2
prs_lim = [0]*2

time = []
mass = []
e    = []
Ek   = []

gamma0 = 5.0/3.0

dx = x[1]-x[0]

for i in range(N):
  path = f'ite_{i}'
  group = f[path]
  rho = np.array(group['rho'])
  vel = np.array(group['vel'])
  prs = np.array(group['prs'])

  if first:
    first = False

    rho_lim[0] = rho.min()
    rho_lim[1] = rho.max()
    vel_lim[0] = vel.min()
    vel_lim[1] = vel.max()
    prs_lim[0] = prs.min()
    prs_lim[1] = prs.max()
  else:
    rho_lim[0] = min(rho_lim[0], rho.min())
    rho_lim[1] = max(rho_lim[1], rho.max())
    vel_lim[0] = min(vel_lim[0], vel.min())
    vel_lim[1] = max(vel_lim[1], vel.max())
    prs_lim[0] = min(prs_lim[0], prs.min())
    prs_lim[1] = max(prs_lim[1], prs.max())

  time.append(group.attrs['time'])
  mass_ = np.sum(rho*dx)
  Ek_ = np.sum(0.5 * rho * vel**2.0)
  e_  = np.sum(prs / (gamma0-1.0))

  mass.append(mass_)
  Ek.append(Ek_)
  e.append(e_)


for lim in (rho_lim, vel_lim, prs_lim):
  lim[0] *= (0.9 if lim[0] > 0.0 else 1.1)
  lim[1] *= 1.1
  
fig, ax = plt.subplots(2, 3, figsize=(14, 10))

ax[0,0].set_ylim(rho_lim[0], rho_lim[1])
ax[0,1].set_ylim(vel_lim[0], vel_lim[1])
ax[0,2].set_ylim(prs_lim[0], prs_lim[1])

ax[1,0].plot(time, mass)
ax[1,1].plot(time, Ek)
ax[1,2].plot(time, e)
mp = ax[1,0].scatter(time[0], mass[0], s=30, color='red')
Ekp = ax[1,1].scatter(time[0], Ek[0], s=30, color='red')
ep = ax[1,2].scatter(time[0], e[0], s=30, color='red')

rho = np.array(f['ite_0/rho'])
vel = np.array(f['ite_0/vel'])
prs = np.array(f['ite_0/prs'])

rho_line = ax[0,0].plot(x, rho)[0]
vel_line = ax[0,1].plot(x, vel)[0]
prs_line = ax[0,2].plot(x, prs)[0]
plt.suptitle(f'Problem = {f.attrs["problem"]}; ite = 0')

def update(frame):
  rho = np.array(f[f'ite_{frame}/rho'])
  vel = np.array(f[f'ite_{frame}/vel'])
  prs = np.array(f[f'ite_{frame}/prs'])

  rho_line.set_ydata(rho)
  vel_line.set_ydata(vel)
  prs_line.set_ydata(prs)

  ax[0,0].set_ylim(rho_lim[0], rho_lim[1])
  ax[0,1].set_ylim(vel_lim[0], vel_lim[1])
  ax[0,2].set_ylim(prs_lim[0], prs_lim[1])

  mp.set_offsets((time[frame], mass[frame]))
  Ekp.set_offsets((time[frame], Ek[frame]))
  ep.set_offsets((time[frame], e[frame]))
  plt.suptitle(f'Problem = {f.attrs["problem"]}; ite = {frame}')

ax[0,0].set_xlabel('x')
ax[0,1].set_xlabel('x')
ax[0,2].set_xlabel('x')
ax[0,0].set_title('Density')
ax[0,1].set_title('Velocity')
ax[0,2].set_title('Pressure')
ax[1,0].set_xlabel('Time')
ax[1,1].set_xlabel('Time')
ax[1,2].set_xlabel('Time')
ax[1,0].set_title('Mass')
ax[1,1].set_title('Kinetic energy')
ax[1,2].set_title('Internal energy')

plt.tight_layout()

ani = animation.FuncAnimation(fig=fig, func=update, frames=N, interval=interval, repeat=loop)

if no_anim:
  plt.show()
else:
  print('Writing animation')
  ani.save(filename="animation.gif", writer="pillow")    

f.close()