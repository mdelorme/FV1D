import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import h5py
import matplotlib

matplotlib.use('GTK3Agg')

if '--interval' in sys.argv:
  i = sys.argv.index('--interval')
  interval = int(sys.argv[i+1])
else:
  interval = 100

xmin = None
xmax = None

if '--xmin' in sys.argv:
  i = sys.argv.index('--xmin')
  xmin = float(sys.argv[i+1])

if '--xmax' in sys.argv:
  i = sys.argv.index('--xmax')
  xmax = float(sys.argv[i+1])

save_movie = '--save-movie' in sys.argv
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
    mask = None

    xmin_v = xmin if xmin else x[0]
    xmax_v = xmax if xmax else x[-1]

    mask = (xmin_v <= x) & (x <= xmax_v)

    rho_lim[0] = rho[mask].min()
    rho_lim[1] = rho[mask].max()
    vel_lim[0] = vel[mask].min()
    vel_lim[1] = vel[mask].max()
    prs_lim[0] = prs[mask].min()
    prs_lim[1] = prs[mask].max()
  else:
    rho_lim[0] = min(rho_lim[0], rho[mask].min())
    rho_lim[1] = max(rho_lim[1], rho[mask].max())
    vel_lim[0] = min(vel_lim[0], vel[mask].min())
    vel_lim[1] = max(vel_lim[1], vel[mask].max())
    prs_lim[0] = min(prs_lim[0], prs[mask].min())
    prs_lim[1] = max(prs_lim[1], prs[mask].max())

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

if no_anim:
  cur_ite = N-1
else:
  cur_ite = 0

mp = ax[1,0].scatter(time[cur_ite], mass[cur_ite], s=30, color='red')
Ekp = ax[1,1].scatter(time[cur_ite], Ek[cur_ite], s=30, color='red')
ep = ax[1,2].scatter(time[cur_ite], e[cur_ite], s=30, color='red')

base_ite = f'ite_{cur_ite}'

rho = np.array(f[f'{base_ite}/rho'])
vel = np.array(f[f'{base_ite}/vel'])
prs = np.array(f[f'{base_ite}/prs'])

rho_line = ax[0,0].plot(x, rho)[0]
vel_line = ax[0,1].plot(x, vel)[0]
prs_line = ax[0,2].plot(x, prs)[0]

xmin_v = xmin if xmin else x[0]
xmax_v = xmax if xmax else x[-1]
ax[0,0].set_xlim(xmin, xmax)
ax[0,1].set_xlim(xmin, xmax)
ax[0,2].set_xlim(xmin, xmax)

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
  ax[0,0].set_xlim(xmin, xmax)
  ax[0,1].set_xlim(xmin, xmax)
  ax[0,2].set_xlim(xmin, xmax)

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

if not no_anim:
  ani = animation.FuncAnimation(fig=fig, func=update, frames=N, interval=interval, repeat=loop)

if not save_movie or no_anim:
  plt.show()
else:
  print('Writing animation')
  ani.save(filename="animation.gif", writer="pillow")    

f.close()