from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
Nxy = np.array([4,
8,
10,
12,
14,
16,
18,
32,
64
])
times_ms = np.array([ 25971,

46545,

5841,

8631,

11696,

15546,

18321,

49847,

192075,
])

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.title("Time taken")
ax1.set_ylim(5000,200000)
ax1.set_ylabel("Time taken /s")
ax1.set_xlabel("$N_{xy}$")
ax1.set_xlim(0,64)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.set_yscale('log')
ax1.plot(Nxy,times_ms,label="Multiscat",color=[0., 0.5, 1.],linestyle=(0,(8,10)),marker='+',markersize=15)
plt.show()
Nxy = np.array([
10,
12,
14,
16,
18,
32,
64
])
intensity_specular = np.array([0.0603913,
0.0604757,
0.0604173,
0.0603834,
0.0603591,
0.0603817,
0.0603818,
])
intensity_fws = np.array([0.0914582,

0.0914012,

0.0913918,

0.0913966,

0.0913918,

0.0913970,

0.0913969
])
intensity_bk = np.array([0.010028,

0.0100523,

0.0100629,

0.0100657,

0.0100639,

0.0100650,

0.0100651,

])
fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.title("Intensity")
ax1.set_ylim(0.06,0.061)
ax1.set_ylabel("$I_{0,0}$")
ax1.set_xlabel("$N_{xy}$")
ax1.set_xlim(8,64)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.plot(Nxy,intensity_specular,label="Multiscat",color=[1., 0.5, 1.],linestyle=(0,(8,10)),marker='x',markersize=10)
plt.show()