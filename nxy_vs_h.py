from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
Nsuper = range(1,6)
times_ms = np.array([6, 283, 2671, 10791, 50317])
times_matlab = np.array([20, 50, 110, 193, 420])

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.title("Time taken")
ax1.set_ylim(6,51000)
ax1.set_ylabel("Time taken /s")
ax1.set_xlabel("Supercell size")
ax1.set_xlim(1,5)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.set_yscale('log')
ax1.plot(Nsuper,times_ms,label="Multiscat",color=[0., 0.5, 1.],linestyle=(0,(8,10)),marker='+',markersize=15)
ax1.plot(Nsuper,times_matlab,label="MATLAB",color=[1., 0.5, 0.],linestyle=(0,(8,10)),marker='x',markersize=10)
ax1.legend()
plt.show()