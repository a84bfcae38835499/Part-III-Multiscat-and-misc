from matplotlib import pyplot as plt
import numpy as np

times = [3.35, 37.35, 408.47]
nch = [170, 671, 1515]
ncell = [1, 2, 3]
ncellcont = np.linspace(1,10,100)
expmodel = 3.35 * np.exp((ncellcont-1)*(np.log(2)/0.29))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(nch,times,color='red')
ax.tick_params(axis='x', labelcolor='red')
ax.set_xlabel("Number of diffraction channels")
ax2 = ax.twiny()
ax2.plot(ncell,times,color='green')
ax2.plot(ncellcont,expmodel,color='green')
ax2.tick_params(axis='x', labelcolor='green')
ax2.set_xlabel("Supercell size")
ax.set_ylabel("Time / s")

plt.show()