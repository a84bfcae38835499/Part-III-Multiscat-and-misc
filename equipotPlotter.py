import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm
import datetime

values = np.array(pd.read_csv('Equipotential.csv'))[:,1:]
print("Values = ")
print(values)

X = np.linspace(0, 1, np.shape(values)[0])
Y = np.linspace(0, 1, np.shape(values)[0])
X, Y = np.meshgrid(X, Y)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.plot_surface(X, Y, values, vmin=-2,vmax=6)

ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[])

plt.show()

hm = sns.heatmap(data = values)
savestr = "Figures/" + datetime.datetime.now().strftime('Potential_%Y-%m-%d_%H-%M') + ".png"
plt.savefig(fname=savestr)
plt.show()