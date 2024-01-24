import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm
import datetime

f = open("Equipotential.csv", "r")
lines = f.readlines()

comment = ""
for line in lines:
    if(line[0] == '#'):
           comment += str(line[1:])
print(comment)
f.close()


values = np.transpose(np.array(pd.read_csv('Equipotential.csv',comment='#'))[:,1:])
#transpose to get x to match with x
print("Values = ")
print(values)


X = np.linspace(0, 1, np.shape(values)[0])
Y = np.linspace(0, 1, np.shape(values)[0])
X, Y = np.meshgrid(X, Y)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.plot_surface(X, Y, values, vmin=-2,vmax=6)
#ax.plot_surface(X, Y, values, vmin=-1,vmax=1)

ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[])

plt.show()

hm = sns.heatmap(data = values)
plt.gca().invert_yaxis()
plt.gca().set_aspect(1/np.sqrt(3))
savestr = "Figures/Potentials/" + datetime.datetime.now().strftime('Potential_%Y-%m-%d_%H-%M') + ".png"
hm.set_title(comment)
plt.savefig(fname=savestr)
plt.show()