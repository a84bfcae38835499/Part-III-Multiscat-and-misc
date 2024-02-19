from matplotlib import pyplot as plt

Nxys = [64, 32, 16, 8, 4]
Entropies = [0.850912701202775, 0.850912701202775,0.850912701202775, 0.850912701202775, 0.850912701202775]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(Nxys,Entropies,linestyle="-",marker='o',c=[0,0,0])
plt.title("Entropy vs Nxy")
plt.xlim(0,64)
plt.xlabel("Nxy")
plt.ylim(0,1)
plt.ylabel("Entropy")
plt.show()