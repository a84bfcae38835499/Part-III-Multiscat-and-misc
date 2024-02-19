from matplotlib import pyplot as plt

DefectDensities = [0,0.25,0.5,0.75,1]
Entropies = [0.6495282977909468,
0.852275076163138,
0.7376804580552242,
0.7510843204920689,
0.4909806615832078]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(DefectDensities,Entropies,linestyle="-",marker='x',c=[0.1,0.5,1])
plt.title("Entropy vs defect density")
plt.xlim(0,1)
plt.xlabel("Fraction of sulphur atoms missing")
plt.ylim(0,1)
plt.ylabel("Entropy")
plt.show()