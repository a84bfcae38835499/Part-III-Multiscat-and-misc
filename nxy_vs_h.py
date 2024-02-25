from matplotlib import pyplot as plt
import numpy as np
n1 = range(-3,0)
intensities = np.array([0.360322E-01,0.520823E-01,0.138636E-01,0.811935E-01])

#     -3     0          0.360322E-01
#     -2     0          0.520823E-01
#     -1     0          0.138636E-01
#      0     0          0.811935E-01

intensities_n = intensities/max(intensities)
print(intensities_n)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(intensities_n,n1)
plt.title("Intensity (normalised)")
plt.ylim(0,1)
plt.xlabel("n1")
plt.xlim(-3,0)
plt.ylabel("Entropy")
plt.show()