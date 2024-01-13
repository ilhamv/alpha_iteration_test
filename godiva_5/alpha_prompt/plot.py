import h5py
import matplotlib.pyplot as plt

with h5py.File("statepoint.450.h5", "r") as f:
    alpha = f["alpha_mode_tallies/alpha_generation"][:]
    k = f["k_generation"][:]
plt.plot(alpha, "b")
ax = plt.gca().twinx()
ax.plot(k, "r")
plt.show()
