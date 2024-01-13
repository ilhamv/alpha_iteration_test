import matplotlib.pyplot as plt
import h5py
import numpy as np

# Load MC and analytical results
with h5py.File("result.h5", "r") as f:
    leak_list = f["leak_list"][:]
    leak_list_mc = f["leak_list_mc"][:]
    alpha0 = f["alpha0"][:]
    k = f["k"][:]
    tr = f["tr"][:]
    alpha_mc = f["alpha0_mc"][:]
    k_mc = f["k_mc"][:]
    tr_mc = f["tr_mc"][:]
    E = f["E"][:]
    dE = E[1:] - E[:-1]
    E_mid = 0.5 * (E[1:] + E[:-1])

plt.plot(leak_list, alpha0, "b-", label="Ref.")
plt.plot(leak_list_mc, alpha_mc, "ro", fillstyle="none", label="OpenMC")
plt.legend()
plt.yscale("symlog")
plt.xlabel("Leakage XS [/cm]")
plt.ylabel(r"$\alpha$ [/s]")
plt.grid()
plt.show()

plt.plot(leak_list, tr, "b-", label="Ref.")
plt.plot(leak_list_mc, tr_mc, "ro", fillstyle="none", label="OpenMC")
plt.legend()
plt.xlabel("Leakage XS [/cm]")
plt.ylabel("Removal time [s]")
plt.grid()
plt.show()

plt.plot(leak_list, k, "b-", label="Ref.")
plt.plot(leak_list_mc, k_mc, "ro", fillstyle="none", label="OpenMC")
plt.legend()
plt.xlabel("Leakage XS [/cm]")
plt.ylabel("Multiplication factor")
plt.grid()
plt.show()
