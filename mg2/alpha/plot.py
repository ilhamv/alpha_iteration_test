import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import h5py
import numpy as np

# Load MC and analytical results
with h5py.File("result.h5", "r") as f:
    delta_list = f["delta_list"][()]
    delta_list_mc = f["delta_list_mc"][()]
    alpha0 = f["alpha0"][()]
    rat12 = f["rat12"][()]
    alpha_mc = f["alpha_mc"][()]
    rat12_mc = f["rat12_mc"][()]

# Plot alpha
ax1 = plt.subplot()
ln1a = ax1.plot(delta_list, alpha0, "b-", label=r"$\alpha$")
ln1b = ax1.plot(
    delta_list_mc, alpha_mc, "bo", fillstyle="none", label=r"$\alpha$ - OpenMC"
)
ax1.set_xlabel(r"$\delta$")
ax1.set_ylabel(r"$\alpha$", color="b")
ax1.tick_params("y", colors="b")
ax1.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
plt.grid()

# Plot flux ratio
ax2 = ax1.twinx()
ln2a = ax2.plot(delta_list, rat12, "r--", label=r"$\phi_2/\phi_1$")
ln2b = ax2.plot(
    delta_list_mc, rat12_mc, "rD", fillstyle="none", label=r"$\phi_2/\phi_1$ - OpenMC"
)
ax2.set_ylabel(r"$\phi_2/\phi_1$", color="r")
ax2.tick_params("y", colors="r")
ln = ln1a + ln1b + ln2a + ln2b
lab = [l.get_label() for l in ln]
plt.legend(ln, lab, loc=2)

plt.show()
