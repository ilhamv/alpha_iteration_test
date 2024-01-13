import openmc
import numpy as np
import h5py
import matplotlib.pyplot as plt

# ===============================================================================
# Data
# ===============================================================================

# Load material data
with np.load("../SHEM-361.npz") as data:
    SigmaT = data["SigmaT"]
    SigmaC = data["SigmaC"]
    SigmaS = data["SigmaS"]
    nuSigmaF_p = data["nuSigmaF_p"]
    SigmaF = data["SigmaF"]
    nu_p = data["nu_p"]
    nu_d = data["nu_d"]
    chi_p = data["chi_p"]
    chi_d = data["chi_d"]
    G = data["G"]
    J = data["J"]
    E = data["E"]
    v = data["v"]
    lamd = data["lamd"]
G = G.item()
J = J.item()
SigmaA = SigmaT - np.sum(SigmaS, axis=0)
chi_p = np.sum(chi_p / np.sum(chi_p), axis=1)
nuSigmaF_p = np.outer(chi_p, SigmaF * nu_p)

# ===============================================================================
# Analytical
# ===============================================================================

N = 100
leak_list = np.linspace(0, 0.005, N)
results_alpha = np.zeros(N)
results_phi = np.zeros([N, G])
results_tr = np.zeros(N)
results_k = np.zeros(N)

for n in range(N):
    leak = leak_list[n]

    # The matrix
    A = np.zeros([G + J, G + J])

    # Top-left [GxG]: phi --> phi
    A[:G, :G] = SigmaS + nuSigmaF_p - np.diag(SigmaT) - np.diag(np.ones(G) * leak)

    # Top-right [GxJ]: C --> phi
    A[:G, G:] = np.multiply(chi_d, lamd)

    # Bottom-left [JxG]: phi --> C
    A[G:, :G] = np.multiply(nu_d, SigmaF)

    # bottom-right [JxJ]: C --> C
    A[G:, G:] = -np.diag(lamd)

    # Multiply with neutron speed
    M = np.copy(A)
    M[:G, :] = np.dot(np.diag(v), A[:G, :])

    M = M[:G, :G]

    alpha, phi = np.linalg.eig(M)
    idx = alpha.argsort()[::-1]
    alpha = alpha[idx][0].real
    phi = phi[:, idx]
    phi = phi[:, 0][:G].real
    phi /= np.sum(phi)

    tr = np.dot(phi, 1.0 / v) / np.dot(phi, SigmaA + leak)
    k = np.dot(phi, nu_p * SigmaF) / np.dot(phi, SigmaA + leak)

    results_alpha[n] = alpha
    results_phi[n, :] = phi
    results_tr[n] = tr
    results_k[n] = k

# ===============================================================================
# Monte Carlo
# ===============================================================================

N = 20
leak_list_mc = np.linspace(0, 0.005, N)
results_alpha_mc = np.zeros(N)
results_phi_mc = np.zeros([N, G])
results_tr_mc = np.zeros(N)
results_k_mc = np.zeros(N)

SigmaS = SigmaS.transpose()
SigmaS = SigmaS[:, :, np.newaxis]

for n in range(N):
    leak = leak_list_mc[n]

    # ===========================================================================
    # Set Library
    # ===========================================================================

    groups = openmc.mgxs.EnergyGroups(openmc.mgxs.GROUP_STRUCTURES["SHEM-361"])

    uo2_xsdata = openmc.XSdata("uo2", groups, num_delayed_groups=J)
    uo2_xsdata.order = 0

    uo2_xsdata.set_inverse_velocity(1 / v, temperature=294.0)

    uo2_xsdata.set_total(SigmaT + leak, temperature=294.0)
    uo2_xsdata.set_absorption(SigmaA + leak, temperature=294.0)

    uo2_xsdata.set_scatter_matrix(SigmaS, temperature=294.0)
    uo2_xsdata.set_decay_rate(lamd, temperature=294.0)

    uo2_xsdata.set_prompt_nu_fission(nu_p * SigmaF, temperature=294.0)
    uo2_xsdata.set_delayed_nu_fission(
        np.multiply(nu_d, SigmaF),
        temperature=294.0,
    )
    uo2_xsdata.set_chi_prompt(chi_p, temperature=294.0)
    uo2_xsdata.set_chi_delayed(chi_d.transpose(), temperature=294.0)
    mg_cross_sections_file = openmc.MGXSLibrary(groups, J)
    mg_cross_sections_file.add_xsdata(uo2_xsdata)
    mg_cross_sections_file.export_to_hdf5("mgxs.h5")

    # ===========================================================================
    # Simulation Input File Parameters
    # ===========================================================================

    # OpenMC simulation parameters
    batches = 50
    inactive = 10
    particles = int(1e6)

    # ===========================================================================
    # Exporting to OpenMC materials.xml file
    # ===========================================================================

    materials = {}
    materials["uo2"] = openmc.Material(name="uo2")
    materials["uo2"].set_density("macro", 1.0)
    materials["uo2"].add_macroscopic("uo2")
    materials_file = openmc.Materials(materials.values())
    materials_file.cross_sections = "mgxs.h5"
    materials_file.export_to_xml()

    # ===========================================================================
    # Exporting to OpenMC geometry.xml file
    # ===========================================================================

    # Instantiate ZCylinder surfaces
    surf_Z1 = openmc.ZPlane(surface_id=1, z0=0.0, boundary_type="reflective")
    surf_Z2 = openmc.ZPlane(surface_id=2, z0=5.0, boundary_type="reflective")

    # Instantiate Cells
    cell_F = openmc.Cell(cell_id=1, name="F")

    # Use surface half-spaces to define regions
    cell_F.region = +surf_Z1 & -surf_Z2

    # Register Materials with Cells
    cell_F.fill = materials["uo2"]

    # Instantiate Universes
    root = openmc.Universe(universe_id=0, name="root universe", cells=[cell_F])

    # Instantiate a Geometry, register the root Universe, and export to XML
    geometry = openmc.Geometry(root)
    geometry.export_to_xml()

    # ===========================================================================
    # Exporting to OpenMC settings.xml file
    # ===========================================================================

    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings_file = openmc.Settings()
    settings_file.batches = batches
    settings_file.inactive = inactive
    settings_file.particles = particles
    settings_file.output = {"tallies": False}
    settings_file.alpha_mode = True
    settings_file.prompt_only = True
    settings_file.energy_mode = "multi-group"

    # Create an initial uniform spatial source distribution over fissionable zones
    bounds = [0.0, 0.0, 0.0, 100.0, 100.0, 5.0]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings_file.source = openmc.source.Source(space=uniform_dist)

    settings_file.export_to_xml()

    # ===========================================================================
    # Exporting to OpenMC tallies.xml file
    # ===========================================================================

    energy_filter = openmc.EnergyFilter(
        openmc.mgxs.GROUP_STRUCTURES["SHEM-361"], filter_id=1
    )

    # Instantiate the first Tally
    first_tally = openmc.Tally(tally_id=1, name="first tally")
    first_tally.scores = ["flux"]
    first_tally.filters = [energy_filter]

    # Instantiate a Tallies collection and export to XML
    tallies_file = openmc.Tallies([first_tally])
    tallies_file.export_to_xml()

    openmc.run(mpi_args=["srun", "-n", "8"])
    # openmc.run(threads=1)

    with openmc.StatePoint("statepoint.50.h5") as osp:
        alpha, err = osp.alpha_eff
        basic_tally = osp.get_tally(name="first tally")
        flux = basic_tally.get_values(scores=["flux"]).reshape(G)
        results_alpha_mc[n] = alpha
        results_phi_mc[n, :] = flux

    with h5py.File("statepoint.50.h5", "r") as f:
        tr = f["alpha_mode_tallies/removal_time"][0]
        k = f["alpha_mode_tallies/k_effective"][0]
        results_tr_mc[n] = tr
        results_k_mc[n] = k


# ===============================================================================
# Plot
# ===============================================================================

with h5py.File("result.h5", "w") as hdf:
    hdf.create_dataset("leak_list", data=leak_list)
    hdf.create_dataset("leak_list_mc", data=leak_list_mc)
    hdf.create_dataset("alpha0", data=results_alpha)
    hdf.create_dataset("phi", data=results_phi)
    hdf.create_dataset("tr", data=results_tr)
    hdf.create_dataset("k", data=results_k)
    hdf.create_dataset("alpha0_mc", data=results_alpha_mc)
    hdf.create_dataset("phi_mc", data=results_phi_mc)
    hdf.create_dataset("tr_mc", data=results_tr_mc)
    hdf.create_dataset("k_mc", data=results_k_mc)
    hdf.create_dataset("E", data=E)
