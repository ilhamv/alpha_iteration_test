import openmc
import numpy as np
import h5py

# ===============================================================================
# Data
# ===============================================================================

G = 2
J = 2
N = G + J

v1 = 10.0
v2 = 5.0

SigmaS11 = 1 / 2
SigmaS12 = 1 / 2
SigmaS21 = 0
SigmaS22 = 1

SigmaF1 = 0
SigmaF2 = 1

SigmaT1 = 2
SigmaT2 = 3

beta1 = 1 / 4
beta2 = 1 / 8
beta = beta1 + beta2

chiP1 = 1.0
chiP2 = 0.0

chiD11 = 3 / 4
chiD12 = 1 / 4

chiD21 = 1 / 2
chiD22 = 1 / 2

lam1 = 0.5
lam2 = 0.07


# ===============================================================================
# Analytical
# ===============================================================================

nu_list = np.linspace(1, 15, 500)
alpha0 = []
rat12 = []

for nu in nu_list:
    M = np.zeros([N, N])

    M[0, 0] = v1 * (SigmaS11 + chiP1 * (1 - beta) * nu * SigmaF1 - SigmaT1)
    M[0, 1] = v1 * (SigmaS21 + chiP1 * (1 - beta) * nu * SigmaF2)
    M[1, 0] = v2 * (SigmaS12 + chiP2 * (1 - beta) * nu * SigmaF1)
    M[1, 1] = v2 * (SigmaS22 + chiP2 * (1 - beta) * nu * SigmaF2 - SigmaT2)

    M[0, G + 0] = v1 * chiD11 * lam1
    M[0, G + 1] = v1 * chiD21 * lam2
    M[1, G + 0] = v2 * chiD12 * lam1
    M[1, G + 1] = v2 * chiD22 * lam2

    M[G + 0, 0] = beta1 * nu * SigmaF1
    M[G + 0, 1] = beta1 * nu * SigmaF2
    M[G + 1, 0] = beta2 * nu * SigmaF1
    M[G + 1, 1] = beta2 * nu * SigmaF2

    M[G + 0, G + 0] = -lam1
    M[G + 1, G + 1] = -lam2

    M = M[:2, :2]

    alpha, phi = np.linalg.eig(M)
    idx = alpha.argsort()[::-1]
    alpha = alpha[idx]
    phi = phi[:, idx]

    alpha0.append(alpha[0])
    phi = phi[:, 0]
    phi21 = phi[1] / phi[0]
    rat12.append(phi21)


# ===============================================================================
# Monte Carlo
# ===============================================================================

nu_list_mc = np.linspace(1, 15, 10)
nu_list_mc = np.append(nu_list_mc, 24 / 5)
alpha_mc = []
rat12_mc = []

for nu in nu_list_mc:
    # ===========================================================================
    # Set Library
    # ===========================================================================

    groups = openmc.mgxs.EnergyGroups([0.0, 1e-5, 2e7])

    uo2_xsdata = openmc.XSdata("uo2", groups, num_delayed_groups=J)
    uo2_xsdata.order = 0

    uo2_xsdata.set_inverse_velocity([1 / v1, 1 / v2], temperature=294.0)

    uo2_xsdata.set_total([SigmaT1, SigmaT2], temperature=294.0)
    uo2_xsdata.set_absorption([SigmaT1 - 1, SigmaT2 - 1], temperature=294.0)
    scatter_matrix = [[[SigmaS11], [SigmaS12]], [[SigmaS21], [SigmaS22]]]
    uo2_xsdata.set_scatter_matrix(scatter_matrix, temperature=294.0)
    uo2_xsdata.set_decay_rate([lam1, lam2], temperature=294.0)

    nu_p = (1.0 - beta) * nu
    nu_d1 = beta1 * nu
    nu_d2 = beta2 * nu
    uo2_xsdata.set_prompt_nu_fission(
        [nu_p * SigmaF1, nu_p * SigmaF2], temperature=294.0
    )
    uo2_xsdata.set_delayed_nu_fission(
        [[nu_d1 * SigmaF1, nu_d1 * SigmaF2], [nu_d2 * SigmaF1, nu_d2 * SigmaF2]],
        temperature=294.0,
    )
    uo2_xsdata.set_chi_prompt([chiP1, chiP2], temperature=294.0)
    uo2_xsdata.set_chi_delayed([[chiD11, chiD12], [chiD21, chiD22]], temperature=294.0)
    mg_cross_sections_file = openmc.MGXSLibrary(groups, 2)
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

    energy_filter = openmc.EnergyFilter([0.0, 1e-5, 2e7], filter_id=1)

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
        rat = flux[0] / flux[1]
        alpha_mc.append(alpha)
        rat12_mc.append(rat)


# ===============================================================================
# Plot
# ===============================================================================

with h5py.File("result.h5", "w") as hdf:
    hdf.create_dataset("nu_list", data=nu_list)
    hdf.create_dataset("nu_list_mc", data=nu_list_mc)
    hdf.create_dataset("alpha0", data=alpha0)
    hdf.create_dataset("rat12", data=rat12)
    hdf.create_dataset("alpha_mc", data=alpha_mc)
    hdf.create_dataset("rat12_mc", data=rat12_mc)
