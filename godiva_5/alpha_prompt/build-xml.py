import openmc


# ===============================================================================
# Simulation Input File Parameters
# ===============================================================================

# OpenMC simulation parameters
batches = 450
inactive = 150
particles = int(1e6)


# ===============================================================================
# Exporting to OpenMC materials.xml file
# ===============================================================================

# Instantiate some Materials and register the appropriate Nuclides
fuel = openmc.Material(material_id=1, name="fuel")
fuel.set_density("g/cc", 9.3999)
fuel.add_nuclide("U235", 0.937695)
fuel.add_nuclide("U238", 0.052053)
fuel.add_nuclide("U234", 0.010252)

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([fuel])
materials_file.export_to_xml()


# ===============================================================================
# Exporting to OpenMC geometry.xml file
# ===============================================================================

# Instantiate ZCylinder surfaces
surf1 = openmc.Sphere(surface_id=1, r=8.7407, boundary_type="vacuum")

# Instantiate Cells
godiva = openmc.Cell(cell_id=1, name="godiva")

# Use surface half-spaces to define regions
godiva.region = -surf1

# Register Materials with Cells
godiva.fill = fuel

# Instantiate Universes
root = openmc.Universe(universe_id=0, name="root universe", cells=[godiva])

# Instantiate a Geometry, register the root Universe, and export to XML
geometry = openmc.Geometry(root)
geometry.export_to_xml()


# ===============================================================================
# Exporting to OpenMC settings.xml file
# ===============================================================================

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.output = {"tallies": False}
settings_file.alpha_mode = True
settings_file.prompt_only = True

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-9.0, -9.0, -9.0, 9.0, 9.0, 9.0]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = (-9, -9, -9)
entropy_mesh.upper_right = (9, 9, 9)
entropy_mesh.dimension = (10, 10, 10)

settings_file.entropy_mesh = entropy_mesh

settings_file.export_to_xml()


# ===============================================================================
# Exporting to OpenMC tallies.xml file
# ===============================================================================

# Instantiate some tally Filters
energy_filter = openmc.EnergyFilter(openmc.mgxs.GROUP_STRUCTURES["UKAEA-1102"])

# Instantiate the first Tally
first_tally = openmc.Tally(tally_id=1, name="first tally")
first_tally.scores = ["nu-fission"]
first_tally.filters = [energy_filter]

# Instantiate a Tallies collection and export to XML
tallies_file = openmc.Tallies([first_tally])
tallies_file.export_to_xml()
