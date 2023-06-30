from stellarutil.simulation import Simulation, get_field, get_field_name
from stellarutil.graph import graph

sim = Simulation()
data = sim.get_field("n_star")
sim.get_stars_in_galaxy(0)
print(data)