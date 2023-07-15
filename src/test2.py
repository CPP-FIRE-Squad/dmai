from stellarutil.simulation import Star, Simulation


star = Star()

print(star.velocity())
print(star.vx)


sim = Simulation()
x_gal, y_gal, z_gal, a_gal, m_gal, v_gal = sim.get_stars_in_halo(0)



