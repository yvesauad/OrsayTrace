import orsaytrace.trace as ot
import numpy
import matplotlib.pyplot as plt

x, y, z, res = 5, 5, 5, 0.008

a = ot.Simu(x, y, z, res)

a.create_chessboard(0.5, [0, 0, -2.5], 10)

a.create_parabolic_surface_element([0.0, 0, 2.0], -1.0, 2, 3.0, 1.0)
#a.create_rectangle_element([-x/2, x/2, 0-0.3, y/2, -z/2, z/2], 1.0, [0, 0, 0])

a.create_analysis_plan([0, 0, 1], 0.5, reflection_count = (0, 0))
a.create_analysis_plan([0, 0, 1], 1.3, reflection_count = (1, 1))
a.create_analysis_plan([0, 1, 0], 0.2, reflection_count = (1, 1))
a.create_analysis_plan([0, 1, 0], -0.2, reflection_count = (1, 1))

#a.show_created_elements('all')
pl = a.run()
#a.show_elements(pl, 'photons-verbose')
a.show_photons2D(pl, mode='-verbose', binning=100)

