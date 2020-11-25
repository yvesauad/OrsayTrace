import trace as ot
import numpy
import matplotlib.pyplot as plt

z_array = numpy.linspace(-4.5, 4.5, 17)
res = 0.05

focus = 2.0
zlens = -2.0
r = 0.25
na = 0.2

# This example simple shows how to apply rotation to itens

a = ot.Simu(5, 5, 10, res)

a.d2_source(r, [0.5, 0.0, -4.5], [-na, 0, 1], 0.0, 1)
a.create_thin_lens([0, 0, zlens], focus, 1.5, 1.43, 'convex-plane')
a.rotate(numpy.arcsin(-na), [0, 1, 0], [0, 0, zlens])

#This next lens will not be rotated as it is instantiated after our rotation method
a.create_thin_lens([1.0, 0, zlens+4.0], focus, 1.5, 1.43, 'convex-plane')

for z in z_array:
    a.create_analysis_plan([-na, 0, 1], z)

photon_lists = a.run()
a.show_elements(photon_lists, 'all')
