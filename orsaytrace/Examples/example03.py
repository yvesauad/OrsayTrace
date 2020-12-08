import orsaytrace.trace as ot
import numpy
import matplotlib.pyplot as plt

z_array = numpy.linspace(-4.5, 4.5, 51)
res = 0.1

focus = 2.0
zlens = -2.0
r = 0.25
na = 0.2

# This example simple shows how to apply rotation to itens

a = ot.Simu(5, 5, 10, res)

a.d2_source(r, [0.5, 0.0, -4.5], [-na, 0, 1], 0.0, 1)
a.create_thin_lens([0, 0, zlens], focus, 1.5, 1.43, 'convex-plane')
a.rotate(numpy.arcsin(-na), [0, 1, 0], [0, 0, zlens])

for z in z_array:
    a.create_analysis_plan([-na, 0, 1], z, refraction_count=(1, 2))

photon_lists = a.run()
a.show_elements(photon_lists, 'all-noplan')
#a.show_elements(photon_lists, 'all')

div_z = numpy.asarray([])
div_tilted = numpy.asarray([])
for photon_list in photon_lists:
    div_z = numpy.append(div_z, photon_list.avg_divergence([0, 0, 1]))
    div_tilted = numpy.append(div_tilted, photon_list.avg_divergence([-na, 0, 1]))

fig, axes = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False)
axes[0].plot(z_array, div_z)
axes[1].plot(z_array, div_tilted)

axes[0].set_ylabel('Beam Divergence')
axes[1].set_ylabel('Beam Divergence')

axes[0].set_xlabel('Z (A.U.)')
axes[1].set_xlabel('Z (A.U.)')

plt.show()
