import orsaytrace.trace as ot
import numpy
import matplotlib.pyplot as plt

z_array = numpy.linspace(-9.0, 9.0, 101)
focus1 = focus2 = 2.0
zlens1 = -6.0
d12 = 6.0
res = 0.08
r_src = 0.25

a = ot.Simu(5, 5, 20, res)

a.d2_source(r_src, [0, 0, -9.5], [0, 0, 1], 0.0, 1)

a.create_thin_lens([0, 0, zlens1], focus1, 1.5, 1.43, 'convex-plane')
a.create_thin_lens([0, 0, zlens1+d12], focus2, 1.5, 1.43, 'plane-convex')

for z in z_array:
    a.create_analysis_plan([0, 0, 1], z)

a.show_created_elements('all-noplan')
photon_lists = a.run()
a.show_elements(photon_lists, 'all-noplan')

vals = numpy.asarray([])
vals_distance = numpy.asarray([])
vals_x = numpy.asarray([])
vals_y = numpy.asarray([])
for photon_list in photon_lists:
    vals = numpy.append(vals, photon_list.avg_divergence([0, 0, 1]))
    vals_distance = numpy.append(vals_distance, photon_list.avg_distance_axis_z([0, 0]))
    vals_x = numpy.append(vals_x, photon_list.avg_position()[0])
    vals_y = numpy.append(vals_y, photon_list.avg_position()[1])


fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, dpi=200)
axes[0, 0].plot(z_array, vals)
axes[0, 1].plot(z_array, vals_distance)
axes[1, 0].plot(z_array, vals_x)
axes[1, 1].plot(z_array, vals_y)

axes[0, 0].set_ylabel('Beam Divergence')
axes[0, 1].set_ylabel('Distance from Optical Axis')
axes[1, 0].set_ylabel('Average X')
axes[1, 1].set_ylabel('Average Y')

axes[0, 0].set_xlabel('Z (A.U.)')
axes[0, 1].set_xlabel('Z (A.U.)')
axes[1, 0].set_xlabel('Z (A.U.)')
axes[1, 1].set_xlabel('Z (A.U.)')

plt.show()
