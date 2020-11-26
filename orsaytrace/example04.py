import trace as ot
import numpy
import matplotlib.pyplot as plt

x, y, z, res = 5, 5, 10, 0.06


focus = 0.3
yvertex = 0.5+0.3
thickness = 1.2
p = 2.0
r = 0.25

y_array = numpy.linspace(-y/2, y/2, 11)
z_array = numpy.linspace(-x/2., +x/2., 11)

a = ot.Simu(x, y, z, res)

a.point_source(r, [0, 0, -z/2], [0, 0, 1])

a.create_parabolic_section_element([0.0, yvertex, 0.0], -1.0, 2*thickness, 3.0, p) 
a.create_rectangle_element([-x/2, x/2, yvertex-focus, y/2, -z/2, z/2], 1.0, [0, 0, 0]) 
#a.rotate(-numpy.pi/16, [0, 1, 0], [0, yvertex, 0.0])


for z in z_array:
    a.create_analysis_plan([0, 0, 1], z)

#a.show_created_elements('all-noplan')
photon_lists = a.run()
a.show_elements(photon_lists, 'all-noplan')

vals = numpy.asarray([])
#vals_distance = numpy.asarray([])
#vals_x = numpy.asarray([])
#vals_z = numpy.asarray([])
for photon_list in photon_lists:
    vals = numpy.append(vals, photon_list.avg_divergence([0, 1, 0]))
#    vals_distance = numpy.append(vals_distance, photon_list.avg_distance_axis_z([0, 0]))
#    vals_x = numpy.append(vals_x, photon_list.avg_position_axis(0))
#    vals_z = numpy.append(vals_z, photon_list.avg_position_axis(2))


#fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False)
#axes[0, 0].plot(z_array, vals)
#axes[0, 1].plot(z_array, vals_distance)
#axes[1, 0].plot(z_array, vals_x)
#axes[1, 1].plot(z_array, vals_z)

#axes[0, 0].set_ylabel('Beam Divergence')
#axes[0, 1].set_ylabel('Distance from Optical Axis')
#axes[1, 0].set_ylabel('Average X')
#axes[1, 1].set_ylabel('Average Y')

#axes[0, 0].set_xlabel('Z (A.U.)')
#axes[0, 1].set_xlabel('Z (A.U.)')
#axes[1, 0].set_xlabel('Z (A.U.)')
#axes[1, 1].set_xlabel('Z (A.U.)')

#plt.show()
