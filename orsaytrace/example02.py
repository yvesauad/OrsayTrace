import trace as ot
import numpy
import matplotlib.pyplot as plt

z_array = numpy.linspace(-4.5, 14.5, 31)
f = 2.0
lens_pos = -2.5
res = 0.05
r = 0.25
na=0.0
ang = 1
thick = 0.3

f_pos = -1.29

#round part is at -2.5 - f/2. = -3.5. Plane is at -3.5+thick = -3.2. Focus is at -1.29, which is at 1.91 from plane side.


### 1st PART ###

a = ot.Simu(5, 5, 30, res)

#Plane convex lens. Source is point source diverging.
#a.d2_source(r, [-2*r, 0, -5.0], [0, 0, 1], na, ang)
#a.d2_source(r, [-2*r, 2*r, -5.0], [0, 0, 1], na, ang)
#a.d2_source(r, [-2*r, -2*r, -5.0], [0, 0, 1], na, ang)
a.d2_source(r, [0, 0, -5.0], [0, 0, 1], na, ang)
#a.d2_source(r, [0, 2*r, -5.0], [0, 0, 1], na, ang)
#a.d2_source(r, [0, -2*r, -5.0], [0, 0, 1], na, ang)
#a.d2_source(r, [2*r, 0, -5.0], [0, 0, 1], na, ang)
#a.d2_source(r, [2*r, 2*r, -5.0], [0, 0, 1], na, ang)
#a.d2_source(r, [2*r, -2*r, -5.0], [0, 0, 1], na, ang)

lens_pos02 = lens_pos+2*f-0.5

#a.create_sphere_element([0.0, 0.0, lens_pos], f/2., 1.43)
#a.create_sphere_element([0.0, 0.0, lens_pos02], f/2., 1.43)

#a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos-(f/2.-thick), lens_pos02+(f/2.-thick)], 1.0, [0, 0, 1], inclusive=False)
#a.create_cylinder_element([0, 0, lens_pos02+(f/2.-thick)], f/2., 0.0, 1.43, [0, 0, 1])

a.create_thin_lens([0, 0, 0.0], 2.0, 1.5, 1.43, 'plane-convex')

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
    vals_x = numpy.append(vals_x, photon_list.avg_position_axis(0))
    vals_y = numpy.append(vals_y, photon_list.avg_position_axis(1))


fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False)
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

