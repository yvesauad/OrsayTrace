import trace as ot
import numpy
import matplotlib.pyplot as plt

x, y, z, res = 5, 5, 10, 0.02


focus = 0.3
yvertex = 0.5+0.3
thickness = 1.2
p = 2.0
r = 0.25

xmax, xmin, zmax, zmin = 0.45, -0.45, -0.25, -2.0

y_array = numpy.linspace(0, y/4, 101)

a = ot.Simu(x, y, z, res)

a.point_source(r, [0, 0, -z/2], [0, 0, 1])

a.create_parabolic_section_element([0.0, yvertex, 0.0], -1.0, 2*thickness, 3.0, p) 
a.create_rectangle_element([-x/2, x/2, yvertex-focus, y/2, -z/2, z/2], 1.0, [0, 0, 0]) 
#a.rotate(numpy.pi/32, [1, 0, 0], [0, yvertex, 0.0])


for y in y_array:
    a.create_analysis_plan([0, 1, 0], y, reflection_count = 1)

#a.show_created_elements('all-noplan')
photon_lists = a.run()
#a.show_elements(photon_lists, 'all-noplan')


pd = numpy.asarray([photon_list.get_average_weighted_inverse() for photon_list in photon_lists])
avg = numpy.asarray([photon_list.avg_position() for photon_list in photon_lists])
pos_max = numpy.asarray([photon_list.max_position() for photon_list in photon_lists])
pos_min = numpy.asarray([photon_list.min_position() for photon_list in photon_lists])
std = numpy.asarray([photon_list.std_position() for photon_list in photon_lists])


fig, axes = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False)
axes[0, 0].plot(avg[:, 0], label='avg(X)')
axes[0, 0].plot(avg[:, 2], label='avg(Z)')
axes[0, 0].set_xlabel('Y (A.U.)')
axes[0, 0].legend()

axes[0, 1].plot(std[:, 0], label='std(X)')
axes[0, 1].plot(std[:, 2], label = 'std(Z)')
axes[0, 1].set_xlabel('Y (A.U.)')
axes[0, 1].legend()

axes[0, 2].plot(pd, label='PD')
axes[0, 2].set_xlabel('Y (A.U.)')
axes[0, 2].legend()

#xmax = numpy.amax(pos_max[:, 0])
#xmin = numpy.amin(pos_min[:, 0])
#zmax = numpy.amax(pos_max[:, 2])
#zmin = numpy.amin(pos_min[:, 2])


list_number = (numpy.where(std[:, 0]==min(std[:, 0])))[0][0]
fac = 16
lists_pos = numpy.asarray([photon_list.get_positions() for photon_list in photon_lists])


axes[1, 0].hist2d(lists_pos[list_number-fac][:, 0], lists_pos[list_number-fac][:, 2], 100, range=[[xmin, xmax], [zmin, zmax]])
axes[1, 1].hist2d(lists_pos[list_number][:, 0], lists_pos[list_number][:, 2], 100, range=[[xmin, xmax], [zmin, zmax]])
axes[1, 2].hist2d(lists_pos[list_number+fac][:, 0], lists_pos[list_number+fac][:, 2], 100, range=[[xmin, xmax], [zmin, zmax]])

axes[1, 0].scatter([0], [-p/2], c='yellow', alpha=0.1)
axes[1, 1].scatter([0], [-p/2], c='yellow', alpha=0.1)
axes[1, 2].scatter([0], [-p/2], c='yellow', alpha=0.1)

axes[1, 0].set_xlabel('X'); axes[1, 0].set_ylabel('Z')
axes[1, 1].set_xlabel('X'); axes[1, 1].set_ylabel('Z')
axes[1, 2].set_xlabel('X'); axes[1, 2].set_ylabel('Z')

plt.show()
