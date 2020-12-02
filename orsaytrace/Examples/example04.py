import orsaytrace.trace as ot
import numpy
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

x, y, z, res = 5, 5, 10, 0.04

focus = 0.3
yvertex = 0.8
thickness = 1.2
p = 2.0
r = 0.2

xmax, xmin, zmax, zmin = 0.35, -0.35, -0.25, -2.0

y_array = numpy.linspace(0, y/4, 201)

a = ot.Simu(x, y, z, res)

a.d2_source(r, [0, 0, -z/2], [0, 0, 1], 0.0, 1)

a.create_parabolic_surface_element([0.0, yvertex, 0.0], -1.0, 2*thickness, 3.0, p) 
a.create_rectangle_element([-x/2, x/2, yvertex-focus, y/2, -z/2, z/2], 1.0, [0, 0, 0])
#a.rotate(numpy.pi/64, [0, 1, 0], [0, yvertex, 0.0], [-0.5, 0.5, -0.5, 0.5, -2.0, 0])

for y in y_array:
    a.create_analysis_plan([0, 1, 0], y, reflection_count = 1)

#a.show_created_elements('all-noplan')
photon_lists = a.run()
#a.show_elements(photon_lists, 'all-noplan')

pd = numpy.asarray([photon_list.get_average_weighted_inverse() for photon_list in photon_lists])
pdvertex = numpy.asarray([photon_list.get_average_weighted_inverse_axis_y([0, -p/2.])for photon_list in photon_lists])
avg = numpy.asarray([photon_list.avg_position() for photon_list in photon_lists])
pos_max = numpy.asarray([photon_list.max_position() for photon_list in photon_lists])
pos_min = numpy.asarray([photon_list.min_position() for photon_list in photon_lists])
std = numpy.asarray([photon_list.std_position() for photon_list in photon_lists])


list_number = (numpy.where(std[:, 0]==min(std[:, 0])))[0][0]
fac = 50

fig, axes = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=False, dpi=200)
axes[0].plot(y_array, avg[:, 0], label='avg(X)')
axes[0].scatter(y_array[list_number], avg[list_number, 0])
axes[0].plot(y_array, avg[:, 2], label='avg(Z)')
axes[0].scatter(y_array[list_number], avg[list_number, 2])

axes[0].scatter(y_array[list_number-fac], avg[list_number-fac, 0], c='black')
axes[0].scatter(y_array[list_number+fac], avg[list_number+fac, 0], c='black')
axes[0].scatter(y_array[list_number-fac], avg[list_number-fac, 2], c='black')
axes[0].scatter(y_array[list_number+fac], avg[list_number+fac, 2], c='black')

axes[0].set_xlabel('Y')
axes[0].legend()

axes[1].plot(y_array, std[:, 0], label='std(X)')
axes[1].scatter(y_array[list_number], std[list_number, 0])
axes[1].plot(y_array, std[:, 2], label = 'std(Z)')
axes[1].scatter(y_array[list_number], std[list_number, 2])
axes[1].set_xlabel('Y')
axes[1].legend()

axes[2].plot(y_array, pd, label='WI CoM')
axes[2].plot(y_array, pdvertex, label='WI FP', c='red')
axes[2].set_xlabel('Y')
axes[2].legend()
plt.show()

lists_pos = numpy.asarray([photon_list.get_positions() for photon_list in photon_lists])

fig, axes = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=False, dpi=200)
axes[0].hist2d(lists_pos[list_number-fac][:, 0], lists_pos[list_number-fac][:, 2], 100, range=[[xmin, xmax], [zmin, zmax]])
axes[1].hist2d(lists_pos[list_number][:, 0], lists_pos[list_number][:, 2], 100, range=[[xmin, xmax], [zmin, zmax]])
axes[2].hist2d(lists_pos[list_number+fac][:, 0], lists_pos[list_number+fac][:, 2], 100, range=[[xmin, xmax], [zmin, zmax]])

axes[0].scatter([0], [-p/2], c='red', alpha=0.3)
axes[1].scatter([0], [-p/2], c='red', alpha=0.3)
axes[2].scatter([0], [-p/2], c='red', alpha=0.3)

axes[0].set_xlabel('X'); axes[0].set_ylabel('Z')
axes[1].set_xlabel('X'); axes[1].set_ylabel('Z')
axes[2].set_xlabel('X'); axes[2].set_ylabel('Z')
plt.show()


fig, ax = plt.subplots()
plt.hist2d(lists_pos[0][:, 0], lists_pos[0][:, 2], 101, range=[[xmin, xmax], [zmin, zmax]])

def init():
    pass

def update(i):
    plt.clf()
    plt.hist2d(lists_pos[int(i)][:, 0], lists_pos[int(i)][:, 2], 101, range=[[xmin, xmax], [zmin, zmax]],
               cmap=plt.cm.jet)
    plt.text(0.18, -0.4, 'y = '+format(y_array[i], '.2f'), fontsize=12, color='white', weight='bold')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title('X-Z Plan')
ani = FuncAnimation(fig, update, range(len(y_array)), init_func=init)

writer = PillowWriter(fps=40)
ani.save("demo_parabolic.gif", writer=writer)
