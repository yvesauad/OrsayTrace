import trace as ot
import numpy
import matplotlib.pyplot as plt
import concurrent.futures
import time

x, y, z, res = 5, 5, 10, 0.1


focus = 0.3
yvertex = 0.5+0.3
thickness = 1.2
p = 2.0
r = 0.2

xmax, xmin, zmax, zmin = 0.45, -0.45, -0.25, -2.0

y_array = numpy.linspace(0, y/4, 101)

if __name__ == "__main__":
    with concurrent.futures.ProcessPoolExecutor(max_workers = 4) as executor:

        a = ot.Simu(x, y, z, res)

        a.point_source(r, [0, 0, -z/2], [0, 0, 1])

        a.create_parabolic_section_element([0.0, yvertex, 0.0], -1.0, 2*thickness, 3.0, p) 
        #a.create_rectangle_element([-x/2, x/2, yvertex-focus, y/2, -z/2, z/2], 1.0, [0, 0, 0]) 
        #a.rotate(numpy.pi/32, [0, 1, 0], [0, yvertex, 0.0], [-1.5, 1.5, -1.5, 1.5, -1.5, 1.5])

        for y in y_array:
            a.create_analysis_plan([0, 1, 0], y, reflection_count = 1)

        #a.show_created_elements('all-noplan')
        
        future_values = {executor.submit(a.run, i, 100, False, False): i for i in [0, 1]}
        #future1 = executor.submit(a.run, 1, 4, False, False)
        #future2 = executor.submit(a.run, 2, 4, False, False)
        #future3 = executor.submit(a.run, 3, 4, False, False)
        
        for future in concurrent.futures.as_completed(future_values):
            print(len(future.result()))


        a.show_elements(future.result(), 'all-noplan')
        #a.show_elements(future1.result(), 'all-noplan')



"""pd = numpy.asarray([photon_list.get_average_weighted_inverse() for photon_list in photon_lists])
pdvertex = numpy.asarray([photon_list.get_weighted_inverse_axis_y([0, -p/2.])for photon_list in photon_lists])
avg = numpy.asarray([photon_list.avg_position() for photon_list in photon_lists])
pos_max = numpy.asarray([photon_list.max_position() for photon_list in photon_lists])
pos_min = numpy.asarray([photon_list.min_position() for photon_list in photon_lists])
std = numpy.asarray([photon_list.std_position() for photon_list in photon_lists])


list_number = (numpy.where(std[:, 0]==min(std[:, 0])))[0][0]
fac = 10

fig, axes = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False)
axes[0, 0].plot(y_array, avg[:, 0], label='avg(X)')
axes[0, 0].scatter(y_array[list_number], avg[list_number, 0])
axes[0, 0].plot(y_array, avg[:, 2], label='avg(Z)')
axes[0, 0].scatter(y_array[list_number], avg[list_number, 2])

axes[0, 0].scatter(y_array[list_number-fac], avg[list_number-fac, 0], c='black')
axes[0, 0].scatter(y_array[list_number+fac], avg[list_number+fac, 0], c='black')
axes[0, 0].scatter(y_array[list_number-fac], avg[list_number-fac, 2], c='black')
axes[0, 0].scatter(y_array[list_number+fac], avg[list_number+fac, 2], c='black')

axes[0, 0].set_xlabel('Y')
axes[0, 0].legend()

axes[0, 1].plot(y_array, std[:, 0], label='std(X)')
axes[0, 1].scatter(y_array[list_number], std[list_number, 0])
axes[0, 1].plot(y_array, std[:, 2], label = 'std(Z)')
axes[0, 1].scatter(y_array[list_number], std[list_number, 2])
axes[0, 1].set_xlabel('Y (A.U.)')
axes[0, 1].legend()

axes[0, 2].plot(y_array, pd, label='Pd Centroid')
axes[0, 2].plot(y_array, pdvertex, label='Pd Vertex', c='red')
axes[0, 2].set_xlabel('Y')
axes[0, 2].legend()

#xmax = numpy.amax(pos_max[:, 0])
#xmin = numpy.amin(pos_min[:, 0])
#zmax = numpy.amax(pos_max[:, 2])
#zmin = numpy.amin(pos_min[:, 2])


#axes[0, 2].scatter(list_number, pd[list_number], label='Pd Centroid')
#axes[0, 2].scatter(list_number, pdvertex[list_number], label='Pd Vertex')
lists_pos = numpy.asarray([photon_list.get_positions() for photon_list in photon_lists])


axes[1, 0].hist2d(lists_pos[list_number-fac][:, 0], lists_pos[list_number-fac][:, 2], 100, range=[[xmin, xmax], [zmin, zmax]])
axes[1, 1].hist2d(lists_pos[list_number][:, 0], lists_pos[list_number][:, 2], 100, range=[[xmin, xmax], [zmin, zmax]])
axes[1, 2].hist2d(lists_pos[list_number+fac][:, 0], lists_pos[list_number+fac][:, 2], 100, range=[[xmin, xmax], [zmin, zmax]])

axes[1, 0].scatter([0], [-p/2], c='red', alpha=0.3)
axes[1, 1].scatter([0], [-p/2], c='red', alpha=0.3)
axes[1, 2].scatter([0], [-p/2], c='red', alpha=0.3)

axes[1, 0].set_xlabel('X'); axes[1, 0].set_ylabel('Z')
axes[1, 1].set_xlabel('X'); axes[1, 1].set_ylabel('Z')
axes[1, 2].set_xlabel('X'); axes[1, 2].set_ylabel('Z')

plt.show()"""
