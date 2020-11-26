import trace as ot
import numpy
import matplotlib.pyplot as plt

x, y, z, res = 5, 5, 10, 0.03


focus = 0.3
yvertex = 0.5+0.3
thickness = 1.2
p = 2.0
r = 0.25

y_array = numpy.linspace(-y/2, y/2, 11)
z_array = numpy.linspace(-z/2., z/2, 31)

a = ot.Simu(x, y, z, res)

a.point_source(r, [0, 0, -z/2], [0, 0, 1])

a.create_parabolic_section_element([0.0, yvertex, 0.0], -1.0, 2*thickness, 3.0, p) 
#a.create_rectangle_element([-x/2, x/2, yvertex-focus, y/2, -z/2, z/2], 1.0, [0, 0, 0]) 
#a.rotate(-numpy.pi/16, [0, 1, 0], [0, yvertex, 0.0])


for z in z_array:
    a.create_analysis_plan([0, 0, 1], z, reflection_count = (1, 1))

#a.show_created_elements('all-noplan')
photon_lists = a.run()
a.show_elements(photon_lists, 'all-noplan')

vals = numpy.asarray([])
for photon_list in photon_lists:
    vals = numpy.append(vals, photon_list.avg_divergence([0, 1, 0]))

fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False)
axes.plot(z_array, vals)
axes.set_ylabel('Beam Divergence')
axes.set_xlabel('Z (A.U.)')

#plt.show()
