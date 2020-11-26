import trace as ot
import numpy

x, y, z, res = 5, 5, 10, 0.05


focus = 0.3
yvertex = 0.5+0.3
thickness = 1.2
p = 2.0
r = 0.5

y_array = numpy.linspace(r+0.1, yvertex, 5)

a = ot.Simu(x, y, z, res)

a.point_source(r, [0, 0, -z/2], [0, 0, 1])

a.create_parabolic_section_element([0.0, yvertex, 0.0], -1.0, 2*thickness, 3.0, p) 
a.create_rectangle_element([-x/2, x/2, yvertex-focus, y/2, -z/2, z/2], 1.0, [0, 0, 0]) 

for y in y_array:
    a.create_analysis_plan([0, 1, 0], y)

a.show_created_elements('all-noplan')
photon_lists = a.run()
a.show_elements(photon_lists, 'all')
