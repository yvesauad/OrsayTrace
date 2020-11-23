import parabolic as ot
import numpy

a = ot.Simu(5, 5, 12, 0.05)

zs = -5.0
z1 = -2.0
f1 = 2.0
f2 = 2.0
d = 4.0

a.d2_source(0.5, -5.0, [0, 0, 1], 0.0, 1)
a.create_sphere_element([0.0, 0.0, z1], f1/2., 1.43) #from 2 to 3.
a.create_sphere_element([0.0, 0.0, z1+d], f2/2., 1.43) #from 2 to 3.
a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, -2.0, 2.0], 1.0, [0, 0, 1]) 



a.create_analysis_plan(z=-2.0)
a.create_analysis_plan(z=0)
a.create_analysis_plan(z=4.0)
a.create_analysis_plan(z=5.0)

a.show_elements(False, 'all')
a.run()
a.show_elements(True, 'photons')
a.show_photons2D('xy')
