import parabolic as ot
import numpy

z_array = numpy.linspace(-4.5, 2.0, 11)
f = 2.0
lens_pos = -2.5
res = 0.05

### 1st PART ###

a = ot.Simu(5, 5, 10, res)

#Plane convex lens. Source is point source diverging.
a.d2_source(0.20, [0, 0, -5.0], [0, 0, 1], 0.0, 1)
a.d2_source(0.20, [0, 1.0, -5.0], [0, 0, 1], 0.55, 3)
a.d2_source(0.0, [1.0, 0.0, -5.0], [0, 0, 1], 0.12, 3)
#a.d2_source(0.0, [1.0, 1.0, -5.0], [0, 0, 1], 0.12, 31)
#a.create_sphere_element([0.0, 0.0, lens_pos], f/2., 1.43) #from 2 to 3.
#a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos-f, lens_pos], 1.0, [0, 0, 1], inclusive=False)
#a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos, lens_pos], 1.43, [0, 0, 1], inclusive=True)

for z in z_array:
    a.create_analysis_plan([0, 0, 1], z)

a.show_elements(False, 'all')
a.run()
a.show_elements(True, 'photons')
