import parabolic as ot
import numpy

z_array = numpy.linspace(-4.5, 2.0, 11)
f = 2.0
lens_pos = -2.5
res = 0.03

a = ot.Simu(5, 5, 10, res)


#Plane convex lens. Source is point source diverging.
a.d2_source(0.0, -5.0, [0, 0, 1], 0.12, 31)
a.create_sphere_element([0.0, 0.0, lens_pos], f/2., 1.43) #from 2 to 3.
a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos-f, lens_pos], 1.0, [0, 0, 1], inclusive=False)
a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos, lens_pos], 1.43, [0, 0, 1], inclusive=True)


for z in z_array:
    a.create_analysis_plan([0, 0, 1], z)


a.show_elements(False, 'all')

oi

a.run()
a.show_elements(True, 'photons')
a.show_photons2D('xy')

### 2ND PART ###

a = ot.Simu(5, 5, 10, res)

#Convex plane lens. Extended source collimated.
a.d2_source(0.5, -5.0, [0, 0, 1], 0.00, 1)
a.create_sphere_element([0.0, 0.0, lens_pos], f/2., 1.43) #from 2 to 3.
a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos, lens_pos+f], 1.0, [0, 0, 1]) 


for z in z_array:
    a.create_analysis_plan(z=z)


a.show_elements(False, 'all')
a.run()
a.show_elements(True, 'photons')
a.show_photons2D('xy')

### 3RD PART ###

a = ot.Simu(5, 5, 10, res)

lens_pos = lens_pos-0.82 #Here we put lens 0.82 backwards

#Plane convex lens. Source is point source diverging.
a.d2_source(0.0, -5.0, [0, 0, 1], 0.12, 31)
a.create_sphere_element([0.0, 0.0, lens_pos], f/2., 1.43) #from 2 to 3.
a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos-f, lens_pos], 1.0, [0, 0, 1], inclusive=False)
a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos, lens_pos], 1.43, [0, 0, 1], inclusive=True)


for z in z_array:
    a.create_analysis_plan(z=z)


a.show_elements(False, 'all')
a.run()
a.show_elements(True, 'photons')
a.show_photons2D('xy')
