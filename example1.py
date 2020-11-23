import parabolic as ot
import numpy

'''a = ot.Simu(5, 5, 10, 0.05)

lens_pos = -2.5
f = 2.0

#Plane convex lens. Source is point source diverging.
a.d2_source(0.0, -5.0, [0, 0, 1], 0.22, 31)
a.create_sphere_element([0.0, 0.0, lens_pos], f/2., 1.43) #from 2 to 3.
a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos-f, lens_pos], 1.0, [0, 0, 1]) 


for z in numpy.linspace(-4.8, 4.8, 11):
    a.create_analysis_plan(z=z)


a.show_elements(False, 'all')
a.run()
a.show_elements(True, 'photons')
a.show_photons2D('xy')'''

############ 2ND PART #############

a = ot.Simu(5, 5, 10, 0.05)

lens_pos = -2.5
f = 2.0

#Convex plane lens. Extended source collimated.
a.d2_source(0.25, -5.0, [0, 0, 1], 0.00, 1)
a.create_sphere_element([0.0, 0.0, lens_pos], f/2., 1.43) #from 2 to 3.
a.create_rectangle_element([-1.1, 1.1, -1.1, 1.1, lens_pos, lens_pos+f], 1.0, [0, 0, 1]) 


for z in numpy.linspace(-4.8, 4.8, 11):
    a.create_analysis_plan(z=z)


a.show_elements(False, 'all')
a.run()
a.show_elements(True, 'photons')
a.show_photons2D('xy')
