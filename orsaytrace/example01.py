import trace as ot
import numpy

z_array = numpy.linspace(-4.5, 2.0, 11)
f = 2.0
lens_pos = -2.5
res = 0.1

### 1st PART ###

a = ot.Simu(5, 5, 10, res)

#Plane convex lens. Source is point source diverging.
a.d2_source(0.0, [0, 0, -5.0], [0, 0, 1], 0.12, 11)
a.create_thin_lens([0, 0, lens_pos], f, 1.75, 1.43, 'plane-convex')

for z in z_array:
    a.create_analysis_plan([0, 0, 1], z)

a.show_created_elements('all-noplan')
my_pl = a.run()
a.show_elements(my_pl, 'photons')


print('Running 2nd part.')
### 2ND PART ###

a = ot.Simu(5, 5, 10, res)

#Convex plane lens. Extended source collimated.
a.d2_source(0.5, [0, 0, -5.0], [0, 0, 1], 0.00, 1)
a.create_thin_lens([0, 0, lens_pos], f, 1.5, 1.43, 'convex-plane')

for z in z_array:
    a.create_analysis_plan([0, 0, 1], z)

a.show_created_elements('all-noplan')
my_pl = a.run(xsym=True, ysym=True)
a.show_elements(my_pl, 'photons')


print('Running 3rd part.')
### 3RD PART ###

a = ot.Simu(5, 5, 10, res)

lens_pos = lens_pos-0.72 #Here we put lens 0.82 backwards

#Plane convex lens. Source is point source diverging.
a.d2_source(0.0, [0, 0, -5.0], [0, 0, 1], 0.12, 11)
a.create_thin_lens([0, 0, lens_pos], f, 1.5, 1.43, 'plane-convex')

for z in z_array:
    a.create_analysis_plan([0, 0, 1], z)


a.show_created_elements('all-noplan')
my_pl = a.run()
a.show_elements(my_pl, 'photons')
