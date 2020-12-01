import os
import sys

currentdir = os.getcwd()
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import trace as ot
import numpy
import matplotlib.pyplot as plt

########## SIMULATION PARAMETERS ###########

f = 2.0
lens_pos = -2.5
res = 0.04
pts = 201
sim_pts = 7

z_plans = numpy.linspace(-4.5, 4.5, pts)
z_lens = numpy.linspace(-2.5, -1.2, sim_pts)
all_my_photons = list()


for lens_pos in z_lens:


    a = ot.Simu(5, 5, 10, res)

    #Plane convex lens. Source is point source diverging.
    a.d2_source(0.0, [0, 0, -4.0], [0, 0, 1], 0.12, 3)
    a.create_thin_lens([0, 0, lens_pos], f, 1.75, 1.43, 'plane-convex')

    for z in z_plans:
        a.create_analysis_plan([0, 0, 1], z)

    #a.show_created_elements('all-noplan')
    all_my_photons.append(a.run())


results = numpy.zeros((sim_pts, pts))

for isim, simu_part in enumerate(all_my_photons):
    #a.show_elements(all_my_photons[isim], 'photons')
    for ilist, photon_list in enumerate(simu_part):
        results[isim, ilist] = (photon_list.std_position()[1])
            
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, dpi=200)

for index, result in enumerate(results):
    axes.plot(z_plans, result, label='z_lens = ' + format(z_lens[index], '.2f'))

axes.set_xlabel('Z')
axes.set_ylabel('stdY')
plt.legend()
plt.show()
