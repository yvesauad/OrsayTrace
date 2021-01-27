import orsaytrace.trace as ot
import numpy
import matplotlib.pyplot as plt

for p in [2.0, 1.0]:
    for tr in [0.5, 3.0, 2.0, 1.0, 0.5, 0.25, 0.125]:
        x, y, z, res = 5., 5., 5., 0.05
        y_dist, na, angles = 0.3, 0.00, 3
        xc, yc, zc = 0, -0.5, -2.0 #Fiber Central Position
        min_x, max_x = -1.0, 1.0
        min_z, max_z = 2.0-p/2.-0.65, 2.0-p/2.+0.35
        hist = 25
        a = ot.Simu(x, y, z, res)

        a.d2_flex_source(tr, [xc, yc, -2.0], [0, 0, 1], na, angles, 'xy')

        a.create_parabolic_surface_element([0.0, 0, 2.0], -1.0, 2, 3.0, p)
        a.create_rectangle_element([-x/2., x/2., 0.0-y_dist, y/2., -z/2., z/2.], 1.0, [0, 0, 0])

        a.create_analysis_plan([0, 1, 0], 0.05, reflection_count=(1, 1))
        a.create_analysis_plan([0, 1, 0], 0.0, reflection_count = (1, 1))
        a.create_analysis_plan([0, 1, 0], -0.05, reflection_count=(1, 1))
        a.create_analysis_plan([0, 0, 1], -1.0, reflection_count=(0, 0))

        a.show_created_elements('photons')
        pl = a.run()
        a.show_elements(pl, mode='-verbose')
        a.show_photons2D(pl, mode='-verbose', binning=100)

        pl = pl[1:]

        x_new_list = [[ photon.init['pos'][0] for photon in pl[i].photons] for i in range(len(pl))]
        y_new_list = [[ photon.init['pos'][1] for photon in pl[i].photons] for i in range(len(pl))]
        z_new_list = [[ photon.init['pos'][2] for photon in pl[i].photons] for i in range(len(pl))]


        fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False)
        theta = numpy.linspace(0, 2* numpy.pi, 180)
        a = tr * numpy.cos(theta) + xc
        b = tr * numpy.sin(theta) + yc
        axes.plot(a, b, linewidth=2)
        for ipt, pt in enumerate(pts):
            a = unitr/1. * numpy.cos(theta) + pt[0]
            b = unitr/1. * numpy.sin(theta) + pt[1]
            axes.text(pt[0], pt[1], str(ipt))
            axes.plot(a, b)
        plt.savefig('map_' + str(p) + '_' + str(tr) + '_' + str(mesh) + '.png')
        plt.clf()

        if it==1:
            fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False)
            plt.hist2d(x_new_list[0], z_new_list[0], bins=hist-1, range=[[min_x, max_x], [min_z, max_z]])
            plt.colorbar()
            plt.savefig(str(p) + '_' + str(tr) + '_' + str(mesh) + '.png')
        else:
            for i in range(it):
                fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False)
                #plt.clf()
                plt.hist2d(x_new_list[i], z_new_list[i], bins=hist-1, range=[[min_x, max_x], [min_z, max_z]])
                plt.colorbar()
                plt.savefig(str(p) + '_' + str(tr) + '_' + str(mesh) + '_' + str(i) + '.png')