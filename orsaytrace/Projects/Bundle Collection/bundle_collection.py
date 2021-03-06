import orsaytrace.trace as ot
import numpy
import matplotlib.pyplot as plt
import os

SAVEPATH = os.path.join(os.path.dirname(__file__), "Data") + "\\"

for p in [2.0]:
    for tr in [2.0, 1.0, 0.5, 0.25]:
        for mesh in [1, 3, 5]:
            x, y, z, res = 6., 7., 5., 0.02 #Grid and Resolution
            y_dist, na, angles = 0.3, 2.0, 13 #y_dist is cut in Z. na is aperture of source and angles is discretization
            xc, yc, zc = 0, -1.65, -2.0 #Fiber Central Position
            min_x, max_x = -2.5, 2.5 #min and max of the source in X
            min_z, max_z = 2.0-p/2.-1.80, 2.0-p/2.+0.80 #min and max of the source in Z
            hist = 50 #Histogram discretization
            a = ot.Simu(x, y, z, res)

            for x_src in numpy.linspace(min_x, max_x, hist):
                for z_src in numpy.linspace(min_z, max_z, hist):
                    a.d2_flex_source(0.0, [x_src, 0, z_src], [0, -1, 0], na, angles, 'point') #Y Focus at 0. Directed to -Y.

            a.create_parabolic_surface_element([0.0, 0, 2.0], -1.0, 6, 5.0, p)
            a.create_rectangle_element([-x/2., x/2., 0.0-y_dist, y/2., -z/2., z/2.], 1.0, [0, 0, 0])

            a.create_analysis_plan([0, 0, 1], -2.0, reflection_count = (1, 1))

            #Lets create a bundle of plans
            pts = list()
            it = 0
            unitr = tr / mesh
            x = numpy.linspace(xc-tr+unitr, xc+tr-unitr, mesh)
            y = numpy.linspace(yc-tr+unitr, yc+tr-unitr, mesh)

            for xpos in x:
                for ypos in y:
                    if (xpos-xc)**2+(ypos-yc)**2<=(tr-unitr)**2:
                        it=it+1
                        pts.append((xpos, ypos))
                        print(f'Creating plans centered at {xpos} and {ypos}.')
                        a.create_analysis_plan([0, 0, 1], zc, reflection_count = (1, 1), pos=([xpos, ypos, zc], (0, unitr)))

            #a.show_created_elements('photons')
            #a.show_created_elements('all')
            pl = a.run()
            #a.show_elements(pl, mode='-verbose')
            #a.show_photons2D(pl, mode='-verbose', binning=100)

            pl = pl[1:] #Discard first plan

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
            plt.savefig(SAVEPATH + 'map_' + str(p) + '_' + str(tr) + '_' + str(mesh) + '.png')
            plt.clf()

            if it==1:
                fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False)
                plt.hist2d(x_new_list[0], z_new_list[0], bins=hist-1, range=[[min_x, max_x], [min_z, max_z]])
                plt.colorbar()
                plt.savefig(SAVEPATH + str(p) + '_' + str(tr) + '_' + str(mesh) + '.png')
            else:
                for i in range(it):
                    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False)
                    plt.hist2d(x_new_list[i], z_new_list[i], bins=hist-1, range=[[min_x, max_x], [min_z, max_z]])
                    plt.colorbar()
                    plt.savefig(SAVEPATH + str(p) + '_' + str(tr) + '_' + str(mesh) + '_' + str(i) + '.png')