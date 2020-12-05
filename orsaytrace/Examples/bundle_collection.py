import orsaytrace.trace as ot
import numpy
import matplotlib.pyplot as plt

for tr in [3.0, 2.0, 1.0, 0.5, 0.25, 0.125]:
    for mesh in [1, 2, 4]:
        x, y, z, p, res = 5., 5., 5., 1.0, 0.02
        min_x, max_x = -1.0, 1.0
        min_z, max_z = 2.0-p/2.-0.65, 2.0-p/2.+0.35
        hist = 50
        a = ot.Simu(x, y, z, res)

        for x in numpy.linspace(min_x, max_x, hist):
            for z in numpy.linspace(min_z, max_z, hist):
                a.d2_flex_source(0.0, [x, 0, z], [0, -1, 0], 1.0, 9, 'point')

        a.create_parabolic_surface_element([0.0, 0, 2.0], -1.0, 2, 3.0, p)
        a.create_rectangle_element([-x/2, x/2, 0-0.3, y/2, -z/2, z/2], 1.0, [0, 0, 0])

        a.create_analysis_plan([0, 0, 1], -2.0, reflection_count = (1, 1))
        #Lets create a bundle of plans
        it = 0
        #tr = 1.0
        xc, yc, zc = 0, -1.0, -2.0
        #mesh = 1
        unitr = tr / mesh
        x = numpy.linspace(xc-tr/2.+unitr/2., xc+tr/2.-unitr/2., mesh)
        y = numpy.linspace(yc-tr/2.+unitr/2., yc+tr/2.-unitr/2., mesh)
        for xpos in x:
            for ypos in y:
                if (xpos-xc)**2+(ypos-yc)**2<=(tr-unitr)**2:
                    it=it+1
                    print(f'Plans centered at {xpos} and {ypos}.')
                    a.create_analysis_plan([0, 0, 1], zc, reflection_count = (1, 1), pos=([xpos, ypos, zc], (0, unitr)))

        #a.show_created_elements('photons')
        pl = a.run()
        #a.show_elements(pl, mode='-verbose')
        #a.show_photons2D(pl, mode='-verbose', binning=100)

        pl = pl[1:]

        x_new_list = [[ photon.init['pos'][0] for photon in pl[i].photons] for i in range(len(pl))]
        y_new_list = [[ photon.init['pos'][1] for photon in pl[i].photons] for i in range(len(pl))]
        z_new_list = [[ photon.init['pos'][2] for photon in pl[i].photons] for i in range(len(pl))]

        fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False)
        if it==1:
            axes.hist2d(x_new_list[0], z_new_list[0], bins=hist, range=[[min_x, max_x], [min_z, max_z]])
            plt.savefig(str(tr) + '_dum' + str(mesh) + '.png')
        else:
            for i in range(it):
                #plt.clf()
                axes.hist2d(x_new_list[i], z_new_list[i], bins=hist, range=[[min_x, max_x], [min_z, max_z]])
                plt.savefig(str(tr) + '_' + str(mesh) + '_dum' + str(i) + '.png')

        #plt.show()
