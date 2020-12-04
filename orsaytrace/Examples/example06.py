import orsaytrace.trace as ot
import numpy
import matplotlib.pyplot as plt

x, y, z, res = 5, 5, 5, 0.05

a = ot.Simu(x, y, z, res)

#a.create_chessboard(1.0, [0, 0, -2.4], 5)
a.d2_source(0.5, [0, 0, -2.5], [0, 0, 1])
a.create_analysis_plan([0, 0, 1], 1.9, pos=([0, 0, 0], (0, 5.0)))
a.show_created_elements('all')
pl = a.run()
#a.show_elements(pl)

#for photon in pl[0].photons:
#    print(photon.pos, photon.__dict__)
