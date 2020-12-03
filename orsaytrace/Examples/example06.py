import orsaytrace.trace as ot
import numpy
import matplotlib.pyplot as plt

x, y, z, res = 5, 5, 5, 0.03

a = ot.Simu(x, y, z, res)

a.d2_source(0.2, [0, 0, 0], [0, 0, 1])

#a.create_chessboard(1.0, [0, 0, 0], 4, [0, 0, 1])
#a.create_fiberbundle(1.0, [0, 0, -2.5], 8, [0, 0, 1])

a.create_analysis_plan([0, 0, 1], 0, pos = 1)

a.run()