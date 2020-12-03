import orsaytrace.trace as ot
import numpy
import matplotlib.pyplot as plt

x, y, z, res = 5, 5, 10, 0.05

a = ot.Simu(x, y, z, res)

a.create_chessboard(1.0, [0, 0, 0], 5)

a.show_created_elements('all-noplan')