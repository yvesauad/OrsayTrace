import matplotlib.pyplot as plt
import os
import numpy

#my computer. 121 angles. 30 plans, cell 5, 5, 5 and res 0.05. cpu_count = 8.
proc = [1, 2, 4, 8, 12]
time = [1355.5, 712.2, 567.4, 605.0, 620.8]
time_linux = [1852.3, 804.2, 502.3, 520.8, 559.1]
time_patrick = [1766.68, 852.75, 606.66, 611.52, 645.72]
pc = [45140, 22570, 11285, 5642, 3762] #photons / core

time = numpy.divide(time, 3600)
time_linux = numpy.divide(time_linux, 3600)
time_patrick = numpy.divide(time_patrick, 3600)

#Patrick running same code as my computer using single core. 30 plans 121 angles
# 45140 photons/core, 1 core, 1165.87s (749s running)


'''fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, dpi=200)

axes.plot(proc, time, color='red', alpha=0.25)
axes.scatter(proc, time, color='red', label='W10 - Intel I7 8665U')
axes.plot(proc, time_linux, color='orange', alpha=0.25)
axes.scatter(proc, time_linux, color='orange', label='U20.04 - Intel I7 8665U')
axes.plot(proc, time_patrick, color='green', alpha=0.25)
axes.scatter(proc, time_patrick, color='green', label='W10 - Xeon Silver 4214')

axes.set_ylabel('Time (hours)')
axes.tick_params(axis='y', labelcolor='red')
axes.axvline(8, linestyle='--', label = 'I7 8665U #core')
axes.axvline(24, linestyle='dashdot', label = 'Xeon Silver 4214 #core')

axes2 = axes.twinx()
axes2.plot(proc, pc, color='blue', alpha=0.15)
axes2.scatter(proc, pc, color='blue')
axes2.set_ylabel('Photons/core', color='blue')
axes2.tick_params(axis='y', labelcolor='blue')

axes.set_xlabel('Number of Processes')
axes.set_yscale('log')
axes.set_xscale('log')
axes2.set_yscale('log')
axes2.set_xscale('log')

axes.legend()
plt.show()'''

#*************************************************************************

##Patrick. 101 angles, 1 plans, cell 1, 1, 320 and res 0.02. cpu_count = 48. 1305 photons/core with
#ncore = 24. Much more paralell programming

proc = [1, 2, 4, 8, 12, 18, 24, 30, 36, 42]
time =[14444.24, 7504.97, 3775.15, 1967.66, 1370.34, 999.13, 836.31, 814.64, 814.11, 802.43]
pc = [31320, 15660, 7830, 3915, 2610, 1740, 1305, 1044, 870, 746]



#Patrick. 221 angles, 30 plans, cell 5, 5, 5 and res 0.05. cpu_count = 48
#proc = [1, 2, 4, 8]
#time = [19510.85, 9761.39, 6500.21, 5637.74]
#run_time = [14206, 4085, 1170, 391]
#pc = [161603, 80802, 40401, 20200]

time = numpy.divide(time, 3600)
#run_time = numpy.divide(run_time, 3600)

fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, dpi=200)

axes.plot(proc, time, color='red', alpha=0.25)
axes.scatter(proc, time, color='red', label='W10 - Xeon Silver 4214')
#axes.plot(proc, run_time, color='orange', alpha=0.25)
#axes.scatter(proc, run_time, color='orange', label='W10 - Xeon Silver 4214 (Run Time)')

axes.set_ylabel('Time (hours)', color='red')
axes.tick_params(axis='y', labelcolor='red')
axes.axvline(24, linestyle='dashdot', label = 'Xeon Silver 4214 #core')

axes2 = axes.twinx()
axes2.plot(proc, pc, color='blue', alpha=0.15)
axes2.scatter(proc, pc, color='blue')
axes2.set_ylabel('Photons/core', color='blue')
axes2.tick_params(axis='y', labelcolor='blue')

axes.set_xlabel('Number of Processes')
axes.set_yscale('log')
axes.set_xscale('log')
axes2.set_yscale('log')
axes2.set_xscale('log')

axes.legend()

plt.show()