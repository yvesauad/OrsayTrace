import orsaytrace.trace as ot
import numpy
import concurrent.futures
import time

x, y, z, res = 5, 5, 5, 0.02
y_array = numpy.linspace(0, y/4, 201)
process = 12

start = time.clock()

if __name__ == "__main__":
    with concurrent.futures.ProcessPoolExecutor(max_workers = process) as executor:

        a = ot.Simu(x, y, z, res)

        a.d2_source(0.2, [0, 0, -z/2], [0, 0, 1], 0.0, 1)

        for y in y_array:
            a.create_analysis_plan([0, 0, 1], y)
        
        future_values = [executor.submit(a.run, i, process) for i in numpy.arange(0, process)]
        
        for index, futures in enumerate(future_values):
            new_photon_list = a.merge_photon_lists(futures.result())

        end = time.clock()
        print(end - start)


