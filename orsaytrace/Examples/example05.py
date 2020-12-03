import orsaytrace.trace as ot
import numpy
#import concurrent.futures
import multiprocessing
import time

x, y, z, res = 5, 5, 5, 0.05
z_array = numpy.linspace(-z/4, +z/4, 30)
nproc = 1

start = time.perf_counter()

if __name__ == "__main__":

    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    a = ot.Simu(x, y, z, res)

    a.d2_source(0.0, [0, 0, -z/4], [0, 0, 1], 0.39)

    for z in z_array:
        a.create_analysis_plan([0, 0, 1], z)

    a.prepare_acquisition(nproc)

    jobs = []
    for i in numpy.arange(0, nproc):
        p = multiprocessing.Process(target=a.run, args=(i, True, return_dict))
        jobs.append(p)
        p.start()

    for index, proc in enumerate(jobs):
        proc.join()
        if index==0: a.merge_photon_lists(return_dict.values()[index])

    #a.show_elements(a.photon_lists)

    end = time.perf_counter()
    print(end - start)


