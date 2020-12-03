import orsaytrace.trace as ot
import numpy
import multiprocessing
import time

for nproc in [1, 2, 4, 8, 12]:
    if __name__ == "__main__":

        x, y, z, res = 5, 5, 5, 0.05
        z_array = numpy.linspace(-z / 4, +z / 4, 30)
        angles = 2

        start = time.perf_counter()

        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        a = ot.Simu(x, y, z, res)

        a.d2_source(0.0, [0, 0, -z/4], [0, 0, 1], 0.39, angles)

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
            a.merge_photon_lists(return_dict.values()[index])

        end = time.perf_counter()

        f = open('mp.txt', 'a+')
        f.write(str(angles) + '_' + str(nproc) + '_' + str(end - start)+ '\n')
        f.close()





