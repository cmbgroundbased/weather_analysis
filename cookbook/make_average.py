import pycal
import datetime
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

weather = pycal.Weather("./weather_STRIP.fits", 1, 1)

if rank == 0:
    nbytes = (24*12*27*120)*MPI.DOUBLE.Get_size()
    
    print("Start to make average from 1900 to 2020...")
    print("You have to free at least: "+str(nbytes/(1024*1024))+" MB" )
    print("Comunicator size: "+str(size))
    print("Years per thread: "+str(int(120/size)))
else:
    nbytes = 0

win = MPI.Win.Allocate_shared(nbytes, MPI.DOUBLE.Get_size(), comm=comm)
buf, itemsize = win.Shared_query(0)

t_atm_40GHz_K = np.ndarray(buffer=buf, dtype='d', shape=(12, 24, 27*120))

threads = size
years_per_threads = int(120/size)

global_index = rank*years_per_threads*27

for y in range(rank*years_per_threads+1901, (rank+1)*years_per_threads+1901):
    print("Thread: "+str(rank)+" has taken in charge the year "+str(y))
    for j in range(1, 28):
        day = j
        for i in range(1, 13):
            year = y
            month = i

            for k in range(0, 24):
                hour = k
                primes = 0
                seconds = 0

                data = datetime.datetime(year, month, day, hour, primes, seconds)
                timestamp = data.timestamp()
                weather.set_time(int(timestamp))

                pwv=weather.pwv
                t0=weather.surface_temperature
                p0=weather.surface_pressure
                t_atm = pycal.atm_atmospheric_loading(2390.0, t0, p0, pwv, 43)

                t_atm_40GHz_K[i-1, k, global_index] += t_atm

        global_index = global_index + 1

comm.Barrier()

print("Thread: "+str(rank)+" is synchronized")
if rank == 0:
    print("Gather and write the result on file...")
    print("Number of samples: "+"27")
    print("Global index: "+str(global_index))
    np.save("t_atm_par", t_atm_40GHz_K)