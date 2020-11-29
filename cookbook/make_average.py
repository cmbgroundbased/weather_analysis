import pycal
import datetime
from mpi4py import MPI
import numpy as np
import yaml
from scipy.integrate import simps

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# All the threads must have the 

N_years = 72 # Limited by the INDACO wall-time

def bandshape(focal_plane, detectors, det):
    id_pol = focal_plane['horns'][det]['polarimeter_id']
    for i in detectors:
        if i['id'] == id_pol:
            bandshape = i['bandpass']['bandshape']
            banderror = i['bandpass']['bandshape_error']
            f_max = i['bandpass']['highest_frequency_hz']
            f_min = i['bandpass']['lowest_frequency_hz']
            samples = i['bandpass']['num_of_frequencies']
            
            freq = np.linspace(f_min, f_max, samples)
            
    return freq, bandshape, banderror


# All the threads have to read weather, focal_plane and detectors files
weather = pycal.Weather("./weather_STRIP.fits", 1, 1)
with open(r'strip_focal_plane.yaml') as file:
    focal_plane = yaml.full_load(file)
    
with open(r'strip_detectors.yaml') as file:
    detectors = yaml.full_load(file)
    
Q_detectors = ["B0", "B1", "B2", "B3", "B4", "B5", "B6",
               "G0", "G1", "G2", "G3", "G4", "G5", "G6",
               "I0", "I1", "I2", "I3", "I4", "I5", "I6",
               "O0", "O1", "O2", "O3", "O4", "O5", "O6",
               "R0", "R1", "R2", "R3", "R4", "R5", "R6",
               "V0", "V1", "V2", "V3", "V4", "V5", "V6",
               "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6"]

# The win array is shared between all the threads but all the results
# have to be gathered in a single array allocated by the "rank == 0"
if rank == 0:
    nbytes = (49*24*12*27*N_years)*MPI.DOUBLE.Get_size()
    
    print("Start to make average from 1948 to 2020...")
    print("You have to free at least: "+str(nbytes/(1024*1024))+" MB" )
    print("Comunicator size: "+str(size))
    print("Years per thread: "+str(int(N_years/size)))
else:
    nbytes = 0

# Allocate the shared array. If rank==0 is an actual allocation in memory else
# it's just istantiate as MPI shared object
win = MPI.Win.Allocate_shared(nbytes, MPI.DOUBLE.Get_size(), comm=comm)
buf, itemsize = win.Shared_query(0)

# The array allocated by the rank = 0 is a numpy array.
t_atm_40GHz_K = np.ndarray(buffer=buf, dtype='d', shape=(49, 12, 24, 27*N_years))

# The number of the years that every thread have to evaluate
threads = size
years_per_threads = int(N_years/size)

# A global index that rely on the rank number
global_index = rank*years_per_threads*27


for year in range(rank*years_per_threads+1949, (rank+1)*years_per_threads+1949):
    
    print("Thread: "+str(rank)+" has taken in charge the year "+str(year))
    t1 = MPI.Wtime()
    
    for day in range(1, 28):
        
        for month in range(1, 13):
            
            for hour in range(0, 24):
                #t3 = MPI.Wtime()
                
                primes = 0
                seconds = 0

                data = datetime.datetime(year, month, day, hour, primes, seconds)
                timestamp = data.timestamp()
                weather.set_time(int(timestamp))

                pwv=weather.pwv
                t0=weather.surface_temperature
                p0=weather.surface_pressure
                
                # Convolve the atmospheric emission with the detectors bandshapes
                det_index = 0
                for det in Q_detectors:
                    
                    index_decimated= np.round(np.linspace(0, 120, 40))
                    index_decimated = np.array(index_decimated, dtype=int)
                    
                    freq_det, band_det, erro_det = bandshape(focal_plane, detectors, det)
                    freq_det = np.array(freq_det)
                    band_det = np.array(band_det)
                    erro_det = np.array(erro_det)
                    
                    freq_det = freq_det[index_decimated]
                    band_det = band_det[index_decimated]
                    erro_det = erro_det[index_decimated]
                    
                    atm_emission = pycal.atm_atmospheric_loading_vec(2390.0, t0, p0, pwv, np.amin(freq_det/1e9), np.amax(freq_det/1e9), len(freq_det))
                    t_atm = simps(band_det*atm_emission, freq_det) / simps(band_det, freq_det)
                    t_atm_40GHz_K[det_index, month-1, hour, global_index] += t_atm
                    det_index = det_index + 1
                # print("Rank: "+str(rank)+" has just finished the whole focal_plane")
                #t4 = MPI.Wtime()
                #print("Rank "+str(rank)+" fnished hour "+str(hour)+" in "+str(t4-t3)+" sec.")
                #print("Days: "+str((t4-t3)*24*12*27*(N_years/size)/(3600*24)))

        global_index = global_index + 1

    t2 = MPI.Wtime()
    print("Rank "+str(rank)+" finished year "+str(year)+" in "+str(t2-t1)+" sec.")
comm.Barrier()

print("Thread: "+str(rank)+" is synchronized")
if rank == 0:
    print("Gather and write the result on file...")
    print("Number of samples: "+"27")
    print("Global index: "+str(global_index))
    np.save("t_atm_par_bandshapes", t_atm_40GHz_K)