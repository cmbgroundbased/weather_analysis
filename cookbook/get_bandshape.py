import argparse
import yaml
import pylab as plt
import numpy as np

def bandshape(focal_plane, detectors, det):
    id_pol = focal_plane['horns'][det]['polarimeter_id']
    for i in detectors:
        if i['id'] == id_pol:
            bandshape = i['bandpass']['bandshape']
            banderror = i['bandpass']['bandshape_error']
            f_max_hz = i['bandpass']['highest_frequency_hz']
            f_min_hz = i['bandpass']['lowest_frequency_hz']
            samples = i['bandpass']['num_of_frequencies']
            
            freq_hz = np.linspace(f_min_hz, f_max_hz, samples)
            
    return freq_hz, bandshape, banderror



if __name__ == "__main__":
    
    with open(r'strip_focal_plane.yaml') as file:
        focal_plane = yaml.full_load(file)
        
    with open(r'strip_detectors.yaml') as file:
        detectors = yaml.full_load(file)
    
    # lista dei detectors in banda Q
    # Q_detectors = ["B0", "B1", "B2", "B3", "B4", "B5", "B6",
    #                "G0", "G1", "G2", "G3", "G4", "G5", "G6",
    #                "I0", "I1", "I2", "I3", "I4", "I5", "I6",
    #                "O0", "O1", "O2", "O3", "O4", "O5", "O6",
    #                "R0", "R1", "R2", "R3", "R4", "R5", "R6",
    #                "V0", "V1", "V2", "V3", "V4", "V5", "V6",
    #                "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6"]

    freq, band, err = bandshape(focal_plane, detectors, "I0")
    
    plt.errorbar(freq/1E9, band, err)
    plt.xlabel("Frequency [GHz]")
    plt.show()



    
    
        
            
            