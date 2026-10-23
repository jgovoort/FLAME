import pandas as pd
import numpy as np
import glob 


def single_channel(path='./*_SPECTRUM.txt') :
    files = glob.glob(path)
    files.sort()
    for i, file in enumerate(files) :
        spectrum = pd.read_csv(file, delimiter='\t', index_col='Timestamp_CH')
        if i == 0 :
            arrspect = spectrum.to_numpy()
            timespect = spectrum.index
            wavelength = spectrum.columns.values
            continue
        arrspect = np.append(arrspect, spectrum.to_numpy(), axis=0)
        timespect = np.append(timespect, spectrum.index)
    timespect = np.char.replace(timespect.astype(str), ':,', '.')
    timespect = np.char.replace(timespect.astype(str), '_CH1', '')
    timespect = timespect.astype(np.datetime64)
    return arrspect, timespect, wavelength.astype(float)