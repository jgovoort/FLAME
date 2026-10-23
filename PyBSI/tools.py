import numpy as np
import xarray as xr


def find_nearest(array, value):
    """
    Find the index of the nearest value in the array.

    Parameters
    ----------
    array : array-like
        A list or numpy array of numerical values.
    value : float or int
        The target value to find the closest element to.

    Returns
    -------
    int
        The index of the element in the array closest to the input value.
    
    Examples
    --------
    >>> find_nearest([1, 3, 7, 10], 8)
    2
    
    >>> find_nearest(np.array([1.5, 2.8, 3.3]), 2)
    1
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def verif_date(absTime) :
    """
    Verify if the day has changed within the given array of timestamps.

    This function checks if there is a change in the calendar day between the 
    first and the last timestamps in the input array.

    Parameters
    ----------
    absTime : array-like
        A list or numpy array of timestamps (e.g., in ISO format or compatible 
        with `np.datetime64`).

    Returns
    -------
    bool
        True if the day has changed between the first and the last timestamp in 
        the array, False otherwise.
    
    Examples
    --------
    >>> verif_date(['2023-10-10T10:15', '2023-10-11T09:45'])
    True

    >>> verif_date(['2023-10-10T10:15', '2023-10-10T23:45'])
    False
    """
    dt = np.datetime64(absTime[-1], 'D') - np.datetime64(absTime[0], 'D') 
    if dt.astype(int) > 0:
        return True
    return False

def export_nc(wav, time, refwav=None, Tempref=None, strainref=None, path='./', mode='LS'): 
    """
    Export BSI FBG data from two dictionaries (data and timing) into a NetCDF file.

    This function takes two dictionaries containing data and timing information, and
    exports them into a NetCDF file. The function supports different sampling modes 
    and can handle additional optional reference parameters for wavelengths, 
    temperature, and strain.

    Parameters
    ----------
    wav : dict or np.ndarray
        Dictionary (PerCH) or np.ndarray (LS) containing the data to be exported. The keys are typically the 
        channels number, and the values are the corresponding data arrays.
    
    time : dict or np.ndarray
        Dictionary (PerCH) or np.ndarray (LS) containing the timing information associated with the data. 
        The keys are the same as in `dict_data`, and the values are arrays of 
        timestamps.
    
    refwav : np.ndarray, optional
        Boolean array indicating reference wavelengths. If not provided, it will be 
        initialized with NaN arrays. This parameter is mainly used for temperature 
        and strain reference calibration (default is None).
    
    Tempref : np.ndarray, optional
        Wavelength corresponding to a temperature of 20°C. Used only for temperature 
        probes. If not provided, it will be initialized as NaN (default is None).
    
    strainref : np.ndarray, optional
        Wavelength corresponding to a reference strain value. If not provided, it 
        will be initialized as NaN (default is None).
    
    path : str, optional
        Path where the output NetCDF file will be saved (default is './').
    
    mode : str, optional
        Sampling mode. Options are:
        - 'HS': High Sampling frequency (only one channel),
        - 'LS': Low Sampling frequency,
        - 'DUAL': Mixed time sampling.
        Default is 'LS'.
    
    Returns
    -------
    None
        This function saves a NetCDF file at the specified `path` but does not 
        return any value.

    Notes
    -----
    - The shape of `dict_data` and `dict_time` must be compatible for exporting 
      into a dataset (e.g., time dimensions should match).
    - The output NetCDF file will be named using the first timestamp in `dict_time`.

    Examples
    --------
    >>> dict_data = {'channel_1': np.random.rand(10, 100), 'channel_2': np.random.rand(10, 100)}
    >>> dict_time = {'time': np.linspace(0, 100, 100)}
    >>> export_nc(dict_data, dict_time, path='/data/exports/', mode='LS')

    #### MUST BE IMPLEMENTED ####
    
    refwav : Bool array
    tempref : Wavelength at 20°C --> Only for the temperature probes
    strainref : Wavelength with reference value of strain 
    
    """
    
    #if refwav is None :
    #    refwav = OutWav[:,:,0]
    #    Tempref = np.full(refwav.shape,np.nan)
    #    strainref = np.full(refwav.shape, np.nan)
    
    if mode == 'LS' : 
        data = xr.Dataset(data_vars=dict(Wavelength_nm = (['Channel',"FBG","Time"], wav),
                                    Calibration_Wavelength = (['Channel', 'FBG'], refwav),
                                    Temperature_calibration = (['Channel', 'FBG'], Tempref),
                                    Strain_calibration = (['Channel', 'FBG'], strainref)),
                      coords = dict(Time = time),
                      attrs=dict(description="FBG data -  Dimension of the data (Channel*FBG*Time)"))
        OutName = str(time[0]).replace('-','.').replace(':','.')
        data.to_netcdf(path+'/FBG_%s_%s.nc' %(mode, OutName))
        
    elif mode == 'PerCH':
        sttime = wav[list(wav.keys())[0]][0]
        sttime = sttime.astype('datetime64[s]')
        data = xr.Dataset(data_vars=wav,
                          coords=time,
                          attrs=dict(description="FBG data -  Dimension of the data (Channel*FBG*Time)"))                
        OutName = str(sttime).replace('-','.').replace(':','.')
        
        print("Exported data to %sFBG_%s_%s.nc" %(path, mode, OutName))
        data.to_netcdf(path+'/FBG_%s_%s.nc' %(mode, OutName))