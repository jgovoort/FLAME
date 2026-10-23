def fbg_converter_extended(path, folder = True, refwav=None, version='1.3', nChannel=8, maxFBGs=12) :
    '''
    IN
    ####
    
    version : Version of the BSI software
                1.3 --> Low sampling frequency
                1.6.1DUAL --> HSV version with an interrogation time different for each channel
    
    '''
    
    if folder == False :
        if version == '1.3' :
            dict_data, dict_time = read_LS_FBG(path, nChannel=nChannel, maxFBGs=maxFBGs, extended=True)
            mode = 'LS'
        elif version == '1.6.1DUAL' :
            dict_data, dict_time = bsi_hsv_dual_161_reader(path, maxFBGs =maxFBGs, nChannel=nChannel)
            mode = 'DUAL'
        export_nc(dict_data, dict_time, mode=mode)  
        return
    
    files = os.listdir(path)
    files.sort()
    
    ji=0
    for file in files:
        if file[-5:] != 'Z.txt':
            continue

        if version == '1.3' :
            dict_data, dict_time = read_LS_FBG(path+file, nChannel=nChannel, maxFBGs=maxFBGs, extended=True)
            mode = 'LS'
        elif version == '1.6.1DUAL' :
            dict_data, dict_time = bsi_hsv_dual_161_reader(path+file, maxFBGs=maxFBGs, nChannel=nChannel)
            mode = 'DUAL'
        if ji == 0 :
            ALL_dict_data, ALL_dict_time = dict_data.copy(), dict_time.copy()
            ji+=1 
            continue
            
        #print(file, TabsTime.shape, TrefTime.shape, TOutWav.shape)
        
        for i, name in enumerate(dict_data) :
            #print('Debug 1', name, ALL_dict_data[name][1].shape, dict_data[name][1].shape)

            ch = i+1
            name_time = name+'_time'
            
            ALL_dict_time[name_time] = np.append(ALL_dict_time[name_time], dict_time[name_time]) 
            out_chann_arr = np.append(ALL_dict_data[name][1], dict_data[name][1], axis=1)
            
            ALL_dict_data[name]  = (['Channel %i : FBGs' %ch, "Channel %i : Over time" %ch], out_chann_arr)  
            
            #print('Debug 2', name, ALL_dict_data[name][1].shape, dict_data[name][1].shape)

            del out_chann_arr
            
        if verif_date(ALL_dict_time['Channel_1_time']) == True :
            absTime = np.array(ALL_dict_time['Channel_1_time'], dtype='datetime64')
            IdxEndTime = np.where(absTime.astype('datetime64[D]')==absTime[-1].astype('datetime64[D]'))[0][0]

            OUT_dict_time = dict()
            
            OUT_dict_time['Channel_1_time'] = ALL_dict_time['Channel_1_time'][:IdxEndTime].copy()
            OUT_dict_data = ALL_dict_data.copy()
            OUT_dict_data['Channel_1']  = (['Channel 1 : FBGs' , "Channel 1 : Over time"], ALL_dict_data['Channel_1'][1][:,:IdxEndTime])
            
            export_nc(OUT_dict_data, OUT_dict_time, mode=mode)
            
            del OUT_dict_data, OUT_dict_time
            
            out_chann_arr = ALL_dict_data['Channel_1'][1][:,IdxEndTime:].copy()
            ALL_dict_time['Channel_1_time'] = ALL_dict_time['Channel_1_time'][IdxEndTime:].copy()
            ALL_dict_data['Channel_1']  = (['Channel 1 : FBGs' , "Channel 1 : Over time"], out_chann_arr) 

            del out_chann_arr
            
            for i, name in enumerate(dict_data) :
                ch = i+1
                name_time = name+'_time'
                if ch == 1 :
                    continue
                ALL_dict_time[name_time] = np.array([])
                ALL_dict_data[name] = (['Channel %i : FBGs'%ch , "Channel 1 : Over time"], np.empty((maxFBGs,0)))             
    try :
        export_nc(ALL_dict_data, ALL_dict_time, mode=mode)#,refwav=refwav)
    except : pass


def read_LS_FBG(filepath, refwav=None, nChannel = 8, maxFBGs=12, extended=False):
    '''
    Read BSI text files with multichannel
    
    return OutWav 3D Array (Channel*FBG*Time)
    '''
    absTime = []
    refTime = []
    with open(filepath, mode='r', newline='') as file:
        reader = csv.reader(file)
        # Iterate over each row in the CSV file
        for i, row in tqdm(enumerate(reader)):
            if i == 0:
                continue
            lst = str(row[0]).split('\t')
            wav = np.array(np.hsplit(np.array(lst[2:]), nChannel))
            wav = wav[:nChannel, :maxFBGs]
            wav = wav.reshape((wav.shape[0],wav.shape[1],1)).astype(float)
            
            if i == 1 and refwav is None :
                refwav = wav ## if no reference, designate the first value as reference. Used to check if no peak is missing
            
            if wav.shape == refwav.shape and (np.abs(wav-refwav).any() > 3.).any() : #threshold in nm or equivalent +/-300°C
                warnings.warn("A shift of more than 3 nm has been found. Trying to solve automatically the shift")
                wav = fbg_verifier(refwav, wav, nFBG=maxSpect)
            wav = np.where(wav==0,np.nan,wav)
            if i == 1 : # initiate the output array
                print('Starting file : '+filepath)
                lst = str(row[0]).split('\t')
                absTime.append(lst[0]), refTime.append(float(lst[1]))
                OutWav = wav                
                continue
                
            try :
                OutWav = np.append(OutWav, wav, axis=2)
                absTime.append(lst[0]), refTime.append(float(lst[1]))
            except :
                print(lst[0], wav.shape)
                
        absTime = np.array(absTime, dtype='datetime64')

    if extended == True :
        dict_data, dict_time = dict(), dict()
        for ch in range(1,nChannel+1) :
            dict_data['Channel_%i' %ch] = (['Channel %i : FBGs' %ch, "Channel %i : Over time" %ch], OutWav[ch-1,:,:])
            dict_time['Channel_%i_time' %ch] = absTime
        return dict_data, dict_time
    return absTime, np.array(refTime), OutWav, refwav


def fbg_converter_extended(path, folder = True, refwav=None, version='1.3', nChannel=8, maxFBGs=12) :
    '''
    IN
    ####
    
    version : Version of the BSI software
                1.3 --> Low sampling frequency
                1.6.1DUAL --> HSV version with an interrogation time different for each channel
    
    '''
    
    if folder == False :
        if version == '1.3' :
            dict_data, dict_time = read_LS_FBG(path, nChannel=nChannel, maxFBGs=maxFBGs, extended=True)
            mode = 'LS'
        elif version == '1.6.1DUAL' :
            dict_data, dict_time = bsi_hsv_dual_161_reader(path, maxFBGs =maxFBGs, nChannel=nChannel)
            mode = 'DUAL'
        export_nc(dict_data, dict_time, mode=mode)  
        return
    else :
        if version == '1.3' :
            folder_reader_LS_FBG(path, nChannel=nChannel, maxFBGs=maxFBGs)
        elif version == '1.6.1DUAL' :
            folder_reader_DUAL_FBG(path, nChannel=nChannel, maxFBGs=maxFBGs)



def folder_reader_LS_FBG(path, nChannel=8, maxFBGs=12) :
    files = os.listdir(path)
    files.sort()
    ji=0
    
    for file in files:
        if file[-5:] != 'Z.txt':
            continue

        dict_data, dict_time = read_LS_FBG(path+file, nChannel=nChannel, maxFBGs=maxFBGs, extended=True)
        mode = 'LS'
        
        if ji == 0 :
            ALL_dict_data, ALL_dict_time = dict_data.copy(), dict_time.copy()
            ji+=1 
            continue
            
        for i, name in enumerate(dict_data) :
            ch = i+1
            name_time = name+'_time'
            
            ALL_dict_time[name_time] = np.append(ALL_dict_time[name_time], dict_time[name_time]) 
            out_chann_arr = np.append(ALL_dict_data[name][1], dict_data[name][1], axis=1)
            
            ALL_dict_data[name]  = (['Channel %i : FBGs' %ch, "Channel %i : Over time" %ch], out_chann_arr)  
            
            del out_chann_arr
            
        if verif_date(ALL_dict_time['Channel_1_time']) == True :
            absTime = np.array(ALL_dict_time['Channel_1_time'], dtype='datetime64')
            IdxEndTime = np.where(absTime.astype('datetime64[D]')==absTime[-1].astype('datetime64[D]'))[0][0]

            OUT_dict_time, OUT_dict_data = dict(),dict()

            for ch in range(1, nChannel+1) :
                OUT_dict_time['Channel_%i_time' %ch] = ALL_dict_time['Channel_%i_time' %ch][:IdxEndTime].copy()
                OUT_dict_data['Channel_%i'%ch]  = (['Channel %i : FBGs' %ch , "Channel %i : Over time" %ch], ALL_dict_data['Channel_%i'%ch][1][:,:IdxEndTime])

                out_chann_arr = ALL_dict_data['Channel_%i' %ch][1][:,IdxEndTime:].copy()

                ALL_dict_time['Channel_%i_time' %ch] = ALL_dict_time['Channel_%i_time' %ch][IdxEndTime:].copy()
                ALL_dict_data['Channel_%i' %ch]  = (['Channel %i : FBGs' %ch , "Channel %i : Over time" %ch], out_chann_arr)
            
            export_nc(OUT_dict_data, OUT_dict_time, mode=mode)
            
            del OUT_dict_data, OUT_dict_time, out_chann_arr
            
    try :
        export_nc(ALL_dict_data, ALL_dict_time, mode=mode)#,refwav=refwav)
    except : pass


def bsi_hsv_dual_161_reader(filename, maxFBGs = 8, nChannel=8) :
    '''
    IN
    #####
    filename : path of the file
    maxFBGs : maximum number of FBGs on a channel
    nChannel : number maximium of channels 
    
    OUT
    ####
    array_data : data (3D array) extracted from the BSI files (Channel, FBG, Time)
    refTime : Time given by the interrogator, not recommanded to use
    absTime : Absolute value of Time in datetime64 array
    
    '''
    
    data_bsi = pd.read_csv(filename, header=None, skiprows=1, sep='\t', decimal=',',names=np.append(np.array(['Time', 'Channel']),np.arange(1,maxFBGs+1)))
    #array_data = np.zeros((nChannel, maxFBGs, data_bsi.shape[0]))
    dict_data, dict_time = dict(), dict()
     
    refTime = data_bsi.Time
    absTime = np.datetime64('%s-%s-%sT%s:%s:%s' %(filename[-24:-20], filename[-20:-18], filename[-18:-16],
                                                  filename[-15:-13],filename[-13:-11], filename[-11:-5].replace(',', '.')),'ms')
    #print(filename, refTime.shape)
    if refTime.shape[0] == 0 :
        warnings.warn('File %s is empty !' %(filename))
        for ch in range(1, nChannel+1) :
            dict_data['Channel_%i' %ch] = (['Channel %i : FBGs' %ch, "Channel %i : Over time" %ch], np.empty((maxFBGs,0)))
            dict_time['Channel_%i_time' %ch] = np.empty((0,))
        return dict_data, dict_time
        
    dtime = (refTime - refTime[0])*1000

    absTime = absTime+np.array(dtime, dtype=int)
    
    for ch in range(1,nChannel+1) :
        cond = data_bsi['Channel'] == ch
        if np.sum(cond) == 0:
            warnings.warn('No data or No existing channel %i' %ch)

        condMask = np.where(cond==True, 1, 0)
        for fbg in range(maxFBGs) :    
            data_ch =  np.extract(condMask, np.array(data_bsi)[:,2+fbg]) 
            medianfbg = np.nanmedian(data_ch)
            
            cleanfbg = np.where((data_ch > medianfbg + 1) | (data_ch < medianfbg - 1), np.nan, data_ch)
            if fbg == 0 :
                out_chann_time  = np.extract(condMask, absTime).astype('object')
                out_chann_arr = cleanfbg.reshape(1, len(cleanfbg))
                continue
            out_chann_arr = np.append(out_chann_arr, cleanfbg.reshape(1, len(cleanfbg)), axis=0)
                   
        dict_data['Channel_%i' %ch] = (['Channel %i : FBGs' %ch, "Channel %i : Over time" %ch], out_chann_arr)
        dict_time['Channel_%i_time' %ch] = out_chann_time
        
    return dict_data, dict_time


def folder_reader_DUAL_FBG(path, nChannel=8, maxFBGs=12) :
    files = os.listdir(path)
    files.sort()
    ji=0

    for file in files:
        if file[-5:] != 'Z.txt':
            continue
            
        dict_data, dict_time = bsi_hsv_dual_161_reader(path+file, maxFBGs=maxFBGs, nChannel=nChannel)
        mode = 'DUAL'
        
        if ji == 0 :
            ALL_dict_data, ALL_dict_time = dict_data.copy(), dict_time.copy()
            ji+=1 
            continue
            
        #print(file, TabsTime.shape, TrefTime.shape, TOutWav.shape)
        
        for i, name in enumerate(dict_data) :
            #print('Debug 1', name, ALL_dict_data[name][1].shape, dict_data[name][1].shape)

            ch = i+1
            name_time = name+'_time'
            
            ALL_dict_time[name_time] = np.append(ALL_dict_time[name_time], dict_time[name_time]) 
            out_chann_arr = np.append(ALL_dict_data[name][1], dict_data[name][1], axis=1)
            
            ALL_dict_data[name]  = (['Channel %i : FBGs' %ch, "Channel %i : Over time" %ch], out_chann_arr)  
            
            #print('Debug 2', name, ALL_dict_data[name][1].shape, dict_data[name][1].shape)

            del out_chann_arr
            
        if verif_date(ALL_dict_time['Channel_1_time']) == True :
            absTime = np.array(ALL_dict_time['Channel_1_time'], dtype='datetime64')
            IdxEndTime = np.where(absTime.astype('datetime64[D]')==absTime[-1].astype('datetime64[D]'))[0][0]

            OUT_dict_time = ALL_dict_time.copy()
            OUT_dict_time['Channel_1_time'] = ALL_dict_time['Channel_1_time'][:IdxEndTime].copy()
            OUT_dict_data = ALL_dict_data.copy()
            OUT_dict_data['Channel_1']  = (['Channel 1 : FBGs' , "Channel 1 : Over time"], ALL_dict_data['Channel_1'][1][:,:IdxEndTime])
            
            export_nc(OUT_dict_data, OUT_dict_time, mode=mode)
            
            del OUT_dict_data, OUT_dict_time
            
            out_chann_arr = ALL_dict_data['Channel_1'][1][:,IdxEndTime:].copy()
            ALL_dict_time['Channel_1_time'] = ALL_dict_time['Channel_1_time'][IdxEndTime:].copy()
            ALL_dict_data['Channel_1']  = (['Channel 1 : FBGs' , "Channel 1 : Over time"], out_chann_arr) 

            del out_chann_arr
            
            for i, name in enumerate(dict_data) :
                ch = i+1
                name_time = name+'_time'
                if ch == 1 :
                    continue
                ALL_dict_time[name_time] = np.array([])
                ALL_dict_data[name] = (['Channel %i : FBGs'%ch , "Channel %i : Over time" %ch], np.empty((maxFBGs,0)))             
    try :
        export_nc(ALL_dict_data, ALL_dict_time, mode=mode)#,refwav=refwav)
    except : pass