#station names
sta0 = 'VDB'
sta1 = 'MBW'
sta2 = 'MULN'
sta3 = 'PASS'
sta4 = 'SAXON'
sta5 = 'SHUK'
#panda dataframe for detections at a station
readsta0 = pd.read_csv('/home/smocz/expand_redpy_new_files/detections/Baker_VDB_2019_clean_detections.csv')
readsta1 = pd.read_csv('/home/smocz/expand_redpy_new_files/detections/Baker_MBW_2019_clean_detections.csv')
readsta2 = pd.read_csv('/home/smocz/expand_redpy_new_files/detections/Baker_MULN_2019_clean_detections.csv')
readsta3 = pd.read_csv('/home/smocz/expand_redpy_new_files/detections/Baker_PASS_2019_clean_detections.csv')
readsta4 = pd.read_csv('/home/smocz/expand_redpy_new_files/detections/Baker_SAXON_2019_clean_detections.csv')
readsta5 = pd.read_csv('/home/smocz/expand_redpy_new_files/detections/Baker_SHUK_2019_clean_detections.csv')
#list of stations for saving (must be in order)
sta_list = ['sta0', 'sta1', 'sta2', 'sta3', 'sta4', 'sta5']
#concat of station dataframes for running through all detections
readsta = pd.concat(['readsta0', 'readsta1', 'readsta2', 'readsta3', 'readsta4', 'readsta5'])
