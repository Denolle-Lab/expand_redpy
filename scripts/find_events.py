import pandas as pd
import yaml
import numpy as np
import obspy
from obspy import UTCDateTime
from time import time
import csv
from glob import glob
import re

with open('/home/smocz/expand_redpy/scripts/config.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

vv = config['vv']
# vv=0
volc_list_names = config['volc_list_names']
volc = volc_list_names[vv]
print(volc)
year = config['year']
years = config['years']
# years = [2013,2014]
print(years)

homedir = config['homedir']
readdir = config['readdir']

minsta = config['minsta']

#time in seconds before and after a detection to check for similar detections
wi = config['wi']

print(minsta,wi)

volc_md = pd.read_csv(readdir+'Volcano_Metadata.csv') #read volcano metadata
volc_md['netsta'] = volc_md['Network'].astype(str)+'.'+volc_md['Station'].astype(str) 
#combine network and station to keep them associated

# for vv,v in enumerate(volc_list_names): #for each volcano
#     volc = v
    
t00 = time() #record time

# for year in years: #for each year
for year in years:
    t0 = time()

    #get detections at this volcano in this year
    sta_list = volc_md[volc_md['Volcano_Name']==volc]['netsta'].values.tolist() #get a list of stations at this volcano
    readstadict = {} #will become a dictionary of dataframes of station detections
    for i in sta_list: #for each station at the volcano
#         print(homedir+f'detections/{volc}_{i}_{year}_detections.csv')
        try: 
            readstadict[i] = pd.read_csv(homedir+f'detections/{volc}_{i}_{year}_detections.csv') 
            print(f'successfully read {i} in {year}')
        except: #if the detection final does NOT exist
            print(f'no detections for {i} in {year}') #say so
            readstadict[i] = pd.DataFrame() #make an empty data frame maintain index

    print('---')

    rsta_list = [] #list version of readstadict for making a dataframe
    for i in readstadict: #making rsta_list
        rsta_list.append(readstadict[i])
    readsta = pd.concat(rsta_list) #make a single dataframe for the whole volcano from the list of station dataframs
    print(type(rsta_list[0]))

    print('----')

    print(readsta)
    
    #for all detections we want to run (concatenated), make cl_list for the volcano
    temp_name_list = readsta['Template_Name'].values.tolist() #make a list of template names
    cl_list_long = [] # make a list of the numbers in each template name
    for i in temp_name_list: 
        if volc=='Baker' or volc=='Hood' or volc=='Newberry' or volc=='Rainier':
            num = i[-3:] #account for zfill, to parameterize, use what make_templates uses to get zfill amount (cid[-1] from volcmd)
        if volc=='St_Helens':
            num = i[-4:] #account for zfill
    #     num = re.findall(r'\d+', i)
        cl_list_long.append(num) #*num for findall but some stations have numbers in their name, so using zfill instead of 
        #just numbers
    cl_list = np.unique(cl_list_long) #get rid of duplicates

    csv_name = homedir+f'events/{volc}_{year}_events.csv'
    with open(csv_name, 'w', newline='') as file: #make a csv to save to
        writer = csv.writer(file)
        writer.writerow(["Earliest_Detection_Time","Cluster_ID","Stations_Found","Stations"]) #,"Stations_Diff"
        #detection time, cluster id, and number of stations with a detection for this event, what stations that is, ?
        file.close()

    for cl in cl_list: #for each cluster
        t1 = time() #record time
        times = [] #get list of datetimes for all templates for this cluster
        stas = [] #get list of station names for all templates for this cluster (same index as times)
        for i in np.unique(temp_name_list): #for each template
            if i.endswith(cl): #if the template ends with the cluster ID
                all_times = readsta[readsta['Template_Name']==i]['Detection_Time'].values.tolist() #find the times for that template
                temp_sta = i[:-10].upper()
                for at in all_times:
                    times.append(at) #append each template time to the list of datetimes
                    stas.append(temp_sta) #append the station to the list
                    
        for ii,i in enumerate(times): #for all detections
            t3 = time()
    #run through each detection time on all stations
            print( )
            print('------------')
            print('trying detection',ii,'cluster',cl,i)

            times_diff = []
            match_list = []
            for ss,s in enumerate(sta_list): #for each station
                
                if stas[ii]==s.split('.')[1]: #if this station is the same one the detection came from
                    match=2
                    print('Auto overlap with station '+s+' detections')
                    
                else:
                
                    statimes = [] #get the times for this station
                    try:
                        readstadict[s]['Template_Name']
                        print(f'successfully read {s} in {year}')
                    except:
                        print(f'no data available for {s} in {year}')
                        match_list.append(0)
                        continue
                    for a in np.unique(readstadict[s]['Template_Name'].values.tolist()):
                        # for every template in this station's df
                        if a.endswith(cl): #if it matches this cluster
                            all_times = readstadict[s][readstadict[s]['Template_Name']==a]['Detection_Time'].values.tolist() 
                            #get a list of detection times at this station for this cluster
                            for at in all_times:
                                statimes.append(at) #append to statimes
                    match=0 #set variable to arbitrary number
                    for tt,t in enumerate(statimes): #for every datetime in the station's detections for this cluster
                        ts = UTCDateTime(t)-wi #find the time wi s before
                        te = UTCDateTime(t)+wi #find the time wi s after
                        if UTCDateTime(i)>ts and UTCDateTime(i)<te:
                            diff = UTCDateTime(i)-UTCDateTime(t)
                            match=2 #if there is an overlap, reset the variable and break out of the loop
                            print('Overlap with station '+s+' detections')
                            break
    #                     times_diff.append(diff)
                match_list.append(match)
            print('match_list:',match_list)
            
    

            t2 = time()
            print(t2-t3,'seconds to test overlap for detection',ii)

            save_list = [] #for saving a list of stations
            for mm,m in enumerate(match_list):
                if m==2:
                    save_list.append(sta_list[mm])
                    
            print(save_list)

            if match_list.count(2) >= minsta: #if at least 4 matches equal 2:
                print('saving...')
                check = pd.read_csv(csv_name) #read the csv we made
                checktimes = check[check['Cluster_ID']==int(cl)]['Earliest_Detection_Time'].values.tolist() 
                #find the datetimes that have already been saved
                checking = 0
                for tt,t in enumerate(checktimes):
                    ts = UTCDateTime(t)-wi
                    te = UTCDateTime(t)+wi
                    if UTCDateTime(i)>ts and UTCDateTime(i)<te: 
                        print('already recorded event')
                        checking=1 #checking still =1 because we don't want the normal save process
                        #if i is earlier (less than) than the recorded event, overwrite the event (we want the earliest detected time)
                        if i<t: #if this time is earlier than the saved time
                            print('replacing...')
                            print('old time:',t,'new time:',i) 
                            check.replace(t,i,inplace=True) #replace t with i in the df
                            check.to_csv(csv_name,index=False)
                            #read panda dataframe - change dataframe, overwrite the csv with this dataframe
                            print('changed csv')
                        break
                if checking==1: continue # because it has already been saved
                row = [i,int(cl),match_list.count(2),' '.join(save_list),] #put together the row, #, ' '.join(times_diff)
                print(row)

                with open(csv_name, 'a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(row) #save the row to csv
                    file.close()
        t4 = time()
        print(t4-t1,'seconds for cluster',cl)
    t5 = time()
    print(t5-t0,'seconds to find all new events on',volc,'in',year)
print(t5-t00,'seconds to find all new events on',volc,'in',years)