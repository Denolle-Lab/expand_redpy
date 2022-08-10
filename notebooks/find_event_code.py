sta0times = []
for a in np.unique(readsta0['Template_Name'].values.tolist()):
	if a.endswith(cl):
		all_times = readsta0[readsta0['Template_Name']==a['Detection_Time'].values.tolist
		for at in all_times:
			sta0times.append(at)
match0=0 #set variable to arbitrary number
for tt,t in enumerate(sta0times):
	ts = UTCDateTime(t)-wi
	te = UTCDateTime(t)+wi
	if UTCDateTime(i)>ts and UTCDateTime(i)<te:
		match0=2 #if there is an overlap, reset the variable and break out of the loop
		print('Overlap with station '+sta0+' detections')
 
sta1times = []
for a in np.unique(readsta1['Template_Name'].values.tolist()):
	if a.endswith(cl):
		all_times = readsta1[readsta1['Template_Name']==a['Detection_Time'].values.tolist
		for at in all_times:
			sta1times.append(at)
match1=0 #set variable to arbitrary number
for tt,t in enumerate(sta1times):
	ts = UTCDateTime(t)-wi
	te = UTCDateTime(t)+wi
	if UTCDateTime(i)>ts and UTCDateTime(i)<te:
		match1=2 #if there is an overlap, reset the variable and break out of the loop
		print('Overlap with station '+sta1+' detections')
 
sta2times = []
for a in np.unique(readsta2['Template_Name'].values.tolist()):
	if a.endswith(cl):
		all_times = readsta2[readsta2['Template_Name']==a['Detection_Time'].values.tolist
		for at in all_times:
			sta2times.append(at)
match2=0 #set variable to arbitrary number
for tt,t in enumerate(sta2times):
	ts = UTCDateTime(t)-wi
	te = UTCDateTime(t)+wi
	if UTCDateTime(i)>ts and UTCDateTime(i)<te:
		match2=2 #if there is an overlap, reset the variable and break out of the loop
		print('Overlap with station '+sta2+' detections')
 
sta3times = []
for a in np.unique(readsta3['Template_Name'].values.tolist()):
	if a.endswith(cl):
		all_times = readsta3[readsta3['Template_Name']==a['Detection_Time'].values.tolist
		for at in all_times:
			sta3times.append(at)
match3=0 #set variable to arbitrary number
for tt,t in enumerate(sta3times):
	ts = UTCDateTime(t)-wi
	te = UTCDateTime(t)+wi
	if UTCDateTime(i)>ts and UTCDateTime(i)<te:
		match3=2 #if there is an overlap, reset the variable and break out of the loop
		print('Overlap with station '+sta3+' detections')
 
sta4times = []
for a in np.unique(readsta4['Template_Name'].values.tolist()):
	if a.endswith(cl):
		all_times = readsta4[readsta4['Template_Name']==a['Detection_Time'].values.tolist
		for at in all_times:
			sta4times.append(at)
match4=0 #set variable to arbitrary number
for tt,t in enumerate(sta4times):
	ts = UTCDateTime(t)-wi
	te = UTCDateTime(t)+wi
	if UTCDateTime(i)>ts and UTCDateTime(i)<te:
		match4=2 #if there is an overlap, reset the variable and break out of the loop
		print('Overlap with station '+sta4+' detections')
 
sta5times = []
for a in np.unique(readsta5['Template_Name'].values.tolist()):
	if a.endswith(cl):
		all_times = readsta5[readsta5['Template_Name']==a['Detection_Time'].values.tolist
		for at in all_times:
			sta5times.append(at)
match5=0 #set variable to arbitrary number
for tt,t in enumerate(sta5times):
	ts = UTCDateTime(t)-wi
	te = UTCDateTime(t)+wi
	if UTCDateTime(i)>ts and UTCDateTime(i)<te:
		match5=2 #if there is an overlap, reset the variable and break out of the loop
		print('Overlap with station '+sta5+' detections')
 
match_list = ['match0', 'match1', 'match2', 'match3', 'match4', 'match5']
