import SVBfunc
VAR=['PHIHYD','VVEL','UVEL','WVEL']
FILT=['no','filt']
#COAST=['smooth','orginal']
#for coast in COAST:
coast='original'
for var in VAR:
	for filt in FILT:
		SVBfunc.linearregressionSave(filt,var,coast)
		print('Done with ' + str(coast) + ' ' + str(var) + ' ' +  str(filt))

