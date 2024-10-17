import SVBfunc
VAR=['ashore','PHIHYD','WVEL']
FILT=['no','filt']
#COAST=['smooth','original']
#for coast in COAST:

coast='original'
for var in VAR:
	for filt in FILT:
		SVBfunc.linearregressionSave(filt,var,coast)
		print('Done with ' + str(coast) + ' ' + str(var) + ' ' +  str(filt))

