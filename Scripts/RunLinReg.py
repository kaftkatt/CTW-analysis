import SVBfunc
VAR=['PHIHYD','WVEL'] # , ashore]
FILT=['no','filt']
#COAST=['smooth','original']
#for coast in COAST:

coast='smooth'
for var in VAR:
	for filt in FILT:
		SVBfunc.linearregressionSave(filt,var,coast)
		print('Done with ' + str(coast) + ' ' + str(var) + ' ' +  str(filt))

