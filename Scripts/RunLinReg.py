import SVBfunc
VAR=['PHIHYD'] #['ashore','WVEL'] #['ashore','PHIHYD','WVEL']
FILT=['no','filt']
COAST=['smooth','original']
for coast in COAST:
	for var in VAR:
		for filt in FILT:
			SVBfunc.linearregressionSave(filt,var,coast)
			print('Done with ' + str(coast) + ' ' + str(var) + ' ' +  str(filt))

