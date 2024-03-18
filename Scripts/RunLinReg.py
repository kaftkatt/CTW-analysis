import SVBfunc
VAR=['PHIHYD','VVEL','UVEL','WVEL']
FILT=['no','filt']
for var in VAR:
  var='WVEL'
  for filt in FILT:
	  SVBfunc.linearregressionSave(filt,var)
