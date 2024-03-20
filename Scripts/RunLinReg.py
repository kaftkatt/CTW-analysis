import SVBfunc
VAR=['PHIHYD','VVEL','UVEL','WVEL']
FILT=['no','filt']
for var in VAR:
  for filt in FILT:
	  SVBfunc.linearregressionSave(filt,var)
