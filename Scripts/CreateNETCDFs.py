from SVBfunc import createNetCDF

#for coast in ['smooth', 'straight']:
coast='originial'
for prefix in ['phiHyd', 'rhoAnoma', 'eta', 'dynVars']:
        createNetCDF(coast,prefix, prefix)
