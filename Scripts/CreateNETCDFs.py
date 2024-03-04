from SVBfunc import createNetCDF

for coast in ['smooth', 'straight']:
  for prefix in ['phiHyd', 'rhoAnoma', 'eta', 'dynVars']:
    createNetCDF(coast,prefix, prefix)
