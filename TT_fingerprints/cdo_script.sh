#!/bin/bash

list="HadISST_sst_masked sst.mnmean HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean_SST HadSST.4.0.1.0_median"
paths=($list)
lis="sst sst tas_mean tos"
vars=($lis)

### cdo weights means by gridcells automatically, so can compute rectangle area fingerprints in cdo, and also create the global mean for the exact SPG fingerprint (_gl.nc)
### create gridarea files to use for weighting SPG area in python script

for i in 0 1 2 3
do
path=${paths[${i}]}
var=${vars[${i}]}
echo ${path}
echo ${var}
cdo -L sub -fldmean -sellonlatbox,-55,-20,46,61 -selvar,${var} ${path}.nc  -fldmean -selvar,${var} ${path}.nc ${path}_M20.nc
cdo -L sub -fldmean -sellonlatbox,-70,-30,45,80 -selvar,${var} ${path}.nc -fldmean -sellonlatbox,-70,-30,-45,0 -selvar,${var} ${path}.nc ${path}_dipole.nc
cdo -L -fldmean -selvar,${var} ${path}.nc ${path}_gl.nc
cdo -L gridarea -selvar,${var} ${path}.nc ${path}_gridarea.nc
done

