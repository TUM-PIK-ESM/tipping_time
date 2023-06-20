#!/bin/bash
git pull
# git add *
git add . :!TT_fingerprints/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean_SST.nc :!TT_fingerprints/HadISST_sst_masked.nc :!HadSST.4.0.1.0_median.nc :!sst.mnmean.nc
git commit -m "upload fingerprints, code and results"
git push
