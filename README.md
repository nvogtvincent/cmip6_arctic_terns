# CMIP6 Arctic Terns (replace with actual title when we choose it...)
This is the github repository for scripts required to reproduce figures and analyses in Morten et al. (202X) [Title].

## Repository structure:
```
cmip6_arctic_terns/
├─ ternsim/                      
│  ├─ PREPROCESSING.sh              | Shell script used to preprocess raw CMIP6 output for particle tracking
│  ├─ TernSim.py                    | Simulates vTern trajectories from preprocessed CMIP6 surface velocity data
│  ├─ flight_path.py 	            | Generates Figure 4
│  ├─ flight_time.py 	            | Generates Supplementary Figure 3
│  ├─ path_comp.py                  | Generates Supplementary Figure 2
│  ├─ ternmethods.py                | Methods for TernSim.py   
│  ├─ FIGURES/
│  │  ├─ vTern_path_future.pdf      | Figure 4
│  │  ├─ vTern_path_validation.pdf  | Supp. Figure 2
│  │  ├─ vTern_time_future.pdf      | Supp. Figure 3
│  ├─ resources/
│  │  ├─ gridsource.txt             | Regridding file used by PREPROCESSING.sh
├─ winds/
│  ├─ model_ensembler.py            | Generates ensemble mean data from CMIP6 catalogue
│  ├─ zonal_winds.ipynb             | Generates Figure 3
│  ├─ DATA/
│  │  ├─ ens_uas_hist.nc            | Historical zonal winds
│  │  ├─ ens_uas_s245.nc            | SSP2-4.5 zonal winds
│  │  ├─ ens_uas_s585.nc            | SSP5-8.5 zonal winds
│  │  ├─ ens_vas_hist.nc            | Historical meridional winds
│  │  ├─ ens_vas_s245.nc            | SSP2-4.5 meridional winds
│  │  ├─ ens_vas_s585.nc            | SSP5-8.5 meridional winds
│  ├─ FIGURES/
│  │  ├─ ens_wind_zonal_nh.pdf      | Figure 3 (left panel)
│  │  ├─ ens_wind_zonal_sh.pdf      | Figure 3 (right panel)
```


