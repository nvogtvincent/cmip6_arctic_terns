# Global Warming and Arctic Terns: Predicting Climate Change Impacts on the World's Longest Migrant
This is the github repository for scripts required to reproduce figures and analyses in Morten et al. (202X) *Global Warming and Arctic Terns: Predicting Climate Change Impacts on the World's Longest Migrant*.

## Repository structure:
```
cmip6_arctic_terns/
├─ seaice/
│  ├─ DATA/
│  │  ├─ Collated migration data.csv | Arctic Tern tracking data
│  │  ├─ ETOPO_siconc_SImon_UKESM1-  | Historical sea ice concentration
│  │  │	 0-LL_historical_r1i1p1f2_
│  │  │	 198201-201411_SH.nc	    
│  │  ├─ ETOPO_siconc_SImon_UKESM1-  | SSP5-8.5 sea ice concentration
│  │  │	 0-LL_ssp585_r1i1p1f2_
│  │  │	 201501-202101_SH.nc	   
│  │  ├─ icec.mnmean.nc 	     | Observed sea ice concentration
│  ├─ FIGURES/
│  │  ├─ fig_7.pdf 		     | Figure 7
│  │  ├─ sup_fig_5.pdf 		     | Supp. Figure 5
│  │  ├─ sup_fig_6.pdf 		     | Supp. Figure 6
│  │  ├─ sup_fig_7.pdf 		     | Supp. Figure 7
│  │  ├─ sup_fig_8.pdf 		     | Supp. Figure 8
│  ├─ SCRIPTS/
│  │  ├─ observed_sie_winter.py	     | Generates Supplementary Figure 7
│  │  ├─ observed_vs_modelled_sie.py | Generates Supplementary Figure 6
│  │  ├─ productivity_change.py      | Generates Figure 7
│  │  ├─ scatter_trends_prod.py	     | Generates Supplmenetary Figure 8
│  │  ├─ terns_at_iceedge.py	     | Generates Supplementary Figure 5
├─ ternsim/                      
│  ├─ PREPROCESSING.sh               | Shell script used to preprocess raw CMIP6 output for particle tracking
│  ├─ TernSim.py                     | Simulates vTern trajectories from preprocessed CMIP6 surface velocity data
│  ├─ flight_path.py 	             | Generates Figure 4
│  ├─ flight_time.py 	             | Generates Supplementary Figure 3
│  ├─ path_comp.py                   | Generates Supplementary Figure 2
│  ├─ ternmethods.py                 | Methods for TernSim.py   
│  ├─ FIGURES/
│  │  ├─ vTern_path_future.pdf       | Figure 4
│  │  ├─ vTern_path_validation.pdf   | Supp. Figure 2
│  │  ├─ vTern_time_future.pdf       | Supp. Figure 3
│  ├─ resources/
│  │  ├─ gridsource.txt              | Regridding file used by PREPROCESSING.sh
├─ winds/
│  ├─ DATA/
│  │  ├─ bird_heatmap.nc		         | Heatmap of Tern locations
│  │  ├─ ens_uas_hist.nc             | Historical zonal winds
│  │  ├─ ens_uas_s245.nc             | SSP2-4.5 zonal winds
│  │  ├─ ens_uas_s585.nc             | SSP5-8.5 zonal winds
│  │  ├─ ens_vas_hist.nc             | Historical meridional winds
│  │  ├─ ens_vas_s245.nc             | SSP2-4.5 meridional winds
│  │  ├─ ens_vas_s585.nc             | SSP5-8.5 meridional winds
│  ├─ FIGURES/
│  │  ├─ ens_wind_zonal_nh.pdf       | Figure 3 (left panel)
│  │  ├─ ens_wind_zonal_sh.pdf       | Figure 3 (right panel)
│  ├─ zonal_winds.py                 | Generates Figure 3
```


