#!/bin/bash -l

### ----------------------------- ###
### Downloading models selected in Kwiatkowski et al 2020 Biogeosciences
### ----------------------------- ###

cd /mnt/lustre/users/pearseb/hackathon_productivity/


#### ACCESS-ESM1-5
#wget -nc http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/intpp/gn/v20191115/intpp_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc intpp_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_185001-201412
#for year in $(seq 1850 1 2014); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_185001-201412.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_185001-201412_yearmonths.nc
#rm intpp_*.nc
#
#wget -nc http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/ScenarioMIP/CSIRO/ACCESS-ESM1-5/ssp585/r1i1p1f1/Omon/intpp/gn/v20210318/intpp_Omon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_201501-210012.nc
#wget -nc http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/ScenarioMIP/CSIRO/ACCESS-ESM1-5/ssp585/r1i1p1f1/Omon/intpp/gn/v20210318/intpp_Omon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_210101-230012.nc
#ncrcat -O intpp_Omon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_201501-210012.nc intpp_Omon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_210101-230012.nc intpp_Omon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_201501-230012.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_201501-230012.nc intpp_Omon_ACCESS-ESM1-5_ssp585_r1i1p1f1_201501-230012
#for year in $(seq 2015 1 2300); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_ACCESS-ESM1-5_ssp585_r1i1p1f1_201501-230012.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_ACCESS-ESM1-5_ssp585_r1i1p1f1_201501-230012_yearmonths.nc
#rm intpp_*.nc
#
#
#### CanESM
#wget -nc http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgA_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Omon/intpp/gn/v20190429/intpp_Omon_CanESM5_historical_r1i1p2f1_gn_185001-201412.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_CanESM5_historical_r1i1p2f1_gn_185001-201412.nc intpp_Omon_CanESM5_historical_r1i1p2f1_185001-201412
#for year in $(seq 1850 1 2014); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_CanESM5_historical_r1i1p2f1_185001-201412.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_CanESM5_historical_r1i1p2f1_185001-201412_yearmonths.nc
#rm intpp_*.nc
#
#
#wget -nc http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgA_dataroot/AR6/CMIP6/ScenarioMIP/CCCma/CanESM5/ssp585/r1i1p2f1/Omon/intpp/gn/v20190429/intpp_Omon_CanESM5_ssp585_r1i1p2f1_gn_201501-210012.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_CanESM5_ssp585_r1i1p2f1_gn_201501-210012.nc intpp_Omon_CanESM5_ssp585_r1i1p2f1_201501-210012
#for year in $(seq 2015 1 2100); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_CanESM5_ssp585_r1i1p2f1_201501-210012.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_CanESM5_ssp585_r1i1p2f1_201501-210012_yearmonths.nc
#rm intpp_*.nc
#
#
#### CESM2
#wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/NCAR/CESM2/historical/r4i1p1f1/Omon/intpp/gn/v20190308/intpp_Omon_CESM2_historical_r4i1p1f1_gn_185001-201412.nc 
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_CESM2_historical_r4i1p1f1_gn_185001-201412.nc intpp_Omon_CESM2_historical_r4i1p1f1_185001-201412
#for year in $(seq 1850 1 2014); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_CESM2_historical_r4i1p1f1_185001-201412.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_CESM2_historical_r4i1p1f1_185001-201412_yearmonths.nc
#rm intpp_*.nc
#
#wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCAR/CESM2/ssp585/r4i1p1f1/Omon/intpp/gn/v20200528/intpp_Omon_CESM2_ssp585_r4i1p1f1_gn_201501-206412.nc 
#wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCAR/CESM2/ssp585/r4i1p1f1/Omon/intpp/gn/v20200528/intpp_Omon_CESM2_ssp585_r4i1p1f1_gn_206501-210012.nc
#ncrcat -O intpp_Omon_CESM2_ssp585_r4i1p1f1_gn_201501-206412.nc intpp_Omon_CESM2_ssp585_r4i1p1f1_gn_206501-210012.nc intpp_Omon_CESM2_ssp585_r4i1p1f1_gn_201501-210012.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_CESM2_ssp585_r4i1p1f1_gn_201501-210012.nc intpp_Omon_CESM2_ssp585_r4i1p1f1_201501-210012
#for year in $(seq 2015 1 2100); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_CESM2_ssp585_r4i1p1f1_201501-210012.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_CESM2_ssp585_r4i1p1f1_201501-210012_yearmonths.nc
#rm intpp_*.nc
#
#
#### CNRM-ESM2
#wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/CNRM-CERFACS/CNRM-ESM2-1/historical/r1i1p1f2/Omon/intpp/gn/v20181206/intpp_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_185001-201412.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_185001-201412.nc intpp_Omon_CNRM-ESM2-1_historical_r1i1p1f2_185001-201412
#for year in $(seq 1850 1 2014); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_CNRM-ESM2-1_historical_r1i1p1f2_185001-201412.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_CNRM-ESM2-1_historical_r1i1p1f2_185001-201412_yearmonths.nc
#rm intpp_*.nc
#
#wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/CNRM-CERFACS/CNRM-ESM2-1/ssp585/r1i1p1f2/Omon/intpp/gn/v20191021/intpp_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_201501-210012.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_201501-210012.nc intpp_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_201501-210012
#for year in $(seq 2015 1 2100); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_201501-210012.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_201501-210012_yearmonths.nc
#rm intpp_*.nc
#
#
#### GFDL-CM4
#wget -nc http://aims3.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_185001-186912.nc
#wget -nc http://aims3.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_187001-188912.nc
#wget -nc http://aims3.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_189001-190912.nc
#wget -nc http://aims3.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_191001-192912.nc
#wget -nc http://aims3.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_193001-194912.nc
#wget -nc http://aims3.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_195001-196912.nc
#wget -nc http://aims3.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_197001-198912.nc
#wget -nc http://aims3.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_199001-200912.nc
#wget -nc http://aims3.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_201001-201412.nc
#files=`ls intpp_*.nc`
#ncrcat -O $files intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_185001-201412.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_GFDL-CM4_historical_r1i1p1f1_gr_185001-201412.nc intpp_Omon_GFDL-CM4_historical_r1i1p1f1_185001-201412
#for year in $(seq 1850 1 2014); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_GFDL-CM4_historical_r1i1p1f1_185001-201412.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O $files ETOPO_intpp_Omon_GFDL-CM4_historical_r1i1p1f1_185001-201412_yearmonths.nc
#rm intpp_*.nc
#
#wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-CM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_201501-203412.nc
#wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-CM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_203501-205412.nc
#wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-CM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_205501-207412.nc
#wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-CM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_207501-209412.nc
#wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-CM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_209501-210012.nc
#ncrcat -O intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_201501-203412.nc intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_203501-205412.nc intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_205501-207412.nc intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_207501-209412.nc intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_209501-210012.nc intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_201501-210012.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_gr_201501-210012.nc intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_201501-210012
#for year in $(seq 2015 1 2100); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_201501-210012.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_GFDL-CM4_ssp585_r1i1p1f1_201501-210012_yearmonths.nc
#rm intpp_*.nc
#
#
#### GFDL-ESM4
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/intpp/gr/v20190726/intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_185001-186912.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/intpp/gr/v20190726/intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_187001-188912.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/intpp/gr/v20190726/intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_189001-190912.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/intpp/gr/v20190726/intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_191001-192912.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/intpp/gr/v20190726/intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_193001-194912.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/intpp/gr/v20190726/intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_195001-196912.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/intpp/gr/v20190726/intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_197001-198912.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/intpp/gr/v20190726/intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_199001-200912.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/intpp/gr/v20190726/intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_201001-201412.nc
#files=`ls intpp_*.nc`
#ncrcat -O $files intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_185001-201412.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_185001-201412.nc intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_185001-201412
#for year in $(seq 1850 1 2014); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_185001-201412.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O $files ETOPO_intpp_Omon_GFDL-ESM4_historical_r1i1p1f1_185001-201412_yearmonths.nc
#rm intpp_*.nc
#
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_201501-203412.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_203501-205412.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_205501-207412.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_207501-209412.nc
#wget -nc http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/intpp/gr/v20180701/intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_209501-210012.nc
#ncrcat -O intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_201501-203412.nc intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_203501-205412.nc intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_205501-207412.nc intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_207501-209412.nc intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_209501-210012.nc intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_201501-210012.nc
#source /users/pearseb/load_cdo.env
#/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_201501-210012.nc intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_201501-210012
#for year in $(seq 2015 1 2100); do
# echo ${year}
# cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_201501-210012.nc intpp_${year}.nc
#done
#source /users/pearseb/purge.env
#files=`ls intpp_????.nc`
#ncecat -O ${files} ETOPO_intpp_Omon_GFDL-ESM4_ssp585_r1i1p1f1_201501-210012_yearmonths.nc
#rm intpp_*.nc


### IPSL-CM6A-LR
wget -nc http://vesg.ipsl.upmc.fr/thredds/fileServer/cmip6/CMIP/IPSL/IPSL-CM6A-LR/historical/r1i1p1f1/Omon/intpp/gn/v20180803/intpp_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc intpp_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_185001-201412
for year in $(seq 1850 1 2014); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_185001-201412.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O $files ETOPO_intpp_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_185001-201412_yearmonths.nc
rm intpp_*.nc

wget -nc http://vesg.ipsl.upmc.fr/thredds/fileServer/cmip6/ScenarioMIP/IPSL/IPSL-CM6A-LR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190903/intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_201501-210012.nc
wget -nc http://vesg.ipsl.upmc.fr/thredds/fileServer/cmip6/ScenarioMIP/IPSL/IPSL-CM6A-LR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190903/intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_210101-230012.nc
ncrcat -O intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_201501-210012.nc intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_210101-230012.nc intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_201501-230012.nc
ncks -O -C -x -v area intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_201501-230012.nc intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_201501-230012.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_201501-230012.nc intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_201501-230012
for year in $(seq 2015 1 2300); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_201501-230012.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O ${files} ETOPO_intpp_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_201501-230012_yearmonths.nc
rm intpp_*.nc


### MIROC-ES2L
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/Omon/intpp/gn/v20191129/intpp_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc intpp_Omon_MIROC-ES2L_historical_r1i1p1f2_185001-201412
for year in $(seq 1850 1 2014); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_MIROC-ES2L_historical_r1i1p1f2_185001-201412.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O $files ETOPO_intpp_Omon_MIROC-ES2L_historical_r1i1p1f2_185001-201412_yearmonths.nc
rm intpp_*.nc

wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/MIROC/MIROC-ES2L/ssp585/r1i1p1f2/Omon/intpp/gn/v20191129/intpp_Omon_MIROC-ES2L_ssp585_r1i1p1f2_gn_201501-210012.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_MIROC-ES2L_ssp585_r1i1p1f2_gn_201501-210012.nc intpp_Omon_MIROC-ES2L_ssp585_r1i1p1f2_201501-210012
for year in $(seq 2015 1 2100); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_MIROC-ES2L_ssp585_r1i1p1f2_201501-210012.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O ${files} ETOPO_intpp_Omon_MIROC-ES2L_ssp585_r1i1p1f2_201501-210012_yearmonths.nc
rm intpp_*.nc


### MPI-ESM1.2-HR
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_185001-185412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_185501-185912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_186001-186412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_186501-186912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_187001-187412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_187501-187912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_188001-188412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_188501-188912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_189001-189412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_189501-189912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_190001-190412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_190501-190912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_191001-191412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_191501-191912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_192001-192412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_192501-192912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_193001-193412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_193501-193912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_194001-194412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_194501-194912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_195001-195412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_195501-195912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_196001-196412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_196501-196912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197001-197412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_197501-197912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_198001-198412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_198501-198912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_199001-199412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_199501-199912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200001-200412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200501-200912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls intpp_*.nc`
ncrcat -O $files intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_185001-201412.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_185001-201412.nc intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_185001-201412 
for year in $(seq 1850 1 2014); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_185001-201412.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O $files ETOPO_intpp_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_185001-201412_yearmonths.nc
rm intpp_*.nc

wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_201501-201912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_202001-202412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_202501-202912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_203001-203412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_203501-203912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_204001-204412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_204501-204912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_205001-205412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_205501-205912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_206001-206412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_206501-206912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_207001-207412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_207501-207912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_208001-208412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_208501-208912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_209001-209412.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_209501-209912.nc
wget -nc http://esgf3.dkrz.de/thredds/fileServer/cmip6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp585/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_210001-210012.nc
files=`ls intpp_*.nc`
ncrcat -O $files intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_201501-210012.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_gn_201501-210012.nc intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_201501-210012
for year in $(seq 2015 1 2100); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_201501-210012.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O ${files} ETOPO_intpp_Omon_MPI-ESM1-2-HR_ssp585_r1i1p1f1_201501-210012_yearmonths.nc
rm intpp_*.nc


### MRI-ESM2
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i2p1f1/Omon/intpp/gn/v20210311/intpp_Omon_MRI-ESM2-0_historical_r1i2p1f1_gn_185001-201412.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_MRI-ESM2-0_historical_r1i2p1f1_gn_185001-201412.nc intpp_Omon_MRI-ESM2-0_historical_r1i2p1f1_185001-201412 
for year in $(seq 1850 1 2014); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_MRI-ESM2-0_historical_r1i2p1f1_185001-201412.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O $files ETOPO_intpp_Omon_MRI-ESM2-0_historical_r1i2p1f1_185001-201412_yearmonths.nc
rm intpp_*.nc

wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i2p1f1/Omon/intpp/gn/v20210329/intpp_Omon_MRI-ESM2-0_ssp585_r1i2p1f1_gn_201501-210012.nc 
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_MRI-ESM2-0_ssp585_r1i2p1f1_gn_201501-210012.nc intpp_Omon_MRI-ESM2-0_ssp585_r1i2p1f1_201501-210012
for year in $(seq 2015 1 2100); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_MRI-ESM2-0_ssp585_r1i2p1f1_201501-210012.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O ${files} ETOPO_intpp_Omon_MRI-ESM2-0_ssp585_r1i2p1f1_201501-210012_yearmonths.nc
rm intpp_*.nc


### NorESM-LM
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_185001-185912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_186001-186912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_187001-187912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_188001-188912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_189001-189912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_190001-190912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_191001-191912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_192001-192912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_193001-193912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_194001-194912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_195001-195912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_196001-196912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_197001-197912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_198001-198912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_199001-199912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_200001-200912.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls intpp_*.nc`
ncrcat -O $files intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc intpp_Omon_NorESM2-LM_historical_r1i1p1f1_185001-201412 
for year in $(seq 1850 1 2014); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_NorESM2-LM_historical_r1i1p1f1_185001-201412.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O $files ETOPO_intpp_Omon_NorESM2-LM_historical_r1i1p1f1_185001-201412_yearmonths.nc
rm intpp_*.nc


wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCC/NorESM2-LM/ssp585/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_201501-202012.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCC/NorESM2-LM/ssp585/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_202101-203012.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCC/NorESM2-LM/ssp585/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_203101-204012.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCC/NorESM2-LM/ssp585/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_204101-205012.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCC/NorESM2-LM/ssp585/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_205101-206012.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCC/NorESM2-LM/ssp585/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_206101-207012.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCC/NorESM2-LM/ssp585/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_207101-208012.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCC/NorESM2-LM/ssp585/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_208101-209012.nc
wget -nc https://esgf-data1.llnl.gov/thredds/fileServer/css03_data/CMIP6/ScenarioMIP/NCC/NorESM2-LM/ssp585/r1i1p1f1/Omon/intpp/gn/v20191108/intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_209101-210012.nc
files=`ls intpp_*.nc`
ncrcat -O $files intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_201501-210012.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_gn_201501-210012.nc intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_201501-210012
for year in $(seq 2015 1 2100); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_201501-210012.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O ${files} ETOPO_intpp_Omon_NorESM2-LM_ssp585_r1i1p1f1_201501-210012_yearmonths.nc
rm intpp_*.nc


### UKESM-1-0-LL
wget -nc http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/intpp/gn/v20190627/intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_185001-194912.nc
wget -nc http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/intpp/gn/v20190627/intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc
files=`ls intpp_*.nc`
ncrcat -O $files intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_185001-201412.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_185001-201412.nc intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_185001-201412 
for year in $(seq 1850 1 2014); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_185001-201412.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O $files ETOPO_intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_185001-201412_yearmonths.nc
rm intpp_*.nc

wget -nc http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/ScenarioMIP/MOHC/UKESM1-0-LL/ssp585/r1i1p1f2/Omon/intpp/gn/v20190726/intpp_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_201501-204912.nc
wget -nc http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/ScenarioMIP/MOHC/UKESM1-0-LL/ssp585/r1i1p1f2/Omon/intpp/gn/v20190726/intpp_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_205001-210012.nc
files=`ls intpp_*.nc`
ncrcat -O $files intpp_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_201501-210012.nc
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/regrid_to_1x1.sh intpp_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_201501-210012.nc intpp_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_201501-210012
for year in $(seq 2015 1 2100); do
 echo ${year}
 cdo -O -f nc -selyear,${year} ETOPO_intpp_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_201501-210012.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncecat -O ${files} ETOPO_intpp_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_201501-210012_yearmonths.nc
rm intpp_*.nc
