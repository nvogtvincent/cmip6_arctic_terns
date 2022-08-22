#!/bin/bash
set -e

# Model name
INST_NAME=$1
MODEL_NAME=$2
INIT=$3

FULL_NAME=${INST_NAME}/${MODEL_NAME}

# Define directories
DIR_HIST=/badc/cmip6/data/CMIP6/CMIP/${FULL_NAME}/historical/${INIT}/day
DIR_SSP245=/badc/cmip6/data/CMIP6/ScenarioMIP/${FULL_NAME}/ssp245/${INIT}/day
DIR_SSP585=/badc/cmip6/data/CMIP6/ScenarioMIP/${FULL_NAME}/ssp585/${INIT}/day

# Define files
vars=('uas' 'vas')

# Make a temporary folder
TEMP_FOLDER=${MODEL_NAME}_processing
mkdir ${TEMP_FOLDER}

# Copy all the relevant files to the folder
module load jaspy

declare -A pids=()
cnt=$((0))

echo 'Copying files...'

cd ${TEMP_FOLDER}
for i in ${vars[@]}
do
    cp -Lr ${DIR_HIST}/${i}/gn/latest/*.nc .
    cp -Lr ${DIR_SSP245}/${i}/gn/latest/*.nc .
    cp -Lr ${DIR_SSP585}/${i}/gn/latest/*.nc .
done

echo 'Complete!'
echo ''

# Merge time series
declare -A pids=()
cnt=$((0))

echo 'Merging time series...'

for i in ${vars[@]}
    do
    cdo mergetime ${i}*historical*.nc ${i}_${MODEL_NAME}_historical.nc &
    cnt=$((${cnt}+1))
    pids[${cnt}]=$!

    cdo mergetime ${i}*ssp245*.nc ${i}_${MODEL_NAME}_ssp245.nc &
    cnt=$((${cnt}+1))
    pids[${cnt}]=$!

    cdo mergetime ${i}*ssp585*.nc ${i}_${MODEL_NAME}_ssp585.nc &
    cnt=$((${cnt}+1))
    pids[${cnt}]=$!
done

wait ${pids}

echo 'Merging complete!'
echo ''



# Extract surface if relevant
if [ ${vars[0]} = 'ua' ]; then
    declare -A pids=()
    cnt=$((0))
     
    echo 'Extracting surface level only...'  
 
    for i in ${vars[@]}
    do
        cdo sellevidx,1 ${i}_${MODEL_NAME}_historical.nc ${i}_historical.nc.1 &
        cnt=$((${cnt}+1))
        pids[${cnt}]=$!

        cdo sellevidx,1 ${i}_${MODEL_NAME}_ssp245.nc ${i}_ssp245.nc.1 &
        cnt=$((${cnt}+1))
        pids[${cnt}]=$!

        cdo sellevidx,1 ${i}_${MODEL_NAME}_ssp585.nc ${i}_ssp585.nc.1 &
        cnt=$((${cnt}+1))
        pids[${cnt}]=$!
    done

    wait ${pids}

    for i in ${vars[@]}
    do
        mv ${i}_historical.nc.1 ${i}_${MODEL_NAME}_historical.nc
        mv ${i}_ssp245.nc.1 ${i}_${MODEL_NAME}_ssp245.nc
        mv ${i}_ssp585.nc.1 ${i}_${MODEL_NAME}_ssp585.nc
    done

   echo 'Extraction complete!'
   echo ''
fi

echo 'Merging complete!'
echo ''

# Carry out remapping
cp ../../resources/gridsource.txt .

declare -A pids=()
for i in ${vars[@]}
do
    cdo remapbil,gridsource.txt ${i}_${MODEL_NAME}_historical.nc ${i}_${MODEL_NAME}_HISTORICAL.nc
    cnt=$((${cnt}+1))
    pids[${cnt}]=$!

    cdo remapbil,gridsource.txt ${i}_${MODEL_NAME}_ssp585.nc ${i}_${MODEL_NAME}_SSP585.nc
    cnt=$((${cnt}+1))
    pids[${cnt}]=$!

    cdo remapbil,gridsource.txt ${i}_${MODEL_NAME}_ssp245.nc ${i}_${MODEL_NAME}_SSP245.nc
    cnt=$((${cnt}+1))
    pids[${cnt}]=$!
done

wait ${pids}
echo 'Regridding complete!'
echo 'Cleaning up...'

# Move to final folder
cd ../
mkdir ${MODEL_NAME}
mv ${TEMP_FOLDER}/*HISTORICAL.nc ${MODEL_NAME}/
mv ${TEMP_FOLDER}/*SSP585.nc ${MODEL_NAME}/
mv ${TEMP_FOLDER}/*SSP245.nc ${MODEL_NAME}/

# Clean up the processing folder
rm -rf  ${TEMP_FOLDER}

echo 'Cleaning complete!'
