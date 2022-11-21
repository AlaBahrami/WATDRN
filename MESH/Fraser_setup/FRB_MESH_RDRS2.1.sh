#!/bin/bash

# foring, config and run dirs
dir_forcing=./Forcing/RDRS2.1/200009
config=./config
dir_mesh=../MESH_code/r1813_modified

ln -sf $dir_forcing/Fraser_RDRS_v2.1_A_PR0_SFC_200009_MESH.nc    ./basin_rain.nc      
ln -sf $dir_forcing/Fraser_RDRS_v2.1_P_UVC_09944_200009_MESH.nc  ./basin_wind.nc
ln -sf $dir_forcing/Fraser_RDRS_v2.1_P_HU_09944_200009_MESH.nc   ./basin_humidity.nc
ln -sf $dir_forcing/Fraser_RDRS_v2.1_P_FB_SFC_200009_MESH.nc     ./basin_shortwave.nc
ln -sf $dir_forcing/Fraser_RDRS_v2.1_P_FI_SFC_200009_MESH.nc     ./basin_longwave.nc
ln -sf $dir_forcing/Fraser_RDRS_v2.1_P_TT_09944_200009_MESH.nc   ./basin_temperature.nc
ln -sf $dir_forcing/Fraser_RDRS_v2.1_P_P0_SFC_200009_MESH.nc     ./basin_pres.nc

# config
ln -sf $config/MESH_input_streamflow.txt 								./MESH_input_streamflow.txt
ln -sf $config/MESH_drainage_database.r2c 							    ./MESH_drainage_database.r2c
ln -sf $config/MESH_input_reservoir.txt								    ./MESH_input_reservoir.txt
ln -sf $config/MESH_input_soil_levels.txt								./MESH_input_soil_levels.txt
ln -sf $config/MESH_input_run_options_RDRS.ini		        		    ./MESH_input_run_options.ini
ln -sf $config/distributed_param/MESH_parameters_Thresh100.r2c		    ./MESH_parameters.r2c
ln -sf $config/class_hydrology/MESH_parameters_CLASS_Liard_edit4.ini	./MESH_parameters_CLASS.ini
ln -sf $config/class_hydrology/MESH_parameters_hydrology_Liard.ini	    ./MESH_parameters_hydrology.ini
  
mpirun -np 1 $dir_mesh/mpi_sa_mesh_r1813      