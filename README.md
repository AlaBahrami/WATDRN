# WATDRN
## Purpose
The purpose of this repository is to provide a stand-alone version of the WATDRN [Soulis et al., 2000]( https://www.tandfonline.com/doi/abs/10.1080/07055900.2000.9649648) algorithm which has been integrated into the MESH/CLASSIC models. 

___
# Folder Structure
Given the file size limitations of GitHub, only smaller files are stored here and the rest are stored on Graham. The files can be synced with the local machine via the respective push/pull bash scripts included in the Data/Raw and Model folders.


## MESH

### MESH_code 
- Includes the modified version of MESH-code to write NetCDF files that you can test WATDRN as a stand-lone for any domain of interest. Here the CLASSW subroutine has been modified to write time-variant and time-invariant outputs, so users can read it later in the WATDRN_SA program.

### Fraser_setup
- Includes MESH Fraser setup and configuration files which are used to test whether the stand-alone version of WATDRN works properly. 

### WATDRN_Standalone
- Includes a stand-alone version of the WATDRN code which is used to find out how the WATDRN and WATROF subroutines are working when input variables and parameters are obtained from the MESH Fraser setup. Both Python and Fortran routines are provided here (WATDRN_SA.f90, WATDRN_SA.py). 
- Includes the WATDRN_SA_simulation.f90 which is a translation of the Python version (WATDRN_SA_simulation.py). The purpose is to test whether any changes between Fortran and Python synthetic simulations will be observed.

## WATDRN_synthetic_tests
- Includes a stand-alone version of WATDRN (watdrn_simulation.py) which is used for running a series of synthetic tests and comparing the performance of WATDRN against the method of characteristics and the TOPMODEL. The entire idea and mathematical background for the synthetic tests have been provided by Prof. [Martyn Clark](https://github.com/martynpclark). A detailed explanation has been provided in the WATDRN critique draft. For more information check laugh tests implemented for the SUMMA model based on test cases defined in Wigmosta and Lettenmaier (1999) [SUMMA LaughTest](https://github.com/CH-Earth/laughTests/tree/master/lt4_wigmosta1994)  
   

