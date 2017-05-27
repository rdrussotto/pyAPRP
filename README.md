# pyAPRP

Python implementation of the Approximate Partial Radiation Perturbation method (Taylor et al., J. Climate, 2007), 
abbreviated as APRP. This method calculates changes in net radiative fluxes at the top of atmosphere in climate models
due to changes in cloud properties, surface albedo, and non-cloud atmospheric absorption and scattering.
It is useful for diagnosing SW cloud feedbacks when only monthly mean model output is available. 

Python code written by Rick Russotto, based on a Matlab implementation of APRP by Yen-Ting Hwang. 
Some sections of code were translated to Python with few changes; others were somewhat more modified.

# Dependencies

Code originally written in Python 2.7, importing 
numpy version 1.8.2 and NetCDF4 version 1.1.0.

# How to run

Import the aprp.py file into another script, and run the aprp_main function. APRP is based on 
comparing two different time periods; the two time periods can be from the same model run or different model runs.
The function takes 6 arguments: 

dataPaths1:  dictionary of paths to the netCDF output for time period 1

firstMonth1: first month (indexed from beginning of output) for time period 1--note Python indices start with 0

lastMonth1:   last month (indexed from beginning of output) for time period 1

dataPaths2:  dictionary of paths to the netCDF output for time period 2 
             (if the two states being compared are different times from the same run, make this the same as dataPaths1)
             
firstMonth2: first month (indexed from beginning of output) for time period 2

lastMonth2:   last month (indexed from beginning of output) for time period 2

Keys in the dataPaths dictionaries should be variable names following the 
standards used in the Coupled Model Intercomparison Project, phase 5 (CMIP5). 
Values should be strings containing paths to the netCDF files. 
The following model output variables must be included:

clt: total cloud fraction

rsds: downwelling SW radiation at surface

rsdt: downwelling SW radiation at top of atmosphere

rsus: upwelling SW radiation at surface

rsut: upwelling SW radiation at top of atmosphere

rsdscs: clear-sky downwelling SW radiation at surface

rsuscs: clear-sky upwelling SW radiation at surface

rsutcs: clear-sky upwelling SW radiation at top of atmosphere

There are alternative versions of the main and loadNetCDF functions at the bottom of the file 
for native Community Earth System Model (CESM) model output. 

# References

Original reference for the APRP method: 

Taylor, K. E., Crucifix, M., Braconnot, P., Hewitt, C. D., Doutriaux, C., 
Broccoli, A. J., Mitchell, J. F. B., and Webb, M. J.: 
Estimating Shortwave Radiative Forcing and Response in Climate Models, Journal of Climate, 20, 2530–2543, 
doi:10.1175/JCLI4143.1, 2007.

Yen-Ting Hwang's Matlab implementation was first used in: 

Hwang, Y.-T. and Frierson, D. M. W.: Increasing atmospheric poleward energy transport with global warming, 
Geophysical Research Letters, 37, L24807, doi:10.1029/2010GL045440, 2010.

The Python version will be used in an upcoming paper by Rick Russotto and Thomas Ackerman.
