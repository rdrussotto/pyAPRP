# APRP.py: module containing functions, etc. related to the Approximate Partial Radiation Perturbation method 
# (Taylor et al., 2007). Based on Matlab code written by Yen-Ting Hwang.
#
# To run, import this module into an outside script and then run "aprp_main" for 
# the desired model 

import numpy as np
import netCDF4 as nc4



#Main function to run, for an individual model (which model is specified via the "dataPaths" arguments).
#This version assumes variable names and dimensions following CMIP convention.
#
#Inputs:
#   dataPaths1:  dictionary of paths to the netCDF output for time period 1 
#   firstMonth1: first month (indexed from beginning of output) for time period 1--note Python indices start with 0
#   lastMonth1:   last month (indexed from beginning of output) for time period 1
#   dataPaths2:  dictionary of paths to the netCDF output for time period 2 
#                (if the two states being compared are different times from the same run, make this the same as dataPaths1)
#   firstMonth2: first month (indexed from beginning of output) for time period 2
#   lastMonth2:   last month (indexed from beginning of output) for time period 2
#
#Outputs: 
#   A dictionary of dictionaries: 
#   returnDict['APRP']: contains results from comparing the two time periods (see "d_albedo" function for list of variables)
#   returnDict['Time1_preliminaries']: contains the relevant model output variables, having been read in from NetCDF
#                                      files, and processed including cloudy-sky calculations and multiannual means
#   returnDict['Time2_preliminaries']: same as above but for time period 2 (see "loadNetCDF" function for list of variables)
#   returnDict['Time1_parameters']: contains tuning parameters for the idealized single-layer radiative transfer model
#   returnDict['Time2_parameters']: same as above but for time period 2 (see "parameters" function for list of variables)
#   Syntax for accessing: e.g. to get the radiative effect of surface albedo changes, type returnDict['APRP']['surface']
#
def aprp_main(dataPaths1, firstMonth1, lastMonth1, dataPaths2, firstMonth2, lastMonth2):
    #Load files and run calculations for first time period
    dict1A = loadNetCDF(dataPaths1, firstMonth1, lastMonth1)
    dict1B = parameters(dict1A)
    #Load files and run calculations for second time period
    dict2A = loadNetCDF(dataPaths2, firstMonth2, lastMonth2)
    dict2B = parameters(dict2A)
    #Run calculations regarding change betweeen 2 time periods
    dictC = d_albedo(dict1A, dict1B, dict2A, dict2B)
    #Nest the dictionaries into an outside dictionary to return 
    returnDict = dict()
    returnDict['APRP'] = dictC
    returnDict['Time1_preliminaries'] = dict1A
    returnDict['Time1_parameters'] = dict1B
    returnDict['Time2_preliminaries'] = dict2A
    returnDict['Time2_parameters'] = dict2B
    return returnDict
    
        
    
    
    

#Load variables from netCDF files (run twice, once for each time) and calculate overcast sky data. 
#Based on "load_nc_coupled.m" in Ting's code.
#Inputs: see aprp_main
#Outputs: a dictionary containing monthly mean SW fluxes at surface and TOA under all-sky, clear-sky and 
#         overcast conditions, as well as the model's latitude and longitude grids. These are monthly mean 
#         data; unlike in Ting's code, not doing multi-annual mean yet. Better to do APRP calculations on 
#         the individual months. 
def loadNetCDF(dataPaths, firstMonth, lastMonth): 
    #Variable names from CMIP convention (dictionary of data paths should have labels corresponding to these)
    #variables = ['rsds', 'rsus', 'rsut', 'rsdt', 'rsutcs', 'rsdscs', 'rsuscs', 'clt']
    
    #For each of the variables, import the netCDF file and extract array from the netCDF Dataset object, subsetted
    #...by the times specified in the arguments (Ting used a loop with "eval" but I'd rather avoid that 
    #...for readability) and mask values greater than 10^10
    Dataset = nc4.Dataset(dataPaths['rsds'])
    rsds = Dataset.variables['rsds'][firstMonth:lastMonth+1, :,:]
    rsds = np.ma.masked_greater(rsds,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['rsus'])
    rsus = Dataset.variables['rsus'][firstMonth:lastMonth+1, :,:]
    rsus = np.ma.masked_greater(rsus,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['rsut'])
    rsut = Dataset.variables['rsut'][firstMonth:lastMonth+1, :,:]
    rsut = np.ma.masked_greater(rsut,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['rsdt'])
    rsdt = Dataset.variables['rsdt'][firstMonth:lastMonth+1, :,:]
    rsdt = np.ma.masked_greater(rsdt,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['rsutcs'])
    rsutcs = Dataset.variables['rsutcs'][firstMonth:lastMonth+1, :,:]
    rsutcs = np.ma.masked_greater(rsutcs,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['rsdscs'])
    rsdscs = Dataset.variables['rsdscs'][firstMonth:lastMonth+1, :,:]
    rsdscs = np.ma.masked_greater(rsdscs,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['rsuscs'])
    rsuscs = Dataset.variables['rsuscs'][firstMonth:lastMonth+1, :,:]
    rsuscs = np.ma.masked_greater(rsuscs,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['clt'])
    clt = Dataset.variables['clt'][firstMonth:lastMonth+1, :,:]
    clt = np.ma.masked_greater(clt,1.e10)
        
    #Alternative to the repetitive code above: shorter but harder to read/debug        
    #for variable in variables: 
    #    Dataset = nc4.Dataset(dataPaths[variable])
    #    eval(variable+'= Dataset.variables['+variable+'][firstMonth:lastMonth+1, :,:]')
    #    eval(variable+'= np.ma.masked_greater('+variable+', 1.e10)') 
       
    #Obtain the latitude and longitude for the model (using last Dataset in the loop which should still be available)
    lat = Dataset.variables['lat'][:]
    lon = Dataset.variables['lon'][:]
    
    #Here Ting calculated multi-year means for individual months. I need to do this too.
    #Dimensions are time, lat, lon
    #Need to average over every 12th time element, leave the lat and lon dependence. 
    #Will end up with a 3D array whose dimensions are month (1-12), lat, lon. 
    #What is best way to do this? 
    #Ting looped over the 12 months.
    #She also saved separate 1-month means, but never used them so I'll skip that for now. 
    numMonths = lastMonth - firstMonth + 1
    m_rsds = np.zeros([12,len(lat),len(lon)])
    m_rsus = np.zeros([12,len(lat),len(lon)])
    m_rsut = np.zeros([12,len(lat),len(lon)])
    m_rsdt = np.zeros([12,len(lat),len(lon)])
    m_rsutcs = np.zeros([12,len(lat),len(lon)])
    m_rsdscs = np.zeros([12,len(lat),len(lon)])
    m_rsuscs = np.zeros([12,len(lat),len(lon)])
    m_clt = np.zeros([12,len(lat),len(lon)])
    for i in range(0,12):
        m_rsds[i,:,:] = np.mean(rsds[i:numMonths:12,:,:], axis=0)
        m_rsus[i,:,:] = np.mean(rsus[i:numMonths:12,:,:], axis=0)
        m_rsut[i,:,:] = np.mean(rsut[i:numMonths:12,:,:], axis=0)
        m_rsdt[i,:,:] = np.mean(rsdt[i:numMonths:12,:,:], axis=0)
        m_rsutcs[i,:,:] = np.mean(rsutcs[i:numMonths:12,:,:], axis=0)
        m_rsdscs[i,:,:] = np.mean(rsdscs[i:numMonths:12,:,:], axis=0)
        m_rsuscs[i,:,:] = np.mean(rsuscs[i:numMonths:12,:,:], axis=0)
        m_clt[i,:,:] = np.mean(clt[i:numMonths:12,:,:], axis=0)
    
        
    #Calculate the overcast versions of rsds, rsus, rsut from the clear-sky and all-sky data
    #First mask zero values of cloud fraction so you don't calculate overcast values in clear-sky pixels
    m_clt = np.ma.masked_values(m_clt, 0)
    c = m_clt/100. #c is cloud fraction. clt was in percentages
    m_rsdsoc = (m_rsds-(1.-c)*(m_rsdscs))/c  #Can derive this algebraically from Taylor et al., 2007, Eq. 3
    m_rsusoc = (m_rsus-(1.-c)*(m_rsuscs))/c
    m_rsutoc = (m_rsut-(1.-c)*(m_rsutcs))/c
    
    #Mask zero values of the downward SW radiation (I assume this means polar night, for monthly mean)
    m_rsds = np.ma.masked_values(m_rsds, 0)
    m_rsdscs = np.ma.masked_values(m_rsdscs, 0)
    m_rsdsoc = np.ma.masked_values(m_rsdsoc, 0)
    m_rsdt = np.ma.masked_values(m_rsdt, 0)
    
    #Return dictionary with all the variables calculated here (called "dictA" because calculated in first function called)
    dictA = dict()
    dictA['rsds'] = m_rsds
    dictA['rsus'] = m_rsus
    dictA['rsut'] = m_rsut
    dictA['rsdt'] = m_rsdt
    dictA['rsutcs'] = m_rsutcs
    dictA['rsdscs'] = m_rsdscs
    dictA['rsuscs'] = m_rsuscs
    dictA['clt'] = m_clt
    dictA['lat'] = lat
    dictA['lon'] = lon
    dictA['rsdsoc'] = m_rsdsoc
    dictA['rsusoc'] = m_rsusoc
    dictA['rsutoc'] = m_rsutoc
    dictA['c'] = c #Cloud fraction as fraction, not %
    
    return dictA
    
    
    
    
    
    
    
#Calculate the tuning parameters for the idealized single-layer radiative transfer model
#for the individual time period (i.e. control or warmed)
#See Figure 1 of Taylor et al., 2007, and other parts of that paper. Equations referenced are from there.
#
#Based on Ting's "parameters.m".
#
#Inputs: the dictionary output by loadNetCDF
#Outputs: a dictionary of additional outputs
def parameters(dictA):
    #Clear-sky parameters
    a_clr = dictA['rsuscs']/dictA['rsdscs'] #Surface albedo   
    Q = dictA['rsdscs']/dictA['rsdt'] #Ratio of incident surface flux to insolation
    mu_clr = dictA['rsutcs']/dictA['rsdt']+Q*(1.-a_clr) #Atmospheric transmittance (Eq. 9)  #"Invalid value in divide"
    ga_clr = (mu_clr-Q)/(mu_clr-a_clr*Q) #Atmospheric scattering coefficient (Eq. 10)

    #Overcast parameters
    a_oc = dictA['rsusoc']/dictA['rsdsoc'] #Surface albedo
    Q = dictA['rsdsoc']/dictA['rsdt'] #Ratio of incident surface flux to insolation
    mu_oc = dictA['rsutoc']/dictA['rsdt']+Q*(1.-a_oc) #Atmospheric transmittance (Eq. 9)
    ga_oc = (mu_oc-Q)/(mu_oc-a_oc*Q) #Atmospheric scattering coefficient (Eq. 10)   
    
    #Calculating cloudy parameters based on clear-sky and overcast ones 
    #Difference between _cld and _oc: _cld is due to the cloud itself, as opposed to 
    #scattering and absorption from all constituents including clouds in overcast skies.
    mu_cld = mu_oc / mu_clr            #Eq. 14
    ga_cld = (ga_oc-1.)/(1.-ga_clr)+1. #Eq. 13  
    
    #Save the relevant variables to a dictionary for later use
    dictB = dict()
    dictB['a_clr'] = a_clr
    dictB['a_oc'] = a_oc
    dictB['mu_clr'] = mu_clr
    dictB['mu_cld'] = mu_cld
    dictB['ga_clr'] = ga_clr
    dictB['ga_cld'] = ga_cld
    
    #Ting saved a cloud fraction variable here--I did this in earlier function instead.
    
    return dictB







#Calculations for the differences between time periods
def d_albedo(dict1A, dict1B, dict2A, dict2B):
    
    #First, Ting set cloud values that were masked in one time period 
    #equal to the value in the other time period, assuming no cloud changes.
    #I'll take these variables out of the dictionary before modifying them. 
    a_oc1 = dict1B['a_oc']
    a_oc2 = dict2B['a_oc']
    a_oc2[a_oc2.mask == True] = a_oc1[a_oc2.mask == True]
    a_oc1[a_oc1.mask == True] = a_oc2[a_oc1.mask == True]
    
    mu_cld1 = dict1B['mu_cld']
    mu_cld2 = dict2B['mu_cld']
    mu_cld2[mu_cld2.mask == True] = mu_cld1[mu_cld2.mask == True]
    mu_cld1[mu_cld1.mask == True] = mu_cld2[mu_cld1.mask == True]
    
    ga_cld1 = dict1B['ga_cld']
    ga_cld2 = dict2B['ga_cld']
    ga_cld2[ga_cld2.mask == True] = ga_cld1[ga_cld2.mask == True]
    ga_cld1[ga_cld1.mask == True] = ga_cld2[ga_cld1.mask == True]
    
    #Now a bunch of calls to the "albedo" function to see how the albedo changes as a result of 
    #...the changes to each of the radiative components. 
    
    #Retrieve other variables from dictionaries to make calls to albedo shorter/more readable
    c1 = dict1A['c']
    c2 = dict2A['c']
    a_clr1 = dict1B['a_clr']
    a_clr2 = dict2B['a_clr']
    mu_clr1 = dict1B['mu_clr']
    mu_clr2 = dict2B['mu_clr']
    ga_clr1 = dict1B['ga_clr']
    ga_clr2 = dict2B['ga_clr']
    
    #Base state albedo
    A1 = albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld1)
    A2 = albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld2)
    
    #Change in albedo due to each component (Taylor et al., 2007, Eq. 12b)
    dA_c =      .5*(albedo(c2, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c1, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld2))    
    
    dA_a_clr =  .5*(albedo(c1, a_clr2, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr1, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld2))
                
    dA_a_oc =   .5*(albedo(c1, a_clr1, a_oc2, mu_clr1, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc1, mu_clr2, mu_cld2, ga_clr2, ga_cld2))
               
    dA_mu_clr = .5*(albedo(c1, a_clr1, a_oc1, mu_clr2, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr1, mu_cld2, ga_clr2, ga_cld2))               
                 
    dA_mu_cld = .5*(albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld2, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld1, ga_clr2, ga_cld2))                 
                 
    dA_ga_clr = .5*(albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr2, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr1, ga_cld2))
                 
    dA_ga_cld = .5*(albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld2)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld1))
                 
    #Set changes due to overcast or cloudy sky parameters, or changes to clouds themselves, to zero
    #...if cloud fraction is less than 3% in either time period
    dA_a_oc[dict1A['c'] < .03] = 0
    dA_a_oc[dict2A['c'] < .03] = 0
    dA_mu_cld[dict1A['c'] < .03] = 0
    dA_mu_cld[dict2A['c'] < .03] = 0
    dA_ga_cld[dict1A['c'] < .03] = 0
    dA_ga_cld[dict2A['c'] < .03] = 0
    dA_c[dict1A['c'] < .03] = 0
    dA_c[dict2A['c'] < .03] = 0
    
    #Combine different components into changes due to surface albedo, atmospheric clear-sky and atmospheric cloudy-sky
    dA_a = dA_a_clr + dA_a_oc                #Eq. 16a
    dA_cld = dA_mu_cld + dA_ga_cld + dA_c    #Eq. 16b
    dA_clr = dA_mu_clr + dA_ga_clr           #Eq. 16c
    
    #Set all planetary albedo changes = zero when incoming solar radaition is zero    
    #(This will replace NaNs with zeros in the polar night--affects annual means)
    dA_a[dict2A['rsdt']<0.1] = 0
    dA_clr[dict2A['rsdt']<0.1] = 0
    dA_cld[dict2A['rsdt']<0.1] = 0
    dA_a_clr[dict2A['rsdt']<0.1] = 0
    dA_a_oc[dict2A['rsdt']<0.1] = 0
    dA_mu_cld[dict2A['rsdt']<0.1] = 0
    dA_ga_cld[dict2A['rsdt']<0.1] = 0
    dA_c[dict2A['rsdt']<0.1] = 0
    dA_mu_clr[dict2A['rsdt']<0.1] = 0
    dA_ga_clr[dict2A['rsdt']<0.1] = 0
    
    
    #Calculate radiative effects in W/m^2 by multiplying negative of planetary albedo changes by downward SW radation
    #(This means positive changes mean more downward SW absorbed)
    surface = -dA_a*dict2A['rsdt']   #Radiative effect of surface albedo changes
    surface[dict2A['rsdt']<0.1] = 0
    surface = np.ma.masked_outside(surface, -100, 100) # Ting called this "boundary for strange output"
    
    cloud = -dA_cld*dict2A['rsdt']   #Radiative effect of cloud changes
    cloud[dict2A['rsdt']<0.1] = 0
    cloud = np.ma.masked_outside(cloud, -100, 100) # Ting called this "boundary for strange output"
    
    noncloud = -dA_clr*dict2A['rsdt'] #Radiative effect of non-cloud SW changes (e.g. SW absorption)
    noncloud[dict2A['rsdt']<0.1] = 0
    
    #Broken down further into the individual terms in Eq. 16
    surface_clr = -dA_a_clr*dict2A['rsdt']    #Effects of surface albedo in clear-sky conditions
    surface_clr[dict2A['rsdt']<0.1] = 0
    
    surface_oc = -dA_a_oc*dict2A['rsdt']      #Effects of surface albedo in overcast conditions
    surface_oc[dict2A['rsdt']<0.1] = 0
    
    cloud_c = -dA_c*dict2A['rsdt']            #Effects of changes in cloud fraction
    cloud_c[dict2A['rsdt']<0.1] = 0
    
    cloud_ga = -dA_ga_cld*dict2A['rsdt']      #Effects of atmospheric scattering in cloudy conditions
    cloud_ga[dict2A['rsdt']<0.1] = 0
    
    cloud_mu = -dA_mu_cld*dict2A['rsdt']      #Effects of atmospheric absorption in cloudy conditions
    cloud_mu[dict2A['rsdt']<0.1] = 0
    
    noncloud_ga = -dA_ga_clr*dict2A['rsdt']   #Effects of atmospheric scattering in clear-sky conditions
    noncloud_ga[dict2A['rsdt']<0.1] = 0
    
    noncloud_mu = -dA_mu_clr*dict2A['rsdt']   #Effects of atmospheric absorption in clear-sky conditions
    noncloud_mu[dict2A['rsdt']<0.1] = 0
    
    #Calculate more useful radiation output
    CRF = dict1A['rsut'] - dict1A['rsutcs'] - dict2A['rsut'] + dict2A['rsutcs'] #Change in cloud radiative effect
    cs = dict1A['rsutcs'] - dict2A['rsutcs']  #Change in clear-sky upward SW flux at TOA
    
    #Define a dictionary to return all the variables calculated here
    dictC = dict()
    dictC['A1'] = A1
    dictC['A2'] = A2
    dictC['dA_c'] = dA_c
    dictC['dA_a_clr'] = dA_a_clr
    dictC['dA_a_oc'] = dA_a_oc
    dictC['dA_mu_clr'] = dA_mu_clr
    dictC['dA_mu_cld'] = dA_mu_cld
    dictC['dA_ga_clr'] = dA_ga_clr
    dictC['dA_ga_cld'] = dA_ga_cld
    dictC['dA_a'] = dA_a
    dictC['dA_cld'] = dA_cld
    dictC['dA_clr'] = dA_clr
    dictC['surface'] = surface
    dictC['cloud'] = cloud
    dictC['noncloud'] = noncloud
    dictC['surface_clr'] = surface_clr
    dictC['surface_oc'] = surface_oc
    dictC['cloud_c'] = cloud_c
    dictC['cloud_ga'] = cloud_ga
    dictC['cloud_mu'] = cloud_mu
    dictC['noncloud_ga'] = noncloud_ga
    dictC['noncloud_mu'] = noncloud_mu
    dictC['CRF'] = CRF
    dictC['cs'] = cs
    
    return dictC
    
    
    
    
    
#Function to calculate the planetary albedo, A.
#Inputs: (see Fig. 1 of Taylor et al., 2007)
#   c: fraction of the region occupied by clouds
#   a_clr: clear sky surface albedo (SW flux up / SW flux down)
#   a_oc: overcast surface albedo
#   mu_clr: clear-sky transmittance of SW radiation
#   mu_cld: cloudy-sky transmittance of SW radiation
#   ga_clr: clear-sky atmospheric scattering coefficient
#   ga_cld: cloudy-sky atmospheric scattering coefficient
def albedo(c, a_clr, a_oc, mu_clr, mu_cld, ga_clr, ga_cld): #Labeled with equation numbers from Taylor et al. 2007
    mu_oc = mu_clr*mu_cld                                                            #Eq. 14
    ga_oc =  1. - (1.-ga_clr)*(1.-ga_cld)                                            #Eq. 13
    A_clr = mu_clr*ga_clr + mu_clr*a_clr*(1.-ga_clr)*(1.-ga_clr)/(1.-a_clr*ga_clr)   #Eq. 7 (clear-sky)
    A_oc = mu_oc*ga_oc + mu_oc*a_oc*(1.-ga_oc)*(1.-ga_oc)/(1.-a_oc*ga_oc)            #Eq. 7 (overcast sky)
    A = (1-c)*A_clr + c*A_oc                                                         #Eq. 15
    return A
    
    
    
    
    
    
    
#### Alternative versions for CESM model runs with different output variable names ####




#Alternative main function to run the different loading function
def aprp_main_cesm(dataPaths1, firstMonth1, lastMonth1, dataPaths2, firstMonth2, lastMonth2):
    #Load files and run calculations for first time period
    dict1A = loadNetCDF_cesm(dataPaths1, firstMonth1, lastMonth1)
    dict1B = parameters(dict1A)
    #Load files and run calculations for second time period
    dict2A = loadNetCDF_cesm(dataPaths2, firstMonth2, lastMonth2)
    dict2B = parameters(dict2A)
    #Run calculations regarding change betweeen 2 time periods
    dictC = d_albedo(dict1A, dict1B, dict2A, dict2B)
    #Nest the dictionaries into an outside dictionary to return 
    returnDict = dict()
    returnDict['APRP'] = dictC
    returnDict['Time1_preliminaries'] = dict1A
    returnDict['Time1_parameters'] = dict1B
    returnDict['Time2_preliminaries'] = dict2A
    returnDict['Time2_parameters'] = dict2B
    return returnDict

#Loading function for CESM output variable names (output the same as the loadNetCDF() function)
def loadNetCDF_cesm(dataPaths, firstMonth, lastMonth): 
    #Variable names from CAM output (dictionary of data paths should have labels corresponding to these)
    #vars_CAM = ['FSDS', 'FSNS', 'FSUTOA', 'FSNTOA', 'FSNTOAC', 'FSDSC', 'FSNSC', 'CLDTOT']
        
    #For each variable, import the netCDF file and extract array from the netCDF Dataset object, subsetted
    #...by the times specified in the arguments 
    Dataset = nc4.Dataset(dataPaths['FSDS'])
    FSDS = Dataset.variables['FSDS'][firstMonth:lastMonth+1, :,:]
    FSDS = np.ma.masked_greater(FSDS,1.e10)    
    
    Dataset = nc4.Dataset(dataPaths['FSNS'])
    FSNS = Dataset.variables['FSNS'][firstMonth:lastMonth+1, :,:]
    FSNS = np.ma.masked_greater(FSNS,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['FSUTOA'])
    FSUTOA = Dataset.variables['FSUTOA'][firstMonth:lastMonth+1, :,:]
    FSUTOA = np.ma.masked_greater(FSUTOA,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['FSNTOA'])
    FSNTOA = Dataset.variables['FSNTOA'][firstMonth:lastMonth+1, :,:]
    FSNTOA = np.ma.masked_greater(FSNTOA,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['FSNTOAC'])
    FSNTOAC = Dataset.variables['FSNTOAC'][firstMonth:lastMonth+1, :,:]
    FSNTOAC = np.ma.masked_greater(FSNTOAC,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['FSDSC'])
    FSDSC = Dataset.variables['FSDSC'][firstMonth:lastMonth+1, :,:]
    FSDSC = np.ma.masked_greater(FSDSC,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['FSNSC'])
    FSNSC = Dataset.variables['FSNSC'][firstMonth:lastMonth+1, :,:]
    FSNSC = np.ma.masked_greater(FSNSC,1.e10)
    
    Dataset = nc4.Dataset(dataPaths['CLDTOT'])
    CLDTOT = Dataset.variables['CLDTOT'][firstMonth:lastMonth+1, :,:]
    CLDTOT = np.ma.masked_greater(CLDTOT,1.e10)
        

    #Variable names from CMIP convention (used in rest of the program)
    #variables = ['rsds', 'rsus', 'rsut', 'rsdt', 'rsutcs', 'rsdscs', 'rsuscs', 'clt']    
        
    #Get the variables into the CMIP format/name convention (some are already fine; some need processing)
    rsds = FSDS
    rsut = FSUTOA
    rsdscs = FSDSC
    clt = CLDTOT
    rsus = FSDS - FSNS   #SW: positive down: net = down - up ---> up = down - net
    rsdt = FSUTOA + FSNTOA    #down = net + up
    rsuscs = FSDSC - FSNSC
    rsutcs = rsdt - FSNTOAC #Downward SW at TOA should be same regardless of clouds.
        
    ####### from here down, same as regular loadNetCDF #######


    #Obtain the latitude and longitude for the model (using last Dataset in the loop which should still be available)
    lat = Dataset.variables['lat'][:]
    lon = Dataset.variables['lon'][:]
    
    #Here Ting calculated multi-year means for individual months. I need to do this too.
    #Dimensions are time, lat, lon
    #Need to average over every 12th time element, leave the lat and lon dependence. 
    #Will end up with a 3D array whose dimensions are month (1-12), lat, lon. 
    #What is best way to do this? 
    #Ting looped over the 12 months.
    #She also saved separate 1-month means, but never used them so I'll skip that for now. 
    numMonths = lastMonth - firstMonth + 1
    m_rsds = np.zeros([12,len(lat),len(lon)])
    m_rsus = np.zeros([12,len(lat),len(lon)])
    m_rsut = np.zeros([12,len(lat),len(lon)])
    m_rsdt = np.zeros([12,len(lat),len(lon)])
    m_rsutcs = np.zeros([12,len(lat),len(lon)])
    m_rsdscs = np.zeros([12,len(lat),len(lon)])
    m_rsuscs = np.zeros([12,len(lat),len(lon)])
    m_clt = np.zeros([12,len(lat),len(lon)])
    for i in range(0,12):
        m_rsds[i,:,:] = np.mean(rsds[i:numMonths:12,:,:], axis=0)
        m_rsus[i,:,:] = np.mean(rsus[i:numMonths:12,:,:], axis=0)
        m_rsut[i,:,:] = np.mean(rsut[i:numMonths:12,:,:], axis=0)
        m_rsdt[i,:,:] = np.mean(rsdt[i:numMonths:12,:,:], axis=0)
        m_rsutcs[i,:,:] = np.mean(rsutcs[i:numMonths:12,:,:], axis=0)
        m_rsdscs[i,:,:] = np.mean(rsdscs[i:numMonths:12,:,:], axis=0)
        m_rsuscs[i,:,:] = np.mean(rsuscs[i:numMonths:12,:,:], axis=0)
        m_clt[i,:,:] = np.mean(clt[i:numMonths:12,:,:], axis=0)
    
        
    #Calculate the overcast versions of rsds, rsus, rsut from the clear-sky and all-sky data
    #First mask zero values of cloud fraction so you don't calculate overcast values in clear-sky pixels
    m_clt = np.ma.masked_values(m_clt, 0)
    # c = m_clt/100. #c is cloud fraction. clt was in percentages  #No-Not true in CESM output
    c = m_clt
    m_rsdsoc = (m_rsds-(1.-c)*(m_rsdscs))/c  #Can derive this algebraically from Taylor et al., 2007, Eq. 3
    m_rsusoc = (m_rsus-(1.-c)*(m_rsuscs))/c
    m_rsutoc = (m_rsut-(1.-c)*(m_rsutcs))/c
    
    #Mask zero values of the downward SW radiation (I assume this means polar night, for monthly mean)
    m_rsds = np.ma.masked_values(m_rsds, 0)
    m_rsdscs = np.ma.masked_values(m_rsdscs, 0)
    m_rsdsoc = np.ma.masked_values(m_rsdsoc, 0)
    m_rsdt = np.ma.masked_values(m_rsdt, 0)
    
    #Return dictionary with all the variables calculated here (called "dictA" because calculated in first function called)
    dictA = dict()
    dictA['rsds'] = m_rsds
    dictA['rsus'] = m_rsus
    dictA['rsut'] = m_rsut
    dictA['rsdt'] = m_rsdt
    dictA['rsutcs'] = m_rsutcs
    dictA['rsdscs'] = m_rsdscs
    dictA['rsuscs'] = m_rsuscs
    dictA['clt'] = m_clt
    dictA['lat'] = lat
    dictA['lon'] = lon
    dictA['rsdsoc'] = m_rsdsoc
    dictA['rsusoc'] = m_rsusoc
    dictA['rsutoc'] = m_rsutoc
    dictA['c'] = c #Cloud fraction as fraction, not %
    
    return dictA
    