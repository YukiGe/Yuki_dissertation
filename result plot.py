# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 19:58:02 2023

@author: admin
"""

##Import standard python modules:
import numpy as np #For array data types and basic math calculations
import matplotlib.pyplot as plt #For making plots
import pandas as pd #For handling 'data frames' (essentially spreadsheets), loading spreadsheets
import geopandas as gpd #For storing GIS data as a data frame. Essentially the same as a pandas data frame but has a column descriping the 'shape'
import datetime  #For handling dates and times e.g. can convert seconds to a dates. 
from calendar import monthrange #to work out how many days are in a month
import xarray as xr #for handling .nc (NetCDF) files which the GIS temperture data is stored as.
from scipy.integrate import cumtrapz #For numerical integration
import os #to do operating system commands

##Import modules for calculating slurry emissions (These are other python files in this folder - look in the .py files to see how they work!)
from BasicParams import * #BasicParams has standard parameters used for working out how much slurry is being generated.  This way of importing means you can access the paramters by using >>MCF_k rather than >>BasicParams.MCF_k
import LS_Energy #Uses IPCC calculations to work out how much energy and food a cow needs given its milk yeild, growth, temperture and how much walking it does.
import LS_Emissions #Uses the output of LS_Energy functions to work out enteric methane production and how much slurry is being generated.
import LS_Slurry_Methane_Temp as MethaneTemp #Uses how much slurry is being generated and temperture to work out methane production from slurry with a model.
import LS_SlurrySpreadTime as SpreadTime #Has a function for working out "voltaile solids" (i.e. slurry that can turn to methane) generation rate Can work out when is the best time to spread slurry on field or use AD to minimize methane emissions for that year.
import os #to do operating system commands


#download dataframe(csv)
solved_gdf_2019=pd.read_csv('../LiveStockCalculations-master/result/result_2019.csv')
FIG_SAVE_DIR='C:/Users/admin/Desktop/LiveStockCalculations-master/plot/'
column_converter={'projection_y_coordinate':'projection','projection_x_coordinate':'projecti_1','month_number':'month_numb'}   #transfer columns's name since saving the shp file in short name                                
#Now we can plot the TM kgCH4, for example, in the UK for each month:
for m in np.arange(1,13): #loop through all months
    month_solved_gdf_2019=solved_gdf_2019[solved_gdf_2019[column_converter['month_number']]==m] #filter to the month of intrest
    month_solved_gdf_2019.plot(column='',cmap='rainbow',markersize=30,legend=True,figsize=(5, 9)) #vmin/vmax sets the limits on the colorbar
    mean_TM=np.mean(month_solved_gdf_2019['TM kgCH4']) #Calculate the mean TM kgCH4 for this month over the UK
    plt.axis('off')
    plt.title('MeanTM kgCH4 {M:1.2f}'.format(M=mean_TM))#Make the title of the plot the mean TM kgCH4 over the UK
    plt.savefig(FIG_SAVE_DIR+'TM kgCH4_2019_'+str(m)+'.png',dpi=300) #Save the figure
    plt.close()
