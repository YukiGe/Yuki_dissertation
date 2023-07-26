# -*- coding: utf-8 -*-
"""
Created on Tue May 16 10:47:15 2023

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
from shapely.geometry import Polygon
import shapely 
def load_temp_nc_data(path,epsg): 
    '''
    Will load the .nc files sourced from Met Office; Hollis, D.; McCarthy, M.; Kendon, M.; Legg, T.; Simpson, I. (2020): HadUK-Grid Gridded Climate Observations on a 25km grid over the UK, v1.0.2.1 (1862-2019).
    -path (string) is the path to the .nc file
    -epsg (int) is a standardised code that defines the 'coordinate reference system' (CRS) of the globe to the plane. Choose the code to define the crs of the output dataframe. Use epsg=3395 for one that has units meters.
    '''
    ar = xr.open_dataset(path) #open netCDF (.nc) file using xarray
    df = ar.to_dataframe() #make a  dataframe
    df=df.reset_index() #reset indexes in the df so it starts from 1
    gdf = gpd.GeoDataFrame(df,geometry=gpd.points_from_xy(df.projection_x_coordinate,df.projection_y_coordinate)) #Make a geopandas dataframe

    gdf=gdf.set_crs(epsg=27700)#this is what  https://doi.org/10.1002/gdj3.78 says the crs is of the data is.
    
    gdf=gdf[gdf['bnds']==0] #Gets rid of the bnds=1 rows. However, these might be good for sorting out the box edges
    
    gdf=gdf[np.isnan(gdf['tas'])==0] ##To remove rows which have a nan in the 'tas' column
    gdf=gdf.reset_index() #reindex the data frame since we just removed a lot of rows.
    gdf=gdf.to_crs(epsg=epsg) # convert to the epsg vaiable
 
    return gdf #the function returns the geopandas data frame
#a CRS with units meters
epsg_m="27700"
cattle_data='../data/scotland-agcensus-2019_5km-grid_5068939/cattle/sc-2019-c103-5km-total-cattle.shp'
temp_data_folder='../data/HadUk-Grid25km-1961-2020_Monthly/'
gdf_temp = load_temp_nc_data (temp_data_folder+'tas_hadukgrid_uk_25km_mon_202001-202012.nc',epsg_m)
gdf=gpd.read_file(cattle_data)
gdf_temp['geometry']=gdf_temp.geometry.buffer(25000,cap_style=3)
gdf_temp.plot()




#transformed coordinate system (CRS: EPSG:27700)
gdf = gdf.to_crs("EPSG:27700")
ax=gdf.plot(column='agcval',cmap='rainbow')
'''
#25km grid
xmin,ymin,xmax,ymax = gdf.total_bounds #bounds: (-9.01, 49.75, 2.01, 61.01)
width = 5000  
height = 5000
rows = int(np.ceil((ymax-ymin) /  height))
cols = int(np.ceil((xmax-xmin) / width))
XleftOrigin = xmin
XrightOrigin = xmin + width
YtopOrigin = ymax
YbottomOrigin = ymax- height
polygons = []
for x0 in np.arange(xmin,xmax+5000,5000):
    for y0 in np.arange(ymin,ymax+5000,5000):
        x1=x0-5000
        y1=y0+5000
        polygons.append(shapely.geometry.box(x0,y0,x1,y1))
'''
'''
for i in range(cols):
    Ytop = YtopOrigin
    Ybottom =YbottomOrigin
for j in range(rows):
    polygons.append(Polygon([(XleftOrigin, Ytop), (XrightOrigin, Ytop), (XrightOrigin, Ybottom), (XleftOrigin, Ybottom)])) 
    Ytop = Ytop - height
    Ybottom = Ybottom - height
    XleftOrigin = XleftOrigin + width
    XrightOrigin = XrightOrigin + width

grid25 = gpd.GeoDataFrame({'geometry':polygons})
grid25.to_excel('25grid.xlsx')
grid25.to_file("grid25.shp")
grid25=grid25.set_crs('EPSG:3395')
grid25.plot()
'''
#Spatial connections
gdf_grid25=gpd.sjoin(gdf_temp,gdf,how="left",predicate="within")
gdf_grid25.head()
gdf_grid25.to_excel('joined_test.xlsx')
gdf.to_excel('gdf.xlsx')
gdf_grid25.plot(column='agcval') 
#The data for "agcval" has not yet been merged

#Combining MCF and cattle number

 