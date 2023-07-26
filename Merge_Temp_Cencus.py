##Import standard python modules:
import numpy as np #For array data types and basic math calculations
import matplotlib.pyplot as plt #For making plots
import pandas as pd #For handling 'data frames' (essentially spreadsheets), loading spreadsheets
import geopandas as gpd #For storing GIS data as a data frame. Essentially the same as a pandas data frame but has a column descriping the 'shape'
import datetime  #For handling dates and times e.g. can convert seconds to a dates. 
from calendar import monthrange #to work out how many days are in a month
import os #to do operating system commands
import shapely
import xarray as xr
def load_temp_nc_data(path,epsg): 
    '''
    Will load the .nc files sourced from Met Office; Hollis, D.; McCarthy, M.; Kendon, M.; Legg, T.; Simpson, I. (2020): HadUK-Grid Gridded Climate Observations on a 25km grid over the UK, v1.0.2.1 (1862-2008).
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
'''
def our_clip(gdf_big,gdf_small):
    polys=list(gdf_small['geometry'])
    scotland_poly=shapely.MultiPolygon(polys)
    scotland_poly_simple=scotland_poly.simplify(12000)
    indices_to_drop=[]
    for indx,row in gdf_big.iterrows():
        square=row.geometry
        if square.disjoint(scotland_poly_simple):
            indices_to_drop.append(indx)
    clipped=gdf_big.drop(indices_to_drop,axis=0)  
    clipped=clipped.reset_index()
    return clipped
'''
def our_clip_latitude(gdf_big,latitude):
    gdf_big_clip=gdf_big[gdf_big['latitude']>latitude]
    clipped=gdf_big_clip.reset_index()
    return clipped

def down_sample(old_grid,new_grid,col_name):
    #will take old_grid (gdf) and resample col_name onto new_grid. We assume old_grid is finer than new grids
    cattle_new=[]
    for indx,row in new_grid.iterrows():
        big_square = row.geometry
        cropped_old_grid=old_grid[old_grid.intersects(big_square)]
        cattle_in_big=0
        for indx_little,row_little in cropped_old_grid.iterrows():
            little_square = row_little.geometry
            inter=shapely.intersection(little_square,big_square)  
            inter_area=inter.area
            cattle_in_big=cattle_in_big+row_little[col_name]*inter_area/(little_square.area)
        cattle_new.append(cattle_in_big)
    new_grid['agcval']=cattle_new
    return new_grid

def df_dt_to_string(df):
    dates=list(df['time'])
    dates_str=[d.strftime('%Y-%m-%d') for d in dates]
    df['time']=dates_str
    return df

if __name__ == '__main__':
            
    epsg_m=27700 #a CRS with units meters

    temp_data_folder='../data/HadUk-Grid5km-1961-2020_Monthly/' 
    gdf_temp_2008 = load_temp_nc_data (temp_data_folder+'tas_hadukgrid_uk_5km_mon_200801-200812.nc',epsg_m) 
    gdf_temp_2008['geometry']=gdf_temp_2008.buffer(5000,cap_style=3)
    cencus_dairy_1='../data/Download_2295053/scotland-agcensus-2008_5km-grid_5103000/cattle/sc-2008-c81-5km-dairy-cows-in-milk.shp'
    cencus_dairy_2='../data/Download_2295053/scotland-agcensus-2008_5km-grid_5103000/cattle/sc-2008-c83-5km-dairy-cows-in-calf.shp'
    cencus_dairy_3='../data/Download_2295053/scotland-agcensus-2008_5km-grid_5103000/cattle/sc-2008-c85-5km-dairy-heifers-2-years-and-over.shp'
    cencus_dairy_4='../data/Download_2295053/scotland-agcensus-2008_5km-grid_5103000/cattle/sc-2008-c87-5km-dairy-heifers-under-2-years.shp'
    
    gdf_cattle_1=gpd.read_file(cencus_dairy_1).to_crs(epsg_m)
    gdf_cattle_2=gpd.read_file(cencus_dairy_2).to_crs(epsg_m)   
    gdf_cattle_3=gpd.read_file(cencus_dairy_3).to_crs(epsg_m)
    gdf_cattle_4=gpd.read_file(cencus_dairy_4).to_crs(epsg_m)
    gdf_cattle=gdf_cattle_1.copy()
    gdf_cattle['agcval']=gdf_cattle_1['agcval']+gdf_cattle_2['agcval']+gdf_cattle_3['agcval']+gdf_cattle_4['agcval']
    SCOTLAND_gdf_temp_2008= our_clip_latitude( gdf_temp_2008, 54.4)
    merged_gdf=down_sample(gdf_cattle,SCOTLAND_gdf_temp_2008,'agcval')
    ax=merged_gdf.plot(column='agcval',cmap='rainbow')
    merged_gdf_nodate=merged_gdf.drop(columns=['time_bnds'])
    merged_gdf_nodate=df_dt_to_string(merged_gdf_nodate)
    merged_gdf_nodate.to_file('C:/Users/admin/Desktop/LiveStockCalculations-master/data/merge_gdf/2008/2008_merged.shp')