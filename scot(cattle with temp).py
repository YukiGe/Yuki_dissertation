# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 17:57:50 2023

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
import shapely
import xarray as xr

##Import modules for calculating slurry emissions (These are other python files in this folder - look in the .py files to see how they work!)
from BasicParams import * #BasicParams has standard parameters used for working out how much slurry is being generated.  This way of importing means you can access the paramters by using >>MCF_k rather than >>BasicParams.MCF_k
import LS_Energy #Uses IPCC calculations to work out how much energy and food a cow needs given its milk yeild, growth, temperture and how much walking it does.
import LS_Emissions #Uses the output of LS_Energy functions to work out enteric methane production and how much slurry is being generated.
import LS_Slurry_Methane_Temp as MethaneTemp #Uses how much slurry is being generated and temperture to work out methane production from slurry with a model.
import LS_SlurrySpreadTime as SpreadTime #Has a function for working out "voltaile solids" (i.e. slurry that can turn to methane) generation rate Can work out when is the best time to spread slurry on field or use AD to minimize methane emissions for that year.

def load_temp_nc_data(path,epsg): 
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

def str_to_dt(df):
    dates=list(df['time'])
    dates_dt=[datetime.datetime.strptime(d,'%Y-%m-%d') for d in dates]
    df['time']=dates_dt
    return df
    

if __name__ == '__main__':

    epsg_m=27700 #a CRS with units meters

    temp_data_folder='C:/Users/admin/Desktop/LiveStockCalculations-master/data/merge_gdf/2019/2019_merged.shp'
    merged_gdf_2019 = gpd.read_file(temp_data_folder) 
    
    scot_cattle_data_folder='../data/merge_gdf/2019/'
    gdf_cattle_2019=gpd.read_file('../data/merge_gdf/2019/2019_merged.shp').to_crs(epsg_m)
    merged_gdf_2019=str_to_dt(merged_gdf_2019 )
    merged_gdf_2019=merged_gdf_2019.drop(columns=['level_0'])
    merged_gdf_2019=merged_gdf_2019.fillna(0)
    column_converter={'projection_y_coordinate':'projection','projection_x_coordinate':'projecti_1','month_number':'month_numb'}   #transfer columns's name since saving the shp file in short name
'''
#plot for January:
ax=merged_gdf_2019.plot(column='agcval',cmap='rainbow')
gdf_2019_January =  merged_gdf_2019[merged_gdf_2019[column_converter['month_number']]==1] #Filter the data frame to only include rows which have 3 in the month_number column
gdf_2019_January.plot(column='tas',cmap='rainbow',markersize=15,legend=True,legend_kwds={'label': "Temperture [C]"}) #the other bits define the color map and add a label to the legend.
plt.show()
'''

#First we define some variables which determine what goes on in the farm:
month_in=11 #what month they are brought inside the barn
month_out=5 #what month they are let out of the barn
summer_prop=4/24 #what proportion of the day are they brought inside to the barn for milking in when out to pasture
winter_prop=24/24 #What proportion of the time they are in the barn. i.e. all the time 
DE_summer=0.7 #The 'Digestible Energy' of the food they are eating when out to pasture i.e. they grass in the field
DE_winter=0.75 #The 'Digestible Energy' of the food they are eating when in the barn over winter i.e. Hay and cereals.
milk_yield=9000/365 #The daily milk yield in kg of Milk/day

#Define how much weight the cows are gaining [kg/day]. assume all the dairy cows are not growing (WG_dict['dairy']=0)
WG_dict={'dairy':0,'dairy_nl':0,'MF':LS_Energy.calc_WG(birth_weight_var,MW_dict['MF'],slaughter_age_dict['MF']),'MM':LS_Energy.calc_WG(birth_weight_var,MW_dict['MM'],slaughter_age_dict['MM'])} 
#'dairy' = dairy cows; 'dairy_nl' = dairy not lactating; 'MF' =meat female; 'MM' = meat male.
#You can look at 'dairy' and 'dairy non lactaing - maybe some of the younger dairy cows are growing?

#extract each months tempertures&cattle for one grid point
y_coord=502500
x_coord=302500
 #these coordinates define the grid point we want to extract (you can choose another one if you want)
gdf_2019_point = merged_gdf_2019[(merged_gdf_2019[column_converter['projection_y_coordinate']] == y_coord) & (merged_gdf_2019[column_converter['projection_x_coordinate']] == x_coord)] #filter down the data frame to only include this point
months =list(gdf_2019_point[column_converter['month_number']]) #extract all the months (will be 1 to 12) and convert to a list
temps_monthly = list(gdf_2019_point['tas']) #extract all the tempertures and convert to a list
cattle_yearly = list(gdf_2019_point['agcval'])
start_date = (2019,1,1) #(year, month, day) 
 
#Just now we have the temperture for each month but the model takes daily resolution. We assume the the temperture on each day is the monthly mean:
temps_daily=np.array([]) #initialise an empty array - we will fill this up with daily tempertures using a loop:
for i in range(len(months)): #loop through each month
    months_daily_temps=np.ones(monthrange(start_date[0],months[i])[1])*temps_monthly[i] #Make an array of all ones with length equal to the number of days in the month then multiply it by the temperture of that month.
    temps_daily=np.concatenate((temps_daily,months_daily_temps)) #update temps_daily by 'concatenating' this months daily tempertures to it.
days=np.arange(0,len(temps_daily),1) #an array representing each day 0 to 365
#Work out how much volatile solids (i.e. slurry) is being generated by ONE cow as well as other metrics(for now, assume as dairy cow)
dates,GEs,VS_slurry,VS_pasture,Enteric_methane=SpreadTime.GE_VS_per_day(days,start_date,month_in,month_out,winter_prop,summer_prop,DE_summer,DE_winter,milk_yield,WG_dict,temp_days=temps_daily,cow_type='dairy') 
#The outputs are: dates in the simulation; Gross Energy required by a cow [MJ/day per cow]; Volatile solids going into slurry [kgVS/day per cow];
                    #Volatile solids going directly to field [kgVS/day per cow]; Enteric methane being generarted by the cow [kgCH4/day per cow]
#The function inputs are: days in the year;the start date of the simulation (year,month,day); month cows are brought in; month they are out to pasture; winter_prop,summer_prop; DE in summer;
                        #DE in winter; milk yield [kg/day]; daily weight dictionarry; the daily tempertures; what type of cow it is.
                        
                    
#Plot Gross Energy requirements for ALL cows per day:
plt.plot(dates,GEs*cattle_yearly[0]) #plot dates vs enteric methane
plt.ylabel('All cows gross energy requiremetns [MJ/day]')
plt.xlabel('Day')
plt.gcf().autofmt_xdate() #To make the dates look nice
plt.show() #You can see when the cows are outside in summer, they need more energy because they are walking more and the temperture can be cold.

#!!Have a go at plotting the volatile solid production rate and enteric methane below!!##


######Now we work out how much methane is being generated from the slurry by using the tempertures and the volatile solids production rate using the model ####
total_slurry_methane,total_VS_slurry,dM_gdt,VS,TM,M_s_at_empty,V_at_empty=SpreadTime.slurry_spread_time(days,MCF_m,B_0_dict['dairy'],temps_daily,VS_slurry*cattle_yearly[0],0.1,gauss_popt,spread_time=[])
#Outputs are: Total slurry methane generated over the simulations time [kgCH4]; Total amount of volatile solid generated over the year [kg VS]; Methane generated per kg of volatile solids [kgCH4/kgVS/day];
               #Amount of volatile solids stored on each day [kgVS]; Total Methane (TM) generated; Concentration of potentially methonisable solids left in slurry at end of simulation [kgCH4/kgVS]; Volume of volatils solids at empty [kgVS]
               
#Inputs are: days on the simulations; a parameter for the model (from BasicParams); daily tempertures; volatile solids production rate [kgVS/day]; initial amount of slurry [kg]; 
            #parameters controlling methogenesis (from BasicParams); spread_time=[] can be a list of days for the slurry store to be emptied. Just now it is never emptied.

#We can calculate the methane conversion factor (MCF) as the total methane emitted (TM) over the total amount of methane (emitted and potential) in the system MCF=TM/(B_0*0.67*VS) (0.67 is to convert from m3 of CH4 to kg CH4)
MCF=TM/(B_0_dict['dairy']*0.67*VS)*100 
 
##Plot slurry methane production rate as well as temperture:
fig,ax=plt.subplots(figsize=(9,6))
ax2=ax.twinx() #so we can have two axes (one for temperture)
ax.plot(dates,dM_gdt,color='green') #on the first axes plot methane production rate
ax.set_xlabel('Day')
ax.set_ylabel('Methane production rate [kgCH4/kgVS/day]',color='green')
plt.gcf().autofmt_xdate()
ax2.plot(dates,temps_daily,color='orange')#on the second axes plot temperture
ax2.set_ylabel('Temperture [C]',color='orange')
plt.show() #Temperture affects methane production rate but so does the concentration of volatile solids that can be converted to methane

##!Have a go at plotting total methane (TM), Volatile solids (VS) and MCF over time below!##

######################################################################ONLY PROCEED ONCE YOU UNDERSTAND THE ABOVE #################

####We now generalise the above process as a function so we can apply it to all grid points in the tempeture GIS data#####

def month_to_DOY_end(month_number,year):#takes a month number and returns the day of year it ends on. We will need this later for reporting results at the end of each month.
    month_days=[monthrange(year,m)[1] for m in range(1,month_number+1)]
    return sum(month_days)


def calc_VS_rate_one_point(point,gdf,month_in,month_out,winter_prop,summer_prop,DE_summer,DE_winter,milk_yield,WG_dict,cow_type='dairy'): #This generalises the first pointnum
    gdf_point = gdf[(gdf[column_converter['projection_y_coordinate']] == point[1]) & (gdf[column_converter['projection_x_coordinate']] == point[0])] #filter down the data frame to only include this point
    gdf_point=gdf_point.reset_index() #reindex the sub gdf
    yearly_cattle=list(gdf_point['agcval'])
    num_cows=yearly_cattle[0]
    months =list(gdf_point[column_converter['month_number']]) #extract all the months (will be 1 to 12) and convert to a list
    temps_monthly = list(gdf_point['tas']) #extract all the tempertures and convert to a list
    start_date = (gdf_point['time'][0].year,1,1) #extract the year from the gdf to define the start date as the first of january
    
    #Just now we have the temperture for each month but the model takes daily resolution. We assume the the temperture on each day is the monthly mean:
    temps_daily=np.array([]) #initialise an empty array - we will fill this up with daily tempertures using a loop:
    for i in range(len(months)): #loop through each month
        months_daily_temps=np.ones(monthrange(start_date[0],months[i])[1])*temps_monthly[i] #Make an array of all ones with length equal to the number of days in the month then multiply it by the temperture of that month.
        temps_daily=np.concatenate((temps_daily,months_daily_temps)) #update temps_daily by 'concatenating' this months daily tempertures to it.
    days=np.arange(0,len(temps_daily),1) #an array representing each day 0 to 365
    
    dates,GEs,VS_slurry,VS_pasture,Enteric_methane=SpreadTime.GE_VS_per_day(days,start_date,month_in,month_out,winter_prop,summer_prop,DE_summer,DE_winter,milk_yield,WG_dict,temp_days=temps_daily,cow_type=cow_type)
    
    #We integrate Enteric_methane [kgCH4/day/cow] to report the total methane at the end of each month in gdf_point [kgCH4]
    Enteric_methane_int=cumtrapz(Enteric_methane,days,initial=0)
    
    #We now have to add new columns to gdf_point representing the Enteric_methane_int at the end of each month. This is a bit technical working with pandas dataframes:
    end_of_month_indexes=[month_to_DOY_end(m,gdf_point['time'][0].year)-1 for m in np.arange(1,13,1)] #To pick out the indexes in the solution where the end of month is. -1 because python indexes from 0.
    Enteric_methane_int_ends=[Enteric_methane_int[i]*num_cows for i in end_of_month_indexes] #Enteric_methane_int at the end of each month
    
    months_strings=[str(m) for m in np.arange(1,13,1)] #each month as a string
    Enteric_methane_int_column=dict(zip(Enteric_methane_int_ends,months_strings))
    gdf_point['Enteric Methane kgCH4']=Enteric_methane_int_ends #add in in new column for Enteric Methane
    
    return dates,GEs*num_cows,VS_slurry*num_cows,VS_pasture*num_cows,Enteric_methane*num_cows,Enteric_methane_int*num_cows,temps_daily,gdf_point #return the output of GE_VS_per_day but we multiple by the number of cows to get the full farm. We also return the daily tempertures and the sub gdf


def calc_slurry_methane_one_point(gdf_point,MCF_shape,B_0,temps_daily,VS_slurry,initial_volume,gauss_params,spread_time=[]): 

    #We use the output of the above function as input to SpreadTime.slurry_spread_time()
    #We then add new columns to gdf_point with the the output of the model at the end of each month

    solved_gdf_point=gdf_point.copy() #This will be our output
    days=np.arange(0,len(temps_daily),1)
    total_slurry_methane,total_VS_slurry,dM_gdt,VS,TM,M_s_at_empty,V_at_empty=SpreadTime.slurry_spread_time(days,MCF_shape,B_0,temps_daily,VS_slurry,initial_volume,gauss_params,spread_time=spread_time) #Note to Dan Can be made faster by only solving on month ends
    #We can also calculate the Methane conversion factor using the output of the model:
    MCF = TM/(VS*B_0_dict['dairy']*0.67)*100 #multiply by 100 to make a percent
    
    #We now have to add new columns to gdf_point representing the output of the model at the end of each month. This is a bit technical working with pandas dataframes:
    end_of_month_indexes=[month_to_DOY_end(m,gdf_point['time'][0].year)-1 for m in np.arange(1,13,1)] #To pick out the indexes in the solution where the end of month is. -1 because python indexes from 0.
    VS_month_ends=[VS[i] for i in end_of_month_indexes]
    TM_month_ends=[TM[i] for i in end_of_month_indexes]
    dM_gdt_month_ends=[dM_gdt[i] for i in end_of_month_indexes]
    MCF_month_ends=[MCF[i] for i in end_of_month_indexes]
    
   
    months_strings=[str(m) for m in np.arange(1,13,1)] #each month as a string
    
    VS_column=dict(zip(VS_month_ends,months_strings)) #these dictionaries make sure its going in the correct month row    
    TM_column=dict(zip(TM_month_ends,months_strings))   
    dM_gdt_column=dict(zip(dM_gdt_month_ends,months_strings))
    
    
    MCF_column = dict(zip(MCF_month_ends,months_strings))
    
    #add columns to gdf_point
    solved_gdf_point['VS kg']=VS_month_ends
    solved_gdf_point['TM kgCH4']=TM_month_ends
    solved_gdf_point['dM_gdt kgCH4/kgVS/day']=dM_gdt_month_ends
    solved_gdf_point['MCF']=MCF_month_ends
   
    
    return total_slurry_methane,total_VS_slurry,dM_gdt,VS,TM,MCF,M_s_at_empty,V_at_empty,solved_gdf_point
    
    
def find_points(gdf): # extracts all grid points in the tempertue GIS data. We will loop over the output of this later.
    #will get a unique list of unique (projection_y_coordinate,projection_x_coordinate)
    all_tuples=[(row[column_converter['projection_x_coordinate']],row[column_converter['projection_y_coordinate']]) for i,row in gdf.iterrows()]#all points
    unique_tuples=list(set(all_tuples)) #list(set()) removes double entries to get a unqiue list of points.
    return unique_tuples


#Lets see how these functions work!

all_points=find_points(merged_gdf_2019) #Find all the points

#Determine the amount of VS for the first point in all_points. The datframe ex_gdf_point will have an column for enteric methane have a look at it:
ex_dates,ex_GEs,ex_VS_slurry,ex_VS_pasture,ex_Enteric_methane,ex_Enteric_methane_int,ex_temps,ex_gdf_point=calc_VS_rate_one_point(all_points[0],merged_gdf_2019,month_in,month_out,winter_prop,summer_prop,DE_summer,DE_winter,milk_yield,WG_dict,cow_type='dairy')
#Determine the slurry emissions and other metrics. solved_gdf_point will have new columns for slurry methane emissions.
ex_total_slurry_methane,ex_total_VS_slurry,ex_dM_gdt,ex_VS,ex_TM,ex_MCF,ex_M_s_at_empty,ex_V_at_empty,ex_solved_gdf_point=calc_slurry_methane_one_point(ex_gdf_point,MCF_m,B_0_dict['dairy'],ex_temps,ex_VS_slurry,0.1,gauss_popt,spread_time=[])

#Have a look at solved_gdf_point. It is a geo data frame which has monthly information on enteric and slurry methane emissions for only one point
#On the grid from the temperure GIS data.

####We now write a function to solve it on all grid points using the above functions#####
def calc_emissions_all_points(gdf,month_in,month_out,winter_prop,summer_prop,DE_summer,DE_winter,milk_yield,WG_dict,
                                MCF_shape,B_0,initial_volume,gauss_params,spread_time=[],cow_type='dairy'):
    
    all_points=find_points(gdf) #Extract all grid points from the Temperture Data frame
    num_points=len(all_points)
    all_solved_gdfs=[] #we will store all the single point geo data frames in this list
    counter=0
    for p in all_points: #loop through all points
        print('Calculating slurry methane missions for whole grid. On gridpoint {gp:d} out of {ap:d}'.format(gp=counter,ap=num_points),end='\r') #to keep track of calculation progress
        dates,GEs,VS_slurry,VS_pasture,Enteric_methane,Enteric_methane_int,temps,gdf_point = calc_VS_rate_one_point(p,gdf,month_in,month_out,
                                                                                            winter_prop,summer_prop,DE_summer,DE_winter,milk_yield,WG_dict,cow_type=cow_type) #Caluclatue slurry production rate for this points
        
        total_slurry_methane,total_VS_slurry,dM_gdt,VS,TM,MCF,M_s_at_empty,V_at_empty,solved_gdf_point=calc_slurry_methane_one_point(gdf_point,MCF_shape,B_0,temps,VS_slurry,
                                                                                                                                initial_volume,gauss_params,spread_time=spread_time) #Calclualte slurry methane emissions and out put in solved_gdf_point
        #print(solved_gdf_point)
        all_solved_gdfs.append(solved_gdf_point)#add the currents solved_gdf_point to the list
        counter=counter+1
    solved_gdf=pd.concat(all_solved_gdfs) #join all the of them together in a big data frame
    return solved_gdf
 

##No we apply the functoin to the 2020 temperture data (it might take a bit of time to run it is commented out just now. delete the ''' to uncomment. Remember to change the FIG_SAVE_DIR to where you want to save them)

FIG_SAVE_DIR='C:/Users/admin/Desktop/LiveStockCalculations-master/plot/'
solved_gdf_2019=calc_emissions_all_points(merged_gdf_2019,month_in,month_out,winter_prop,summer_prop,DE_summer,DE_winter,milk_yield,WG_dict,
                                MCF_m,B_0_dict['dairy'],0.1,gauss_popt,spread_time=[],cow_type='dairy')
def df_dt_to_string(df):
    dates=list(df['time'])
    dates_str=[d.strftime('%Y-%m-%d') for d in dates]
    df['time']=dates_str
    return df

solved_gdf_2019=df_dt_to_string(solved_gdf_2019)
solved_gdf_2019.to_file('C://Users//admin//Desktop//LiveStockCalculations-master//LiveStockCalculations-master//result//solved_2019.shp')
#Now we can plot the MCF, for example, in the UK for each month:
for m in np.arange(1,13): #loop through all months
    month_solved_gdf_2019=solved_gdf_2019[solved_gdf_2019[column_converter['month_number']]==m] #filter to the month of intrest
    month_solved_gdf_2019.plot(column='TM kgCH4',cmap='rainbow',markersize=30,legend=True,legend_kwds={'label': "TM kgCH4"},figsize=(5, 9)) #vmin/vmax sets the limits on the colorbar
    mean_TM=np.mean(month_solved_gdf_2019['TM kgCH4']) #Calculate the mean MCF for this month over the UK
    plt.axis('off')
    plt.title('Mean TM kgCH4 {M:1.2f}'.format(M=mean_TM))#Make the title of the plot the mean MCF over the UK
    plt.savefig(FIG_SAVE_DIR+'TM kgCH4_2019_'+str(m)+'.png',dpi=300) #Save the figure
    plt.close()


'''
#CHOOSE the years you want to solve for:
start_year=1884
end_year=2020   

NC_files=os.listdir(temp_data_folder) #get the names of all the GIS data in the folder
NC_files = [f for f in NC_files if f.split('.')[-1]=='nc'] #only include files in the foler that end in '.nc' incase other files are in here.

NC_files=[f for f in NC_files if float(f.split('_')[5][0:4])<=end_year and float(f.split('_')[5][0:4])>=start_year and float(f.split('_')[5][0:4])!=1907 and float(f.split('_')[5][0:4])!=1920] #filter by years. We also avoid some years were the data is corrupted.
years=[int(f.split('_')[5][0:4]) for f in NC_files] #to keep track of years

for i in range(len(years)): #loop through all the indexes of the years data
    print('On year {y:d} out of {total:d}'.format(y=i+1, total=len(years)))
    gdf_year = load_temp_nc_data (temp_data_folder+NC_files[i],epsg_m) #load the i'th years GIS data up.
    solved_gdf_year=calc_emissions_all_points(gdf_year,num_cows,month_in,month_out,winter_prop,summer_prop,DE_summer,DE_winter,milk_yield,WG_dict,
                                MCF_m,B_0_dict['dairy'],0.1,gauss_popt,spread_time=[],cow_type='dairy') #Solve the model on the grid.
    #make a folder for us to save the maps into for this year if it doesnt exist:
    save_folder=FIG_SAVE_DIR+str(years[i])+'//'
    if not os.path.isdir(save_folder): 
        os.makedirs(save_folder) #make save_folder if doesnt exist
    
    #you can also save the solved geodata frame for use later.
    solved_gdf_year=solved_gdf_year.drop(columns=['time','time_bnds']) # we drop the time column before saving since shape files dont like python date time objects. The month_number column serves the same purpose
    solved_gdf_year.to_file(save_folder+str(years[i])+'SolvedWithSource.shp') #we save it as a shape file - you could load this in ArcGIS if you want.
    
    #Now plot each month
    for m in np.arange(1,13): #loop through all months
        month_solved_gdf_year=solved_gdf_year[solved_gdf_year[column_converter['month_number']]==m] #filter to the month of intrest
        month_solved_gdf_year.plot(column='MCF',cmap='rainbow',markersize=30,vmin=0,vmax=40,legend=True,legend_kwds={'label': "MCF [%]"},figsize=(5, 9)) #vmin/vmax sets the limits on the colorbar
        mean_MCF=np.mean(month_solved_gdf_year['MCF'])
        plt.axis('off')
        plt.title('Mean MCF {M:1.2f}'.format(M=mean_MCF))#Make the title of the plot the mean MCF over the UK
        plt.savefig(save_folder+str(m)+'.png',dpi=300) #Save the figure
        plt.close()
'''         