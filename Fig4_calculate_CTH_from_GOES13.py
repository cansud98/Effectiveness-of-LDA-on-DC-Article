#******************************************************
# imports
#******************************************************
from netCDF4 import Dataset
import numpy as np
import sys
import wrf
from tqdm import tqdm

wrffile = sys.argv[1]
goesfile = sys.argv[2]

wrffile = Dataset(wrffile,'r')
goesfile = Dataset(goesfile,'r')

lats = goesfile.variables['XLAT']
lats = np.asarray(lats)
lons = goesfile.variables['XLONG']
lons = np.asarray(lons)
brightnessT = goesfile['ch4'] #['calibratedData']


cldfra = np.asarray(wrffile.variables['CLDFRA'])
height = np.zeros_like(cldfra[0,:,:,:])

ph = np.asarray(wrffile.variables['PH'])
phb = np.asarray(wrffile.variables['PHB'])
geop = (ph + phb)/9.81

wrf_temp_p = wrf.getvar(wrffile,'tc')+273.15

# Convert z-staggered variables to mass coordinates
height = (geop[0,0:-1,:,:] + geop[0,1:,:,:]) / 2.0

#******************************************************
# GOES-13 CTH calculations
#******************************************************

#Required function
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

print('Calculating CTH...')

cldtophgt_G13 = np.zeros((height.shape[1],height.shape[2]))

for i in tqdm(range(height.shape[1])):
    for j in range(height.shape[2]):

        zindex_G13 = 0
        layer_number = 68

        while (layer_number>0):
            a = float(wrf_temp_p[layer_number,i,j])
            b = float(wrf_temp_p[layer_number-1,i,j])


            if a<brightnessT[0,i,j]<b:
                zindex_G13 = layer_number
                break
            layer_number = layer_number-1


        if (zindex_G13==0):

            c = find_nearest(wrf_temp_p[:,i,j],brightnessT[0,i,j])

            #print('(',i,',',j,')',lats[i],lons[j],'bt: ',brightnessT[0,i,j],'|','c: ',c )

            zindex_G13_mv   = 0
            layer_number_mv = 68

            while (layer_number_mv>0):
                if c == float(wrf_temp_p[layer_number_mv,i,j]):
                    zindex_G13_mv = layer_number_mv
                    break

                layer_number_mv = layer_number_mv-1

                cldtophgt_G13[i,j]=float(height[layer_number_mv,i,j])

        else:
            [x1,y1] = a,height[layer_number,i,j]
            [x,y] = float(brightnessT[0,i,j]),'nan'
            [x0,y0] = b,height[layer_number-1,i,j]

            y = (((y1-y0)/(x1-x0))*(x-x0)+y0)

            cldtophgt_G13[i,j] = y

print('Complete.')
'''
#******************************************************
# Save visualized GOES-13 CTH data in a .png format
#******************************************************

lat_min = 33
lat_max = 38
lon_min = -103
lon_max = -94

# specify your colormap and projection
mylevs = [-2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0]
cmap = plt.get_cmap('jet')
cmap.set_under(color='white')
crs = ccrs.PlateCarree()

# plot
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, facecolor='None', projection=crs)
ax.coastlines(resolution='10m', alpha=0.5)
ax.add_feature(cfeature.STATES)
plot_rain = ax.pcolormesh(lons, lats, cldtophgt_G13[:,:]/1000, vmin=3, vmax=20.0, cmap=cmap)
cbar = fig.colorbar(plot_rain)
cbar.ax.set_ylabel('LDA CLD TOP HGT (km)')

# set other plot parameters
plt.xlim((lon_min,lon_max))
plt.ylim((lat_min, lat_max))
plt.savefig(goesfile+"_CTH.png")
plt.show()
'''

#******************************************************
# Save calculated GOES-13 CTH data as a new NetCDF file
#******************************************************

try: ncfile.close()  # just to be safe, make sure dataset is not already open.
except: pass
goescthfile = "GOES13_CTH.nc"
ncfile = Dataset(goescthfile,mode='w',format='NETCDF4_CLASSIC')
print(ncfile)

lat_dim = ncfile.createDimension('lati', len(lats[:,0]))    # latitude axis
lon_dim = ncfile.createDimension('longi', len(lons[0,:]))    # longitude axis

for dim in ncfile.dimensions.items():
    print(dim)

ncfile.title='GOES-13 Cloud Top Height Data'
print(ncfile.title)

ncfile.subtitle="Calculated with wrfinput_d02"
print(ncfile.subtitle)
print(ncfile)

# Define two variables with the same names as dimensions,
# a conventional way to define "coordinate variables".
lati = ncfile.createVariable('lati',np.float64,('lati','longi'))
lati.units = 'degrees_north'
lati.long_name = 'lati'
longi = ncfile.createVariable('longi',np.float64,('lati','longi'))
longi.units = 'degrees_east'
longi.long_name = 'longi'

# Define a 2D variable to hold the data
GOES13_CTH = ncfile.createVariable('GOES13_CTH',np.float64,('lati','longi')) 
GOES13_CTH.units = 'm'
GOES13_CTH.standard_name = 'GOES-13 calculated CTH'
print(GOES13_CTH)

print("-- Some pre-defined attributes for variable GOES13_CTH:")
print("CTH.dimensions:", GOES13_CTH.dimensions)
print("CTH.shape:", GOES13_CTH.shape)
print("CTH.dtype:", GOES13_CTH.dtype)
print("CTH.ndim:", GOES13_CTH.ndim)

#Define lat, lon.
lati[:] = lats[:]
longi[:] = lons[:]

print("-- Wrote data, CTH.shape is now ", GOES13_CTH.shape)

#Fill CTH values with calculated ones
GOES13_CTH[:,:] = cldtophgt_G13[:,:]

# first print the Dataset object to see what we've got
print(ncfile)
# close the Dataset.
ncfile.close(); print('Dataset is closed!')
