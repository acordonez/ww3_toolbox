# this script is to combine NH and SH sea ice inputs
# into one file to use as an input for WW3

from netCDF4 import Dataset
import numpy as np

casename = 'b.e11.BRCP85C5CNBDRD.f09_g16.016.cice.h1.'
datadir = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/'
workdir = '/glade/scratch/aordonez/tmp/'
yyyy = 0 # year 0 in CESM lens file, corresponding to 2006
ryyyy = 2006
yyyy_str= '20060101-20801231'
ryyyy_str= '2006'
# output directory
outdir='.'
# output file name
outdata = casename + '_cice_' + yyyy_str + '.nc'
outdata2 = 'cice.nc'
outdata3 = 'wind.nc'

nj_tmp = 384
nj_sh = 76
nj_nh = 104
ni_tmp = 320
ndays = 365

myvar = 'aice_d'
# combine northern and southern hemispheres
fname1 = datadir + myvar + '/' + casename + myvar + '_nh.' +yyyy_str + '.nc'
fname2 = datadir + myvar + '/' + casename + myvar + '_sh.' +yyyy_str + '.nc'
infile1 = Dataset(fname1)
infile2 = Dataset(fname2)

aice1 = infile1.variables[myvar][0:ndays,:,:]
aice2 = infile2.variables[myvar][0:ndays,:,:]
fill_value = aice1.fill_value

aice_d = np.ones((ndays,nj_tmp,ni_tmp)) * fill_value
aice_d[:,384-nj_nh:384,:] = aice1
aice_d[:,0:nj_sh,:] = aice2

myvar = 'hi_d'
# combine northern and southern hemispheres
fname3 = datadir + myvar + '/' + casename + myvar + '_nh.' +yyyy_str + '.nc'
fname4 = datadir + myvar + '/' + casename + myvar + '_sh.' +yyyy_str + '.nc'
infile3 = Dataset(fname3)
infile4 = Dataset(fname4)

hi1 = infile3.variables[myvar][0:ndays,:,:]
hi2 = infile4.variables[myvar][0:ndays,:,:]

hi_d = np.ones((ndays,nj_tmp,ni_tmp)) * fill_value
hi_d[:,384-nj_nh:384,:] = hi1
hi_d[:,0:nj_sh,:,] = hi2

# ocean data file for grids
fname5 = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ocn/proc/tseries/daily/SST/b.e11.BRCP85C5CNBDRD.f09_g16.016.pop.h.nday1.SST.20060102-20801231.nc'
infile5 = Dataset(fname5)

# write the new netcdf
dataset = Dataset(outdata2,'w',format='NETCDF4_CLASSIC')
ni = dataset.createDimension('ni',ni_tmp)
nj = dataset.createDimension('nj',nj_tmp)
time = dataset.createDimension('time',ndays)
times = dataset.createVariable('time',np.float32,('time',))
TLON = dataset.createVariable('TLON',np.float32,('nj','ni',))
TLAT = dataset.createVariable('TLAT',np.float32,('nj','ni'))
ULON = dataset.createVariable('ULON',np.float32,('nj','ni',))
ULAT = dataset.createVariable('ULAT',np.float32,('nj','ni'))
ANGLE = dataset.createVariable('ANGLE',np.float32,('nj','ni'))
hi = dataset.createVariable('hi',np.float64,('time','nj','ni',),fill_value = fill_value)
aice = dataset.createVariable('aice',np.float64,('time','nj','ni',), fill_value = fill_value)

times[:] = infile1.variables['time'][0:ndays] - 2000 * 365
TLON[:] = infile5.variables['TLONG'][:]
TLAT[:] = infile5.variables['TLAT'][:]
ULON[:] = infile5.variables['ULONG'][:]
ULAT[:] = infile5.variables['ULAT'][:]
ANGLE[:] = infile5.variables['ANGLE'][:]
hi[:] = hi_d[:]
aice[:] = aice_d[:]


times.calendar = 'gregorian'
times.long_name = infile1.variables['time'].long_name
times.units = "days since 2000-01-01 00:00:00"
dataset.description = 'Test input file generated from LENS for WW3'
dataset.history = 'Created Nov 2017'
dataset.source = fname1 + '\n' + fname2 + '\n' + fname3 + '\n' + fname4 + '\n' + fname5
dataset.close()

"""# uniform floe size file
floedataset = Dataset('floesize.nc','w',format='NETCDF4_CLASSIC')
ni = floedataset.createDimension('ni',ni_tmp)
nj = floedataset.createDimension('nj',nj_tmp)
time = floedataset.createDimension('time',1)
times = floedataset.createVariable('time',np.float32,('time',))
TLON = floedataset.createVariable('TLON',np.float32,('nj','ni',))
TLAT = floedataset.createVariable('TLAT',np.float32,('nj','ni'))
ULON = floedataset.createVariable('ULON',np.float32,('nj','ni',))
ULAT = floedataset.createVariable('ULAT',np.float32,('nj','ni'))
ANGLE = floedataset.createVariable('ANGLE',np.float32,('nj','ni'))
floe = floedataset.createVariable('floe',np.float64,('time','nj','ni',),fill_value = fill_value)
times[:] = infile1.variables['time'][0:1]
TLON[:] = infile5.variables['TLONG'][:]
TLAT[:] = infile5.variables['TLAT'][:]
ULON[:] = infile5.variables['ULONG'][:]
ULAT[:] = infile5.variables['ULAT'][:]
ANGLE[:] = infile5.variables['ANGLE'][:]
# 10 meter ice floes
floesize = np.zeros(aice_d[0:1,:,:].shape)
floesize[(aice_d[0:1,:,:] > 0) & (aice_d[0:1,:,:] < 101)] = 1000.
floe[:] = floesize
floe.units = 'm'


times.calendar = 'gregorian'
times.long_name = infile1.variables['time'].long_name
times.units = "days since 2000-01-01 00:00:00"
floedataset.description = 'Test input file generated from LENS for WW3'
floedataset.history = 'Created Nov 2017'
floedataset.source = fname2 + '\n' + fname3 + '\n' + fname4 + '\n' + fname5
floedataset.close()"""

# non-uniform floe size file
floedataset = Dataset('floesize.nc','w',format='NETCDF4_CLASSIC')
ni = floedataset.createDimension('ni',ni_tmp)
nj = floedataset.createDimension('nj',nj_tmp)
time = floedataset.createDimension('time',1)
times = floedataset.createVariable('time',np.float32,('time',))
TLON = floedataset.createVariable('TLON',np.float32,('nj','ni',))
TLAT = floedataset.createVariable('TLAT',np.float32,('nj','ni'))
ULON = floedataset.createVariable('ULON',np.float32,('nj','ni',))
ULAT = floedataset.createVariable('ULAT',np.float32,('nj','ni'))
ANGLE = floedataset.createVariable('ANGLE',np.float32,('nj','ni'))
floe = floedataset.createVariable('floe',np.float64,('time','nj','ni',),fill_value = fill_value)
times[:] = infile1.variables['time'][0:1]
TLON[:] = infile5.variables['TLONG'][:]
TLAT[:] = infile5.variables['TLAT'][:]
ULON[:] = infile5.variables['ULONG'][:]
ULAT[:] = infile5.variables['ULAT'][:]
ANGLE[:] = infile5.variables['ANGLE'][:]
# 10 meter ice floes
#floesize = np.zeros(aice_d[0:1,:,:].shape)
floesize = aice_d[0:1,:,:] #* np.random.rand(aice_d[0:1,:,:].shape[0],aice_d[0:1,:,:].shape[1],aice_d[0:1,:,:].shape[2])
floesize[(aice_d[0:1,:,:] <= 0.001)] = 0
floesize[np.isnan(aice_d[0:1,:,:])] = np.nan
floe[:] = floesize
floe.units = 'm'


times.calendar = 'gregorian'
times.long_name = infile1.variables['time'].long_name
times.units = "days since 2000-01-01 00:00:00"
floedataset.description = 'Test input file generated from LENS for WW3'
floedataset.history = 'Created Nov 2017'
floedataset.source = fname2 + '\n' + fname3 + '\n' + fname4 + '\n' + fname5
floedataset.close()

#infile1.close()
infile2.close()
infile3.close()
infile4.close()
infile5.close()

#using wind stress instead of wind to test if grid is problem
fname6 = '/glade/scratch/aordonez/b.e11.BRCP85C5CNBDRD.f09_g16.016.cam.h1.UBOT.20060101-20801231.nc'
fname7 = '/glade/scratch/aordonez/b.e11.BRCP85C5CNBDRD.f09_g16.016.cam.h1.VBOT.20060101-20801231.nc'
infile6 = Dataset(fname6)
infile7  =Dataset(fname7)
u = infile6.variables['UBOT'][0:ndays,:,:]
v = infile7.variables['VBOT'][0:ndays,:,:]

lat = 192
lon = 288

winddataset = Dataset(outdata3,'w',format='NETCDF4_CLASSIC')
nlat = winddataset.createDimension('lat',lat)
nlon = winddataset.createDimension('lon',lon)
time = winddataset.createDimension('time',ndays)
times = winddataset.createVariable('time',np.float32,('time',))
LON = winddataset.createVariable('lon',np.float32,('lon',))
LAT = winddataset.createVariable('lat',np.float32,('lat',))
vbot = winddataset.createVariable('UBOT',np.float32,('time','lat','lon',))
ubot = winddataset.createVariable('VBOT',np.float32,('time','lat','lon',))

times[:] = infile1.variables['time'][0:ndays] - 2000 * 365
LON[:] = infile6.variables['lon'][:]
LAT[:] = infile6.variables['lat'][:]
ubot[:] = u[:]
vbot[:] = v[:]

ubot.units = infile6.variables['UBOT'].units
vbot.units = infile7.variables['VBOT'].units
ubot.long_name = infile6.variables['UBOT'].long_name
vbot.long_name = infile7.variables['VBOT'].long_name
ubot.cell_methods = infile6.variables['UBOT'].cell_methods
vbot.cell_methods = infile7.variables['VBOT'].cell_methods
LON.units = infile6.variables['lon'].units
LAT.units = infile6.variables['lat'].units

times.calendar = 'gregorian'
times.long_name = infile1.variables['time'].long_name
times.units = "days since 2000-01-01 00:00:00"
winddataset.description = 'Test input file generated from LENS for WW3'
winddataset.history = 'Created Nov 2017'
winddataset.source = fname1 + '\n' + fname2 + '\n' + fname3 + '\n' + fname4 + '\n' + fname5
winddataset.close()

infile6.close()
infile7.close()




