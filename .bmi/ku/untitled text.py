# Sampling soil data:
#
BD_file = 'soil/BD5min2.nc';
#
fh = Dataset(BD_file, mode='r')
lon_grid_BD = fh.variables['lon'][:]; 
lat_grid_BD = fh.variables['lat'][:]; 

# Make sure the latitude is strict monotonic increasing sequence
# lat_grid_BD = np.flipud(lat_grid_BD) 

x_coords, y_coords = xxxx.Mapping_Grid_Coords(lat_list, lon_list, lat_grid_BD, lon_grid_BD);

BD1 = fh.variables['BD'][:]

BD1 = BD1[y_coords, x_coords]
