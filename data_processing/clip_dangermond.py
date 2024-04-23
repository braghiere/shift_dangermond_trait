import xarray
import rioxarray
import geopandas
from shapely.geometry import mapping

#xds = xarray.open_dataset("output.nc", drop_variables='time_bnds')
xds = xarray.open_dataset("output.nc")
xds.rio.write_crs("EPSG:4326", inplace=True)
geodf = geopandas.read_file("/net/fluo/data1/students/renato/aviris_dangermond/dangermond_bound/dan_bound.shx")

clipped = xds.rio.clip(geodf.geometry.apply(mapping), geodf.crs,drop=False)
clipped.to_netcdf('output_clipped.nc')
