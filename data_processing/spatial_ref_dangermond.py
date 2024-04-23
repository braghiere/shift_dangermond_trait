import xarray
import geopandas
from shapely.geometry import mapping

xds = xarray.open_dataset("TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll.nc")
geodf = geopandas.read_file("/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/dangermond_bound/dan_bound.shx")

# Set the spatial dimensions for the sif variable
xds.rio.set_spatial_dims("lon", "lat", inplace=True)

# Set the CRS information for the sif variable
xds.rio.write_crs("EPSG:4326", inplace=True)

#clipped = xds.rio.clip(geodf.geometry.apply(mapping), geodf.crs, drop=False)
xds.to_netcdf('TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_spatial_ref.nc')

