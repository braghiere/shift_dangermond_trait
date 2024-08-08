using NetcdfIO: append_nc!, create_nc!, read_nc
# using Emerald.EmeraldVisualization: canvas, save_canvas!

using Distributed: @everywhere, addprocs, pmap, workers
using ProgressMeter: @showprogress

using Base.GC

using Printf

# you can change the number of threads here
# Fluo: within 48
# Tofu: within 48
# Curry: within 256
if length(workers()) == 1
    #addprocs(1; exeflags = "--project");
    addprocs(48; exeflags = "--project");
end;

@everywhere include("target_function_v3.jl");

function fit_shift_traits!(datafile::String, ncresult::String)
    _wlen = read_nc(FT, datafile, "wavelength")
    _refl = read_nc(FT, datafile, "reflectance")
    _lat_values = read_nc(FT, datafile, "latitude")  # Assuming 'latitude' is the variable name
    _lon_values = read_nc(FT, datafile, "longitude")  # Assuming 'longitude' is the variable name

    # Clip reflectance data to be within the range [0, 1]
    _refl = clamp.(_refl, 0.0, 1.0)

    # create a vector of parameters
    _params = []
    #for _i in axes(_refl,1)[1600:1610], _j in axes(_refl,2)[1600:1610] # change the index number of the lat and lon
    for _i in axes(_refl, 1), _j in axes(_refl, 2) # change the index number of the lat and lon
        push!(_params, (deepcopy(_wlen), _refl[_i, _j, :, 1], _i, _j))
    end

    _fittings = (@showprogress pmap(fit_shift_traits, _params))
    _chls = [_f[1] for _f in _fittings]
    _lais = [_f[2] for _f in _fittings]
    _lmas = [_f[3] for _f in _fittings]
    _lwcs = [_f[4] for _f in _fittings]
    _cbcs = [_f[5] for _f in _fittings]
    _pros = [_f[6] for _f in _fittings]

    _len_lat = size(_refl, 2)
    _len_lon = size(_refl, 1)
    #_len_lat = length(axes(_refl,2)[1600:1610]);
    #_len_lon = length(axes(_refl,1)[1600:1610]);
    _mat_chls = reshape(_chls, _len_lat, _len_lon)
    _mat_lais = reshape(_lais, _len_lat, _len_lon)
    _mat_lmas = reshape(_cis, _len_lat, _len_lon)
    _mat_lwcs = reshape(_lwcs, _len_lat, _len_lon)
    _mat_cbcs = reshape(_cbcs, _len_lat, _len_lon)
    _mat_pros = reshape(_pros, _len_lat, _len_lon)

    # save the matrices to a file
    create_nc!(ncresult, String["lon", "lat"], [_len_lon, _len_lat])  
    append_nc!(ncresult, "lat", _lat_values, Dict("latitude" => "latitude"), ["lat"])
    append_nc!(ncresult, "lon", _lon_values, Dict("longitude" => "longitude"), ["lon"])
    #append_nc!(ncresult, "lat", _lat_values[axes(_refl,2)[1600:1610]], Dict("latitude" => "latitude"), ["lat"]);
    #append_nc!(ncresult, "lon", _lon_values[axes(_refl,1)[1600:1610]], Dict("longitude" => "longitude"), ["lon"]);
    append_nc!(ncresult, "chl", _mat_chls, Dict("chl" => "chl"), ["lon", "lat"])
    append_nc!(ncresult, "lai", _mat_lais, Dict("lai" => "lai"), ["lon", "lat"])
    append_nc!(ncresult, "lma", _mat_cbcs + _mat_pros, Dict("lma" => "lma"), ["lon", "lat"])
    append_nc!(ncresult, "lwc", _mat_lwcs, Dict("lwc" => "lwc"), ["lon", "lat"])
    append_nc!(ncresult, "cbc", _mat_cbcs, Dict("cbc" => "cbc"), ["lon", "lat"])
    append_nc!(ncresult, "pro", _mat_pros, Dict("pro" => "pro"), ["lon", "lat"])
    #append_nc!(ncresult, "ci", _mat_cis, Dict("ci" => "ci"), ["lon", "lat"])
    
    GC.gc()

    return _fittings
end

# Loop over files from time_00 to time_12
function process_all_times()
    base_datafile = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/aviris_dangermond_time_%02d/output_clipped.nc"
    base_ncresult = "shift_fittings_lwc_time_%02d_test.nc"
    
    for time in 0:6
        #datafile = @sprintf(base_datafile, time)
        datafile = @sprintf("/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/aviris_dangermond_time_%02d/output_clipped.nc", time)
        #ncresult = @sprintf(base_ncresult, time)
        ncresult = @sprintf("shift_fittings_lwc_time_%02d_test.nc", time)
        fit_shift_traits!(datafile, ncresult)
    end
end

# Run the function to process all times
process_all_times()

