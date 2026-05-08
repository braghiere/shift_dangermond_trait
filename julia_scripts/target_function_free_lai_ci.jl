using Emerald.EmeraldLand.Namespace: BulkSPAC, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize_spac!, prescribe_traits!, read_spectrum, soil_plant_air_continuum!
using Emerald.EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak
using Emerald.EmeraldMath.Stats: rmse

using Base.GC

# Set up the global variables (disable SIF as it is not needed)
FT = Float64;
CONFIG = SPACConfiguration(Float64);
CONFIG.ENABLE_SIF = false;
SHIFT = BulkSPAC(CONFIG);
initialize_spac!(CONFIG, SHIFT);
SHIFT_BAK = deepcopy(SHIFT);


# Get the target curve
function target_curve(ref_x::Vector{FT}, params::Dict{String,FT}) where {FT}
    SHIFT = deepcopy(SHIFT_BAK);

    # update the SPAC parameters
    _keys = keys(params);

    # cab
    if "cab" in _keys
        for _leaf in SHIFT.plant.leaves
            _leaf.bio.trait.cab = params["cab"];
        end;
    end;

    if !("car" in _keys) && "cab" in _keys
        for _leaf in SHIFT.plant.leaves
            _leaf.bio.trait.car = _leaf.bio.trait.cab / 7;
        end;
    end;

    # car
    if "car" in _keys
        for _leaf in SHIFT.plant.leaves
            _leaf.bio.trait.car = params["car"];
        end;
    end;

    # cbc
    if "cbc" in _keys
        for _leaf in SHIFT.plant.leaves
            _leaf.bio.trait.cbc = params["cbc"];
            _leaf.bio.trait.lma = _leaf.bio.trait.pro + _leaf.bio.trait.cbc;
        end;
    end;

    # ci
    if "ci" in _keys
        prescribe_traits!(CONFIG, SHIFT; ci = params["ci"]);
    end;

    # lai
    if "lai" in _keys
        prescribe_traits!(CONFIG, SHIFT; lai = params["lai"]);
    end;

    # lma
    if "lma" in _keys
        for _leaf in SHIFT.plant.leaves
            _leaf.bio.trait.lma = params["lma"];
        end;
    end;

    # lwc
    if "lwc" in _keys
        for _leaf in SHIFT.plant.leaves
            _leaf.capacitor.trait.v_max = params["lwc"];
            _leaf.capacitor.state.v_storage = params["lwc"];
        end;
    end;

    # n (leaf mesophyll structural parameter, meso_n field in Emerald biophysics)
    if "n" in _keys
        for _leaf in SHIFT.plant.leaves
            _leaf.bio.trait.meso_n = params["n"];
        end;
    end;

    # pro
    if "pro" in _keys
        for _leaf in SHIFT.plant.leaves
            _leaf.bio.trait.pro = params["pro"];
            _leaf.bio.trait.lma = _leaf.bio.trait.pro + _leaf.bio.trait.cbc;
        end;
    end;

    # sc: discrete soil colour class 1-20; round to nearest integer on every call
    if "sc" in _keys
        SHIFT.soil_bulk.trait.color = round(Int, params["sc"]);
    end;

    # top soil moisture
    if "tsm" in _keys
        SHIFT.soils[1].state.θ = params["tsm"];
    end;

    # rerun the leaf_spectra!
    initialize_spac!(CONFIG, SHIFT);

    # generate the canopy level spectrum
    soil_plant_air_continuum!(CONFIG, SHIFT, 0);

    # match the reference and target curves
    _tar_ys = similar(ref_x);
    _min_wl = minimum(CONFIG.SPECTRA.Λ);
    _max_wl = maximum(CONFIG.SPECTRA.Λ);
    for _i in eachindex(_tar_ys)
        _mask = (_min_wl <= ref_x[_i] <= _max_wl) && !(1790 <= ref_x[_i] <= 1920) && !(1345 <= ref_x[_i] <= 1415);
        _tar_ys[_i] = _mask ? read_spectrum(CONFIG.SPECTRA.Λ, SHIFT.canopy.sensor_geometry.auxil.reflectance, ref_x[_i]) : FT(NaN);
    end;

    return _tar_ys
end


# Define the target function
function target_function(ref_xy::Tuple, params::Dict{String,FT}) where {FT}
    @assert length(ref_xy[1]) == length(ref_xy[2]) "Length of reference spectrum does not match!";

    _ref_ys = ref_xy[2];
    _tar_ys = target_curve(ref_xy[1], params)

    return rmse(_ref_ys, _tar_ys)
end


# Function to return fitting required solver parameters by the input vars to fit.
# lai and ci keyword args supply prescribed values used when "lai"/"ci" are NOT in vars.
# Supported vars: "cab", "car", "cbc", "ci", "lai", "lma", "lwc", "n", "pro"
function solver_params(ref_xy::Tuple, vars::Vector{String} = ["cab"]; soil_color::Int = 20, top_soil_moisture::Number = 0.3, lai::Number = 1.0, ci::Number = 1.0)
    # returns [x_min, x_max, x_ini, Δ_ini, tol]
    @inline switch_var(var::String) = (
        return if var == "cab"
            FT[0.01, 80, 10, 5, 0.1]
        elseif var == "car"
            FT[0.01, 80, 10, 5, 0.1]
        elseif var == "cbc"
            FT[0.0005, 0.035, 0.01, 0.005, 0.0005]
        elseif var == "ci"
            FT[0.2, 1.0, 0.7, 0.1, 0.01]
        elseif var == "lai"
            FT[0.1, 10.0, 2.0, 0.5, 0.05]
        elseif var == "lma"
            FT[1e-6, 0.05, 0.012, 0.01, 0.001]
        elseif var == "lwc"
            FT[0.1, 20, 5, 1, 0.1]
        elseif var == "n"
            FT[1.0, 2.0, 1.4, 0.2, 0.05]
        elseif var == "pro"
            FT[1e-6, 0.015, 0.005, 0.002, 0.0005]
        elseif var == "sc"
            FT[1.0, 20.0, 10.0, 2.0, 1.0]
        else
            error("$(var) is not supported by switch_var function!")
        end;
    );

    _x_mins = FT[];
    _x_maxs = FT[];
    _x_inis = FT[];
    _Δ_inis = FT[];
    _tols   = FT[];
    for _var in vars
        _params = switch_var(_var);
        push!(_x_mins, _params[1]);
        push!(_x_maxs, _params[2]);
        push!(_x_inis, _params[3]);
        push!(_Δ_inis, _params[4]);
        push!(_tols  , _params[5]);
    end;

    _ms = ReduceStepMethodND{FT}(x_mins = _x_mins, x_maxs = _x_maxs, x_inis = _x_inis, Δ_inis = _Δ_inis);
    _st = SolutionToleranceND{FT}(_tols, 50);

    @inline _dict_func(vals::Vector{FT}) where {FT} = (
        # Prescribed lai/ci are pre-filled; any entry in vars will overwrite its slot
        _params = Dict{String,FT}("sc" => soil_color, "tsm" => top_soil_moisture, "lai" => lai, "ci" => ci);
        for _i in eachindex(vars)
            push!(_params, vars[_i] => vals[_i]);
        end;
        return _params
    );

    @inline _fit_func(vals::Vector{FT}) where {FT} = (
        _params = _dict_func(vals);
        return -1 * target_function(ref_xy, _params)
    );

    return _dict_func, _fit_func, _ms, _st
end;


# Flexible fit function.
# vars  : free parameters, e.g. ["cab", "lma", "lwc", "cbc", "pro", "n"]
# lai   : prescribed LAI (ignored when "lai" is in vars)
# ci    : prescribed CI  (ignored when "ci"  is in vars)
# Returns a vector of fitted values in the same order as vars.
function fit_shift_traits(ref_xy::Tuple; vars::Vector{String} = ["cab", "lma", "lwc", "cbc", "pro"], lai::Number = 1.0, ci::Number = 1.0)
    if all(isnan, ref_xy[2][1:10])
        return fill(NaN, length(vars))
    end;

    (_, _fit_func, _ms, _st) = solver_params(ref_xy, vars; soil_color = 20, top_soil_moisture = 0.3, lai = lai, ci = ci);
    try
        return find_peak(_fit_func, _ms, _st)
    catch
        @info "Error encounter at site $(ref_xy[3]) and $(ref_xy[4])!"
        return fill(NaN, length(vars))
    end;
end;
