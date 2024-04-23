using Emerald.EmeraldLand.CanopyOptics: canopy_optical_properties!, read_spectrum, soil_albedo!, shortwave_radiation!
using Emerald.EmeraldLand.LeafOptics: leaf_spectra!
using Emerald.EmeraldLand.Namespace: MultiLayerSPAC, SPACConfiguration
using Emerald.EmeraldLand.SPAC: initialize!, soil_plant_air_continuum!, update!
using Emerald.EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak
using Emerald.EmeraldMath.Stats: rmse

using Base.GC

# Set up the global variables
FT = Float64;
CONFIG = SPACConfiguration{FT}();
SHIFT = MultiLayerSPAC(CONFIG);
initialize!(CONFIG, SHIFT);
SHIFT_BAK = deepcopy(SHIFT);


# Get the target curve
function target_curve(ref_x::Vector{FT}, params::Dict{String,FT}) where {FT}
    SHIFT = deepcopy(SHIFT_BAK);

    # update the SPAC parameters
    _keys = keys(params);

    # cab
    if "cab" in _keys
        for _leaf in SHIFT.LEAVES
            _leaf.BIO.cab = params["cab"];
        end;
    end;

    if !("car" in _keys) && "cab" in _keys
        for _leaf in SHIFT.LEAVES
            _leaf.BIO.car = _leaf.BIO.cab / 7;
        end;
    end;

    # car
    if "car" in _keys
        for _leaf in SHIFT.LEAVES
            _leaf.BIO.car = params["car"];
        end;
    end;

    # cbc
    if "cbc" in _keys
        for _leaf in SHIFT.LEAVES
            _leaf.BIO.cbc = params["cbc"];
            _leaf.BIO.lma = _leaf.BIO.pro + _leaf.BIO.cbc;
        end;
    end;

    # ci
    if "ci" in _keys
        SHIFT.CANOPY.ci = params["ci"];
        SHIFT.CANOPY.Ω_A = params["ci"];
    end;

    # lai
    if "lai" in _keys
        update!(CONFIG, SHIFT; lai = params["lai"]);
    end;

    # lma
    if "lma" in _keys
        for _leaf in SHIFT.LEAVES
            _leaf.BIO.lma = params["lma"];
        end;
    end;

    # pro
    if "pro" in _keys
        for _leaf in SHIFT.LEAVES
            _leaf.BIO.pro = params["pro"];
            _leaf.BIO.lma = _leaf.BIO.pro + _leaf.BIO.cbc;
        end;
    end;

    # sc
    if "sc" in _keys
        SHIFT.SOIL.COLOR = Int(params["sc"]);
    end;

    # top soil mositure
    if "tsm" in _keys
        SHIFT.SOIL.LAYERS[1].θ = params["tsm"];
    end;

    # rerun the leaf_spectra!
    for _leaf in SHIFT.LEAVES
        _leaf.BIO._v_storage = 0;
    end;
    leaf_spectra!(CONFIG, SHIFT);

    # generate the canopy level spectrum
    # soil_plant_air_continuum!(SHIFT, CONFIG, FT(0));
    (; ANGLES, CANOPY, LEAVES, METEO, SOIL) = SHIFT;
    soil_albedo!(CONFIG, SOIL);
    canopy_optical_properties!(CONFIG, CANOPY, ANGLES);
    canopy_optical_properties!(CONFIG, CANOPY, LEAVES, SOIL);
    shortwave_radiation!(CONFIG, CANOPY, LEAVES, METEO.rad_sw, SOIL);

    # match the reference and target curves
    _tar_ys = similar(ref_x);
    _min_wl = minimum(CONFIG.WLSET.Λ);
    _max_wl = maximum(CONFIG.WLSET.Λ);
    for _i in eachindex(_tar_ys)
        _mask = (_min_wl <= ref_x[_i] <= _max_wl) && !(1790 <= ref_x[_i] <= 1920) && !(1345 <= ref_x[_i] <= 1415);
        _tar_ys[_i] = _mask ? read_spectrum(CONFIG.WLSET.Λ, SHIFT.CANOPY.RADIATION.albedo, ref_x[_i]) : FT(NaN);
    end;
    
    #GC.gc()
    return _tar_ys
end


# Define the target function
function target_function(ref_xy::Tuple{Vector{FT},Vector{FT}}, params::Dict{String,FT}) where {FT}
    @assert length(ref_xy[1]) == length(ref_xy[2]) "Length of reference spectrum does not match!";

    # match the reference and target curves
    _ref_ys = ref_xy[2];
    _tar_ys = target_curve(ref_xy[1], params)
    
    #GC.gc()
    return rmse(_ref_ys, _tar_ys)
end


# Function to return fitting required solver parameters by the input vars to fit
function solver_params(ref_xy::Tuple{Vector{FT},Vector{FT}}, vars::Vector{String} = ["cab"]; soil_color::Int = 1, top_soil_moisture::Number = 0.3)
    # define the var switch
    @inline switch_var(var::String) = (
        # return the x_min, x_max, x_ini, Δ_ini, and tol
        return if var == "cab"
            FT[0.01, 80, 10, 10, 0.1]
        elseif var == "car"
            FT[0.01, 80, 10, 10, 0.1]
        elseif var == "cbc"
            FT[0.001, 1, 0.012, 0.01, 0.001]
        elseif var == "ci"
            FT[0.01, 1, 1, 0.1, 0.01]
        elseif var == "lai"
            FT[0.01, 10, 3, 1, 0.01]
        elseif var == "lma"
            FT[1e-6, 1, 0.012, 0.01, 0.001]
        elseif var == "pro"
            FT[1e-6, 1, 0.012, 0.01, 0.001]
        else
            error("$(var) is not supported by switch_var function!")
        end;
    );

    # iterate the vars to get the x_mins, x_maxs, x_inis, and Δ_inis
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
        _params = Dict{String,FT}("sc" => soil_color, "tsm" => top_soil_moisture);
        for _i in eachindex(vars)
            push!(_params, vars[_i] => vals[_i]);
        end;
        
        #GC.gc()
        
        return _params
    );

    @inline _fit_func(vals::Vector{FT}) where {FT} = (
        _params = _dict_func(vals);
        
        #GC.gc()
        return -1 * target_function(ref_xy, _params)
    );
    
    #GC.gc()
    return _dict_func, _fit_func, _ms, _st
end;


# Function to fit the shift traits
function fit_shift_traits(ref_xy::Tuple{Vector{FT},Vector{FT}})
    if all(isnan, ref_xy[2][1:10])
        #GC.gc()
        return [NaN, NaN, NaN]
    end;

    (_, _fit_func, _ms, _st) = solver_params(ref_xy, ["cab", "lai", "lma"]; soil_color = 1, top_soil_moisture = 0.3);
    #GC.gc()
    return find_peak(_fit_func, _ms, _st)
end;
