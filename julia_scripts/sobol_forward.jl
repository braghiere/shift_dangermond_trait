
include("/home/renatob/data/FluoData1/aviris_dangermond/review/test_field_data/target_function_free_lai_ci.jl")

# Extended forward model that also sets ant and brown pigments
function target_curve_full(ref_x::Vector{FT}, params::Dict{String,FT}) where {FT}
    SHIFT_loc = deepcopy(SHIFT_BAK)
    _keys = keys(params)

    for _leaf in SHIFT_loc.plant.leaves
        if "cab"   in _keys; _leaf.bio.trait.cab    = params["cab"];   end
        if "car"   in _keys; _leaf.bio.trait.car    = params["car"];   end
        if "ant"   in _keys; _leaf.bio.trait.ant    = params["ant"];   end
        if "brown" in _keys; _leaf.bio.trait.brown  = params["brown"]; end
        if "cbc"   in _keys; _leaf.bio.trait.cbc    = params["cbc"];   end
        if "pro"   in _keys; _leaf.bio.trait.pro    = params["pro"];   end
        if "n"     in _keys; _leaf.bio.trait.meso_n = params["n"];     end
        if "lwc"   in _keys
            _leaf.capacitor.trait.v_max     = params["lwc"]
            _leaf.capacitor.state.v_storage = params["lwc"]
        end
        # lma is always pro + cbc
        _leaf.bio.trait.lma = _leaf.bio.trait.pro + _leaf.bio.trait.cbc
    end
    if "lai" in _keys; prescribe_traits!(CONFIG, SHIFT_loc; lai=params["lai"]); end
    if "ci"  in _keys; prescribe_traits!(CONFIG, SHIFT_loc; ci=params["ci"]);   end
    if "sc"  in _keys; SHIFT_loc.soil_bulk.trait.color = Int(params["sc"]);      end
    if "tsm" in _keys; SHIFT_loc.soils[1].state.θ      = params["tsm"];          end

    initialize_spac!(CONFIG, SHIFT_loc)
    soil_plant_air_continuum!(CONFIG, SHIFT_loc, 0)

    _tar_ys = similar(ref_x)
    _min_wl = minimum(CONFIG.SPECTRA.Λ)
    _max_wl = maximum(CONFIG.SPECTRA.Λ)
    for _i in eachindex(_tar_ys)
        _mask = (_min_wl <= ref_x[_i] <= _max_wl) &&
                !(1790 <= ref_x[_i] <= 1920) && !(1345 <= ref_x[_i] <= 1415)
        _tar_ys[_i] = _mask ? read_spectrum(CONFIG.SPECTRA.Λ,
                               SHIFT_loc.canopy.sensor_geometry.auxil.reflectance,
                               ref_x[_i]) : FT(NaN)
    end
    return _tar_ys
end

using DelimitedFiles

input_csv  = ARGS[1]
output_csv = ARGS[2]
wl_path    = ARGS[3]

wl_raw = open(wl_path) do io
    inner = strip(strip(readline(io)), ['[',']'])
    parse.(Float64, split(inner, ','))
end

data, header = readdlm(input_csv, ',', header=true)
param_names = vec(header)
n_samples   = size(data, 1)

out_rows  = Vector{Vector{Float64}}()
wl_header = ["wl_$(round(Int, w))" for w in wl_raw]

for i in 1:n_samples
    row = data[i, :]
    params = Dict{String,FT}(
        n => Float64(row[findfirst(==(n), param_names)]) for n in param_names
    )
    params["sc"]  = Float64(round(Int, params["sc"]))
    params["tsm"] = FT(0.3)

    local rfl = fill(NaN, length(wl_raw))
    try
        rfl = target_curve_full(wl_raw, params)
    catch e
        @warn "Sample $i failed: $e"
    end
    push!(out_rows, rfl)
    i % 200 == 0 && @info "Sobol: $i/$n_samples"
end

open(output_csv, "w") do io
    println(io, join(wl_header, ","))
    for r in out_rows; println(io, join(r, ",")); end
end
@info "Sobol done -- $n_samples rows -> $output_csv"
