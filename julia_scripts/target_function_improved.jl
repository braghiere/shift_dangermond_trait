include("/home/renatob/data/FluoData1/aviris_dangermond/review/test_field_data/target_function_free_lai_ci.jl")

# ─────────────────────────────────────────────────────────────────────────
# IMPROVED INVERSION — fixes for CHL ceiling (cab hitting 80 upper bound)
# Root cause: degeneracy between cab, sc (soil colour), and lai
#
# Idea 1 [J]: cab init 40 µg/cm² (field mean) instead of default 10
# Idea 2 [K]: cab search step 2 instead of 5 — finer, less overshoot
# Idea 3 [L]: soft Gaussian prior on cab centred at 40 µg/cm² (σ=40)
# Idea 4 [M]: red-edge bands 700–740 nm upweighted ×5 in spectral RMSE
# Config N: all 4 ideas combined
# ─────────────────────────────────────────────────────────────────────────

# Idea 4: weighted spectral RMSE — upweight 700–740 nm (cab red-edge fingerprint)
function target_function_rededge(ref_xy::Tuple, params::Dict{String,FT};
                                  band_wt::FT = FT(5.0)) where {FT}
    _ref_xs = ref_xy[1]
    _ref_ys = ref_xy[2]
    _tar_ys = target_curve(_ref_xs, params)
    valid   = .!isnan.(_ref_ys) .& .!isnan.(_tar_ys)
    sum(valid) < 2 && return FT(Inf)
    ws  = FT[(700.0 <= x <= 740.0) ? band_wt : FT(1.0) for x in _ref_xs]
    sw  = sum(ws[valid])
    return sqrt(sum(ws[valid] .* (_ref_ys[valid] .- _tar_ys[valid]) .^ 2) / sw)
end

# Parameterised improved solver — ideas 1–4 as Bool keyword flags
function fit_shift_traits_improved(ref_xy::Tuple;
    vars       = ["cab","cbc","lwc","lai","sc","n","ci"],
    lai        = 1.0,
    ci         = 1.0,
    idea1      = false,
    idea2      = false,
    idea3      = false,
    idea4      = false,
    reg_lambda = 0.01,
    band_wt    = 5.0,
)
    all(isnan, ref_xy[2][1:10]) && return fill(NaN, length(vars))

    @inline sw_var(v::String) = begin
        v == "cab" && return FT[0.01, 80.0, (idea1 ? 40.0 : 10.0), (idea2 ? 2.0 : 5.0), 0.1]
        v == "cbc" && return FT[0.0005, 0.035, 0.01, 0.005, 0.0005]
        v == "ci"  && return FT[0.2, 1.0, 0.7, 0.1, 0.01]
        v == "lai" && return FT[0.1, 10.0, 2.0, 0.5, 0.05]
        v == "lwc" && return FT[0.1, 20.0, 5.0, 1.0, 0.1]
        v == "n"   && return FT[1.0, 2.0, 1.4, 0.2, 0.05]
        v == "sc"  && return FT[1.0, 20.0, 10.0, 2.0, 1.0]
        error("$(v) not supported in sw_var")
    end

    _x_mins=FT[]; _x_maxs=FT[]; _x_inis=FT[]; _Δ_inis=FT[]; _tols=FT[]
    for v in vars
        p = sw_var(v)
        push!(_x_mins, p[1]); push!(_x_maxs, p[2]); push!(_x_inis, p[3])
        push!(_Δ_inis, p[4]); push!(_tols, p[5])
    end

    _ms = ReduceStepMethodND{FT}(x_mins=_x_mins, x_maxs=_x_maxs, x_inis=_x_inis, Δ_inis=_Δ_inis)
    _st = SolutionToleranceND{FT}(_tols, 50)

    cab_idx = findfirst(==("cab"), vars)

    @inline function _fit_func(vals::Vector{FT})
        p = Dict{String,FT}("sc"=>FT(10.0), "tsm"=>FT(0.3), "lai"=>FT(lai), "ci"=>FT(ci))
        for i in eachindex(vars); p[vars[i]] = vals[i]; end
        obj = if idea4
            target_function_rededge(ref_xy, p; band_wt=FT(band_wt))
        else
            target_function(ref_xy, p)
        end
        if idea3 && !isnothing(cab_idx)
            cab_v = vals[cab_idx]
            obj  += FT(reg_lambda) * max(FT(0), cab_v - FT(40))^2 / FT(1600)
        end
        return -obj
    end

    try
        return find_peak(_fit_func, _ms, _st)
    catch e
        @info "Error at site $(ref_xy[3]): $e"
        return fill(NaN, length(vars))
    end
end
