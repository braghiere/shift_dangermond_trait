
include("/home/renatob/data/FluoData1/aviris_dangermond/review/test_field_data/target_function_free_lai_ci.jl")

config     = ARGS[1]
input_csv  = ARGS[2]
output_csv = ARGS[3]

function parse_float_vec(s::AbstractString)
    inner = strip(strip(s), ['[', ']'])
    return parse.(Float64, split(inner, ','))
end

function parse_float_field(s::AbstractString)
    s2 = strip(s)
    return isempty(s2) ? NaN : parse(Float64, s2)
end

# lma is NOT a free variable: in PROSPECT lma = pro + cbc (they are not independent).
# We fit pro and cbc separately and compute lma = pro + cbc as output.
CONFIG_VARS = Dict(
    "A" => ["cab","lwc","cbc","pro"],
    "B" => ["cab","lwc","cbc","pro"],
    "C" => ["cab","lwc","cbc","pro","n"],
    "D" => ["cab","lwc","cbc","pro","lai"],
    "E" => ["cab","lwc","cbc","pro","n","lai"],
)

vars_free = CONFIG_VARS[config]

rows   = readlines(input_csv)
header = split(rows[1], ",")
idx_col = findfirst(h -> h == "sample_idx", header)
lai_col = findfirst(h -> h == "lai_sr",     header)
ci_col  = findfirst(h -> h == "ci_map",     header)

# Output: free params + derived lma (pro+cbc) + rmse
out_header = vcat(["sample_idx"],
                  ["$(config)_$(v)" for v in vars_free],
                  ["$(config)_lma", "$(config)_rmse"])
out_rows = Vector{Vector{Any}}()

pro_idx = findfirst(==("pro"), vars_free)
cbc_idx = findfirst(==("cbc"), vars_free)

for (ri, row_str) in enumerate(rows[2:end])
    m_wl  = match(r"(\[.+?\])", row_str)
    m_rfl = match(r"(\[.+?\])", row_str, m_wl.offset + length(m_wl.match))
    wl_str  = m_wl.match
    rfl_str = m_rfl.match
    stub    = replace(replace(row_str, wl_str => "WL__", count=1),
                               rfl_str => "RFL__", count=1)
    fields  = split(stub, ",")

    sample_idx = parse(Int,     strip(fields[idx_col]))
    lai_val    = parse_float_field(fields[lai_col])
    ci_val     = parse_float_field(fields[ci_col])

    wl  = Vector{FT}(parse_float_vec(wl_str))
    rfl = Vector{FT}(parse_float_vec(rfl_str))

    p_lai = (config == "A") ? 1.0 : (isnan(lai_val) ? 2.0 : lai_val)
    p_ci  = (config == "A") ? 1.0 : (isnan(ci_val)  ? 0.5 : ci_val)

    ref_xy = (wl, rfl, Float64(sample_idx), 0.0)
    fitted = fit_shift_traits(ref_xy; vars=vars_free, lai=p_lai, ci=p_ci)

    # Derived lma = pro + cbc (physically: lma is always pro+cbc in PROSPECT)
    lma_val = (!isnothing(pro_idx) && !isnothing(cbc_idx)) ? fitted[pro_idx] + fitted[cbc_idx] : NaN

    # Spectral RMSE
    local fit_rmse = NaN
    if !any(isnan, fitted)
        fp = Dict{String,FT}("sc"=>20.0,"tsm"=>0.3,"lai"=>p_lai,"ci"=>p_ci)
        for (k,fv) in zip(vars_free, fitted); fp[k] = fv; end
        try
            rfl_fit  = target_curve(wl, fp)
            valid    = .!isnan.(rfl) .& .!isnan.(rfl_fit)
            fit_rmse = sum(valid) < 2 ? NaN : sqrt(sum((rfl[valid].-rfl_fit[valid]).^2)/sum(valid))
        catch e
            @warn "RMSE failed for sample $(sample_idx): $e"
        end
    end

    push!(out_rows, vcat([sample_idx], fitted, [lma_val, fit_rmse]))
    ri % 50 == 0 && @info "Config $(config): $(ri)/$(length(rows)-1) done"
end

open(output_csv, "w") do io
    println(io, join(out_header, ","))
    for r in out_rows; println(io, join(r, ",")); end
end
@info "Config $(config) done -- $(length(out_rows)) rows"
