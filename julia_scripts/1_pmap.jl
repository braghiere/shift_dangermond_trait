using Emerald.EmeraldFrontier: simulation!
using Emerald.EmeraldLand.SPAC: update!
using Base.GC


function run_shift_simulation!(param::Vector)
    if any(isnothing, param)
        return NaN, NaN, NaN, NaN, NaN, NaN
    end;

    config = param[1];
    #spac = param[2];
    spac = deepcopy(param[2]);
    df = param[3];
    simulation!(config, spac, df; initialial_state = true);
    
    spac = nothing;
    
    GC.gc()

    return df.F_GPP[1], df.SIF740[1], df.SIF683[1], df.SIF757[1], df.SIF771[1], df.F_H2O[1]
end;
