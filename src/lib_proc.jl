#= Data processing functions =#
using Arrow, YAML
export codf, standard, wo_s

"""
Standard long form processor
"""
function standard(it)::NamedTuple

    stm, flm, adj_mort, ocm = it

    # Extract state and flow
    stdict = munge_state(stm)
    fldict = munge_flow(flm)

    # Exception - extract treatments
    tx = stm[5, :, :]

    # Output dictionary
    ost::Dict{String,DataFrame} = Dict()
    ofl::Dict{String,DataFrame} = Dict()

    # Total population
    tot_pop::DataFrame =
        codf(permutedims(dropdims(sum(stm, dims = 1), dims = 1)), "Pop", opt = "pi")

    mean_pop::DataFrame =
        codf(permutedims(dropdims(sum(stm, dims = 1), dims = 1)), "Pop", opt = "m")

    @inbounds for (sn, mt) in stdict
        ost[sn] = innerjoin(codf(mt, sn, opt = "m"), mean_pop, on = [:Year, :AgeGrp])
    end

    @inbounds for (sn, mt) in fldict
        ofl[sn] = innerjoin(codf(mt, sn, opt = "s"), mean_pop, on = [:Year, :AgeGrp])
    end

    ofl["TX"] =
        innerjoin(codf(permutedims(tx), "TX", opt = "s"), mean_pop, on = [:Year, :AgeGrp])
    ost["Pop"] = tot_pop

    rocm = reshape(ocm, 9, 400)

    return (; ost = ost, ofl = ofl, adj_mu = adj_mort, rocm = rocm)
end

"""
Extract States
"""
@inline function munge_state(dtm)
    opd = Dict{String,Matrix{Float64}}()
    dl = length(axes(dtm, 1))
    if dl == 6
        opd["S"] = permutedims(dtm[1, :, :])
        opd["L"] = permutedims(dtm[2, :, :])
        opd["I"] = permutedims(dtm[3, :, :])
        opd["N"] = permutedims(dtm[4, :, :])
        opd["T"] = permutedims(dtm[5, :, :])
        opd["R"] = permutedims(dtm[6, :, :])
    elseif dl == 12
        opd["S"] = permutedims(dropdims(sum(dtm[[1, 7], :, :], dims = 1), dims = 1))
        opd["L"] = permutedims(dropdims(sum(dtm[[2, 8], :, :], dims = 1), dims = 1))
        opd["I"] = permutedims(dropdims(sum(dtm[[3, 9], :, :], dims = 1), dims = 1))
        opd["N"] = permutedims(dropdims(sum(dtm[[4, 10], :, :], dims = 1), dims = 1))
        opd["T"] = permutedims(dropdims(sum(dtm[[5, 11], :, :], dims = 1), dims = 1))
        opd["R"] = permutedims(dropdims(sum(dtm[[6, 12], :, :], dims = 1), dims = 1))
    end
    return opd
end

"""
Extract flows
"""
@inline function munge_flow(dtm)
    opd = Dict{String,Matrix{Float64}}()
    dl = length(axes(dtm, 1))
    if dl == 7
        opd["INC"] = permutedims(dtm[1, :, :])
        opd["MORT"] = permutedims(dtm[2, :, :])
        opd["IfSN"] = permutedims(dtm[3, :, :])
        opd["IfLN"] = permutedims(dtm[4, :, :])
        opd["IfLR"] = permutedims(dtm[5, :, :])
        opd["IfRN"] = permutedims(dtm[6, :, :])
        opd["IfRR"] = permutedims(dtm[7, :, :])
    elseif dl == 14
        opd["INC"] = permutedims(dropdims(sum(dtm[[1, 8], :, :], dims = 1), dims = 1))
        opd["MORT"] = permutedims(dropdims(sum(dtm[[2, 9], :, :], dims = 1), dims = 1))
        opd["IfSN"] = permutedims(dropdims(sum(dtm[[3, 10], :, :], dims = 1), dims = 1))
        opd["IfLN"] = permutedims(dropdims(sum(dtm[[4, 11], :, :], dims = 1), dims = 1))
        opd["IfLR"] = permutedims(dropdims(sum(dtm[[5, 12], :, :], dims = 1), dims = 1))
        opd["IfRN"] = permutedims(dropdims(sum(dtm[[6, 13], :, :], dims = 1), dims = 1))
        opd["IfRR"] = permutedims(dropdims(sum(dtm[[7, 14], :, :], dims = 1), dims = 1))
    end
    return opd
end

"""
Long DataFrame
"""
function ldf(dtm, vn)
    tpdf = DataFrame(dtm)
    transform!(tpdf, [:x1, :x2, :x3] => (+) => :x4)
    rename!(tpdf, ["[0,15)", "[15,65)", "[65,99]", "[0,99]"])
    tpdf.Year = 1900:2099
    return stack(tpdf, Not(:Year), variable_name = "AgeGrp", value_name = vn)
end

"""
Convert to dataframe and sum or mean
"""
function codf(mt, nm; opt = "m")::DataFrame
    df = DataFrame(mt, :auto)
    transform!(df, [:x2, :x3] => (+) => :x4)
    transform!(df, [:x1, :x2, :x3] => (+) => :x5)
    df.Year = repeat(1900:2099, inner = 2)
    gdf = groupby(df, :Year)

    if opt == "m"
        opdf = combine(gdf, [:x1, :x2, :x3, :x4, :x5] .=> mean)
    elseif opt == "s"
        opdf = combine(gdf, [:x1, :x2, :x3, :x4, :x5] .=> sum)
    elseif opt == "pi"
        opdf = df[collect(1:2:400), :]
        opdf = opdf[!, [:Year, :x1, :x2, :x3, :x4, :x5]]
    end

    rename!(opdf, ["Year", "[0,15)", "[15,65)", "[65,99]", "[15,99]", "[0,99]"])
    rv = stack(opdf, Not(:Year), variable_name = "AgeGrp", value_name = nm)
    return rv
end

"""
CSV writer
"""
function wo_s(
    dt,
    PSID,
    jn = "testing",
    cma = "M2",
    profile = (;
        vxtype = 0,
        effI = 0.0,
        effD = 0.0,
        dur = 0,
        ageG = [0],
        cov = 0.0,
        init_step = 0,
    ),
)

    vi_h = PSID

    @inbounds for (k, v) in dt[1]
        dt[1][k][!, :ph] .= vi_h
        opp = joinpath(
            "output",
            cma,
            "raw_output",
            "traj",
            jn,
            "state",
            pgen(; PSID = PSID, FN = k, profile = profile),
        )
        mkpath(dirname(opp))
        Arrow.write(opp, dt[1][k])
    end

    @inbounds for (k, v) in dt[2]
        dt[2][k][!, :ph] .= vi_h
        opp = joinpath(
            "output",
            cma,
            "raw_output",
            "traj",
            jn,
            "flow",
            pgen(; PSID = PSID, FN = k, profile = profile),
        )
        mkpath(dirname(opp))
        Arrow.write(opp, dt[2][k])
    end

    mort_path = mkpath(joinpath("output", cma, "params", "adj_mort", jn))
    CSV.write(
        joinpath(mort_path, "Adj_mort_$(vi_h).csv"),
        Tables.table(permutedims(dt[3]));
        writeheader = false,
    )

    cm_path = joinpath(
        "output",
        cma,
        "raw_output",
        "cm",
        jn,
        pgen(; PSID = PSID, FN = "cm", profile = profile, ext = "csv"),
    )
    mkpath(dirname(cm_path))
    CSV.write(cm_path, Tables.table(dt[4]), writeheader = false)

    return nothing

end

"""
Path generator
"""
function pgen(;
    PSID,
    FN,
    profile = (;
        vxtype = 0,
        effI = 0.0,
        effD = 0.0,
        dur = 0,
        ageG = [0],
        cov = 0.0,
        init_step = 0,
    ),
    ext = "arrow",
)

    nms = String.(collect(keys(profile)))
    kvp = ["$(k)=$(v)" for (k, v) in zip(nms, values(profile))]
    return joinpath(FN, kvp..., "$(FN)_$(PSID).$(ext)")
end
