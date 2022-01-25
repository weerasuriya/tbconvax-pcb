using Revise
using tbconvax, CSV, Tables, DataFrames, ProgressBars, Infiltrator, Debugger, YAML

# Job name for parameter sets
CMA = YAML.load_file("data/CMA_indicator.yml")["CMA"];
JN = YAML.load_file(joinpath("output", CMA, "metadata", "latest.yml"))["JN"];

# Fully processed (post-PST) fitted parameters
hm_fp = Tables.matrix(CSV.File("output/$(CMA)/params/sets/$JN.csv"));

@inbounds for PS in ProgressBar(axes(hm_fp, 1))
    params = hm_fp[PS, 1:25]
    PSID = hm_fp[PS, end]

    # Baseline simulation
    bl_output = standard(main(fixed_input = fi, variable_input = params))
    wo_s(bl_output, PSID, JN, CMA)

end

@inbounds for PS in ProgressBar(axes(hm_fp, 1))
    params = hm_fp[PS, 1:25]
    PSID = hm_fp[PS, end]

    # Update mortality
    adj_mort =
        CSV.File(
            joinpath("output", CMA, "params", "adj_mort", JN, "Adj_mort_$PSID.csv"),
            header = false,
        ) |>
        Tables.matrix |>
        permutedims
    ufi = merge(fi, (; demog_mort = adj_mort))

    # Vaccine simulations
    for vxt in 3,
        effI in [0.0],
        effD in [0.5],
        dur in 10,
        ageG in [[1], [2], [3], [4]],
        cov in [0.7]

        vx_profile = (;
            vxtype = vxt,
            effI = effI,
            effD = effD,
            dur = dur,
            ageG = ageG,
            cov = cov,
            init_step = 255,
        )
        vx_output = standard(
            main_sim(;
                vxt = true,
                vx = vx_profile,
                fixed_input = ufi,
                variable_input = params,
            ),
        )
        wo_s(vx_output, PSID, JN, CMA, vx_profile)
    end

end
