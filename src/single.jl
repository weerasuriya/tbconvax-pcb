using Revise
using tbconvax, CSV, Tables, DataFrames, ProgressBars, Infiltrator, Debugger, YAML

# Job name for parameter sets
JN = YAML.load_file("data/params_latest.yml")["params_latest"]

# Fully processed (post-PST) fitted parameters
hm_fp = Tables.matrix(CSV.File("data/proc_data/fitted_params/$JN.csv"))

PS = 300

vx_profile =
    (; vxtype = 2, effI = 0.0, effD = 0.0, dur = 10, ageG = [1], cov = 0.7, init_step = 255)

params = hm_fp[PS, 1:25]
PSID = hm_fp[PS, end]

bl_output = standard(main(fixed_input = fi, variable_input = params))
ufi = merge(fi, (; demog_mort = bl_output[3]))

vx_output = standard(
    main_sim(vxt = true, vx = vx_profile, fixed_input = ufi, variable_input = params),
)
