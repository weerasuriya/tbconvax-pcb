# Single run - entry point
using Revise
using tbconvax, CSV, Tables, DataFrames, JuliaFormatter, BenchmarkTools

# ts = [
#     0.9993643361416588,
#     0.9991462696024104,
#     0.5401289935659376,
#     0.9991010242185427,
#     0.000416498715921727,
#     0.5087320854236573,
#     0.8595264456824183,
#     0.9985674720756037,
#     0.6279166512619861,
#     0.8157484724312818,
#     0.016630091361792197,
#     0.9950984943034717,
#     0.08208689283501618,
#     0.005190867826695006,
#     0.8778960189544066,
#     0.3796269274363202,
#     0.9854182346675122,
#     0.9901533584205806,
#     0.9990488232414719,
#     0.9935613162777513,
#     0.8271003235533392,
#     0.10297870559508938,
#     0.0002326836112379125,
#     0.7325121498654061,
#     0.22126306601275733,
# ]

ts = [0.9999999999999788, 0.9999999999999793, 0.9999, 0.9999999999999964, 1.0156082655814866e-7, 0.17144348827487524, 0.8783046517976585, 0.9999999999999944, 0.6001762059183694, 0.6288759892106566, 4.221965261338382e-5, 0.9999999999995639, 0.040851291556046554, 0.0002767046119145985, 0.08754525055688847, 0.5209273264612461, 0.9999999999997914, 0.9999726332513945, 0.9999999999999218, 0.9999640094460772, 0.3211357736544128, 0.030551856554940946, 3.519994928996584e-8, 0.9999999999999996, 0.2948059468527162]

# st, fl = main(fixed_input = fi, variable_input = pst(ts));

# wo_s(standard(main(fixed_input = fi, variable_input = pst(ts))), pst(ts))
target_extract(main(fixed_input = fi, variable_input = pst(ts)))
btc(target_extract(main(fixed_input = fi, variable_input = pst(ts))))