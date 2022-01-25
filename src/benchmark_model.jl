# Benchmark
using Revise
using tbconvax, BenchmarkTools

@code_warntype main(fixed_input = fi, variable_input = pst(zeros(18)))

@benchmark main(fixed_input = fi, variable_input = pst(zeros(18)))

vxc =
    (; effI = 0.0, effD = 0.0, dur = 0, vxtype = 1, ageG = [1], cov = 0.0, init_step = 255)
vx_tog = true

@benchmark main_sim(
    vxt = vx_tog,
    vx = vxc,
    fixed_input = fi,
    variable_input = zeros(Float64, 18),
)
