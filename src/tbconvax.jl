module tbconvax
# Imports
using CSV,
    Tables,
    LinearAlgebra,
    Interpolations,
    SparseArrays,
    Debugger,
    DataFrames,
    Statistics,
    Infiltrator

export main, init_mat, fi, munge_flow, munge_state

"""
Initialise fixed inputs
"""
function init_mat(datasrc)

    # Demography
    # Master demography fila
    demog = CSV.File(joinpath(datasrc, "proc_data", "3ac_manual.csv")) |> Tables.matrix

    # Mort
    mort = transpose(demog[:, 6:8])

    # Births
    births = demog[:, 12]

    # Ageout
    ag = transpose(demog[:, 9:11])

    # Year zero
    yearzero = demog[1, 3:5]

    # Total
    total = demog[:, 3:5]

    return (;
        demog_total = total,
        demog_births = births,
        demog_mort = mort,
        demog_yearzero = yearzero,
        demog_ageout = ag,
    )

end
const fi = init_mat("data")

"""
Risk/rate scaler
"""
@inline function rs(x)
    1.0 - exp(log(1.0 - x) * 0.5)
end

"""
Population reset
"""
@inline function popreset(sm, yzo)
    op = sm ./ sum(sm, dims = 1) .* transpose(yzo)
    return op
end

"""
Treatment initiation rate generator
"""
@inline function tir_gen(scdr)
    # Treatment initiation rate
    sig = scdr ./ (1.0 .+ exp.(-0.25 .* (collect(1960:2020) .- 1990)))
    tir = interpolate(
        ([1900; 1959; collect(1960:2020); 2021; 2100],),
        [0; 0; sig; scdr; scdr],
        Gridded(Linear()),
    )
    return tir
end

"""
Treatment success rate generator
"""
macro tsr_gen()
    tsr = interpolate(
        ([1900, 1960, 2000, 2020, 2100],),
        [0.0, 0.0, 0.35, 0.85, 0.85],
        Gridded(Linear()),
    )
    return tsr
end

"""
Population fraction by age group
Returns 1 x n matrix of age-proportions in population
"""
@inline function afrac(sm)
    sum(sm, dims = 1) ./ sum(sm)
end

"""
Scale assortativity matrix
Argument: 1xn row vector of age-proportions
Value: n x n daily contact matrix
"""
@inline function scm(afrac)

    asm = [
        48.4 9.4 3.3
        9.4 17.1 6.0
        3.3 6.0 11.6
    ]

    return asm .* afrac
end

# M1 pair-wise correction
function pwc(pop)
    tpop = vec(sum(pop, dims = 1))
    rcm = [
        7.63717251606525 6.2944241384838 0.565957647241523
        1.480300832048 11.4901621621622 1.02108069401742
        0.521315766320647 3.99929106156634 1.99090909090909
    ]
    tcon = rcm .* tpop
    return ((tcon + permutedims(tcon)) / 2) ./ tpop
end

"""
Main iterator, using 3 age class matrix method
"""
@inline function main(; fixed_input, variable_input)
    # Set up conditions
    start_yr = 1900
    end_yr = 2100
    steps_per_year = 2
    dt = 0.5
    steps = 400
    init_dist = [0.39 0.5 0.02 0.02 0.0 0.07]
    # init_dist = [1.0 0.0 0.0 0.0 0.0 0.0]

    # Assign invariant matrix variables and rate scale
    e::Float64 = variable_input[1]
    f::Vector{Float64} = variable_input[[2, 3, 4]]
    n::Vector{Float64} = rs.(variable_input[[5, 5, 6]])
    tr::Float64 = variable_input[7]
    w::Float64 = rs(variable_input[8])
    p::Vector{Float64} = rs.(variable_input[[9, 10, 11]])
    r::Vector{Float64} = rs.(variable_input[[12, 13, 14]])
    uI::Vector{Float64} = rs.(variable_input[[15, 16, 17]])
    uN::Vector{Float64} = rs.(variable_input[[18, 19, 20]])
    v::Vector{Float64} = rs.(variable_input[[21, 22, 23]])
    x::Float64 = variable_input[24]
    c2020::Float64 = variable_input[25]
    uT = 0.0

    tir = tir_gen(c2020)
    κI = repeat(tir(1900:2100), inner = 2)
    κN = repeat(e .* tir(1900:2100), inner = 2)

    tsr = @tsr_gen
    ψ = repeat(tsr(1900:2100), inner = 2)
    ϕ = 1.0 .- ψ

    # Demography - rate scale and construct step-wise matrices
    ag = fixed_input.demog_ageout
    births = fixed_input.demog_births
    u_raw = fixed_input.demog_mort
    u = zeros(Float64, 3, 400)

    # Set up main state container
    state = zeros(Float64, 6, 3, 400)

    # views
    S = @view state[1, :, :]
    L = @view state[2, :, :]
    I = @view state[3, :, :]
    N = @view state[4, :, :]
    T = @view state[5, :, :]
    R = @view state[6, :, :]

    # Initialise state matrix
    state[:, :, 1] .= transpose(fixed_input.demog_yearzero .* init_dist)

    # Set up main flow container and views
    flow = zeros(Float64, 7, 3, 400)
    INC = @view flow[1, :, :]
    MORT = @view flow[2, :, :]
    IfRR = @view flow[3, :, :]
    IfRN = @view flow[4, :, :]
    IfLR = @view flow[5, :, :]
    IfLN = @view flow[6, :, :]
    IfSN = @view flow[7, :, :]

    # Set up container for contact matrices
    cmc = Array{Float64}(undef, 3, 3, 400)
    cmc[:, :, 400] .= 0.0

    # Main iterator
    @inbounds for i = 2:400

        cm = pwc(state[:, :, i-1])
        cmc[:, :, i-1] = cm
        # Age-wise infection prevalence
        inf_prev = I[:, i-1] ./ vec(sum(state[:, :, i-1], dims = 1))

        # S_{t+1} = S_{t} - λS_{t}
        # λ = 1 - exp(-βI)
        # I = inf_prev
        ttr = exp10(-tr)
        β = -180.0 * cm * log(1 - ttr)
        βI = β * inf_prev

        #TODO: Incorporate -kln(1-c) into β
        λ = @. 1 - exp(-βI)

        # λ = βI
        # I = lm
        # β = -k*ln(1-c)
        # λp = tr * lm
        # λp = 1.0 - exp(-tr * lm)
        # λ = fill(λp, 3)
        @assert all(λ .<= 1.0)

        # Mortality correction
        # @bp
        @. MORT[:, i-1] = (I[:, i-1] * uI) + (N[:, i-1] * uN)
        tpop = vec(sum(state[:, :, i-1], dims = 1))
        expected_deaths = tpop .* u_raw[:, i-1]
        actual_deaths = max.(0, (expected_deaths .- MORT[:, i-1]))
        u[:, i-1] = actual_deaths ./ tpop

        ## Model equations
        @inbounds for ac = 1:3
            S[ac, i] = S[ac, i-1] * (1.0 - λ[ac])

            L[ac, i] =
                (L[ac, i-1] * (1.0 - (λ[ac] * x * p[ac]) - v[ac])) +
                (S[ac, i-1] * λ[ac] * (1.0 - p[ac]))

            I[ac, i] =
                I[ac, i-1] * (1.0 - n[ac] - κI[i-1] - uI[ac]) +
                (S[ac, i-1] * λ[ac] * p[ac] * f[ac]) +
                (L[ac, i-1] * ((λ[ac] * x * p[ac]) + v[ac]) * f[ac]) +
                (N[ac, i-1] * w) +
                (R[ac, i-1] * ((λ[ac] * x * p[ac]) + r[ac]) * f[ac])

            N[ac, i] =
                N[ac, i-1] * (1.0 - n[ac] - w - κN[i-1] - uN[ac]) +
                (S[ac, i-1] * λ[ac] * p[ac] * (1.0 - f[ac])) +
                (L[ac, i-1] * ((λ[ac] * x * p[ac]) + v[ac]) * (1.0 - f[ac])) +
                (T[ac, i-1]) * (1.0 - uT) * ϕ[i-1] +
                (R[ac, i-1] * ((λ[ac] * x * p[ac]) + r[ac]) * (1.0 - f[ac]))

            T[ac, i] = (I[ac, i-1] * κI[i-1]) + (N[ac, i-1] * κN[i-1])

            R[ac, i] =
                (R[ac, i-1] * (1.0 - ((λ[ac] * x * p[ac]) + r[ac]))) +
                (I[ac, i-1] * n[ac]) +
                (N[ac, i-1] * n[ac]) +
                (T[ac, i-1] * (1.0 - uT) * ψ[i-1])

            INC[ac, i-1] =
                (S[ac, i-1] * λ[ac] * p[ac]) +
                (L[ac, i-1] * ((λ[ac] * x * p[ac]) + v[ac])) +
                (R[ac, i-1] * ((λ[ac] * x * p[ac]) + r[ac]))

            IfSN[ac, i-1] = (S[ac, i-1] * λ[ac] * p[ac])
            IfLN[ac, i-1] = L[ac, i-1] * (λ[ac] * x * p[ac])
            IfLR[ac, i-1] = L[ac, i-1] * v[ac]
            IfRN[ac, i-1] = R[ac, i-1] * (λ[ac] * x * p[ac])
            IfRR[ac, i-1] = R[ac, i-1] * r[ac]

        end


        # Demography
        @. S[:, i] += -(u[:, i-1] + ag[:, i-1]) * S[:, i-1]
        @. L[:, i] += -(u[:, i-1] + ag[:, i-1]) * L[:, i-1]
        @. I[:, i] += -(u[:, i-1] + ag[:, i-1]) * I[:, i-1]
        @. N[:, i] +=
            -(u[:, i-1] + ag[:, i-1]) * N[:, i-1] - (T[:, i-1] * u[:, i-1] * ϕ[i-1]) -
            ((ag[:, i-1]) * T[:, i-1] * ϕ[i-1])
        @. R[:, i] +=
            -(u[:, i-1] + ag[:, i-1]) * R[:, i-1] - (T[:, i-1] * u[:, i-1] * ψ[i-1]) -
            ((ag[:, i-1]) * T[:, i-1] * ψ[i-1])

        @. S[2:3, i] += (ag[1:2, i-1]) * S[1:2, i-1]
        @. L[2:3, i] += (ag[1:2, i-1]) * L[1:2, i-1]
        @. I[2:3, i] += (ag[1:2, i-1]) * I[1:2, i-1]
        @. N[2:3, i] +=
            (ag[1:2, i-1]) * N[1:2, i-1] + ((ag[1:2, i-1]) * T[1:2, i-1] * ϕ[i-1])
        @. R[2:3, i] +=
            (ag[1:2, i-1]) * R[1:2, i-1] + ((ag[1:2, i-1]) * T[1:2, i-1] * ψ[i-1])

        # Add births
        if i % 2 == 1
            S[1, i] += births[i]
        end

        # 1950 reset
        if i == 101
            state[:, :, i] = popreset(state[:, :, i], fixed_input.demog_yearzero)
        end

    end

    @assert all(state .>= 0.0) "Negative state"
    @assert all(flow .>= 0.0) "Negative flow"
    # @assert all(u .>= 0.0) "Negative background mortality: $(@show [uI, uN]) and $(u)"

    return state, flow, u, cmc

end

# Output processing functions
include("lib_proc.jl")

# Calibration related functions
include("lib_calib.jl")

# Internal vaccine model related functions
include("lib_vx.jl")

# Not otherwise specified
include("lib_nos.jl")

end # module
