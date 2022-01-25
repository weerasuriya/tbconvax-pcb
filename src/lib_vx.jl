"""
Vaccine set up - set up vaccine introduction, campaigns, waning etc
Duration in years
"""

export vxsetup, vxswp!, main_sim

@inline function vxsetup(;
    effI,
    effD,
    dur,
    vxtype,
    ageG,
    cov,
    init_step = 255,
)::Vector{Matrix{Float64}}
    if ageG == [4]
        ageG = [1, 2, 3]
    end
    step_dur::Int64 = dur * 2

    vxnull = zeros(Float64, 400, 3)
    vxi = zeros(Float64, 400, 3)
    vxm = zeros(Float64, 400, 3)
    wm = zeros(Float64, 400, 3)

    vx_steps = [255, 275, 295]
    waning_steps = min.(400, (vx_steps .+ step_dur))

    vxm[vx_steps, ageG] .= cov
    wm[waning_steps, :] .= 1.0
    vxi[vx_steps, ageG] .= 1.0

    if vxtype == 1
        vxmS = vxmL = vxmR = vxm
    elseif vxtype == 2
        vxmS = vxm
        vxmL = vxmR = vxnull
    elseif vxtype == 3
        vxmL = vxmR = vxm
        vxmS = vxnull
    end

    return [vxmS, vxmL, vxmR, vxi, wm]

end

"""
Vaccine, non-vaccine matrix swapper
"""
@inline function vxswp!(st, vxmS, vxmL, vxmR, wm, hld, i)::Array{Float64}
    @inbounds for ag = 1:3
        # NVx -> Vx
        # S -> vS
        hld[1, ag] =
            st[1, ag, i] * (1.0 - vxmS[i, ag]) +
            ((1.0 - vxmS[i, ag]) * wm[i, ag] * st[7, ag, i])
        # L -> vL
        hld[2, ag] =
            st[2, ag, i] * (1.0 - vxmL[i, ag]) +
            ((1.0 - vxmL[i, ag]) * wm[i, ag] * st[8, ag, i])
        # R -> vR
        hld[3, ag] =
            st[6, ag, i] * (1.0 - vxmR[i, ag]) +
            ((1.0 - vxmR[i, ag]) * wm[i, ag] * st[12, ag, i])

        # Vx -> NVx
        # vS -> S
        hld[4, ag] =
            (st[7, ag, i] + (vxmS[i, ag] * st[1, ag, i])) -
            ((1.0 - vxmS[i, ag]) * wm[i, ag] * st[7, ag, i])
        # vL -> L
        hld[5, ag] =
            (st[8, ag, i] + (vxmL[i, ag] * st[2, ag, i])) -
            ((1.0 - vxmL[i, ag]) * wm[i, ag] * st[8, ag, i])
        # vR -> R
        hld[6, ag] =
            (st[12, ag, i] + (vxmR[i, ag] * st[6, ag, i])) -
            ((1.0 - vxmR[i, ag]) * wm[i, ag] * st[12, ag, i])

        st[[1, 2, 6], ag, i] = hld[[1, 2, 3], ag]
        st[[7, 8, 12], ag, i] = hld[[4, 5, 6], ag]
    end

    return st
end

"""
Main iterator, using 3 age class matrix method
"""
@inline function main_sim(; vxt, vx, fixed_input, variable_input)
    # Set up conditions
    start_yr = 1900
    end_yr = 2100
    steps_per_year = 2
    dt = 0.5
    steps = 400
    init_dist = [0.39 0.5 0.02 0.02 0.0 0.07 0.0 0.0 0.0 0.0 0.0 0.0]
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

    # Vaccine parameters
    pv = p .* (1.0 - vx.effD)
    rv = r .* (1.0 - vx.effD)
    vv = v .* (1.0 - vx.effD)

    tir = tir_gen(c2020)
    κI = repeat(tir(1900:2100), inner = 2)
    κN = repeat(e .* tir(1900:2100), inner = 2)

    tsr = @tsr_gen
    ψ = repeat(tsr(1900:2100), inner = 2)
    ϕ = 1.0 .- ψ

    # Demography - rate scale and construct step-wise matrices
    ag = fixed_input.demog_ageout
    births = fixed_input.demog_births
    u = fixed_input.demog_mort

    # Set up main state container
    state = zeros(Float64, 12, 3, 400)

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
    flow = zeros(Float64, 14, 3, 400)
    INC = @view flow[1, :, :]
    MORT = @view flow[2, :, :]
    IfRR = @view flow[3, :, :]
    IfRN = @view flow[4, :, :]
    IfLR = @view flow[5, :, :]
    IfLN = @view flow[6, :, :]
    IfSN = @view flow[7, :, :]

    # Vaccine views and matrices
    if vxt
        vS = @view state[7, :, :]
        vL = @view state[8, :, :]
        vI = @view state[9, :, :]
        vN = @view state[10, :, :]
        vT = @view state[11, :, :]
        vR = @view state[12, :, :]
        vINC = @view flow[8, :, :]
        vMORT = @view flow[9, :, :]
        vIfRR = @view flow[10, :, :]
        vIfRN = @view flow[11, :, :]
        vIfLR = @view flow[12, :, :]
        vIfLN = @view flow[13, :, :]
        vIfSN = @view flow[14, :, :]

        vxmat = vxsetup(; vx...)
        hldm = zeros(Float64, 12, 3)

    end

    # Set up container for contact matrices
    cmc = Array{Float64}(undef, 3, 3, 400)
    cmc[:, :, 400] .= 0.0

    # Main iterator
    @inbounds for i = 2:400

        # Transmission
        # Population age fraction
        paf = afrac(state[:, :, i-1])
        # Scale assortativity matrix to current population
        cm = scm(paf)
        cmc[:, :, i-1] = cm
        # Age-wise infection prevalence
        inf_prev = (I[:, i-1] + vI[:, i-1]) ./ vec(sum(state[:, :, i-1], dims = 1))

        # S_{t+1} = S_{t} - λS_{t}
        # λ = 1 - exp(-βI)
        # I = inf_prev
        ttr = exp10(-tr)
        β = -180.0 * cm * log(1 - ttr)
        βI = β * inf_prev

        #TODO: Incorporate -kln(1-c) into β
        λ = @. 1 - exp(-βI)

        vλ = λ * (1.0 - vx.effI)
        @assert all(λ .<= 1.0)

        # @. MORT[:, i-1] = (I[:, i-1] * uI) + (N[:, i-1] * uN)
        # if vxt && i >= vx.init_step
        #     @. vMORT[:, i-1] = (vI[:, i-1] * uI) + (vN[:, i-1] * uN)
        # end
        # combo_mort = MORT[:, i-1] + vMORT[:, i-1]
        # tpop = vec(sum(state[:, :, i-1], dims = 1))
        # expected_deaths = tpop .* u_raw[:, i-1]
        # actual_deaths = max.(0, (expected_deaths .- combo_mort))
        # u[:, i-1] = actual_deaths ./ tpop

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
            MORT[ac, i-1] = (I[ac, i-1] * uI[ac]) + (N[ac, i-1] * uN[ac])

            if vxt && i >= vx.init_step
                # Vaccine equations
                vS[ac, i] = vS[ac, i-1] * (1.0 - vλ[ac])

                vL[ac, i] =
                    (vL[ac, i-1] * (1.0 - (vλ[ac] * x * pv[ac]) - vv[ac])) +
                    (vS[ac, i-1] * vλ[ac] * (1.0 - pv[ac]))

                vI[ac, i] =
                    vI[ac, i-1] * (1.0 - n[ac] - κI[i-1] - uI[ac]) +
                    (vS[ac, i-1] * vλ[ac] * pv[ac] * f[ac]) +
                    (vL[ac, i-1] * ((vλ[ac] * x * pv[ac]) + vv[ac]) * f[ac]) +
                    (vN[ac, i-1] * w) +
                    (vR[ac, i-1] * ((vλ[ac] * x * pv[ac]) + rv[ac]) * f[ac])

                vN[ac, i] =
                    vN[ac, i-1] * (1.0 - n[ac] - w - κN[i-1] - uN[ac]) +
                    (vS[ac, i-1] * vλ[ac] * pv[ac] * (1.0 - f[ac])) +
                    (vL[ac, i-1] * ((vλ[ac] * x * pv[ac]) + vv[ac]) * (1.0 - f[ac])) +
                    (vT[ac, i-1]) * (1.0 - uT) * ϕ[i-1] +
                    (vR[ac, i-1] * ((vλ[ac] * x * pv[ac]) + rv[ac]) * (1.0 - f[ac]))

                vT[ac, i] = (vI[ac, i-1] * κI[i-1]) + (vN[ac, i-1] * κN[i-1])

                vR[ac, i] =
                    (vR[ac, i-1] * (1.0 - ((vλ[ac] * x * pv[ac]) + rv[ac]))) +
                    (vI[ac, i-1] * n[ac]) +
                    (vN[ac, i-1] * n[ac]) +
                    (vT[ac, i-1] * (1.0 - uT) * ψ[i-1])

                vINC[ac, i-1] =
                    (vS[ac, i-1] * vλ[ac] * pv[ac]) +
                    (vL[ac, i-1] * ((vλ[ac] * x * pv[ac]) + vv[ac])) +
                    (vR[ac, i-1] * ((vλ[ac] * x * pv[ac]) + rv[ac]))

                vIfSN[ac, i-1] = (vS[ac, i-1] * vλ[ac] * pv[ac])
                vIfLN[ac, i-1] = vL[ac, i-1] * (vλ[ac] * x * pv[ac])
                vIfLR[ac, i-1] = vL[ac, i-1] * vv[ac]
                vIfRN[ac, i-1] = vR[ac, i-1] * (vλ[ac] * x * pv[ac])
                vIfRR[ac, i-1] = vR[ac, i-1] * rv[ac]
                vMORT[ac, i-1] = (vI[ac, i-1] * uI[ac]) + (vN[ac, i-1] * uN[ac])
            end

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

        if vxt && i >= vx.init_step
            @. vS[:, i] += -(u[:, i-1] + ag[:, i-1]) * vS[:, i-1]
            @. vL[:, i] += -(u[:, i-1] + ag[:, i-1]) * vL[:, i-1]
            @. vI[:, i] += -(u[:, i-1] + ag[:, i-1]) * vI[:, i-1]
            @. vN[:, i] +=
                -(u[:, i-1] + ag[:, i-1]) * vN[:, i-1] - (vT[:, i-1] * u[:, i-1] * ϕ[i-1]) -
                ((ag[:, i-1]) * vT[:, i-1] * ϕ[i-1])
            @. vR[:, i] +=
                -(u[:, i-1] + ag[:, i-1]) * vR[:, i-1] - (vT[:, i-1] * u[:, i-1] * ψ[i-1]) -
                ((ag[:, i-1]) * vT[:, i-1] * ψ[i-1])

            @. vS[2:3, i] += (ag[1:2, i-1]) * vS[1:2, i-1]
            @. vL[2:3, i] += (ag[1:2, i-1]) * vL[1:2, i-1]
            @. vI[2:3, i] += (ag[1:2, i-1]) * vI[1:2, i-1]
            @. vN[2:3, i] +=
                (ag[1:2, i-1]) * vN[1:2, i-1] + ((ag[1:2, i-1]) * vT[1:2, i-1] * ϕ[i-1])
            @. vR[2:3, i] +=
                (ag[1:2, i-1]) * vR[1:2, i-1] + ((ag[1:2, i-1]) * vT[1:2, i-1] * ψ[i-1])
        end

        # Add births
        if i % 2 == 1
            S[1, i] += births[i]
        end

        # 1950 reset
        if i == 101
            state[:, :, i] = popreset(state[:, :, i], fixed_input.demog_yearzero)
        end

        if vxt && i >= vx.init_step
            vxswp!(state, vxmat[1], vxmat[2], vxmat[3], vxmat[5], hldm, i)
        end

    end
    # @bp
    @assert all(state .>= 0.0) "$variable_input"

    return state, flow, u, cmc

end
