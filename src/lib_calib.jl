using YAML
export target_extract, target_opt, sobgen, pobgen, vobgen, mobgen, wam, pst, btc

#=
Get latest 2050 targets
=#
inc2050hm = YAML.load_file("output/HM/targets/inc2050/latest.yml")
const i2050_hi = inc2050hm["i2050_hi"]
const i2050_lo = inc2050hm["i2050_lo"]
const i2050_med = inc2050hm["i2050_med"]

"""
Target extractor: extract targets from output arrays
"""
@inline function target_extract(optpl)
    st = optpl[1]
    fl = optpl[2]

    # Years
    p2010 = @view st[:, :, [221, 222]]
    p2015 = @view st[:, :, [231, 232]]
    p2019 = @view st[:, :, [239, 240]]
    p2050 = @view st[:, :, [301, 302]]
    # Denominators

    pop_all_2010 = sum(p2010) / 2.0
    pop_all_2015 = sum(p2015) / 2.0
    pop_all_2019 = sum(p2019) / 2.0
    pop_all_2050 = sum(p2050) / 2.0

    pop_014_2019 = sum(p2019[:, 1, 1:2]) / 2.0
    # pop_1564_2019 = sum(p2019[:, 2, 1:2]) / 2.0
    pop_65p_2019 = sum(p2019[:, 3, 1:2]) / 2.0
    pop_1599_2019 = sum(p2019[:, 2:3, 1:2]) / 2.0

    prevalence_all_2015 = (sum(st[3, :, [231, 232]]) / 2.0) / pop_all_2015 * 1e5
    incidence_all_2010 = sum(fl[1, :, [221, 222]]) / pop_all_2010 * 1e5
    incidence_all_2019 = sum(fl[1, :, [239, 240]]) / pop_all_2019 * 1e5
    incidence_014_2019 = sum(fl[1, 1, [239, 240]]) / pop_014_2019 * 1e5
    incidence_1599_2019 = sum(fl[1, 2:3, [239, 240]]) / pop_1599_2019 * 1e5
    incidence_65p_2019 = sum(fl[1, 3, [239, 240]]) / pop_65p_2019 * 1e5
    incidence_all_2050 = sum(fl[1, :, [301, 302]]) / pop_all_2050 * 1e5
    # mort_all_2000 = sum(fl[2, :, [221, 222]]) / pop_all_2000 * 1e5
    mort_all_2019 = sum(fl[2, :, [239, 240]]) / pop_all_2019 * 1e5
    notif_all_2019 = sum(st[5, :, [239, 240]]) / pop_all_2019 * 1e5

    return [
        prevalence_all_2015,
        # mort_all_2010,
        mort_all_2019,
        incidence_all_2010,
        incidence_all_2019,
        incidence_014_2019,
        incidence_1599_2019,
        incidence_65p_2019,
        notif_all_2019,
    ]

end

"""
Target distance calculator
Generates the unsquared distance from target
"""
@inline function target_opt(tv, mid, lo, hi)
    ΔH = hi - mid
    ΔL = mid - lo
    tv >= mid ? ((tv - mid) / ΔH) : ((mid - tv) / ΔL)
end

"""
Base vector objective function
"""
@inline function vobgen(ov)
    hov = similar(ov)
    hov[1] = target_opt(ov[1], 253.0, 195.0, 312.0)
    # hov[2] = target_opt(ov[2], 67.0, 57.0, 79.0)
    hov[2] = target_opt(ov[2], 33.0, 30.0, 35.0)
    hov[3] = target_opt(ov[3], 247.0, 128.0, 405.0)
    hov[4] = target_opt(ov[4], 193.0, 132.0, 266.0)
    hov[5] = target_opt(ov[5], 90.1, 54.9, 125.6)
    hov[6] = target_opt(ov[6], 230.2, 140.1, 321.3)
    hov[7] = target_opt(ov[7], 276.7, 0.0, 622.3)
    hov[8] = target_opt(ov[8], 190.0, 152.0, 228.0)
    return hov
end

# vobgen = @vobgen_mac()

"""
Vector objective function with out of bounds penalty
"""
@inline function bounds_penalty!(hov)
    @inbounds @simd for i in eachindex(hov)
        if hov[i] > 1.0
            hov[i] *= 8
        end
    end
    return hov
end

@inline pobgen(ov) = (bounds_penalty! ∘ vobgen)(ov)

"""
Scalar objective function with bounds penalty
    Generates sum of pobgen
"""
@inline sobgen(ov) = sum(pobgen(ov))

"""
Scalar objective function with maximum value output
"""
@inline mobgen(ov) = (maximum ∘ vobgen)(ov)

"""
Dichotomous target check function
"""
@inline function btc(ov)
    opbv = BitVector(undef, 9)
    opbv[1] = ov[1] >= 195.0 && ov[1] <= 312.0
    opbv[2] = ov[2] >= 30.0 && ov[2] <= 35.0
    opbv[3] = ov[3] >= 149.0 && ov[3] <= 473.0
    opbv[4] = ov[4] >= 132.0 && ov[4] <= 266.0
    opbv[5] = ov[5] >= 54.9 && ov[5] <= 125.6
    opbv[6] = ov[6] >= 140.1 && ov[6] <= 321.3
    opbv[7] = ov[7] >= 0.0 && ov[7] <= 622.3
    opbv[8] = ov[8] >= 152.0 && ov[9] <= 228.0
    return opbv
end

# btc = @btc_mac()

"""
Wrapped for ABC-MCMC
"""
function wam(x)
    pobgen(target_extract(main(variable_input = pst(x), fixed_input = fi)))
end

"""
Parameter set transformer: transform uniform input into constrained model params
"""
@inline function pst(pv)::Vector{Float64}
    ov = similar(pv)
    # [1]
    # e
    ov[1] = 0.4 + (pv[1] * 0.4)

    # [2,3,4]
    # fc, fa, fe
    ov[2] = pv[2] * 0.15
    ov[3] = 0.25 + (pv[3] * 0.5)
    ov[4] = 0.19 + (ov[3] - 0.19) * pv[4]

    # [5, 6]
    # nc, na, ne
    ov[5] = 0.1 + (pv[5] * 0.15)
    ov[6] = 0.1 + (ov[5] - 0.1) * pv[6]

    # [7] tr
    ov[7] = 1 + pv[7] * 3

    # [8] w
    ov[8] = 0.007 + (pv[8] * 0.013)

    # [9,10,11]
    # pc, pa, pe
    ov[9] = 0.01 + (pv[9] * 0.05)
    ov[10] = 0.08 + (pv[10] * 0.12)
    ov[11] = ov[10] + (0.36 - ov[10]) * pv[11]

    # [12, 13, 14]
    # rc, ra, re
    ov[12] = 0.005 + (pv[12] * 0.01)
    ov[13] = 0.005 + (pv[13] * 0.01)
    ov[14] = ov[13] + (0.015 - ov[13]) * pv[14]

    # [15, 16, 17]
    # uIc, uIa, uIe - note order change
    ov[16] = pv[16] * 0.178
    ov[15] = ov[16] + (0.178 - ov[16]) * pv[15]
    ov[17] = ov[16] + (0.178 - ov[16]) * pv[17]

    # [18, 19, 20]
    # uNc, uNa, uNe - note order change
    ov[19] = pv[19] * 0.034
    ov[18] = ov[19] + (0.034 - ov[19]) * pv[18]
    ov[20] = ov[19] + (0.034 - ov[19]) * pv[20]

    # [21, 22, 23]
    # vc, va, ve
    ov[21] = 0.0001 + (pv[21] * 0.0002)
    ov[22] = 0.0001 + (pv[22] * 0.0002)
    ov[23] = ov[22] + (0.04 - ov[22]) * pv[23]

    # [24]
    # x
    ov[24] = 0.25 + (pv[24] * 0.17)

    # [25]
    # c2000
    ov[25] = (1.0 - maximum(ov[15:17]) - ov[5] - ov[8]) * pv[25]

    # @assert all(ov .>= 0.0) "$pv -> $ov"

    return ov
end
