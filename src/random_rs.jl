#= Pure rejection sampler =#
using Revise
using tbconvax
attempts = 10000

# Pre-allocate holding matrix
hm = Matrix{Float64}(undef, attempts, 18)

for i = 1:attempts
    try
        ov = rand(18)
        res = btc(target_extract(main(fixed_input = fi; variable_input = pst(ov))))
        if sum(res) >= 5
            @show res
            hm[i, :] = ov
        end
    catch
        nothing
    end
end
