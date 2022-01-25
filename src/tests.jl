# Test suite
using Test, tbconvax, CSV, Tables

@testset "pst" begin
    @test all(pst(rand(25)) .>= 0.0)
end;

@testset verbose = true "main-zero" begin
    om = main(fixed_input = fi, variable_input = zeros(Float64, 25))
    @test all(om[2] .== 0.0)
    @test all(om[1][5, :, :] .== 0.0)
    @test all(om[1] .>= 0.0)
end;

@testset verbose = true "main-vi-max" begin
    om = main(fixed_input = fi, variable_input = pst(ones(Float64, 25)))
    @test all(om[1] .>= 0.0)
    @test all(om[2] .>= 0.0)
end;

@testset verbose = true "main-vi-min" begin
    om = main(fixed_input = fi, variable_input = pst(zeros(Float64, 25)))
    @test all(om[1] .>= 0.0)
    @test all(om[2] .>= 0.0)
end;

# @testset verbose = true "demography" begin
#     demog = CSV.File("data/proc_data/3ac_total.csv") |> Tables.matrix
#     omz = main(fixed_input=fi, variable_input=pst(zeros(Float64,25)))
#     omfz = main(fixed_input=fi, variable_input=zeros(Float64,25))
#     @test isapprox(vec(sum(omz[1][:, :, 399], dims=1)), demog[150, 2:4])
#     @test isapprox(vec(sum(omfz[1][:, :, 399], dims=1)), demog[150, 2:4])

#     # Zero TB mort with rand other
#     iv = pst(rand(18))
#     iv[11:14] .= 0.0
#     omrz = main(fixed_input=fi, variable_input=iv)
#     @test isapprox(vec(sum(omrz[1][:, :, 399], dims=1)), demog[150, 2:4], atol=1e-9)
# end;

@testset verbose = true "vxsetup" begin
    hld = Matrix{Float64}(undef, 6, 3)
    # Blank vaccine characteristics
    vxc = (;
        effI = 0.0,
        effD = 0.5,
        dur = 10,
        vxtype = 1,
        ageG = [1],
        cov = 1,
        init_step = 255,
    )
    vx_matrices = vxsetup(; vxc...)
    @test vx_matrices[1][vxc.init_step, :] == [1.0, 0.0, 0.0]
    @test vx_matrices[1][275, :] == [1.0, 0.0, 0.0]
    @test vx_matrices[1][295, :] == [1.0, 0.0, 0.0]
    @test vx_matrices[5][275, :] == [1.0, 1.0, 1.0]
    @test vx_matrices[5][295, :] == [1.0, 1.0, 1.0]
end;

@testset verbose = true "vxswap" begin

    hld = Matrix{Float64}(undef, 6, 3)
    # Blank vaccine characteristics
    vxc = (;
        effI = 0.0,
        effD = 0.5,
        dur = 10,
        vxtype = 1,
        ageG = [1],
        cov = 1,
        init_step = 255,
    )
    vx_matrices = vxsetup(; vxc...)
    st = zeros(Float64, 12, 3, 400)
    hldm = zeros(Float64, 6, 3)
    st[1:6, :, :] .= 1.0
    for step = 2:400
        vxswp!(
            st,
            vx_matrices[1],
            vx_matrices[2],
            vx_matrices[3],
            vx_matrices[5],
            hldm,
            step,
        )
    end
    @test st[[7, 8, 12], :, 255] == [1.0 0.0 0.0; 1.0 0.0 0.0; 1.0 0.0 0.0]
end;

@testset verbose = true "null_vx" begin
    vxc = (;
        effI = 0.0,
        effD = 0.0,
        dur = 0,
        vxtype = 1,
        ageG = [1],
        cov = 0.0,
        init_step = 255,
    )
    vx_tog = true
    vx_output = main_sim(
        vxt = vx_tog,
        vx = vxc,
        fixed_input = fi,
        variable_input = zeros(Float64, 25),
    )
    bl_output = main(fixed_input = fi, variable_input = zeros(Float64, 25))
    @test all(vx_output[2] .== 0.0)
    @test all(vx_output[1][5, :, :] .== 0.0)
    @test all(vx_output[1] .>= 0.0)
    @test all(bl_output[2] .== 0.0)
    @test all(bl_output[1][5, :, :] .== 0.0)
    @test all(bl_output[1] .>= 0.0)

    @test isapprox(bl_output[1], vx_output[1][1:6, :, :])
    @test isapprox(bl_output[2], vx_output[2][1:7, :, :])
    @test isapprox(sum(vx_output[1], dims = 1), sum(bl_output[1], dims = 1))

end;

@testset verbose = true "active_vx" begin
    vxc = (;
        effI = 0.5,
        effD = 0.5,
        dur = 10,
        vxtype = 1,
        ageG = [4],
        cov = 0.7,
        init_step = 255,
    )
    vx_tog = true
    vx_output = main_sim(
        vxt = vx_tog,
        vx = vxc,
        fixed_input = fi,
        variable_input = zeros(Float64, 25),
    )
    bl_output = main(fixed_input = fi, variable_input = zeros(Float64, 25))
    @test isapprox(sum(vx_output[1], dims = 1), sum(bl_output[1], dims = 1))
end;
