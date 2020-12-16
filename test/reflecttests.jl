using Jeff, Test, DelimitedFiles

curdir = pwd()
FWHM = 2.0 * sqrt(2.0 * log(2.0))


LAYERS = Array{Any}(undef, (3, 4))
LAYERS[1, 1] = 0.0
LAYERS[1, 2] = 0.0
LAYERS[1, 3] = 0.0
LAYERS[1, 4] = 0.0
LAYERS[2, 1] = 100.
LAYERS[2, 2] = 3.47
LAYERS[2, 3] = 0.
LAYERS[2, 4] = 2.
LAYERS[3, 1] = 0.0
LAYERS[3, 2] = 2.07
LAYERS[3, 3] = 0.0
LAYERS[3, 4] = 3.


@testset "abeles" begin
    data = readdlm(joinpath(curdir, "data", "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]
    result = Jeff.abeles(q_test, LAYERS)
    @test all(isapprox.(result, r_test, atol=1e-9))

    # some of the ORSO tests
    for i in [0, 1, 2, 3, 6]
        layers = readdlm(joinpath(curdir, "data", "test$i.layers"))
        data = readdlm(joinpath(curdir, "data", "test$i.dat"))
        q_test = data[:, 1]
        r_test = data[:, 2]
        result = Jeff.abeles(q_test, layers)
        # same tolerances as used in ORSO test for refnx
        @test all(isapprox.(result, r_test, rtol=8e-5))
    end
end


@testset "constant_smearing" begin
    data = readdlm(joinpath(curdir, "data", "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]

    result = Jeff.constant_smearing(q_test, LAYERS, 0.0, 1.0, 0.0)
    @test all(isapprox.(result, r_test, atol=1e-9))
end


@testset "pointwise_smearing" begin
    # ORSO tests
    data = readdlm(joinpath(curdir, "data", "test4.dat"))
    q_test = data[:, 1]
    r_test = data[:, 2]
    dq_test = data[:, 4] * FWHM
    layers = readdlm(joinpath(curdir, "data", "test0.layers"))
    result = Jeff.pointwise_smearing(q_test, layers, dq_test)
    # same tolerance as used in ORSO test for refnx
    # should look into whether this can be tightened.
    @test all(isapprox.(result, r_test, rtol=0.03))
end
