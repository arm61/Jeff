using Jeff, Test, DelimitedFiles, Distributions

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
    data = readdlm(string(pwd(), Base.Filesystem.path_separator, "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]

    result = Jeff.abeles(q_test, LAYERS)

    @test all(isapprox.(result, r_test, atol=1e-9))
end

@testset "constant_smearing" begin
    data = readdlm(string(pwd(), Base.Filesystem.path_separator, "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]

    result = Jeff.constant_smearing(q_test, LAYERS, 0.0, 1.0, 0.0)

    @test all(isapprox.(result, r_test, atol=1e-9))
end
