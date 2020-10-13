using Jeff, Test, DelimitedFiles, Distributions

LAYERS = Array{Any}(undef, (3, 4))
LAYERS[1, 1] = 0.0
LAYERS[1, 2] = 0.0
LAYERS[1, 3] = 0.0
LAYERS[1, 4] = 0.0
LAYERS[2, 1] = 100
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

    @test all(isapprox.(Jeff.abeles(q_test, LAYERS), r_test, atol=1e-9))
end

@testset "constant_smearing" begin
    data = readdlm(string(pwd(), Base.Filesystem.path_separator, "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]

    @test all(isapprox.(Jeff.constant_smearing(q_test, LAYERS, resolution=5.), r_test, atol=1e-1))
end

@testset "arbitrary_smearing" begin
    data = readdlm(string(pwd(), Base.Filesystem.path_separator, "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]
    res = []
    for i = 1:size(q_test, 1)
        push!(res, Normal(q_test[i], q_test[i] * 0.02))
    end

    @test all(isapprox.(Jeff.arbitrary_smearing(q_test, LAYERS, res), r_test, atol=1e-1))
end
