using Jeff, Test, DelimitedFiles

@testset "abeles" begin
    layers = zeros(Float64, (3, 4))
    layers[2, 1] = 100
    layers[2, 2] = 3.47
    layers[2, 4] = 2
    layers[3, 2] = 2.07
    layers[3, 4] = 3

    data = readdlm(string(pwd(), Base.Filesystem.path_separator, "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]

    @test all(isapprox.(Jeff.abeles(q_test, layers), r_test, atol=1e-9))
end

@testset "constant_smearing" begin
    layers = zeros(Float64, (3, 4))
    layers[2, 1] = 100
    layers[2, 2] = 3.47
    layers[2, 4] = 2
    layers[3, 2] = 2.07
    layers[3, 4] = 3

    data = readdlm(string(pwd(), Base.Filesystem.path_separator, "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]

    @test all(isapprox.(Jeff.constant_smearing(q_test, layers, 5.), r_test, atol=1e-1))
end