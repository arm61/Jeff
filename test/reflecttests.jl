using Jeff, Test, DelimitedFiles

@testset "abeles" begin
    layers = Array(zeros(Float64, 12))
    layers[2] = 0.
    layers[3] = 0.
    layers[5] = 100
    layers[6] = 3.47
    layers[8] = 2
    layers[10] = 2.07
    layers[12] = 3

    l = transpose(reshape(layers, (4, floor(Int64, size(layers, 1) / 4))))

    data = readdlm(string(pwd(), Base.Filesystem.path_separator, "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]

    @test all(isapprox.(Jeff.abeles(q_test, l), r_test, atol=1e-9))
end

@testset "constant_smearing" begin
    layers = Array(zeros(Float64, 12))
    layers[2] = 0.
    layers[3] = 0.
    layers[5] = 100
    layers[6] = 3.47
    layers[8] = 2
    layers[10] = 2.07
    layers[12] = 3

    l = transpose(reshape(layers, (4, floor(Int64, size(layers, 1) / 4))))

    data = readdlm(string(pwd(), Base.Filesystem.path_separator, "theoretical.txt"))
    q_test = data[:, 1]
    r_test = data[:, 2]

    @test all(isapprox.(Jeff.constant_smearing(q_test, l, 5), r_test, atol=1e-1))
end