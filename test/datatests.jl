using Jeff, Test, Measurements, Distributions, LinearAlgebra

X = range(0.01, 0.3, length=100)
Y =range(1, 1e-6, length=100)
DY = Y .* 0.001
DX = X .* 0.4

@testset "data_struct" begin
    q = zeros(Measurement, 100)
    R = zeros(Measurement, 100)
    for i in range(1, 100, step=1)
        q[i] = measurement(X[i], DX[i])
        R[i] = measurement(Y[i], DY[i])
    end
    distribution = MvNormal(Measurements.value.(R), Diagonal(Measurements.uncertainty.(R)))
    data = Jeff.Data(q, R, distribution, "filename.dat")
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), DX, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
    @test all(isapprox.(mean(data.distribution), Y, atol=1e-9)) 
    @test isequal(data.filepath, "filename.dat")
end

@testset "data_four_col" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "four_col.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), DX, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
    @test all(isapprox.(mean(data.distribution), Y, atol=1e-9)) 
    @test isequal(data.filepath, string(pwd(), Base.Filesystem.path_separator, "four_col.dat"))
end

@testset "data_three_col" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "three_col.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), X .* 0.05, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
    @test all(isapprox.(mean(data.distribution), Y, atol=1e-9)) 
    @test isequal(data.filepath, string(pwd(), Base.Filesystem.path_separator, "three_col.dat"))
end

@testset "data_two_col" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "two_col.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), X .* 0.05, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), Y .* 0.1, atol=1e-9))
    @test all(isapprox.(mean(data.distribution), Y, atol=1e-9)) 
    @test isequal(data.filepath, string(pwd(), Base.Filesystem.path_separator, "two_col.dat"))
end

@testset "data_one_col" begin
    @test_throws ArgumentError Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "one_col.dat"))
end