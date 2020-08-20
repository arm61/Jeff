using Jeff, Test, Measurements

X = range(0.01, 0.3, length=100)
Y = range(1, 1e-6, length=100)
DY = Y .* 0.001
DX = X .* 0.4

@testset "data_four_col" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "four_col.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), DX, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
end

@testset "data_three_col" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "three_col.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), X .* 0.05, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
end

@testset "data_two_col" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "two_col.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), X .* 0.05, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), Y .* 0.1, atol=1e-9))
end

@testset "data_one_col" begin
    @test_throws ArgumentError Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "one_col.dat"))
end