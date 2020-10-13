using Jeff, Test, Measurements, Distributions, LinearAlgebra

X = range(0.01, 0.3, length=100)
Y = range(1, 1e-6, length=100)
DY = Y .* 0.001
DX = X .* 0.4
RES = []
for i = 1:100
    push!(RES, Normal(X[i], DX[i] * 0.5))
end

@testset "data_struct" begin
    q = zeros(Measurement, 100)
    R = zeros(Measurement, 100)
    res = []
    for i in range(1, 100, step=1)
        q[i] = measurement(X[i], DX[i])
        R[i] = measurement(Y[i], DY[i])
        push!(res, Normal(X[i], DX[i]))
    end
    data = Jeff.Data(q, R, res, "filename.dat")
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), DX, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
    @test all(isapprox.(mean.(data.resolution), X, atol=1e-9))
    @test isequal(data.name, "filename.dat")
end

@testset "data_four_col" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "four_col.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), DX, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
    @test isequal(data.name, string(pwd(), Base.Filesystem.path_separator, "four_col.dat"))
end

@testset "data_four_col_comma" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "four_col_comma.csv"); delim=',')
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), DX, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
    @test isequal(data.name, string(pwd(), Base.Filesystem.path_separator, "four_col_comma.csv"))
end

@testset "data_three_col" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "three_col.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), X .* 0.05, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
    @test isequal(data.name, string(pwd(), Base.Filesystem.path_separator, "three_col.dat"))
end

@testset "data_two_col" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "two_col.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), X .* 0.05, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), Y .* 0.1, atol=1e-9))
    @test isequal(data.name, string(pwd(), Base.Filesystem.path_separator, "two_col.dat"))
end

@testset "data_with_res" begin
    data = Jeff.read_data(string(pwd(), Base.Filesystem.path_separator, "with_res.dat"))
    @test all(isapprox.(Measurements.value.(data.q), X, atol=1e-9))
    @test all(isapprox.(Measurements.value.(data.R), Y, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.q), X .* 0.05, atol=1e-9))
    @test all(isapprox.(Measurements.uncertainty.(data.R), DY, atol=1e-9))
    @test isequal(data.name, string(pwd(), Base.Filesystem.path_separator, "with_res.dat"))
end

@testset "data_distribution" begin
    q = zeros(Measurement, 100)
    R = zeros(Measurement, 100)
    res = []
    for i in range(1, 100, step=1)
        q[i] = measurement(X[i], DX[i])
        R[i] = measurement(Y[i], DY[i])
        push!(res, Normal(X[i], DX[i]))
    end
    data = Jeff.Data(q, R, res, "test")
    @test all(isapprox.(hcat(mean.(Jeff.get_distribution(data.q, data.R))...)[1, :], X, atol=1e-9))
    @test all(isapprox.(hcat(mean.(Jeff.get_distribution(data.q, data.R))...)[2, :], Y, atol=1e-9))
end

@testset "data_transform_default" begin
    q = zeros(Measurement, 100)
    R = zeros(Measurement, 100)
    res = []
    for i in range(1, 100, step=1)
        q[i] = measurement(X[i], DX[i])
        R[i] = measurement(Y[i], DY[i])
        push!(res, Normal(X[i], DX[i]))
    end
    data = Jeff.Data(q, R, res, "test")
    @test all(isapprox.(Jeff.transform(data.R), log.(Y), atol=1e-9))
end

@testset "data_transform_rq4" begin
    q = zeros(Measurement, 100)
    R = zeros(Measurement, 100)
    res = []
    for i in range(1, 100, step=1)
        q[i] = measurement(X[i], DX[i])
        R[i] = measurement(Y[i], DY[i])
        push!(res, Normal(X[i], DX[i]))
    end
    data = Jeff.Data(q, R, res, "test")
    @test all(isapprox.(Jeff.transform(data.R, f=x->(log.(x) .* data.q .^ 4)), log.(Y) .* data.q .^ 4, atol=1e-9))
end
