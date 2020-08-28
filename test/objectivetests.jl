using Jeff, Test, Distributions, LinearAlgebra

@testset "loglikelihood" begin
    x = [0, 0]
    y = [0, 0]
    dx = [1, 1]
    dy = [1, 1]
    distributions = []
    for i = 1:size(x, 1)
        position = [x[i], y[i]]
        uncertainty = [dx[i], dy[i]]
        push!(distributions, MvNormal(position, Diagonal(uncertainty)))
    end
    data = zeros(Float64, (2, 2))
    @test isapprox(Jeff.loglikelihood(distributions, data), -3.675754132818691)
end

@testset "negativeloglikelihood" begin
    x = [0, 0]
    y = [0, 0]
    dx = [1, 1]
    dy = [1, 1]
    distributions = []
    for i = 1:size(x, 1)
        position = [x[i], y[i]]
        uncertainty = [dx[i], dy[i]]
        push!(distributions, MvNormal(position, Diagonal(uncertainty)))
    end
    data = zeros(Float64, (2, 2))
    @test isapprox(Jeff.negativeloglikelihood(distributions, data), 3.675754132818691)
end
