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

@testset "variablesandbounds" begin
    l1 = Jeff.Layer(Jeff.Parameter(0, false, nothing),
                    Jeff.Parameter(0.000, false, nothing),
                    Jeff.Parameter(0.000, false, nothing),
                    Jeff.Parameter(0.000, false, nothing))
    l2 = Jeff.Layer(Jeff.Parameter(50, true, Uniform(40, 150)),
                    Jeff.Parameter(3.47, true, Normal(3.5, 5)),
                    Jeff.Parameter(0.000, false, nothing),
                    Jeff.Parameter(2., true, truncated(Normal(2, 5), 1.5, 4)))
    l3 = Jeff.Layer(Jeff.Parameter(0, false, nothing),
                    Jeff.Parameter(1.5, true, Normal(2.07, 5)),
                    Jeff.Parameter(0.000, false, nothing),
                    Jeff.Parameter(3., false, nothing))
    scale = Jeff.Parameter(1.5, true, Normal(1, 0.1))
    bkg = Jeff.Parameter(1e-9, true, Uniform(1e-10, 1e-6))
    model = Jeff.Model(scale, bkg, [l1, l2, l3])
    expected_variables = [1.5, 1e-9, 50, 3.47, 2., 1.5]
    expected_bounds = transpose(hcat([[-Inf, Inf], [1e-10, 1e-6], [40, 150], [-Inf, Inf], [1.5, 4.], [-Inf, Inf]]...))
    variables, bounds = Jeff.variablesandbounds(model)
    @test all(isapprox.(variables, expected_variables))
    @test all(isapprox.(bounds, expected_bounds))
end

@testset "populate_model" begin
    l1 = Jeff.Layer(Jeff.Parameter(0, false, nothing),
                    Jeff.Parameter(0.000, false, nothing),
                    Jeff.Parameter(0.000, false, nothing),
                    Jeff.Parameter(0.000, false, nothing))
    l2 = Jeff.Layer(Jeff.Parameter(50, true, Uniform(40, 150)),
                    Jeff.Parameter(3.47, true, Normal(3.5, 5)),
                    Jeff.Parameter(0.000, false, nothing),
                    Jeff.Parameter(2., true, truncated(Normal(2, 5), 1.5, 4)))
    l3 = Jeff.Layer(Jeff.Parameter(0, false, nothing),
                    Jeff.Parameter(1.5, true, Normal(2.07, 5)),
                    Jeff.Parameter(0.000, false, nothing),
                    Jeff.Parameter(3., false, nothing))
    scale = Jeff.Parameter(1.5, true, Normal(1, 0.1))
    bkg = Jeff.Parameter(1e-9, true, Uniform(1e-10, 1e-6))
    model = Jeff.Model(scale, bkg, [l1, l2, l3])
    variables, bounds = Jeff.variablesandbounds(model)
    variables[2] = 1e-10
    variables[4] = 3.
    model = Jeff.populate_model(variables, model)
    @test isapprox.(model.bkg.value, 1e-10)
    @test isapprox.(model.layers[2].sld.value, 3.)
end
