using Jeff, Distributions

@testset "parameter_struct" begin
    parameter = Jeff.Parameter(42.0, true, Normal(42.0, 5.0))
    @test isapprox(parameter.value, 42.0)
    @test isequal(parameter.vary, true)
    @test isapprox(mean(parameter.prior), 42.0)
    @test isapprox(std(parameter.prior), 5.0)
end

@testset "parameter_struct_no_vary" begin
    parameter = Jeff.Parameter(42.0, false, nothing)
    @test isapprox(parameter.value, 42.0)
    @test isequal(parameter.vary, false)
    @test isequal(parameter.prior, nothing)
end

@testset "layer_struct" begin
    sld = Jeff.Parameter(6.335, true, Uniform(6.3, 6.4))
    isld = Jeff.Parameter(0.000, false, nothing)
    thick = Jeff.Parameter(50.0, false, nothing)
    rough = Jeff.Parameter(5, true, LogNormal(5, 1))
    layer = Jeff.Layer(thick, sld, isld, rough)
    @test isapprox(layer.thick.value, 50.0)
    @test isequal(layer.thick.vary, false)
    @test isequal(layer.thick.prior, nothing)
    @test isapprox(layer.sld.value, 6.335)
    @test isequal(layer.sld.vary, true)
    @test isapprox(minimum(layer.sld.prior), 6.3)
    @test isapprox(maximum(layer.sld.prior), 6.4)
    @test isapprox(layer.isld.value, 0.000)
    @test isequal(layer.isld.vary, false)
    @test isequal(layer.isld.prior, nothing)
    @test isequal(layer.rough.vary, true)
    @test isapprox(meanlogx(layer.rough.prior), 5.)
    @test isapprox(stdlogx(layer.rough.prior), 1.)
end