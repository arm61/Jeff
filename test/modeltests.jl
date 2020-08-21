using Jeff, Distributions

SLD = Jeff.Parameter(6.335, true, Uniform(6.3, 6.4))
ISLD = Jeff.Parameter(0.000, false, nothing)
THICK = Jeff.Parameter(50.0, false, nothing)
ROUGH = Jeff.Parameter(5, true, LogNormal(5, 1))
LAYER = Jeff.Layer(THICK, SLD, ISLD, ROUGH)

@testset "parameter_struct" begin
    @test isapprox(SLD.value, 6.335)
    @test isequal(SLD.vary, true)
    @test isapprox(minimum(SLD.prior), 6.3)
    @test isapprox(maximum(SLD.prior), 6.4)
end

@testset "parameter_struct_no_vary" begin
    @test isapprox(THICK.value, 50.0)
    @test isequal(THICK.vary, false)
    @test isequal(THICK.prior, nothing)
end

@testset "layer_struct" begin
    @test isapprox(LAYER.thick.value, 50.0)
    @test isequal(LAYER.thick.vary, false)
    @test isequal(LAYER.thick.prior, nothing)
    @test isapprox(LAYER.sld.value, 6.335)
    @test isequal(LAYER.sld.vary, true)
    @test isapprox(minimum(LAYER.sld.prior), 6.3)
    @test isapprox(maximum(LAYER.sld.prior), 6.4)
    @test isapprox(LAYER.isld.value, 0.000)
    @test isequal(LAYER.isld.vary, false)
    @test isequal(LAYER.isld.prior, nothing)
    @test isequal(LAYER.rough.vary, true)
    @test isapprox(meanlogx(LAYER.rough.prior), 5.)
    @test isapprox(stdlogx(LAYER.rough.prior), 1.)
end

@testset "convert_to_array" begin
    layer2 = Jeff.Layer(Jeff.Parameter(20.0, false, nothing), SLD, ISLD, ROUGH)
    layers = [LAYER, layer2]
    array = Jeff.convert_to_array(layers)
    @test isapprox(array[1, 1], 50.0)
    @test isapprox(array[1, 2], 6.335)
    @test isapprox(array[1, 3], 0.000)
    @test isapprox(array[1, 4], 5)
    @test isapprox(array[2, 1], 20.0)
    @test isapprox(array[2, 2], 6.335)
    @test isapprox(array[2, 3], 0.000)
    @test isapprox(array[2, 4], 5)
end