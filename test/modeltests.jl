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

function test_layer(layer::Jeff.Layer)
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

@testset "layer_struct" begin
    test_layer(LAYER)
end

@testset "model_struct" begin
    scale = Jeff.Parameter(1, false, nothing)
    bkg = Jeff.Parameter(1e-7, true, Uniform(1e-8, 1e-6))
    model = Jeff.Model(scale, bkg, [LAYER, LAYER])
    @test isapprox(model.scale.value, 1)
    @test isequal(model.scale.vary, false)
    @test isequal(model.scale.prior, nothing)
    @test isapprox(model.bkg.value, 1e-7)
    @test isequal(model.bkg.vary, true)
    @test isapprox(minimum(model.bkg.prior), 1e-8)
    @test isapprox(maximum(model.bkg.prior), 1e-6)
    test_layer(model.layers[1])
    test_layer(model.layers[2])
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