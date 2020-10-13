using Jeff, Distributions

SLD = Jeff.Parameter(Uniform(6.3, 6.4), "SLD")
ISLD = Jeff.Parameter(0.000, "iSLD")
THICK = Jeff.Parameter(50.0, "thick")
ROUGH = Jeff.Parameter(LogNormal(5, 1), "rough")
LAYER = Jeff.Layer(THICK, SLD, ISLD, ROUGH, "layer")

@testset "parameter_struct" begin
    @test isapprox(minimum(SLD.value), 6.3)
    @test isapprox(maximum(SLD.value), 6.4)
    @test isequal(SLD.name, "SLD")
end

@testset "parameter_struct_no_vary" begin
    @test isapprox(THICK.value, 50.0)
    @test isequal(THICK.name, "thick")
end

function test_layer(layer::Jeff.Layer)
    @test isapprox(layer.thick.value, 50.0)
    @test isapprox(minimum(layer.sld.value), 6.3)
    @test isapprox(maximum(layer.sld.value), 6.4)
    @test isapprox(layer.isld.value, 0.000)
    @test isapprox(meanlogx(layer.rough.value), 5.)
    @test isapprox(stdlogx(layer.rough.value), 1.)
    @test isequal(layer.name, "layer")
end

@testset "layer_struct" begin
    test_layer(LAYER)
end

@testset "model_struct" begin
    scale = Jeff.Parameter(1, "scale")
    bkg = Jeff.Parameter(Uniform(1e-8, 1e-6), "bkg")
    model = Jeff.Model(scale, bkg, [LAYER, LAYER], "model")
    @test isapprox(model.scale.value, 1)
    @test isapprox(minimum(model.bkg.value), 1e-8)
    @test isapprox(maximum(model.bkg.value), 1e-6)
    @test isequal(model.name, "model")
    test_layer(model.layers[1])
    test_layer(model.layers[2])
end
