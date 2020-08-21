using Jeff, Test, Distributions

@testset "nll" begin
    distribution = MvNormal([0., 0.], Diagonal([1., 1.]))
    data = [0., 0.]
    @test isapprox(Jeff.nll(distribution, data), 1.8378770664093456)
end


