using Jeff, MCMCChains, Turing, Distributions, DelimitedFiles, DataFrames, Optim, TuringCallbacks, Turkie, Makie

layers = Array{Any}(undef, (4, 4))
layers[1, 1] = 0.0
layers[1, 2] = 2.07
layers[1, 3] = 0.0
layers[1, 4] = 0
layers[2, 1] = 38
layers[2, 2] = 3.47
layers[2, 3] = 0.
layers[2, 4] = 5.86
layers[3, 1] = 258.827
layers[3, 2] = 2.4
layers[3, 3] = 0.0
layers[3, 4] = 9.0
layers[4, 1] = 0
layers[4, 2] = 6.36
layers[4, 3] = 0.0
layers[4, 4] = 3.68

scale = 0.8797
bkg = 4.4265e-7

data = readdlm("/Users/andrewmccluskey/work/analysis/Jeff.jl/dev/c_PLP0011859_q.txt")

x = data[:, 1]
y = data[:, 2]
yerr = data[:, 3]

function log_likelihood(ll, scale, bkg, x, y, yerr)
  layers[2, 1] = ll
	model = Jeff.constant_smearing(x, layers, 5., scale, bkg)
  sigma2 = yerr.^2.0
  return -0.5 .* sum((y .- model).^2 ./ sigma2 .+ log.(sigma2))
end

my_loglike = (ll, scale, bkg) -> log_likelihood(ll, scale, bkg, x, y, yerr)

@model function my_model()
  #layers[2, 1] ~ Uniform(15, 50)
  #layers[2, 4] ~ Uniform(1., 15.)
  #layers[3, 1] ~ Uniform(200., 300.)
  #layers[3, 2] ~ Uniform(1.5, 3.5)
  #layers[3, 4] ~ Uniform(1., 15.)
  #layers[4, 4] ~ Uniform(1., 15.)
  ll ~ Uniform(15, 50)
  scale ~ Uniform(0.6, 1.2)
  bkg ~ Uniform(1e-9, 9e-6)
  Turing.@addlogprob! my_loglike(ll, scale, bkg)
end

num_adapts = 250
mod = my_model()

callback = TensorBoardCallback("tensorboard_logs/run")
alg = NUTS(num_adapts, 0.65)
chain = sample(mod, alg, 2500; callback=callback)

display(chain)
#plot(chain)
#savefig("distribtions.pdf")

df = DataFrame(chain)

#p0 = plot(x, y, yerror=yerr, label="data", yaxis=:log10, markerstrokecolor=:auto, ylim=(2e-7, 2))
#for val in eachrow(df[1:5:end, [:"bkg", :"layers[2,1]", :"layers[2,4]", :"layers[3,1]", :"layers[3,2]", :"layers[3,4]", :"layers[4,4]", :"scale"]])
#     layers[2, 1] = val[2]
#     layers[2, 4] = val[3]
#     layers[3, 1] = val[4]
#     layers[3, 2] = val[5]
#     layers[3, 4] = val[6]
#     layers[4, 4] = val[7]
#     plot!(p0, x, Jeff.constant_smearing(x, layers, 5., val[8], val[1]), label="", linealpha=0.01, linecolor=:black)
#end
#savefig("reflectometry.pdf")
