"""
    nll(distribution::Distributions.MultivariateDistribution, model::Array{Float64, 1})

Return the negative log-likelihood between the measured data and some model data. 

### Parameters
- `distribution::Distributions.MultivariateDistribution` : the measured data. 
- `model::Array{Float64, 1}` : the model data. 

### Returns 
- `::Float64` : the negative log-likelihood.
"""
function nll(distribution::Distributions.MultivariateDistribution, model::Array{Float64, 1})
    return -1. * logpdf(distribution, model)
end