using Distributions

"""
    loglikelihood(distributions, model::Array{Float64, 1})

Return the log-likelihood between the measured data and some model data.

### Parameters
- `distributions::Array{Distributions.MultivariateDistribution}` : the measured data distributions.
- `model::Array{Float64, 2}` : the model data as an Nx2 array of abscissa and ordinate.

### Returns
- `::Float64` : the log-likelihood.
"""
function loglikelihood(distributions, model::Array{Float64, 2})
    ll = 0
    for i = 1:size(distributions, 1)
        ll += logpdf(distributions[i], model[i, :])
    end
    return ll
end

"""
    negativeloglikelihood(distributions, model::Array{Float64, 1})

Return the negative log-likelihood between the measured data and some model data.

### Parameters
- `distributions::Array{Distributions.MultivariateDistribution}` : the measured data distributions.
- `model::Array{Float64, 1}` : the model data as an Nx2 array of abscissa and ordinate.

### Returns
- `::Float64` : the negative log-likelihood.
"""
function negativeloglikelihood(distributions, model::Array{Float64, 2})
    return -1 * loglikelihood(distributions, model)
end
