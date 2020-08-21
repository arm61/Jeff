using Distributions

"""
    nll(data::Jeff.Data, model::Array{Float64, 1})

Return the negative log-likelihood between the measured data and some model data. 

### Parameters
- `data::Jeff.Data` : the measured data. 
- `model::Array{Float64, 1}` : the model data. 

### Returns 
- `::Float64` : the negative log-likelihood.
"""
function nll(data::Jeff.Data, model::Array{Float64, 1})
    return -1. * logpdf(data.distribution, model)
end