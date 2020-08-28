using Distributions

"""
    log_likelihood(distributions, model::Array{Float64, 1})

Return the log-likelihood between the measured data and some model data.

### Parameters
- `distributions::Array{Distributions.MultivariateDistribution}` : the measured data distributions.
- `model::Array{Float64, 2}` : the model data as an Nx2 array of abscissa and ordinate.

### Returns
- `::Float64` : the log-likelihood.
"""
function log_likelihood(distributions, model::Array{Float64, 2})
    ll = 0
    for i = 1:size(distributions, 1)
        ll += logpdf(distributions[i], model[i, :])
    end
    return ll
end

"""
    negative_log_likelihood(distributions, model::Array{Float64, 1})

Return the negative log-likelihood between the measured data and some model data.

### Parameters
- `distributions::Array{Distributions.MultivariateDistribution}` : the measured data distributions.
- `model::Array{Float64, 1}` : the model data as an Nx2 array of abscissa and ordinate.

### Returns
- `::Float64` : the negative log-likelihood.
"""
function negative_log_likelihood(distributions, model::Array{Float64, 2})
    return -1 * loglikelihood(distributions, model)
end


"""
    variables_and_bounds(model::Jeff.Model)

Get the variable values and bounds for a given model.

### Parameters
- `model::Jeff.Model` : reflectometry model.

### Returns
- `::Tuple{Array{Float64, 1}, Array{Float64, 2}}` : variable values and upper and lower bounds for each variable
"""
function variables_and_bounds(model::Jeff.Model)
    nlayers = size(model.layers, 1)
    variables = zeros(Float64, 0)
    bounds = []
    if model.scale.vary
        push!(variables, model.scale.value)
        push!(bounds, [minimum(model.scale.prior), maximum(model.scale.prior)])
    end
    if model.bkg.vary
        push!(variables, model.bkg.value)
        push!(bounds, [minimum(model.bkg.prior), maximum(model.bkg.prior)])
    end
    for i = 1:nlayers
        if model.layers[i].thick.vary
            push!(variables, model.layers[i].thick.value)
            push!(bounds, [minimum(model.layers[i].thick.prior), maximum(model.layers[i].thick.prior)])
        end
        if model.layers[i].sld.vary
            push!(variables, model.layers[i].sld.value)
            push!(bounds, [minimum(model.layers[i].sld.prior), maximum(model.layers[i].sld.prior)])
        end
        if model.layers[i].isld.vary
            push!(variables, model.layers[i].isld.value)
            push!(bounds, [minimum(model.layers[i].isld.prior), maximum(model.layers[i].isld.prior)])
        end
        if model.layers[i].rough.vary
            push!(variables, model.layers[i].rough.value)
            push!(bounds, [minimum(model.layers[i].rough.prior), maximum(model.layers[i].rough.prior)])
        end
    end
    return variables, transpose(hcat(bounds...))
end

"""
    populate_model(variables::Array{Float64, 1}, model::Jeff.Model)

Populate the model with a set of variables.

### Parameters
- `variables::Array{Float64, 1}` : variable to populate.
- `model::Jeff.Model` : model to be populated.

### Returns
- `::Jeff.Model` : populated model.
"""
function populate_model(variables::Array{Float64, 1}, model::Jeff.Model)
    j = 1
    if model.scale.vary
        model.scale.value = variables[j]
        j += 1
    end
    if model.bkg.vary
        model.bkg.value = variables[j]
        j += 1
    end
    nlayers = size(model.layers, 1)
    for i = 1:nlayers
        if model.layers[i].thick.vary
            model.layers[i].thick.value = variables[j]
            j += 1
        end
        if model.layers[i].sld.vary
            model.layers[i].sld.value = variables[j]
            j += 1
        end
        if model.layers[i].isld.vary
            model.layers[i].isld.value = variables[j]
            j += 1
        end
        if model.layers[i].rough.vary
            model.layers[i].rough.value = variables[j]
            j += 1
        end
    end
    return model
end
