"""
    Parameter(value::Float64, vary::Bool, prior::Distrbutions.UnivariateDistribution)

A parameter that can be optimised in the analysis procress.

### Parameters
- `value::Float64` : the value for the given parameter.
- `vary::Bool` : should the parameter be varied in optimisation. 
- `prior::Distrbutions.UnivariateDistribution` : the prior probability distribution for the parameter. If `vary` is `false`, then `nothing` can be passed as the `prior`.
"""
mutable struct Parameter
    value::Float64
    vary::Bool
    prior
end

"""
    Layer(thick::Jeff.Parameter, sld::Jeff.Parameter, isld::Jeff.Parameter, rough::Jeff.Parameter)

A description of a layer in a system. 

### Parameters
- `thick::Jeff.Parameter` : the layer thickness.
- `sld::Jeff.Parameter` : layer real scattering length density.
- `isld::Jeff.Parameter` : layer imaginary scattering length density.
- `rough::Jeff.Parameter` : roughness with layer above.
"""
mutable struct Layer
    thick::Jeff.Parameter
    sld::Jeff.Parameter
    isld::Jeff.Parameter
    rough::Jeff.Parameter
end

"""
    Model(scale::Jeff.Parmeter, bkg::Jeff.Parameter, layers::Array{Jeff.Layer})

The model from which the reflectometry should be calculated. 

### Parameters
- `scale::Jeff.Parameter` : scale factor for reflectometry.
- `bkg::Jeff.Parameter` : uniform background. 
- `layers::Array{Jeff.Layer}` : the array of [`Jeff.Layer`] objects that describe the model.
"""
mutable struct Model
    scale::Jeff.Parameter
    bkg::Jeff.Parameter
    layers::Array{Jeff.Layer}
end

"""
    convert_to_array(layers::Array{Jeff.Layer})

Convert from a array of [`Jeff.Layer`](@ref) objects to an array that is compatible with [`reflect.jl`](@ref) functions.

### Parameters
- `layers::Array{Jeff.Layer}` : an array of [`Jeff.Layer`](@ref) objects.

### Returns
- `::Array{Float64, 2}` : a two-dimensional array that is compatible with the [`reflect.jl`](@ref) functions.
"""
function convert_to_array(layers::Array{Jeff.Layer})
    nlayers = size(layers, 1)
    array = zeros(Float64, (nlayers, 4))
    for i in range(1, nlayers, step=1)
        array[i, 1] = layers[i].thick.value
        array[i, 2] = layers[i].sld.value
        array[i, 3] = layers[i].isld.value
        array[i, 4] = layers[i].rough.value
    end
    return array
end