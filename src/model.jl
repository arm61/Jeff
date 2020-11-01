"""
    Parameter(value::Float64, vary::Bool, prior::Distrbutions.UnivariateDistribution)

A parameter that can be optimised in the analysis procress.

### Parameters
- `value::Float64` or `Distributions.UnivariateDistribution` : the value for the given parameter.
"""
mutable struct Parameter
    value
end

"""
    make_layer
"""
function make_layer(thick, sld, isld, rough)
    return Layer(Parameter(thick), Parameter(sld), Parameter(isld), Parameter(rough))
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
    layers_to_array(layers::Array{Jeff.Layer})

Convert from an array of N Jeff.Layer objects to an Nx4 array.

### Parameters
- `layers::Array{Jeff.Layer}`: layers to be converted.

### Returns
- `::Array{Any, 2}` : an array describing the layers.
"""
function layers_to_array(layers::Array{Jeff.Layer})
    nlayers = size(layers, 1)
    array = Array{Any}(undef, nlayers, 4)
    for i in 1:nlayers
        array[i, 1] = layers[i].thick.value
        array[i, 2] = layers[i].sld.value
        array[i, 3] = layers[i].isld.value
        array[i, 4] = layers[i].rough.value
    end
    return array
end
