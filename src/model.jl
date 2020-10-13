"""
    Parameter(value::Float64, vary::Bool, prior::Distrbutions.UnivariateDistribution)

A parameter that can be optimised in the analysis procress.

### Parameters
- `value::Float64` : the value for the given parameter.
- `vary::Bool` : should the parameter be varied in optimisation.
- `prior::Distrbutions.UnivariateDistribution` : the prior probability distribution for the parameter. If `vary` is `false`, then `nothing` can be passed as the `prior`.
"""
mutable struct Parameter
    value
    name::String
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
    name::String
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
    name::String
end
