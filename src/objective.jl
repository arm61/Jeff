using Turing, Measurements


"""
    t_model(q::Array{Float64, 1}, R::Array{Measurement, 1}, array::Array{Any, 2}, scale::Any, bkg::Any, resolution::Any)

The Turing model for the sampling of the reflectometry profile.

### Parameters
- `q::Array{Float}` : the q-vector values.
- `R::Array{Measurement}` : the reflected intensity values, including uncertainty.
- `array::Array{Any, 2}` : the layers from which the reflectometry should be found.
- `scale::Any` : the scale factor.
- `bkg::Any` : the background level.
- `resolution::Any` : the constant resolution width.
"""
@model t_model(q::Array{Float64, 1}, R::Array{Measurement, 1}, layers, scale, bkg, resolution) = begin
    layers_inp = Array{Any}(undef, size(layers))
    for i in 1:size(layers, 1)
        for j in 1:size(layers, 2)
            if isa(layers[i, j].value, AbstractFloat)
                layers_inp[i, j] = layers[i, j].value
            else
                layers_inp[i, j] ~ layers[i, j].value
            end
        end
    end
    if isa(scale.value, AbstractFloat)
        scale_inp = scale.value
    else
        scale_inp ~ scale.value
    end
    if isa(bkg.value, AbstractFloat)
        bkg_inp = bkg.value
    else
        bkg_inp ~ bkg.value
    end
    if isa(resolution.value, AbstractFloat)
        resolution_inp = resolution.value
    else
        resolution_inp ~ resolution.value
    end

    r = zeros(length(q))
    dr = zeros(length(q))
    for i in 1:length(q)
        r[i] = Measurements.value(R[i])
        dr[i] = Measurements.uncertainty(R[i])
    end
    mu = Jeff.constant_smearing(q, layers_inp, resolution_inp, scale_inp, bkg_inp)
    for i in 1:length(r)
        r[i] ~ Normal(mu[i], dr[i])
    end
end
