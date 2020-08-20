using Interpolations

const TINY = 1e-30
const _FWHM = 2 * sqrt(2 * log(2.0))

"""
    abeles(q::Array{Float64, 1}, layers::Array{Float64, 2})

Performs the Abeles optical matrix calculation.

Parameters
----------
- `q::Array{Float64, 1}` : q-vector values.
- `layers::Array{Float64, 2}` : an Nx4 array, where N is the number of layers in the system, 1st item in a given row is the thickness, the 2nd the SLD, the 3rd the imaginary SLD, and the 4th the roughness with the layer above.

Returns
-------
- `::Array{Float64, 1}` : unsmeared reflectometry values for the given q-vectors. 
"""
function abeles(q::Array{Float64, 1}, layers::Array{Float64, 2})
    nlayers = size(layers, 1) - 2
    q = hcat(fill(q, size(layers, 1))...)
    npnts = size(q, 1)
    
    mi00 = ones(ComplexF64, npnts, nlayers + 1)
    
    sld = zeros(ComplexF64, 1, nlayers + 2)
    
    sld[1, 2:end] += ((layers[2:end, 2] .- layers[1, 2]) + 1im * (abs.(layers[2:end, 3]) .+ TINY)) * 1e-6
            
    kn = sqrt.(q .^ 2 ./ 4 .- 4 .* pi .* sld)::Array{ComplexF64, 2}
        
    rj = kn[:, 1:end-1] .- kn[:, 2:end]
    rj ./= (kn[:, 1:end-1] + kn[:, 2:end])
    rj .*= exp.(-2 .* kn[:, 1:end-1] .* kn[:, 2:end] .* transpose(layers[2:end, 4]) .^ 2)
    
    if nlayers > 0
        mi00[:, 2:end] = exp.(kn[:, 2:end-1] * 1im * abs.(layers[2:end-1, 1]))
    end
    mi11 = 1. ./ mi00
    mi10 = rj .* mi00
    mi01 = rj .* mi11    
    
    mrtot00 = mi00[:, 1]
    mrtot01 = mi01[:, 1]
    mrtot10 = mi10[:, 1]
    mrtot11 = mi11[:, 1]
    
    for idx in range(2, nlayers + 1, step=1)
        p0 = mrtot00 .* mi00[:, idx] .+ mrtot10 .* mi01[:, idx]
        p1 = mrtot00 .* mi10[:, idx] .+ mrtot10 .* mi11[:, idx]
        mrtot00 = p0
        mrtot10 = p1

        p0 = mrtot01 .* mi00[:, idx] .+ mrtot11 .* mi01[:, idx]
        p1 = mrtot01 .* mi10[:, idx] .+ mrtot11 .* mi11[:, idx]

        mrtot01 = p0
        mrtot11 = p1
    end
            
    r = (mrtot01 ./ mrtot00)
    
    reflectivity = r .* conj(r)
    reflectivity = real(reflectivity)
    return reflectivity
end

"""
    same_convolution(a::Array{Float64, 1}, b::Array{Float64, 1})

Performs a convolution of one-dimensional arrays equivalent to the [`np.convolve`](https://numpy.org/doc/stable/reference/generated/numpy.convolve.html), where the `mode` is `'same'`. 

Parameters
----------
- `a::Array{Float64, 1}` : first array to convolve.
- `b::Array{Float64, 1}` : second array to convolve.

Returns
-------
- `::Array{Float64, 1}` : discrete, linear convolution of `a` and `b`. 
"""
function same_convolution(a::Array{Float64, 1}, b::Array{Float64, 1})

    m = length(a)
    n = length(b)
    
    padleft = ceil(Int32, n/2) - 1
    padright = floor(Int32, n/2)
    
    a = append!(zeros(padleft), a)
    a = append!(a, zeros(padright))
    
    output = zeros(m)
    for i in 1:m
        for j in 1:n
            output[i] += a[i+j-1] * b[n-j+1]
        end
    end
    return output
end

"""
    gauss(x::Array{Float64, 1}, s::Float64) 

A Gaussian kernel for resolution smearing. 

Parameters
----------
- `x::Array{Float64, 1}` : the kernal positions.
- `s::Float64` : the width of the kernel.

Returns
-------
- `::Array{Float64, 1}`: probabilities for `x`. 
"""
function gauss(x::Array{Float64, 1}, s::Float64) 
    g = 1. / s / sqrt(2 * pi) * exp.(-0.5 * x .^ 2 / s / s)
    return g
end

"""
    constant_smearing(q::Array{Float64, 1}, w::Array{Float64, 2}, resolution::Float64=0.5, scale::Float64=1., bkg::Float64=0.)

Perform the reflectometry calculation with a constant convolutional smearing. 

Parameters
----------
- `q::Array{Float64, 1}` : q-vector values.
- `layers::Array{Float64, 2}` : an Nx4 array, where N is the number of layers in the system, 1st item in a given row is the thickness, the 2nd the SLD, the 3rd the imaginary SLD, and the 4th the roughness with the layer above.
- `resolution::Float64` : the percentage resolution (dq/q) to be used. Defaults to `5.`.
- `scale::Float64` : the multiplicative scale factor assocaited with the reflectometry profile. Defaults to `1.`
- `bkg::Float64` : the uniform background to add to the profile. Defaults to `0.`. 

Returns
-------
- `::Array{Float64, 1}` : smeared reflectometry values for the given q-vectors. 
"""
function constant_smearing(q::Array{Float64, 1}, w::Array{Float64, 2}, resolution::Float64=5., scale::Float64=1., bkg::Float64=0.)
    if resolution < 0.5
        return abeles(q, w)
    end

    resolution /= 100
    gaussnum = 51
    gaussgpoint = (gaussnum - 1) / 2
    

    lowq = minimum(q[:])
    highq = maximum(q[:])
    if lowq <= 0.
        lowq = 1e-6
    end

    start = log10(lowq) - 6 * resolution / _FWHM
    finish = log10(highq * (1 + 6 * resolution / _FWHM))
    interpnum = round(abs(1 * (abs(start - finish)) / (1.7 * resolution / _FWHM / gaussgpoint)))
    xtemp = range(start, finish, length=Int(interpnum))
    xlin = 10. .^ xtemp
    
    gauss_x = collect(range(-1.7 * resolution, 1.7 * resolution, length=gaussnum))
    gauss_y = gauss(gauss_x, resolution / _FWHM)

    rvals = abeles(xlin, w)
    smeared_rvals = same_convolution(rvals, gauss_y)
    smeared_rvals *= gauss_x[2] - gauss_x[1]        
 
    itp = LinearInterpolation(xlin, smeared_rvals)
    smeared_output = itp(q)
    smeared_output *= scale
    smeared_output .+= bkg
    return smeared_output
end