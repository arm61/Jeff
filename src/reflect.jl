using Interpolations, Distributions, ForwardDiff

const TINY = 1e-30
const _FWHM = 2 * sqrt(2 * log(2.0))
const PI4 = 4e-6 * pi


Base.float(d::ForwardDiff.Dual{T}) where T = ForwardDiff.Dual{T}(float(d.value), d.partials)
Base.prevfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(prevfloat(float(d.value)), d.partials)
Base.nextfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(nextfloat(float(d.value)), d.partials)
function Base.ldexp(x::T, e::Integer) where T<:ForwardDiff.Dual
    if e >=0
        x * (1<<e)
    else
        x / (1<<-e)
    end
end

"""
    abeles(q::Array{Float64, 1}, layers::Array{Float64, 2})

Performs the Abeles optical matrix calculation.

### Parameters
- `q::Array{Float64, 1}` : q-vector values.
- `layers::Array{Float64, 2}` : an Nx4 array, where N is the number of layers in the system, 1st item in a given row is the thickness, the 2nd the SLD, the 3rd the imaginary SLD, and the 4th the roughness with the layer above.

### Returns
- `::Array{Float64, 1}` : unsmeared reflectometry values for the given q-vectors.
"""
function abeles(q::Array{Float64, 1}, layers)
    nlayers = size(layers, 1) - 2
    npnts = length(q)
    reflectivity = Array{Any}(undef, npnts)

    sld = Array{Any}(undef, nlayers + 2)
    thick = layers[:, 1] * 1im
    rough = Array{Any}(undef, nlayers + 1)

    for i = 2:nlayers+2
        sld[i] = PI4 * (layers[i, 2] - layers[1, 2] + 1im * (
            abs(layers[i, 3]) + TINY))
        rough[i - 1] = -2 * layers[i, 4] ^ 2
    end

    Threads.@threads for j = 1:npnts
        mr = Array{Any}(undef, (2, 2))
        mi = Array{Any}(undef, (2, 2))
        
        q2 = q[j] ^ 2. / 4. + 0im
        kn = q[j] / 2. + 0im

        for i = 1:nlayers + 1
            knn = sqrt(q2 - sld[i + 1])
            rj = (kn - knn) / (kn + knn) * exp(kn * knn * rough[i])

            if i == 1
                mr[1, 1] = 1. + 0im
                mr[1, 2] = rj
                mr[2, 2] = 1. + 0im
                mr[2, 1] = rj
            else
                beta = exp(kn * thick[i])
                mi[1, 1] = beta
                mi[2, 2] = (1. + 0im) / beta
                mi[2, 1] = rj * mi[1, 1]
                mi[1, 2] = rj * mi[2, 2]
                mr = mi * mr
            end
            kn = knn
        end

        r = mr[1, 2] / mr[1, 1]
        reflectivity[j] = real(r * conj(r))
    end
    return reflectivity
end

"""
    same_convolution(a::Array{Float64, 1}, b::Array{Float64, 1})

Performs a convolution of one-dimensional arrays equivalent to the [`np.convolve`](https://numpy.org/doc/stable/reference/generated/numpy.convolve.html), where the `mode` is `'same'`.

### Parameters
- `a::Array{Float64, 1}` : first array to convolve.
- `b::Array{Float64, 1}` : second array to convolve.

### Returns
- `::Array{Float64, 1}` : discrete, linear convolution of `a` and `b`.
"""
function same_convolution(a::Array{Any, 1}, b::Array{Float64, 1})

    m = length(a)
    n = length(b)

    padleft = ceil(Int32, n/2) - 1
    padright = floor(Int32, n/2)

    zero_padleft = Array{Any}(undef, padleft)
    for i in 1:padleft
        zero_padleft[i] = 0.
    end
    zero_padright = Array{Any}(undef, padleft)
    for i in 1:padleft
        zero_padright[i] = 0.
    end
    a = append!(zero_padleft, a)
    a = append!(a, zero_padright)

    output = Array{Any}(undef, m)
    for i in 1:m
        output[i] = 0.
    end
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

### Parameters
- `x::Array{Float64, 1}` : the kernal positions.
- `s::Float64` : the width of the kernel.

### Returns
- `::Array{Float64, 1}`: probabilities for `x`.
"""
function gauss(x::Array{Float64, 1}, s::Float64)
    g = zeros(size(x, 1))
    for i = 1:size(x, 1)
        g[i] = 1. / s / sqrt(2 * pi) * exp(-0.5 * x[i] ^ 2 / s / s)
    end
    return g
end

"""
    pointwise_smearing(q::Array{Float64, 1}, w::Array{Any, 2}, dq::Array{Float64, 1}; quad_order::Int32=17)

Perform the reflectometry calculation with a pointwise smearing.

### Parmaeters
- `q::Array{Float64, 1}` : q-vector values.
- `w::Array{Any, 2}` : an Nx4 array, where N is the number of layers in the system, 1st item in a given row is the thickness, the 2nd the SLD, the 3rd the imaginary SLD, and the 4th the roughness with the layer above.
- `dq::Array{Float64, 1}` : `dq` values for each `q` value.
- `quad_order{Int32, 1}` : The order of the Gaussian quadrature integration for the resolution integration.

### Returns
- `::Array{Float64, 1}` : smeared reflectometry values for the given q-vectors, with a pointise smearing.
"""
function pointwise_smearing(q::Array{Float64, 1}, w::Array{Any, 2}, dq::Array{Float64, 1}; quad_order::Int32=17)
end

"""
    constant_smearing(q::Array{Float64, 1}, w::Array{Any, 2}, resolution::Any, scale::Any, bkg::Any)

Perform the reflectometry calculation with a constant convolutional smearing.

### Parameters
- `q::Array{Float64, 1}` : q-vector values.
- `w::Array{Any, 2}` : an Nx4 array, where N is the number of layers in the system, 1st item in a given row is the thickness, the 2nd the SLD, the 3rd the imaginary SLD, and the 4th the roughness with the layer above.
- `resolution::Any` : the percentage resolution (dq/q) to be used.
- `scale::Any` : the multiplicative scale factor assocaited with the reflectometry profile.
- `bkg::Any` : the uniform background to add to the profile.

### Returns
- `::Array{Float64, 1}` : smeared reflectometry values for the given q-vectors.
"""
function constant_smearing(q::Array{Float64, 1}, w::Array{Any, 2}, resolution, scale, bkg)
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
    xlin = zeros(size(xtemp, 1))
    for i = 1:size(xtemp, 1)
        xlin[i] = 10. ^ xtemp[i]
    end

    gauss_x = collect(range(-1.7 * resolution, 1.7 * resolution, length=gaussnum))
    gauss_y = gauss(gauss_x, resolution / _FWHM)

    rvals = abeles(xlin, w)
    smeared_rvals = same_convolution(rvals, gauss_y)
    smeared_rvals *= gauss_x[2] - gauss_x[1]

    itp = LinearInterpolation(xlin, smeared_rvals)
    smeared_output = itp(q)
    for i = 1:size(smeared_output, 1)
        smeared_output[i] *= scale
        smeared_output[i] += bkg
    end
    return smeared_output
end
