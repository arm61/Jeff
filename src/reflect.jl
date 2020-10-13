using Interpolations, Distributions, ForwardDiff

const TINY = 1e-30
const _FWHM = 2 * sqrt(2 * log(2.0))

function Base.nextfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N}
     ForwardDiff.Dual{T}(nextfloat(d.value), d.partials)
end

function lssqs(x::T, y::T) where T<:Real
    k::Int = 0
    ρ = x*x + y*y
    if !isfinite(ρ) && (isinf(x) || isinf(y))
        ρ = convert(T, Inf)
    elseif isinf(ρ) || (ρ==0 && (x!=0 || y!=0)) || ρ<nextfloat(zero(T))/(2*eps(T)^2)
        m::T = max(abs(x), abs(y))
        k = m==0 ? m : exponent(m)
        xk, yk = ldexp(x,-k), ldexp(y,-k)
        ρ = xk*xk + yk*yk
    end
    ρ, k
end

lcopysign(x::Float64, y::Real) = copysign(x, Float64(y.value))

function lsqrt(z::Complex)
    z = float(z)
    x, y = reim(z)
    if x==y==0
        return Complex(zero(x),y)
    end
    ρ, k::Int = lssqs(x, y)
    if isa(x, AbstractFloat)
        if isfinite(x) ρ=ldexp(abs(x),-k)+sqrt(ρ) end
    else
        if isfinite(x) ρ=ldexp(abs(x).value,-k)+sqrt(ρ) end
    end
    if isodd(k)
        k = div(k-1,2)
    else
        k = div(k,2)-1
        ρ += ρ
    end
    if isa(ρ, AbstractFloat)
        ρ = ldexp(sqrt(ρ),k) #sqrt((abs(z)+abs(x))/2) without over/underflow
    else
        ρ = ldexp(sqrt(ρ).value,k) #sqrt((abs(z)+abs(x))/2) without over/underflow
    end
    ξ = ρ
    η = y
    if ρ != 0
        if isfinite(η) η=(η/ρ)/2 end
        if x<0
            ξ = abs(η)
            if isa(y, AbstractFloat)
                η = copysign(ρ,y)
            else
                η = lcopysign(ρ,y)
            end
        end
    end
    Complex(ξ,η)
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
function abeles(q::Array{Float64, 1}, layers::Array{Any, 2})
    nlayers = size(layers, 2) - 2
    npnts = size(q, 1)
    reflectivity = Array{Any}(undef, npnts)

    sub = layers[end, 2] + (1im * layers[end, 3] + TINY)
    super = layers[1, 2] + 0im
    sld = Array{Any}(undef, nlayers + 2)
    thick = Array{Any}(undef, nlayers)
    rough = Array{Any}(undef, nlayers)

    for i = 2:nlayers
        rho = layers[i, 2] + (1im * layers[i, 3] + TINY)
        sld[i] = 4e-6 * pi * (rho - super)
        thick[i - 1] = 0 + 1im * layers[i, 1]
        rough[i - 1] = -2 * layers[i, 4] ^ 2
    end
    sld[1] = 0 #+ 0im
    sld[nlayers+1] = 4e-6 * pi * (sub - super)
    rough[nlayers] = -2 * layers[end, 4] ^ 2

    mr = Array{Any}(undef, (2, 2))
    mi = Array{Any}(undef, (2, 2))

    for j = 1:npnts
        q2 = q[j] ^ 2 / 4. + 0im
        kn = q[j] / 2. + 0im
        for i = 1:nlayers
            knn = lsqrt(q2 - sld[i + 1])
            rj = (kn - knn) / (kn + knn) * exp(kn * knn * rough[i])

            if i == 1
                mr[1, 1] = 1. + 0im
                mr[1, 2] = rj
                mr[2, 2] = 1. + 0im
                mr[2, 1] = rj
            else
                beta = exp(kn * thick[i - 1])
                mi[1, 1] = beta
                mi[1, 2] = rj * beta
                mi[2, 2] = (1. + 0im) / beta
                mi[2, 1] = rj * mi[2, 2]
                mr = mr * mi
            end
            kn = knn
        end

        r = mr[2, 1] / mr[1, 1]
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
    constant_smearing(q::Array{Float64, 1}, w::Array{Float64, 2}; resolution::Float64=0.5, scale::Float64=1., bkg::Float64=0.)

Perform the reflectometry calculation with a constant convolutional smearing.

### Parameters
- `q::Array{Float64, 1}` : q-vector values.
- `w::Array{Float64, 2}` : an Nx4 array, where N is the number of layers in the system, 1st item in a given row is the thickness, the 2nd the SLD, the 3rd the imaginary SLD, and the 4th the roughness with the layer above.
- `resolution::Float64` : the percentage resolution (dq/q) to be used. Defaults to `5.`.
- `scale::Float64` : the multiplicative scale factor assocaited with the reflectometry profile. Defaults to `1.`
- `bkg::Float64` : the uniform background to add to the profile. Defaults to `0.`.

### Returns
- `::Array{Float64, 1}` : smeared reflectometry values for the given q-vectors.
"""
function constant_smearing(q::Array{Float64, 1}, w::Array{Any, 2}; resolution::Float64=5., scale::Float64=1., bkg::Float64=0.)
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

"""
    arbitrary_smearing(q::Array{Float64, 1}, w::Array{Float64, 2}, resolution::Array{Distributions.UnivariateDistribution, 1}; nsamples::Int64=500, scale::Float64=1., bkg::Float64=0.)

Perform resolution smearing based on an arbitrary (q-dependent) resolution function.

### Parameters
- `q::Array{Float64, 1}` : q-vector values.
- `w::Array{Float64, 2}` : an Nx4 array, where N is the number of layers in the system, 1st item in a given row is the thickness, the 2nd the SLD, the 3rd the imaginary SLD, and the 4th the roughness with the layer above.
- `resolution::Array{Distributions.UnivariateDistribution, 1}` : the q-dependent resolution function as a `Distributions` object.
- `nsamples::Int64=500` : Number of Monte Carlo samples to use in the smearing. Defaults to `500`.
- `scale::Float64` : the multiplicative scale factor assocaited with the reflectometry profile. Defaults to `1.`
- `bkg::Float64` : the uniform background to add to the profile. Defaults to `0.`.

### Returns
- `::Array{Float64, 1}` : arbitrarily smeared reflectometry values for the given q-vectors.
"""
function arbitrary_smearing(q::Array{Float64, 1}, w::Array{Any, 2}, resolution; nsamples::Int64=500, scale::Any=1., bkg::Any=0.)
    smeared_output = Array{Any}(undef, size(q, 1))
    for i = 1:size(q, 1)
        smeared_output[i] = mean(abeles(rand(resolution[i], nsamples), w))
        smeared_output[i] *= scale
        smeared_output[i] += bkg
    end
    return smeared_output
end
