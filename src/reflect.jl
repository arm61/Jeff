using Interpolations, Distributions, ForwardDiff, LinearAlgebra, FastGaussQuadrature

const TINY = 1e-30
const _FWHM = 2 * sqrt(2 * log(2.0))
const PI4 = 4e-6 * pi
const _INTLIMIT = 3.5


# Base.float(d::ForwardDiff.Dual{T}) where T = ForwardDiff.Dual{T}(float(d.value), d.partials)
# Base.prevfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(prevfloat(float(d.value)), d.partials)
# Base.nextfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(nextfloat(float(d.value)), d.partials)


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
function abeles(q, w)
    nlayers = size(w, 1) - 2
    npnts = length(q)

    reflectivity = Vector{Any}(undef, (npnts))
    oneC = Complex(1.0)

    for j in eachindex(q)
        qq2 = (q[j] * q[j] / 4.0) + 0.0im
        kn = (q[j] / 2.) + 0.0im

        # variables are local to if blocks
        local MRtotal11, MRtotal12, MRtotal21, MRtotal22

        for i = 1:nlayers+1
            # wavevector in the layer
            sld_next = ((w[i+1, 2] - w[1, 2]) + ((abs(w[i+1, 3]) + TINY))im) * pi * 4.0e-6
            kn_next = sqrt(qq2 - sld_next)

            # reflectance of the interface
            rj = (kn - kn_next)/(kn + kn_next) * exp(kn * kn_next * (-2.0 * w[i+1, 4]^2))

            if i == 1
                # characteristic matrix for first interface
                MRtotal11 = oneC
                MRtotal12 = rj
                MRtotal21 = rj
                MRtotal22 = oneC
            else
                # work out the beta for the layer
                beta = exp(kn * (abs(w[i, 1]) * 1im))

                # this is the characteristic matrix of a layer
                MI11 = beta
                MI12 = rj * beta
                MI22 = oneC / beta
                MI21 = rj * MI22

                # propagate optical matrix by matmul
                p11 = MRtotal11 * MI11 + MRtotal12 * MI21
                p12 = MRtotal11 * MI12 + MRtotal12 * MI22
                p21 = MRtotal21 * MI11 + MRtotal22 * MI21
                p22 = MRtotal21 * MI12 + MRtotal22 * MI22

                MRtotal11 = p11
                MRtotal12 = p12
                MRtotal21 = p21
                MRtotal22 = p22

            end
            kn = kn_next;

        end
        reflectivity[j] = MRtotal21 / MRtotal11
    end
    return real(reflectivity .* conj(reflectivity))
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
function same_convolution(a, b)

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

### Parameters
- `q::Array{Float64, 1}` : q-vector values.
- `w::Array{Any, 2}` : an Nx4 array, where N is the number of layers in the system, 1st item in a given row is the thickness, the 2nd the SLD, the 3rd the imaginary SLD, and the 4th the roughness with the layer above.
- `dq::Array{Float64, 1}` : `dq` values for each `q` value.
- `quad_order{Int32, 1}` : The order of the Gaussian quadrature integration for the resolution integration.

### Returns
- `::Array{Float64, 1}` : smeared reflectometry values for the given q-vectors, with a pointise smearing.
"""
function pointwise_smearing(q, w, dq, quad_order::Int=17)
    # not sure how one checks that q, dq have same size.
    npnts = length(q)

    # get the gauss-legendre weights and abscissae
    # TODO: is it possible to use a LRU cache?
    abscissa, weights = gausslegendre(quad_order)

    # get the normal distribution at that point
    prefactor = 1.0 / sqrt(2 * pi)

    function gauss(x)
        return exp(-0.5 * x * x)
    end

    # TODO: is it possible to use a LRU cache?
    gaussvals = prefactor * gauss.(abscissa * _INTLIMIT) .* weights

    # integration between -3.5 and 3.5 sigma
    va = q .- _INTLIMIT .* dq ./ _FWHM
    vb = q .+ _INTLIMIT .* dq ./ _FWHM

    # (quad_order, npnts)
    qvals_for_res = [(ab * (li[2] - li[1]) + li[2] + li[1])/2.0 for ab=abscissa, li=zip(va, vb)]
    smeared_rvals = abeles(qvals_for_res, w)

    # abeles flattens the q vector
    smeared_rvals = reshape(smeared_rvals, size(qvals_for_res))
    m = reshape(gaussvals, (quad_order, 1))

    smeared_rvals .*= m
    rvals = sum(smeared_rvals, dims=1) .* _INTLIMIT
    return reshape(rvals, size(q))
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
function constant_smearing(q, w, resolution, scale, bkg)
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
