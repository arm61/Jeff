using Interpolations

const TINY = 1e-30
const _FWHM = 2 * sqrt(2 * log(2.0))

"""
    matmul(a::Array{ComplexF64, 2}, b::Array{ComplexF64, 2})

Perform matrix multiplication.

### Parameters
- `a::Array{ComplexF64, 2}` : first matrix.
- `b::Array{ComplexF64, 2}` : second matrix.

### Returns 
- `::Array{ComplexF64, 2}` : result of multiplication.
"""
function matmul(a::Array{ComplexF64, 2}, b::Array{ComplexF64, 2})
    c = zeros(ComplexF64, (2, 2))
    c[1, 1] = a[1, 1] * b[1, 1] + a[1, 2] * b[2, 1]
    c[1, 2] = a[1, 1] * b[1, 2] + a[1, 2] * b[2, 2]
    c[2, 1] = a[2, 1] * b[1, 1] + a[2, 2] * b[2, 1]
    c[2, 2] = a[2, 1] * b[1, 2] + a[2, 2] * b[2, 2]
    return c
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
function abeles(q::Array{Float64, 1}, layers::Array{Float64, 2})
    nlayers = size(layers, 2) - 2
    npnts = size(q, 1)
    reflectivity = zeros(Float64, npnts)
    
    sub = layers[end, 2] + (1im * layers[end, 3] + TINY)
    super = layers[1, 2] + 0im
    sld = zeros(ComplexF64, nlayers + 2)
    thick = zeros(ComplexF64, nlayers)
    rough = zeros(ComplexF64, nlayers)

    for ii = 2:nlayers
        _t = layers[ii, 2] + (1im * layers[ii, 3] + TINY)
        sld[ii] = 4e-6 * pi * (_t - super)
        thick[ii - 1] = 0 + 1im * layers[ii, 1]
        rough[ii - 1] = -2 * layers[ii, 4] ^ 2
    end
    sld[1] = 0 + 0im
    sld[nlayers+1] = 4e-6 * pi * (sub - super)
    rough[nlayers] = -2 * layers[end, 4] ^ 2
    
    mr = zeros(ComplexF64, (2, 2))
    mi = zeros(ComplexF64, (2, 2))

    for j = 1:npnts
        q2 = q[j] ^ 2 / 4. + 0im
        kn = q[j] / 2. + 0im
        for ii = 1:nlayers # maybe nlayers + 1
            knn = sqrt(q2 - sld[ii + 1])
            rj = (kn - knn) / (kn + knn) * exp(kn * knn * rough[ii])
        
            if ii == 1
                mr[1, 1] = 1. + 0im
                mr[1, 2] = rj
                mr[2, 2] = 1. + 0im
                mr[2, 1] = rj
            else
                beta = exp(kn * thick[ii - 1])
                mi[1, 1] = beta
                mi[1, 2] = rj * beta
                mi[2, 2] = (1. + 0im) / beta
                mi[2, 1] = rj * mi[2, 2]
                mr = matmul(mr, mi)
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

### Parameters
- `x::Array{Float64, 1}` : the kernal positions.
- `s::Float64` : the width of the kernel.

### Returns
- `::Array{Float64, 1}`: probabilities for `x`. 
"""
function gauss(x::Array{Float64, 1}, s::Float64) 
    g = 1. / s / sqrt(2 * pi) * exp.(-0.5 * x .^ 2 / s / s)
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
function constant_smearing(q::Array{Float64, 1}, w::Array{Float64, 2}; resolution::Float64=5., scale::Float64=1., bkg::Float64=0.)
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