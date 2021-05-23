"""
    hammingWeight(x)

Return the number of ones in the binary expansion of the integer x.
"""
function hammingWeight(x)
  y = 0
  while x > 0
    y += x & 1
    x = x >>> 1
  end
  return y
end

"""
    getkthbit(x, k)

Return the kth bit in the binary expansion of the integer k.
"""
function getkthBit(x, k)
    return (x >>> (k-1)) & 1
end

"""
    int2bits(nQubits, num)

Return the binary expansion of the integer x as a bitstring of length 2*nQubits.
"""
function int2bits(nQubits, x)
    vec = zeros(Int8, 2*nQubits)
    for k in 1:2*nQubits
        vec[end-k+1] = (x >>> (k-1)) & 1
    end
    return vec
end

"""
    symplecticForm(nQubits, a, b)

Returns the symplectic form of a,b ∈ ℤ₂ⁿ×ℤ₂ⁿ defined by
∑_{k=1}^{n} a[k]*b[n+k] + a[n+k]*b[k] mod 2.
"""
function symplecticForm(nQubits, a, b)
    J1 = reduce(⊻, 2^(2*k) for k in 0:(nQubits-1))
    J2 = (2^(2*nQubits)-1) ⊻ J1
    return hammingWeight(a & (((b >>> 1) & J1) ⊻ ((b << 1) & J2))) % 2
end
