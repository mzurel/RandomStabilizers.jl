module RandomStabilizers

export randomStabilizerState, randomSymplecticMatrix, symplecticGroupOrder

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
    symplecticForm(nQubits, a, b)

Returns the symplectic form of a,b ∈ ℤ₂ⁿ×ℤ₂ⁿ defined by
∑_{k=1}^{n} a[k]*b[n+k] + a[n+k]*b[k] mod 2.
"""
function symplecticForm(nQubits,a,b)
    J1 = reduce(⊻, 2^(2*k) for k in 0:(nQubits-1))
    J2 = (2^(2*nQubits)-1) ⊻ J1
    return hammingWeight(a & (((b >>> 1) & J1) ⊻ ((b << 1) & J2))) % 2
end

"""
    transvection(nQubits, k, v)

k and v are integers whose binary expansions represent vectors in ℤ₂ⁿ×ℤ₂ⁿ.
This function returns the symplectic transvection Zₖ[v]:=v+<k,v>k, represented
as the binary expansion of an integer.
"""
function transvection(nQubits, k, v)
    return v ⊻ (symplecticForm(nQubits, k, v) * k)
end

"""
    multiTransvection(nQubits, h, v)

Apply all of the transvections in the Array h to all of the vectors in the Array v.
"""
function multiTransvection(nQubits, h, v)
    output = v
    for k in h
        output = [transvection(nQubits, k, u) for u in output]
    end
    return output
end

"""
    findTransvection(nQubits, x, y)

This function finds vectors h₁,h₂ ∈ ℤ₂ⁿ×ℤ₂ⁿ such that y=Z_h₁ Z_h₂ x.  This procedure
is described in the proof of Lemma 2 in J. Math. Phys. 55, 122202 (2014).
"""
function findTransvection(nQubits, x, y)
    if x == y
        return [0,0]
    elseif symplecticForm(nQubits, x, y) == 1
        return [x⊻y,0]
    end
    for j in 1:nQubits
        x1 = getkthBit(x, 2*j-1)
        x2 = getkthBit(x, 2*j)
        y1 = getkthBit(y, 2*j-1)
        y2 = getkthBit(y, 2*j)
        if ((x1 | x2) & (y1 | y2)) == 1
            for v in 0:3
                if (((v&1)*x1 ⊻ ((v>>>1)&1)*x2) == 1) && ((v&1)*y1 + ((v>>>1)&1)*y2 == 1)
                    z = x ⊻ (((v>>>1)&1)<<(2*j-2)) ⊻ ((v&1)<<(2*j-1))
                    return [x⊻z,z⊻y]
                end
            end
        end
    end
    for j in 1:nQubits
        x1 = getkthBit(x, 2*j-1)
        x2 = getkthBit(x, 2*j)
        if (x1 | x2) == 1
            for k in 1:nQubits
                y1 = getkthBit(y, 2*k-1)
                y2 = getkthBit(y, 2*k)
                if (y1 | y2) == 1
                    for v in 0:15
                        if (((v&1)*x1 + ((v>>>1)&1)*x2) == 1) && (((v>>>2)&1)*y1 + ((v>>>3)&1)*y2== 1)
                            z = x ⊻ (((v>>>1)&1)*2^(2*j-2)) ⊻ ((v&1)*2^(2*j-1)) ⊻ (((v>>>3)&1)*2^(2*k-2)) ⊻ (((v>>>2)&1)*2^(2*k-1))
                            return [x⊻z,z⊻y]
                        end
                    end
                end
            end
        end
    end
end

"""
    symplecticGroupOrder(nQubits)

Returns the order of the symplectic group Sp(2n,ℤ₂) according to
https://groupprops.subwiki.org/wiki/Order_formulas_for_symplectic_groups
"""
function symplecticGroupOrder(nQubits)
    num = BigInt(2)^(nQubits^2)
    for k in 1:nQubits
        num *= BigInt(4)^k - 1
    end
    return num
end

"""
    symplectic(nQubits, i)

Returns the symplectic group element uniquely identified with the integer
1≤i≤|Sp(2n,ℤ₂)|, according to the algorithm SYMPLECTICImproved from
J. Math. Phys. 55, 122202 (2014).
"""
function symplectic(nQubits::BigInt, i::BigInt)
    s = (BigInt(1)<<(2*nQubits)) - 1; k = (i % s) + 1; i = floor(BigInt, i/s)
    e1 = BigInt(1) << (2*nQubits-1)
    T = findTransvection(nQubits, e1, k)

    e = e1 ⊻ ((i >>> 1) & (BigInt(2)^(2*nQubits-2)-1))
    h0 = multiTransvection(nQubits, T, [e])[1]

    if (i & 1)  == 1
        Tprime = [h0]
    else
        Tprime = [h0,k]
    end

    if nQubits == 1
        return [k, multiTransvection(nQubits, [T;Tprime], [(BigInt(1) << (2*nQubits-2))])[1]]  # [f1,f2]
    else
        return multiTransvection(nQubits, [T;Tprime], [(BigInt(1)<<(2*nQubits-1)), (BigInt(1) << (2*nQubits-2)), bigSymplectic(nQubits-1, i >>> (2*nQubits - 1))...])
    end
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
    randomSymplecticMatrix(nQubits)

Return a random element of the symplectic group Sp(2n,ℤ₂).  This is obtained by
generating a random integer i between 1 and |Sp(2n,ℤ₂)| and returning symplectic(nQubits, i).
"""
function randomSymplecticMatrix(nQubits)
    nQubits = BigInt(nQubits)
    i = rand(1:symplecticGroupOrder)
    A = symplectic(nQubits, i)
    A = [int2bit(nQubits, a) for a in A]
    return A
end

"""
    randomStabilizerState(nQubits)

Returns the check matrix of a random n-qubit stabilizer state.  This is obtained
by generating a random element T of Sp(2n,ℤ₂) with the function randomSymplecticMatrix
and returning the Matrix with Te₁,Te₂,... as rows where e₁,e₂,... are n pairwise
orthogonal standard symplectic basis vectors.
"""
function randomStabilizerState(nQubits)
    nQubits = BigInt(nQubits)
    i = rand(1:symplecticGroupOrder(nQubits))
    A = symplectic(nQubits, i)[1:2:2*nQubits]
    A = [int2bits(nQubits, a) for a in A]
    A = [[a[1:2:2*nQubits]; a[2:2:2*nQubits]] for a in A]
    return A
end

end
