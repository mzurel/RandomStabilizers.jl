function hammingWeight(x)
  y = 0
  while x > 0
    y += x & 1
    x = x >>> 1
  end
  return y
end

function getkthBit(x, k)
    return (x >>> (k-1)) & 1
end

function symplecticForm(nQubits,a,b)
    J1 = reduce(⊻, 2^(2*k) for k in 0:(nQubits-1))
    J2 = (2^(2*nQubits)-1) ⊻ J1
    return hammingWeight(a & (((b >>> 1) & J1) ⊻ ((b << 1) & J2))) % 2
end

function transvection(nQubits, k, v)  # Applies the transvection Z_k to v
    return v ⊻ (symplecticForm(nQubits, k, v) * k)
end

function multiTransvection(nQubits, h, v)
    # Apply all of the transvections in the Array h to all of the vectors in the Array v
    output = v
    for k in h
        output = [transvection(nQubits, k, u) for u in output]
    end
    return output
end

function findTransvection(nQubits, x, y)  # Finds h1, h2 such that y = Z_h1 Z_h2 x
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

function symplecticGroupOrder(nQubits)
    num = BigInt(2)^(nQubits^2)
    for k in 1:nQubits
        num *= BigInt(4)^k - 1
    end
    return num
end

function symplectic(nQubits, i)
    s = (1<<(2*nQubits)) - 1; k = (i % s) + 1; i = floor(BigInt, i/s)
    e1 = 1 << (2*nQubits-1)
    T = findTransvection(nQubits, e1, k)

    e = e1 ⊻ ((i >>> 1) & (2^(2*nQubits-2)-1))
    h0 = multiTransvection(nQubits, T, [e])[1]

    if (i & 1)  == 1
        Tprime = [h0]
    else
        Tprime = [h0,k]
    end

    if nQubits == 1
        return [k, multiTransvection(nQubits, [T;Tprime], [(1 << (2*nQubits-2))])[1]]  # [f1,f2]
    else
        return multiTransvection(nQubits, [T;Tprime], [(1<<(2*nQubits-1)), (1 << (2*nQubits-2)), symplectic(nQubits-1, i >>> (2*nQubits - 1))...])
    end
end

function bigSymplectic(nQubits::BigInt, i::BigInt)
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

function int2bits(nQubits, num)
    vec = zeros(Int8, 2*nQubits)
    for k in 1:2*nQubits
        vec[end-k+1] = (num >>> (k-1)) & 1
    end
    return vec
end

function randomSymplecticGroupElement(nQubits)
    nQubits = BigInt(nQubits)
    i = rand(1:symplecticGroupOrder)
    A = bigSymplectic(nQubits, i)
    A = [int2bit(nQubits, a) for a in A]
    return A
end

function randomStabilizerState(nQubits)
    nQubits = BigInt(nQubits)
    i = rand(1:symplecticGroupOrder(nQubits))
    A = bigSymplectic(nQubits, i)[1:2:2*nQubits]
    A = [int2bits(nQubits, a) for a in A]
    A = [[a[1:2:2*nQubits]; a[2:2:2*nQubits]] for a in A]
    return A
end
