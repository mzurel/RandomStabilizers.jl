module RandomStabilizers

include("utils.jl")
include("symplectic.jl")
include("Paulis.jl")

export randomStabilizerState, randomSymplecticMatrix, symplecticGroupOrder

"""
    randomSymplecticMatrix(nQubits)

Return a random element of the symplectic group Sp(2n,ℤ₂).  This is obtained by
generating a random integer i between 1 and |Sp(2n,ℤ₂)| and returning symplectic(nQubits, i).
"""
function randomSymplecticMatrix(nQubits)
    nQubits = BigInt(nQubits)
    i = rand(1:symplecticGroupOrder(nQubits))
    A = symplectic(nQubits, i)
    A = [int2bits(nQubits, a) for a in A]
    return A
end

"""
    randomStabilizerState(nQubits)

Returns the check matrix of a random n-qubit stabilizer state.  This is obtained
by generating a random element T of Sp(2n,ℤ₂) with the function randomSymplecticMatrix
and returning the Matrix with Te₁,Te₂,... as rows where e₁,e₂,... are n pairwise
orthogonal standard symplectic basis vectors.
"""
function randomCheckMatrix(nQubits)
    nQubits = BigInt(nQubits)
    i = rand(1:symplecticGroupOrder(nQubits))
    A = symplectic(nQubits, i)[1:2:2*nQubits]
    A = [int2bits(nQubits, a) for a in A]
    A = [[a[1:2:2*nQubits]; a[2:2:2*nQubits]] for a in A]
    return A
end

"""
    randomStabilizerState(nQubits)

Returns a random stabilizer state represented by a density matrix.
"""
function randomStabilizerState(nQubits)
    nQubits = BigInt(nQubits)
    i = rand(1:symplecticGroupOrder(nQubits))
    generators = symplectic(nQubits, i)[1:2:2*nQubits]
    signs = rand([0,1], nQubits)
    return getStabilizerState(nQubits, generators, signs)
end

end
