using LinearAlgebra, Kronecker, Combinatorics

include("utils.jl")

X = [0 1; 1 0]
Y = [0 -1im; 1im 0]
Z = [1 0; 0 -1]

function phase(n, a)
    ϕ = 0
    for k in 1:n
        ϕ -= getkthBit(a, 2*k-1) * getkthBit(a, 2*k)
    end
    return float(1im)^ϕ
end

"""
    Pauli(n, a)

Return the Pauli operator on n qubits associated with the vector a.
"""
function pauli(n, a)
    operator = Z^getkthBit(a, 1) * X^getkthBit(a, 2)
    for k in 2:n
        localOperator = Z^getkthBit(a, 2*k-1) * X^getkthBit(a, 2*k)
        operator = operator ⊗ localOperator
    end
    return phase(n, a) * operator
end

"""
    beta(n, a, b)

The function beta tracks how commuting Pauli operators compose according to
T_aT_b=(-1)^betaT_{a+b}.
"""
function beta(n, a, b)
    A = pauli(n, a)
    B = pauli(n, b)
    if A*B == pauli(n, a⊻b)
        return 0
    end
    return 1
end

function isClosed(n, set)
    for k in with_replacement_combinations(set, 2)
        if symplecticForm(n, k[1], k[2]) == 0
            if k[1] ⊻ k[2] ∉ set
                return false
            end
        end
    end
    return true
end

function inferProduct(n, set)
    if size(set)[1] == 2
        x1 = set[1]
        x2 = set[2]
    else
        x1 = set[1]
        x2 = inferProduct(set[2:end])
    end
    return [x1[1] ⊻ x2[1], x1[2] ⊻ x2[2] ⊻ beta(n, x1[1], x2[1])]
end

"""
    getStabilizerState(n, generators, signs)

Given a set of generators of an isotropic subspace in ℤ₂ⁿ×ℤ₂ⁿ with signs, return
the corresponding stabilizer state represented as a density matrix.
"""
function getStabilizerState(n, generators, generatorSigns)
    subspace = vcat([0], generators)
    signs = vcat([0], generatorSigns)
    while !isClosed(n, subspace)
        for k in with_replacement_combinations(1:size(subspace)[1], 2)
            if subspace[k[1]] ⊻ subspace[k[2]] ∉ subspace
                newVector, newSign = inferProduct(n, [[subspace[k[1]], signs[k[1]]], [subspace[k[2]], signs[k[2]]]])
                push!(subspace, newVector)
                push!(signs, newSign)
            end
        end
    end
    dm = (-1)^(signs[1]) * pauli(n, subspace[1])
    for k in 2:size(subspace)[1]
        dm += (-1)^(signs[k]) * pauli(n, subspace[k])
    end
    return dm / 2^n
end
