using TensorToolbox, LinearAlgebra, Statistics

struct MultiPCA1{T,N}
    μ::Array{T}
    core::Array{T, N}
    factors::Vector{Matrix{T}}
end

MultiPCA = MultiPCA1
rank(m::MultiPCA) = size(m.core)
Base.show(io::IO, m::MultiPCA) = print(io, "MultiPCA on data $(first.(size.(m.factors))) with rank ", rank(m))

firstrT(r, A) = Matrix(A[:, 1:r]')
firstr(r, A) = A[:, 1:r]

# last dimension is the sampling dimension
function MultiPCA1(A::AbstractArray{T,N}; rank=80, eps=0.01) where {T,N}
    μ = mean(A, dims=N)
    H = hosvd(A .- μ; reqrank=rank, eps_rel=eps)
    return MultiPCA1{T,N}(μ, H.cten, H.fmat)
    # V = ttm(H.cten[1:r,1:r,1:r], [firstr(f) for f in H.fmat])
    # return MultiPCA(V, [firstr(f) for f in H.fmat[1:2]])
end

function transform(mpca::MultiPCA{T,N}, V::AbstractArray, r=rank(mpca)) where {T,N}
    B = V .- mpca.μ
    return ttm(B, firstrT.(r[1:end-1], mpca.factors[1:end-1]))
end

function inverse(mpca::MultiPCA{T,N}, C::AbstractArray) where {T,N}
    r = size(C, 1)
    B = ttm(C, firstr.(Ref(r), mpca.factors[1:end-1]))
    return B .+ mpca.μ
end
