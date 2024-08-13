using Images

function registration(X::AbstractArray{<:Number, 3}, center, radius, search_radius=div(radius, 3))
    nframes = size(X, 3)
    w = -search_radius:search_radius
    δ = CartesianIndices((w, w))
    window = CartesianIndices((-radius:radius, -radius:radius))  # viewing window
    centers = fill(Tuple(center), nframes)
    indices = fill(window .+ CartesianIndex(Tuple(center)), nframes)
    for frame in nframes-1:-1:1
        center = CartesianIndex(centers[frame + 1])
        X0 = view(X, window .+ center, frame+1)  # reference to be matched
        dif(shift) = view(X, window .+ (center + shift), frame) - X0  # difference from reference
        idx = argmin( sum(abs, dif(shift)) for shift in δ )
        centers[frame] = Tuple(center + δ[idx])
        indices[frame] = window .+ (center + δ[idx])
    end
    return indices, centers
end
registration(X::AbstractArray{<:Gray}, args...) = registration(Float16.(X), args...)
registration(X::AbstractArray{<:RGB}, args...) = registration(Gray.(X), args...)
registration(X::AbstractVector{<:AbstractMatrix}, args...) = registration(cat(X..., dims=3), args...)


## Demo for me
X = load.(readdir("/Users/driscoll/Downloads/t1_extracted", join=true));
X = X[17:93];  # 77 usable frames
data, centers = registration(X, (585, 460), 30, 6);
X[77][data[77]]  # last frame
X[67][data[67]]
