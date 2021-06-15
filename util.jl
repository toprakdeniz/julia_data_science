euclidean_distance(p1::AbstractVector{<:AbstractFloat},p2::AbstractVector{<:AbstractFloat}) = sqrt(sum((p1 .- p2).^2))

function gaussian_kernel(distance, bandwidth)
 euclidian_dist = sum(distance.^2)
 (1/ (bandwidth* sqrt(2*π))) * exp(-0.5*(euclidian_dist/bandwidth).^2)
end
