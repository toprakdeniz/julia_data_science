# kmean
include("util.jl")

const min_change = 0.000001
const group_distance_gap = 0.1
function shifting(point, points,bandwidth, kernel = gaussian_kernel)
    weights = kernel.(point.-points, bandwidth)
    weighted_points = weights .* points
    numerator = sum(weighted_points, dim=1)
    denominator = sum(weights)
    numerator / denominator
end
function mean_shift(points::AbstractVector{:<AbstractFloat}, bandwidth::AbstractFloat, kernel = gaussian_kernel)
    centers = collect(eachrow(points))
    m = lenght(centers)
    is_shifting = trues(m)
    while change_exists
        change_exists = false
        for i in eachindex(centers)
            if @inbounds is_shifting[i] == false
                continue
            end
            new_center = shifting(centers[i], centers, kernel)
            if euclidean_distance(new_center,centers[i]) < min_change
                @inbounds is_shifting[i] = false
            else
                change_exists = true
            end
            centers[i] = new_center
        end
    end
    groups = reduce_close_centers(centers)
    return centers, groups
end


function reduce_centers_and_group_assignments(centers, distance_tolerance = group_distance_gap)
    m = length(centers)
    groups = Vector{Vector{Int64}}()
    group_index = 0
    index_of_rest = collect(1:m)
    while !isempty(index_of_rest)
        group_index += 1
        push(groups, [])
        groups[group_index] = current_group_indexes = [popfirst!(index_of_rest)]
        while !isempty(current_group_indexes) && !isempty(index_of_rest)
            cur = popfirst!(current_group_indexes)
            close_centers =  euclidean_distance.(centers[cur], centers[index_of_rest]) < distance_tolerance
            index_of_rest = index_of_rest[findall(.!close_centers)]
            close_center_indexes = index_of_rest[findall(close_centers)]
            append!(current_group_indexes, close_center_indexes)
            append!(groups[group_index], close_center_indexes)
        end
    end
    return groups
end
