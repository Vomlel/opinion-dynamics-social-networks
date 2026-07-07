module StaticCharacteristics

using DataFrames
using Graphs
using Random
using Statistics

export average_degree,
       degree_distribution,
       density,
       average_shortest_path_length


_infer_n(edges::DataFrame) = max(maximum(edges.source), maximum(edges.target))


function _adjlist_undirected(edges::DataFrame, n::Int)
    adj = [Int[] for _ in 1:n]
    @inbounds for i in 1:nrow(edges)
        u = edges.source[i]
        v = edges.target[i]
        if u == v
            continue
        end
        push!(adj[u], v)
        push!(adj[v], u)
    end
    return adj
end


function _degrees(adj::Vector{Vector{Int}})
    deg = Vector{Int}(undef, length(adj))
    @inbounds for i in eachindex(adj)
        deg[i] = length(adj[i])
    end
    return deg
end


function average_degree(g::AbstractGraph)
    n = nv(g)
    m = ne(g)
    n <= 1 && return 0.0
    return 2.0 * m / n
end


function degree_distribution(g::AbstractGraph)
    n = nv(g)
    deg = degree(g)

    counts = Dict{Int,Int}()
    @inbounds for k in deg
        counts[k] = get(counts, k, 0) + 1
    end

    ks = sort(collect(keys(counts)))
    cnt = [counts[k] for k in ks]
    frac = Float64.(cnt) ./ n

    return DataFrame(degree = ks, count = cnt, fraction = frac)
end

"""
    density(edges::DataFrame; n::Union{Int,Nothing}=nothing) -> Float64

Hustota neorientované sítě:

    ρ = 2m / (n(n-1))

kde `m` je počet hran a `n` počet uzlů.
"""
function density(g::AbstractGraph)
    n = nv(g)
    m = ne(g)
    n <= 1 && return 0.0
    return 2.0 * m / (n * (n - 1))
end

function average_shortest_path_length(g::AbstractGraph)
    n = nv(g)

    total = 0.0
    pairs = 0

    for s in 1:n
        dist = gdistances(g, s)

        @inbounds for v in 1:n
            if v != s && isfinite(dist[v])
                total += dist[v]
                pairs += 1
            end
        end
    end

    return total / pairs
end

end # module
