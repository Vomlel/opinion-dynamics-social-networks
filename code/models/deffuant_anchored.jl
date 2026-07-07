module DeffuantAnchored

using Random
using Graphs

export anchored_step!, run_deffuant_anchored, run_deffuant_anchored_series

"""
Ukotveny socialni Deffuantuv model (ASD).

Oproti klasickemu modelu pouziva (1) hladkou pravdepodobnost prijeti nazoru,
(2) socialni kontext v podobe prumeru okoli, (3) heterogenni ovlivnitelnost
podle stupne uzlu a (4) slaby navrat k puvodnimu nazoru aktéra.
"""
function _susceptibility(g::AbstractGraph, degree_inertia::Float64)
    maxdeg = maximum(degree(g); init=0)
    maxdeg == 0 && return ones(Float64, nv(g))
    return [1.0 / (1.0 + degree_inertia * degree(g, i) / maxdeg) for i in vertices(g)]
end

function _local_mean(x::AbstractVector{<:Real}, g::AbstractGraph, i::Int)
    ns = neighbors(g, i)
    isempty(ns) && return Float64(x[i])
    return sum(Float64(x[j]) for j in ns) / length(ns)
end

"""
    anchored_step!(x, anchors, g, edge_list, susceptibility; epsilon, mu,
                   anchor_strength=0.02, social_weight=0.25,
                   confidence_shape=2.0, rng=Random.default_rng())

Vybere nahodnou hranu. Interakce se uskutecni s pravdepodobnosti
`exp(-(abs(xu-xv)/epsilon)^confidence_shape)`; epsilon tedy neni tvrda hranice.
Kazdy uzel se posune k mixu nazoru partnera a lokalniho prumeru a soucasne je
slabe pritahovan ke svemu pocatecnimu nazoru (`anchors`).
"""
function anchored_step!(
    x::Vector{Float64}, anchors::Vector{Float64}, g::AbstractGraph, edge_list,
    susceptibility::Vector{Float64}; epsilon::Real, mu::Real,
    anchor_strength::Real=0.02, social_weight::Real=0.25,
    confidence_shape::Real=2.0, rng::AbstractRNG=Random.default_rng(),
    clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0),
)
    isempty(edge_list) && throw(ArgumentError("graf nesmi byt bez hran"))
    eps = Float64(epsilon); muf = Float64(mu); a = Float64(anchor_strength)
    w = Float64(social_weight); q = Float64(confidence_shape)
    eps > 0 || throw(ArgumentError("epsilon musi byt > 0"))
    0 < muf <= 0.5 || throw(ArgumentError("mu musi byt v (0, 0.5]"))
    0 <= a <= 1 || throw(ArgumentError("anchor_strength musi byt v [0,1]"))
    0 <= w <= 1 || throw(ArgumentError("social_weight musi byt v [0,1]"))
    q > 0 || throw(ArgumentError("confidence_shape musi byt > 0"))

    e = edge_list[rand(rng, eachindex(edge_list))]
    u, v = src(e), dst(e)
    xu, xv = x[u], x[v]
    p_accept = exp(-((abs(xu - xv) / eps)^q))
    rand(rng) <= p_accept || return (u, v, false)

    context_u = _local_mean(x, g, u)
    context_v = _local_mean(x, g, v)
    target_u = (1 - w) * xv + w * context_u
    target_v = (1 - w) * xu + w * context_v
    new_u = xu + muf * susceptibility[u] * (target_u - xu) + a * (anchors[u] - xu)
    new_v = xv + muf * susceptibility[v] * (target_v - xv) + a * (anchors[v] - xv)

    if clamp_range !== nothing
        low, high = Float64.(clamp_range)
        low < high || throw(ArgumentError("clamp_range musi mit low < high"))
        new_u, new_v = clamp(new_u, low, high), clamp(new_v, low, high)
    end
    x[u], x[v] = new_u, new_v
    return (u, v, new_u != xu || new_v != xv)
end

function run_deffuant_anchored(g::AbstractGraph, x0::AbstractVector{<:Real};
    iterations::Integer, epsilon::Real, mu::Real=0.2, anchor_strength::Real=0.02,
    social_weight::Real=0.25, degree_inertia::Real=1.0,
    confidence_shape::Real=2.0, rng::AbstractRNG=Random.default_rng(),
    clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0))
    iterations >= 0 || throw(ArgumentError("iterations musi byt >= 0"))
    length(x0) == nv(g) || throw(ArgumentError("delka x0 musi byt rovna nv(g)"))
    degree_inertia >= 0 || throw(ArgumentError("degree_inertia musi byt >= 0"))
    x = Float64.(x0); anchors = copy(x); edge_list = collect(edges(g))
    susceptibility = _susceptibility(g, Float64(degree_inertia))
    for _ in 1:iterations
        anchored_step!(x, anchors, g, edge_list, susceptibility; epsilon=epsilon,
            mu=mu, anchor_strength=anchor_strength, social_weight=social_weight,
            confidence_shape=confidence_shape, rng=rng, clamp_range=clamp_range)
    end
    return x
end

function run_deffuant_anchored_series(g::AbstractGraph, x0::AbstractVector{<:Real};
    iterations::Integer, epsilon::Real, mu::Real=0.2, anchor_strength::Real=0.02,
    social_weight::Real=0.25, degree_inertia::Real=1.0,
    confidence_shape::Real=2.0, rng::AbstractRNG=Random.default_rng(),
    clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0), snapshots::Integer=101)
    iterations >= 0 || throw(ArgumentError("iterations musi byt >= 0"))
    snapshots >= 2 || throw(ArgumentError("snapshots musi byt >= 2"))
    length(x0) == nv(g) || throw(ArgumentError("delka x0 musi byt rovna nv(g)"))
    x = Float64.(x0); anchors = copy(x); edge_list = collect(edges(g))
    susceptibility = _susceptibility(g, Float64(degree_inertia))
    times = round.(Int, range(0, Int(iterations); length=Int(snapshots)))
    series = [copy(x)]; next_idx = 2
    for t in 1:Int(iterations)
        anchored_step!(x, anchors, g, edge_list, susceptibility; epsilon=epsilon,
            mu=mu, anchor_strength=anchor_strength, social_weight=social_weight,
            confidence_shape=confidence_shape, rng=rng, clamp_range=clamp_range)
        while next_idx <= snapshots && t >= times[next_idx]
            push!(series, copy(x)); next_idx += 1
        end
    end
    while length(series) < snapshots; push!(series, copy(x)); end
    return series
end

end
