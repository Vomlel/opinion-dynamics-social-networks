module DeffuantPartisan

using Random
using Statistics
using Graphs

export partisan_step!, run_deffuant_partisan, run_deffuant_partisan_series

"""
Partisan Identity Deffuant (PID): model stranickeho trideni.

Identita aktéra je na zacatku urcena jeho pozici vuci medianu. Kontakty uvnitr
tabora vedou ke konvergenci, blizky kontakt s opacnym taborem ke kontrastu.
Slaby tah ke stranickemu centru modeluje stranicke trideni; tah k pocatecnimu
nazoru zachovava individualni heterogenitu.
"""
function _setup(g::AbstractGraph, x0::Vector{Float64}, degree_inertia::Float64)
    split = median(x0)
    camp = [x <= split ? -1.0 : 1.0 for x in x0]
    maxdeg = maximum(degree(g); init=0)
    susceptibility = maxdeg == 0 ? ones(length(x0)) :
        [1 / (1 + degree_inertia * degree(g, i) / maxdeg) for i in vertices(g)]
    return camp, susceptibility
end

function partisan_step!(x::Vector{Float64}, anchors::Vector{Float64}, camp::Vector{Float64},
    susceptibility::Vector{Float64}, edge_list; epsilon::Real, mu::Real,
    cross_epsilon::Real=0.35, repulsion::Real=0.35, identity_strength::Real=0.015,
    anchor_strength::Real=0.01, left_pole::Real=0.25, right_pole::Real=0.75,
    rng::AbstractRNG=Random.default_rng(), clamp_range::Tuple{Real,Real}=(0.0, 1.0))

    eps, muf, ce = Float64(epsilon), Float64(mu), Float64(cross_epsilon)
    0 < eps || throw(ArgumentError("epsilon musi byt > 0"))
    0 < muf <= 0.5 || throw(ArgumentError("mu musi byt v (0,0.5]"))
    ce > 0 || throw(ArgumentError("cross_epsilon musi byt > 0"))
    repulsion >= 0 || throw(ArgumentError("repulsion musi byt >= 0"))
    identity_strength >= 0 || throw(ArgumentError("identity_strength musi byt >= 0"))
    anchor_strength >= 0 || throw(ArgumentError("anchor_strength musi byt >= 0"))
    isempty(edge_list) && throw(ArgumentError("graf nesmi byt bez hran"))

    e = edge_list[rand(rng, eachindex(edge_list))]
    u, v = src(e), dst(e)
    xu, xv = x[u], x[v]
    peer_u = 0.0; peer_v = 0.0
    distance = abs(xu - xv)

    if camp[u] == camp[v]
        if distance <= eps
            peer_u = muf * (xv - xu)
            peer_v = muf * (xu - xv)
        end
    elseif distance < ce
        # Diferenciace: nejvetsi je, kdyz se protistrany nazorove priblizi.
        force = muf * Float64(repulsion) * (ce - distance)
        peer_u = camp[u] * force
        peer_v = camp[v] * force
    end

    pole_u = camp[u] < 0 ? Float64(left_pole) : Float64(right_pole)
    pole_v = camp[v] < 0 ? Float64(left_pole) : Float64(right_pole)
    new_u = xu + susceptibility[u] * (peer_u + identity_strength * (pole_u - xu)) +
            anchor_strength * (anchors[u] - xu)
    new_v = xv + susceptibility[v] * (peer_v + identity_strength * (pole_v - xv)) +
            anchor_strength * (anchors[v] - xv)
    low, high = Float64.(clamp_range)
    x[u], x[v] = clamp(new_u, low, high), clamp(new_v, low, high)
    return (u, v, x[u] != xu || x[v] != xv)
end

function _run(g, x0; iterations, epsilon, mu, cross_epsilon, repulsion,
    identity_strength, anchor_strength, degree_inertia, left_pole, right_pole,
    rng, clamp_range, snapshots)
    iterations >= 0 || throw(ArgumentError("iterations musi byt >= 0"))
    length(x0) == nv(g) || throw(ArgumentError("delka x0 musi byt rovna nv(g)"))
    0 <= left_pole < right_pole <= 1 || throw(ArgumentError("musi platit 0 <= left_pole < right_pole <= 1"))
    x = Float64.(x0); anchors = copy(x); edge_list = collect(edges(g))
    camp, susceptibility = _setup(g, x, Float64(degree_inertia))
    times = snapshots === nothing ? Int[] : round.(Int, range(0, iterations; length=snapshots))
    series = snapshots === nothing ? Vector{Vector{Float64}}() : [copy(x)]
    next_idx = 2
    for t in 1:Int(iterations)
        partisan_step!(x, anchors, camp, susceptibility, edge_list; epsilon=epsilon, mu=mu,
            cross_epsilon=cross_epsilon, repulsion=repulsion, identity_strength=identity_strength,
            anchor_strength=anchor_strength, left_pole=left_pole, right_pole=right_pole,
            rng=rng, clamp_range=clamp_range)
        while snapshots !== nothing && next_idx <= snapshots && t >= times[next_idx]
            push!(series, copy(x)); next_idx += 1
        end
    end
    snapshots === nothing && return x
    while length(series) < snapshots; push!(series, copy(x)); end
    return series
end

function run_deffuant_partisan(g::AbstractGraph, x0::AbstractVector{<:Real};
    iterations::Integer, epsilon::Real, mu::Real=0.2, cross_epsilon::Real=0.35,
    repulsion::Real=0.35, identity_strength::Real=0.015, anchor_strength::Real=0.01,
    degree_inertia::Real=0.5, left_pole::Real=0.25, right_pole::Real=0.75,
    rng::AbstractRNG=Random.default_rng(), clamp_range::Tuple{Real,Real}=(0.0,1.0))
    _run(g, x0; iterations=iterations, epsilon=epsilon, mu=mu, cross_epsilon=cross_epsilon,
        repulsion=repulsion, identity_strength=identity_strength, anchor_strength=anchor_strength,
        degree_inertia=degree_inertia, left_pole=left_pole, right_pole=right_pole,
        rng=rng, clamp_range=clamp_range, snapshots=nothing)
end

function run_deffuant_partisan_series(g::AbstractGraph, x0::AbstractVector{<:Real};
    iterations::Integer, epsilon::Real, mu::Real=0.2, cross_epsilon::Real=0.35,
    repulsion::Real=0.35, identity_strength::Real=0.015, anchor_strength::Real=0.01,
    degree_inertia::Real=0.5, left_pole::Real=0.25, right_pole::Real=0.75,
    rng::AbstractRNG=Random.default_rng(), clamp_range::Tuple{Real,Real}=(0.0,1.0), snapshots::Integer=101)
    snapshots >= 2 || throw(ArgumentError("snapshots musi byt >= 2"))
    _run(g, x0; iterations=iterations, epsilon=epsilon, mu=mu, cross_epsilon=cross_epsilon,
        repulsion=repulsion, identity_strength=identity_strength, anchor_strength=anchor_strength,
        degree_inertia=degree_inertia, left_pole=left_pole, right_pole=right_pole,
        rng=rng, clamp_range=clamp_range, snapshots=Int(snapshots))
end

end
