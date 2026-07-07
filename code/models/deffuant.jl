module Deffuant

using Random
using DataFrames
using Graphs

# Modul řeší čistě simulaci nad grafem. IO a generování počátečních dat zůstává
# mimo něj, aby nevznikaly zbytečné závislosti mezi moduly.

export deffuant_step!,
       run_deffuant,
       run_deffuant_series

function _clamp_value(x::Float64, low::Float64, high::Float64)
    return x < low ? low : (x > high ? high : x)
end

function _deffuant_step_edges!(
    x::Vector{Float64},
    edge_list;
    ϵ::Float64,
    μ::Float64,
    rng::AbstractRNG,
    do_clamp::Bool,
    low::Float64,
    high::Float64,
)
    e = edge_list[rand(rng, eachindex(edge_list))]
    u = src(e)
    v = dst(e)

    u == v && return false

    xu = x[u]
    xv = x[v]

    if abs(xu - xv) <= ϵ
        new_u = xu + μ * (xv - xu)
        new_v = xv + μ * (xu - xv)
        if do_clamp
            new_u = _clamp_value(new_u, low, high)
            new_v = _clamp_value(new_v, low, high)
        end
        x[u] = new_u
        x[v] = new_v
        return true
    end

    return false
end

"""
    deffuant_step!(x::AbstractVector{<:Real}, g::AbstractGraph;
        ϵ::Real,
        μ::Real=0.5,
        rng::AbstractRNG=Random.default_rng()
    ) -> Tuple{Int,Int,Bool}

Provede 1 náhodnou interakci v Deffuantově modelu:
- vybere náhodnou hranu (u,v)
- pokud |x[u]-x[v]| ≤ ϵ, provede aktualizaci:

    x[u] ← x[u] + μ (x[v] - x[u])
    x[v] ← x[v] + μ (x[u]_old - x[v])

Vrací `(u, v, changed)`, kde `changed == true` značí, že se názory
po vybrané interakci skutečně změnily.

Pozn.: `x` typicky drž jako `Vector{Float64}`.
"""
function deffuant_step!(
    x::AbstractVector{<:Real},
    g::AbstractGraph;
    ϵ::Real,
    μ::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
)
    μ = Float64(μ)
    ϵ = Float64(ϵ)
    0.0 < μ <= 0.5 || throw(ArgumentError("μ musí být v intervalu (0, 0.5]"))
    ϵ >= 0.0 || throw(ArgumentError("ϵ musí být >= 0"))

    ne(g) > 0 || throw(ArgumentError("graf nesmí být bez hran"))
    length(x) == nv(g) || throw(ArgumentError("Délka x ($(length(x))) musí být rovna počtu vrcholů v grafu ($(nv(g)))"))

    e = rand(rng, edges(g))
    u = src(e)
    v = dst(e)

    u == v && return (u, v, false)

    xu = Float64(x[u])
    xv = Float64(x[v])

    if abs(xu - xv) <= ϵ
        new_u = xu + μ * (xv - xu)
        new_v = xv + μ * (xu - xv)
        x[u] = new_u
        x[v] = new_v
        return (u, v, true)
    else
        return (u, v, false)
    end
end

"""
    run_deffuant(g::AbstractGraph, x0::AbstractVector{<:Real};
        iterations::Integer,
        ϵ::Real,
        μ::Real=0.5,
        rng::AbstractRNG=Random.default_rng(),
        clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0),
        snapshot_every::Integer=0,
        on_snapshot::Union{Nothing,Function}=nothing
    ) -> Vector{Float64}

Spustí Deffuantův model po `iterations` interakcí a vrátí finální vektor názorů.

Parametry:
- `x0` je počáteční stav délky `nv(g)`.
- `clamp_range`: pokud není `nothing`, po každém kroku ořízne hodnoty do [low, high].
- `snapshot_every`: pokud > 0, volá `on_snapshot(step, x)` pro step=0 a pak každých `snapshot_every` kroků.
- `on_snapshot`: callback, např. ` (step,x)->save_opinion_snapshot("out/op", step, x)`.

Pozn.: `iterations` jsou *interakce* (tj. náhodné vybrané hrany), ne „synchronní iterace přes všechny hrany“.
"""
function run_deffuant(
    g::AbstractGraph,
    x0::AbstractVector{<:Real};
    iterations::Integer,
    ϵ::Real,
    μ::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
    clamp_range::Union{Nothing,Tuple{Real,Real}} = (0.0, 1.0),
    snapshot_every::Integer = 0,
    on_snapshot::Union{Nothing,Function} = nothing,
)
    iterations = Int(iterations)
    iterations >= 0 || throw(ArgumentError("iterations musí být >= 0"))

    n = nv(g)

    length(x0) == n || throw(ArgumentError("Délka x0 ($(length(x0))) musí být rovna n ($n)"))

    x = Float64.(x0)

    # Počáteční stav vystavíme jako snapshot s časem 0.
    if snapshot_every > 0 && on_snapshot !== nothing
        on_snapshot(0, x)
    end

    low = 0.0
    high = 1.0
    do_clamp = clamp_range !== nothing
    if do_clamp
        low = Float64(clamp_range[1])
        high = Float64(clamp_range[2])
        low < high || throw(ArgumentError("clamp_range musí mít low < high"))
    end

    edge_list = collect(edges(g))
    isempty(edge_list) && throw(ArgumentError("graf nesmí být bez hran"))
    μf = Float64(μ)
    ϵf = Float64(ϵ)
    0.0 < μf <= 0.5 || throw(ArgumentError("μ musí být v intervalu (0, 0.5]"))
    ϵf >= 0.0 || throw(ArgumentError("ϵ musí být >= 0"))

    for s in 1:iterations
        _deffuant_step_edges!(x, edge_list; ϵ=ϵf, μ=μf, rng=rng, do_clamp=do_clamp, low=low, high=high)

        if snapshot_every > 0 && on_snapshot !== nothing
            if s % snapshot_every == 0
                on_snapshot(s, x)
            end
        end
    end

    return x
end

"""
    run_deffuant_series(g::AbstractGraph, x0::AbstractVector{<:Real};
        iterations::Integer,
        ϵ::Real,
        μ::Real=0.5,
        rng::AbstractRNG=Random.default_rng(),
        clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0),
        snapshots::Integer=101
    ) -> Vector{Vector{Float64}}

Spustí Deffuantův model a vrátí přímo časovou řadu stavů (včetně počátečního).

- `snapshots` říká kolik stavů chceš uložit včetně t=0 a t=iterations.
  Typicky 101 => 0 + 99 mezistavů + finální.

Interně ukládá přibližně rovnoměrně rozložené snapshoty v čase.
"""
function run_deffuant_series(
    g::AbstractGraph,
    x0::AbstractVector{<:Real};
    iterations::Integer,
    ϵ::Real,
    μ::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
    clamp_range::Union{Nothing,Tuple{Real,Real}} = (0.0, 1.0),
    snapshots::Integer = 101,
)
    iterations = Int(iterations)
    iterations >= 0 || throw(ArgumentError("iterations musí být >= 0"))

    snapshots = Int(snapshots)
    snapshots >= 2 || throw(ArgumentError("snapshots musí být >= 2"))

    n = nv(g)
    length(x0) == n || throw(ArgumentError("Délka x0 ($(length(x0))) musí být rovna počtu vrcholů v grafu ($n)"))

    x = Float64.(x0)
    series = Vector{Vector{Float64}}()
    sizehint!(series, snapshots)
    push!(series, copy(x))

    low = 0.0
    high = 1.0
    do_clamp = clamp_range !== nothing
    if do_clamp
        low = Float64(clamp_range[1])
        high = Float64(clamp_range[2])
        low < high || throw(ArgumentError("clamp_range musí mít low < high"))
    end

    edge_list = collect(edges(g))
    isempty(edge_list) && throw(ArgumentError("graf nesmí být bez hran"))
    μf = Float64(μ)
    ϵf = Float64(ϵ)
    0.0 < μf <= 0.5 || throw(ArgumentError("μ musí být v intervalu (0, 0.5]"))
    ϵf >= 0.0 || throw(ArgumentError("ϵ musí být >= 0"))

    # Předem si určíme časy snapshotů, aby série měla přesně danou délku.
    snap_iterations = round.(Int, range(0, iterations; length=snapshots))
    next_idx = 2
    next_snap = snap_iterations[next_idx]

    for s in 1:iterations
        _deffuant_step_edges!(x, edge_list; ϵ=ϵf, μ=μf, rng=rng, do_clamp=do_clamp, low=low, high=high)

        if s == next_snap
            push!(series, copy(x))
            next_idx += 1
            if next_idx <= snapshots
                next_snap = snap_iterations[next_idx]
            end
        end
    end

    while length(series) < snapshots
        push!(series, copy(x))
    end

    return series
end

end # module
