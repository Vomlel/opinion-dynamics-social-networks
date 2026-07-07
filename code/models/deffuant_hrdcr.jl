module DeffuantHRDCR

# Strukturálně rozšířená varianta Deffuantova modelu (HRDCR).
using Random
using DataFrames

export deffuant_hrdcr_step!,
       run_deffuant_hrdcr,
       run_deffuant_hrdcr_series

"""Interní helper: odhad počtu uzlů z edge listu (1-based indexace)."""
_infer_n(edges::DataFrame) = max(maximum(edges.source), maximum(edges.target))

"""Interní helper: seznam sousedů pro neorientovaný graf."""
function _adjlist_undirected(edges::DataFrame, n::Int)
    adj = [Int[] for _ in 1:n]
    @inbounds for i in 1:nrow(edges)
        u = Int(edges.source[i])
        v = Int(edges.target[i])
        u == v && continue
        push!(adj[u], v)
        push!(adj[v], u)
    end
    return adj
end

"""Interní helper: Jaccard mezi sousedstvími uzlů u a v."""
function _jaccard_neighbors(u::Int, v::Int, adj)::Float64
    Nu = adj[u]
    Nv = adj[v]
    isempty(Nu) && isempty(Nv) && return 0.0
    Su = Set(Nu)
    Sv = Set(Nv)
    inter = length(intersect(Su, Sv))
    uni = length(union(Su, Sv))
    return uni == 0 ? 0.0 : inter / uni
end

"""Interní helper: předpočítá míru setrvačnosti uzlů a Jaccard pro každou hranu."""
function _precompute_hrdcr(edges::DataFrame, n::Int; s_max::Real, s_beta::Real)
    adj = _adjlist_undirected(edges, n)

    deg = [length(adj[i]) for i in 1:n]
    maxdeg = maximum(deg)

    stubbornness = zeros(Float64, n)
    if maxdeg > 0
        @inbounds for i in 1:n
            x = (deg[i] / maxdeg)^Float64(s_beta)
            stubbornness[i] = min(Float64(s_max), x * Float64(s_max))
        end
    end

    m = nrow(edges)
    jacc = Vector{Float64}(undef, m)
    @inbounds for i in 1:m
        u = Int(edges.source[i])
        v = Int(edges.target[i])
        jacc[i] = _jaccard_neighbors(u, v, adj)
    end

    return stubbornness, jacc
end

"""
    deffuant_hrdcr_step!(x, edges, stubbornness, jacc;
        d::Real,
        μ::Real=0.5,
        eps_j::Real=0.2,
        gamma_j::Real=0.5,
        rng::AbstractRNG=Random.default_rng()
    ) -> Tuple{Int,Int,Bool}

Provede jednu náhodnou interakci modelu Deffuant-HRDCR:
1. Vybere náhodnou hranu (u, v).
2. Spočítá prah důvěry pro danou hranu:
   `d_uv = d * (eps_j + (1 - eps_j) * J_uv)`.
3. Pokud `|x[u] - x[v]| <= d_uv`, provede aktualizaci s komunitním posílením a
   uzlovou setrvačností (hub nodes se hýbou méně).
"""
function deffuant_hrdcr_step!(
    x::AbstractVector{<:Real},
    edges::DataFrame,
    stubbornness::AbstractVector{<:Real},
    jacc::AbstractVector{<:Real};
    d::Real,
    μ::Real = 0.5,
    eps_j::Real = 0.2,
    gamma_j::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
)
    m = nrow(edges)
    m > 0 || throw(ArgumentError("edges nesmí být prázdné"))

    d = Float64(d)
    μ = Float64(μ)
    eps_j = Float64(eps_j)
    gamma_j = Float64(gamma_j)
    0.0 <= d || throw(ArgumentError("d musí být >= 0"))
    0.0 < μ <= 0.5 || throw(ArgumentError("μ musí být v intervalu (0, 0.5]"))
    0.0 < eps_j <= 1.0 || throw(ArgumentError("eps_j musí být v intervalu (0, 1]"))

    k = rand(rng, 1:m)
    u = Int(edges.source[k])
    v = Int(edges.target[k])
    u == v && return (u, v, false)

    xu = Float64(x[u])
    xv = Float64(x[v])
    Δ = xv - xu

    j = Float64(jacc[k])
    d_uv = d * (eps_j + (1.0 - eps_j) * j)

    if abs(Δ) <= d_uv
        boost = 1.0 + gamma_j * j
        baseδ = μ * Δ * boost
        su = Float64(stubbornness[u])
        sv = Float64(stubbornness[v])
        x[u] = xu + (1.0 - su) * baseδ
        x[v] = xv - (1.0 - sv) * baseδ
        return (u, v, true)
    end

    return (u, v, false)
end

"""
    run_deffuant_hrdcr(edges, x0; kwargs...) -> Vector{Float64}

Spustí Deffuant-HRDCR model po `steps` interakcí a vrátí finální názory.

Klíčové parametry:
- `d`: základní confidence bound.
- `s_max`, `s_beta`: síla a tvar uzlové setrvačnosti.
- `eps_j`, `gamma_j`: vliv komunitního překryvu (Jaccard) na práh a krok.
"""
function run_deffuant_hrdcr(
    edges::DataFrame,
    x0::AbstractVector{<:Real};
    steps::Integer,
    d::Real,
    μ::Real = 0.5,
    s_max::Real = 0.7,
    s_beta::Real = 1.0,
    eps_j::Real = 0.2,
    gamma_j::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
    n::Union{Int,Nothing} = nothing,
    clamp_range::Union{Nothing,Tuple{Real,Real}} = (0.0, 1.0),
)
    steps = Int(steps)
    steps >= 0 || throw(ArgumentError("steps musí být >= 0"))

    n === nothing && (n = _infer_n(edges))
    n = Int(n)
    length(x0) == n || throw(ArgumentError("Délka x0 ($(length(x0))) musí být rovna n ($n)"))

    x = Float64.(x0)
    stubbornness, jacc = _precompute_hrdcr(edges, n; s_max=s_max, s_beta=s_beta)

    do_clamp = clamp_range !== nothing
    low = 0.0
    high = 1.0
    if do_clamp
        low = Float64(clamp_range[1])
        high = Float64(clamp_range[2])
        low < high || throw(ArgumentError("clamp_range musí mít low < high"))
    end

    for _ in 1:steps
        deffuant_hrdcr_step!(x, edges, stubbornness, jacc; d=d, μ=μ, eps_j=eps_j, gamma_j=gamma_j, rng=rng)
        if do_clamp
            @inbounds for i in 1:n
                xi = x[i]
                x[i] = xi < low ? low : (xi > high ? high : xi)
            end
        end
    end

    return x
end

"""
    run_deffuant_hrdcr_series(edges, x0; kwargs...) -> Vector{Vector{Float64}}

Stejné jako `run_deffuant_hrdcr`, ale vrací časovou řadu stavů
(včetně počátečního stavu) v `snapshots` rovnoměrných časech.
"""
function run_deffuant_hrdcr_series(
    edges::DataFrame,
    x0::AbstractVector{<:Real};
    steps::Integer,
    d::Real,
    μ::Real = 0.5,
    s_max::Real = 0.7,
    s_beta::Real = 1.0,
    eps_j::Real = 0.2,
    gamma_j::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
    n::Union{Int,Nothing} = nothing,
    clamp_range::Union{Nothing,Tuple{Real,Real}} = (0.0, 1.0),
    snapshots::Integer = 101,
)
    steps = Int(steps)
    steps >= 0 || throw(ArgumentError("steps musí být >= 0"))
    snapshots = Int(snapshots)
    snapshots >= 2 || throw(ArgumentError("snapshots musí být >= 2"))

    n === nothing && (n = _infer_n(edges))
    n = Int(n)
    length(x0) == n || throw(ArgumentError("Délka x0 ($(length(x0))) musí být rovna n ($n)"))

    x = Float64.(x0)
    stubbornness, jacc = _precompute_hrdcr(edges, n; s_max=s_max, s_beta=s_beta)

    series = Vector{Vector{Float64}}()
    sizehint!(series, snapshots)
    push!(series, copy(x))

    do_clamp = clamp_range !== nothing
    low = 0.0
    high = 1.0
    if do_clamp
        low = Float64(clamp_range[1])
        high = Float64(clamp_range[2])
        low < high || throw(ArgumentError("clamp_range musí mít low < high"))
    end

    snap_steps = round.(Int, range(0, steps; length=snapshots))
    next_idx = 2
    next_snap = snap_steps[next_idx]

    for s in 1:steps
        deffuant_hrdcr_step!(x, edges, stubbornness, jacc; d=d, μ=μ, eps_j=eps_j, gamma_j=gamma_j, rng=rng)
        if do_clamp
            @inbounds for i in 1:n
                xi = x[i]
                x[i] = xi < low ? low : (xi > high ? high : xi)
            end
        end

        if s == next_snap
            push!(series, copy(x))
            next_idx += 1
            if next_idx <= snapshots
                next_snap = snap_steps[next_idx]
            end
        end
    end

    while length(series) < snapshots
        push!(series, copy(x))
    end

    return series
end

end # module
