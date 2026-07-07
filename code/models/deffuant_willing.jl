module DeffuantWilling

using Random
using DataFrames

# Optional snapshot support (only used if user passes callbacks)
# We keep this module independent of GenerateOpinion to avoid circular includes.

export deffuant_willing_step!,
       run_deffuant_willing,
       run_deffuant_willing_series

"""Internal: infer number of vertices from edge list (assumes 1-based ids)."""
_infer_n(edges::DataFrame) = max(maximum(edges.source), maximum(edges.target))

"""Internal: fast random edge index."""
@inline function _rand_edge_index(rng::AbstractRNG, m::Int)
    return rand(rng, 1:m)
end

"""Internal: build fixed mask of willing participants (percentage of population)."""
function _build_willing_mask(n::Int, willing::Real, rng::AbstractRNG)
    p = Float64(willing)
    0.0 <= p <= 1.0 || throw(ArgumentError("willing musí být v intervalu [0,1]"))

    k = round(Int, p * n)
    k = clamp(k, 0, n)

    mask = falses(n)
    if k > 0
        idx = randperm(rng, n)[1:k]
        @inbounds for i in idx
            mask[i] = true
        end
    end
    return mask
end

"""
    deffuant_willing_step!(x::AbstractVector{<:Real}, edges::DataFrame, willing_mask::AbstractVector{Bool};
        ϵ::Real,
        μ::Real=0.5,
        rng::AbstractRNG=Random.default_rng()
    ) -> Tuple{Int,Int,Bool}

Provede 1 náhodnou interakci v asymetrickém „willing“ modelu:
- vybere náhodnou hranu (u,v),
- náhodně zvolí směr interakce (kdo je learner a kdo partner),
- learner změní názor vždy (pokud patří mezi willing účastníky),
- partner změní názor pouze pokud |x[learner]-x[partner]| ≤ ϵ

Aktualizace learnera:
    x[learner] ← x[learner] + μ (x[partner] - x[learner])

Aktualizace partnera (jen v limitu ϵ):
    x[partner] ← x[partner] + μ (x[learner]_old - x[partner])

Vrací `(learner, partner, changed)`.

Pozn.: `x` typicky drž jako `Vector{Float64}`.
"""
function deffuant_willing_step!(
    x::AbstractVector{<:Real},
    edges::DataFrame;
    willing_mask::AbstractVector{Bool},
    ϵ::Real,
    μ::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
)
    μ = Float64(μ)
    ϵ = Float64(ϵ)
    0.0 < μ <= 0.5 || throw(ArgumentError("μ musí být v intervalu (0, 0.5]"))
    ϵ >= 0.0 || throw(ArgumentError("ϵ musí být >= 0"))

    m = nrow(edges)
    m > 0 || throw(ArgumentError("edges nesmí být prázdné"))
    length(willing_mask) == length(x) || throw(ArgumentError("willing_mask musí mít stejnou délku jako x"))

    i = _rand_edge_index(rng, m)
    u = Int(edges.source[i])
    v = Int(edges.target[i])
    u == v && return (u, v, false)

    learner, partner = rand(rng, Bool) ? (u, v) : (v, u)
    willing_mask[learner] || return (learner, partner, false)

    xl = Float64(x[learner])
    xp = Float64(x[partner])
    Δ = xp - xl

    # willing účastník se posune vždy směrem k partnerovi
    new_l = xl + μ * Δ
    x[learner] = new_l

    # partner se posune jen pokud je rozdíl v confidence bound
    changed_partner = false
    if abs(Δ) <= ϵ
        new_p = xp + μ * (xl - xp)
        x[partner] = new_p
        changed_partner = (new_p != xp)
    end

    changed_learner = (new_l != xl)
    return (learner, partner, changed_learner || changed_partner)
end

"""
    run_deffuant_willing(edges::DataFrame, x0::AbstractVector{<:Real};
        steps::Integer,
        ϵ::Real,
        μ::Real=0.5,
        rng::AbstractRNG=Random.default_rng(),
        n::Union{Int,Nothing}=nothing,
        clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0),
        snapshot_every::Integer=0,
        on_snapshot::Union{Nothing,Function}=nothing
    ) -> Vector{Float64}

Spustí Deffuantův model po `steps` interakcí a vrátí finální vektor názorů.

Parametry:
- `x0` je počáteční stav (délka n).
- `n`: pokud není zadáno, odvodí se z `edges` přes maximum id.
- `willing`: podíl uzlů v [0,1], které jsou „ochotné“ měnit svůj názor.
  Tento podíl se vybere fixně na začátku běhu.
- `clamp_range`: pokud není `nothing`, po každém kroku ořízne hodnoty do [low, high].
- `snapshot_every`: pokud > 0, volá `on_snapshot(step, x)` pro step=0 a pak každých `snapshot_every` kroků.
- `on_snapshot`: callback, např. ` (step,x)->save_opinion_snapshot("out/op", step, x)`.

Pozn.: `steps` jsou *interakce* (tj. náhodné vybrané hrany), ne „synchronní iterace přes všechny hrany“.
"""
function run_deffuant_willing(
    edges::DataFrame,
    x0::AbstractVector{<:Real};
    steps::Integer,
    ϵ::Real,
    μ::Real = 0.5,
    willing::Real = 0.1,
    rng::AbstractRNG = Random.default_rng(),
    n::Union{Int,Nothing} = nothing,
    clamp_range::Union{Nothing,Tuple{Real,Real}} = (0.0, 1.0),
    snapshot_every::Integer = 0,
    on_snapshot::Union{Nothing,Function} = nothing,
)
    steps = Int(steps)
    steps >= 0 || throw(ArgumentError("steps musí být >= 0"))

    n === nothing && (n = _infer_n(edges))
    n = Int(n)

    length(x0) == n || throw(ArgumentError("Délka x0 ($(length(x0))) musí být rovna n ($n)"))

    x = Float64.(x0)
    willing_mask = _build_willing_mask(n, willing, rng)

    # snapshot step 0
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

    for s in 1:steps
        deffuant_willing_step!(x, edges; willing_mask=willing_mask, ϵ=ϵ, μ=μ, rng=rng)
        if do_clamp
            @inbounds for i in 1:n
                xi = x[i]
                x[i] = xi < low ? low : (xi > high ? high : xi)
            end
        end

        if snapshot_every > 0 && on_snapshot !== nothing
            if s % snapshot_every == 0
                on_snapshot(s, x)
            end
        end
    end

    return x
end

"""
    run_deffuant_willing_series(edges::DataFrame, x0::AbstractVector{<:Real};
        steps::Integer,
        ϵ::Real,
        μ::Real=0.5,
        rng::AbstractRNG=Random.default_rng(),
        n::Union{Int,Nothing}=nothing,
        clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0),
        snapshots::Integer=101
    ) -> Vector{Vector{Float64}}

Spustí Deffuantův model a vrátí přímo časovou řadu stavů (včetně počátečního).

V této variantě `willing` značí podíl uzlů, které v interakci mohou měnit
svůj názor vždy; druhý uzel v páru se změní jen pokud je v limitu `ϵ`.

- `snapshots` říká kolik stavů chceš uložit včetně t=0 a t=steps.
  Typicky 101 => 0 + 99 mezistavů + finální.

Interně ukládá rovnoměrně rozložené snapshoty v čase.
"""
function run_deffuant_willing_series(
    edges::DataFrame,
    x0::AbstractVector{<:Real};
    steps::Integer,
    ϵ::Real,
    μ::Real = 0.5,
    willing::Real = 0.1,
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
    willing_mask = _build_willing_mask(n, willing, rng)
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

    # kroky, kdy snapshotovat (včetně finále)
    snap_steps = round.(Int, range(0, steps; length=snapshots))
    next_idx = 2
    next_snap = snap_steps[next_idx]

    for s in 1:steps
        deffuant_willing_step!(x, edges; willing_mask=willing_mask, ϵ=ϵ, μ=μ, rng=rng)
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

    # když steps == 0, range vrátí samé nuly a série už má 1 stav; doplníme finále
    while length(series) < snapshots
        push!(series, copy(x))
    end

    return series
end

end # module
