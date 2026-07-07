module HegselmannKrause

using Random
using Graphs

export hk_step!,
       run_hk,
       run_hk_series

"""
    hk_step!(x::AbstractVector{<:Real}, g::AbstractGraph;
        ϵ::Real,
        rng::AbstractRNG=Random.default_rng(),
        include_self::Bool=true,
        clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0)
    ) -> Vector{Float64}

Provede 1 synchronní iteraci Hegselmann–Krause (HK) modelu na dané síti.

Definice (síťová varianta):
- pro každý uzel i se vezmou jeho sousedé v grafu (a volitelně i i samo)
- do množiny Nᵋ(i) se zahrnou pouze ti j, pro které |x[i] - x[j]| ≤ ϵ
- aktualizace je průměr:

    xᵢ(t+1) = average( xⱼ(t) for j ∈ Nᵋ(i) )

Poznámky:
- HK je synchronní: používá se původní x(t) pro výpočet všech x(t+1).
- Pokud Nᵋ(i) vyjde prázdné (např. include_self=false a žádný soused nevyhoví),
  nechá se xᵢ beze změny.
- Funkce vrací nově vytvořený vektor `x_next` a zároveň jím přepíše původní `x`.
"""
function hk_step!(
    x::AbstractVector{<:Real},
    g::AbstractGraph;
    ϵ::Real,
    rng::AbstractRNG = Random.default_rng(),
    include_self::Bool = true,
    clamp_range::Union{Nothing,Tuple{Real,Real}} = (0.0, 1.0),
)
    ϵ = Float64(ϵ)
    ϵ >= 0.0 || throw(ArgumentError("ϵ musí být >= 0"))

    n = nv(g)
    length(x) == n || throw(ArgumentError("Délka x ($(length(x))) musí být rovna počtu vrcholů v grafu ($n)"))
    n == 0 && return Float64[]

    # HK je synchronní, takže všechny uzly musí vycházet ze stejného stavu x(t).
    x0 = Float64.(x)
    x_next = similar(x0)

    do_clamp = clamp_range !== nothing
    low = 0.0
    high = 1.0
    if do_clamp
        low = Float64(clamp_range[1])
        high = Float64(clamp_range[2])
        low < high || throw(ArgumentError("clamp_range musí mít low < high"))
    end

    @inbounds for i in 1:n
        xi = x0[i]
        sumv = 0.0
        cnt = 0

        if include_self
            sumv += xi
            cnt += 1
        end

        for j in neighbors(g, i)
            xj = x0[j]
            if abs(xi - xj) <= ϵ
                sumv += xj
                cnt += 1
            end
        end

        newxi = (cnt > 0) ? (sumv / cnt) : xi
        if do_clamp
            newxi = newxi < low ? low : (newxi > high ? high : newxi)
        end
        x_next[i] = newxi
    end

    @inbounds for i in 1:n
        x[i] = x_next[i]
    end

    return x_next
end

"""
    run_hk(g::AbstractGraph, x0::AbstractVector{<:Real};
        iterations::Integer,
        ϵ::Real,
        rng::AbstractRNG=Random.default_rng(),
        include_self::Bool=true,
        clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0),
        snapshot_every::Integer=0,
        on_snapshot::Union{Nothing,Function}=nothing
    ) -> Vector{Float64}

Spustí HK model po `iterations` synchronních iterací a vrátí finální názory.

Snapshoty:
- pokud `snapshot_every > 0` a `on_snapshot` není `nothing`, volá `on_snapshot(step, x)` pro step=0
  a pak každých `snapshot_every` iterací.
"""
function run_hk(
    g::AbstractGraph,
    x0::AbstractVector{<:Real};
    iterations::Integer,
    ϵ::Real,
    rng::AbstractRNG = Random.default_rng(),
    include_self::Bool = true,
    clamp_range::Union{Nothing,Tuple{Real,Real}} = (0.0, 1.0),
    snapshot_every::Integer = 0,
    on_snapshot::Union{Nothing,Function} = nothing,
)
    iterations = Int(iterations)
    iterations >= 0 || throw(ArgumentError("iterations musí být >= 0"))

    n = nv(g)
    length(x0) == n || throw(ArgumentError("Délka x0 ($(length(x0))) musí být rovna počtu vrcholů v grafu ($n)"))

    x = Float64.(x0)

    if snapshot_every > 0 && on_snapshot !== nothing
        on_snapshot(0, x)
    end

    for t in 1:iterations
        hk_step!(x, g; ϵ=ϵ, rng=rng, include_self=include_self, clamp_range=clamp_range)

        if snapshot_every > 0 && on_snapshot !== nothing
            if t % snapshot_every == 0
                on_snapshot(t, x)
            end
        end
    end

    return x
end

"""
    run_hk_series(g::AbstractGraph, x0::AbstractVector{<:Real};
        iterations::Integer,
        ϵ::Real,
        rng::AbstractRNG=Random.default_rng(),
        include_self::Bool=true,
        clamp_range::Union{Nothing,Tuple{Real,Real}}=(0.0, 1.0),
        snapshots::Integer=101
    ) -> Vector{Vector{Float64}}

Spustí HK model a vrátí časovou řadu stavů (včetně počátečního).

- `snapshots` je počet uložených stavů včetně t=0 a t=iterations.
- Snapshoty jsou rovnoměrně rozložené v čase.
"""
function run_hk_series(
    g::AbstractGraph,
    x0::AbstractVector{<:Real};
    iterations::Integer,
    ϵ::Real,
    rng::AbstractRNG = Random.default_rng(),
    include_self::Bool = true,
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

    # Snapshoty rozložíme do času předem, aby výstup měl stabilní délku.
    snap_steps = round.(Int, range(0, iterations; length=snapshots))
    next_idx = 2
    next_snap = snap_steps[next_idx]

    for t in 1:iterations
        hk_step!(x, g; ϵ=ϵ, rng=rng, include_self=include_self, clamp_range=clamp_range)

        if t == next_snap
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
