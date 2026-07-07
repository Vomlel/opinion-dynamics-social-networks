

module GenerateOpinion

using Random
using DataFrames
using CSV
using Statistics
using Printf
using Graphs

export generate_opinions,
       opinions_to_dataframe,
       save_opinions_csv,
       save_opinion_snapshot,
       load_opinion_snapshot,
       save_opinion_series_long_csv,
       save_opinion_series_wide_csv

"""
    generate_opinions(g::AbstractGraph;
        model::Symbol=:continuous,
        distribution::Symbol=:uniform,
        rng::AbstractRNG=Random.default_rng(),
        μ::Float64=0.0,
        σ::Float64=0.2,
        low::Float64=0.0,
        high::Float64=1.0,
        polarized_ratio::Float64=0.5,
        eps::Float64=0.05,
        seed::Union{Int,Nothing}=nothing,
    ) -> Vector{Float64}

Vygeneruje počáteční názory pro vrcholy grafu `g`.

Parametry:
- `model`:
  - `:continuous` => názory jsou reálná čísla v intervalu [`low`,`high`] (typicky [0,1]).
  - `:binary`     => názory jsou 0.0 nebo 1.0.
- `distribution` (pro `:continuous`):
  - `:uniform`    => U(low, high)
  - `:normal`     => N(μ, σ) a ořízne do [low, high]
  - `:polarized`  => směs dvou shluků u krajů: část u `low+eps`, část u `high-eps`
- `seed`: pokud zadáš, vytvoří se lokální `MersenneTwister(seed)` (přebije `rng`).

Vrací vektor `x::Vector{Float64}` délky `n`, kde `x[i]` je názor uzlu i.
"""
function generate_opinions(
    g::AbstractGraph;
    model::Symbol=:continuous,
    distribution::Symbol=:uniform,
    rng::AbstractRNG=Random.default_rng(),
    μ::Float64=0.0,
    σ::Float64=0.2,
    low::Float64=0.0,
    high::Float64=1.0,
    polarized_ratio::Float64=0.5,
    eps::Float64=0.05,
    seed::Union{Int,Nothing}=nothing,
)
    local_rng = seed === nothing ? rng : MersenneTwister(seed)

    n = nv(g)
    n >= 1 || throw(ArgumentError("n musí být >= 1"))
    low < high || throw(ArgumentError("low musí být < high"))

    x = Vector{Float64}(undef, n)

    if model == :binary
        @inbounds for i in 1:n
            x[i] = rand(local_rng) < 0.5 ? 0.0 : 1.0
        end
        return x
    elseif model != :continuous
        throw(ArgumentError("Neznámý model: $model (povoleno :continuous nebo :binary)"))
    end

    if distribution == :uniform
        span = high - low
        @inbounds for i in 1:n
            x[i] = low + span * rand(local_rng)
        end

    elseif distribution == :normal
        # Box-Muller transformace, aby nebyla potřeba další statistická knihovna.
        @inbounds for i in 1:n
            u1 = max(rand(local_rng), eps(Float64))
            u2 = rand(local_rng)
            z  = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2)
            val = μ + σ * z
            x[i] = clamp(val, low, high)
        end

    elseif distribution == :polarized
        0.0 <= polarized_ratio <= 1.0 || throw(ArgumentError("polarized_ratio musí být v [0,1]"))
        0.0 <= eps < (high-low)/2 || throw(ArgumentError("eps musí být v [0, (high-low)/2)"))

        left  = low + eps
        right = high - eps
        @inbounds for i in 1:n
            if rand(local_rng) < polarized_ratio
                # Levý pól s malým šumem kolem `low + eps`.
                x[i] = clamp(left + eps * (2rand(local_rng) - 1), low, high)
            else
                # Pravý pól s malým šumem kolem `high - eps`.
                x[i] = clamp(right + eps * (2rand(local_rng) - 1), low, high)
            end
        end

    else
        throw(ArgumentError("Neznámá distribuce: $distribution (povoleno :uniform, :normal, :polarized)"))
    end

    return x
end

"""
    opinions_to_dataframe(x::AbstractVector{<:Real}) -> DataFrame

Převede vektor názorů na tabulku s 1-based `node` a `opinion`.
"""
function opinions_to_dataframe(x::AbstractVector{<:Real})
    n = length(x)
    return DataFrame(node = collect(1:n), opinion = Float64.(x))
end

"""
    save_opinions_csv(path::AbstractString, x::AbstractVector{<:Real})

Uloží názory do CSV se sloupci `node,opinion`.
"""
function save_opinions_csv(path::AbstractString, x::AbstractVector{<:Real})
    df = opinions_to_dataframe(x)
    CSV.write(path, df)
    return nothing
end

"""
    save_opinion_snapshot(dir::AbstractString, step::Integer, x::AbstractVector{<:Real};
        prefix::AbstractString="opinions",
        width::Int=3,
        ext::AbstractString="csv"
    ) -> String

Uloží jeden stav názorů jako samostatný soubor do adresáře `dir`.

- `step` je index kroku (např. 0 = počáteční stav, 100 = finální stav).
- Vytváří soubor ve formátu: `{prefix}_{step_padded}.{ext}` (např. `opinions_000.csv`).
- CSV má sloupce `node,opinion`.

Vrací plnou cestu k uloženému souboru.
"""
function save_opinion_snapshot(
    dir::AbstractString,
    step::Integer,
    x::AbstractVector{<:Real};
    prefix::AbstractString = "opinions",
    width::Int = 3,
    ext::AbstractString = "csv",
)
    mkpath(dir)
    filename = @sprintf("%s_%0*d.%s", prefix, width, Int(step), ext)
    path = joinpath(dir, filename)
    save_opinions_csv(path, x)
    return path
end

"""
    load_opinion_snapshot(path::AbstractString) -> Vector{Float64}

Načte uložený snapshot z CSV (`node,opinion`) a vrátí vektor názorů.

Pozn.: Předpokládá, že `node` jsou 1..n. Pokud tam budou díry, vyhodí chybu.
"""
function load_opinion_snapshot(path::AbstractString)
    df = CSV.read(path, DataFrame)
    (:node in propertynames(df) && :opinion in propertynames(df)) ||
        throw(ArgumentError("CSV musí mít sloupce 'node' a 'opinion'"))

    nodes = Int.(df.node)
    ops = Float64.(df.opinion)

    n = length(nodes)
    n == 0 && return Float64[]

    # Seřadíme řádky podle indexu uzlu a ověříme souvislou indexaci 1..n.
    ord = sortperm(nodes)
    nodes = nodes[ord]
    ops = ops[ord]
    nodes[1] == 1 || throw(ArgumentError("Node index musí začínat na 1"))
    @inbounds for i in 1:n
        nodes[i] == i || throw(ArgumentError("Node indexy musí být 1..n bez děr (chyba u $i)"))
    end

    return ops
end

"""
    save_opinion_series_long_csv(path::AbstractString, series::Vector{<:AbstractVector{<:Real}})

Uloží celou časovou řadu názorů do jednoho CSV v "long" tvaru:

    step,node,opinion
    0,1,0.12
    0,2,0.77
    ...
    100,1,0.34

Výhoda: snadná analýza v pandas/R/SQL, dobře škáluje i pro velké T.
"""
function save_opinion_series_long_csv(path::AbstractString, series::Vector{<:AbstractVector{<:Real}})
    T = length(series)
    T == 0 && throw(ArgumentError("series nesmí být prázdné"))

    n = length(series[1])
    @inbounds for t in 1:T
        length(series[t]) == n || throw(ArgumentError("Všechny stavy musí mít stejnou délku (n)"))
    end

    steps = Int[]
    nodes = Int[]
    opinions = Float64[]

    sizehint!(steps, T * n)
    sizehint!(nodes, T * n)
    sizehint!(opinions, T * n)

    @inbounds for t in 1:T
        x = series[t]
        # V exportu držíme kroky 0-based, takže 0 značí počáteční stav.
        step = t - 1
        for i in 1:n
            push!(steps, step)
            push!(nodes, i)
            push!(opinions, Float64(x[i]))
        end
    end

    df = DataFrame(step = steps, node = nodes, opinion = opinions)
    CSV.write(path, df)
    return nothing
end

"""
    save_opinion_series_wide_csv(path::AbstractString, series::Vector{<:AbstractVector{<:Real}};
        prefix::AbstractString="step"
    )

Uloží časovou řadu do jednoho CSV v "wide" tvaru:

    node,step0,step1,...,step100

Výhoda: jednoduché pro heatmapy / matice. Nevýhoda: hodně sloupců.
"""
function save_opinion_series_wide_csv(
    path::AbstractString,
    series::Vector{<:AbstractVector{<:Real}};
    prefix::AbstractString = "step",
)
    T = length(series)
    T == 0 && throw(ArgumentError("series nesmí být prázdné"))

    n = length(series[1])
    @inbounds for t in 1:T
        length(series[t]) == n || throw(ArgumentError("Všechny stavy musí mít stejnou délku (n)"))
    end

    df = DataFrame(node = collect(1:n))
    for t in 1:T
        col = Symbol("$(prefix)$(t-1)")
        df[!, col] = Float64.(series[t])
    end

    CSV.write(path, df)
    return nothing
end

end # module
