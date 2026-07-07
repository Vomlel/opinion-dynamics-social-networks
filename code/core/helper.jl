module Helper
export parse_args, ensure_dirs, dynamic_metrics_series, static_metrics, opinion_vector, edge_dataframe

using Random
using DataFrames
using CSV
using Graphs

using ..StaticCharacteristics
using ..DynamicCharacteristics
using ..GenerateOpinion

"""
    parse_args(args::Vector{String}) -> String

Zpracuje argumenty příkazové řádky a vrátí cestu ke konfiguračnímu souboru.

Podporovaný formát je `-c path/to/config.toml`.
"""
function parse_args(args::Vector{String})
    config_path = "config.toml"
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "-c"
            i += 1
            i <= length(args) || error("Za přepínačem -c chybí cesta ke konfiguraci")
            config_path = args[i]
        else
            error("Neznámý argument: $a")
        end
        i += 1
    end
    return config_path
end

"""
    ensure_dirs(out_dir::String)

Vytvoří výstupní adresáře používané během experimentu.
"""
function ensure_dirs(out_dir::String)
    mkpath(out_dir)
    mkpath(joinpath(out_dir, "data"))
    mkpath(joinpath(out_dir, "aggregated"))
    mkpath(joinpath(out_dir, "plots"))
end

"""
    dynamic_metrics_series(series::Vector{Vector{Float64}}) -> DataFrame

Spočítá dynamické charakteristiky pro všechny snapshoty v časové řadě.

Každý snapshot se před výpočtem seřadí, protože downstream metriky
očekávají seřazený vektor názorů.
"""
function dynamic_metrics_series(series::Vector{Vector{Float64}})
    T = length(series)
    t = collect(0:(T-1))

    clusters = Vector{Int}(undef, T)
    pad = Vector{Float64}(undef, T)
    pol = Vector{Float64}(undef, T)

    @inbounds for i in 1:T
        x = sort(series[i])
        clusters[i] = count_of_clusters_dbscan_auto(x)
        pad[i] = pairwise_average_distance(x)
        pol[i] = polarization_der(x)
    end

    return DataFrame(t=t, clusters=clusters, pad=pad, polarization=pol)
end

"""
    edge_dataframe(g::AbstractGraph) -> DataFrame

Převede graf na edge list se sloupci `source` a `target`.

Tato funkce slouží hlavně pro export do CSV; uvnitř simulací pracujeme přímo
nad `Graphs.jl` reprezentací.
"""
function edge_dataframe(g::AbstractGraph)
    src_nodes = Int[]
    dst_nodes = Int[]

    sizehint!(src_nodes, ne(g))
    sizehint!(dst_nodes, ne(g))

    for e in edges(g)
        push!(src_nodes, src(e))
        push!(dst_nodes, dst(e))
    end

    return DataFrame(source=src_nodes, target=dst_nodes)
end

"""
    static_metrics(g::AbstractGraph) -> NamedTuple

Spočítá základní statické charakteristiky grafu.
"""
function static_metrics(g::AbstractGraph)
    avgk = StaticCharacteristics.average_degree(g)
    rho  = StaticCharacteristics.density(g)
    apl  = StaticCharacteristics.average_shortest_path_length(g)
    return (avg_degree=avgk, density=rho, apl=apl)
end

"""
    opinion_vector(cfg::OpinionsConfig, g::AbstractGraph; seed::Int) -> Vector{Float64}

Vygeneruje počáteční názory podle konfigurace experimentu.
"""
function opinion_vector(cfg, g::AbstractGraph; seed::Int)
    model = Symbol(cfg.model)
    dist  = Symbol(cfg.distribution)

    return generate_opinions(
        g;
        model=model,
        distribution=dist,
        low=cfg.low,
        high=cfg.high,
        μ=cfg.mu,
        σ=cfg.sigma,
        polarized_ratio=cfg.polarized_ratio,
        eps=cfg.polarized_eps,
        seed=seed,
    )
end

end # module
