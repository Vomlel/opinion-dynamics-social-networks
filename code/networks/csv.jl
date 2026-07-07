module CSVGraph

using CSV
using DataFrames
using Graphs

export loadCSVGraph, loadCSVGraph2

"""
    loadCSVGraph(path::AbstractString) -> DataFrame

Načte hrany ze souboru CSV ve tvaru:

source,target
0,1972
0,5111
...

Vrací edge list jako `DataFrame(source, target)` s 1-based indexací, bez
self-loopů a bez duplicit.
"""
function loadCSVGraph(path::AbstractString)
    df = CSV.read(path, DataFrame)

    # Očekáváme standardní dvojici sloupců se zdrojem a cílem hrany.
    (:source in propertynames(df) && :target in propertynames(df)) ||
        throw(ArgumentError("CSV musí mít sloupce 'source' a 'target'"))

    # CSV může vrátit obecnější numerické typy, proto vše převedeme na `Int`.
    df.source = Int.(df.source)
    df.target = Int.(df.target)

    # Vstupní data jsou 0-based, interně v projektu používáme 1-based indexaci.
    df.source .+= 1
    df.target .+= 1

    # Pro jednoduchý neorientovaný graf nechceme self-loop hrany.
    df = df[df.source .!= df.target, :]

    # Hrany normalizujeme do tvaru `u < v`, aby šly dobře deduplikovat.
    u = similar(df.source)
    v = similar(df.target)
    @inbounds for i in eachindex(df.source)
        a = df.source[i]
        b = df.target[i]
        if a < b
            u[i] = a
            v[i] = b
        else
            u[i] = b
            v[i] = a
        end
    end
    df = DataFrame(source=u, target=v)

    # Paralelní nebo obousměrný export může obsahovat duplicitní hrany.
    unique!(df)

    return df
end

"""
    loadCSVGraph2(path::AbstractString) -> SimpleGraph

Načte graf z adresáře s `edges.csv` a `nodes.csv` a vrátí ho jako
`SimpleGraph`.

Předpokládá, že `reference_id` v `nodes.csv` po převodu na 1-based indexaci
tvoří přesně množinu `1:N`.
"""
function loadCSVGraph2(path::AbstractString)
    csv_edges = CSV.read(joinpath(path, "edges.csv"), DataFrame)
    csv_nodes = CSV.read(joinpath(path, "nodes.csv"), DataFrame)
    nodes = csv_nodes.reference_id

    csv_edges.source .+= 1
    csv_edges.target .+= 1
    nodes .+= 1
    N = length(nodes)
    E = length(csv_edges.source)

    _check_indexes_1_to_n(nodes)

    g = SimpleGraph(N)
    for i in 1:E
        add_edge!(g, csv_edges.source[i], csv_edges.target[i])
    end
    return g
end

"""
    _check_indexes_1_to_n(nodes::Vector{Int}) -> Bool

Ověří, že identifikátory uzlů tvoří souvislý interval `1:n`.
"""
function _check_indexes_1_to_n(nodes::Vector{Int})
    n = length(nodes)
    if sort(nodes) != collect(1:n)
        error("Reference IDs nejsou 1..N bez mezer")
    end
    return true
end
end # module
