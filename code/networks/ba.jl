module BA

export generateBA

using Random
using Graphs

"""
    generateBA(N::Integer, M::Integer; rng=Random.default_rng()) -> SimpleGraph

Vygeneruje neorientovanou Barabási–Albert síť G(N, M) bez self-loopů.
- Začíná z úplného grafu na (M+1) uzlech.
- Každý nový uzel připojí k M různým existujícím uzlům preferenčním připojováním (pravděpodobnost ∝ stupni).

Vrací `SimpleGraph` s 1-based indexy uzlů.
"""
function generateBA(N::Integer, M::Integer; rng::AbstractRNG = Random.default_rng())
    N = Int(N)
    M = Int(M)
    N >= 2 || throw(ArgumentError("N musí být >= 2"))
    M >= 1 || throw(ArgumentError("M musí být >= 1"))
    M < N  || throw(ArgumentError("M musí být < N"))

    # `m0` je velikost počátečního úplného grafu, ze kterého BA startuje.
    m0 = M + 1
    m0 <= N || throw(ArgumentError("Pro BA musí platit N >= M+1"))

    g = SimpleGraph(N)

    # Nejprve postavíme úplný graf na prvních `m0` vrcholech.
    degrees = zeros(Int, N)
    @inbounds for i in 1:(m0 - 1)
        for j in (i + 1):m0
            add_edge!(g, i, j)
            degrees[i] += 1
            degrees[j] += 1
        end
    end

    # "Urna" reprezentuje preferenční připojování: uzel se v ní vyskytuje
    # tolikrát, kolik má aktuálně incidentních hran.
    urn = Int[]
    @inbounds for u in 1:m0
        for _ in 1:degrees[u]
            push!(urn, u)
        end
    end
    isempty(urn) && throw(ArgumentError("Interní chyba: prázdná urna"))

    # Každý další uzel se připojí k `M` různým existujícím uzlům.
    for v in (m0 + 1):N
        chosen = Set{Int}()
        # Cíle vybíráme z urny, takže pravděpodobnost je úměrná stupni.
        while length(chosen) < M
            u = urn[rand(rng, 1:length(urn))]
            u != v && push!(chosen, u)
        end

        # Po přidání hrany musíme aktualizovat stupně i obsah urny.
        @inbounds for u in chosen
            add_edge!(g, u, v)

            degrees[u] += 1
            degrees[v] += 1

            # Každá nová hrana zvýší stupeň obou konců o 1, takže oba vrcholy
            # přidáme do urny právě jednou.
            push!(urn, u)
            push!(urn, v)
        end
    end

    return g
end

end # module
