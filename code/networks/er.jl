module ER

using Random
using Graphs

export generateER

"""
    generateER(N::Integer, p::Real; rng=Random.default_rng()) -> SimpleGraph

Vygeneruje neorientovaný Erdős-Rényi graf `G(N, p)` bez self-loopů a duplicit.
- Uzly jsou 1..N (1-based).

Vrací `SimpleGraph` s vrcholy `1..N`.
"""
function generateER(N::Integer, p::Real; rng::AbstractRNG = Random.default_rng())
    N = Int(N)
    p = Float64(p)

    N >= 2 || throw(ArgumentError("N musí být >= 2"))
    0.0 <= p <= 1.0 || throw(ArgumentError("p musí být v intervalu [0,1]"))

    g = SimpleGraph(N)

    @inbounds for u in 1:(N-1)
        for v in (u+1):N
            if rand(rng) < p
                add_edge!(g, u, v)
            end
        end
    end

    return g
end

end # module
