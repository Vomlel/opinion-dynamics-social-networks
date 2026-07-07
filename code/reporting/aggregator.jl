

module Aggregator

using DataFrames
using Statistics

export aggregate_static,
       aggregate_dynamic

"""
    _safe_quantile(v, p) -> Float64

Vrátí kvantil pro neprázdný vektor. Prázdný vstup mapuje na `NaN`, aby se
agregační pipeline nerozbila na chybějících nebo neplatných hodnotách.
"""
_safe_quantile(v::AbstractVector{<:Real}, p::Real) = isempty(v) ? NaN : quantile(v, p)

"""
    aggregate_static(static_rows::DataFrame) -> DataFrame

Agreguje statické charakteristiky přes běhy (runs) pro každou kombinaci
(network, model, epsilon, mu).

Očekávané sloupce ve vstupu:
- `network::String`
- `model::String`
- `epsilon::Float64`
- `mu::Float64`
- `run::Int`
- `avg_degree::Float64`
- `density::Float64`
- `apl::Float64` (může obsahovat NaN/Inf; agregace je robustní)

Výstup obsahuje agregace přes běhy pro každou parametrickou kombinaci.
"""
function aggregate_static(static_rows::DataFrame)
    required = (:network, :model, :epsilon, :mu, :run, :avg_degree, :density, :apl)
    for c in required
        c in propertynames(static_rows) || throw(ArgumentError("aggregate_static: chybí sloupec $(c)"))
    end

    g = groupby(static_rows, [:network, :model, :epsilon, :mu])

    return combine(g) do sdf
        apl_vals = Float64.(sdf.apl)
        apl_vals = apl_vals[isfinite.(apl_vals)]

        DataFrame(
            runs = length(unique(sdf.run)),
            avg_degree_mean = mean(sdf.avg_degree),
            avg_degree_std  = std(sdf.avg_degree),
            density_mean    = mean(sdf.density),
            density_std     = std(sdf.density),
            apl_mean        = isempty(apl_vals) ? NaN : mean(apl_vals),
            apl_std         = isempty(apl_vals) ? NaN : std(apl_vals),
        )
    end
end

"""
    aggregate_dynamic(dynamic_rows::DataFrame) -> DataFrame

Agreguje dynamické charakteristiky přes běhy (runs) pro každou kombinaci
(network, model, epsilon, mu, t).

Očekávané sloupce ve vstupu:
- `network::String`
- `model::String`
- `epsilon::Float64`
- `mu::Float64`
- `run::Int`
- `t::Int`
- `clusters::Int`
- `pad::Float64`
- `polarization::Float64`

Výstup obsahuje agregace pro každý časový krok `t`, včetně kvantilů.
"""
function aggregate_dynamic(dynamic_rows::DataFrame)
    required = (:network, :model, :epsilon, :mu, :run, :t, :clusters, :pad, :polarization)
    for c in required
        c in propertynames(dynamic_rows) || throw(ArgumentError("aggregate_dynamic: chybí sloupec $(c)"))
    end

    g = groupby(dynamic_rows, [:network, :model, :epsilon, :mu, :t])

    return combine(g) do ddf
        pad_vals = Float64.(ddf.pad)
        pol_vals = Float64.(ddf.polarization)

        # Neplatné hodnoty odfiltrujeme, aby statistiky a kvantily nepadaly.
        pad_vals = sort(pad_vals[isfinite.(pad_vals)])
        pol_vals = sort(pol_vals[isfinite.(pol_vals)])

        clusters_f = sort(Float64.(ddf.clusters))

        DataFrame(
            runs = length(unique(ddf.run)),
            clusters_mean = mean(clusters_f),
            clusters_std  = std(clusters_f),
            pad_mean      = isempty(pad_vals) ? NaN : mean(pad_vals),
            pad_std       = isempty(pad_vals) ? NaN : std(pad_vals),
            polarization_mean = isempty(pol_vals) ? NaN : mean(pol_vals),
            polarization_std  = isempty(pol_vals) ? NaN : std(pol_vals),
            clusters_q25 = quantile(clusters_f, 0.25),
            clusters_q50 = quantile(clusters_f, 0.5),
            clusters_q75 = quantile(clusters_f, 0.75),
            pad_q25 = _safe_quantile(pad_vals, 0.25),
            pad_q50 = _safe_quantile(pad_vals, 0.5),
            pad_q75 = _safe_quantile(pad_vals, 0.75),
            polarization_q25 = _safe_quantile(pol_vals, 0.25),
            polarization_q50 = _safe_quantile(pol_vals, 0.5),
            polarization_q75 = _safe_quantile(pol_vals, 0.75),
        )
    end
end

end # module
