module DynamicCharacteristics

using Statistics

export count_of_clusters_auto, avg, polarization, pairwise_average_distance, pad_inner, pad_outer, group_values, count_of_group_clusters, polarization_der, count_of_clusters_dbscan_auto

#=
Definuje 3 dynymické veličiny pro pole názorů.
Vstupní pole musí být seřazené.
PAD, Polarization a Count of clusters
=#

#=
Průměrná hodnota názoru
=#
function avg(x::AbstractVector{<:Real})
    n = length(x)
    s = 0.0
    for value in x
        s += value
    end
    return s / n
end


#=
Průměrná vzdálenost dvou názorů od sebe.
Součet přes všechny dvojice |x[i]-x[j]| / počet dvojic.
Ale lze počítat optimalizovaně.
Vstupní pole musí být seřazené!!!
=#
function pairwise_average_distance(x::AbstractVector{<:Real})
    n = length(x)
    s = 0.0
    for k in 1:n
        s += x[k] * (2k - n - 1)
    end
    return 2 * s / n / (n - 1)
end


#=
Rozdělíme názory na dvě poloviny, horní a dolní.
Spočteme průměrnou vzdálenost dvou názorů od sebe, ale jen pro dvojice uvnitř skupin.
Součet přes všechny dvojice |x[i]-x[j]| / počet dvojic, jen pro dvojice, kde jsou oba prvky ve stejné skupině.
Ale lze počítat optimalizovaně.
Vstupní pole musí být seřazené!!!
=#
function pad_inner(x::AbstractVector{<:Real})
    n = length(x)
    half = n ÷ 2
    lower = view(x, 1:half)
    upper = view(x, n-half+1:n)
    return (pairwise_average_distance(lower)+pairwise_average_distance(upper)) / 2
end


#=
Rozdělíme názory na dvě poloviny, horní a dolní.
Spočteme průměrnou vzdálenost dvou názorů od sebe, ale jen pro dvojice mezi skupinami.
Součet přes všechny dvojice |x[i]-x[j]| / počet dvojic, jen pro dvojice, kde je každý prvek v jiné skupině.
Ale lze počítat optimalizovaně.
Vstupní pole musí být seřazené!!!
=#
function pad_outer(x::AbstractVector{<:Real})
    n = length(x)
    half = n ÷ 2
    lower = view(x, 1:half)
    upper = view(x, n-half+1:n)
    return avg(upper)-avg(lower)
end


#=
Polorizaci chápeme jako pad_outer/pad_inner.
Kolikrát jsou si názory mezi oběma skupinami dál než názory uvnitř skupin.
Vstupní pole musí být seřazené!!!
=#
function polarization(x::AbstractVector{<:Real})
    return pad_outer(x) / pad_inner(x)
end

"""
    polarization_der(x::AbstractVector{<:Real};
        α::Real = 1.3,
        gridsize::Integer = 256,
        bandwidth::Union{Nothing,Real} = nothing,
        range_min::Real = 0.0,
        range_max::Real = 1.0
    ) -> Float64

Numerická aproximace spojitého Duclos–Esteban–Ray (2004) indexu polarizace:

    P_α(f) = ∬ f(x)^(1+α) f(y) |y-x| dy dx

kde `f` je hustota odhadnutá z dat `x` pomocí gaussovského KDE.

Parametry:
- `α`: parametr polarizační citlivosti, obvykle v intervalu (0, 1.6]
- `gridsize`: počet bodů mřížky pro numerickou integraci
- `bandwidth`: šířka jádra pro KDE; pokud je `nothing`, použije se Silvermanovo pravidlo
- `range_min`, `range_max`: interval názorové škály

Poznámka:
- Výstup není automaticky normalizovaný do [0,1].
- Hodnoty mají smysl hlavně pro porovnávání mezi běhy při stejné metodice.
ZDROJ: https://osf.io/preprints/socarxiv/e5vp8_v1 - kapitola 6.2.11.
"""

function polarization_der(
    x::AbstractVector{<:Real};
    α::Real = 1.3,
    gridsize::Integer = 256,
    bandwidth::Union{Nothing,Real} = nothing,
    range_min::Real = 0.0,
    range_max::Real = 1.0,
)
    n = length(x)
    n == 0 && throw(ArgumentError("x nesmí být prázdné"))
    α = Float64(α)
    0.0 < α <= 1.6 || throw(ArgumentError("α musí být v intervalu (0, 1.6]"))

    a = Float64(range_min)
    b = Float64(range_max)
    a < b || throw(ArgumentError("Musí platit range_min < range_max"))

    xs = Float64.(x)
    any(isnan, xs) && throw(ArgumentError("x nesmí obsahovat NaN"))
    any(v -> v < a || v > b, xs) &&
        throw(ArgumentError("Všechny názory musí ležet v intervalu [$a, $b]"))

    # bandwidth: Silvermanovo pravidlo
    h = if bandwidth === nothing
        s = std(xs)
        iqr = quantile(xs, 0.75) - quantile(xs, 0.25)
        σ = min(s, iqr / 1.34)
        σ = σ == 0.0 ? max((b - a) / 100, eps()) : σ
        0.9 * σ * n^(-1/5)
    else
        Float64(bandwidth)
    end
    h > 0 || throw(ArgumentError("bandwidth musí být > 0"))

    m = Int(gridsize)
    m >= 2 || throw(ArgumentError("gridsize musí být >= 2"))

    grid = collect(range(a, b; length=m))
    Δ = (b - a) / (m - 1)

    # Gaussovské KDE
    invh = 1.0 / h
    normconst = 1.0 / (n * h * sqrt(2π))
    f = Vector{Float64}(undef, m)

    @inbounds for i in 1:m
        gx = grid[i]
        s = 0.0
        for v in xs
            z = (gx - v) * invh
            s += exp(-0.5 * z * z)
        end
        f[i] = normconst * s
    end

    # Renormalizace na [a,b]
    mass = sum(f) * Δ
    mass > 0 || return 0.0
    @inbounds for i in 1:m
        f[i] /= mass
    end

    # Prefixové sumy pro rychlý výpočet:
    # row_i = sum_j f[j] * |grid[j] - grid[i]|
    pref_f = Vector{Float64}(undef, m)
    pref_fx = Vector{Float64}(undef, m)

    s_f = 0.0
    s_fx = 0.0
    @inbounds for i in 1:m
        s_f += f[i]
        s_fx += f[i] * grid[i]
        pref_f[i] = s_f
        pref_fx[i] = s_fx
    end

    total_f = pref_f[m]
    total_fx = pref_fx[m]

    total = 0.0
    @inbounds for i in 1:m
        xi = grid[i]

        left_f  = i > 1 ? pref_f[i - 1]  : 0.0
        left_fx = i > 1 ? pref_fx[i - 1] : 0.0

        right_f  = total_f - pref_f[i]
        right_fx = total_fx - pref_fx[i]

        # sum_{j<i} f[j] * (xi - xj) + sum_{j>i} f[j] * (xj - xi)
        row = (xi * left_f - left_fx) + (right_fx - xi * right_f)

        total += f[i]^(1 + α) * row
    end

    return total * Δ * Δ
end

"""
Clusterovani
zdroj: https://medium.com/data-science/the-5-clustering-algorithms-data-scientists-need-to-know-a36d136ef68
kapitola: Density-Based Spatial Clustering of Applications with Noise (DBSCAN)
    dbscan_1d(x::AbstractVector{<:Real}; eps::Real, min_pts::Int=2)

DBSCAN pro jednorozměrná data.

Parametry:
- `x`: vstupní hodnoty
- `eps`: maximální vzdálenost dvou bodů, aby byli sousedé
- `min_pts`: minimální počet bodů v eps-okolí (včetně samotného bodu),
             aby byl bod považován za core point

Vrací:
- `labels::Vector{Int}`:
    - `-1` = noise
    - `1,2,3,...` = ID clusteru
"""
function dbscan_1d(x::AbstractVector{<:Real}; eps::Real, min_pts::Int=2)
    n = length(x)
    n == 0 && return Int[]
    eps >= 0 || throw(ArgumentError("eps musí být >= 0"))
    min_pts >= 1 || throw(ArgumentError("min_pts musí být >= 1"))

    xs = Float64.(x)

    # left[i]  = nejmenší index v eps-okolí bodu i
    # right[i] = největší index v eps-okolí bodu i
    left  = Vector{Int}(undef, n)
    right = Vector{Int}(undef, n)

    # Spočítání levých hranic pomocí sliding window
    l = 1
    @inbounds for i in 1:n
        while xs[i] - xs[l] > eps
            l += 1
        end
        left[i] = l
    end

    # Spočítání pravých hranic pomocí sliding window
    r = 1
    @inbounds for i in 1:n
        while r < n && xs[r + 1] - xs[i] <= eps
            r += 1
        end
        right[i] = r
    end

    # Core body
    core = Vector{Bool}(undef, n)
    @inbounds for i in 1:n
        core[i] = (right[i] - left[i] + 1) >= min_pts
    end

    labels = fill(-1, n)   # defaultně noise
    cluster_id = 0
    i = 1

    while i <= n
        if !core[i]
            i += 1
            continue
        end

        # Začátek nového clusteru
        cluster_id += 1
        cluster_start = left[i]
        cluster_end = right[i]

        # Rozšiřujeme cluster přes core body, které do něj ještě spadají
        j = i
        while j <= n && j <= cluster_end
            if core[j]
                cluster_end = max(cluster_end, right[j])
            end
            j += 1
        end

        # Označ celý interval clusteru
        @inbounds for k in cluster_start:cluster_end
            labels[k] = cluster_id
        end

        i = j
    end

    return labels
end

function count_of_clusters_dbscan_auto(x::AbstractVector{<:Real}; eps::Real=0.01, min_pts::Int=-1)
    isempty(x) && return 0
    if min_pts === -1
        min_pts = max(2, ceil(Int, length(x) / 100)) # zaokrouhleno na cela cisla nahoru; minimalne ale 2 velikost clusteru alespon 1% celkove populace
    end

    labels = dbscan_1d(x; eps=eps, min_pts=min_pts)

    cluster_count = maximum(labels)

    cluster_count <= 0 && return Int[]

    sizes = zeros(Int, cluster_count)
    @inbounds for label in labels
        if label > 0
            sizes[label] += 1
        end
    end
    return length(Set(l for l in labels if l > 0))
end


#=
Názory převede na skupiny, kde u každé skupiny je počet názorů, které do ní patří.
=#
function group_values(x::AbstractVector{<:Real}, from::Real, to::Real, count::Int)
    groups =  zeros(Int, count)
    for value in x
        if from <= value < to
            index = floor(Int, (value - from) / (to - from) * count) + 1
            groups[index] += 1
        elseif value == to
            groups[count] += 1
        end
    end
    return groups
end


#=
Spočítá počet clustrů.
Clustery musí být odděleny skupinou menší než separationLimit.
Cluster musí být alespoň sizeLimit velký. Velikost clusteru je součet velikostí jeho skupin.
=#
function count_of_group_clusters(groups::AbstractVector{<:Int}, separationLimit::Real, sizeLimit::Real)
    clusters = 0
    separationOccured = true
    sizeLimitOccured = false
    size = 0
    for g in groups
        if g >= separationLimit
            size += g
        end
        if separationOccured && !sizeLimitOccured && size >= sizeLimit
            clusters += 1
            sizeLimitOccured = true
            separationOccured = false
        end
        if g <= separationLimit
            separationOccured = true
            size = 0.0
            sizeLimitOccured = false
        end
    end

    return clusters
end


#=
Spočítá počet clustrů.
Rozdělí názory na skupiny a spočte počet clusterů.
=#
function count_of_clusters(x::AbstractVector{<:Real}, from::Real, to::Real, count::Int, separationLimit::Real, sizeLimit::Real)
    groups = group_values(x, from, to, count)
    #println("groups=$groups, sizeLimit=$sizeLimit, separationLimit=$separationLimit")
    return count_of_group_clusters(groups, separationLimit, sizeLimit)
end

#=
Spočítá počet clustrů.
Rozdělí názory na skupiny a spočte počet clusterů.
Parametry rozdělení sám nastaví "smysluplně"
Vstupní pole musí být seřazené!!!
=#
function count_of_clusters_auto(x::AbstractVector{<:Real})
    n = length(x)
    if n == 0
        return 0
    end

    #Od do si vezme z krajních názorů
    from = x[1]
    to = x[n]

    #Počet slupin tak abychom měli hodně skupin i hodně názorů ve skupinách
    count = floor(Int, sqrt(n))

    # Cluster musí mít alespoň n/count prvků, což je průměrný počet prvků ve skupině
    sizeLimit = n / count

    #Separační limit je o dost menší než průměrný počet ve skupině n^0.5
    separationLimit = sqrt(n/count)

    return count_of_clusters(x, from, to, count, separationLimit, sizeLimit)
end


end # module
