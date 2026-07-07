module ThesisFigures

ENV["GKSwstype"] = "100"

import Plots

using CSV
using DataFrames
using Statistics
using Base.Filesystem: dirname, mkpath
using Plots

export generate_thesis_figures

const METRICS = [:clusters, :pad, :polarization]
const MODEL_ORDER = ["deffuant", "hk", "deffuant_hrdcr", "deffuant_anchored", "deffuant_partisan"]
const NETWORK_ORDER = ["er", "ba", "real"]

const MODEL_LABELS = Dict(
    "deffuant" => "Deffuantův model",
    "hk" => "HK",
    "deffuant_hrdcr" => "HRDCR",
    "deffuant_anchored" => "Anchored Deffuant",
    "deffuant_partisan" => "Partisan Deffuant",
)

const NETWORK_LABELS = Dict(
    "er" => "ER",
    "ba" => "BA",
    "real" => "reálná síť",
)

const METRIC_LABELS = Dict(
    :clusters => "počet clusterů",
    :pad => "průměrná párová vzdálenost",
    :polarization => "polarizace",
)

const MODEL_COLORS = Dict(
    "deffuant" => :steelblue,
    "hk" => :darkorange,
    "deffuant_hrdcr" => :seagreen,
    "deffuant_anchored" => :darkorange,
    "deffuant_partisan" => :purple,
)

const MODEL_LINESTYLES = Dict(
    "deffuant" => :solid,
    "hk" => :dash,
    "deffuant_hrdcr" => :dot,
    "deffuant_anchored" => :dashdot,
    "deffuant_partisan" => :dashdotdot,
)

const MODEL_MARKERS = Dict(
    "deffuant" => :circle,
    "hk" => :diamond,
    "deffuant_hrdcr" => :utriangle,
    "deffuant_anchored" => :square,
    "deffuant_partisan" => :star5,
)

const SERIES_LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]
const SERIES_MARKERS = [:circle, :diamond, :utriangle, :square, :star5]

function _metric_mean_col(metric::Symbol)
    metric == :clusters && return :clusters_mean
    metric == :pad && return :pad_mean
    metric == :polarization && return :polarization_mean
    throw(ArgumentError("Unknown metric: $metric"))
end

function _model_label(model::AbstractString)
    return get(MODEL_LABELS, String(model), String(model))
end

function _network_label(network::AbstractString)
    return get(NETWORK_LABELS, String(network), String(network))
end


function _network_locative(network::AbstractString)
    network in ("csv", "real") && return "reálné síti"
    network == "er" && return "síti ER"
    network == "ba" && return "síti BA"
    return "síti $(uppercase(String(network)))"
end

function _metric_label(metric::Symbol)
    return get(METRIC_LABELS, metric, String(metric))
end

function _normalize_network_names!(df::DataFrame)
    :network in propertynames(df) || return df
    df.network = [String(n) == "csv" ? "real" : String(n) for n in df.network]
    return df
end

function _read_senate_reference(path::AbstractString)
    if isfile(path)
        df = CSV.read(path, DataFrame)
        senate = df[df.chamber .== "Senate", :]
        if !isempty(senate)
            last = maximum(senate.congress)
            ref = senate[senate.congress .== last, :]
            if nrow(ref) == 1
                return (
                    congress=Int(ref.congress[1]),
                    n=Int(ref.n[1]),
                    clusters=Float64(ref.clusters[1]),
                    pad=Float64(ref.pad[1]),
                    polarization=Float64(ref.polarization[1]),
                )
            end
        end
    end

    # Fallback z poslední ruční analýzy Senate 119.
    return (congress=119, n=103, clusters=10.0, pad=0.2719649723967256, polarization=0.4131908000162851)
end

function _read_senate_series(path::AbstractString)
    if isfile(path)
        df = CSV.read(path, DataFrame)
        senate = df[df.chamber .== "Senate", :]
        if !isempty(senate)
            sort!(senate, :congress)
            first_c = minimum(senate.congress)
            last_c = maximum(senate.congress)
            span = max(last_c - first_c, 1)
            senate.relative_t = (Float64.(senate.congress) .- first_c) ./ span
            return senate
        end
    end
    return DataFrame()
end

function _add_relative_time!(df::DataFrame)
    isempty(df) && return df
    sort!(df, :congress)
    first_c = minimum(df.congress)
    last_c = maximum(df.congress)
    span = max(last_c - first_c, 1)
    df.relative_t = (Float64.(df.congress) .- first_c) ./ span
    return df
end

function _read_chamber_series(path::AbstractString)
    if !isfile(path)
        return Dict{String,DataFrame}()
    end

    df = CSV.read(path, DataFrame)
    :chamber in propertynames(df) || return Dict{String,DataFrame}()

    series = Dict{String,DataFrame}()
    for chamber in sort(unique(String.(df.chamber)))
        sub = df[df.chamber .== chamber, :]
        isempty(sub) && continue
        series[chamber] = _add_relative_time!(copy(sub))
    end
    return series
end

function _slug(value::AbstractString)
    return lowercase(replace(String(value), " " => "_", "/" => "_"))
end

function _final_rows(dyn::DataFrame)
    final_t = maximum(dyn.t)
    return dyn[dyn.t .== final_t, :]
end

function _final_by_epsilon(dyn::DataFrame)
    final = _final_rows(dyn)
    return combine(
        groupby(final, [:network, :model, :epsilon]),
        :clusters_mean => mean => :clusters,
        :pad_mean => mean => :pad,
        :polarization_mean => mean => :polarization,
    )
end

function _final_by_epsilon_mu(dyn::DataFrame, network::String, model::String, metric::Symbol)
    final = _final_rows(dyn)
    sub = final[(final.network .== network) .& (final.model .== model), :]
    isempty(sub) && return Float64[], Float64[], Matrix{Float64}(undef, 0, 0)

    col = _metric_mean_col(metric)
    eps_values = sort(unique(Float64.(sub.epsilon)))
    mu_values = sort(unique(Float64.(sub.mu)))
    z = fill(NaN, length(eps_values), length(mu_values))

    for r in eachrow(sub)
        i = findfirst(==(Float64(r.epsilon)), eps_values)
        j = findfirst(==(Float64(r.mu)), mu_values)
        if i !== nothing && j !== nothing
            z[i, j] = Float64(r[col])
        end
    end

    return eps_values, mu_values, z
end

function _safe_mkpath(path::AbstractString)
    mkpath(dirname(path))
end

function _heatmap_with_senate_reference(dyn::DataFrame, ref, network::String, model::String, metric::Symbol, out_path::AbstractString)
    eps_values, mu_values, z = _final_by_epsilon_mu(dyn, network, model, metric)
    isempty(eps_values) && return nothing

    ref_value = Float64(getproperty(ref, metric))
    metric_label = _metric_label(metric)
    _safe_mkpath(out_path)

    p = Plots.heatmap(
        mu_values,
        eps_values,
        z;
        xlabel="μ",
        ylabel="ε",
        title="$(_model_label(model)) na $(_network_locative(network)): $(metric_label)",
        colorbar_title="",
        color=:viridis,
        size=(900, 650),
        right_margin=12Plots.mm,
    )
    Plots.annotate!(
        p,
        minimum(mu_values),
        maximum(eps_values),
        Plots.text("Senát $(ref.congress): $(round(ref_value, digits=3))", 9, :left, :white),
    )

    if minimum(skipmissing(vec(z))) <= ref_value <= maximum(skipmissing(vec(z)))
        Plots.contour!(
            p,
            mu_values,
            eps_values,
            z;
            levels=[ref_value],
            color=:white,
            linewidth=3,
            label="Senát $(ref.congress)",
        )
    end

    Plots.savefig(p, out_path)
    return nothing
end

function _hk_line_with_reference(final_eps::DataFrame, ref, network::String, metric::Symbol, out_path::AbstractString)
    sub = final_eps[(final_eps.network .== network) .& (final_eps.model .== "hk"), :]
    isempty(sub) && return nothing
    sort!(sub, :epsilon)

    ref_value = Float64(getproperty(ref, metric))
    metric_label = _metric_label(metric)
    _safe_mkpath(out_path)

    p = Plots.plot(
        sub.epsilon,
        sub[!, metric];
        marker=:circle,
        linewidth=3,
        color=MODEL_COLORS["hk"],
        label="HK",
        xlabel="ε",
        ylabel=metric_label,
        title="HK na $(_network_locative(network)): finální $(metric_label)",
        size=(900, 550),
        legend=:topright,
    )
    Plots.hline!(p, [ref_value]; color=:black, linestyle=:dash, linewidth=2.5, label="Senát $(ref.congress)")
    Plots.savefig(p, out_path)
    return nothing
end

function _model_metric_space(final_eps::DataFrame, ref, network::String, out_path::AbstractString)
    sub = final_eps[final_eps.network .== network, :]
    isempty(sub) && return nothing

    _safe_mkpath(out_path)
    p = Plots.plot(
        xlabel="průměrná párová vzdálenost",
        ylabel="polarizace",
        title="Metrický prostor na $(_network_locative(network)) (detail)",
        legend=:outertopright,
        size=(950, 650),
        ylims=(0, 1.2),
    )

    for model in MODEL_ORDER
        ms = sub[sub.model .== model, :]
        isempty(ms) && continue
        Plots.scatter!(
            p,
            ms.pad,
            ms.polarization;
            color=get(MODEL_COLORS, model, :steelblue),
            markerstrokewidth=0,
            markersize=5 .+ 0.7 .* ms.clusters,
            label=_model_label(model),
            alpha=0.85,
        )
    end

    Plots.scatter!(
        p,
        [ref.pad],
        [ref.polarization];
        color=:red,
        marker=:star5,
        markersize=13,
        markerstrokecolor=:black,
        label="Senát $(ref.congress)",
    )

    Plots.savefig(p, out_path)
    return nothing
end

function _comparison_distance_table(final_eps::DataFrame, ref)
    rows = DataFrame(
        network=String[],
        model=String[],
        epsilon=Float64[],
        clusters=Float64[],
        pad=Float64[],
        polarization=Float64[],
        relative_distance=Float64[],
    )

    for r in eachrow(final_eps)
        dc = abs(Float64(r.clusters) - ref.clusters) / max(ref.clusters, eps())
        dp = abs(Float64(r.pad) - ref.pad) / max(ref.pad, eps())
        dpol = abs(Float64(r.polarization) - ref.polarization) / max(ref.polarization, eps())
        push!(rows, (
            String(r.network),
            String(r.model),
            Float64(r.epsilon),
            Float64(r.clusters),
            Float64(r.pad),
            Float64(r.polarization),
            mean([dc, dp, dpol]),
        ))
    end

    sort!(rows, :relative_distance)
    return rows
end

function _interp_linear(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, x::Real)
    n = length(xs)
    n == 0 && return NaN
    n == 1 && return Float64(ys[1])
    x <= xs[1] && return Float64(ys[1])
    x >= xs[end] && return Float64(ys[end])

    i = searchsortedlast(xs, x)
    i = clamp(i, 1, n - 1)
    x0 = Float64(xs[i])
    x1 = Float64(xs[i + 1])
    y0 = Float64(ys[i])
    y1 = Float64(ys[i + 1])
    x1 == x0 && return y0
    w = (Float64(x) - x0) / (x1 - x0)
    return y0 + w * (y1 - y0)
end

function _senate_time_distance_table(dyn::DataFrame, senate_series::DataFrame)
    isempty(senate_series) && return DataFrame()
    rows = DataFrame(
        network=String[],
        model=String[],
        epsilon=Float64[],
        clusters_rmse=Float64[],
        pad_rmse=Float64[],
        polarization_rmse=Float64[],
        combined_relative_rmse=Float64[],
    )

    senate_x = Float64.(senate_series.relative_t)
    senate_clusters = Float64.(senate_series.clusters)
    senate_pad = Float64.(senate_series.pad)
    senate_pol = Float64.(senate_series.polarization)
    scale_clusters = max(mean(abs.(senate_clusters)), eps())
    scale_pad = max(mean(abs.(senate_pad)), eps())
    scale_pol = max(mean(abs.(senate_pol)), eps())

    for sdf in groupby(dyn, [:network, :model, :epsilon])
        model_name = String(first(sdf.model))
        net_name = String(first(sdf.network))
        eps_value = Float64(first(sdf.epsilon))

        avg = combine(
            groupby(sdf, :t),
            :clusters_mean => mean => :clusters,
            :pad_mean => mean => :pad,
            :polarization_mean => mean => :polarization,
        )
        sort!(avg, :t)
        tmin = minimum(avg.t)
        tmax = maximum(avg.t)
        span = max(tmax - tmin, 1)
        avg.relative_t = (Float64.(avg.t) .- tmin) ./ span

        model_x = Float64.(avg.relative_t)
        model_clusters = Float64.(avg.clusters)
        model_pad = Float64.(avg.pad)
        model_pol = Float64.(avg.polarization)

        dc = Float64[]
        dp = Float64[]
        dpol = Float64[]
        for i in eachindex(senate_x)
            x = senate_x[i]
            push!(dc, _interp_linear(model_x, model_clusters, x) - senate_clusters[i])
            push!(dp, _interp_linear(model_x, model_pad, x) - senate_pad[i])
            push!(dpol, _interp_linear(model_x, model_pol, x) - senate_pol[i])
        end

        clusters_rmse = sqrt(mean(dc .^ 2))
        pad_rmse = sqrt(mean(dp .^ 2))
        pol_rmse = sqrt(mean(dpol .^ 2))
        combined = mean([clusters_rmse / scale_clusters, pad_rmse / scale_pad, pol_rmse / scale_pol])

        push!(rows, (net_name, model_name, eps_value, clusters_rmse, pad_rmse, pol_rmse, combined))
    end

    sort!(rows, :combined_relative_rmse)
    return rows
end

function _distance_lollipop(distance::DataFrame, out_path::AbstractString; top_n::Int=12)
    top = first(distance, min(top_n, nrow(distance)))
    y = reverse(1:nrow(top))
    labels = ["$(_network_label(r.network)) | $(_model_label(r.model)) | ε=$(r.epsilon)" for r in eachrow(top)]

    _safe_mkpath(out_path)
    p = Plots.plot(
        xlabel="průměrná relativní odchylka od Senátu",
        ylabel="",
        title="Parametry nejbližší reálnému Senátu",
        yticks=(y, labels),
        legend=false,
        size=(1300, 760),
        left_margin=22Plots.mm,
        right_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )

    for (xi, yi, model) in zip(top.relative_distance, y, top.model)
        Plots.plot!(p, [0.0, xi], [yi, yi]; color=:gray70, linewidth=2, label=false)
        Plots.scatter!(
            p,
            [xi],
            [yi];
            color=get(MODEL_COLORS, String(model), :steelblue),
            markersize=8,
            markerstrokewidth=0,
            label=false,
        )
        Plots.annotate!(p, xi, yi, Plots.text("  $(round(xi, digits=3))", 9, :left))
    end

    Plots.xlims!(p, 0, maximum(top.relative_distance) * 1.25)
    Plots.savefig(p, out_path)
    return nothing
end

function _dynamic_selected(dyn::DataFrame, distance::DataFrame, ref, metric::Symbol, out_path::AbstractString; top_n::Int=5)
    selected = first(distance, min(top_n, nrow(distance)))
    col = _metric_mean_col(metric)
    metric_label = _metric_label(metric)
    ref_value = Float64(getproperty(ref, metric))

    _safe_mkpath(out_path)
    p = Plots.plot(
        xlabel="časový krok uloženého snapshotu",
        ylabel=metric_label,
        title="Dynamika nejbližších nastavení: $(metric_label)",
        legend=:outertopright,
        size=(1100, 650),
    )

    for (i, r) in enumerate(eachrow(selected))
        sub = dyn[
            (dyn.network .== r.network) .&
            (dyn.model .== r.model) .&
            (dyn.epsilon .== r.epsilon),
            :,
        ]
        isempty(sub) && continue
        avg = combine(groupby(sub, :t), col => mean => :value)
        sort!(avg, :t)
        Plots.plot!(
            p,
            avg.t,
            avg.value;
            linewidth=2.2,
            linestyle=SERIES_LINESTYLES[mod1(i, length(SERIES_LINESTYLES))],
            label="$(_network_label(r.network)) | $(_model_label(r.model)) | ε=$(r.epsilon)",
        )
    end

    Plots.hline!(p, [ref_value]; color=:black, linestyle=:dash, linewidth=2.5, label="Senát $(ref.congress)")
    Plots.savefig(p, out_path)
    return nothing
end

function _plot_senate_history(senate_series::DataFrame, out_path::AbstractString)
    isempty(senate_series) && return nothing
    _safe_mkpath(out_path)

    p1 = Plots.plot(
        senate_series.congress,
        senate_series.clusters;
        marker=:circle,
        linewidth=2.5,
        color=:firebrick,
        label=false,
        ylabel="počet clusterů",
        title="Historický vývoj Senátu podle NOMINATE",
    )
    p2 = Plots.plot(
        senate_series.congress,
        senate_series.pad;
        marker=:circle,
        linewidth=2.5,
        color=:steelblue,
        label=false,
        ylabel="PAD",
    )
    p3 = Plots.plot(
        senate_series.congress,
        senate_series.polarization;
        marker=:circle,
        linewidth=2.5,
        color=:seagreen,
        label=false,
        xlabel="kongres",
        ylabel="polarizace",
    )

    p = Plots.plot(p1, p2, p3; layout=(3, 1), size=(1100, 900))
    Plots.savefig(p, out_path)
    return nothing
end

function _plot_all_chambers_history(chamber_series::Dict{String,DataFrame}, metric::Symbol, out_path::AbstractString)
    isempty(chamber_series) && return nothing
    metric_label = _metric_label(metric)
    _safe_mkpath(out_path)

    colors = Dict("All" => :black, "House" => :steelblue, "Senate" => :firebrick)
    styles = Dict("All" => :solid, "House" => :dash, "Senate" => :dot)
    markers = Dict("All" => :circle, "House" => :diamond, "Senate" => :utriangle)
    p = Plots.plot(
        xlabel="kongres",
        ylabel=metric_label,
        title="Historický vývoj podle komory: $(metric_label)",
        legend=:outertopright,
        size=(1100, 650),
    )

    for chamber in ["All", "House", "Senate"]
        haskey(chamber_series, chamber) || continue
        sub = chamber_series[chamber]
        Plots.plot!(
            p,
            sub.congress,
            sub[!, metric];
            marker=get(markers, chamber, :circle),
            linestyle=get(styles, chamber, :solid),
            linewidth=2.3,
            color=get(colors, chamber, :gray40),
            label=chamber,
        )
    end

    Plots.savefig(p, out_path)
    return nothing
end

function _plot_chamber_history(series::DataFrame, chamber::String, out_path::AbstractString)
    isempty(series) && return nothing
    _safe_mkpath(out_path)

    p1 = Plots.plot(
        series.congress,
        series.clusters;
        marker=:circle,
        linewidth=2.5,
        color=:firebrick,
        label=false,
        ylabel="počet clusterů",
        title="Historický vývoj: $(chamber)",
    )
    p2 = Plots.plot(
        series.congress,
        series.pad;
        marker=:circle,
        linewidth=2.5,
        color=:steelblue,
        label=false,
        ylabel="PAD",
    )
    p3 = Plots.plot(
        series.congress,
        series.polarization;
        marker=:circle,
        linewidth=2.5,
        color=:seagreen,
        label=false,
        xlabel="kongres",
        ylabel="polarizace",
    )

    p = Plots.plot(p1, p2, p3; layout=(3, 1), size=(1100, 900))
    Plots.savefig(p, out_path)
    return nothing
end

function _senate_time_lollipop(distance::DataFrame, out_path::AbstractString; top_n::Int=12)
    isempty(distance) && return nothing
    top = first(distance, min(top_n, nrow(distance)))
    y = reverse(1:nrow(top))
    labels = ["$(_network_label(r.network)) | $(_model_label(r.model)) | ε=$(r.epsilon)" for r in eachrow(top)]

    _safe_mkpath(out_path)
    p = Plots.plot(
        xlabel="průměrná relativní RMSE vůči průběhu Senátu",
        ylabel="",
        title="Modelové průběhy nejbližší historickému vývoji Senátu",
        yticks=(y, labels),
        legend=false,
        size=(1350, 760),
        left_margin=24Plots.mm,
        right_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )

    for (xi, yi, model) in zip(top.combined_relative_rmse, y, top.model)
        Plots.plot!(p, [0.0, xi], [yi, yi]; color=:gray70, linewidth=2, label=false)
        Plots.scatter!(
            p,
            [xi],
            [yi];
            color=get(MODEL_COLORS, String(model), :steelblue),
            markersize=8,
            markerstrokewidth=0,
            label=false,
        )
        Plots.annotate!(p, xi, yi, Plots.text("  $(round(xi, digits=3))", 9, :left))
    end

    Plots.xlims!(p, 0, maximum(top.combined_relative_rmse) * 1.25)
    Plots.savefig(p, out_path)
    return nothing
end

function _plot_senate_vs_model_time(dyn::DataFrame, senate_series::DataFrame, distance::DataFrame, metric::Symbol, out_path::AbstractString; top_n::Int=5)
    isempty(senate_series) && return nothing
    isempty(distance) && return nothing

    selected = first(distance, min(top_n, nrow(distance)))
    col = _metric_mean_col(metric)
    metric_label = _metric_label(metric)

    _safe_mkpath(out_path)
    p = Plots.plot(
        xlabel="relativní čas",
        ylabel=metric_label,
        title="Historický vývoj Senátu vs modely: $(metric_label)",
        legend=:outertopright,
        size=(1150, 680),
    )

    Plots.plot!(
        p,
        senate_series.relative_t,
        senate_series[!, metric];
        color=:black,
        linewidth=3.2,
        marker=:circle,
        label="Senát (všechny kongresy)",
    )

    for (i, r) in enumerate(eachrow(selected))
        sub = dyn[
            (dyn.network .== r.network) .&
            (dyn.model .== r.model) .&
            (dyn.epsilon .== r.epsilon),
            :,
        ]
        isempty(sub) && continue
        avg = combine(groupby(sub, :t), col => mean => :value)
        sort!(avg, :t)
        tmin = minimum(avg.t)
        tmax = maximum(avg.t)
        span = max(tmax - tmin, 1)
        xs = (Float64.(avg.t) .- tmin) ./ span
        Plots.plot!(
            p,
            xs,
            avg.value;
            linewidth=2.2,
            color=get(MODEL_COLORS, String(r.model), :steelblue),
            linestyle=SERIES_LINESTYLES[mod1(i, length(SERIES_LINESTYLES))],
            alpha=0.75,
            label="$(_network_label(r.network)) | $(_model_label(r.model)) | ε=$(r.epsilon)",
        )
    end

    Plots.savefig(p, out_path)
    return nothing
end

function _plot_chamber_vs_model_time(dyn::DataFrame, series::DataFrame, chamber::String, distance::DataFrame, metric::Symbol, out_path::AbstractString; top_n::Int=5)
    isempty(series) && return nothing
    isempty(distance) && return nothing

    selected = first(distance, min(top_n, nrow(distance)))
    col = _metric_mean_col(metric)
    metric_label = _metric_label(metric)

    _safe_mkpath(out_path)
    p = Plots.plot(
        xlabel="relativní čas",
        ylabel=metric_label,
        title="$(chamber) vs modely: $(metric_label)",
        legend=:outertopright,
        size=(1150, 680),
    )

    Plots.plot!(
        p,
        series.relative_t,
        series[!, metric];
        color=:black,
        linewidth=3.2,
        marker=:circle,
        label="$(chamber) (kongresy)",
    )

    for (i, r) in enumerate(eachrow(selected))
        sub = dyn[
            (dyn.network .== r.network) .&
            (dyn.model .== r.model) .&
            (dyn.epsilon .== r.epsilon),
            :,
        ]
        isempty(sub) && continue
        avg = combine(groupby(sub, :t), col => mean => :value)
        sort!(avg, :t)
        tmin = minimum(avg.t)
        tmax = maximum(avg.t)
        span = max(tmax - tmin, 1)
        xs = (Float64.(avg.t) .- tmin) ./ span
        Plots.plot!(
            p,
            xs,
            avg.value;
            linewidth=2.2,
            color=get(MODEL_COLORS, String(r.model), :steelblue),
            linestyle=SERIES_LINESTYLES[mod1(i, length(SERIES_LINESTYLES))],
            alpha=0.75,
            label="$(_network_label(r.network)) | $(_model_label(r.model)) | ε=$(r.epsilon)",
        )
    end

    Plots.savefig(p, out_path)
    return nothing
end

function _small_multiple_final_metric(final_eps::DataFrame, ref, metric::Symbol, out_path::AbstractString; include_hrdcr::Bool=true)
    models = include_hrdcr ? MODEL_ORDER : ["deffuant", "hk"]
    metric_label = _metric_label(metric)
    ref_value = Float64(getproperty(ref, metric))
    panels = Plots.Plot[]

    for network in NETWORK_ORDER
        sub = final_eps[final_eps.network .== network, :]
        isempty(sub) && continue
        p = Plots.plot(
            xlabel="ε",
            ylabel=metric_label,
            title=_network_label(network),
            legend=:topright,
            size=(850, 420),
        )
        for model in models
            ms = sub[sub.model .== model, :]
            isempty(ms) && continue
            sort!(ms, :epsilon)
            Plots.plot!(
                p,
                ms.epsilon,
                ms[!, metric];
                marker=get(MODEL_MARKERS, model, :circle),
                linestyle=get(MODEL_LINESTYLES, model, :solid),
                linewidth=2.4,
                color=get(MODEL_COLORS, model, :steelblue),
                label=_model_label(model),
            )
        end
        Plots.hline!(p, [ref_value]; color=:black, linestyle=:dash, linewidth=2.2, label="Senát $(ref.congress)")
        push!(panels, p)
    end

    isempty(panels) && return nothing
    _safe_mkpath(out_path)
    plt = Plots.plot(
        panels...;
        layout=(length(panels), 1),
        size=(1000, 360 * length(panels)),
        plot_title="Finální $(metric_label): porovnání modelů se Senátem",
    )
    Plots.savefig(plt, out_path)
    return nothing
end

function _write_markdown(out_dir::AbstractString, distance::DataFrame, ref, time_distance::DataFrame)
    mkpath(out_dir)
    best = first(distance, min(5, nrow(distance)))
    best_time = isempty(time_distance) ? DataFrame() : first(time_distance, min(5, nrow(time_distance)))

    chapter7 = """
# Poznámky ke kapitole 7: Experimentální analýza

## 7.1 Deffuant

Do části o Deffuantově modelu je vhodné vložit grafy `figures/01_deffuant/deffuant_real_clusters_heatmap_senate.png`, `figures/01_deffuant/deffuant_real_pad_heatmap_senate.png` a `figures/01_deffuant/deffuant_real_polarization_heatmap_senate.png`. Grafy ukazují finální hodnoty sledovaných charakteristik na reálné síti v závislosti na parametrech ε a μ. Přerušovaná/kontrastní referenční linie odpovídá hodnotě vypočtené z rozložení `nominate_dim1` pro 119. kongres Senátu.

Z těchto grafů lze popsat, že Deffuantův model vykazuje plynulou parametrickou závislost. Nejvýraznější fragmentace vzniká při nízkém prahu důvěry, zatímco vyšší hodnota ε podporuje sbližování názorů. Parametr μ mění rychlost a intenzitu párového přiblížení, ale samotný počet výsledných clusterů je v datech nejcitlivější především na ε.

## 7.2 Hegselmann-Krause

Do části o HK modelu vlož grafy `figures/02_hk/hk_real_clusters_line_senate.png`, `figures/02_hk/hk_real_pad_line_senate.png` a `figures/02_hk/hk_real_polarization_line_senate.png`. Tyto grafy porovnávají finální hodnoty HK modelu s referenční hodnotou Senátu.

HK model má pouze parametr ε, proto jsou výsledky zobrazeny jako jednorozměrné křivky. Ve srovnání s Deffuantovým modelem je přechod k nízkému počtu clusterů ostřejší. To odpovídá synchronní aktualizaci názorů, při které se lokální konsenzus šíří rychleji než v párovém Deffuantově mechanismu.

## 7.3 Porovnání Deffuantova modelu a HK

Pro porovnání základních modelů použij `figures/03_comparison/base_models_clusters_senate.png`, `figures/03_comparison/base_models_pad_senate.png`, `figures/03_comparison/base_models_polarization_senate.png` a `figures/03_comparison/real_metric_space_base_models_senate.png`.

Grafy ukazují, že oba modely reprodukují obecný vztah mezi prahem důvěry a mírou fragmentace, ale liší se rychlostí přechodu mezi fragmentovaným stavem a konsenzem. Deffuantův model dává jemnější parametrickou škálu výsledků, zatímco HK model se u vyšších hodnot ε rychle přibližuje konsenzu.

Jako doplňkové srovnání lze použít také historický průběh Senátu přes všechny kongresy. K tomu slouží grafy `figures/06_senate_time_comparison/senate_historical_dynamic_characteristics.png`, `figures/06_senate_time_comparison/senate_vs_models_clusters.png`, `figures/06_senate_time_comparison/senate_vs_models_pad.png` a `figures/06_senate_time_comparison/senate_vs_models_polarization.png`. Časová osa je v tomto případě normalizovaná na interval 0 až 1, protože modelová iterace a historická posloupnost kongresů nejsou stejný typ času.

Pro širší pohled lze doplnit také grafy ze složky `figures/07_chambers_time_comparison`. Ty porovnávají modelové průběhy nejen se Senátem, ale také se Sněmovnou reprezentantů (`House`) a s celkovým rozložením všech členů Kongresu (`All`).
"""

    chapter8 = """
# Poznámky ke kapitole 8: Rozšířený model HRDCR

## Experimentální ověření HRDCR

Do části experimentálního ověření HRDCR vlož grafy `figures/04_hrdcr/hrdcr_real_clusters_heatmap_senate.png`, `figures/04_hrdcr/hrdcr_real_pad_heatmap_senate.png` a `figures/04_hrdcr/hrdcr_real_polarization_heatmap_senate.png`. Tyto grafy zachycují, jak se rozšířený model chová při změně ε a μ na reálné síti.

HRDCR je vhodné interpretovat jako strukturálně citlivou variantu Deffuantova modelu. Oproti původnímu modelu zohledňuje lokální překryv sousedství a rozdílnou setrvačnost uzlů. Výsledky proto nejsou dány pouze vzdáleností názorů, ale také tím, v jaké části sítě interakce probíhá.

## Porovnání s Deffuantovým modelem a HK

Pro společné porovnání všech modelů použij grafy `figures/05_all_models/all_models_clusters_senate.png`, `figures/05_all_models/all_models_pad_senate.png`, `figures/05_all_models/all_models_polarization_senate.png`, `figures/05_all_models/real_metric_space_all_models_senate.png` a `figures/05_all_models/closest_to_senate_ranked.png`.

Referenční hodnoty posledního Senátu jsou: počet clusterů = $(round(ref.clusters, digits=3)), průměrná párová vzdálenost = $(round(ref.pad, digits=3)) a polarizace = $(round(ref.polarization, digits=3)). Podle průměrné relativní odchylky se k této referenci nejvíce přiblížila tato nastavení:

$(join(["- $(_network_label(r.network)), $(_model_label(r.model)), ε = $(r.epsilon), relativní odchylka = $(round(r.relative_distance, digits=3))" for r in eachrow(best)], "\n"))

Tyto výsledky je vhodné chápat opatrně. Referenční Senát není síťová simulace, ale empirické rozložení ideologických skóre. Srovnání proto neslouží jako přímá validace modelů, ale jako orientační měřítko toho, zda modely dokážou vytvořit podobnou míru fragmentace, rozptylu a polarizace.

Pro porovnání celého průběhu metrik v čase lze vložit také `figures/06_senate_time_comparison/senate_time_distance_ranked.png`. Tento graf neporovnává pouze finální stav, ale podobnost celé časové řady vůči historickému vývoji Senátu. Nejbližší průběhy podle průměrné relativní RMSE jsou:

$(isempty(best_time) ? "- časové porovnání nebylo vytvořeno" : join(["- $(_network_label(r.network)), $(_model_label(r.model)), ε = $(r.epsilon), relativní RMSE = $(round(r.combined_relative_rmse, digits=3))" for r in eachrow(best_time)], "\n"))

Pokud je cílem vyhnout se příliš úzké interpretaci pouze přes Senát, je možné použít i složku `figures/07_chambers_time_comparison`. Obsahuje samostatná porovnání pro `All`, `House` a `Senate`, takže lze posoudit, zda je podobnost modelů stabilní napříč různými částmi kongresových dat.
"""

    readme = """
# Thesis Figures

Tato složka obsahuje prezentační grafy vytvořené z `prod_real_run_1` a referenčních hodnot 119. kongresu Senátu.

Struktura:

- `figures/01_deffuant`: grafy pro Deffuantův model
- `figures/02_hk`: grafy pro HK model
- `figures/03_comparison`: porovnání Deffuantova modelu a HK
- `figures/04_hrdcr`: experimentální analýza HRDCR
- `figures/05_all_models`: společné porovnání všech modelů
- `figures/06_senate_time_comparison`: porovnání modelové dynamiky s historickým vývojem Senátu přes všechny kongresy
- `figures/07_chambers_time_comparison`: porovnání modelové dynamiky s historickým vývojem komor All, House a Senate
- `tables`: CSV tabulky s finálními hodnotami a vzdáleností od Senátu
- `chapter_7_text.md`: návrh textu pro kapitolu 7
- `chapter_8_text.md`: návrh textu pro kapitolu 8
"""

    write(joinpath(out_dir, "chapter_7_text.md"), chapter7)
    write(joinpath(out_dir, "chapter_8_text.md"), chapter8)
    write(joinpath(out_dir, "README.md"), readme)
end

function generate_thesis_figures(;
    prod_summary_path::AbstractString,
    senate_summary_path::AbstractString,
    out_dir::AbstractString,
)
    dyn = CSV.read(prod_summary_path, DataFrame)
    _normalize_network_names!(dyn)
    ref = _read_senate_reference(senate_summary_path)
    senate_series = _read_senate_series(senate_summary_path)
    chamber_series = _read_chamber_series(senate_summary_path)
    final_eps = _final_by_epsilon(dyn)
    distance = _comparison_distance_table(final_eps, ref)
    time_distance = _senate_time_distance_table(dyn, senate_series)

    tables_dir = joinpath(out_dir, "tables")
    mkpath(tables_dir)
    CSV.write(joinpath(tables_dir, "final_by_epsilon.csv"), final_eps)
    CSV.write(joinpath(tables_dir, "distance_to_senate_ranked.csv"), distance)
    CSV.write(joinpath(tables_dir, "senate_reference.csv"), DataFrame([ref]))
    if !isempty(senate_series)
        CSV.write(joinpath(tables_dir, "senate_historical_series.csv"), senate_series)
    end
    if !isempty(time_distance)
        CSV.write(joinpath(tables_dir, "time_distance_to_senate_ranked.csv"), time_distance)
    end
    for chamber in sort(collect(keys(chamber_series)))
        chamber_distance = _senate_time_distance_table(dyn, chamber_series[chamber])
        isempty(chamber_distance) && continue
        CSV.write(joinpath(tables_dir, "time_distance_to_$(lowercase(chamber))_ranked.csv"), chamber_distance)
    end

    for metric in METRICS
        _heatmap_with_senate_reference(
            dyn,
            ref,
            "real",
            "deffuant",
            metric,
            joinpath(out_dir, "figures", "01_deffuant", "deffuant_real_$(metric)_heatmap_senate.png"),
        )
        _hk_line_with_reference(
            final_eps,
            ref,
            "real",
            metric,
            joinpath(out_dir, "figures", "02_hk", "hk_real_$(metric)_line_senate.png"),
        )
        _heatmap_with_senate_reference(
            dyn,
            ref,
            "real",
            "deffuant_hrdcr",
            metric,
            joinpath(out_dir, "figures", "04_hrdcr", "hrdcr_real_$(metric)_heatmap_senate.png"),
        )
        _small_multiple_final_metric(
            final_eps,
            ref,
            metric,
            joinpath(out_dir, "figures", "03_comparison", "base_models_$(metric)_senate.png");
            include_hrdcr=false,
        )
        _small_multiple_final_metric(
            final_eps,
            ref,
            metric,
            joinpath(out_dir, "figures", "05_all_models", "all_models_$(metric)_senate.png");
            include_hrdcr=true,
        )
        _dynamic_selected(
            dyn,
            distance,
            ref,
            metric,
            joinpath(out_dir, "figures", "05_all_models", "closest_settings_dynamic_$(metric).png"),
        )
        _plot_senate_vs_model_time(
            dyn,
            senate_series,
            time_distance,
            metric,
            joinpath(out_dir, "figures", "06_senate_time_comparison", "senate_vs_models_$(metric).png"),
        )
        _plot_all_chambers_history(
            chamber_series,
            metric,
            joinpath(out_dir, "figures", "07_chambers_time_comparison", "all_chambers_history_$(metric).png"),
        )
    end

    _model_metric_space(
        final_eps[final_eps.model .!= "deffuant_hrdcr", :],
        ref,
        "real",
        joinpath(out_dir, "figures", "03_comparison", "real_metric_space_base_models_senate.png"),
    )
    _model_metric_space(
        final_eps,
        ref,
        "real",
        joinpath(out_dir, "figures", "05_all_models", "real_metric_space_all_models_senate.png"),
    )
    _distance_lollipop(
        distance,
        joinpath(out_dir, "figures", "05_all_models", "closest_to_senate_ranked.png"),
    )
    _plot_senate_history(
        senate_series,
        joinpath(out_dir, "figures", "06_senate_time_comparison", "senate_historical_dynamic_characteristics.png"),
    )
    _senate_time_lollipop(
        time_distance,
        joinpath(out_dir, "figures", "06_senate_time_comparison", "senate_time_distance_ranked.png"),
    )
    for chamber in sort(collect(keys(chamber_series)))
        series = chamber_series[chamber]
        chamber_distance = _senate_time_distance_table(dyn, series)
        isempty(chamber_distance) && continue
        chamber_slug = _slug(chamber)
        _plot_chamber_history(
            series,
            chamber,
            joinpath(out_dir, "figures", "07_chambers_time_comparison", "$(chamber_slug)_historical_dynamic_characteristics.png"),
        )
        _senate_time_lollipop(
            chamber_distance,
            joinpath(out_dir, "figures", "07_chambers_time_comparison", "$(chamber_slug)_time_distance_ranked.png"),
        )
        for metric in METRICS
            _plot_chamber_vs_model_time(
                dyn,
                series,
                chamber,
                chamber_distance,
                metric,
                joinpath(out_dir, "figures", "07_chambers_time_comparison", "$(chamber_slug)_vs_models_$(metric).png"),
            )
        end
    end

    _write_markdown(out_dir, distance, ref, time_distance)
    println("Thesis figures saved to: $out_dir")
    return (final_by_epsilon=final_eps, distance=distance, time_distance=time_distance, senate_reference=ref)
end

end # module
