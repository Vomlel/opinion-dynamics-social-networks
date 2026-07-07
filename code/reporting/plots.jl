module Plotting

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

ENV["GKSwstype"] = "100"

import Plots

using CSV
using DataFrames
using Statistics
using Base.Filesystem: dirname, mkpath
using Plots

const SERIES_LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]
const SERIES_MARKERS = [:circle, :diamond, :utriangle, :square, :star5]

_cz_number(value::Real; digits::Int=3) = replace(string(round(Float64(value), digits=digits)), "." => ",")

function _metric_label(metric::Symbol)
    metric == :clusters && return "počet clusterů"
    metric == :pad && return "průměrná párová vzdálenost"
    metric == :polarization && return "polarizace"
    return String(metric)
end

export plot_dynamic_metric,
       plot_dynamic_metric_params,
       plot_parameter_heatmap,
       plot_parameter_surface,
       plot_parameter_contour,
       plot_parameter_slices,
       plot_parameter_line,
       plot_models_parameter_surfaces,
       plot_deffuant_hk_final_comparison,
       plot_deffuant_hk_dynamic_comparison,
       plot_models_final_comparison,
       plot_models_dynamic_comparison,
       plot_opinion_series,
       plot_nominate_dim1_distribution

"""
    _subset(dyn::DataFrame, network::AbstractString, model::AbstractString) -> DataFrame

Vrátí seřazený výřez agregovaných dat pro zvolenou síť a model.
"""
function _subset(dyn::DataFrame, network::AbstractString, model::AbstractString)
    sub = dyn[(dyn.network .== network) .& (dyn.model .== model), :]
    isempty(sub) && return sub
    sort!(sub, :t)
    return sub
end

function _metric_columns(metric::Symbol)
    if metric == :clusters
        return (:clusters_mean, :clusters_q25, :clusters_q50, :clusters_q75)
    elseif metric == :pad
        return (:pad_mean, :pad_q25, :pad_q50, :pad_q75)
    elseif metric == :polarization
        return (:polarization_mean, :polarization_q25, :polarization_q50, :polarization_q75)
    else
        throw(ArgumentError("Unknown metric: $metric"))
    end
end

"""
    plot_dynamic_metric(dyn_summary::DataFrame; ...) -> Nothing

Vykreslí jednu dynamickou metriku pro zvolenou dvojici `(network, model)`.

Do grafu kreslí průměr a volitelně i kvantily 25/50/75.
"""
function plot_dynamic_metric(
    dyn_summary::DataFrame;
    network::String,
    metric::Symbol,
    model::String,
    out_path::String,
    title::String = "",
    xlabel::String = "t",
    ylabel::String = "",
    quantiles::Bool = true,
)

    mean_col, q25_col, q50_col, q75_col = _metric_columns(metric)

    sub = _subset(dyn_summary, network, model)
    isempty(sub) && return nothing

    plt = Plots.plot()
    mkpath(dirname(out_path))
    ys = sub[!, mean_col]
    Plots.plot!(plt, sub.t, ys; label="$(model)_avg")

    if quantiles
        s_q25 = sub[!, q25_col]
        Plots.plot!(plt, sub.t, s_q25; label="$(model)_q25", linestyle=:dash)
        s_q50 = sub[!, q50_col]
        Plots.plot!(plt, sub.t, s_q50; label="$(model)_q50", linestyle=:dashdot)
        s_q75 = sub[!, q75_col]
        Plots.plot!(plt, sub.t, s_q75; label="$(model)_q75", linestyle=:dash)
    end
    Plots.xlabel!(plt, xlabel)
    Plots.ylabel!(plt, ylabel == "" ? String(metric) : ylabel)
    if title != ""
        Plots.title!(plt, title)
    end

    current_ylims = Plots.ylims(plt)
    Plots.ylims!(plt, 0, max(1, current_ylims[2])) # takto abychom mohli nastavit pouze spodni hranici na 0

    Plots.savefig(plt, out_path)

    return nothing
end

function _param_label(model::AbstractString, epsilon::Real, mu::Real)
    if model == "hk" || isnan(Float64(mu))
        return "ε=$(_cz_number(epsilon))"
    end
    return "ε=$(_cz_number(epsilon)), μ=$(_cz_number(mu))"
end

"""
    plot_dynamic_metric_params(dyn_summary::DataFrame; ...) -> Nothing

Vykreslí vývoj metriky pro všechny parametrické kombinace v jednom obrázku.
"""
function plot_dynamic_metric_params(
    dyn_summary::DataFrame;
    network::String,
    metric::Symbol,
    model::String,
    out_path::String,
    title::String = "",
    xlabel::String = "t",
    ylabel::String = "",
)
    mean_col, _, _, _ = _metric_columns(metric)
    sub = _subset(dyn_summary, network, model)
    isempty(sub) && return nothing

    sort!(sub, [:epsilon, :mu, :t])
    mkpath(dirname(out_path))

    plt = Plots.plot(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel == "" ? String(metric) : ylabel,
        legend=:outertopright,
        size=(1100, 650),
        left_margin=16Plots.mm,
        right_margin=8Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )

    for (i, sdf) in enumerate(groupby(sub, [:epsilon, :mu]))
        eps = first(sdf.epsilon)
        mu = first(sdf.mu)
        Plots.plot!(
            plt,
            sdf.t,
            sdf[!, mean_col];
            label=_param_label(model, eps, mu),
            linewidth=1.8,
            linestyle=SERIES_LINESTYLES[mod1(i, length(SERIES_LINESTYLES))],
        )
    end

    Plots.savefig(plt, out_path)
    return nothing
end

function _final_parameter_slice(dyn_summary::DataFrame, network::String, model::String)
    sub = _subset(dyn_summary, network, model)
    isempty(sub) && return sub
    final_t = maximum(sub.t)
    return sub[sub.t .== final_t, :]
end

function _parameter_matrix(dyn_summary::DataFrame, network::String, model::String, metric::Symbol)
    mean_col, _, _, _ = _metric_columns(metric)
    sub = _final_parameter_slice(dyn_summary, network, model)
    isempty(sub) && return Float64[], Float64[], Matrix{Float64}(undef, 0, 0)

    eps_values = sort(unique(Float64.(sub.epsilon)))
    mu_values = sort(unique(Float64.(sub.mu)))
    z = fill(NaN, length(eps_values), length(mu_values))

    for r in eachrow(sub)
        i = findfirst(==(Float64(r.epsilon)), eps_values)
        j = findfirst(==(Float64(r.mu)), mu_values)
        if i !== nothing && j !== nothing
            z[i, j] = Float64(r[mean_col])
        end
    end

    return eps_values, mu_values, z
end

"""
    plot_parameter_heatmap(dyn_summary::DataFrame; ...) -> Nothing

Pro Deffuantův model vykreslí finální hodnotu metriky jako heatmapu
`epsilon × mu`.
"""
function plot_parameter_heatmap(
    dyn_summary::DataFrame;
    network::String,
    model::String,
    metric::Symbol,
    out_path::String,
    title::String = "",
    xlabel::String = "μ",
    ylabel::String = "ε",
)
    eps_values, mu_values, z = _parameter_matrix(dyn_summary, network, model, metric)
    isempty(eps_values) && return nothing

    mkpath(dirname(out_path))
    plt = heatmap(
        mu_values,
        eps_values,
        z;
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        colorbar_title="",
        size=(900, 650),
        left_margin=16Plots.mm,
        right_margin=12Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )
    Plots.savefig(plt, out_path)
    return nothing
end

"""
    plot_parameter_surface(dyn_summary::DataFrame; ...) -> Nothing

Vykreslí finální hodnotu metriky jako 3D povrch nad parametry `mu` a
`epsilon`.
"""
function plot_parameter_surface(
    dyn_summary::DataFrame;
    network::String,
    model::String,
    metric::Symbol,
    out_path::String,
    title::String = "",
    xlabel::String = "μ",
    ylabel::String = "ε",
    zlabel::String = "",
)
    eps_values, mu_values, z = _parameter_matrix(dyn_summary, network, model, metric)
    isempty(eps_values) && return nothing

    mkpath(dirname(out_path))
    plt = surface(
        mu_values,
        eps_values,
        z;
        xlabel=xlabel,
        ylabel=ylabel,
        zlabel=zlabel == "" ? String(metric) : zlabel,
        title=title,
        colorbar_title="",
        size=(900, 650),
        left_margin=16Plots.mm,
        right_margin=12Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
        camera=(45, 30),
    )
    Plots.savefig(plt, out_path)
    return nothing
end

"""
    plot_parameter_contour(dyn_summary::DataFrame; ...) -> Nothing

Vykreslí finální hodnotu metriky jako vrstevnicový graf nad parametry.
"""
function plot_parameter_contour(
    dyn_summary::DataFrame;
    network::String,
    model::String,
    metric::Symbol,
    out_path::String,
    title::String = "",
    xlabel::String = "μ",
    ylabel::String = "ε",
)
    eps_values, mu_values, z = _parameter_matrix(dyn_summary, network, model, metric)
    isempty(eps_values) && return nothing

    mkpath(dirname(out_path))
    plt = contourf(
        mu_values,
        eps_values,
        z;
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        colorbar_title="",
        size=(900, 650),
        left_margin=16Plots.mm,
        right_margin=12Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )
    Plots.savefig(plt, out_path)
    return nothing
end

"""
    plot_parameter_slices(dyn_summary::DataFrame; ...) -> Nothing

Vykreslí finální hodnotu metriky podle `mu`; každá křivka odpovídá jedné
hodnotě `epsilon`.
"""
function plot_parameter_slices(
    dyn_summary::DataFrame;
    network::String,
    model::String,
    metric::Symbol,
    out_path::String,
    title::String = "",
    group_label::String = "ε",
    xlabel::String = "μ",
    ylabel::String = "",
)
    mean_col, _, _, _ = _metric_columns(metric)
    sub = _final_parameter_slice(dyn_summary, network, model)
    isempty(sub) && return nothing

    sort!(sub, [:epsilon, :mu])
    mkpath(dirname(out_path))
    plt = plot(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel == "" ? String(metric) : ylabel,
        legend=:outertopright,
        size=(950, 600),
        left_margin=16Plots.mm,
        right_margin=8Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )

    for (i, sdf) in enumerate(groupby(sub, :epsilon))
        eps = first(sdf.epsilon)
        Plots.plot!(
            plt,
            sdf.mu,
            sdf[!, mean_col];
            marker=SERIES_MARKERS[mod1(i, length(SERIES_MARKERS))],
            linestyle=SERIES_LINESTYLES[mod1(i, length(SERIES_LINESTYLES))],
            linewidth=2,
            label="$(group_label)=$(_cz_number(eps))",
        )
    end

    Plots.savefig(plt, out_path)
    return nothing
end

"""
    plot_parameter_line(dyn_summary::DataFrame; ...) -> Nothing

Pro modely s jedním parametrem vykreslí finální hodnotu metriky podle
`epsilon`.
"""
function plot_parameter_line(
    dyn_summary::DataFrame;
    network::String,
    model::String,
    metric::Symbol,
    out_path::String,
    title::String = "",
    xlabel::String = "ε",
    ylabel::String = "",
)
    mean_col, _, _, _ = _metric_columns(metric)
    sub = _final_parameter_slice(dyn_summary, network, model)
    isempty(sub) && return nothing

    sort!(sub, :epsilon)
    mkpath(dirname(out_path))
    plt = plot(
        sub.epsilon,
        sub[!, mean_col];
        marker=:circle,
        linewidth=2,
        legend=false,
        xlabel=xlabel,
        ylabel=ylabel == "" ? String(metric) : ylabel,
        title=title,
        size=(900, 550),
        left_margin=16Plots.mm,
        right_margin=8Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )
    Plots.savefig(plt, out_path)
    return nothing
end

function _final_metric_by_epsilon(sub::DataFrame, metric::Symbol)
    mean_col, _, _, _ = _metric_columns(metric)
    isempty(sub) && return sub

    final_t = maximum(sub.t)
    final = sub[sub.t .== final_t, :]
    return combine(groupby(final, :epsilon), mean_col => mean => :value)
end

"""
    plot_deffuant_hk_final_comparison(dyn_summary::DataFrame; ...) -> Nothing

Porovná finální hodnotu metriky mezi Deffuantovým modelem a HK modelem.
HK má jednu křivku podle `epsilon`; Deffuant je zprůměrovaný přes hodnoty
`mu`, aby bylo porovnání jednorozměrné.
"""
function plot_deffuant_hk_final_comparison(
    dyn_summary::DataFrame;
    network::String,
    metric::Symbol,
    out_path::String,
    title::String = "",
    xlabel::String = "ε",
    ylabel::String = "",
)
    deffuant = _final_metric_by_epsilon(_subset(dyn_summary, network, "deffuant"), metric)
    hk = _final_metric_by_epsilon(_subset(dyn_summary, network, "hk"), metric)
    isempty(deffuant) && isempty(hk) && return nothing

    mkpath(dirname(out_path))
    plt = plot(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel == "" ? String(metric) : ylabel,
        legend=:outertopright,
        size=(1050, 650),
        left_margin=16Plots.mm,
        right_margin=8Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )

    if !isempty(deffuant)
        sort!(deffuant, :epsilon)
        plot!(
            plt,
            deffuant.epsilon,
            deffuant.value;
            label="Deffuantův model (průměr přes μ)",
            marker=:circle,
            linestyle=:solid,
            linewidth=2.4,
        )
    end

    if !isempty(hk)
        sort!(hk, :epsilon)
        plot!(
            plt,
            hk.epsilon,
            hk.value;
            label="HK",
            marker=:diamond,
            linestyle=:dash,
            linewidth=2.4,
        )
    end

    savefig(plt, out_path)
    return nothing
end

function _dynamic_metric_by_epsilon(sub::DataFrame, metric::Symbol)
    mean_col, _, _, _ = _metric_columns(metric)
    isempty(sub) && return sub
    return combine(groupby(sub, [:epsilon, :t]), mean_col => mean => :value)
end

"""
    plot_deffuant_hk_dynamic_comparison(dyn_summary::DataFrame; ...) -> Nothing

Porovná časový průběh metriky mezi Deffuantovým modelem a HK modelem.
Pro každou hodnotu `epsilon` kreslí dvojici křivek; Deffuant je
zprůměrovaný přes hodnoty `mu`.
"""
function plot_deffuant_hk_dynamic_comparison(
    dyn_summary::DataFrame;
    network::String,
    metric::Symbol,
    out_path::String,
    title::String = "",
    xlabel::String = "t",
    ylabel::String = "",
)
    deffuant = _dynamic_metric_by_epsilon(_subset(dyn_summary, network, "deffuant"), metric)
    hk = _dynamic_metric_by_epsilon(_subset(dyn_summary, network, "hk"), metric)
    isempty(deffuant) && isempty(hk) && return nothing

    eps_values = sort(unique(vcat(
        isempty(deffuant) ? Float64[] : Float64.(deffuant.epsilon),
        isempty(hk) ? Float64[] : Float64.(hk.epsilon),
    )))

    mkpath(dirname(out_path))
    plt = plot(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel == "" ? String(metric) : ylabel,
        legend=:outertopright,
        size=(1100, 650),
        left_margin=16Plots.mm,
        right_margin=8Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )

    for (i, eps) in enumerate(eps_values)
        dsub = deffuant[Float64.(deffuant.epsilon) .== eps, :]
        hsub = hk[Float64.(hk.epsilon) .== eps, :]

        if !isempty(dsub)
            sort!(dsub, :t)
            plot!(
                plt,
                dsub.t,
                dsub.value;
                label="Deffuantův model ε=$(_cz_number(eps))",
                linewidth=2,
                linestyle=SERIES_LINESTYLES[mod1(i, length(SERIES_LINESTYLES))],
            )
        end

        if !isempty(hsub)
            sort!(hsub, :t)
            plot!(
                plt,
                hsub.t,
                hsub.value;
                label="HK ε=$(_cz_number(eps))",
                linewidth=2,
                linestyle=SERIES_LINESTYLES[mod1(i + 1, length(SERIES_LINESTYLES))],
            )
        end
    end

    savefig(plt, out_path)
    return nothing
end

"""
    plot_models_parameter_surfaces(dyn_summary; ...)

Porovná finální parametrické plochy dvou modelů nad společnými osami
`ε × μ`. Plochy používají odlišné barevné škály a průhlednost, aby bylo
vidět jejich vzájemnou polohu i v místech, kde se překrývají.
"""
function plot_models_parameter_surfaces(
    dyn_summary::DataFrame;
    network::String,
    metric::Symbol,
    first_model::String="deffuant",
    second_model::String="deffuant_hrdcr",
    out_path::String,
    title::String="",
    xlabel::String="μ",
    ylabel::String="ε",
    zlabel::String="",
)
    eps1, mu1, z1 = _parameter_matrix(dyn_summary, network, first_model, metric)
    eps2, mu2, z2 = _parameter_matrix(dyn_summary, network, second_model, metric)
    (isempty(eps1) || isempty(eps2)) && return nothing

    mkpath(dirname(out_path))
    plt = Plots.surface(
        mu1,
        eps1,
        z1;
        color=:blues,
        alpha=0.72,
        colorbar=false,
        label=_model_label(first_model),
        xlabel=xlabel,
        ylabel=ylabel,
        zlabel=zlabel == "" ? _metric_label(metric) : zlabel,
        title=title,
        camera=(42, 28),
        size=(1050, 750),
        left_margin=12Plots.mm,
        right_margin=10Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )
    Plots.surface!(
        plt,
        mu2,
        eps2,
        z2;
        color=:reds,
        alpha=0.62,
        colorbar=false,
        label=_model_label(second_model),
    )

    # Samostatné pomocné série zajistí čitelnou legendu i u 3D povrchů.
    Plots.plot!(plt, [NaN], [NaN], [NaN]; color=:steelblue, linewidth=5, label=_model_label(first_model))
    Plots.plot!(plt, [NaN], [NaN], [NaN]; color=:firebrick, linewidth=5, label=_model_label(second_model))

    Plots.savefig(plt, out_path)
    return nothing
end

function _model_label(model::AbstractString)
    if model == "deffuant"
        return "Deffuantův model"
    elseif model == "hk"
        return "HK"
    elseif model == "deffuant_hrdcr"
        return "HRDCR"
    elseif model == "deffuant_anchored"
        return "Anchored Deffuant"
    elseif model == "deffuant_partisan"
        return "Partisan Deffuant"
    end
    return String(model)
end

"""
    plot_models_final_comparison(dyn_summary::DataFrame; ...) -> Nothing

Porovná finální hodnoty metriky pro zadané modely v jednom grafu. Modely
s parametrem `mu` jsou pro každou hodnotu `epsilon` zprůměrované přes `mu`.
"""
function plot_models_final_comparison(
    dyn_summary::DataFrame;
    network::String,
    metric::Symbol,
    models::Vector{String},
    out_path::String,
    title::String = "",
    xlabel::String = "ε",
    ylabel::String = "",
)
    mkpath(dirname(out_path))
    plt = plot(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel == "" ? String(metric) : ylabel,
        legend=:outertopright,
        size=(1000, 600),
        left_margin=16Plots.mm,
        right_margin=8Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )

    drew_anything = false
    for (i, model) in enumerate(models)
        model_data = _final_metric_by_epsilon(_subset(dyn_summary, network, model), metric)
        isempty(model_data) && continue

        sort!(model_data, :epsilon)
        plot!(
            plt,
            model_data.epsilon,
            model_data.value;
            label=_model_label(model),
            marker=SERIES_MARKERS[mod1(i, length(SERIES_MARKERS))],
            linestyle=SERIES_LINESTYLES[mod1(i, length(SERIES_LINESTYLES))],
            linewidth=2.4,
        )
        drew_anything = true
    end

    drew_anything || return nothing
    savefig(plt, out_path)
    return nothing
end

"""
    plot_models_dynamic_comparison(dyn_summary::DataFrame; ...) -> Nothing

Porovná časový vývoj metriky pro zadané modely v jednom grafu. Aby graf
zůstal čitelný, kreslí pouze zadanou hodnotu `epsilon`; modely s parametrem
`mu` jsou pro tuto hodnotu zprůměrované přes `mu`.
"""
function plot_models_dynamic_comparison(
    dyn_summary::DataFrame;
    network::String,
    metric::Symbol,
    models::Vector{String},
    epsilon::Real,
    out_path::String,
    title::String = "",
    xlabel::String = "t",
    ylabel::String = "",
)
    mkpath(dirname(out_path))
    plt = plot(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel == "" ? String(metric) : ylabel,
        legend=:outertopright,
        size=(1100, 650),
        left_margin=16Plots.mm,
        right_margin=8Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )

    drew_anything = false
    for (i, model) in enumerate(models)
        model_data = _dynamic_metric_by_epsilon(_subset(dyn_summary, network, model), metric)
        isempty(model_data) && continue

        sub = model_data[Float64.(model_data.epsilon) .== Float64(epsilon), :]
        isempty(sub) && continue

        sort!(sub, :t)
        plot!(
            plt,
            sub.t,
            sub.value;
            label=_model_label(model),
            linewidth=2.4,
            linestyle=SERIES_LINESTYLES[mod1(i, length(SERIES_LINESTYLES))],
        )
        drew_anything = true
    end

    drew_anything || return nothing
    savefig(plt, out_path)
    return nothing
end

function plot_opinion_series(
    series::Vector;
    out_path::String,
    title::String = "",
    xlabel::String = "time",
    ylabel::String = "opinion",
    iterations::Int,
)
    mkpath(dirname(out_path))

    steps = length(series)

    # osa x odpovídá skutečným iteracím, ne jen indexům snapshotů
    xs = range(0, iterations; length=steps)

    p = plot(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        legend=false,
    )

    n = length(series[1])

    for i in 1:n
        ys = [series[t][i] for t in 1:steps]

        plot!(
            p,
            xs,
            ys,
            linewidth=0.4,
            alpha=0.15,
        )
    end

    savefig(p, out_path)
    return nothing
end

"""
    plot_nominate_dim1_distribution(; csv_path, out_path, bins=50) -> Nothing

Načte `HSall_members.csv` a vykreslí histogram hodnot ve sloupci
`nominate_dim1`. Prázdné hodnoty v CSV ignoruje.
"""
function plot_nominate_dim1_distribution(;
    csv_path::AbstractString = joinpath(PROJECT_ROOT, "data", "opinions", "HSall_members.csv"),
    out_path::AbstractString = joinpath(PROJECT_ROOT, "out", "plots", "nominate_dim1_distribution.png"),
    bins::Int = 50,
)
    df = CSV.read(csv_path, DataFrame)

    :nominate_dim1 in propertynames(df) ||
        throw(ArgumentError("CSV musí obsahovat sloupec 'nominate_dim1'"))

    values = collect(skipmissing(df.nominate_dim1))
    isempty(values) && throw(ArgumentError("Sloupec 'nominate_dim1' neobsahuje žádné hodnoty k vykreslení"))

    mkpath(dirname(out_path))

    p = histogram(
        values;
        bins=bins,
        normalize=:pdf,
        legend=false,
        xlabel="nominate_dim1",
        ylabel="hustota",
        title="Rozložení hodnot nominate_dim1",
        color=:steelblue,
        alpha=0.75,
    )

    savefig(p, out_path)
    return nothing
end

end # module
