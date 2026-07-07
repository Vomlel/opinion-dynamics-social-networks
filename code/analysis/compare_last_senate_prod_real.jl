ENV["GKSwstype"] = "100"

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "code", "reporting", "plots.jl"))

using CSV
using DataFrames
using Statistics
using .Plotting

const Plots = Plotting.Plots

const MODEL_LABELS = Dict(
    "deffuant" => "Deffuant",
    "hk" => "HK",
    "deffuant_hrdcr" => "HRDCR",
)

const MODEL_LINESTYLES = Dict("deffuant" => :solid, "hk" => :dash, "deffuant_hrdcr" => :dot)
const MODEL_MARKERS = Dict("deffuant" => :circle, "hk" => :diamond, "deffuant_hrdcr" => :utriangle)

const METRIC_LABELS = Dict(
    :clusters => "počet clusterů",
    :pad => "průměrná párová vzdálenost",
    :polarization => "polarizace",
)

function _metric_mean_column(metric::Symbol)
    if metric == :clusters
        return :clusters_mean
    elseif metric == :pad
        return :pad_mean
    elseif metric == :polarization
        return :polarization_mean
    end
    throw(ArgumentError("Unknown metric: $metric"))
end

function _senate_reference(nominate_summary::DataFrame)
    senate = nominate_summary[nominate_summary.chamber .== "Senate", :]
    isempty(senate) && throw(ArgumentError("V nominate summary chybi chamber=Senate"))
    last_congress = maximum(senate.congress)
    ref = senate[senate.congress .== last_congress, :]
    nrow(ref) == 1 || throw(ArgumentError("Pro posledni Senate congress cekam prave jeden radek"))
    return ref[1, :]
end

function _final_prod_by_epsilon(prod_summary::DataFrame)
    final_t = maximum(prod_summary.t)
    final = prod_summary[prod_summary.t .== final_t, :]

    return combine(
        groupby(final, [:network, :model, :epsilon]),
        :clusters_mean => mean => :clusters,
        :pad_mean => mean => :pad,
        :polarization_mean => mean => :polarization,
    )
end

function _comparison_table(prod_final::DataFrame, senate_ref)
    rows = DataFrame(
        network=String[],
        model=String[],
        epsilon=Float64[],
        metric=String[],
        prod_value=Float64[],
        senate_value=Float64[],
        abs_diff=Float64[],
        rel_diff=Float64[],
    )

    for r in eachrow(prod_final)
        for metric in [:clusters, :pad, :polarization]
            prod_value = Float64(r[metric])
            senate_value = Float64(senate_ref[metric])
            abs_diff = abs(prod_value - senate_value)
            rel_diff = senate_value == 0.0 ? abs_diff : abs_diff / abs(senate_value)
            push!(rows, (
                String(r.network),
                String(r.model),
                Float64(r.epsilon),
                String(metric),
                prod_value,
                senate_value,
                abs_diff,
                rel_diff,
            ))
        end
    end

    sort!(rows, [:metric, :network, :model, :epsilon])
    return rows
end

function _distance_table(comparison::DataFrame)
    return combine(
        groupby(comparison, [:network, :model, :epsilon]),
        :rel_diff => mean => :mean_relative_distance,
        :abs_diff => mean => :mean_absolute_distance,
    ) |> df -> sort(df, :mean_relative_distance)
end

function _plot_metric(prod_final::DataFrame, senate_ref, metric::Symbol, out_path::AbstractString)
    networks = sort(unique(String.(prod_final.network)))
    models = ["deffuant", "hk", "deffuant_hrdcr"]
    senate_value = Float64(senate_ref[metric])

    panels = []
    for network in networks
        sub_network = prod_final[prod_final.network .== network, :]
        p = Plots.plot(
            xlabel="ε",
            ylabel=get(METRIC_LABELS, metric, String(metric)),
            title=network == "csv" ? "reálná síť" : uppercase(network),
            legend=:outertopright,
            left_margin=14Plots.mm,
            right_margin=8Plots.mm,
            top_margin=6Plots.mm,
            bottom_margin=8Plots.mm,
        )

        for model in models
            sub = sub_network[sub_network.model .== model, :]
            isempty(sub) && continue
            sort!(sub, :epsilon)
            Plots.plot!(
                p,
                sub.epsilon,
                sub[!, metric];
                marker=get(MODEL_MARKERS, model, :circle),
                linestyle=get(MODEL_LINESTYLES, model, :solid),
                linewidth=2.2,
                label=get(MODEL_LABELS, model, model),
            )
        end

        eps_values = sort(unique(Float64.(sub_network.epsilon)))
        if !isempty(eps_values)
            Plots.plot!(
                p,
                eps_values,
                fill(senate_value, length(eps_values));
                color=:black,
                linestyle=:dash,
                linewidth=2.4,
                label="Senát 119",
            )
        end

        push!(panels, p)
    end

    plt = Plots.plot(
        panels...;
        layout=(length(panels), 1),
        size=(1100, 330 * length(panels)),
        plot_title="Finální výsledky modelů a reference 119. Senátu: $(get(METRIC_LABELS, metric, String(metric)))",
    )
    Plots.savefig(plt, out_path)
end

function _plot_distance(distance::DataFrame, out_path::AbstractString; top_n::Int=12)
    top = first(distance, min(top_n, nrow(distance)))
    labels = ["$(r.network) | $(get(MODEL_LABELS, String(r.model), String(r.model))) | eps=$(r.epsilon)" for r in eachrow(top)]
    y = reverse(1:nrow(top))
    x = top.mean_relative_distance

    p = Plots.scatter(
        x,
        y;
        legend=false,
        xlabel="prumerna relativni odchylka",
        ylabel="",
        title="Finální výsledky nejbližší 119. Senátu",
        markersize=7,
        markerstrokewidth=0,
        color=:steelblue,
        yticks=(y, labels),
        yflip=false,
        left_margin=18Plots.mm,
        right_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
        size=(1300, 760),
    )
    for (xi, yi) in zip(x, y)
        Plots.annotate!(p, xi, yi, Plots.text("  $(round(xi, digits=3))", 9, :left))
    end
    Plots.xlims!(p, 0, maximum(x) * 1.2)
    Plots.savefig(p, out_path)
end

function run_last_senate_prod_real_comparison(;
    prod_summary_path::AbstractString = joinpath(PROJECT_ROOT, "out", "prod_real_run_1", "aggregated", "dynamic_summary.csv"),
    nominate_summary_path::AbstractString = joinpath(PROJECT_ROOT, "out", "nominate_dim1_download", "nominate_dim1_dynamic_summary.csv"),
    out_dir::AbstractString = joinpath(PROJECT_ROOT, "out", "last_senate_vs_prod_real_run_1"),
)
    mkpath(out_dir)
    plots_dir = joinpath(out_dir, "plots")
    mkpath(plots_dir)

    prod_summary = CSV.read(prod_summary_path, DataFrame)
    nominate_summary = CSV.read(nominate_summary_path, DataFrame)
    senate_ref = _senate_reference(nominate_summary)
    prod_final = _final_prod_by_epsilon(prod_summary)
    comparison = _comparison_table(prod_final, senate_ref)
    distance = _distance_table(comparison)

    CSV.write(joinpath(out_dir, "prod_real_final_by_epsilon.csv"), prod_final)
    CSV.write(joinpath(out_dir, "last_senate_reference.csv"), DataFrame(senate_ref))
    CSV.write(joinpath(out_dir, "comparison_metrics.csv"), comparison)
    CSV.write(joinpath(out_dir, "comparison_distance_ranked.csv"), distance)

    for metric in [:clusters, :pad, :polarization]
        _plot_metric(prod_final, senate_ref, metric, joinpath(plots_dir, "prod_real_vs_last_senate_$(metric).png"))
    end
    _plot_distance(distance, joinpath(plots_dir, "prod_real_vs_last_senate_closest_settings.png"))

    println("Saved outputs to: $out_dir")
    return (prod_final=prod_final, senate_ref=senate_ref, comparison=comparison, distance=distance)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_last_senate_prod_real_comparison()
end
