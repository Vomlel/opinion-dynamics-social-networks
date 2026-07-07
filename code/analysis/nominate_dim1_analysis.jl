ENV["GKSwstype"] = "100"

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CODE_DIR = joinpath(PROJECT_ROOT, "code")

include(joinpath(CODE_DIR, "metrics", "dynamic_characteristics.jl"))
include(joinpath(CODE_DIR, "reporting", "plots.jl"))

using CSV
using DataFrames
using Statistics
using Base.Filesystem: dirname, mkpath
using .DynamicCharacteristics
using .Plotting

const Plots = Plotting.Plots

function _normalize_nominate_dim1(values)
    return (Float64.(values) .+ 1.0) ./ 2.0
end

function _valid_nominate_rows(df::DataFrame)
    :nominate_dim1 in propertynames(df) ||
        throw(ArgumentError("CSV musi obsahovat sloupec nominate_dim1"))
    :congress in propertynames(df) ||
        throw(ArgumentError("CSV musi obsahovat sloupec congress"))
    :chamber in propertynames(df) ||
        throw(ArgumentError("CSV musi obsahovat sloupec chamber"))

    clean = dropmissing(df, [:nominate_dim1, :congress, :chamber])
    clean.nominate_dim1 = Float64.(clean.nominate_dim1)
    clean.congress = Int.(clean.congress)
    return clean
end

function _metrics_for_values(values)
    x = sort(_normalize_nominate_dim1(values))
    return (
        n = length(x),
        clusters = count_of_clusters_dbscan_auto(x),
        pad = pairwise_average_distance(x),
        polarization = polarization_der(x),
        mean_nominate_dim1 = mean(Float64.(values)),
        median_nominate_dim1 = median(Float64.(values)),
    )
end

function nominate_dim1_dynamic_summary(df::DataFrame; min_group_size::Int=20)
    clean = _valid_nominate_rows(df)

    rows = DataFrame(
        congress=Int[],
        chamber=String[],
        n=Int[],
        clusters=Int[],
        pad=Float64[],
        polarization=Float64[],
        mean_nominate_dim1=Float64[],
        median_nominate_dim1=Float64[],
    )

    for sdf in groupby(clean, [:congress, :chamber])
        nrow(sdf) < min_group_size && continue
        metrics = _metrics_for_values(sdf.nominate_dim1)
        push!(rows, (
            first(sdf.congress),
            String(first(sdf.chamber)),
            metrics.n,
            metrics.clusters,
            metrics.pad,
            metrics.polarization,
            metrics.mean_nominate_dim1,
            metrics.median_nominate_dim1,
        ))
    end

    for sdf in groupby(clean, :congress)
        nrow(sdf) < min_group_size && continue
        metrics = _metrics_for_values(sdf.nominate_dim1)
        push!(rows, (
            first(sdf.congress),
            "All",
            metrics.n,
            metrics.clusters,
            metrics.pad,
            metrics.polarization,
            metrics.mean_nominate_dim1,
            metrics.median_nominate_dim1,
        ))
    end

    sort!(rows, [:chamber, :congress])
    return rows
end

function _plot_distribution(clean::DataFrame, out_dir::AbstractString)
    mkpath(out_dir)

    values = Float64.(clean.nominate_dim1)
    p = Plots.histogram(
        values;
        bins=60,
        normalize=:pdf,
        legend=false,
        xlabel="nominate_dim1",
        ylabel="hustota",
        title="Rozlozeni nominate_dim1",
        color=:steelblue,
        alpha=0.75,
    )
    Plots.savefig(p, joinpath(out_dir, "nominate_dim1_distribution.png"))

    p_chamber = Plots.plot(
        xlabel="nominate_dim1",
        ylabel="hustota",
        title="Rozlozeni nominate_dim1 podle komory",
        legend=:topright,
    )
    for chamber in ["House", "Senate"]
        sub = clean[clean.chamber .== chamber, :]
        isempty(sub) && continue
        Plots.histogram!(
            p_chamber,
            Float64.(sub.nominate_dim1);
            bins=60,
            normalize=:pdf,
            alpha=0.45,
            label=chamber,
        )
    end
    Plots.savefig(p_chamber, joinpath(out_dir, "nominate_dim1_distribution_by_chamber.png"))
end

function _plot_metric(summary::DataFrame, metric::Symbol, out_path::AbstractString; ylabel::String=String(metric))
    p = Plots.plot(
        xlabel="congress",
        ylabel=ylabel,
        title="$(ylabel) podle kongresu",
        legend=:outertopright,
    )

    for chamber in ["All", "House", "Senate"]
        sub = summary[summary.chamber .== chamber, :]
        isempty(sub) && continue
        sort!(sub, :congress)
        Plots.plot!(
            p,
            sub.congress,
            sub[!, metric];
            marker=:circle,
            linewidth=2,
            label=chamber,
        )
    end

    Plots.savefig(p, out_path)
end

function _plot_dynamic_characteristics(summary::DataFrame, out_dir::AbstractString)
    mkpath(out_dir)

    _plot_metric(summary, :clusters, joinpath(out_dir, "nominate_dim1_clusters_by_congress.png"); ylabel="pocet clusteru")
    _plot_metric(summary, :pad, joinpath(out_dir, "nominate_dim1_pad_by_congress.png"); ylabel="PAD")
    _plot_metric(summary, :polarization, joinpath(out_dir, "nominate_dim1_polarization_by_congress.png"); ylabel="polarizace")
    _plot_metric(summary, :mean_nominate_dim1, joinpath(out_dir, "nominate_dim1_mean_by_congress.png"); ylabel="prumer nominate_dim1")

    all_rows = summary[summary.chamber .== "All", :]
    sort!(all_rows, :congress)
    if !isempty(all_rows)
        p1 = Plots.plot(all_rows.congress, all_rows.clusters; marker=:circle, linewidth=2, label=false, ylabel="pocet clusteru")
        p2 = Plots.plot(all_rows.congress, all_rows.pad; marker=:circle, linewidth=2, label=false, ylabel="PAD")
        p3 = Plots.plot(all_rows.congress, all_rows.polarization; marker=:circle, linewidth=2, label=false, ylabel="polarizace", xlabel="congress")
        p = Plots.plot(p1, p2, p3; layout=(3, 1), size=(1000, 900), title="Dynamicke charakteristiky nominate_dim1")
        Plots.savefig(p, joinpath(out_dir, "nominate_dim1_dynamic_characteristics.png"))
    end
end

function _plot_last_congress_distribution(clean::DataFrame, out_dir::AbstractString; bins::Int=180)
    mkpath(out_dir)

    last_congress = maximum(clean.congress)
    last = clean[(clean.congress .== last_congress) .& (clean.chamber .!= "President"), :]
    isempty(last) && return nothing

    values = Float64.(last.nominate_dim1)
    clusters = count_of_clusters_dbscan_auto(sort(_normalize_nominate_dim1(values)))

    p = Plots.histogram(
        values;
        bins=bins,
        normalize=false,
        legend=false,
        xlabel="nominate_dim1",
        ylabel="pocet",
        title="Rozlozeni nominate_dim1, congress $(last_congress) (DBSCAN clusters=$(clusters))",
        color=:steelblue,
        alpha=0.78,
        size=(1200, 650),
    )
    Plots.savefig(p, joinpath(out_dir, "nominate_dim1_last_congress_$(last_congress)_fine_histogram.png"))

    p_chamber = Plots.plot(
        xlabel="nominate_dim1",
        ylabel="pocet",
        title="Rozlozeni nominate_dim1 podle komory, congress $(last_congress)",
        legend=:topright,
        size=(1200, 650),
    )
    for chamber in ["House", "Senate"]
        sub = last[last.chamber .== chamber, :]
        isempty(sub) && continue
        Plots.histogram!(
            p_chamber,
            Float64.(sub.nominate_dim1);
            bins=bins,
            normalize=false,
            alpha=0.45,
            label=chamber,
        )
    end
    Plots.savefig(p_chamber, joinpath(out_dir, "nominate_dim1_last_congress_$(last_congress)_fine_histogram_by_chamber.png"))

    sorted = sort(last, :nominate_dim1)
    p_sorted = Plots.scatter(
        1:nrow(sorted),
        sorted.nominate_dim1;
        group=sorted.chamber,
        xlabel="poradi podle nominate_dim1",
        ylabel="nominate_dim1",
        title="Serazene nazory, congress $(last_congress)",
        markerstrokewidth=0,
        markersize=3.0,
        alpha=0.8,
        legend=:topleft,
        size=(1200, 650),
    )
    Plots.savefig(p_sorted, joinpath(out_dir, "nominate_dim1_last_congress_$(last_congress)_sorted_values.png"))

    return last_congress
end

function _plot_last_congress_chamber_distribution(
    clean::DataFrame,
    out_dir::AbstractString,
    chamber::AbstractString;
    bins::Int=120,
)
    mkpath(out_dir)

    last_congress = maximum(clean.congress)
    last = clean[(clean.congress .== last_congress) .& (clean.chamber .== chamber), :]
    isempty(last) && return nothing

    values = Float64.(last.nominate_dim1)
    clusters = count_of_clusters_dbscan_auto(sort(_normalize_nominate_dim1(values)))
    chamber_label = lowercase(String(chamber))

    p = Plots.histogram(
        values;
        bins=bins,
        normalize=false,
        legend=false,
        xlabel="nominate_dim1",
        ylabel="pocet",
        title="$(chamber), congress $(last_congress): rozlozeni nominate_dim1 (clusters=$(clusters))",
        color=:darkorange,
        alpha=0.8,
        size=(1200, 650),
    )
    Plots.savefig(p, joinpath(out_dir, "nominate_dim1_last_congress_$(last_congress)_$(chamber_label)_fine_histogram.png"))

    sorted_values = sort(values)
    p_sorted = Plots.scatter(
        1:length(sorted_values),
        sorted_values;
        xlabel="poradi podle nominate_dim1",
        ylabel="nominate_dim1",
        title="$(chamber), congress $(last_congress): serazene nazory",
        markerstrokewidth=0,
        markersize=4.0,
        alpha=0.85,
        legend=false,
        size=(1200, 650),
    )
    Plots.savefig(p_sorted, joinpath(out_dir, "nominate_dim1_last_congress_$(last_congress)_$(chamber_label)_sorted_values.png"))

    return last_congress
end

function plot_senate_119_static_opinion_series(;
    csv_path::AbstractString = joinpath(PROJECT_ROOT, "data", "konkresy", "HSall_members.csv"),
    out_path::AbstractString = joinpath(PROJECT_ROOT, "out", "konkresy_senate_119", "plots", "senate_119_nominate_dim1_static_opinion_series.png"),
    congress::Int = 119,
    fake_time_steps::Int = 100,
    normalize::Bool = false,
    color_clusters::Bool = false,
)
    df = CSV.read(csv_path, DataFrame)
    clean = _valid_nominate_rows(df)
    senate = clean[(clean.congress .== congress) .& (clean.chamber .== "Senate"), :]
    isempty(senate) && throw(ArgumentError("Nenalezeny zadne zaznamy pro Senate v congress=$congress"))

    mkpath(dirname(out_path))

    values = sort(Float64.(senate.nominate_dim1))
    if normalize
        values = _normalize_nominate_dim1(values)
    end
    xs = collect(0:fake_time_steps)
    ylabel = normalize ? "normalizovany nominate_dim1" : "nominate_dim1"
    ylims = normalize ? (0.0, 1.0) : (-1.0, 1.0)
    title_suffix = normalize ? "normalizovane" : "surove"
    cluster_labels = color_clusters ? DynamicCharacteristics.dbscan_1d(values; eps=0.01, min_pts=max(2, ceil(Int, length(values) / 100))) : fill(1, length(values))
    cluster_count = maximum(cluster_labels)
    palette = color_clusters ? Plots.palette(:tab10, max(cluster_count, 1)) : [:steelblue]

    p = Plots.plot(
        title=color_clusters ?
              "Senate $(congress): nominate_dim1 jako staticke nazory ($(title_suffix), clusters=$(cluster_count))" :
              "Senate $(congress): nominate_dim1 jako staticke nazory ($(title_suffix))",
        xlabel="cas",
        ylabel=ylabel,
        legend=color_clusters ? :outertopright : false,
        size=(1200, 700),
        ylims=ylims,
    )

    seen_clusters = Set{Int}()
    for (value, label) in zip(values, cluster_labels)
        color = label > 0 ? palette[label] : :gray
        line_label = if color_clusters && label > 0 && !(label in seen_clusters)
            push!(seen_clusters, label)
            "cluster $(label)"
        elseif color_clusters && label <= 0 && !(label in seen_clusters)
            push!(seen_clusters, label)
            "noise"
        else
            false
        end

        Plots.plot!(
            p,
            xs,
            fill(value, length(xs));
            linewidth=1.2,
            alpha=0.75,
            color=color,
            label=line_label,
        )
    end

    Plots.savefig(p, out_path)
    return out_path
end

function run_nominate_dim1_analysis(;
    csv_path::AbstractString = joinpath(PROJECT_ROOT, "data", "konkresy", "HSall_members.csv"),
    out_dir::AbstractString = joinpath(PROJECT_ROOT, "out", "nominate_dim1_download"),
    min_group_size::Int = 20,
)
    df = CSV.read(csv_path, DataFrame)
    clean = _valid_nominate_rows(df)

    plots_dir = joinpath(out_dir, "plots")
    mkpath(out_dir)
    mkpath(plots_dir)

    summary = nominate_dim1_dynamic_summary(clean; min_group_size=min_group_size)
    CSV.write(joinpath(out_dir, "nominate_dim1_dynamic_summary.csv"), summary)

    _plot_distribution(clean, plots_dir)
    _plot_dynamic_characteristics(summary, plots_dir)
    _plot_last_congress_distribution(clean, plots_dir)
    _plot_last_congress_chamber_distribution(clean, plots_dir, "Senate")
    plot_senate_119_static_opinion_series(
        csv_path=csv_path,
        out_path=joinpath(plots_dir, "senate_119_nominate_dim1_static_opinion_series.png"),
    )

    println("Saved outputs to: $out_dir")
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    csv_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(PROJECT_ROOT, "data", "konkresy", "HSall_members.csv")
    out_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(PROJECT_ROOT, "out", "nominate_dim1_download")
    run_nominate_dim1_analysis(csv_path=csv_path, out_dir=out_dir)
end
