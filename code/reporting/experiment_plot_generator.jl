module ExperimentPlotGenerator

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using Base.Filesystem: mkpath
using DataFrames
using ..Plotting
using ..ThesisFigures

export generate_first_run_opinion_plot,
       generate_standard_experiment_plots,
       generate_thesis_experiment_plots

function _metric_symbol(metric::AbstractString)
    metric == "clusters" && return :clusters
    metric == "pad" && return :pad
    metric == "polarization" && return :polarization
    return nothing
end

function _metric_label(metric::Symbol)
    metric == :clusters && return "počet clusterů"
    metric == :pad && return "průměrná párová vzdálenost"
    metric == :polarization && return "polarizace"
    return String(metric)
end

function _model_label(model::AbstractString)
    model == "deffuant" && return "Deffuantův model"
    model == "hk" && return "HK"
    model == "deffuant_hrdcr" && return "HRDCR"
    model == "deffuant_anchored" && return "ukotvený Deffuant"
    model == "deffuant_partisan" && return "stranický Deffuant"
    return String(model)
end

function _network_label(network::AbstractString)
    network in ("csv", "real") && return "reálná síť"
    return uppercase(String(network))
end

function _network_locative(network::AbstractString)
    network in ("csv", "real") && return "reálné síti"
    network == "er" && return "síti ER"
    network == "ba" && return "síti BA"
    return "síti $(uppercase(String(network)))"
end

function _data_network_name(dynamic_summary::DataFrame, configured_network::AbstractString)
    configured = String(configured_network)
    available = Set(String.(unique(dynamic_summary.network)))
    configured in available && return configured
    configured == "real" && "csv" in available && return "csv"
    configured == "csv" && "real" in available && return "real"
    return configured
end

function _model_iterations(cfg, model::AbstractString)
    if model == "deffuant"
        return cfg.models.deffuant.iterations
    elseif model == "hk"
        return cfg.models.hk.iterations
    elseif model == "deffuant_hrdcr"
        return cfg.models.deffuant_hrdcr.iterations
    elseif model == "deffuant_anchored"
        return cfg.models.deffuant_anchored.iterations
    elseif model == "deffuant_partisan"
        return cfg.models.deffuant_partisan.iterations
    end
    throw(ArgumentError("$model musí mít svou podmínku zde"))
end

function _epsilon_label(epsilon::Real)
    return replace(string(Float64(epsilon)), "." => "p", "-" => "m")
end

_cz_number(value::Real; digits::Int=3) = replace(string(round(Float64(value), digits=digits)), "." => ",")

function generate_first_run_opinion_plot(
    series::Vector;
    cfg,
    model::String,
    network::String,
    suffix::String,
)
    getproperty(cfg.outputs, :plot_first_run_opinions) || return nothing

    plots_dir = joinpath(cfg.experiment.out_dir, "plots", "first_runs")
    mkpath(plots_dir)

    outp = joinpath(plots_dir, "$(model)_$(network)_$(suffix)_first_run_opinions.png")
    iterations = _model_iterations(cfg, model)

    Plotting.plot_opinion_series(
        series;
        out_path=outp,
        title="$(_model_label(model)) na $(_network_locative(network)): vývoj názorů",
        xlabel="iterace",
        ylabel="názor",
        iterations=iterations,
    )

    return outp
end

function generate_standard_experiment_plots(dynamic_summary::DataFrame, cfg)
    isempty(cfg.outputs.plots) && return nothing

    plots_dir = joinpath(cfg.experiment.out_dir, "plots")

    for network in cfg.networks.enabled
        data_network = _data_network_name(dynamic_summary, network)
        for metric_s in cfg.outputs.plots
            metric = _metric_symbol(metric_s)
            metric === nothing && continue
            metric_label = _metric_label(metric)

            for model in cfg.models.enabled
                title = "$(_model_label(model)) na $(_network_locative(network)): $(metric_label)"
                outp = joinpath(plots_dir, "$(model)_$(data_network)_$(metric_s)_params.png")
                iterations = _model_iterations(cfg, model)

                Plotting.plot_dynamic_metric_params(
                    dynamic_summary;
                    network=data_network,
                    metric=metric,
                    model=model,
                    out_path=outp,
                    title=title,
                    xlabel="podíl celkového počtu iterací (celkem $iterations)",
                    ylabel=metric_label,
                )

                if model in ("deffuant", "deffuant_hrdcr", "deffuant_anchored", "deffuant_partisan")
                    heatmap_out = joinpath(plots_dir, "$(model)_$(data_network)_$(metric_s)_parameter_heatmap.png")
                    Plotting.plot_parameter_heatmap(
                        dynamic_summary;
                        network=data_network,
                        model=model,
                        metric=metric,
                        out_path=heatmap_out,
                        title="$(_model_label(model)) na $(_network_locative(network)): finální $(metric_label)",
                        ylabel="ε",
                        xlabel="μ",
                    )

                    surface_out = joinpath(plots_dir, "$(model)_$(data_network)_$(metric_s)_parameter_surface.png")
                    Plotting.plot_parameter_surface(
                        dynamic_summary;
                        network=data_network,
                        model=model,
                        metric=metric,
                        out_path=surface_out,
                        title="$(_model_label(model)) na $(_network_locative(network)): finální $(metric_label)",
                        ylabel="ε",
                        xlabel="μ",
                        zlabel=metric_label,
                    )

                    contour_out = joinpath(plots_dir, "$(model)_$(data_network)_$(metric_s)_parameter_contour.png")
                    Plotting.plot_parameter_contour(
                        dynamic_summary;
                        network=data_network,
                        model=model,
                        metric=metric,
                        out_path=contour_out,
                        title="$(_model_label(model)) na $(_network_locative(network)): finální $(metric_label)",
                        ylabel="ε",
                        xlabel="μ",
                    )

                    slices_out = joinpath(plots_dir, "$(model)_$(data_network)_$(metric_s)_parameter_slices.png")
                    Plotting.plot_parameter_slices(
                        dynamic_summary;
                        network=data_network,
                        model=model,
                        metric=metric,
                        out_path=slices_out,
                        title="$(_model_label(model)) na $(_network_locative(network)): finální $(metric_label)",
                        group_label="ε",
                        xlabel="μ",
                        ylabel=metric_label,
                    )
                elseif model == "hk"
                    line_out = joinpath(plots_dir, "$(model)_$(data_network)_$(metric_s)_parameter_line.png")
                    Plotting.plot_parameter_line(
                        dynamic_summary;
                        network=data_network,
                        model=model,
                        metric=metric,
                        out_path=line_out,
                        title="$(_model_label(model)) na $(_network_locative(network)): finální $(metric_label)",
                        xlabel="ε",
                        ylabel=metric_label,
                    )
                end
            end

            if ("deffuant" in cfg.models.enabled) && ("hk" in cfg.models.enabled)
                comparison_final_out = joinpath(plots_dir, "comparison_deffuant_hk_$(data_network)_$(metric_s)_final.png")
                Plotting.plot_deffuant_hk_final_comparison(
                    dynamic_summary;
                    network=data_network,
                    metric=metric,
                    out_path=comparison_final_out,
                    title="Deffuantův model a HK na $(_network_locative(network)): finální $(metric_label)",
                    ylabel=metric_label,
                )

                comparison_dynamic_out = joinpath(plots_dir, "comparison_deffuant_hk_$(data_network)_$(metric_s)_dynamic.png")
                Plotting.plot_deffuant_hk_dynamic_comparison(
                    dynamic_summary;
                    network=data_network,
                    metric=metric,
                    out_path=comparison_dynamic_out,
                    title="Deffuantův model a HK na $(_network_locative(network)): vývoj $(metric_label)",
                    xlabel="časový krok",
                    ylabel=metric_label,
                )
            end

            if length(cfg.models.enabled) > 1
                all_models_final_out = joinpath(plots_dir, "comparison_all_models_$(data_network)_$(metric_s)_final.png")
                Plotting.plot_models_final_comparison(
                    dynamic_summary;
                    network=data_network,
                    metric=metric,
                    models=String.(cfg.models.enabled),
                    out_path=all_models_final_out,
                    title="Porovnání modelů na $(_network_locative(network)): finální $(metric_label)",
                    ylabel=metric_label,
                )

                eps_values = sort(unique(Float64.(dynamic_summary[
                    (dynamic_summary.network .== data_network) .&
                    in.(dynamic_summary.model, Ref(cfg.models.enabled)),
                    :epsilon,
                ])))

                for epsilon in eps_values
                    eps_label = _epsilon_label(epsilon)
                    all_models_dynamic_out = joinpath(plots_dir, "comparison_all_models_$(data_network)_$(metric_s)_eps_$(eps_label)_dynamic.png")
                    Plotting.plot_models_dynamic_comparison(
                        dynamic_summary;
                        network=data_network,
                        metric=metric,
                        models=String.(cfg.models.enabled),
                        epsilon=epsilon,
                        out_path=all_models_dynamic_out,
                        title="Porovnání modelů na $(_network_locative(network)): $(metric_label), ε=$(_cz_number(epsilon))",
                        xlabel="časový krok",
                        ylabel=metric_label,
                    )
                end
            end

            if ("deffuant" in cfg.models.enabled) && ("deffuant_hrdcr" in cfg.models.enabled)
                surfaces_out = joinpath(
                    plots_dir,
                    "comparison_deffuant_hrdcr_$(data_network)_$(metric_s)_parameter_surfaces.png",
                )
                Plotting.plot_models_parameter_surfaces(
                    dynamic_summary;
                    network=data_network,
                    metric=metric,
                    first_model="deffuant",
                    second_model="deffuant_hrdcr",
                    out_path=surfaces_out,
                    title="Deffuantův model a HRDCR na $(_network_locative(network)): finální $(metric_label)",
                    xlabel="μ",
                    ylabel="ε",
                    zlabel=metric_label,
                )
            end
        end
    end

    return plots_dir
end

function generate_thesis_experiment_plots(cfg)
    cfg.outputs.save_aggregated_csv || return nothing

    prod_summary_path = joinpath(cfg.experiment.out_dir, "aggregated", "dynamic_summary.csv")
    senate_summary_path = joinpath(PROJECT_ROOT, "out", "nominate_dim1_download", "nominate_dim1_dynamic_summary.csv")
    thesis_out_dir = joinpath(cfg.experiment.out_dir, "thesis_figures")

    try
        ThesisFigures.generate_thesis_figures(
            prod_summary_path=prod_summary_path,
            senate_summary_path=senate_summary_path,
            out_dir=thesis_out_dir,
        )
    catch e
        @warn "Nepodařilo se vygenerovat bakalářkové grafy" exception=(e, catch_backtrace())
    end

    return thesis_out_dir
end

end # module
