# Orchestrator
# Z korene projektu: julia --project=code code/main.jl -c config.toml

ENV["GKSwstype"] = "100"

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "core", "config.jl"))
include(joinpath(@__DIR__, "networks", "er.jl"))
include(joinpath(@__DIR__, "networks", "ba.jl"))
include(joinpath(@__DIR__, "networks", "csv.jl"))
include(joinpath(@__DIR__, "opinions", "generate_opinion.jl"))
include(joinpath(@__DIR__, "metrics", "static_characteristics.jl"))
include(joinpath(@__DIR__, "metrics", "dynamic_characteristics.jl"))
include(joinpath(@__DIR__, "models", "deffuant.jl"))
include(joinpath(@__DIR__, "models", "hk.jl"))
include(joinpath(@__DIR__, "models", "deffuant_hrdcr.jl"))
include(joinpath(@__DIR__, "models", "deffuant_anchored.jl"))
include(joinpath(@__DIR__, "models", "deffuant_partisan.jl"))
include(joinpath(@__DIR__, "reporting", "aggregator.jl"))
include(joinpath(@__DIR__, "reporting", "plots.jl"))
include(joinpath(@__DIR__, "reporting", "thesis_figures.jl"))
include(joinpath(@__DIR__, "reporting", "experiment_plot_generator.jl"))
include(joinpath(@__DIR__, "core", "helper.jl"))

using Random
using Graphs
using CSV
using DataFrames
using .Config
using .ER
using .BA
using .CSVGraph
using .GenerateOpinion
using .StaticCharacteristics
using .DynamicCharacteristics
using .Deffuant
using .HegselmannKrause
using .DeffuantHRDCR
using .DeffuantAnchored
using .DeffuantPartisan
using .Aggregator
using .ExperimentPlotGenerator
using .Helper

function _float_label(x::Real)
    s = replace(string(Float64(x)), "." => "p", "-" => "m")
    return s
end

function _param_suffix(model::String, epsilon::Float64, mu::Float64)
    if model == "hk"
        return "eps_$(_float_label(epsilon))"
    end
    return "eps_$(_float_label(epsilon))_mu_$(_float_label(mu))"
end

function _model_params(cfg, model::String)
    if model == "deffuant"
        d = cfg.models.deffuant
        d === nothing && error("Deffuant povolen ale chybi [models.deffuant]")
        return [(epsilon=eps, mu=mu) for eps in d.epsilon for mu in d.mu]
    elseif model == "hk"
        h = cfg.models.hk
        h === nothing && error("HK povolen ale chybi [models.hk]")
        return [(epsilon=eps, mu=NaN) for eps in h.epsilon]
    elseif model == "deffuant_hrdcr"
        h = cfg.models.deffuant_hrdcr
        h === nothing && error("Deffuant-HRDCR povolen ale chybi [models.deffuant_hrdcr]")
        return [(epsilon=eps, mu=mu) for eps in h.epsilon for mu in h.mu]
    elseif model == "deffuant_anchored"
        a = cfg.models.deffuant_anchored
        a === nothing && error("Deffuant-Anchored povolen, ale chybi konfigurace")
        return [(epsilon=eps, mu=mu) for eps in a.epsilon for mu in a.mu]
    elseif model == "deffuant_partisan"
        p = cfg.models.deffuant_partisan
        return [(epsilon=eps, mu=mu) for eps in p.epsilon for mu in p.mu]
    else
        error("Unknown model: $model")
    end
end

function _format_duration(seconds::Real)
    total_seconds = max(0, round(Int, Float64(seconds)))
    hours = total_seconds ÷ 3600
    minutes = (total_seconds % 3600) ÷ 60
    secs = total_seconds % 60
    return lpad(string(hours), 2, '0') * ":" *
           lpad(string(minutes), 2, '0') * ":" *
           lpad(string(secs), 2, '0')
end

function main()
    # Vsechny relativni vystupni cesty jsou stabilne vztažene ke koreni projektu.
    cd(PROJECT_ROOT)
    # konfigurace
    config_arg = parse_args(ARGS)
    config_path = isabspath(config_arg) ? config_arg : joinpath(@__DIR__, config_arg)
    cfg = load_config(config_path)
    ensure_dirs(cfg.experiment.out_dir) # vytvor out slozky

    # uloz kopii konfiguracniho souboru
    cp = joinpath(cfg.experiment.out_dir, "config_used.toml")
    try
        write(cp, read(config_path, String))
    catch
        # ignore
    end

    # nacti real graf ze zdrojovych CSV souboru, pokud je nastaven v konfiguraci
    real_graph_fixed = nothing
    if "real" in cfg.networks.enabled
        real_path = isabspath(cfg.networks.csv.path) ? cfg.networks.csv.path : joinpath(PROJECT_ROOT, cfg.networks.csv.path)
        real_graph_fixed = loadCSVGraph2(real_path)
    end

    # inicializuj staticke vlastnosti (nad grafy) a dynamicke vlastnosti (nad modely)
    static_rows = DataFrame(
        network=String[],
        model=String[],
        epsilon=Float64[],
        mu=Float64[],
        run=Int[],
        avg_degree=Float64[],
        density=Float64[],
        apl=Float64[],
    )
    dynamic_rows = DataFrame(
        network=String[],
        model=String[],
        epsilon=Float64[],
        mu=Float64[],
        run=Int[],
        t=Int[],
        clusters=Int[],
        pad=Float64[],
        polarization=Float64[],
    )

    function _save_data_tables(network::String, model::String, epsilon::Float64, mu::Float64, run::Int, stat::NamedTuple, dyn_df::DataFrame)
        cfg.outputs.save_data_csv || return
        suffix = _param_suffix(model, epsilon, mu)
        base = joinpath(cfg.experiment.out_dir, "data", network, model, suffix)
        mkpath(base)

        stat_df = DataFrame(run=[run], epsilon=[epsilon], mu=[mu], avg_degree=[stat.avg_degree], density=[stat.density], apl=[stat.apl])
        CSV.write(joinpath(base, "run_$(lpad(string(run), 3, '0'))_static.csv"), stat_df)
        CSV.write(joinpath(base, "run_$(lpad(string(run), 3, '0'))_dynamic.csv"), dyn_df)
    end

    # pro kazdou sit spust kazdy algoritmus a kazdou kombinaci jeho parametru
    done = 0
    total = sum(length(_model_params(cfg, model)) for model in cfg.models.enabled) *
            length(cfg.networks.enabled) * cfg.experiment.runs
    simulation_started_at = time()
    for network in cfg.networks.enabled
        for model in cfg.models.enabled
            for params in _model_params(cfg, model)
                epsilon = Float64(params.epsilon)
                mu = Float64(params.mu)
                suffix = _param_suffix(model, epsilon, mu)

                for run in 1:cfg.experiment.runs
                done += 1
                # deterministické seedy (jine pro kazdy beh)
                seeds = derive_seeds(cfg, run, network, model)

                graph = nothing

                if network == "er"
                    cfg.networks.er === nothing && error("ER enabled but missing [networks.er]")
                    erc = cfg.networks.er
                    graph = generateER(erc.N, erc.p; rng=MersenneTwister(seeds.graph_seed))
                elseif network == "ba"
                    cfg.networks.ba === nothing && error("BA enabled but missing [networks.ba]")
                    bac = cfg.networks.ba
                    graph = generateBA(bac.N, bac.M; rng=MersenneTwister(seeds.graph_seed))
                elseif network == "real"
                    c = cfg.networks.csv
                    c === nothing && error("real enabled but missing [networks.real]")
                    graph = real_graph_fixed
                else
                    error("Unknown network: $network")
                end

                graph isa AbstractGraph || error("Nepodařilo se vytvořit graf pro síť $network")
                n = nv(graph)

                if cfg.outputs.save_edges_csv
                    base = joinpath(cfg.experiment.out_dir, "data", network, model, suffix)
                    mkpath(base)
                    CSV.write(joinpath(base, "run_$(lpad(string(run), 3, '0'))_edges.csv"), edge_dataframe(graph))
                end

                x0 = opinion_vector(cfg.opinions, graph; seed=seeds.opinions_seed)

                rng_sim = MersenneTwister(seeds.sim_seed)

                series = if model == "deffuant"
                    d = cfg.models.deffuant
                    d === nothing && error("Deffuant povolen ale chybi [models.deffuant]")
                    Deffuant.run_deffuant_series(graph, x0; iterations=d.iterations, ϵ=epsilon, μ=mu, rng=rng_sim, clamp_range=(d.clamp_low, d.clamp_high), snapshots=cfg.experiment.snapshots)
                elseif model == "hk"
                    h = cfg.models.hk
                    h === nothing && error("HK povolen ale chybi [models.hk]")
                    HegselmannKrause.run_hk_series(graph, x0; iterations=h.iterations, ϵ=epsilon, rng=rng_sim, include_self=h.include_self, clamp_range=(h.clamp_low, h.clamp_high), snapshots=cfg.experiment.snapshots)
                elseif model == "deffuant_hrdcr"
                    h = cfg.models.deffuant_hrdcr
                    h === nothing && error("Deffuant-HRDCR povolen ale chybi [models.deffuant_hrdcr]")
                    edges_df = edge_dataframe(graph)
                    DeffuantHRDCR.run_deffuant_hrdcr_series(
                        edges_df,
                        x0;
                        steps=h.iterations,
                        d=epsilon,
                        μ=mu,
                        s_max=h.s_max,
                        s_beta=h.s_beta,
                        eps_j=h.eps_j,
                        gamma_j=h.gamma_j,
                        rng=rng_sim,
                        n=n,
                        clamp_range=(h.clamp_low, h.clamp_high),
                        snapshots=cfg.experiment.snapshots,
                    )
                elseif model == "deffuant_anchored"
                    a = cfg.models.deffuant_anchored
                    DeffuantAnchored.run_deffuant_anchored_series(graph, x0;
                        iterations=a.iterations, epsilon=epsilon, mu=mu,
                        anchor_strength=a.anchor_strength, social_weight=a.social_weight,
                        degree_inertia=a.degree_inertia, confidence_shape=a.confidence_shape,
                        rng=rng_sim, clamp_range=(a.clamp_low, a.clamp_high),
                        snapshots=cfg.experiment.snapshots)
                elseif model == "deffuant_partisan"
                    p = cfg.models.deffuant_partisan
                    DeffuantPartisan.run_deffuant_partisan_series(graph, x0;
                        iterations=p.iterations, epsilon=epsilon, mu=mu,
                        cross_epsilon=p.cross_epsilon, repulsion=p.repulsion,
                        identity_strength=p.identity_strength, anchor_strength=p.anchor_strength,
                        degree_inertia=p.degree_inertia, left_pole=p.left_pole, right_pole=p.right_pole,
                        rng=rng_sim, clamp_range=(p.clamp_low,p.clamp_high), snapshots=cfg.experiment.snapshots)
                else
                    error("Unknown model: $model")
                end

                if cfg.outputs.save_final_opinion_snapshot
                    base = joinpath(cfg.experiment.out_dir, "data", network, model, suffix, "final_opinion_run_$(lpad(string(run), 3, '0'))")
                    mkpath(base)
                    last_idx = length(series) - 1
                    last_x = sort(series[end])
                    save_opinion_snapshot(base, last_idx, last_x; width=3)
                end

                if run == 1
                    ExperimentPlotGenerator.generate_first_run_opinion_plot(
                        series;
                        cfg=cfg,
                        model=model,
                        network=network,
                        suffix=suffix,
                    )
                end

                if cfg.outputs.save_opinion_snapshots
                    base = joinpath(cfg.experiment.out_dir, "data", network, model, suffix, "opinions_run_$(lpad(string(run), 3, '0'))")
                    mkpath(base)
                    for (idx, x) in enumerate(series)
                        save_opinion_snapshot(base, idx - 1, x; width=3)
                    end
                end

                stat = static_metrics(graph)
                dyn_df = dynamic_metrics_series(series)
                insertcols!(dyn_df, 1, :epsilon => fill(epsilon, nrow(dyn_df)), :mu => fill(mu, nrow(dyn_df)))

                push!(static_rows, (network, model, epsilon, mu, run, stat.avg_degree, stat.density, stat.apl))
                for r in eachrow(dyn_df)
                    push!(dynamic_rows, (network, model, epsilon, mu, run, Int(r.t), Int(r.clusters), Float64(r.pad), Float64(r.polarization)))
                end

                _save_data_tables(network, model, epsilon, mu, run, stat, dyn_df)

                elapsed = time() - simulation_started_at
                avg_per_iteration = elapsed / done
                remaining = avg_per_iteration * (total - done)
                percent = round(done / total * 100, digits=2)
                println(
                    "$(percent)% | $(done)/$(total) | elapsed $(_format_duration(elapsed)) | ETA $(_format_duration(remaining))"
                )
                end
            end
        end
    end

    dynamic_summary = Aggregator.aggregate_dynamic(dynamic_rows)

    if cfg.outputs.save_aggregated_csv
        agg_dir = joinpath(cfg.experiment.out_dir, "aggregated")

        static_summary = Aggregator.aggregate_static(static_rows)
        CSV.write(joinpath(agg_dir, "static_summary.csv"), static_summary)

        CSV.write(joinpath(agg_dir, "dynamic_summary.csv"), dynamic_summary)
    end

    ExperimentPlotGenerator.generate_standard_experiment_plots(dynamic_summary, cfg)

    if cfg.outputs.save_data_csv
        CSV.write(joinpath(cfg.experiment.out_dir, "data", "_all_static.csv"), static_rows)
        CSV.write(joinpath(cfg.experiment.out_dir, "data", "_all_dynamic.csv"), dynamic_rows)
    end

    ExperimentPlotGenerator.generate_thesis_experiment_plots(cfg)

    println("All done. Output: $(cfg.experiment.out_dir)")
end

main()
