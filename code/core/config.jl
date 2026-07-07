

module Config

using TOML

export ExperimentConfig,
       NetworksConfig,
       OpinionsConfig,
       DeffuantConfig,
       DeffuantHRDCRConfig,
       DeffuantAnchoredConfig,
       DeffuantPartisanConfig,
       HKConfig,
       DynamicsConfig,
       OutputsConfig,
       FullConfig,
       load_config,
       derive_seeds

# -------------------------------
# Typed config structs
# -------------------------------

struct ExperimentConfig
    name::String
    out_dir::String
    runs::Int
    seed_base::Int
    snapshots::Int
end

struct OutputsConfig
    save_aggregated_csv::Bool
    save_data_csv::Bool
    save_edges_csv::Bool
    save_opinion_snapshots::Bool
    save_final_opinion_snapshot::Bool
    plot_first_run_opinions::Bool
    plots::Vector{String}
end

struct ERConfig
    N::Int
    p::Float64
end

struct BAConfig
    N::Int
    M::Int
end

struct CSVConfig
    path::String
    zero_based::Bool
    undirected::Bool
    n::Union{Nothing,Int}
end

struct NetworksConfig
    enabled::Vector{String}
    er::Union{Nothing,ERConfig}
    ba::Union{Nothing,BAConfig}
    csv::Union{Nothing,CSVConfig}
end

struct OpinionsConfig
    model::String            # "continuous" | "binary"
    distribution::String     # "uniform" | "normal" | "polarized"
    low::Float64
    high::Float64
    mu::Float64
    sigma::Float64
    polarized_ratio::Float64
    polarized_eps::Float64
end

struct DeffuantConfig
    iterations::Int
    epsilon::Vector{Float64}
    mu::Vector{Float64}
    clamp_low::Float64
    clamp_high::Float64
end

struct HKConfig
    iterations::Int
    epsilon::Vector{Float64}
    include_self::Bool
    clamp_low::Float64
    clamp_high::Float64
end

struct DeffuantHRDCRConfig
    iterations::Int
    epsilon::Vector{Float64}
    mu::Vector{Float64}
    s_max::Float64
    s_beta::Float64
    eps_j::Float64
    gamma_j::Float64
    clamp_low::Float64
    clamp_high::Float64
end

struct DeffuantAnchoredConfig
    iterations::Int
    epsilon::Vector{Float64}
    mu::Vector{Float64}
    anchor_strength::Float64
    social_weight::Float64
    degree_inertia::Float64
    confidence_shape::Float64
    clamp_low::Float64
    clamp_high::Float64
end

struct DeffuantPartisanConfig
    iterations::Int; epsilon::Vector{Float64}; mu::Vector{Float64}
    cross_epsilon::Float64; repulsion::Float64; identity_strength::Float64
    anchor_strength::Float64; degree_inertia::Float64
    left_pole::Float64; right_pole::Float64; clamp_low::Float64; clamp_high::Float64
end

struct ModelsConfig
    enabled::Vector{String}
    deffuant::Union{Nothing,DeffuantConfig}
    hk::Union{Nothing,HKConfig}
    deffuant_hrdcr::Union{Nothing,DeffuantHRDCRConfig}
    deffuant_anchored::Union{Nothing,DeffuantAnchoredConfig}
    deffuant_partisan::Union{Nothing,DeffuantPartisanConfig}
end

struct DynamicsConfig
    cluster_delta::Float64
    polarization_threshold::Float64
end

struct FullConfig
    experiment::ExperimentConfig
    outputs::OutputsConfig
    networks::NetworksConfig
    opinions::OpinionsConfig
    models::ModelsConfig
    dynamics::DynamicsConfig
end

# -------------------------------
# Internal helpers
# -------------------------------

_get(tbl::AbstractDict, key::AbstractString) = haskey(tbl, key) ? tbl[key] : throw(ArgumentError("Chybí klíč '$key'"))

function _get_opt(tbl::AbstractDict, key::AbstractString)
    return haskey(tbl, key) ? tbl[key] : nothing
end

function _as_int(x, name::AbstractString)
    x isa Integer || throw(ArgumentError("$name musí být Integer"))
    return Int(x)
end

function _as_float(x, name::AbstractString)
    x isa Real || throw(ArgumentError("$name musí být číslo"))
    return Float64(x)
end

function _as_bool(x, name::AbstractString)
    x isa Bool || throw(ArgumentError("$name musí být Bool"))
    return x
end

function _as_string(x, name::AbstractString)
    x isa AbstractString || throw(ArgumentError("$name musí být String"))
    return String(x)
end

function _as_vec_string(x, name::AbstractString)
    x isa AbstractVector || throw(ArgumentError("$name musí být seznam (array)"))
    out = String[]
    for v in x
        push!(out, _as_string(v, name))
    end
    return out
end

function _as_vec_float(x, name::AbstractString)
    if x isa AbstractVector
        out = Float64[]
        for v in x
            push!(out, _as_float(v, name))
        end
        isempty(out) && throw(ArgumentError("$name nesmí být prázdný seznam"))
        return out
    end
    return [_as_float(x, name)]
end

function _validate_in(value::String, allowed::Vector{String}, name::AbstractString)
    value in allowed || throw(ArgumentError("$name='$value' není povoleno. Povoleno: $(join(allowed, ", "))"))
    return value
end

function _validate_positive_int(v::Int, name::AbstractString; allow_zero::Bool=false)
    if allow_zero
        v >= 0 || throw(ArgumentError("$name musí být >= 0"))
    else
        v >= 1 || throw(ArgumentError("$name musí být >= 1"))
    end
    return v
end

function _validate_range(low::Float64, high::Float64, name::AbstractString)
    low < high || throw(ArgumentError("$name: low musí být < high"))
end

# -------------------------------
# Public API
# -------------------------------

"""
    load_config(path::AbstractString) -> FullConfig

Načte `config.toml` a vrátí typovaný `FullConfig`.

Vyhazuje `ArgumentError` při chybějících klíčích nebo neplatných hodnotách.
"""
function load_config(path::AbstractString)::FullConfig
    raw = TOML.parsefile(path)

    # -------- experiment --------
    exp = _get(raw, "experiment")
    exp isa AbstractDict || throw(ArgumentError("[experiment] musí být tabulka"))

    name = _as_string(_get(exp, "name"), "experiment.name")
    out_dir = _as_string(_get(exp, "out_dir"), "experiment.out_dir")
    runs = _validate_positive_int(_as_int(_get(exp, "runs"), "experiment.runs"), "experiment.runs")
    seed_base = _as_int(_get(exp, "seed_base"), "experiment.seed_base")
    snapshots = _validate_positive_int(_as_int(_get(exp, "snapshots"), "experiment.snapshots"), "experiment.snapshots")

    experiment = ExperimentConfig(name, out_dir, runs, seed_base, snapshots)

    # -------- outputs --------
    out = _get(raw, "outputs")
    out isa AbstractDict || throw(ArgumentError("[outputs] musí být tabulka"))

    save_aggregated_csv = _as_bool(_get(out, "save_aggregated_csv"), "outputs.save_aggregated_csv")
    save_data_csv = _as_bool(_get(out, "save_data_csv"), "outputs.save_data_csv")
    save_edges_csv = _as_bool(_get(out, "save_edges_csv"), "outputs.save_edges_csv")
    save_opinion_snapshots = _as_bool(_get(out, "save_opinion_snapshots"), "outputs.save_opinion_snapshots")
    save_final_opinion_snapshot = _as_bool(_get(out, "save_final_opinion_snapshot"), "outputs.save_final_opinion_snapshot")
    plot_first_run_opinions = _as_bool(_get(out, "plot_first_run_opinions"), "outputs.plot_first_run_opinions")
    plots = _as_vec_string(_get(out, "plots"), "outputs.plots")

    outputs = OutputsConfig(save_aggregated_csv, save_data_csv, save_edges_csv, save_opinion_snapshots, save_final_opinion_snapshot, plot_first_run_opinions, plots)

    # -------- networks --------
    nets = _get(raw, "networks")
    nets isa AbstractDict || throw(ArgumentError("[networks] musí být tabulka"))

    net_enabled = _as_vec_string(_get(nets, "enabled"), "networks.enabled")
    net_enabled = [n == "csv" ? "real" : n for n in net_enabled]

    er_cfg = nothing
    if "er" in net_enabled
        er = _get(nets, "er")
        er isa AbstractDict || throw(ArgumentError("[networks.er] musí být tabulka"))
        N = _validate_positive_int(_as_int(_get(er, "N"), "networks.er.N"), "networks.er.N")
        p = _as_float(_get(er, "p"), "networks.er.p")
        (0.0 <= p <= 1.0) || throw(ArgumentError("networks.er.p musí být v [0,1]"))
        er_cfg = ERConfig(N, p)
    end

    ba_cfg = nothing
    if "ba" in net_enabled
        ba = _get(nets, "ba")
        ba isa AbstractDict || throw(ArgumentError("[networks.ba] musí být tabulka"))
        N = _validate_positive_int(_as_int(_get(ba, "N"), "networks.ba.N"), "networks.ba.N")
        M = _validate_positive_int(_as_int(_get(ba, "M"), "networks.ba.M"), "networks.ba.M")
        M < N || throw(ArgumentError("networks.ba.M musí být < networks.ba.N"))
        ba_cfg = BAConfig(N, M)
    end

    csv_cfg = nothing
    if "real" in net_enabled
        real = haskey(nets, "real") ? _get(nets, "real") : _get(nets, "csv")
        real isa AbstractDict || throw(ArgumentError("[networks.real] musí být tabulka"))

        path_s = _as_string(_get(real, "path"), "networks.real.path")
        zero_based = _as_bool(_get(real, "zero_based"), "networks.real.zero_based")
        undirected = _as_bool(_get(real, "undirected"), "networks.real.undirected")

        n_opt = nothing
        if haskey(real, "n")
            n_raw = real["n"]
            n_opt = _validate_positive_int(_as_int(n_raw, "networks.real.n"), "networks.real.n")
        end

        csv_cfg = CSVConfig(path_s, zero_based, undirected, n_opt)
    end

    networks = NetworksConfig(net_enabled, er_cfg, ba_cfg, csv_cfg)

    # -------- opinions --------
    op = _get(raw, "opinions")
    op isa AbstractDict || throw(ArgumentError("[opinions] musí být tabulka"))

    model = _validate_in(_as_string(_get(op, "model"), "opinions.model"), ["continuous", "binary"], "opinions.model")
    distribution = _validate_in(_as_string(_get(op, "distribution"), "opinions.distribution"), ["uniform", "normal", "polarized"], "opinions.distribution")

    low = _as_float(_get(op, "low"), "opinions.low")
    high = _as_float(_get(op, "high"), "opinions.high")
    _validate_range(low, high, "opinions")

    mu = _as_float(_get(op, "mu"), "opinions.mu")
    sigma = _as_float(_get(op, "sigma"), "opinions.sigma")
    sigma >= 0.0 || throw(ArgumentError("opinions.sigma musí být >= 0"))

    polarized_ratio = _as_float(_get(op, "polarized_ratio"), "opinions.polarized_ratio")
    0.0 <= polarized_ratio <= 1.0 || throw(ArgumentError("opinions.polarized_ratio musí být v [0,1]"))

    polarized_eps = _as_float(_get(op, "polarized_eps"), "opinions.polarized_eps")
    polarized_eps >= 0.0 || throw(ArgumentError("opinions.polarized_eps musí být >= 0"))

    opinions = OpinionsConfig(model, distribution, low, high, mu, sigma, polarized_ratio, polarized_eps)

    # -------- models --------
    mods = _get(raw, "models")
    mods isa AbstractDict || throw(ArgumentError("[models] musí být tabulka"))

    model_enabled = _as_vec_string(_get(mods, "enabled"), "models.enabled")

    deff_cfg = nothing
    if "deffuant" in model_enabled
        d = _get(mods, "deffuant")
        d isa AbstractDict || throw(ArgumentError("[models.deffuant] musí být tabulka"))

        iterations = _validate_positive_int(_as_int(_get(d, "iterations"), "models.deffuant.iterations"), "models.deffuant.iterations")
        eps = _as_vec_float(_get(d, "epsilon"), "models.deffuant.epsilon")
        all(v -> v >= 0.0, eps) || throw(ArgumentError("models.deffuant.epsilon musí obsahovat hodnoty >= 0"))
        mu_d = _as_vec_float(_get(d, "mu"), "models.deffuant.mu")
        all(v -> 0.0 < v <= 0.5, mu_d) || throw(ArgumentError("models.deffuant.mu musí obsahovat hodnoty v (0,0.5]"))

        clamp_low = _as_float(_get(d, "clamp_low"), "models.deffuant.clamp_low")
        clamp_high = _as_float(_get(d, "clamp_high"), "models.deffuant.clamp_high")
        _validate_range(clamp_low, clamp_high, "models.deffuant.clamp_range")

        deff_cfg = DeffuantConfig(iterations, eps, mu_d, clamp_low, clamp_high)
    end

    hk_cfg = nothing
    if "hk" in model_enabled
        h = _get(mods, "hk")
        h isa AbstractDict || throw(ArgumentError("[models.hk] musí být tabulka"))

        iterations = _validate_positive_int(_as_int(_get(h, "iterations"), "models.hk.iterations"), "models.hk.iterations")
        eps = _as_vec_float(_get(h, "epsilon"), "models.hk.epsilon")
        all(v -> v >= 0.0, eps) || throw(ArgumentError("models.hk.epsilon musí obsahovat hodnoty >= 0"))
        include_self = _as_bool(_get(h, "include_self"), "models.hk.include_self")

        clamp_low = _as_float(_get(h, "clamp_low"), "models.hk.clamp_low")
        clamp_high = _as_float(_get(h, "clamp_high"), "models.hk.clamp_high")
        _validate_range(clamp_low, clamp_high, "models.hk.clamp_range")

        hk_cfg = HKConfig(iterations, eps, include_self, clamp_low, clamp_high)
    end

    hrdcr_cfg = nothing
    if "deffuant_hrdcr" in model_enabled
        h = _get(mods, "deffuant_hrdcr")
        h isa AbstractDict || throw(ArgumentError("[models.deffuant_hrdcr] musí být tabulka"))

        iterations = _validate_positive_int(_as_int(_get(h, "iterations"), "models.deffuant_hrdcr.iterations"), "models.deffuant_hrdcr.iterations")
        eps = _as_vec_float(_get(h, "epsilon"), "models.deffuant_hrdcr.epsilon")
        all(v -> v >= 0.0, eps) || throw(ArgumentError("models.deffuant_hrdcr.epsilon musí obsahovat hodnoty >= 0"))

        mu_h = _as_vec_float(_get(h, "mu"), "models.deffuant_hrdcr.mu")
        all(v -> 0.0 < v <= 0.5, mu_h) || throw(ArgumentError("models.deffuant_hrdcr.mu musí obsahovat hodnoty v (0,0.5]"))

        s_max = _as_float(_get(h, "s_max"), "models.deffuant_hrdcr.s_max")
        (0.0 <= s_max <= 1.0) || throw(ArgumentError("models.deffuant_hrdcr.s_max musí být v [0,1]"))
        s_beta = _as_float(_get(h, "s_beta"), "models.deffuant_hrdcr.s_beta")
        s_beta >= 0.0 || throw(ArgumentError("models.deffuant_hrdcr.s_beta musí být >= 0"))

        eps_j = _as_float(_get(h, "eps_j"), "models.deffuant_hrdcr.eps_j")
        (0.0 < eps_j <= 1.0) || throw(ArgumentError("models.deffuant_hrdcr.eps_j musí být v (0,1]"))
        gamma_j = _as_float(_get(h, "gamma_j"), "models.deffuant_hrdcr.gamma_j")
        gamma_j >= 0.0 || throw(ArgumentError("models.deffuant_hrdcr.gamma_j musí být >= 0"))

        clamp_low = _as_float(_get(h, "clamp_low"), "models.deffuant_hrdcr.clamp_low")
        clamp_high = _as_float(_get(h, "clamp_high"), "models.deffuant_hrdcr.clamp_high")
        _validate_range(clamp_low, clamp_high, "models.deffuant_hrdcr.clamp_range")

        hrdcr_cfg = DeffuantHRDCRConfig(iterations, eps, mu_h, s_max, s_beta, eps_j, gamma_j, clamp_low, clamp_high)
    end

    anchored_cfg = nothing
    if "deffuant_anchored" in model_enabled
        a = _get(mods, "deffuant_anchored")
        a isa AbstractDict || throw(ArgumentError("[models.deffuant_anchored] musi byt tabulka"))
        iterations = _validate_positive_int(_as_int(_get(a, "iterations"), "models.deffuant_anchored.iterations"), "models.deffuant_anchored.iterations")
        eps = _as_vec_float(_get(a, "epsilon"), "models.deffuant_anchored.epsilon")
        all(>(0.0), eps) || throw(ArgumentError("models.deffuant_anchored.epsilon musi byt > 0"))
        mus = _as_vec_float(_get(a, "mu"), "models.deffuant_anchored.mu")
        all(v -> 0 < v <= 0.5, mus) || throw(ArgumentError("models.deffuant_anchored.mu musi byt v (0,0.5]"))
        anchor = _as_float(_get(a, "anchor_strength"), "models.deffuant_anchored.anchor_strength")
        social = _as_float(_get(a, "social_weight"), "models.deffuant_anchored.social_weight")
        inertia = _as_float(_get(a, "degree_inertia"), "models.deffuant_anchored.degree_inertia")
        shape = _as_float(_get(a, "confidence_shape"), "models.deffuant_anchored.confidence_shape")
        0 <= anchor <= 1 || throw(ArgumentError("anchor_strength musi byt v [0,1]"))
        0 <= social <= 1 || throw(ArgumentError("social_weight musi byt v [0,1]"))
        inertia >= 0 || throw(ArgumentError("degree_inertia musi byt >= 0"))
        shape > 0 || throw(ArgumentError("confidence_shape musi byt > 0"))
        low = _as_float(_get(a, "clamp_low"), "models.deffuant_anchored.clamp_low")
        high = _as_float(_get(a, "clamp_high"), "models.deffuant_anchored.clamp_high")
        _validate_range(low, high, "models.deffuant_anchored.clamp_range")
        anchored_cfg = DeffuantAnchoredConfig(iterations, eps, mus, anchor, social, inertia, shape, low, high)
    end

    partisan_cfg = nothing
    if "deffuant_partisan" in model_enabled
        p = _get(mods, "deffuant_partisan")
        iterations = _validate_positive_int(_as_int(_get(p,"iterations"), "iterations"), "iterations")
        eps = _as_vec_float(_get(p,"epsilon"), "epsilon"); mus = _as_vec_float(_get(p,"mu"), "mu")
        ce = _as_float(_get(p,"cross_epsilon"), "cross_epsilon")
        rep = _as_float(_get(p,"repulsion"), "repulsion")
        ids = _as_float(_get(p,"identity_strength"), "identity_strength")
        anc = _as_float(_get(p,"anchor_strength"), "anchor_strength")
        di = _as_float(_get(p,"degree_inertia"), "degree_inertia")
        lp = _as_float(_get(p,"left_pole"), "left_pole"); rp = _as_float(_get(p,"right_pole"), "right_pole")
        low = _as_float(_get(p,"clamp_low"), "clamp_low"); high = _as_float(_get(p,"clamp_high"), "clamp_high")
        all(>(0), eps) || throw(ArgumentError("epsilon musi byt > 0")); all(v->0<v<=0.5,mus) || throw(ArgumentError("mu mimo rozsah"))
        ce > 0 || throw(ArgumentError("cross_epsilon musi byt > 0")); minimum((rep,ids,anc,di)) >= 0 || throw(ArgumentError("sily musi byt >= 0"))
        0 <= lp < rp <= 1 || throw(ArgumentError("neplatne stranicke poly")); _validate_range(low,high,"clamp")
        partisan_cfg = DeffuantPartisanConfig(iterations,eps,mus,ce,rep,ids,anc,di,lp,rp,low,high)
    end

    models = ModelsConfig(model_enabled, deff_cfg, hk_cfg, hrdcr_cfg, anchored_cfg, partisan_cfg)

    # -------- dynamics --------
    dyn = _get(raw, "dynamics")
    dyn isa AbstractDict || throw(ArgumentError("[dynamics] musí být tabulka"))

    cluster_delta = _as_float(_get(dyn, "cluster_delta"), "dynamics.cluster_delta")
    cluster_delta >= 0.0 || throw(ArgumentError("dynamics.cluster_delta musí být >= 0"))

    pol_thr = _as_float(_get(dyn, "polarization_threshold"), "dynamics.polarization_threshold")

    dynamics = DynamicsConfig(cluster_delta, pol_thr)

    return FullConfig(experiment, outputs, networks, opinions, models, dynamics)
end

"""
    derive_seeds(cfg::FullConfig, run_id::Integer, network::AbstractString, model::AbstractString)

Vygeneruje deterministické seedy pro:
- generování grafu (`graph_seed`)
- generování počátečních názorů (`opinions_seed`)
- samotnou simulaci (`sim_seed`)

Použití: `seeds = derive_seeds(cfg, run, "er", "deffuant")`.
"""
function derive_seeds(cfg::FullConfig, run_id::Integer, network::AbstractString, model::AbstractString)
    base = cfg.experiment.seed_base

    # jednoduché stabilní mapování string -> int
    function s2i(s::AbstractString)
        h = UInt(0)
        for c in codeunits(s)
            h = h * UInt(131) + UInt(c)
        end
        return Int(mod(h, UInt(1_000_000_000)))
    end

    r = Int(run_id)
    r >= 1 || throw(ArgumentError("run_id musí být >= 1"))

    netv = s2i(network)
    modv = s2i(model)

    graph_seed    = base + 10_000_000 * r + 1_000 * netv + 1
    opinions_seed = base + 10_000_000 * r + 1_000 * netv + 2
    sim_seed      = base + 10_000_000 * r + 1_000 * netv + 10 + modv % 900

    return (graph_seed=graph_seed, opinions_seed=opinions_seed, sim_seed=sim_seed)
end

end # module
