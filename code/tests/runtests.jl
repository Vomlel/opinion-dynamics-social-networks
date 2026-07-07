module Tests

const CODE_DIR = normpath(joinpath(@__DIR__, ".."))

include(joinpath(CODE_DIR, "metrics", "dynamic_characteristics.jl"))
include(joinpath(CODE_DIR, "metrics", "static_characteristics.jl"))
include(joinpath(CODE_DIR, "core", "config.jl"))
include(joinpath(CODE_DIR, "reporting", "aggregator.jl"))
include(joinpath(CODE_DIR, "models", "deffuant.jl"))
include(joinpath(CODE_DIR, "models", "deffuant_hrdcr.jl"))
include(joinpath(CODE_DIR, "models", "deffuant_anchored.jl"))
include(joinpath(CODE_DIR, "models", "deffuant_partisan.jl"))
include(joinpath(CODE_DIR, "opinions", "generate_opinion.jl"))
include(joinpath(CODE_DIR, "models", "hk.jl"))
include(joinpath(CODE_DIR, "core", "helper.jl"))

using Test
using .DynamicCharacteristics
using .StaticCharacteristics
using .Config
using .Aggregator
using .GenerateOpinion
using CSV, DataFrames
using Graphs
using Random
using .Deffuant
using .DeffuantHRDCR
using .DeffuantAnchored
using .DeffuantPartisan
using .HegselmannKrause
using .Helper

@testset "avg" begin
    @test isnan(avg(Real[]))
    @test avg([1]) == 1
    @test avg([0,1]) == 0.5
    @test avg([0,1,1]) == 2/3
end

@testset "pairwise_average_distance" begin
    @test isnan(pairwise_average_distance(Real[]))
    @test isnan(pairwise_average_distance([1]))
    @test pairwise_average_distance([0,1]) == 1
    @test pairwise_average_distance([0,1,1]) == 2/3
end

@testset "pad_inner" begin
    @test isnan(pad_inner(Real[]))
    @test isnan(pad_inner([1]))
    @test isnan(pad_inner([1, 2]))
    @test isnan(pad_inner([1, 2, 3]))
    @test pad_inner([0,1,2,3]) == 1
    @test pad_inner([0,0,1,1,1]) == 0
end

@testset "pad_outer" begin
    @test isnan(pad_outer(Real[]))
    @test isnan(pad_outer([1]))
    @test pad_outer([1, 2]) == 1
    @test pad_outer([1, 2, 3]) == 2
    @test pad_outer([0,1,2,3]) == 8/4
    @test pad_outer([0,0,1,1,1]) == 4/4
end

@testset "polarization" begin
    @test isnan(polarization(Real[]))
    @test isnan(polarization([1]))
    @test isnan(polarization([1, 2]))
    @test isnan(polarization([1, 2, 3]))
    @test polarization([0,1,2,3]) == 2
    @test isinf(polarization([0,0,1,1,1]))
    @test polarization([0,1,1,1,1]) == 1
end

@testset "group_values" begin
    @test group_values(Real[], 0, 1, 3) == [0, 0, 0]
    @test group_values(Real[0, 0.1, 0.1, 0.5, 0.9, 1], 0, 1, 3) == [3, 1, 2]
end

@testset "count_of_group_clusters" begin
    @test count_of_group_clusters(Int[], 2, 4) == 0
    @test count_of_group_clusters([5, 1, 6, 1, 7], 2, 4) == 3
    #Mezi prvními dvema clustery není dostatečná malá skupins separationLimit=2, proto jsou jen 2
    @test count_of_group_clusters([5, 3, 6, 2, 7], 2, 4) == 2
    #Prostřední cluister je příliš malý sizeLimit=7
    @test count_of_group_clusters([7, 1, 3, 3, 1, 7], 2, 7) == 2
end

@testset "count_of_clusters_auto" begin
    @test count_of_clusters_auto(Real[]) == 0
    @test count_of_clusters_auto([1]) == 1
    @test count_of_clusters_auto([1, 2]) == 1
    @test count_of_clusters_auto([1, 1, 2, 2]) == 1
    @test count_of_clusters_auto([1, 1, 2, 2, 3, 3]) == 1
    c1 = [1, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13]
    c2 = [2, 2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.10, 2.11, 2.12, 2.13]
    c3 = [3, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.10, 3.11, 3.12, 3.13]
    @test count_of_clusters_auto([c1; c2; c3]) == 3
end

@testset "dyn_chars_1" begin
    graph = SimpleGraph(6)
    add_edge!(graph, 1, 2)
    add_edge!(graph, 2, 3)
    add_edge!(graph, 1, 4)
    add_edge!(graph, 2, 5)
    add_edge!(graph, 3, 4)
    add_edge!(graph, 4, 6)
    x0 = Vector{Float64}(undef, 6)
    x0[1] = 0.1
    x0[2] = 0.3
    x0[3] = 0.5
    x0[4] = 0.5
    x0[5] = 0.7
    x0[6] = 0.9
    series = Deffuant.run_deffuant_series(graph, x0; iterations=3, ϵ=0.3, μ=0.5, rng=Random.Xoshiro(123), clamp_range=(0, 1), snapshots=4)
    dyn_df = dynamic_metrics_series(series)
    @test dyn_df[!, :t] == [0, 1, 2, 3]
    @test dyn_df[!, :clusters] == [1, 1, 1, 1]
    @test dyn_df[!, :pad] ≈ [0.3466666666666666, 0.3466666666666666, 0.3466666666666666, 0.3466666666666666] atol=1e-10
    @test dyn_df[!, :polarization] ≈ [0.29520600208549064, 0.29520600208549064, 0.29520600208549064, 0.29520600208549064] atol=1e-10
end

@testset "dyn_chars_2" begin
    graph = SimpleGraph(6)
    add_edge!(graph, 1, 2)
    add_edge!(graph, 2, 3)
    add_edge!(graph, 1, 3)
    add_edge!(graph, 1, 4)
    add_edge!(graph, 4, 5)
    add_edge!(graph, 5, 6)
    add_edge!(graph, 4, 6)
    x0 = Vector{Float64}(undef, 6)
    x0[1] = 0.5
    x0[2] = 0.4
    x0[3] = 0.4
    x0[4] = 0.5
    x0[5] = 0.6
    x0[6] = 0.6
    series = Deffuant.run_deffuant_series(graph, x0; iterations=3, ϵ=0.1, μ=0.5, rng=Random.Xoshiro(123), clamp_range=(0, 1), snapshots=4)
    dyn_df = dynamic_metrics_series(series)
    @test dyn_df[!, :t] == [0, 1, 2, 3]
    @test dyn_df[!, :clusters] == [3, 3, 2, 2]
    @test dyn_df[!, :pad] ≈ [0.10666666666666665, 0.10666666666666665, 0.09999999999999999, 0.09666666666666664] atol=1e-10
    @test dyn_df[!, :polarization] ≈ [0.385749640681452, 0.385749640681452, 0.4020333696944648, 0.41340195033273014] atol=1e-10
end

@testset "dyn_chars_3" begin
    graph = SimpleGraph(6)
    add_edge!(graph, 1, 2)
    add_edge!(graph, 2, 3)
    add_edge!(graph, 1, 3)
    add_edge!(graph, 1, 4)
    add_edge!(graph, 4, 5)
    add_edge!(graph, 5, 6)
    add_edge!(graph, 4, 6)
    x0 = Vector{Float64}(undef, 6)
    x0[1] = 0.1
    x0[2] = 0.1
    x0[3] = 0.1
    x0[4] = 0.9
    x0[5] = 0.9
    x0[6] = 0.9
    series = Deffuant.run_deffuant_series(graph, x0; iterations=3, ϵ=0.1, μ=0.5, rng=Random.Xoshiro(123), clamp_range=(0, 1), snapshots=4)
    dyn_df = dynamic_metrics_series(series)
    @test dyn_df[!, :t] == [0, 1, 2, 3]
    @test dyn_df[!, :clusters] == [2, 2, 2, 2]
    @test dyn_df[!, :pad] ≈ [0.48, 0.48, 0.48, 0.48] atol=1e-10
    @test dyn_df[!, :polarization] ≈ [0.36546926110234157, 0.36546926110234157, 0.36546926110234157, 0.36546926110234157] atol=1e-10
end

@testset "config parameter grids" begin
    path = tempname() * ".toml"
    write(path, """
    [experiment]
    name = "TEST"
    out_dir = "out/test"
    runs = 1
    seed_base = 1
    snapshots = 3

    [outputs]
    save_aggregated_csv = true
    save_data_csv = false
    save_edges_csv = false
    save_opinion_snapshots = false
    save_final_opinion_snapshot = false
    plot_first_run_opinions = false
    plots = ["clusters"]

    [networks]
    enabled = ["er"]

    [networks.er]
    N = 10
    p = 0.2

    [opinions]
    model = "continuous"
    distribution = "uniform"
    low = 0.0
    high = 1.0
    mu = 0.5
    sigma = 0.1
    polarized_ratio = 0.5
    polarized_eps = 0.05

    [models]
    enabled = ["deffuant", "hk", "deffuant_hrdcr"]

    [models.deffuant]
    iterations = 5
    epsilon = [0.1, 0.2]
    mu = [0.1, 0.2]
    clamp_low = 0.0
    clamp_high = 1.0

    [models.hk]
    iterations = 2
    epsilon = [0.1, 0.2]
    include_self = true
    clamp_low = 0.0
    clamp_high = 1.0

    [models.deffuant_hrdcr]
    iterations = 5
    epsilon = [0.1, 0.2]
    mu = [0.1, 0.2]
    s_max = 0.7
    s_beta = 1.0
    eps_j = 0.2
    gamma_j = 0.5
    clamp_low = 0.0
    clamp_high = 1.0

    [dynamics]
    cluster_delta = 0.01
    polarization_threshold = 0.5
    """)

    cfg = Config.load_config(path)
    @test cfg.models.deffuant.epsilon == [0.1, 0.2]
    @test cfg.models.deffuant.mu == [0.1, 0.2]
    @test cfg.models.hk.epsilon == [0.1, 0.2]
    @test cfg.models.deffuant_hrdcr.epsilon == [0.1, 0.2]
    @test cfg.models.deffuant_hrdcr.mu == [0.1, 0.2]
end

@testset "aggregate_dynamic keeps parameter combinations separate" begin
    rows = DataFrame(
        network=["er", "er", "er", "er"],
        model=["deffuant", "deffuant", "deffuant", "deffuant"],
        epsilon=[0.1, 0.1, 0.2, 0.2],
        mu=[0.1, 0.1, 0.1, 0.1],
        run=[1, 2, 1, 2],
        t=[0, 0, 0, 0],
        clusters=[2, 4, 10, 12],
        pad=[0.1, 0.3, 0.5, 0.7],
        polarization=[0.2, 0.4, 0.6, 0.8],
    )

    summary = Aggregator.aggregate_dynamic(rows)
    @test nrow(summary) == 2
    @test sort(summary.clusters_mean) == [3.0, 11.0]
    @test Set(summary.epsilon) == Set([0.1, 0.2])
end

@testset "deffuant_hrdcr series" begin
    edges = DataFrame(source=[1], target=[2])
    x0 = [0.0, 1.0]

    unchanged = DeffuantHRDCR.run_deffuant_hrdcr_series(
        edges,
        x0;
        steps=3,
        d=0.0,
        μ=0.5,
        eps_j=1.0,
        gamma_j=0.0,
        s_max=0.0,
        s_beta=1.0,
        rng=Random.Xoshiro(1),
        n=2,
        snapshots=4,
    )
    @test unchanged[end] == x0

    converged = DeffuantHRDCR.run_deffuant_hrdcr_series(
        edges,
        x0;
        steps=1,
        d=2.0,
        μ=0.5,
        eps_j=1.0,
        gamma_j=0.0,
        s_max=0.0,
        s_beta=1.0,
        rng=Random.Xoshiro(1),
        n=2,
        snapshots=2,
    )
    @test converged[end] ≈ [0.5, 0.5]
end

@testset "deffuant anchored" begin
    g = path_graph(2)
    x0 = [0.0, 1.0]
    series = DeffuantAnchored.run_deffuant_anchored_series(g, x0;
        iterations=2, epsilon=100.0, mu=0.5, anchor_strength=0.0,
        social_weight=0.0, degree_inertia=0.0, rng=Random.Xoshiro(7), snapshots=3)
    @test length(series) == 3
    @test series[2] ≈ [0.5, 0.5]
    @test all(0 .<= series[end] .<= 1)
end

@testset "deffuant partisan separates camps" begin
    g = complete_graph(4); x0 = [0.45,0.49,0.51,0.55]
    out = DeffuantPartisan.run_deffuant_partisan(g,x0; iterations=2000, epsilon=0.2,
        mu=0.2, identity_strength=0.02, anchor_strength=0.01, rng=Random.Xoshiro(3))
    @test sum(out[3:4])/2 - sum(out[1:2])/2 > 0.2
    @test all(0 .<= out .<= 1)
end

end # module
