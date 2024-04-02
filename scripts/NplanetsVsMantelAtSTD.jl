## NplanetsVsMantelAtSTD.jl
## Recreates Fig. A3 in Smith & Sinapayen, 2024 
## https://arxiv.org/abs/2403.14195
## Written by Harrison B. Smith, 2024
## Instructions:
## > julia 
## > include("NplanetsVsMantelAtSTD.jl")
## > main()

using Pkg
println(@__DIR__)
Pkg.activate(@__DIR__)

using Plots
using Random
# Pkg.develop(path=".") #../TerraformingAgents" ## Don't need this if I start the Julia REPL from the package itself
using TerraformingAgents
using Distances
using LinearAlgebra: diag, issymmetric, tril!
using Statistics

"""
    PermuteDistanceMatrix(D; rng::AbstractRNG = Random.default_rng())

Randomly permute the rows and columns of a distance matrix `D`.
"""
function PermuteDistanceMatrix(D; rng::AbstractRNG = Random.default_rng())
    order = shuffle(rng, collect(1:size(D)[1]))
    return D[order,:][:,order]
end

"""
    triangleFlatten(D)

Flatten the lower triangle of a square matrix `D` into a 1D array.
"""
triangleFlatten(D) = D[tril!(trues(size(D)), -1)]


"""
    MantelTestPermutations(
        x, 
        y; 
        rng::AbstractRNG=Random.default_rng(), 
        dist_metric=Euclidean(), 
        method=:pearson, 
        permutations=999, 
        alternative=:twosided)

Perform a Mantel test between two distance matrices `x` and `y` with `permutations` 
    permutations. Based off of: 
    https://github.com/biocore/scikit-bio/blob/0.1.3/skbio/math/stats/distance/_mantel.py#L23

Returns the original statistic (the mantel coefficient), the p-value, and the 
    permuted statistics (the mantel coefficients of the permutations).


"""
function MantelTestPermutations(
    x, 
    y; 
    rng::AbstractRNG=Random.default_rng(), 
    dist_metric=Euclidean(), 
    method=:pearson, 
    permutations=999, 
    alternative=:twosided
    )

    ## Check validity of arguments
    method == :pearson ? corr_func = cor : throw(ArgumentError("Not yet implemented"))
    permutations < 0 && throw(ArgumentError("Number of permutations must be greater than or equal to zero."))
    alternative ∉ [:twosided, :greater, :less] && throw(ArgumentError("Invalid alternative hypothesis $alternative."))

    ## Check x and y to verify they're distance matrices
    size(x) != size(y) && throw(ArgumentError("Distance matrices must have the same shape."))
    size(x)[1] < 3 && throw(ArgumentError("Distance matrices must be at least 3x3 in size."))
    sum(abs.(diag(x))) + sum(abs.(diag(y))) != 0 && throw(ArgumentError("Distance matrices must be hollow."))
    ~issymmetric(x) | ~issymmetric(y) && throw(ArgumentError("Distance matrices must be symmetric."))

    ## This part just needs to get a flattened version of the diagonal of a hollow, square, symmetric matrix
    x_flat = triangleFlatten(x)
    y_flat = triangleFlatten(y)

    orig_stat = corr_func(x_flat, y_flat)

    ## Permutation tests
    if (permutations == 0) | isnan(orig_stat)
        p_value = NaN
        permuted_stats = NaN
    else
        perm_gen = (cor(triangleFlatten(PermuteDistanceMatrix(x, rng=rng)), y_flat) for _ in 1:permutations)
        permuted_stats = collect(Iterators.flatten(perm_gen))

        if alternative == :twosided
            count_better = sum(abs.(permuted_stats) .>= abs(orig_stat))
        elseif alternative == :greater
            count_better = sum(permuted_stats .>= orig_stat)
        else
            count_better = sum(permuted_stats .<= orig_stat)
        end

        p_value = (count_better + 1) / (permutations + 1)

    end

    return orig_stat, p_value, permuted_stats

end

"""
    MantelTestCompositionsPositions(
        comp, 
        pos; 
        dist_metric=Euclidean(), 
        rng::AbstractRNG=Random.default_rng(), 
        method=:pearson, 
        permutations=999, 
        alternative=:twosided)

Perform a Mantel test specifically between the composition matrix and position matrix as 
    they are formatting in the TerraformingAgents initialization.

Returns MantelTestPermutations ouputs.
"""
function MantelTestCompositionsPositions(
    comp, 
    pos; 
    dist_metric=Euclidean(), 
    rng::AbstractRNG=Random.default_rng(), 
    method=:pearson, 
    permutations=999, 
    alternative=:twosided
    )

    x = pairwise(dist_metric, comp, dims=2)
    y = pairwise(dist_metric, hcat(collect.(pos)...), dims=2)

    MantelTestPermutations(x, y; rng=rng, dist_metric=dist_metric, method=method, permutations=permutations, alternative=alternative)

end

"""
    MantelTestTerraformingAgents(
        rng::AbstractRNG=Random.default_rng(), 
        extent=(100, 100, 100), 
        nplanets=1000, 
        maxcomp=1, 
        compsize=10, 
        permutations=999, 
        alternative=:twosided)

Initialize from the TerraformingAgents simulation the composition and position matrices, 
    and perform a Mantel test between them.

Returns MantelTestPermutations ouputs.
"""
function MantelTestTerraformingAgents(;
    rng = MersenneTwister(3141),
    extent = (100, 100, 100),
    nplanets = 1000,
    maxcomp = 1,
    compsize = 10,
    permutations = 999,
    alternative=:twosided
    )
    
    pos = TerraformingAgents.random_positions(rng, extent, nplanets)
    comp = TerraformingAgents.random_compositions(rng, maxcomp, compsize, nplanets)
    orig_stat, p_value, permuted_stats = MantelTestCompositionsPositions(comp, pos; dist_metric=Euclidean(), rng=rng, method=:pearson, permutations=permutations, alternative=alternative)
    
end

"""
    NplanetsVsMantelAtSTD(nplanet_range, std_vals; kwargs...)

Calculate what Mantel coefficient corresponds to the sigma anomalies provided in 
    `std_vals`, for a range of nplanets present and observed in the simulation 
    `nplanet_range`.

Returns a 2D matrix of Mantel coefficients for each sigma value in `std_vals` and 
    each nplanet in `nplanet_range`.
    
"""
function NplanetsVsMantelAtSTD(nplanet_range, std_vals; kwargs...)
#permutations = 500 should be good
    mantels = []
    for nplanets in nplanet_range
        @show nplanets
        o, p, ps = MantelTestTerraformingAgents(nplanets=nplanets; kwargs...)
        std1 = Statistics.stdm(ps, mean(ps))
        std2345 = std1 .* std_vals
        push!(mantels, std2345)
    end
    
    hcat(mantels...)

end

"""
    PlotNplanetsVsMantelAtSTD(mantels, nplanet_range, std_vals; kwargs...)

Returns a plot of the Mantel coefficients for each sigma value in `std_vals` and each nplanet in
    `nplanet_range`, calcuated from `NplanetsVsMantelAtSTD`.
"""
function PlotNplanetsVsMantelAtSTD(mantels, nplanet_range, std_vals; kwargs...)
    labels = collect(string(i)*"σ" for i in std_vals)
    
    ## Plot
    P = plot(
        xlabel = "Number of planets",
        ylabel = "log(Mantel coefficient)",
    )
    for i in 1:length(std_vals)
        plot!(
            P, 
            nplanet_range, 
            mantels[i,:], 
            label=labels[i], 
            palette=:darkrainbow, 
            marker=(:circle,5), 
            yaxis=:log10,
            xticks=[10,100,250,500,1000,2000],
            # xminorticks=nplanet_range,
            minorgrid=true)
    end
    P

end

## Creates plot seen in the appendix of the paper
function main()
    std_vals = [2.5, 5]
    nplanet_range = [5, 10, 25, 75, 100, 250, 500, 1000, 2000]
    mantels = NplanetsVsMantelAtSTD(nplanet_range, std_vals; permutations=1000)
    P = PlotNplanetsVsMantelAtSTD(mantels, nplanet_range, std_vals) 
    gui(P)
    # savefig(PlotNplanetsVsMantelAtSTD(mantels, nplanet_range, std_vals), "../../nplanets_vs_mantel_at_sigma_test.pdf")
end
