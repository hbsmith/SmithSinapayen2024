## run_simulation_parallel.jl
## Collect data from a run with 1000 planets for 1000 timesteps
## Uses parallel processing
## Generates data found here: [FIGSHARE]
## Written by Harrison B. Smith, 2024
## Instructions (will output files to TerraformingAgents_Projects/output/):
## > julia
## > include("run_simulation_parallel.jl")

## Set number of cores to use
using Distributed
desired_n_workers = 7

## Below code makes sure you don't accidentally spawn too many procs when calling code 
##  using include()
worker_diff = desired_n_workers - nworkers()
if worker_diff > 0
    addprocs(worker_diff)
elseif worker_diff < 0
    rmprocs(filter!(x -> x > desired_n_workers + 1, workers()))
end
@show workers()

## Must activate independently of imports
@everywhere begin
    using Pkg
    println(@__DIR__)
    Pkg.activate(@__DIR__)
end
@everywhere begin
    using Agents, Plots, Random
    using TerraformingAgents
    using CSVFiles
    using DataFrames
    using Dates
    using DrWatson: @dict
    using JSON
    using YAML
    using OrderedCollections

    n = 1000 # number of steps to run simulation for
    parameters = Dict(
        :nplanets => 1000,
        :extent => (100,100,100),
        :nool => [1,10],
        :dt => 10,
        :maxcomp => 1,
        :compsize => 10,
        :spawn_rate => 0.01,
        :agent_step! => galaxy_agent_step_spawn_at_rate!,
        :compmix_func => horizontal_gene_transfer,
        :mutation_rate => [0,0.1],
        :n_idxs_to_keep_from_destination => [0,1,2,5],
        :compatibility_func => planets_in_range,
        :r => [
            9,
            10,
            11,
            12,
            20,
            50,
            174],
        :destination_func => most_similar_planet,
        :rng => MersenneTwister(3141)
        )

end

"""
    write_parameters_to_json(parameters, n)

Write the parameters to a JSON and YAML file in the output directory

# Arguments
- `parameters::Dict`: A dictionary of parameters to write to file
- `n::Int`: The number of steps to run the simulation for

# Returns
- `pathname::String`: The path to the directory where the parameters were written
"""
function write_parameters_to_json(parameters, n)
    pkg_path = joinpath(splitpath(pathof(TerraformingAgents))[1:end-2])
    curr_time = Dates.format(now(),"YYYY-mm-dd-HHMMSS")
    pathname = joinpath(pkg_path,"output",curr_time)
    mkpath(pathname)
    parameters_copy = deepcopy(OrderedDict(parameters))
    parameters_copy[:compmix_func] = string(parameters_copy[:compmix_func])
    parameters_copy[:compatibility_func] = string(parameters_copy[:compatibility_func])
    parameters_copy[:destination_func] = string(parameters_copy[:destination_func])
    parameters_copy[:agent_step!] = string(parameters_copy[:agent_step!])
    parameters_copy[:n] = n
    parameters_copy = sort(parameters_copy)
    # Write JSON of parameters to file
    open(joinpath(pathname,"parameters.json"),"w") do f
        JSON.print(f, parameters_copy, 4) 
    end
    YAML.write_file(joinpath(pathname,"parameters.yml"), parameters_copy)
    return pathname 
end
pathname = write_parameters_to_json(parameters, n)

@everywhere begin
    """
        hgt_paramscan_setup(;kwargs...)

    Setup function for the parameter scan. This function is called by the `paramscan` 
        function in Agents.jl

    NOTE: It doesn't matter what these default values are because they will all be 
          overwritten by the parameters argument when paramscan is called

    # Returns
    - `model::GalaxyModel`: The model to run the simulation on
    """
    function hgt_paramscan_setup(;
        nplanets = 100,
        extent = (100,100,100),
        nool = 1,
        dt = 10,
        maxcomp = 1,
        compsize = 4,
        spawn_rate = 0.01,
        agent_step! = galaxy_agent_step_spawn_at_rate!,
        compmix_func = horizontal_gene_transfer,
        mutation_rate = 0,
        n_idxs_to_keep_from_destination = 2,
        compatibility_func = compositionally_similar_planets,
        allowed_diff = nothing,
        r = nothing,
        destination_func = nearest_planet,
        rng = MersenneTwister(3141)
        )
        ## 

        compmix_kwargs = Dict(
            :mutation_rate=>mutation_rate, 
            :n_idxs_to_keep_from_destination=>n_idxs_to_keep_from_destination)

        if compatibility_func == compositionally_similar_planets
            compatibility_kwargs = Dict(:allowed_diff=>allowed_diff)
        elseif compatibility_func == planets_in_range
            compatibility_kwargs = Dict(:r=>r)
        end

        galaxyparams = GalaxyParameters(
            rng,
            nplanets,
            extent = extent,
            nool = nool,
            dt = dt,
            maxcomp = maxcomp,
            compsize = compsize,
            spawn_rate = spawn_rate,
            compmix_func = compmix_func,
            compmix_kwargs = compmix_kwargs,
            compatibility_func = compatibility_func,
            compatibility_kwargs = compatibility_kwargs,
            destination_func = nearest_planet)
        
        model = galaxy_model_setup(galaxyparams)

        model

    end

    ## Functions specific to agent data collection. 
    ##  The main purpose is to turn references to objects into references to object ids
    destination_id(A)::Int = A.destination.id
    ancestor_ids(A)::Vector{Int} = [i.id for i in A.ancestors]
    adata =  [:pos,
            :composition, # property of Planet and Life
            :alive,
            :claimed,
            :destination_distance, # Life
            destination_id, # Life
            ancestor_ids] # Life
    adata = filter(x -> Symbol(x) ∉ keys(parameters), adata)

    ## Functions specific to model data collection.
    mantel_corr_coeff(model) = TerraformingAgents.PlanetMantelTest(model; permutations=0)[1] ## Permutations set to 0 because mdata won't be able to show this p-value anyways
    mantel_p_value(model) = TerraformingAgents.PlanetMantelTest(model; permutations=99)[2] ## 99 permutations is relatively low, but will significantly speed up the simulation
    mutation_rate(model) = model.properties[:compmix_kwargs][:mutation_rate]
    n_idxs_to_keep_from_destination(model) = model.properties[:compmix_kwargs][:n_idxs_to_keep_from_destination]
    allowed_diff(model) = model.properties[:compmix_kwargs][:allowed_diff]
    r(model) = model.properties[:compmix_kwargs][:r]
    mdata = [:dt,
            :lifespeed, 
            :interaction_radius, 
            :nplanets, 
            :maxcomp, 
            :compsize,
            :spawn_rate,
            :compmix_func,
            :n_living_planets,
            :terraformed_on_step,
            :n_terraformed_on_step,
            :compatibility_func,
            :destination_func,
            mutation_rate,
            n_idxs_to_keep_from_destination,
            mantel_corr_coeff,
            mantel_p_value]
    if parameters[:compatibility_func] == compositionally_similar_planets
        push!(mdata, allowed_diff)
    elseif parameters[:compatibility_func] == planets_in_range
        push!(mdata, r)
    end
    mdata = filter(x -> Symbol(x) ∉ keys(parameters), mdata)

    ## the paramscan function requires the step functions to be named exactly agent_step! and model_step!
    agent_step! = parameters[:agent_step!]
    model_step! = galaxy_model_step!

    ## only collect data on steps where planets are terraformed
    terraformed_on_step(model,s) = model.terraformed_on_step

end

df_agent, df_model = paramscan(
                        parameters, 
                        hgt_paramscan_setup; 
                        agent_step!, 
                        model_step!, 
                        n,
                        adata=adata, 
                        mdata=mdata,
                        when=terraformed_on_step,
                        when_model=terraformed_on_step,
                        showprogress=true,
                        obtainer=deepcopy,
                        include_constants = false,
                        parallel=true)

## Agents.jl combines all all agent data in a single dataframe, but we want to separate 
##  "planets" from "life", as they have many non-overlapping and empty properties that 
##  makes data processing cumbersome
df_planets, df_lifes = split_df_agent(df_agent, length(parameters[:extent]), parameters[:compsize])

save(File(format"CSV", joinpath(pathname, "df_planets.csv.gz")), clean_df(df_planets), delim=';')
save(File(format"CSV", joinpath(pathname, "df_lifes.csv.gz")), clean_df(df_lifes), delim=';')
save(File(format"CSV", joinpath(pathname, "df_model.csv.gz")), clean_df(df_model), delim=';')

