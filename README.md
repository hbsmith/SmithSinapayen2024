This repository contains the R and Julia scripts used to generate the results and create the plots found in Smith & Sinapayen, 2024 (https://arxiv.org/abs/2403.14195).

The simulation data is available without rerunning on figshare: https://figshare.com/articles/dataset/2024-03-12-175828_zip/25515067

## Running Julia `\scripts`

- Install Julia (https://julialang.org/downloads/)
- Clone this repo (or a forked version) to your machine
- Navigate to the repo directory, and into `\scripts`
- Start the Julia REPL from your terminal (`julia`)


```julia
julia> 
```

- Access Julia's package mode by hitting the `]` key. 

  ```julia
  julia> ]
  ```

  > *Note:* You will notice your promt change from `julia>` to `(@v1.6) pkg> `, or whatever your version of Julia is.

  ```julia
  (@v1.6) pkg> 
  ```

  > *Note*: This is Julia's powerful Pkg manager, you can learn more about it [here](https://docs.julialang.org/en/v1/stdlib/Pkg/).

  > *Note*: You can also access `shell` from within Julia by hitting `;` instead of `]`.

- Next we will `activate` our environment. This is kind of like a virtual environment within python--this allows for us to install packages at specific versions for this particular project without affecting other projects.

  ```julia
  (@v1.6) pkg> activate .
  ```

  > *Note:* You'll notice your prompt switch to the directory/repo name, to let you know you're in this environment. 

  ```julia
  (scripts) pkg>
  ```

- `add` the `TerraformingAgents` package via it's GitHub URL and verson number. 

  ```julia
  (scripts) pkg> add https://github.com/hbsmith/TerraformingAgents#0.1.0
  ```

  > *Note:* This is the only package that must be manually installed, because it is not registered in the Julia package registry. 

- Then we `instantiate`. This means that Julia looks at the packages required by the project (inside the `Project.toml` file), and installs all of them.

  ```julia
  (scripts) pkg> instantiate
  ```

Now you're ready to run the two scripts:

### `run_simulation_parallel.jl`

This generates the data used in all analyses in Smith & Sinapayen, 2024. By default, the data will create a directory called `output` in a directory where the `TerraformingAgents` package is. On macOS, this will be somewhere like, `/Users/me/.julia/packages/TerraformingAgents/WyHKo/output`. Note that it can take a few hours to run depending on the number of cores allocated, and the total data size will be ~2.3GB. The already generated data is available on figshare: https://figshare.com/articles/dataset/2024-03-12-175828_zip/25515067

To run, simply navigate to the `/scripts` directory and:

``` julia
% julia
julia> include("run_simulation_parallel.jl")
```

### `NplanetsVsMantelAtSTD.jl`

This generates the Figure A3 in Smith & Sinapayen, 2024. By default, the plot will be output to the system GUI.

```julia
% julia
julia> include("NplanetsVsMantelAtSTD.jl")
julia> main()
```

## Running R `\analysis`

...
