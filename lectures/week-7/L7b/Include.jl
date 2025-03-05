# setup paths -
const _ROOT = @__DIR__
const _PATH_TO_SRC = joinpath(_ROOT, "src");
const _PATH_TO_DATA = joinpath(_ROOT, "data");

# check: do we have a manifest file?
using Pkg
if (isfile(joinpath(_ROOT, "Manifest.toml")) == false) # have manifest file, we are good. Otherwise, we need to instantiate the environment
    Pkg.activate("."); Pkg.resolve(); Pkg.instantiate(); Pkg.update();
end

# load external packages
using Sundials
using Plots
using Colors
using DataFrames
using CSV
using FileIO
using JLD2

# load my codes -
include(joinpath(_PATH_TO_SRC, "DataFile.jl"));
include(joinpath(_PATH_TO_SRC, "Kinetics.jl"));
include(joinpath(_PATH_TO_SRC, "Flow.jl"));
include(joinpath(_PATH_TO_SRC, "Balances.jl"));
include(joinpath(_PATH_TO_SRC, "SolveBalances.jl"));