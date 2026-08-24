# To run this script:
# julia --project=.. run_single_catchment.jl settings/calibration_251120.toml 2.11 calibration
using BenchmarkTools
using Random
include("calibration.jl")

# Shell arguments
@assert length(ARGS) == 3
path_settings, id, period = ARGS
# Load input
input = inputSingleCatchmentRun(path_settings, id, period)
# Wrapper to make sure seed is set
function wrappersinglerun(input)
    Random.seed!(0)
    ddd(input...)
end
# Timing
@btime wrappersinglerun(input)
