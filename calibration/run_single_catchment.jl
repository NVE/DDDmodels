# To run this script:
# julia --project=.. run_single_catchment.jl settings/calibration_251120.toml 2.11 calibration
include("calibration.jl")

# Shell arguments
@assert length(ARGS) == 3
path_settings, id, period = ARGS
# Load input
input = inputSingleCatchmentRun(path_settings, id, period)
#runSingleCatchment(path_settings, id, period)

#using Profile
#Profile.clear()
#@profile runSingleCatchment(path_settings, id, period)
#Profile.print(format=:tree, sortedby=:count, maxdepth=10, mincount=5)

#using InteractiveUtils
#@code_warntype runSingleCatchment(path_settings, id, period)

using BenchmarkTools
@btime ddd(input...)

# ProfileView displays works in REPL and not script
#include("calibration.jl")
#using ProfileView
#using Cthulhu
#@profview runSingleCatchment("settings/calibration_251120.toml", "2.11", "calibration")

# REPL
#using Profile, OwnTime
#include("calibration.jl")
#Profile.clear()
#@profile runSingleCatchment("settings/calibration_251120.toml", "2.11", "calibration")
#owntime()
