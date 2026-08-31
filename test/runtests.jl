include(joinpath(@__DIR__, "regression_testing.jl"))

# julia --project=.. --threads 5 runtests.jl ../calibration/settings/calibration_251120.toml data/catchments.txt

# Shell arguments
@assert length(ARGS) == 2
path_settings, path_catchments = ARGS
# Catchment lists
catchments = readCatchmentList(path_catchments)
# Tests
tests = RegressionTest.(catchments, Ref(path_settings))
run(tests)
