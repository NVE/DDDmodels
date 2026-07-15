using ProgressMeter
include(joinpath(@__DIR__, "regression_testing.jl"))

# julia --project=.. --threads 15 create_tests.jl ../calibration/settings/calibration_251120.toml data/catchments.txt

# Shell arguments
@assert length(ARGS) == 2
path_settings, path_catchments = ARGS
# Catchment lists
catchments = readCatchmentList(path_catchments)
# Create tests
println("Create test files for $(length(catchments)) catchments:")
@showprogress Threads.@threads for id in catchments
    test = RegressionTest(id, path_settings)
    create(test)
end
