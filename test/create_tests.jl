using ProgressMeter
include(joinpath(@__DIR__, "regression_testing.jl"))

# julia --project=.. --threads 10 create_tests.jl ../calibration/settings/calibration_251120.toml 2.11 247.3

# Shell arguments
@assert length(ARGS) > 1
path_settings = ARGS[1]
catchments = ARGS[2:end]
# Create tests
println("Create test files for $(length(catchments)) catchments:")
@showprogress Threads.@threads for id in catchments
    test = RegressionTest(id, path_settings)
    create(test)
end
