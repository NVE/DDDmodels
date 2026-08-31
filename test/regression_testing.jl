using Random
using Test
using ProgressMeter
include(joinpath(dirname(@__DIR__), "calibration", "calibration.jl"))

struct RegressionTest
    id::String
    settings::SettingsCalibration
    path_parameters::String
    path_ptq::String
    path_out_series::String
    path_out_r2::String
    path_benchmark::String
end

function RegressionTest(id::String, path_settings::String)
    settings = from_toml(SettingsCalibration, path_settings)
    dir_data = joinpath(@__DIR__, "data")
    dir_in = joinpath(dir_data, "input")
    dir_out = joinpath(dir_data, "output")
    path_parameters = joinpath(dir_in, "parameters_$(id).csv")
    path_ptq = joinpath(dir_in, "ptq_$(id).csv")
    path_out_series = joinpath(dir_out, "series_$(id).csv")
    path_out_r2 = joinpath(dir_out, "r2_$(id).csv")
    path_benchmark = joinpath(dir_data, "benchmark", "series_$(id).csv")
    RegressionTest(id, settings, path_parameters, path_ptq, path_out_series, path_out_r2, path_benchmark)
end

function create(rt::RegressionTest)
    # Forcing
    path_ptq = pathsPTQ(rt.id, rt.settings)["calibration"]
    timesteps, precipitation, temperature, discharge = loadPTQ(path_ptq)
    # Parameters
    path_par = pathIniPar(rt.id, rt.settings)
    parameters = ParameterSet(path_par)
    # Run DDD
    Random.seed!(0)
    ddd(timesteps, precipitation, temperature, discharge, 1, getHydrologicParameters(parameters),
        parameters.values, rt.path_benchmark, rt.path_out_r2, 0, 0, 0, rt.settings.spinup, true)
    # Delete r2 file (not needed for testing)
    rm(rt.path_out_r2)
    # Copy input files to test data folders
    cp(path_par, rt.path_parameters, force=true)
    cp(path_ptq, rt.path_ptq, force=true)
end

function run(tests::Vector{RegressionTest})
    # Run DDD for each catchment
    println("Running DDD on $(length(tests)) catchments:")
    @showprogress Threads.@threads for rt in tests
        # Load input
        parameters = ParameterSet(rt.path_parameters)
        timesteps, precipitation, temperature, discharge = loadPTQ(rt.path_ptq)
        # Run DDD
        Random.seed!(0)
        ddd(timesteps, precipitation, temperature, discharge, 1, getHydrologicParameters(parameters),
            parameters.values, rt.path_out_series, rt.path_out_r2, 0, 0, 0, rt.settings.spinup, true)
        # Delete r2 file (not needed for testing)
        rm(rt.path_out_r2)
    end
    # Test for each catchment
    @testset "Check simulated output against benchmark for $(length(tests)) catchments" verbose=true begin
        for rt in tests
            benchmark = CSV.read(rt.path_benchmark, DataFrame)
            output = CSV.read(rt.path_out_series, DataFrame)
            @test isapprox(output, benchmark) || "Failed for $(rt.id)"
        end
    end
end
