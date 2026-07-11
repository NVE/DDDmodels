using Test
using Random
include(joinpath(dirname(@__DIR__), "calibration", "calibration.jl"))

# Input: load parameters and path to PTQ file
dir_data = joinpath(@__DIR__, "data")
parameters = ParameterSet(joinpath(dir_data, "parameters_12.70.csv"))
path_ptq = joinpath(dir_data, "ptq_12.70_2009-2024.csv")
# Output paths
dir_out = mkpath(joinpath(dir_data, "output"))
println("Test output in ", dir_out)
path_out_series = joinpath(dir_out, "output_series_12.70.csv")
path_out_r2 = joinpath(dir_out, "output_r2_12.70.csv")
# Load PTQ
timesteps, precipitation, temperature, discharge = loadPTQ(path_ptq)
# Load output and benchmark time series
output = CSV.read(path_out_series, DataFrame)
benchmark = CSV.read(joinpath(dir_data, "benchmark_series_12.70.csv"), DataFrame)
# Test
@testset "Check simulated output against benchmark for 12.70" verbose=true begin
    Random.seed!(0)
    ddd(timesteps, precipitation, temperature, discharge, 1, getHydrologicParameters(parameters),
        parameters.values, path_out_series, path_out_r2, 0, 0, 0, 365, true)
    @test isapprox(output, benchmark)
end
