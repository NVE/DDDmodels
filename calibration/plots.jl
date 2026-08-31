using TOML
using ProgressMeter
include("calibration.jl")

# Run this script with
# julia --project=.. plots.jl 50 settings/calibration_251120.toml

# Check shell arguments
@assert length(ARGS) == 2
# Percentage of evaluations to keep
pct_keep = parse(Int, ARGS[1])
# Load settings from TOML file
settings = from_toml(SettingsCalibration, ARGS[2])
# Load catchment list and keep only those not done yet
catchments = readCatchmentList(settings.path_catchments_list)
# Plots
dir_fig = mkpath(joinpath(settings.root_output, "figures"))
println("Plotting parameter and KGE plots in $(dir_fig):")
@showprogress for id in catchments
    plotParameters(id, pct_keep, dir_fig, settings)
end
