include("calibration.jl")

@assert length(ARGS) == 1
# Load settings and catchments from input TOML file
settings = from_toml(SettingsCalibration, ARGS[1])
catchments = readCatchmentList(settings.path_catchments_list)
# Loop through catchments
for (n, id) in enumerate(catchments)
    folder = joinpath(dirCatchment(id, settings), "calibrated")
    parameters = ParameterSet(joinpath(folder, "parameters.csv"))
    for k in ("Timeresinsec", "GscInt")
        setValue!(parameters, k, getValue(parameters, k) / 3.)
    end
    toCSV(parameters, joinpath(folder, "rescaled_parameters_1h.csv"))
end
