# Calibration tools for DDD


## 1. Requirements

Julia version and libraries are defined in `Project.toml` and `Manifest.toml`:
See [Pkg documentation](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project)
for how to recreate the Julia environment to calibrate DDD.


## 2. How to set up and run a calibration

To calibrate DDD on a set of catchments using 10 parallel threads, run

`julia --project=.. --threads 11 calibrate_multiple_catchments.jl <path/to/settings.toml>`

where the last argument is the path to a TOML file (for example `settings/calibration_251120.toml`)
containing the following entries:

- `root_output`: absolute path to the folder under which all calibration output will be stored.
- `path_catchments_list`: path to the text file with the IDs of the catchments to be calibrated,
  where lines starting with `#` are ignored (example: `settings/catchments/catchments_251120.txt`).
- `path_parameter_ranges`: path to the TOML file defining upper and lower bounds for the hydrologic
  parameters to be calibrated, where the `[default]` section is compulsory; and sections for
  individual catchments can be added by using catchment IDs for headers, where custom ranges for
  individual parameters can be set
  (example: `settings/parameter_ranges/parameter_ranges_251120.toml`).
- `template_path_ptq`: template path to the PTQ (precipitation, temperature and discharge) CSV file
  for each catchment, where `<CATCHMENT>` and `<PERIOD>` are replaced by the catchment ID and
  the abbreviation for each period (e.g. `kal` for calibration) in the function `pathsPTQ`.
- `template_path_inipar`:template path to the files with initial values for the full parameter set
  defining each catchment (example: files in `/hdata/hmdata01/DDD_calibration/BestPrimo2025/`).
- `spinup`: model spinup time (days).
- `steps_max`: maximum number of steps used by the optimisation algorithm.
- `periods`: section defining names (arbitrary) and abbreviations (must match the tokens used to
  replace `<PERIOD>` in `template_path_ptq`) for each period the model is run on.


## 3. Output

For each catchment, results are written in folder `root_output/<CATCHMENT>/`, where the `initial/`
and `calibrated/` subfolders contains output from DDD runs with initial and calibrated parameter
values, respectively.

In each subfolder, the following files are written by
`../DDDFunctions/DDDAllTerrain22012024.jl`:

- `r2_<PERIOD>.csv`: row file with NSE, KGE, bias and hydrologic parameter values.
- `series_<PERIOD>.csv`: DDD output time series.

Additional items in the `calibrated/` subfolder:

- `parameters.csv`: new full parameter set defining the catchment, including the newly calibrated
  hydrologic parameters.
- `log/`: all parameter sets, and corresponding scores, tried by the optimisation algorithm and
  logged in separate files (one file per computational thread to avoid concurrent access).
  Each row has the same format as the `r2_<PERIOD>.csv` files described above.


## 4. Restart

The presence of the empty file `root_output/<CATCHMENT>/done` means that calibration ran
until the stopping criterion for that catchment: If the calibration script is re-run, this
catchment will be skipped.

## 5. Possible improvements

- avoid reversing STempvec and optimise related calculations
- Check that functions are type-stable to avoid performance losses.
- improve upon `ptqinn` data structure for faster data access in time loop (now time steps are rows)
- Test alternatives to DataFrames to ensure type stability (e.g. DimensionalData, TypedTables, TimeArrays)
  or is it enough with function barrier on loading?
- Specify which parameters should be calibrated and which have a fixed value in settings,
  instead of using collapsed ranges, to reduce the formal number of parameters.
- Sensitivity analysis (e.g. [GlobalSensitivity.jl](https://docs.sciml.ai/GlobalSensitivity/stable/))
- Use [workspaces](https://pkgdocs.julialang.org/dev/creating-packages/#Test-specific-dependencies)
   to define test- and calibration-specific dependencies.

## 6. Questions

- NewSnowDensityEB.jl: conversion to Fahrenheit (overwriting) propagates to DensityAge.jl via NewSnowSDEB.jl:
  is it intentional?
- `tempstart` is not necessary to run the model as it can be replaced by `temperature[:,(i-len+1):i]` when
  computing `skintempsnow`. Is it needed to restart the model from a state file? Can it be avoided by loading
  temperature time series further back in time?
- LayerEvap.jl: reset evapUH to 0 where evapUH >= Layers. Why not Layers set to 0?
