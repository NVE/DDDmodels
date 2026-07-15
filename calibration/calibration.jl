using Infiltrator
using Configurations
using TOML
using CSV
using Random
using Dates
using DataFrames
using Glob
using StaticArrays
using OrderedCollections
using BlackBoxOptim
using Plots
default(legend=false)
gr() # remove to plot interactively
ENV["GKSwstype"] = "100" # remove to plot interactively
include(joinpath(dirname(@__DIR__), "DDDFunctions", "DDDAllTerrain22012024.jl"))

@option struct SettingsCalibration
    root_output::String
    path_catchments_list::String
    path_parameter_ranges::String
    template_path_ptq::String
    template_path_inipar::String
    spinup::Int
    steps_max::Int
    periods::Dict{String,String}
end

mutable struct ParameterSet
    names::Vector{String}
    values::Vector{Float64}
    positions_hyd::Vector{UInt8}
    function ParameterSet(path_in::String)
        raw = CSV.read(path_in, DataFrame, header=["Name", "val"], delim=';')
        positions_hyd::Vector{UInt8} = [20, 21, 22, 18, 19, 33, 34, 35, 36, 37]
        new(raw[:,"Name"], raw[:,"val"], positions_hyd)
    end
end

function getHydrologicParameters(parameters::ParameterSet)::Vector{Float64}
    parameters.values[parameters.positions_hyd]
end

function setHydrologicParameters!(parameters::ParameterSet, hydrologic::Vector{Float64})
    parameters.values[parameters.positions_hyd] .= hydrologic
end

function toCSV(parameters::ParameterSet, path::String)
   CSV.write(path, zip(parameters.names, parameters.values), delim=';', writeheader=false)
end

function pathsPTQ(id::String, settings::SettingsCalibration)
    Dict(k => replace(settings.template_path_ptq, "<CATCHMENT>" => id, "<PERIOD>" => v) for (k, v) in settings.periods)
end

function pathIniPar(id::String, settings::SettingsCalibration)
    replace(settings.template_path_inipar, "<CATCHMENT>" => id)
end

function dirCatchment(id::String, settings::SettingsCalibration)
    joinpath(settings.root_output, id)
end

function pathDone(id::String, settings::SettingsCalibration)
    joinpath(dirCatchment(id, settings), "done")
end

function checkInput(id::String, settings::SettingsCalibration)
    paths = pathsPTQ(id, settings)
    paths["inipar"] = pathIniPar(id, settings)
    Dict(vcat(["id" => id],[k => isfile(p) for (k, p) in paths]))
end

function loadParameterRanges(id::String, settings::SettingsCalibration)::OrderedDict{String,Tuple{Float64,Float64}}
    names = ("u", "pro", "TX", "pkorr", "skorr", "GscInt", "OFVP", "OFVIP", "Lv", "Rv")
    raw = TOML.parsefile(settings.path_parameter_ranges)
    bounds = Dict(k => convert(Dict{String,SVector{2,Float64}}, d) for (k, d) in raw)
    if haskey(bounds, id)
        return OrderedDict(k => Tuple(k in keys(bounds[id]) ? bounds[id][k] : bounds["default"][k]) for k in names)
    else
        return OrderedDict(k => Tuple(bounds["default"][k]) for k in names)
    end
end

function runDDD(paths_ptq::Dict{String,String}, params_hyd::Vector{Float64}, params_all::Vector{Float64}, spinup::Int, dir_out::String, txt::String)
    print("\tRuns with ", txt)
    kge_score = Dict(k => NaN for k in eachindex(paths_ptq))
    Threads.@threads for p in collect(eachindex(paths_ptq))
        path_out_series = joinpath(dir_out, "series_$(p).csv")
        path_out_r2 = joinpath(dir_out, "r2_$(p).csv")
        timesteps, precipitation, temperature, discharge = loadPTQ(paths_ptq[p])
        Random.seed!(0)
        kge_score[p] = ddd(timesteps, precipitation, temperature, discharge, 1, params_hyd,
                           params_all, path_out_series, path_out_r2, 0, 0, 0, spinup, true)[3]
    end
    println(" -> KGE (period): ", join(["$(round(v, digits=4)) ($k)" for (k, v) in kge_score], ", "))
end

function makeEvaluator(timesteps::Vector{DateTime}, precipitation::Matrix{Float64}, temperature::Matrix{Float64},
                       discharge::Vector{Float64}, parameters_all::Vector{Float64}, spinup::Int, path_r2::String)
    function wrapper(hydpar::Vector{Float64})
        Random.seed!(0)
        kge_score = ddd(timesteps, precipitation, temperature, discharge,
                        1, hydpar, parameters_all, "", path_r2, 0, 0, 1, spinup, true)[3]
        return 1. - kge_score
    end
end

function readCatchmentList(path)
    catchments = readlines(path) .|> strip .|> String
    filter(l -> !isempty(l) && !startswith(l, "#"), catchments)
end

function plotParameters(id::String, pct_keep::Int, dir_fig::String, settings::SettingsCalibration)
    # 1. Load parameter ranges
    bounds = loadParameterRanges(id, settings)
    # 2. Identify actually calibrated parameters
    names_par = [k for (k, r) in bounds if r[1] != r[2]]
    num_par = length(names_par)
    # 3. Load parameter sets and KGE values, and keep the best as per pct_keep
    paths_r2 = glob(joinpath(joinpath(settings.root_output, id, "calibrated"), "log", "r2_*.csv"))
    filter!(p -> filesize(p) > 0, paths_r2)
    header = vcat(["NSE", "KGE", "Bias"], ["u", "pro", "TX", "pkorr", "skorr", "GscInt", "OFVP", "OFVIP", "Lv", "Rv"])
    df = sort(CSV.read(paths_r2, DataFrame, header=header)[:,vcat(["KGE"], names_par)], "KGE")
    df = df[end-Int(floor(nrow(df) * pct_keep/100)):end,:]
    mm_KGE = round.(extrema(df[:,"KGE"]), digits=4)
    # 4. Matrix plots
    plots = Matrix{Plots.Plot}(undef, num_par, num_par);
    for c in 1:num_par
        p_x = names_par[c]
        for r in 1:num_par
            xt = r==num_par ? bounds[p_x] : false
            xl = r==num_par ? p_x : ""
            if c == r
                plots[r,c] = scatter(df[1:end-1,p_x], df[1:end-1,"KGE"], color="black", msw=0, grid=false,
                                     xticks=xt, yticks=false, xlim=bounds[p_x], xlabel=xl, ymirror=true);
                scatter!(plots[r,c], [df[end,p_x]], [df[end,"KGE"]], color="red", msw=0);
            elseif c < r
                p_y = names_par[r]
                yt = c==1 ? bounds[p_y] : false
                yl = c==1 ? p_y : ""
                plots[r,c] = scatter(df[:,p_x], df[:,p_y], marker_z=df[:,"KGE"], color=:viridis, msw=0, grid=false,
                                     xticks=xt, yticks=yt, xlim=bounds[p_x], ylim=bounds[p_y],
                                     xlabel=xl, ylabel=yl);
                scatter!(plots[r,c], [df[end,p_x]], [df[end,p_y]], color="red", marker=:cross, msw=0.5);
            else
                plots[r,c] = plot(showaxis=false, grid=false, framestyle=:none, background_color=:white, border=:none);
            end
        end
    end
    fig = plot(permutedims(plots, (2, 1))..., layout=(num_par, num_par), size=(1000, 1000),
               plot_title="$(id): best $(pct_keep)% KGE ∈ [$(mm_KGE[1]), $(mm_KGE[2])]");
    savefig(fig, joinpath(dir_fig, "KGE_parameters_$(id).pdf")) # to display interactively: gui(fig)
end

function calibrateMultipleCatchments(path_toml::String)
    # Load settings from TOML file
    settings = from_toml(SettingsCalibration, path_toml)
    num_threads = max(Threads.nthreads() - 1, 1)
    # Load catchment list and keep only those not done yet
    catchments = readCatchmentList(settings.path_catchments_list)
    num_tot = length(catchments)
    filter!(id -> !isfile(pathDone(id, settings)), catchments)
    if num_tot > length(catchments)
        println("Of ", num_tot, " catchments, ", num_tot - length(catchments), " are already calibrated")
    end
    # Check that input files are valid
    ok_input = DataFrame(checkInput.(catchments, Ref(settings)))
    is_bad = (!all).(eachrow(ok_input[:,setdiff(names(ok_input), ["id"])]))
    if any(is_bad)
        error(println("Some input files do not exist or are not valid:\n", ok_input[is_bad,:]))
    end
    # Loop through catchments
    for (n, id) in enumerate(catchments)
        ## Root folder for catchment output
        dir_out = mkpath(dirCatchment(id, settings))
        println("\nCATCHMENT ", id, " (", n, " of ", length(catchments), "): output in ", dir_out)
        ## Paths to PTQ input for each period
        paths_ptq = pathsPTQ(id, settings)
        ## Load initial parameters and run DDD
        path_inipar = pathIniPar(id, settings)
        parameters = ParameterSet(path_inipar)
        dir_out_ini = mkpath(joinpath(dir_out, "initial"))
        runDDD(paths_ptq, getHydrologicParameters(parameters), parameters.values, settings.spinup, dir_out_ini, "initial parameters")
        ## Load parameter ranges
        ranges = collect(values(loadParameterRanges(id, settings)))
        ## Calibrate
        dir_out_cal = mkpath(joinpath(dir_out, "calibrated"))
        dir_log_cal = mkpath(joinpath(dir_out_cal, "log"))
        template_path_r2 = joinpath(dir_log_cal, "r2.csv")
        timesteps, precipitation, temperature, discharge = loadPTQ(paths_ptq["calibration"])
        evaluator = makeEvaluator(timesteps, precipitation, temperature, discharge,
                                  parameters.values, settings.spinup, template_path_r2)
        print("\tCalibration started on ", Dates.format(now(), "yyyy-mm-dd HH:MM"))
        res = redirect_stdio(stdout=devnull, stderr=devnull) do
            bboptimize(evaluator; SearchRange=ranges, MaxSteps=settings.steps_max, TraceMode=:silent, SaveTrace=true, NThreads=num_threads)
        end
        print(" and ended after ", canonicalize(Second(Int(round(res.elapsed_time)))))
        println(" -> KGE: ", round(1 - best_fitness(res), digits=4))
        ## Write calibrated parameters to file
        setHydrologicParameters!(parameters, best_candidate(res))
        toCSV(parameters, joinpath(dir_out_cal, "parameters.csv"))
        ## Parameter plots
        plotParameters(id, 50, dir_out, settings)
        ## Run DDD with calibrated parameters
        runDDD(paths_ptq, getHydrologicParameters(parameters), parameters.values, settings.spinup, dir_out_cal, "calibrated parameters")
        ## Create empty file ("done") to be used in case of restart to skip this catchment
        touch(pathDone(id, settings))
    end
end

function inputSingleCatchmentRun(path_toml::String, id::String, period::String)
    # Load settings from TOML file
    settings = from_toml(SettingsCalibration, path_toml)
    # PTQ input
    path_ptq = pathsPTQ(id, settings)[period]
    timesteps, precipitation, temperature, discharge = loadPTQ(path_ptq)
    # Load initial parameters
    path_inipar = pathIniPar(id, settings)
    parameters = ParameterSet(path_inipar)
    # Root folder for catchment output
    dir_out = mkpath(joinpath(settings.root_output, "single_runs", id))
    path_out_series = joinpath(dir_out, "series_$(id)_$(period).csv")
    path_out_r2 = joinpath(dir_out, "r2_$(id)_$(period).csv")
    println("Input prepared for single run with output in ", dir_out)
    # Return input required for running DDD
    return (timesteps, precipitation, temperature, discharge, 1, getHydrologicParameters(parameters),
            parameters.values, path_out_series, path_out_r2, 0, 0, 0, settings.spinup, true)
end
