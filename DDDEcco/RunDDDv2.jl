"""
Runs DDD model either in single run model or for calibration
The model itself is called as a function which calls on several functions
"""
#using CSV
using Distributions
using LsqFit
using Statistics
using Dates
using DataFrames
using Plots
using CSV
using BlackBoxOptim
using JLD2

prefix = "../"

##Preprocessing routines
include(prefix * "DDDFunctions/Big2SmallLambda.jl")
include(prefix * "DDDFunctions/CeleritySubSurface.jl")
include(prefix * "DDDFunctions/SingleUH.jl")
include(prefix * "DDDFunctions/SingleNormalUH.jl")
include(prefix * "DDDFunctions/LayerEstimation.jl")
include(prefix * "DDDFunctions/PyrAreas.jl")
include(prefix * "DDDFunctions/GrWPoint.jl")
include(prefix * "DDDFunctions/RiverPoint.jl")
include(prefix * "DDDFunctions/TemperatureVector.jl")

##EB and Snow Routines
include(prefix * "DDDFunctions/NedbEBGlac_debug04072022.jl")
include(prefix * "DDDFunctions/SnowpackTemp.jl")
include(prefix * "DDDFunctions/TempstartUpdate.jl")
include(prefix * "DDDFunctions/SmeltEBGlac_debug04072022.jl")
include(prefix * "DDDFunctions/CloudCoverGlac_debug04072022.jl")
include(prefix * "DDDFunctions/TssDewpoint.jl")
include(prefix * "DDDFunctions/SolradTransAlbedoper_hrs_debug04072022.jl")
include(prefix * "DDDFunctions/LongWaveRad_debug04072022.jl")
include(prefix * "DDDFunctions/SensibleLatHeat_debug04072022.jl")
include(prefix * "DDDFunctions/AlbedoUEB_debug04072022.jl")
include(prefix * "DDDFunctions/GroundPrecCC.jl")
include(prefix * "DDDFunctions/SnowGamma.jl")
include(prefix * "DDDFunctions/Varc.jl")
include(prefix * "DDDFunctions/NewSnowDensityEB.jl")
include(prefix * "DDDFunctions/NewSnowSDEB.jl")
include(prefix * "DDDFunctions/DensityAge.jl")

#Subsurface and Evaporation routines
include(prefix * "DDDFunctions/LayerCapacityUpdate.jl")
include(prefix * "DDDFunctions/PotentialEvapPT.jl")
include(prefix * "DDDFunctions/UnsaturatedEvapEB.jl")
include(prefix * "DDDFunctions/LayerEvap.jl")
include(prefix * "DDDFunctions/UnsaturatedExEvap.jl")
include(prefix * "DDDFunctions/WetlandsEB.jl")
include(prefix * "DDDFunctions/GrvInputDistributionICap2022.jl")
include(prefix * "DDDFunctions/OFICap.jl")
include(prefix * "DDDFunctions/LayerUpdate.jl")
include(prefix * "DDDFunctions/BogLayerUpdate.jl")
include(prefix * "DDDFunctions/RiverUpdate.jl")
## Overland Flow routine
include(prefix * "DDDFunctions/OverlandFlowDynamicDD.jl")
## Efficiency criteria
include(prefix * "DDDFunctions/NSEJM.jl")
include(prefix * "DDDFunctions/KGEJM.jl")
# Model Module
#include("F:/HB/HB-modellering/DDDtestbenk/DDD_Julia/DDDFunctions/DDDUrbanFunc.jl")
include(prefix * "DDDFunctions/DDDAllTerrain22012024.jl")
########################################################################################

catchment = "56.1"  # stationnumber
inputPrefix = "./input/"
outputPrefix = "./output/"
parameterPrefix = "./parameter/"

TR = "5min"         # this is just a marker for naming files, does NOT set the temporal resolution

ptqfile = string(inputPrefix, "S1_", catchment, "_", TR, "_ptq_kal.csv")

r2fil = string(outputPrefix, "r2_", catchment, "_", TR, "_kal.csv")

utfile  = string(outputPrefix, "simres_", catchment, "_", TR, "_kal_DDDv2.csv")

paramfile = string(parameterPrefix, "ParDDDv2_", catchment, "_", TR, ".csv")

println(ptqfile, "\n", r2fil, "\n", utfile, "\n", paramfile)

spinup = (31*4) #days used to spin up the model. 

prm = CSV.read(paramfile,DataFrame,header=["Name", "val"], delim=';')
#prm = CSV.read(paramfile,header=["Name", "val"], delim=';')

#            u,          pro           TX,         Pkorr        skorr,     GscInt      OVP          
#         OVIP       Lv            rv        
tprm = [prm.val[20], prm.val[21], prm.val[22], prm.val[18], prm.val[19],prm.val[33], prm.val[34], 
    prm.val[35],prm.val[36],prm.val[37]]
println(tprm)


Gshape, Gscale = Big2SmallLambda(prm.val[32], prm.val[33]) # Coverting integrated celerity to layers takes too long in calibration: preprocessing
Gpar = [Gshape, Gscale]


startsim = 1 
kal = 0
modstate = 0
savestate = 0

t1= time_ns()

function calib_wrapper_model(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil, modstate, savestate, kal, spinup)
 qobs, qberegn, KGE, NSE, bias = DDDAllTerrain(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil, modstate, savestate,
        kal, spinup)  
 return qobs,qberegn, KGE,NSE,bias 
end

function calib_single_wsh(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil, modstate, savestate, kal, spinup)
 qobs, qberegn, KGE, NSE, bias = DDDAllTerrain(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil, modstate, savestate,
        kal, spinup)    
 return (1.0 - KGE)
end

if(kal == 0)
    qobs,qberegn,KGE,NSE, bias = calib_wrapper_model(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil,
        modstate, savestate,kal, spinup) # a single run 
    
    println(catchment)
    println("KGE=",round(KGE,digits=3))
    println("NSE=",round(NSE,digits=3))
    println("bias=",round(bias,digits=3))
end

if(kal == 1) # calibrate
    #                   u,        pro,         TX,        Pkorr,    skorr,          GscInt,         OVP     OVIP 
    param_range = [(1.0,3.0), (0.05,0.05), (-0.5, 0.5), (0.5, 2.0), (0.5,2.0), (0.065,0.075), (tprm[7],tprm[7]),
        (tprm[8],tprm[8]), (tprm[9],tprm[9]),(tprm[10],tprm[10])] # 
    
    println(param_range)
    calib_single_wsh_tmp(param) = calib_single_wsh(Gpar,startsim, param, prm, ptqfile, utfile, r2fil,
                                           modstate, savestate, kal, spinup)
    res = bboptimize(calib_single_wsh_tmp; SearchRange = param_range, MaxSteps = 1000, TraceMode = :verbose)
    param_hydro = best_candidate(res)
    println(param_hydro)
end

t2 = time_ns()
println("Pkorr=", round(tprm[4],digits=3))
println("Time elapsed[s]= ",(t2-t1)/1.0e9)

if(kal==0)
 plot(qobs[15240:16240], color="black",label = "Observed",lw =1)
 plot!(qberegn[15240:16240],label = "Simulated", lw = 2)
end 
