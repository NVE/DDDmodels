#    Function LayerEstimation
#
#------------------------------------------------------------------------
#     Description:  Calculates the total capacity of groundwater reservoir and for its Layers  
#
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#--------------------------------------------------------------------------

function LayerEstimationDesignEmpSM(GshRes,GscRes, NoL, gtcel)

# using Distributions
# include("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Julia\\DDDFunctions\\SingleUH.jl")

 # GshInt: scalar float = 1.5
 # GscInt: scalar float = 0.22 
 # midDL: scalar float = 129.84
 # maxDl: scalar float = 538
 # MAD: scalar float = 2.35
 # Timeresinsec: scalar float = 86400
 # NoL: scalar integer = 5
 # area2: scalar float = 2300000 
 # gtcel: scalar float = 0.99
 


MLev = [1/(NoL-1):1/(NoL-1):1.0;]              # (sequence)Quantiles  to calculate reservoir levels [0.1:0.1:0.9;]

MLev[NoL-1] = gtcel                            # quantile for start overland flow
Res_prob = zeros(Float64,(NoL-1))
Magkap = zeros(Float64,NoL)
g = Gamma(GshRes,GscRes) 
#calculates the reservoir levels associated with quantiles. Mean is GshRes*GscRes
Res_prob .= quantile.(g,MLev)

#Capasity of Layers
ssRes1 = zeros(Float64,NoL)
ssRes1[1] = 2000                               # capacity of overland flow level

for i in 2:(NoL-1)
  ssRes1[i] = Res_prob[NoL-i+1]-Res_prob[(NoL-i)]
end

ssRes1[NoL] = Res_prob[1]                     # capasity for the first slowest level         
                       
Magkap = ssRes1                                 # capasity for Layers
M = Res_prob[(NoL-1)]                           # Total groundwater reservoir
#println(GshRes) 
#println(GscRes)
println("M fra Subrutine ", M)
return Magkap, M
end
