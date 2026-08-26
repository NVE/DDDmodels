#    Function CeleritySubSurface
#
#------------------------------------------------------------------------
#     Description:  Calculates Sub surface celerities from level of saturation
#
#     Author: Thomas Skaugen
#     Revised: 07.05.2024
#--------------------------------------------------------------------------

function CeleritySubSurface1UH(shLam,scLam,prob,midDL,Timeresinsec) 
#function CeleritySubSurface1UH(Gshape, Gscale, midDL, Timeresinsec, prob) 
#probvec[1] = 0.99                            #Quantile in celerity distribution for overland flow fixed at 0.99

g = Gamma(shLam, scLam)        #Gamma distributed groundwater fluctuations 
#k = quantile.(g,prob)*midDL/Timeresinsec

k = quantile.(g,prob) *midDL/Timeresinsec #Celerity [m/s]
return k                        #Celerities [m/s] 
end
