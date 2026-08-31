#    Function OverlandFlowDynamicDD
#
#------------------------------------------------------------------------
#     Description:  Estimates a dynamic distance distribution for overland flow. 
#                   ddist is a vector of weights, distributing outx to the layers
#     Author: Thomas Skaugen and Anne E. Stavang
#     Revised: 08.03.2019
#--------------------------------------------------------------------------


function OverlandFlowDynamicDD!(k, ddist, outx, layerUH, nodaysvector, NoL, midDL, CritFlux, Timeresinsec)

# using Distributions
    
# k: 1 dim array float
# ddist: 1 dim array float
# layerUH: 1 dim array float
# nodaysvector: 1 dim array integer
# NoL: scalar, integer    
# midDL: scalar, integer
# Critflux: scalar, float
# Timeresinsec: scalar, integer    

#include("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Julia\\DDDFunctions\\SingleUH.jl") 
#NoL = 2
#midDL = 156.0
#maxDL = 928.0
#dmean = 82.3
#zS = 0.07
#Timeresinsec = 86400
#k = [0.001,0.0002]
#nodaysvector = [11, 54]
#UH1 = SingleUH(k[1],Timeresinsec, midDL, maxDL, zS)
#UH2 =SingleUH(k[2],Timeresinsec, midDL, maxDL, zS)
#layerUH =zeros(2,54)
#layerUH[1,1:nodaysvector[1]] = UH1
#layerUH[2,1:nodaysvector[2]] = UH2
    
  A2 = CritFlux/(outx*ddist[1]/1000)       # Area in m2 that generates critical volume
  #dmean <- A2^0.25# (A2^0.5)^0.5          # doesn't work
  #dmean <- ((A2+A2)^0.5)/2                # half of the diagonal Works better
  #dmean =(A2/(pi*1.25))^0.5*(1.5/2)       # approx ovoid shaped catchment. Two circles, with one double the radius of the other.dmean  is the mean of the two diameters: Works pretty good
  dmean = 0.5*A2^0.5        #The mean describes half of the critical area, as it does for the entire catchment. Kind of logical

  if (dmean >= midDL) 
    dmean = midDL                                     #degenerates to natural river network
  else                                  # Dynamic OF RN
   eksp = Exponential(dmean)     
   maxdistOF = quantile(eksp, 0.99)                   #such that 0.99 of the area contributes
   celerityOF = k[1]                                  #overland flow celerity
   zOF =  0.5165*dmean^-0.948                         # areal fraction of river network Anne's regression
      if (zOF > 1)
        zOF = 0.99
      end
      OFUH = SingleUH(celerityOF,Timeresinsec, dmean, maxdistOF, zOF)
      if length(OFUH) < nodaysvector[1] # ensures DNR is more dense than original RN
         layerUH[1:length(OFUH),1] .= OFUH
         layerUH[(length(OFUH)+1):nodaysvector[NoL],1] .= 0.0
      end
  end
end
