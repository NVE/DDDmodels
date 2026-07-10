#    Function GrW_point
#
#------------------------------------------------------------------------
#     Description: Estimates moisture (groundwater) at a point in the hillslope
#
#     Author:  Thomas Skaugen
#     Revised: 18.08.2023
#--------------------------------------------------------------------------

function GrWPoint(nodaysvector, Tpkt, NoL, Areas, innLayers, totarea, distvec, delta_d)

# actual mm at the site is Layers* totarea/areas, see Eq 27 in Kreboel, 2016
# Note that Layers is a scaled verion of input, scaled according to the fractional area each j represents
# Layers_mm will be the input (current for actual area and from previous upstream areas)  
# delta_d is how wide are each box (for ecah celerity)
# distvec the distance away from the RN we want estimated

  Layersmm = zeros(nodaysvector[NoL], NoL)
  
  for i in 1: NoL
    for j in 1: nodaysvector[i]   
	   Layersmm[j,i] = innLayers[j,i]*(totarea/Areas[j,i])
    end        
  end
  
ant = length(distvec)       # number of groundwater bore hole locations
grw = zeros(ant)            # vector for storing groundwater values for the current Layer 

for i in 1:(ant-1)  # to be certain we do not go beyond the matrix
  tall = distvec[i]
  grw[i] = 0.0
  for j in 1: NoL
   grw[i] += round(Layersmm[(Int(trunc(tall/delta_d[j]))+1),j],digits=3) # Integer division.  Sums up Layer_mm for each layer
   #println("grw: ",grw[i])
  end
end
             
return grw

end
