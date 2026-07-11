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

ant = length(distvec)       # number of groundwater bore hole locations
grw = zeros(ant)            # vector for storing groundwater values for the current Layer 

for i in 1:(ant-1)  # to be certain we do not go beyond the matrix
  tall = distvec[i]
  for j in 1:NoL
   h = Int(trunc(tall/delta_d[j])) + 1
   layersmm = innLayers[h,j] * (totarea/Areas[h,j])
   grw[i] += round(layersmm, digits=3) # Sums up Layer_mm for each layer
  end
end
             
return grw

end
