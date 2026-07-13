#    Function Layer_capacity
#
#------------------------------------------------------------------------
#     Description:  Calculates current capacity in groundwater layers. 
#
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#--------------------------------------------------------------------------

function LayerCapacityUpdate(Layers, nodaysvector, Magkap, NoL)

ddistx = Vector{Float64}(undef, NoL)

#Below are the states (in mm) for each saturation level
@inbounds for j in 1:NoL
                                      #state after this timesteps' water is gone. amount of water  in mm, minus current timestep
  aktMag = sum(@view(Layers[2:nodaysvector[j],j]); init=0.0)
  ddistx[j] = max(Magkap[j] - aktMag, 0.0)
        
end

return ddistx  # ddistx informs on current capacity for each level in mm.
end
