#    Function Layer_capacity
#
#------------------------------------------------------------------------
#     Description:  Calculates current capacity in groundwater layers. 
#
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#--------------------------------------------------------------------------

function LayerCapacityUpdate!(ddistx, Layers, nodaysvector, Magkap, NoL)

#Below are the states (in mm) for each saturation level
@inbounds for j in 1:NoL
                                      #state after this timesteps' water is gone. amount of water  in mm, minus current timestep
  aktMag = 0.0
  @simd for i in 2:nodaysvector[j]
      aktMag += Layers[i,j]
  end

  ddistx[j] = max(Magkap[j] - aktMag, 0.0)
        
end

# ddistx informs on current capacity for each level in mm.
end
