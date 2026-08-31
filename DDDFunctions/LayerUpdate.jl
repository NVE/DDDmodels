#    Function Layer_update
#
#------------------------------------------------------------------------
#     Description:  Updates Layers by todays event and shifting to one timestepahead. 
#
#     Author: Thomas Skaugen
#     Revised: 17.12.2019
#--------------------------------------------------------------------------

"""
LayerUpdate!(ddist, outx, Layers, layerUH, nodaysvector, NoL)

Update the subsurface water storage `Layers`.

# Arguments
- `ddist`: weights distributing soil moisture to the layers [-]
- `outx`: moisture to be distributed [mm]
- `Layers`: subsurface water storage (columns: time steps; rows: layers) [mm]
- `layerUH`: weights distributing moisture in time (columns: time steps; rows: layers) [-]
- `nodaysvector`: number of time steps for each layer [-]
- `NoL`: number of layers [-]
"""
function LayerUpdate!(ddist, outx, Layers, layerUH, nodaysvector, NoL)
  @inbounds for j in 1:NoL
    multiplier = ddist[j] * outx
    Layers[1,j] = multiplier * layerUH[1,j]
    @simd for h in 2:nodaysvector[j]
      Layers[h-1,j] = Layers[h,j] + multiplier * layerUH[h,j]
    end
  end
end
