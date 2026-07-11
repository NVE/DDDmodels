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
- `ddist::Vector{Float64}`: weights distributing soil moisture to the layers [-]
- `outx::Float64`: moisture to be distributed [mm]
- `Layers::Matrix{Float64}`: subsurface water storage (columns: time steps; rows: layers) [mm]
- `layerUH::Matrix{Float64}`: weights distributing moisture in time (columns: time steps; rows: layers) [-]
- `nodaysvector::Vector{Int}`: number of time steps for each layer [-]
- `NoL::Int`: number of layers [-]
"""
function LayerUpdate!(ddist::Vector{Float64}, outx::Float64, Layers::Matrix{Float64}, layerUH::Matrix{Float64}, nodaysvector::Vector{Int}, NoL::Int)
  @simd for j in 1:NoL
    multiplier = ddist[j] * outx
    Layers[1,j] = multiplier * layerUH[1,j]
    @inbounds for h in 2:nodaysvector[j]
      Layers[h-1,j] = Layers[h,j] + multiplier * layerUH[h,j]
    end
  end
end
