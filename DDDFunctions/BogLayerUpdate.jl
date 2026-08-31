#    Function Layer_update
#
#------------------------------------------------------------------------
#     Description:  Updates Layers by todays event and shifting to one timestepahead. 
#
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#--------------------------------------------------------------------------

function BogLayerUpdate!(outbog, BogLayers, UHBog, nodaysvector)
    if nodaysvector == 1
        BogLayers .= 0.
    else
        @inbounds @simd for t in 2:nodaysvector
            BogLayers[t-1] = BogLayers[t] + outbog * UHBog[t] # finds response (mm) and shifts the matrix one timestep ahead
        end
    end
end            


