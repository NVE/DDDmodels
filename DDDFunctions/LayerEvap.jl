#    Function Layer_evap
#
#------------------------------------------------------------------------
#     Description:  Calculates evapotranspiration directly from Layers and updates Layers. 
#
#     Author: Thomas Skaugen and Anne Stavang
#     Revised: 4.7.2019
#--------------------------------------------------------------------------

function LayerEvap!(Layers, nodaysvector, ea_S, layerUH, NoL)
  
  #Layers = [1.0 2.0 0.0 0.0 0.0; 2.0 3.0 4.0 0.0 0.0; 3.0 4.0 5.0 6.0 0.0; 2.0 2.0 2.0 2.0 0.35]
  #NoL = 4
  #nodaysvector = [2 3 4 5]
  #ea_S = 4.5
  #layerUH = [0.6 0.4 0.0 0.0 0.0; 0.4 0.35 0.25 0.0 0.0; 0.35 0.25 0.25 0.15 0.0; 0.3 0.25 0.2 0.15 0.1]

  (ea_S == 0) && return 0.0 # return zero evaporation if input ea_S is zero
   
  total_layers_last = sum(Layers)
  
  #Below are the states (in mm) for each saturation level
  redea = ea_S

  
  @inbounds for j in 1:NoL # 1 is the top(fastest) Layer NoL is the bottom layer 

        (redea <= 0) && break # exit the layers loop if there is no more water to evaporate
        
        aktMag = sum(@view(Layers[1:nodaysvector[j],j])) # this is correct because ea is a nonintegrated with a continuum variable as opposed to discharge

        (aktMag == 0) && continue # skip to next layer if it is empty
        
        differ = aktMag - redea 
  
        if differ > 0 # the Layer has more water than is to be evaporated > ea_S                
            ea_excess = 0.0
            tot_layers = 0.0
            @simd for i in 1:nodaysvector[j]
              evapUH = redea * layerUH[i,j]
              if Layers[i,j] < evapUH
                  ea_excess += evapUH
              else
                  Layers[i,j] -= evapUH
              end
              tot_layers += Layers[i,j]
            end
            avvik = aktMag - tot_layers - redea + ea_excess
            if abs(avvik) > 1e-8
              println("Hei, feil i fordampning", avvik)
            end  
            redea = 0.0
        else
          Layers[1:nodaysvector[j],j] .= 0.
          redea -= aktMag 
        end
  end
   
  total_layers_current = sum(Layers)
  ea = total_layers_last - total_layers_current
   
#  if (ea > 0) && abs(total_layers_last - (ea + total_layers_current)) > 1e-8 # COMMENTED OUT AS IT IS A CIRCULAR CHECK
#      println("Layers ut not OK")
#  end

  return ea
end
