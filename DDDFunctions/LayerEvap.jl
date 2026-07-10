#    Function Layer_evap
#
#------------------------------------------------------------------------
#     Description:  Calculates evapotranspiration directly from Layers and updates Layers. 
#
#     Author: Thomas Skaugen and Anne Stavang
#     Revised: 4.7.2019
#--------------------------------------------------------------------------

function LayerEvap(Layers, nodaysvector, ea_S, layerUH, NoL)
  
  #Layers = [1.0 2.0 0.0 0.0 0.0; 2.0 3.0 4.0 0.0 0.0; 3.0 4.0 5.0 6.0 0.0; 2.0 2.0 2.0 2.0 0.35]
  #NoL = 4
  #nodaysvector = [2 3 4 5]
  #ea_S = 4.5
  #layerUH = [0.6 0.4 0.0 0.0 0.0; 0.4 0.35 0.25 0.0 0.0; 0.35 0.25 0.25 0.15 0.0; 0.3 0.25 0.2 0.15 0.1]
   
  total_layers_last = sum(Layers)
  
  #Below are the states (in mm) for each saturation level
  redea = ea_S
   #println(total_layers_last)
   #println(sum(Layers))

  
  for j in 1:NoL # 1 is the top(fastest) Layer NoL is the bottom layer 

    if redea > 0

      if sum(Layers) > 0
        aktMag = sum(Layers[1:nodaysvector[j],j]) # this is correct because ea is a nonintegrated with a continuum variable as opposed to discharge
        differ = aktMag - redea 
        # println(aktMag)
        # println(differ)
        # println(newLayer)
        # println(total_layers_last)
        # println(sum(Layers))

  
        if differ > 0 # the Layer has more water than is to be evaporated > ea_S                
            ea_excess = 0.0
            @views evapUH = redea .* layerUH[1:nodaysvector[j],j]
            @views newLayer = Layers[1:nodaysvector[j],j] .- evapUH[1:nodaysvector[j]]  

            tull = findall(newLayer .< 0.0)
            if length(tull) > 0
              x = tull    # locate which boxes have not enough water to evaporate
              ea_excess = sum(evapUH[x]) # the amount which is not evaporated
              evapUH[x] .= 0.0            # the identified boxes have zero instead of negative values
              @views newLayer = Layers[1:nodaysvector[j],j] .- evapUH[1:nodaysvector[j]] # updates so that everthing is positive
            end   

            if round(aktMag - sum(newLayer) - redea + ea_excess; digits=8) != 0
              avvik = round(aktMag - sum(newLayer) - redea + ea_excess; digits= 8)
              println("Hei, feil i fordampning", avvik)
            end  
            redea = 0.0
        else
          newLayer = zeros(nodaysvector[j])
          if j == NoL 
              redea = 0.0
          else
              redea -= aktMag                    
          end
          #println(newLayer)
        end
        #println(total_layers_last)       
        Layers[1:nodaysvector[j],j] .= newLayer  # updating the actual layer
        #println(sum(Layers))
        #println(total_layers_last)    
      else                           # sumLayers == 0    
        ea = 0.0
      end                             # sumLayers ==0      
    end                   # redea > 0.0 
  end                          # for loop
   
  total_layers_current = sum(Layers)
  ea = total_layers_last - total_layers_current
   
  if ea > 0
     if round(total_layers_last - (ea + total_layers_current); digits=8) != 0
       println("Layers ut not OK")
     end
  end
#println(ea)
#println(total_layers_last)
#println(sum(Layers))                                   
  return Layers, ea
end
