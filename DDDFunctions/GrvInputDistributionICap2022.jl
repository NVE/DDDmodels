# Function for distributing input to the differnt groundwater layers according to their individual
# level of saturation, after overland flow and evapotranspiration are taken into account.
#
# Author: Thomas Skaugen
# Revised: 30.03.2022
#################################################################################################
# ddist provides the weights to distribute outx to different layers
# The sum of ddist[1:NoL] is always 1 and the higher ddist the more water is received by the layer
################################################################################################


function GrvInputDistributionICap!(outx, NoL,ddistx, ddist, ICap)
    
if outx > 0

 @views defi = sum(ddistx[2:NoL])             # deficit in subsurface

 if defi < outx                   # the deficit is less than water input (outx)
   OF = outx-defi                  # In addition, saturation from below also gives OF
   redoutx = defi
   #println("OFBelow=",OF)
   #println("outx_OF from below=",outx)
 else
   OF = 0.
   redoutx = outx
 end
     
# outx is reduced sucessively by the layers, starting from the slowest, NoL
    
   ddist[1] = OF > 0 ? OF/outx : 0.                  # assigning distribution of outx to overland flow layer layer 1

   for j in NoL:-1:2                   # Remember NoL is the slowest layer, 1 is the fastest and already accounted for
   
    if redoutx > 0
    
       differ = ddistx[j]-redoutx
       if differ < 0                #i.e the deficit i layer j is less than input, input(outx (redoutx)) > deficit       
         ddist[j] = ddistx[j]/outx    # divide by outx informs us what frafction of outx needed for this layer
         redoutx -= ddistx[j] # must reduce  the input correspondingly       
       else               # i.e. deficit in layer j is more than input
           if j < NoL
             ddist[j] = 1.0 - sum(ddist[(j+1):NoL])-ddist[1]  #ddist[1] is already assigned
           else
             ddist[j] = 1.0 - ddist[1]                              #ddist[1] is already assigned
           end
           redoutx = 0.0
       end
    else
        ddist[2:j] .= 0.
        break
    end
   end
   #println("outx_ordinary=",outx)
else
  ddist .= 0.
end #if outx >0.0    

#if(ddist[NoL] > 0.0)
#   println("ddist=",ddist)
#   println("outx=",outx)
#end

#println(ddist*outx)
#println(sum(ddist))

end #end of func
