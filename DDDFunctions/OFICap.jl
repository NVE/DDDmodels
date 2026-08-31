# Function for calculating the Overland Flow (OF only) due to limited infiltration capacity. Distribution of water to 
# differnt groundwater layers according to their individual level of saturation and additional overland flow that may
# come from saturation from below is calculayted later and comes after evaotranpiration
# and updating the soilmoisture.  
#
# Author: Thomas Skaugen
# Revised: 30.03.2022
#################################################################################################
                                       
# ddist provides the weights to distribute outx to different layers
# The sum of ddist[1:NoL] is always 1 and the higher ddist the more water is received by the layer
# This subroutine takes into account Infiltration capacity ICap

################################################################################################


function OFICap(outx, ICap)
 
if ICap < outx             # OF due to exceedance of ICap, In addition we may have OF from below
  OF = outx-ICap            # what cannot be infiltrated goes to overland flow. 
  nyoutx = ICap 
else
  OF = 0.0
  nyoutx = outx
end
    
return OF, nyoutx  

end #end of func

