#    Function UnsaturatedEXEvap
#
#     Description: Calculates soil moisture (unsaturated) 
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#
#     Parameters:
#     sm       Soil moisture content, mm
#     D        Saturation deficit (D>Z) mm
#     OUTX     Outflow form soil moisture reservoir mm/d#     
#--------------------------------------------------------------------------

function UnsaturatedExEvap(toSoil,sm,R,D)

#--------------------------------------------------------------------------------------------------------
#     Estimating retention of moistureinput to sm and runoff
#     The theory is that saturated and unsaturated sone share the same volume, and are hence the complement
#     of each other 
#--------------------------------------------------------------------------------------------------------            
      if(D > 0)
        rat = (sm + toSoil)/D
      else
         rat = 1.0 #taking into account that D can be zero or even negative (overland flow), complete saturation.
      end
      
      if (rat > R)      
        outx = (sm+toSoil) - R*D
        sm = R*D      
      else
        outx = 0.0
        sm += toSoil
      end        
return outx,sm
end        


 
