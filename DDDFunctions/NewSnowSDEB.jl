#calls NewSnowDensityEB

function NewSnowSDEB(PRX,PSX,T,TX,xmw,swe,swex,snowdepth)
MaxDensity = 0.7
#snowdepth is current snowdepth prior to updating

if(T > TX)
  newsnow = 0.0    
else
  newsnow = PSX
end
    
if(swex > 0 && swe > 0)
    density = swex/snowdepth 
    if(xmw > 0) 
       snowdepth = snowdepth - xmw/density
    end
    if(snowdepth < 0)
       snowdepth = 0.0
    end
else   
   density = 0.1
   snowdepth = 0.0
end 
    
density_new = NewSnowDensityEB(T)
    
if(swe > 0)
    if(snowdepth==0)
      delta_depth = 0
    else
     delta_depth = 25.4*(newsnow/25.4)*(snowdepth/25.4)/((swe-newsnow)/25.4)* (snowdepth/25.4/10.0)^0.35
                                                          #From Bras reduction in depth due to weigth of newsnow
     if(delta_depth > snowdepth-swe/MaxDensity)
       delta_depth = snowdepth-swe/MaxDensity
     end
    end                                                                 
    snowdepth = snowdepth + (newsnow/density_new) - delta_depth # updatet snowdepth due to newsnow
end
 
DensityAge!(snowdepth,swe,T)  #snowdepth due to ageing 

return snowdepth, density
end
