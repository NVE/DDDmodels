function TssDewpoint(Ta, RH)
#Fra Raleigh et al. Wrr 2013
# Using estimated Dew point temperaure as proxy for Tss  
#  Author Thomas Skaugen
#  Last revised : 2.12 2019 
       
	
	if Ta  > 0.0 
       b=17.625
       c=243.04 # Celscius
    else
      b=22.587
      c=273.86
    end

    Tss = c*(log(RH) +(b*Ta)/(c+Ta))/(b-log(RH)-(b*Ta)/(c+Ta))
    
return Tss

end
