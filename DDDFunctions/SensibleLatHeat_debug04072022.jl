function SensibleLatHeat(regn,Ta,Timeresinsec,u,Tss,Pa, RH)
#regn [PR) is not used        

ca = 1.01            # heat capacity of air kJ/kgK, p197 i Dingmann
rhoa = 1.29          # kg/m3
zu = 10#             # Height [m] of wind measurements
#zt = 2               # Height[m] of temperature measurements
zm = 0.0015#         # roughness parameter for  snow Also called z0
d = 0                # Zeroplane displacement for snow
#zh = 0.0002          # 
k = 0.4              # von Karmans konstant
common = k^2/(log((zu-d)/zm))^2
H = ca*rhoa*common*u*(Ta-Tss) * Timeresinsec # positive if Tss < Ta (same for latent Heat) 

Wa =(100.0-5.0*(Ta-Tss))/100.0 # M.G. Lawrence BAMS 2005
if Wa > 1
   Wa = 1.0
elseif Wa <= 0
   Wa = 0.1
end

#Latent heat                     
ea = 0.611*exp((17.3*Ta)/(Ta+237.3))*Wa    # always positive. D7, D10 page 587 Dingman actual  vapor pressure, (because of RH) 
es =  0.611*exp((17.3*Tss)/(Tss+237.3))    # always positive
LaV = 2470                                 # kJm^2, latent heat from evaporation
LaF = 334                                  # kJ,    latent heat from fusion

if Tss < 0
 LE = (LaV+LaF)*0.622*(rhoa/Pa)*common*u*(ea-es)
else
 LE = (LaV)*0.622*(rhoa/Pa)*common*u*(ea-es)
end

LE *= Timeresinsec

#if((ea-es) > -100.0)
#  println("Tss=",Tss)
#  println("Ta=",Ta)
#  println("RH=",RH)
#  println("LE=",LE)
#  println("PA=",Pa)
#  println("ea=",ea)
#  println("es=",es)
#  println("ea-es=",ea-es)
#end  
 
return H, LE
end
