function LongWaveRad(regn,Ta,Timeresinsec,Tss, Cl)

#Input: Precipitation and temperature, 
#Tidsoppløningsbestemte parametere Boltzmanns konstant
#Emissivity is without dimensions
#regn is not used

#New 24.10.2016 derived from Filefjell data
#epsa = (1.0+0.0025*Ta)-(1-Cl)*(0.25*(1.0+0.0025*Ta))# Quite good, use in RH modelling?
#epsa = (1.0+0.0015*Ta)-(1-Cl)*(0.25*(1.0+0.0020*Ta))# Quite good, use in RH modelling?

#Wa = (Cl+1.667)/2.667

#Wa =(100.0-5.0*(Ta-Tss))/100.0 # M.G. Lawrence BAMS 2005
#if(Wa > 1.0) 
#   Wa = 1.0
#end

#if(Wa <= 0.0)
#   Wa = 0.1
#end

ea = 0.611*exp((17.3*Ta)/(Ta+237.3))# equals mean RH from Filefjell How to Wiki page said not to have Wa multiplied with ea 
epsa = 1.72*((ea/(Ta+273.2))^0.1428)* (1+0.22*Cl^2) #Brutsaert, P.196 in Dingman
epsa = min(epsa, 1.0)

#sigma = 4.89*10^-6                  # Boltzman's constant per day!
sigma2 = (5.67*10^-11)*Timeresinsec # kJ/m^2*s*K^4, sigma per Timeresinsec

LA = epsa*sigma2*(Ta+273.2)^4              #Atmospherical radiation. Ta in degrees Celsisus. 

ess = Tss >= 0 ? 0.9 : 0.97

LT = ess*sigma2*(Tss+273.2)^4 + (1 - ess)*LA #Terrestrial LW radiation. Emmisivity for snow, 0.97, also used for bare ground,
                                           #see Dingman p.583
return LA, LT
        
end
