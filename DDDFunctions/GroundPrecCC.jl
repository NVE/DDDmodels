function GroundPrecCC(SWE,Ta,Rp,Timeresinsec,snittT,PS)

# CC is Cold content as a fucntion of SWE and snowpack temperature (Dingman)
# GH Ground Heat
# PH Precipitation Heat
# Rp is precipitation (RP) comes in as mm but is used as[m], see p. 199 Dingmann

# Rp = 2.0
# Ta = -5.0
# Timeresinsec = 10800   
# snittT = -2.3  
# PS = 2.0
# SWE = 2.00
    
 GH = 173.0/86400.0         #kj/secund
 GH *= Timeresinsec   # desired resolution
    
#precipitation heat
rhow = 1000         #kg/m3
CW = 4.19           #kJ/kg K Heatcapacity of water, see p.199 Dingmann
if(Rp > 0.0)
 PH = rhow*CW*(Rp/1000)*Ta
else
 PH = 0.0
end

#Energy stored in snowpack, Cold content
cc = 2.102                    #kJ/kgK  varmekapasitet is, se Dingman p. 189
                              #snittT is the mean temperature of snowpack
CC = rhow*cc*(SWE+PS/1000)*(snittT) #is negativ or zero, SWE is in [m] 

return GH, PH, CC

end
