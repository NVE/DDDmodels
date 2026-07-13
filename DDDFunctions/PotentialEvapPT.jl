#    Function Potential_evap_PM(isoilx,sncovx,smx,outx,ceax,tempx,lp,ep)
#
#C------------------------------------------------------------------------
#
#     Description: Estimates the Potential Evapotranspiration using Priestly-Taylor (Penman-Monteith model start is also here (Dingman, p.299) 
#     Author: Thomas Skaugen
#     Revised: 02.03.2018
#
#     Parameters:
#     
#--------------------------------------------------------------------------

function PotentialEvapPT(eatemp,SWrad,LA,LT,Pa)
 
ca = 1.01                                                          # heat capacity air kJ/kgK, p197 in Dingmann
LaV = 2470                                                         # [kJm^2] latent varme fra fordampning
rhow = 1000                                                        # density of water [kg/m3]
gam = ca*Pa/(0.622*LaV)                                            # pyscrometric constant (p.275 in Dingman)
delta = (2508.3/(eatemp+237.3)^2)*exp((17.3*eatemp)/(eatemp+237.3)) 
alfa = 1.26                            # constant in Priestly Taylor for humid areas, Weiss and Menzel (2008) Adv. Geosci.
 
ET = alfa*(delta/(delta+gam))*(SWrad+LA-LT)*(1000/(LaV*rhow))      # [mm/dt] Priestly-Taylor potential evapotranspiration

return max(ET, 0.)
end
