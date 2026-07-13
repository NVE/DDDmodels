function DensityAge!(snowdepth,swe,T)

#compaction due to ageing

MaxDensity = 0.7
MaxChange = 0.9
kcomp = 0.5
G = 9.81
ETA0 = 3.6e6
C5 = 0.08
C6 = 0.021
rowater = 1000 #[kg/m3]
secperday = 86400

if(snowdepth > 0)
  density = swe/snowdepth
  overburden = kcomp*G*rowater*swe/1000
  viscosity = ETA0*exp(-C5*min(T, 0.0) + C6*density*1000)
  delta_depth = min((overburden/viscosity*snowdepth/1000*secperday)*1000, snowdepth-swe/MaxDensity) 

    if snowdepth < delta_depth
      snowdepth = 0.0
    else
      snowdepth = snowdepth - delta_depth
    end
end
    
end

