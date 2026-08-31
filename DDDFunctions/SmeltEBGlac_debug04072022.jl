function SmeltEB(DN,PR,PS,Ta,SWE,Timeresinsec,thr,u,taux,Pa,snittT, phi, thi)

#import Statistics
#using Distributions 
    
SWE /= 1000.0           #[m]
regn = PR + PS

Cl, RH = CloudCover(regn, Timeresinsec) # Gives a vector(Cl, Wa)

Tss = TssDewpoint(Ta,RH) # 

if (Tss > 0.0) && (SWE > 0.0) 
  Tss = 0.0
end

#ShortWave
SWrad, A, taux =SolradTransAlbedo(DN,Ta,taux,SWE,regn,thr,Timeresinsec,Tss,PS,phi,thi,Cl)# Walter and UEB albedo 
       #Net radiation appears to be OK. sr,ss Trans,zen_ang,theta, theat2 checked at 1h res and found to be OK

#LongWave, Atmospheric, Terrestrial
                 
LWradA, LWradT = LongWaveRad(regn,Ta,Timeresinsec,Tss,Cl)# uses Snow surface temperature, TSS

#Heat
SH, LE = SensibleLatHeat(PR,Ta,Timeresinsec,u,Tss,Pa,RH)    #uses Snow surface temperature, TSS

#Ground_prec_CC
GH, PH, CC = GroundPrecCC(SWE,Ta,PR,Timeresinsec,snittT,PS)    

    rhow = 1000 #[kgm^-3]
lam = 334 #kJkg^-1 latent heat of fusion    se Dingman p.190

if CC == 0
   melt = 1000* (SWrad+LWradA-LWradT+SH+LE+GH+PH+CC)/(lam*rhow) #potential melt in mm
elseif CC > 0
   melt = 0.0
end

if (snittT < 0.0) || (SWE == 0.0)
   melt = 0.0 
end

return melt, SWrad, LWradA, LWradT, SH,LE,GH,PH,CC,A,taux,Tss,Cl, RH 

end
