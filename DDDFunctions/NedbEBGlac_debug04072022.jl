function NedbEB(CFR,DN,P,Ta,SWE,Timeresinsec,TX,prkorr,pskorr,thr,u,Pa,CGLAC,CX,TS,taux,snittT,phi,thi)
 
#Rain or snow?
    if(Ta > TX)
      PR = P*prkorr
      PS = 0.0
    else     
      PR = 0.0
	  PS = P*pskorr
    end
        
#merk at hvis albedoUEB brukes er Aprim egentlig forriges tidsskritts snøalder, taux
     FraSmeltEB = SmeltEB(DN,PR,PS,Ta,SWE,Timeresinsec,thr,u,taux,Pa,snittT,phi,thi)
     MW = FraSmeltEB[1]           #potential melt
     SWrad = FraSmeltEB[2]       #net SW radiation to get W/m2 multiply with 1000/Timresinsec
     LA = FraSmeltEB[3]
     LT = FraSmeltEB[4]
     SH = FraSmeltEB[5]
     LE = FraSmeltEB[6]
     GH = FraSmeltEB[7]
     PH = FraSmeltEB[8]
     CC = FraSmeltEB[9]
     A = FraSmeltEB[10]
     taux = FraSmeltEB[11]
     Tss = FraSmeltEB[12]
     Cl = FraSmeltEB[13]
     RH = FraSmeltEB[14] 

    AGlac = max(0.3, 0.5*A) # Glac Albedo is minimum 0.35 or 0.5*Snow albedo. Accounts for snowfree conditions.Range from Hock,2005
    rhow = 1000 #[kgm^-3]
    lam = 334 #kJkg^-1 latent heat of fusion    se Dingman p.190

    #CC == 0.0
    #TS = 0.0                # Temperature threshold for melt (glaciers)
    Δ = Ta - TS
	MWCX = CX * Δ		# potential snowmelt using degreeday factor
    
    if Δ < 0              #melt (MW) returns as negative    
      MW *= CFR            #refreezes the melt, returns as negative and as mm refreezed water
      MWGLAC = 0.0  
      MWGLAC1 = 0.0    
    else
      MWGLAC = CGLAC * Δ  #Glacial melt
      MWGLAC1 = 1000.0*(SWrad*((1-AGlac)/(1-A))+LA-LT+SH+LE+GH+PH)/(lam*rhow) #potential melt in mm Glacial melt. 
	  #Note that SWrad already has been corrected for Albedo. We thus have to correct it back by 1/(1-A). We multiply with 1000 to go from [m] to [mm] 
    end      

return PS,PR,MW,MWGLAC,SWrad, LA,LT,SH,LE,GH,PH,CC,A,taux,snittT, Tss,Cl, RH, MWGLAC1, MWCX
end
