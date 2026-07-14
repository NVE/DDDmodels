function SolradTransAlbedo(DN,Ta,taux,SWE,regn,thr,Timeresinsec,TSS,PS,phi,thi,Cl)

#phi = 59.94970*pi/180  #converts Lat degrees to radians degrees lat
#thi = 10.79776*pi/180  #converts Lon degrees to radians Lon
# regn is not used in this subroutine 07.07.2022
#PS is assumed in mm (used in Albedo routine to determine taux)

    deltaSWE = PS/1000.0
    thrx = (Timeresinsec == 86400) || (thr == 0) ? 24 : thr # time step for this calculation (in hours)

    ddphi = phi*180/pi #latitude decimal degrees
    ddthi = thi*180/pi #longitude decimal degrees
    gamma = (2*pi/365)*(DN-1)
    theta = 0.4092*cos((2*pi/365.25)*(DN-173)) #solar declinati0on angleday angle Liston 1995 
    theta2 = (2*pi/365.25)*(DN-80)
  
    r = 149598000     #distance from the sun [km]
    R = 6378          #Radius of earth [km]
    timezone = -4*(mod1(thi,15))*sign(thi)#ang lengdegrad ikke 
    epsilon = 0.4092 #rad(23.45)
    z_s = r*sin(theta2)*sin(epsilon)
    r_p = sqrt(r^2-z_s^2)   
    nevner = (R-z_s*sin(phi))/(r_p*cos(phi)) 
  
  if -1 < nevner < 1   
    t0 = 1440/(2*pi)*acos((R-z_s*sin(phi))/(r_p*cos(phi)))
    that = t0+5
    n =720-10*sin(4*pi*(DN-80)/365.25)+8*sin(2*pi*DN/365.25)
    sr =(n-that+timezone)/60        #sunrise
    ss = (n+that+timezone)/60       #sunset
  elseif nevner <= -1 # Midnightsun
    sr = 0.0
    ss = 24.0
  else # Dark all day
    sr = 12.0
    ss = 12.0
  end

   TTList = zeros(Float64,24)  #vector for Transmissivity
   dingom = zeros(Float64,24)  #vector for solar zenith angle
   zenangtid = zeros(Float64,24)
   SW = zeros(Float64,24)  #vector for short wave radiation per hour
    
   S0 = (118.1*10^3)/86400                  #Solar constant kJ/m2*s. Dingman 's number'
    
  for tid in 1:24
    if tid > sr && tid < ss
     tom =  -(12-(tid))         #number of hours from solar noon, which gives tom=0.0
     cosarg = 0.2618 *tom #radianer pr time
     dingom[tid] = acos(sin(phi)*sin(theta)+cos(phi)*cos(theta)*cos(cosarg))  #Dingmans zenith angle
     #zenangtid[tid] =(pi/2)-dingom[tid]  
     
     #TTList[tid] = (0.5 + 0.3*cos(dingom[tid]))*(0.20 +0.50*(1-Cl)) #Inspired by Bacellar 2008 transmissivitet and 
                                                                   #Dingman2000 p. 193 (CloudsSkyer)        
     TTList[tid] = (0.5 + 0.3*cos(dingom[tid]))*(0.355 +0.68*(1-Cl)) #Inspired by Bacellar 2008 transmissivitet and 
                                                                   #Dingman2000 (CloudsSkyer)
     #TTList[tid] = (0.55 + 0.35*cos(dingom[tid]))*(0.355 +0.68*(1-Cl)) #Inspired by Bacellar 2008 transmissivitet and 
                                                                   #Dingman2000 (CloudsSkyer) 																   
     SW[tid] = TTList[tid]*sin((pi/2)-dingom[tid])*S0*3600         # Solar constant multiplied with the seconds        
                                                                   # see eq. Liston 1995, eq. 26 (NettoSWstråling)
    end 
    if tid < sr || tid > ss  #If time is outside sunrise sunset
     TTList[tid] = 0.0 #Transmissivity, equals zero if sun below horizon
     dingom[tid] = pi*0.5 # pi/2 minus zenith angle
     #zenangtid[tid] = 0.0
     SW[tid] = 0.0 # No SW at night       
    end
  end

  interv = max(trunc(Int, Timeresinsec/3600), 1) #how many hours? we don not go below 1 hours for these calculations
  Sinn = Statistics.mean(@view(SW[(thrx-interv+1):thrx]))
  zenang = Statistics.mean(@view(dingom[(thrx-interv+1):thrx])) 
  #zenang =Statistics.mean(zenangtid[(thrx-interv+1):thrx])   
    
#Albedo: needs to include todays snow for calc. albedo Må ta med dagens snø for å regne ut albeo, eller smelter den bare vekk
SD = (1000/300)*(SWE+deltaSWE) #[m] We use 300 kg/m3 as standard density for calc albedo in UEB SDi [m]

#UEB albedo
 taux, A_UEB = AlbedoUEB(PS,SD,TSS,zenang,taux,Timeresinsec)
  #A_UEB= albUEB[2]
  #taux = albUEB[1]

A = A_UEB

#Solar radiation
S = (1-A)*Sinn*(Timeresinsec/3600) # accumulate to number of hours
    
return S, A, taux # Sinn, theta, theta2, sr, ss, Trans, zen_ang

end
