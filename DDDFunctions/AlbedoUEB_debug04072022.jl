function AlbedoUEB(PS,SD,TSS,zenang,taux,Timeresinsec)

#constants from Tarboton and Luce, Utah energy balance Snow Accumulation and Melt Model UEB, 1996
A_bg = 0.20         #Bare ground albedo
C_v = 0.2
C_ir = 0.5
alfa_v0 = 0.85# 0.95
alfa_ir0 = 0.65 #0.65
tau0 = 1000000
h_sd = 0.01 # 0.1# meter
b = 2 #Dickinson et al 1993

##TEST VARIABLE som er input
#PS = 0.005 #meter snøfall
#TSS = -3.2# Snøoverflate temperatur
#zenang = 1.2 #zenith vinkel i radianer
#SD = 0.02# snø i snømagasin (med PS?)
#taux = 0.7# lader fra forrige tidsskritt 
#Timeresinsec = 10800
    
# if zenang < 0.5
#    println("zenang_albedo= ", zenang)
# end      
####
    
# with newly fallen snow, tau equals zero
if PS >= 1.0#  1.0 #0.01
    taux = 0.0
end
    
r1 = exp(5000*((1/273.16)-(1/(TSS+273.16))))
r2 = min(r1^10,1)
r3 = 0.03
d_tau =  ((r1+r2+r3)/tau0)*Timeresinsec

taux += d_tau
F_age = taux/(1+taux)

alfa_vd = (1-C_v*F_age)*alfa_v0
alfa_ird = (1-C_ir*F_age)*alfa_ir0

if cos(zenang) < 0.5
   f_omg = (1/b)*(((b+1)/(1+2*b*cos(zenang)))-1)
end
if cos(zenang) > 0.5
   f_omg = 0.0
end

alfa_v = alfa_vd+0.4*f_omg*(1-alfa_vd)
alfa_ir = alfa_ird+0.4*f_omg*(1-alfa_ird)

albedo_init = 0.9*alfa_v+0.1*alfa_ir#vekter det mot visible <0.7 mum, alfa_vd dette gir max albedo på 0.92

R = (1-(SD/h_sd))*exp(-(SD/(2*h_sd)))

if SD < h_sd 
    albedo = R*A_bg+(1-R)*albedo_init
  elseif SD >= h_sd 
    albedo = albedo_init #interpolates albedo to bare ground albedo
end    

return taux, albedo                    

end
