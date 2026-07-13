#Function SnowGamma
#This function has input:
#MW              :potenial melt/refreeze(negative) mm
#PR              :liquid precipitation
#PS              :solid precipitation
#ny0             :ny parameter in precip gamma distribution (shape)
#alfa0           :alfa parameter in precip gamma distribution (scale)
#cc              :temporal correlasjon between events
#pro             :max % liquid water in snow
#sca             :Snow covered area (for every elvation zone)array(idim))
#idim            : #no elevation zone 

#Recall that SPD, WCD are unconditional values, must be multiplied with PPA for arealvalues!!!
#Output:
#ISOIL            :runoff to subsurface- or overland flow (unconditional value, UCV)
#SPD              :swe, conditional value (CV)
#WCD              :fritt vann i snøpakke CV
#nsno             :number of events
#sca              :updated snow coverage 
#hc and ac        :describes the relations between areal mean and standard deviation of SWE (see Skaugen and Weltzien, 2016)

#using Distributions

function SnowGamma(PRX,PSX,MWX,scax,spdx,wcdx,prox,nsnox,alfax,nyx,alfa0x,ny0x,ac,hc)
 
 if(scax*spdx < 0.1 && PSX == 0.0)

    scax =0.0
    ppa = 0.0
    spt = 0.0
    isoil = PRX
    nsnox = 0.0
    nn = 0.0
    alfadyn = alfa0x
    nydyn = ny0x
    wct = 0.0  
    xmw = 0.0
    
 else                                 #scax*spdx > 0.1 eller PSX == 0.0   
    redsca = 0.0
    na = 0.0              # number of units for accumulation event
    u = 0.0               # number of units for melt event
    totps = PSX           # rpecipitation as snow UCV
    spt = spdx            # SWE CV kept as CV 
    wct = wcdx            # liquid water in snow CV
    xmw = MWX             # UCV melt for areas with snow. Must multiply with ppa for actual runoff CV                                                             
    ppa = scax            # actual coverage for catchment, not relative to coverage!
    ppaold = ppa
    nn = round(nsnox)     # number of events spt/(ny0/alfa)
                          # CV conditional on snow. UCV, unconditional on snow (areal values)       
    alfadyn =  alfax      # Transferred variable   
    nydyn = nyx           # Transferred variable   
                    
    if(nn > 0.0)                            
      alfadyn = alfadyn
      nydyn = round(nn) * nydyn # nydyn exits the routine as per unit, enters as the product of nn* nydyn
                                # her blir den jekket opp
    else
      nydyn = ny0x
      alfadyn = alfa0x
    end
        
#Recall that spt is the ice part of SWE, whereas wct er liquid water in SWE. Thes are treated separately
    totsn = (spt+wct)            # SWE and liquid water in snowpack, CV before actual avant associated with scax
    
    #Updating the water content in snowpack 
    if (xmw < 0.0)               # Refreezing liquid water in snow, the sum of spt+wct and SCA will not change
      if -xmw < wct
           wct = wct + xmw       # update wct, reduces water content in snow, NB xmv is negative
           spt = spt - xmw       # update, increases spt, reduces wct 
      end
      if -xmw >= wct
        spt = spt + wct          # update spt, increases spt and wct is set to 0
        wct = 0.0
      end
    end          
# spt is conditional value. For areal values needs to multiply with ppa
      
# Proceeds  with gamma distributed accumulation and melt 
# ACCUMULATION
     if(totps >= 0.1)                        # totps equals snø (plus poss. refreezing) for event
        na = round(totps/(ny0x/alfa0x),digits = 0)        #number of units in accumulation event 
            
        if(na == 0.0) 
            na = 1.0                         # round up if there is actually some snow
        end
     
        if(nn == 0.0)                        # no inital snow, i.e number of initial units equals zero
        
            ppa = 1.0                        # default setting after a snowfall event
            ppaold = ppa
          
            if (na == 1.0)                   # accumulation events equals zero                                         
              alfadyn = alfa0x
              nydyn = ny0x        
            else                              # na > 1.0
             fraVarc = Varc(ppa,ac,hc,nydyn,alfadyn,ny0x,alfa0x,nn,u,na,redsca) 
             nydyn = fraVarc[1]
             alfadyn = fraVarc[2]
             nnVarc = fraVarc[3]
             if(alfadyn < 0.0)
              println("jeg er i snogamma 4.1 ppa=",ppa,ac,hc,nydyn,alfadyn,ny0x,alfa0x,nn,u,na, redsca, PSX, PRX)
             end
             if(nydyn/alfadyn/na > 0.10000000001)
               println("diff mellom nnVarc og nn",nnVarc,nn, nnVarc-nn)
               println("nn er 0 og na er >0, akkumulasjon",nydyn/alfadyn/n)
               sleep(2)
             end
            end
            spt = round(na,digits=0)*(ny0x/alfa0x)# Cond mean = UCond mean(ppa = 1.0) and is updatetd due to snowfall (na >0, 1.0)                                
        else #nn and na > 0.0               
          fraVarc = Varc(ppa,ac,hc,nydyn,alfadyn,ny0x,alfa0x,nn,u,na,redsca) 
          #return Varc nudyn, alphadyn, nnn
          nydyn = fraVarc[1]
          alfadyn = fraVarc[2]
          nnVarc = fraVarc[3]
          nn = round(nn*ppa)+ round(na) # # adjusting the mean for initial coverage < 1.0
          spt = nn*(ny0x/alfa0x)# Conditional mean equals unconditional mean (ppa = 1.0) and is updatett due to snowfall (na >0, 1.0)  
          ppa = 1.0
          ppaold = ppa
         if(nydyn/alfadyn/nn >0.10000000001)
           println("diff mellom nnVarc og nn",nnVarc,nn, nnVarc-nn)
           println("nn er >0 og na er >0, akkumulasjon",nydyn/alfadyn/nn)
           sleep(2)
         end
        end # slutt else nn > 0.0
       
     nn = round(spt*alfa0x/ny0x, digits=0) #oppdaterer etter akkumulasjon  
     nsnox = nn  #oppdaterer etter akkumulasjon
                                    
     end # if sentence for  totps > 0.1
     #for accumulation.  New coverage is set to 1.0  (wise?)
     
# ABLATION xmw
     if (xmw >= 0.1) # xmw reduserer SWE (is) bidrar til wct og avrenning 
      u = Int(round(xmw/(ny0x/alfa0x), digits=0))  # UCV potentially over the entire catchment, but only for SCA                                   
                                             # uu is potential melt over SCA
                                             # nsno() becomes the new n after melting event CV
                                             # nn is original n before melt event
                                             # snømag == 0.0 if snømag < 0.2 of event mean       
                                             # u == 1 if < 0.5
   
     if((nn - u) < 2)
       nn = 0.0
       alfadyn = alfa0x
       nydyn = ny0x 
       spt = 0.0
       wct = 0.0
       ppa = 0.0
     
     else # nn-u is greater than 2 units
         
         function corrvec(nnn,drange) # Correlation function
          corrv = exp(-nnn/drange)
          return(corrv)
         end

        vars = (ny0x/alfa0x^2)*(u+u*(u-1)*corrvec(u,hc))   # varianca of melt
        ms = u*(ny0x/alfa0x)                               # mean of melt
        nys = ms^2/vars                                    # melt ny, shape parameter
        alfas = ms/vars                                    # melt alfa, scale paramater
        alfaa = alfadyn                                    # accumulation alfa
        nya = nydyn                                        # accumulation ny teste med nn*nydyn
   
        antall = Int(round(nya/alfaa) +1)                  #  How far we need to go (in mm) for the pfd.
        diff = zeros(Float64,antall)  
        xkrit = zeros(Int,2)
                                       
         smelt = Gamma(nys,1/alfas)                                        
         acc = Gamma(nya,1/alfaa) # Note that R uses alfa whereas Julia uses 1/alfa
         
       for i in 1:antall
        s = pdf(smelt,i) # PDF; wants the crossing point bewtween dist (s=melt)
        a = pdf(acc,i)   # PDF; wants the crossing point bewtween dist (a=accumulation)
        diff[i] = a - s
       end        
                                
        #find crossing
        krysspos = findall(diff .> 0)
        if(length(krysspos) >= 1)        
          xkrit[2] = Int(krysspos[1])
        else
          xkrit[2] = 1
        end
                        
        if(xkrit[2]==0)
          println("No crossing point!") 
          sleep(2)                     #pause for 2 seconds
        end

#start exact method 
        pa = cdf(acc,xkrit[2])    # CDF; wants the area
        ps = cdf(smelt, xkrit[2]) # CDF; wants the area
        redsca = round((pa + (1-ps)),digits = 4)

# end exact method
        if (redsca <= 0.0)
           redsca = 0.0001
        end
        if (redsca >=1.0)
          redsca = 0.9999
        end
                       
        newsca = 1-redsca                          # reduction from previous melting; actual coverage is ppa <- ppa *newsca  
        spt = (1/(1-redsca))*((nn-u)*(ny0x/alfa0x))# New CV mean; mean for snowcovered area
        ppa = ppaold*newsca                        # updating actual SCA. NOTE that ppa is the coverage for the elevation zone
       
        if(spt <= 0)
             ppa = 0.0 # 
        end
        wcmax = spt*prox # current max level of free water in snow
                                                                
        if((xmw + PRX) >= wcmax) 
             wct = wcmax 
        end
                                                                
        if((xmw + PRX) <  wcmax) 
             wct = xmw +PRX                                           
        end
       # NOTE it might be a problem that the potential melting XMW can never be actual, even for the areas with snow
       # because the areas left snow free are areas with snow less than XMW      
       
       if(ppa < 0.02)                          #Initializing
         nn = 0.0
         alfadyn = alfa0x
         nydyn = ny0x
         spt = 0.0
         wct = 0.0
         ppa = 0.0       
       else  
        newnn = spt/(ny0x/alfa0x)             # accumulation updated after ablation
        nn = nsnox                            # old accumulation 
        fraVarc = Varc(ppa,ac,hc,nydyn,alfadyn,ny0x,alfa0x,nn,u,na,redsca) 
        nydyn = fraVarc[1]
        alfadyn = fraVarc[2]            
        nsnox = newnn                         # update units etafter melting
        nn = newnn
       end                                     # if ppa > 0.05
             
       end     #if (nsnox-uu < 0.1)
     end       #if(xmw >0.1)
          
# ABLATION is finished

sptgml = spt     # this TS's hydrology beregnes ut i fra ikke oppdatet spt. The effect of updating in next TS
ppagml = ppa       
nnn = nn

isoil = (scax*totsn)+PSX+PRX -(ppagml*(sptgml+wct))# This is moisture according to model and input NOT due to updating. 
                                                   # The next  time step will take updating into account
                                                   # isoil <- (scax*totsn)+PSX+PRX-(ppa*(spt+wct))
                                                   # pr and ps are already areal values
                                                   # spt and wct must be multiplied with new ppa 
                                                   # meaning, not updated, but this timesteps meteorologi is taken into account
if(isoil < 0.0001)
    isoil = 0.0
end

end #if setning som er if(sca*spd < 1.0 && PSX == 0.0)

if(nn <= 0.0) # 
    nydyn = ny0x
else
   nydyn = nydyn/nn # parameter exits as not for units
end
#print(paste("ut av snogamma 2:nydyn=",nydyn, "alfadyn:",alfadyn,"nn=", nn))
#variabler som skal overføres
spdx = spt
wcdx = wct
scax = ppa
nsnox = nn
    
return isoil, spdx, wcdx, scax,nsnox,alfadyn, nydyn
end



