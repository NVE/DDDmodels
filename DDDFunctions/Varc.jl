#Function varc
#This function updates and returns dynamical values of nu and alfa, of the snowroutine.
# The correlation coefficient is evaluated as the sum of the correlation N*N, 2*N*u and u*u matrices
# Should be a decreasing value with n, large n small c
# All variables must be conditional (i.e. calculated on snowcovered area)
# Note the distinction between conditional and unconditional units of melt
#######################################################
#This function has input:
#alfanull (ac) and drange      : parameters of relationship mean standarddeviation of conditional precip.
#sdnux,sdalphax,nudyn,alfadyn   :initial distr parameters, current distr parameters
#ppa,redsca,nnn,n,u   :SCA,relative reduction of SCA, previously accumulated units,current accu units and melt units
#########################################################

# sdnux og sdalphax = ny0 og alfa0

function Varc(ppa,ac,drange,nudyn,alphadyn,sdnux,sdalphax,nnn,u,n,redsca)
    
 #ppa = 1.0
 #nudyn = 0.118124339507407
 #alphadyn = 1.18124339507407
 #sdnux = 4.077 #Nynull
 #sdalphax = 40.77 #alfanull
 #nnn = 10
 #drange = 113 #range of correlation function
 #n = 2
 #u = 0
 #n = 1
 #redsca = 0 

nnn = round(nnn) #previously accumulated units
n = round(n)     # new accumulated units
u = round(u)     # melting units
krysscorr = 0.0

#variance in: remember that nudyn inkludes nnn
vardyn = nudyn/alphadyn^2

#correlation =1 for index (1,1)
function corrvec(nnn,drange) # Correlation function
  corrv = exp(-nnn/drange)
  return corrv
end

# ACCUMULATION ####################
if(n > 0)
#if nnn = 0 then the totalvariance is equal to autocovariance calculated for n (fresh start)
#vardyn and kryssvar equals 0 and the mean value equals n*unitmean (unitmean = sdnux/sdalphax).
   if nnn == 0
     vartot =  n*(sdnux/sdalphax^2)*(1 + (n-1)*corrvec(n,drange))
     meantot = n*(sdnux/sdalphax)
   elseif nnn > 0 # must scale the corrvec
     varoldcov = vardyn + n*(sdnux/sdalphax^2)*(1 + (n-1)*corrvec(n,drange)) # 08.01.2016  antar uavhengighet mellom gammelt og nytt, hensikten er å øke variansen
     # calculating the n matrix
     varnewcov = n*(sdnux/sdalphax^2)*(1+(n-1)*corrvec(n,drange)) #denne er OK
     #total variance: (nnn+n)*ppa +n*(1-ppa)
     vartot = (varoldcov*ppa^2)+(varnewcov*(1-ppa)^2)
     #new condtional mean for 100%coverage
     meantot = ((ppa*(nnn+n))+((1-ppa)*n))*sdnux/sdalphax #For alle dekningsgrader
   end
end # for n >0

# ABLATION ########################
if(u > 0.0)
   faktor = (((nudyn*sdalphax^2)/(nnn*sdnux*alphadyn^2))+1+(nnn-1)*corrvec(nnn,drange))/(2*nnn)
# We must ensure that u is conditional since u enters the routine as an unconditional value (melts over the entire area)
    meanubet = (nnn-u)*(sdnux/sdalphax)
    meantot = meanubet/(1-redsca)            # Is the CONDITIONAL VALUE
    betu = round((1-redsca)*nnn-(nnn-u))     # according to the Skaugen and Randen, 2013
    if(betu <= 0)                   # taking into account that new conditional mean is higher than previous conditional mean, can't have that
        betu = 0.0 
        autovar = 0.0
        kryssvar = 0.0
    else
        autovar = (sdnux/sdalphax^2)*(betu+betu*(betu-1)*corrvec(betu,drange))
        corrun = faktor*(betu/nnn)
        kryssvar = (sdnux/sdalphax^2)*(2*nnn*betu*corrun)
    end
    vartot = vardyn+autovar-kryssvar #NB Kryssvar should be reduced. Vartot is the CONDITIONAL VALUE
end #for u > 0.0 END ABLATION
        
#Calculating the updated parmeters for the spatial distribution
#NB note that nudyn and alphadyn includes the nnn
      if(meantot==0)# the reservoir (nnn) and melt (u/betu) match eaxtly       
        nudyn = sdnux      #This feature gives nnn as one, i.e. never zero, but that is taken care of in other subroutines 
        alphadyn = sdalphax
      else
       nudyn = (meantot^2/vartot)
       alphadyn = meantot/vartot  # Conditional distribution
      end
      if isnan(vartot)
       println("vartot=NaN (ny_varc) nudyn= ",nudyn," alphadyn= ",alphadyn, " meantot= ", meantot," vartot= ",vartot," nnn= ",nnn," n= ",n," u= ",u)
       sleep(5) #stops for 5 seconds
      end
      nnn = round((nudyn/alphadyn)/(sdnux/sdalphax))  #ensures correspondence between nnn, nydyn and alfadyn
      alphadyn = nudyn/((sdnux/sdalphax)*nnn)

return nudyn, alphadyn, nnn
end
