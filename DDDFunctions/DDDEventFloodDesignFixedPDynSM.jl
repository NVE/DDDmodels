# DDD-Flood stochastic simulation of floods in urban areas
#
#----------------------------------------------------------------------------------
#     Script name:  DDDurbanFloodDesign.jl
#     Author:     Thomas Skaugen and Ina Cecilie Storteig
#     Revised: 04.07.2023
#              
# Purpose   : DDD-Flood (Distance Distribution Dynamics) hydrological model simulates
#               1) saturated soil water flow, states of saturated (2-D)) and unsaturates zones
# 
#               3) urban runoff for a given catchment.
#               4) 
#               5) Uses EB for snowmelt and evapotranspiration
#               6) estimates overland flow rivernetwork 
#               7) we will erform simulations for three reservoirs: 1)IP impermeable ares 2) P, permeable areas and 3) water in conduits (we need to know capacity)
#
# Model input: 1) catchment physical/geographical description, distance distributions for landscape features and river network
#              2) weather data, temperature,precipitation 
#              
#
# Model output: simulated 1)timeseries of runoff and subsurface states
#                         2)temporal distribution of extreme precipitation
#                         3)flood values
#                         4) ++
#                         
#----------------------------------------------------------------------------------

function DDDEventFloodDesign(tprm, prm ,utfile, utfile2, NumSim, ext_precip_max)

    
Lty = 2                      # number of landscape types in urban areas, for now Lty = 2, where Permeable (P) = 1 and ImPermeable (IP) =2, does not include wetlands

#-------------------------------------preprocess start------------------------------------------
#-----------------------------------------------------------------------------------------------

#####################VARIABLES FOR 10 ELEVATION ZONES and for Parts ###########################################
isoil = zeros(Lty)# precipitation and snowmelt from the elevation zones
PR = zeros(Lty)
PS = zeros(Lty)

#################################Parameters as vectors#############################################
MAD = zeros(Lty)          # really just one value, total for the entire catchment
Ltyfrac = zeros(Lty) # areal fraction landscape type permeable area, includes forest and parks,fields, impermeable are built up areas and streets
Ltymax = zeros(Lty) 
Ltymid = zeros(Lty) 
Ltyz = zeros(Lty)   
   
#################################States as vectors#################################################
M = zeros(Lty)
GshRes = zeros(Lty) # shape parameters og subsurface reservoir [resevoir in mm]
GscRes = zeros(Lty) # scale parameters og subsurface reservoir [reservoir in mm]
ICap = zeros(Lty)       # Infiltration capacity foor P[1] and IP[2] in [mm/s]
ICapdummy = zeros(Lty)
lyrs = zeros(Lty)
subsurface = zeros(Lty)
area = zeros(Lty)     # includes the area for wetlands
totdef = zeros(Lty)
Def = zeros(Lty)
PotEvap = zeros(Lty)
tosoil = zeros(Lty)
outx = zeros(Lty)
gmloutx = zeros(Lty)
OF = zeros(Lty)
########################Reading parameters#########################################################

# Misc.
Timeresinsec = prm.val[3]       # in seconds
MAD[1:Lty] .= prm.val[4]     # mean annual discharge
totarea = prm.val[5]            # Total area of catchment in m2
NoL = Int(prm.val[6])                # Number of layers in subsurface including OF level
R = prm.val[7]                  # soil moisture content/100 for field capacity of soils 
    
#Distance distributions Lty
Ltyfrac[1:(Lty)] .= prm.val[8:9]  # areal fraction permeable area:forest and parks,fields, impermeable area:roads and built up areas, wetlands 
Ltymax[1:(Lty)] .= prm.val[10:11]  # Maximum distances Landscapetypes
Ltymid[1:(Lty)] .= prm.val[12:13] 
Ltyz[1:(Lty)] .= prm.val[14:15]     # frac zero dist permable, impermeable area, 

#Distance distribution conduit network
midFl = prm.val[16]        # mean of distance in river/conduits network
stdFl = prm.val[17]        # standard deviation of distances in river/conduits network
maxFl = prm.val[18]        # maximum distance in river/conduits network

#Celerities, subsurface and conduits
OFVP  = tprm[1]     #prm.val[19]       # Overland flow velocity P
OFVIP = prm.val[20]       # Overland flow velocity IP
GshInt = prm.val[21]        # shapeparameter saturation Capital lambda
GscInt = tprm[2] #prm.val[22]        # scaleparameter saturation Capital lambda
Gshape, Gscale = Big2SmallLambda(GshInt, GscInt) # Coverting integrated celerity to layer celerity    

rv = prm.val[23]    # celerity of water in river/conduits network, we'll probably want a celeritydistribution associated with levels here, have data!

persons = tprm[3] # prm.val[24]

#Precipitation 
Precresinsec = prm.val[25]  # Duration of extreme value precipitation
GPloc = prm.val[26] # location parameters of Extremevalue precipitation (Generalised Pareto)
GPsc = prm.val[27] # Scale parameters of Extremevalue precipitation (Generalised Pareto)
GPsh = tprm[4] # prm.val[28] # Shape parameters of Extremevalue precipitation (Generalised Pareto)
SMShConst =prm.val[29] # We estimate the Shape parameter of gamma distributed SM as a function of precip through linear regression, the constant 
SMShSlope =prm.val[30] # We estimate the Shape parameter of gamma distributed SM as a function of precipt hrough linear regression, the slope 
SMScConst = prm.val[31] # We estimate the Scale parameter of gamma distributed SM as a function of precip through linear regression, the constant 
SMScSlope = prm.val[32] # We estimate the Scale parameter of gamma distributed SM as a function of precip through linear regression, the slope
    
ICap[1] = tprm[5]/ 3600.0  # prm.val[31]/3600.0 # Infiltration capacity for Permeable[1] surfaces. Input from paramfile mm/hour. Previously (1000)*(GshInt*GscInt*Ltymid[1]/Timeresinsec)
ICap[2] = prm.val[34]/3600.0# # Infiltration capacity for ImPermeable[2]  previously as ICap[1]*0.002  
GshRes[1:(Lty)] .= prm.val[35]
GscRes[1:(Lty)] .= prm.val[36]
#Done reading parameters
gtcel = 0.99 # Quantile decribing full subsurface 

days = Int(Precresinsec/Timeresinsec)     # Length of Time series simulated
days2 = 10     #2*days                           # Length of we want to see the recession . 10 fits the current outputformat

#vectors and matrixes
ddist = zeros(Lty,NoL)          # distribution of water in layers
ddistx = zeros(Lty,NoL)         # water present in layers
ddistxtemp = zeros(Lty,NoL)     # temporay of the above
Magkap = zeros(Lty,NoL)
aktMag = zeros(Lty,NoL)         # current water in Layers
k = zeros(Lty,NoL)              # subsurface celerities (including that of overland flow) 

qberegn = zeros(days2)            # runoff 
lyrberegn = zeros(days2)          # for monitorng the SS
SSberegn = zeros(days2)   # SS contribution of runoff
Pberegn = zeros(days2)    # Precip contribution of runoff
Initberegn = zeros(days2) # runoff solely due to initial conditions

wgtprec = zeros(Float64,days)
ext_precipvec = zeros(Float64,days2) # extra days to see recession (no precip for these days)


# Areas and weights
area = Ltyfrac.* totarea                   #areas of differnt Lty, a vector
    
#Infiltration capacity in mm/s. Must mutiply with Timeresinsec
for Lst in 1:Lty 
    ICap[Lst] = ICap[Lst] .* Timeresinsec
end
#println("ICap ",ICap[1]," ",ICap[2]) 
                                       
#UH for Conduits, normally distributed 
#number of time units in river routing, rounded up. If set equal to one day, time is actually less
UHriver, noDT, nodaysRiv = SingleNormalUH(rv,Timeresinsec,midFl,stdFl,maxFl) #This subroutine extracts a UH normally distributed.

#println("noDT =", noDT)

#Estimating reseroir capacity and  (gamma) parameters of subsurface state (in mm) ONhold for the moment, 04.06.2025
#for Lst in 1: Lty
#  Magkap[Lst,1:NoL], M[Lst], GshRsrv[Lst], GscRsrv[Lst] = LayerEstimationDesign(GshInt,GscInt,Timeresinsec,Ltymax[Lst],Ltymid[Lst], 
#        Ltyfrac[Lst]*MAD[Lst], Ltyfrac[Lst]*totarea, NoL, gtcel) 
#end
#We try a new formulation where Sm is specified from the simulated SM of the DDDv2 model 04.06.2025
for Lst in 1: Lty
  Magkap[Lst,1:NoL], M[Lst] = LayerEstimationDesignEmpSM(GshRes[Lst],GscRes[Lst], NoL, gtcel) 
end
#println("Magkap ", Magkap[1,1]," ",Magkap[1,2] )
#println("M fra Mainruitne ", M[1])
#println("Shape Sm= ",GshRsrv[1])
#println("Scale Sm = ", GscRsrv[1]) 
#Note that MAD is for the entire catchment. It is fractioned in the above subroutine
numwgt = Timeresinsec/Precresinsec # How does the optimum temporal runoff resolution "chop" up the duration of extreme rainfall 

precip_rounds = Int(trunc(ext_precip_max*2))  # Number of different precip values we are going to simulate for, up to  ext_precip_max by 0.5 mm
#Output matrices ##################################################################

    simresult = zeros(days2,24)        # matrix which into the results are written
    simresult2 = zeros(NumSim*precip_rounds,24)        # matrix which into the results are written
    simresult3 = zeros(NumSim*precip_rounds,(5*days2)) # matrix which into the results are written, 50 because the initial precip is evetually multiplied with 50 
 
###################################################################################

for l in 1: precip_rounds
   ext_precip = 1 + (l-1)*0.5 # Fixed (extreme precipitation value) for critical duration to be dragged through all possible combinations of conditional soilmoisture and precip sequence

   println("Precip= ", ext_precip)
 # ext_precip = 50.0  #Now supplied by the runscript. fixed extreme value: Risvollan X years, 15 minutes duration. IVF KSS
  
#############################################################################################
#Stochastic simulations start here
#############################################################################################
   for j in 1:NumSim
    
##Drawing stochastic varaibles from groundwater state and precipitation    
#grvstate = rand(Gamma(GshRsrv[1], GscRsrv[1]),1) # Draw groundwater state, ordinary distributions (not extreme)
#Groundwater layers; 2dim levels, 1 fastest, NoL slowest#
    
#Initializing (everything is zero) everything for current round 
qsimutx = zeros(noDT)      # common vector for runoff from all landscapetypes [m3/s]
qsimutxP = zeros(noDT)     # vector for runoff permeabel i m3/s
qsimutxIP = zeros(noDT)    # vector for runoff impermeable i m3/s
qsimutxOF = zeros(noDT)    # vector for Overland Flow
qsimutxOFP = zeros(noDT)    # vector for Overland Flow P
qsimutxOFIP = zeros(noDT)    # vector for Overland Flow IP
qsimutxSSFP = zeros(noDT)    # vector for Overland Flow P
qsimutxSSFIP = zeros(noDT)    # vector for Overland Flow IP 
qsimutxSS_contribP = zeros(noDT)
qsimutxSS_contribIP = zeros(noDT)
qsimutxP_contribP = zeros(noDT)
qsimutxP_contribIP = zeros(noDT)
qsimutxSS_contrib = zeros(noDT) # vector for SS contribution 
qsimutxP_contrib = zeros(noDT)   # vector for precip contribution 
SS_contribP = zeros(noDT)
SS_contribIP = zeros(noDT)
P_contribP = zeros(noDT)
P_contribIP = zeros(noDT)
qsimutxInit = zeros(noDT)

QRivx = zeros(noDT)        # all qsimutvectors are stored UHriver timesteps back. 
QRivxP = zeros(noDT)
QRivxIP = zeros(noDT)
QRivxOF = zeros(noDT)
QRivxOFP = zeros(noDT)
QRivxOFIP = zeros(noDT)
QRivxSSF_P = zeros(noDT)
QRivxSSF_IP = zeros(noDT)
QRivxSS_contrib = zeros(noDT)
QRivxPrecip_contrib = zeros(noDT)
QRivxInit = zeros(noDT)

#Random element 1, Sunbsurface state
# we have that soilmoisture (SM) distribution is gamma distributed
# and a function of precipitation
# Distribution parameters SMSh and SMSc are regressed with precipitation (ext_precip)
SMSh = SMShConst + SMShSlope*ext_precip # Shape parameter 
SMSc = SMScConst + SMScSlope*ext_precip # Scale parameter

 grvstate = rand(Gamma(SMSh ,SMSc),1) # Draw groundwater state [mm], ordinary distributions (not extreme),Gamma(shape ,scale) 
 #println("grvstate, ext_precip, M[1] =", grvstate, " ", ext_precip, " ",M[1])
####################################################################################
# This is for variabel celerity conditioned by saturation 
#satur = ext_precip + grvstate[1]
satur = grvstate[1]               #This is tthe drawn groundwater state from which dist is conditioned on precip 
g = Gamma(GshRes[1],GscRes[1])  # This should be correct Gamma distributed soilmoisture. We need the quantile to estimate the celerity
#g = Gamma(SMSh ,SMSc)  # Gamma distributed soilmoisture. We need the probability to estimate the celerity this is wrong should drawn from original SM
prob = cdf(g, (satur + ext_precip)) # we extract the probability of the saturation must be higher than the linear version and is, besides, more correct 
# the celerity should be due to drawn saturation plus the precip
#println("prob=", prob)

probmax = 0.99
probmin = 0.01
if (prob > probmax)
 prob = probmax
end
if (prob < probmin)
 prob = probmin
end
mVel = GshInt*GscInt*Ltymid[1]/Timeresinsec
#println("mVel= ", mVel)

#Dynamic Subsurface Celerities
k[1,2] = CeleritySubSurface1UH(GshInt, GscInt,prob,Ltymid[1],Timeresinsec ) # Gammadistributed Celerity of subsurface flow
k[2,2] = k[1,2]
#println("k[1,2]= ", k[1,2])


# Fixed Subsurface celerities 
#for Lst in 1:Lty
#  k[Lst,1:NoL] = CeleritySubSurface(NoL, Gshape, Gscale, Ltymid[Lst], Timeresinsec) # Celerity of subsurface (and overland) flow
#end

cel = round(k[1,2], digits = 6)
# println("cel= ",round(k[1,2], digits = 6)) 
#Overlandflow celerities

k[1,1] = 2*CeleritySubSurface1UH(GshInt, GscInt,0.99,Ltymid[1],Timeresinsec ) # Overland flow celerites P. Fixed- non dynamical
k[2,1] = k[1,1]                                                                # Overland flow celerites IP Fixed- non dynamical
########################################################################################

#We need to place all the estimation of vectorlengths etc after the celerity has been estimated
#Unit hydrographs for P based on celerities and distance distributions,exponentially distributed
antHorlag = zeros(Int64,Lty) # the absolute slowest layer
for Lst in 1:Lty
    antHorlag[Lst] = Int(trunc(Ltymax[Lst]/(k[Lst,NoL]*Timeresinsec))+1)       # +1 includes the final day.  NB only in soils part
end

nodaysvector = zeros(Int64,Lty, NoL)                        # Number of timesteps to drain each layer
for Lst in 1:Lty 
    nodaysvector[Lst,1:NoL] = Int.(trunc.(Ltymax[Lst]./k[Lst,1:NoL]./Timeresinsec).+1)   # integer time steps
end

layerUH_P = zeros(NoL,antHorlag[1])
layerUH_IP = zeros(NoL,antHorlag[2])  

for i in 1: NoL  
  for Lst in 1:1  
   layerUH_P[i,1:nodaysvector[Lst,i]] .= SingleUH(k[Lst,i], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
  end
  for Lst in 2:2 
   layerUH_IP[i,1:nodaysvector[Lst,i]] .= SingleUH(k[Lst,i], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
  end
end
    
#Groundwater layers; 2dim levels, 1 fastest, NoL slowest#
LayersP = zeros(NoL,antHorlag[1])
LayersIP = zeros(NoL,antHorlag[2])
LayersInit = zeros(NoL,antHorlag[1])

#Estimating current capacity in Layers before Evapotranspiration
ddistx[1,1:NoL] = LayerCapacityUpdate(LayersP, nodaysvector[1,1:NoL], Magkap[1,1:NoL], NoL)
ddistx[2,1:NoL] = LayerCapacityUpdate(LayersIP, nodaysvector[2,1:NoL], Magkap[2,1:NoL], NoL)    
 
ICapdummy[1:2] .= 1000.0 #infinite infiltration capacity to establish subsurface states Does not apply, ICap is not 
    # used in GrvInputDistributionICap 
#println("distx1_1 ",ddistx[1,1], " ",ddistx[1,2]) #OK 13.06
#println("distx2_1 ",ddistx[2,1], " ",ddistx[2,2]) #OK 13.06   
#outx[1] = grvstate[1] # state for subsurface for P areas
#outx[2] = grvstate[1] # 0.0 # subsurface state for IP areas, assuming precipitation on IP areas goes directly to runoff

# Initial Subsurface state from completely dry conditions before rainfall. 
# Estimating input to Layers according to capacity, Note ICap dummy for getting the water in, only for this time
  for Lst in 1:Lty  
   ddist[Lst,1:NoL] = GrvInputDistributionICap(grvstate[1], NoL, ddistx[Lst,1:NoL], ddist[Lst,1:NoL], ICapdummy[Lst])
  end
#println("ddist1_1 ",ddist[1,1], " ",ddist[1,2])
#println("ddist2_1 ",ddist[2,1], " ",ddist[2,2])

#initializing the saturation Layers, putting water in - groundwater state, needs ddist
    for Lst in 1:Lty
      if(Lst==1)
         LayersP = LayerInitUrbanDesign(ddist[Lst,1:NoL],grvstate[1], LayersP, layerUH_P, nodaysvector[Lst,1:NoL], NoL)
    end       
      if(Lst==2)
        LayersIP = LayerInitUrbanDesign(ddist[Lst,1:NoL],grvstate[1], LayersIP, layerUH_IP, nodaysvector[Lst,1:NoL], NoL)       
      end
    end  
# We are set, stauration wise, for  extrem precipitation input

LayersInit = LayersP   # Routing LayersInit will illustrate the runoff solely due to the initial moisture conditions 

#println("layerUH_P, i er 1= ",round(layerUH_P[1,1],digits=3), " ", round(layerUH_P[1,2],digits=3), " ", round(layerUH_P[1,3],digits=3)," ", round(layerUH_P[1,4],digits=3),
#" ", round(layerUH_P[1,5],digits=3)," ", round(layerUH_P[1,6],digits=3)," ", round(layerUH_P[1,7],digits=3)," ", round(layerUH_P[1,8],digits=3)," ", round(layerUH_P[1,9],digits=3)) 

#println("layerUH_P, i er NOL= ",round(layerUH_P[2,1],digits=3), " ", round(layerUH_P[2,2],digits=3), " ", round(layerUH_P[2,3],digits=3)," ", round(layerUH_P[2,4],digits=3),
#" ", round(layerUH_P[2,5],digits=3)," ", round(layerUH_P[2,6],digits=3)," ", round(layerUH_P[2,7],digits=3)," ", round(layerUH_P[2,8],digits=3)," ", round(layerUH_P[2,9],digits=3)) 

#println("LayersP etter grvstate= ",sum(LayersP), " grvstate  ", grvstate[1])  # OK! 13.06
#println("LayersP etter grvstate= ",round(LayersP[2,1],digits=2), " ", round(LayersP[2,2],digits=2), " ", round(LayersP[2,3],digits=2)," ", round(LayersP[2,4],digits=2),
#" ", round(LayersP[2,5],digits=2)," ", round(LayersP[2,6],digits=2)," ", round(LayersP[2,7],digits=2)," ", round(LayersP[2,8],digits=2)," ", round(LayersP[2,9],digits=2)) 


#println("layersIP ",sum(LayersIP), " grvstate ", grvstate[1]) # OK 13.06
    
#Random element 2, temporal distribution of extreme precipitation
# Creating weights for precipitation disaggregation, drawn from a beta distribution
shp1 = rand(1:4) # 3 as maximum secure that none of the values are too large/small 
shp2 = rand(1:4) #
Precseq = [0: numwgt: 1;] # sequence of weights, how many weights for the duration
precipdist = Beta(shp1, shp2)          # Draw precipitatin distribution
wgtprec = zeros(length(Precseq)-1) #the actual weights is the average of the subssequebt betadistributed density one less than the number of drawn values
PBeta = pdf.(precipdist,Precseq) # draws the density values 
for k in 2: (length(wgtprec)+1)    # estimates the acutal weights, mean of two consequitive drawn values  
   wgtprec[k-1] = (PBeta[k-1]+PBeta[k])*0.5
end
wgtprec .= wgtprec/sum(wgtprec) # regularize to make proper weights
    
#println("Sum of prec weights = ", sum(wgtprec))    #   OK! 13.06
    
# End for Creating weights for precipitation disaggregation
    
#Random element 3, value of extreme precipitation      
#ext_precip = rand(GeneralizedPareto(μ, σ, ξ),1) # shape parameter ξ, scale σ and location μ draw extreme precipitaton
#ext_precip =  rand(GeneralizedPareto(GPloc, GPsc, GPsh),1) 
ext_precipvec[1:days] = wgtprec .* ext_precip # precipitation distributed over Temp. resol for runoff

############################################################################
#                                                                          #
#Runoff simulation starts here.............. days is no. time steps               #
############################################################################
for i in 1:days2 #Length of precip timeseries +5 (recession)
    
  for Lst in 1:Lty     # landscape types, one snow regime for each landscape type
      gmloutx[Lst] = ext_precipvec[i]    #precipitation and snowmelt     
  end #For landscape types snow and rain
  
    #if(i == 1)  Shall not shift contents of Layers
     ddistx[1,1:NoL] = LayerCapacityInit(LayersP, nodaysvector[1,1:NoL], Magkap[1,1:NoL], NoL)
     ddistx[2,1:NoL] = LayerCapacityInit(LayersIP, nodaysvector[2,1:NoL], Magkap[2,1:NoL], NoL)       
    #end        
    ##Estimating current capacity in Layers after initial drawn saturation
    #if(i > 1)
    #  ddistx[1,1:NoL] = LayerCapacityUpdate(LayersP, nodaysvector[1,1:NoL], Magkap[1,1:NoL], NoL)
    #  ddistx[2,1:NoL] = LayerCapacityUpdate(LayersIP, nodaysvector[2,1:NoL], Magkap[2,1:NoL], NoL)
    #end
    
    #println("after init ss ddistx1_2 ",ddistx[1,1], " ",ddistx[1,2]) #OK 13.06
    #println("after Init ss ddistx2_2 ",ddistx[2,1], " ",ddistx[2,2]) #OK 13.06 
   
    #Updating the deficit (for all sub surface layers, NOT overland flow layer) (before Evapotranspiration, NOT applicable)
    for Lst in 1:Lty
      totdef[Lst] = sum(ddistx[Lst,2:NoL])  # This may be negative if input outstrips deficits
      Def[Lst] = max(totdef[Lst],0)
    end  
  
  # Assigning overland flow if infiltration capapcity is exceeded and reduces outx in the process
  # OF is only overland flow due to infiltration excess      
  for Lst in 1:Lty
      OF[Lst],outx[Lst] = OFICap(gmloutx[Lst],ICap[Lst])
  end       
   #println("First OFdue to ICap, j ",OF[1]," ",OF[2], " ", j)
   #println("gmloutx ",gmloutx[1]," ",gmloutx[2])
        
  #Estimating input to Layers according to capacity. ddist always initially zero. Estimates possible OF due to Sat Excess
  #Test of OF due to SatExcess is if ddist[1 or 2 ,1] > 0.0
  for Lst in 1:Lty  
   ddist[Lst,1:NoL] = GrvInputDistributionICap(outx[Lst], NoL, ddistx[Lst,1:NoL], ddist[Lst,1:NoL], ICap[Lst])
  end       
  #println("1 outx ",outx[1]," ",outx[2])
  #println("1 ddistx ",ddistx[1,1]," ",ddistx[1,2]," ",ddistx[2,1]," ",ddistx[2,2])      
  #println("1 ddist Sat excess ",ddist[1,1]," ",ddist[1,2]," ",ddist[2,1]," ",ddist[2,2])
        
  #if(ddist[1,1] >0)
  #    println("Hei ", j)
  #end     
        
  #We need to include OF due to exceeeded infiltration capacity and UPDTATE dist and outx
  for Lst in 1:Lty  
   if(OF[Lst] > 0.0) # must have this condition or maybe division by zero
     ddist[Lst,1] = ddist[Lst,1] + (OF[Lst]/(OF[Lst]+outx[Lst])) # will include OF from Inf_Ex
     ddist[Lst,2:NoL] = (1-ddist[Lst,1]).*ddist[Lst,2:NoL]  #updating the rest of ddist if case of OF
     outx[Lst]= outx[Lst]+ OF[Lst]    
   end
  end      
  #println("2 outx ",outx[1]," ",outx[2])
  #println("2 ddist corrected due to Inf_excess ",ddist[1,1]," ",ddist[1,2]," ",ddist[2,1]," ",ddist[2,2])
        
#waterbalance calculations Initial soilmoisture pluss precip
#At this point nothing is put in Layers yet, not overland flow or otherwise
  #RinnP =sum(QRivxP[2:noDT])*(Timeresinsec*1000/area[1]) # P water stored in RN from last TimeSep
  #RinnIP= sum(QRivxIP[2:noDT])*(Timeresinsec*1000/area[2]) # IP water stored in RN from last TimeStep
        
  #SPinn1 = sum(LayersP[1,1:nodaysvector[1,1]]) + sum(ddist[1,1] .* outx[1] .* layerUH_P[1,1:nodaysvector[1,1]])#lag1
  #SPinn2 = sum(LayersP[2,1:nodaysvector[1,2]]) + sum(ddist[1,2] .* outx[1] .* layerUH_P[2,1:nodaysvector[1,2]])#lag2 
  #SIPinn1 = sum(LayersIP[1,1:nodaysvector[2,1]]) + sum(ddist[2,1] .* outx[2] .* layerUH_IP[1,1:nodaysvector[2,1]])
  #SIPinn2  = sum(LayersIP[2,1:nodaysvector[2,2]]) + sum(ddist[2,2] .* outx[2] .* layerUH_IP[2,1:nodaysvector[2,2]])           
        
  # Water to runoff this timestep in [mm]
  GDTInit = sum(LayersInit[1:NoL,i])                                                      # runoff this timestep solely due to  initial moisture conditions, goes for both Lsts     
  GDT_P = sum(LayersP[1:NoL,1]) + sum(ddist[1,1:NoL] .* outx[1] .* layerUH_P[1:NoL,1])    # groundwater to be discharged into the rivernetwork + this timesteps contribution
  GDT_IP = sum(LayersIP[1:NoL,1]) + sum(ddist[2,1:NoL] .* outx[2] .* layerUH_IP[1:NoL,1]) # groundwater to be discharged into the rivernetwork 
  SSF_P = sum(LayersP[2:NoL,1]) + sum(ddist[1,2:NoL] .* outx[1] .* layerUH_P[2:NoL,1])    # subsurface contribution (Ex OF)- to assess how important is the infiltrated contributions
  SSF_IP = sum(LayersIP[2:NoL,1]) + sum(ddist[2,2:NoL] .* outx[2] .* layerUH_IP[2:NoL,1]) # subsurface contribution (Ex OF)- to assess how important is the infiltrated contributions
  SS_contribP = sum(LayersP[2:NoL,1])                                                     # subsurface contribution- to assess how important is the SS contrib? Not sure abot this
  SS_contribIP = sum(LayersIP[2:NoL,1])                                                   # subsurface contribution- to assess how important is the SS contrib? Not sure abot this
  P_contribP = sum(ddist[1,1:NoL] .* outx[1] .* layerUH_P[1:NoL,1])                       # precip contribution - to assess how important is the precip contrib? 
  P_contribIP = sum(ddist[2,1:NoL] .* outx[2] .* layerUH_IP[1:NoL,1])                     # precip contribution - to assess how important is the precip contrib? 
        
  GDT_OF = Ltyfrac[1]*(LayersP[1,1]+(ddist[1,1]*outx[1]*layerUH_P[1,1]))+
           Ltyfrac[2]*(LayersIP[1,1]+(ddist[2,1]*outx[2]*layerUH_IP[1,1]))                 # Overland flow part, Satur and infil. excess, contrib from both LSTs
        
  GDT_OFP = (LayersP[1,1]+(ddist[1,1]*outx[1]*layerUH_P[1,1]))                             # Overland flow P
  GDT_OFIP = (LayersIP[1,1]+(ddist[2,1]*outx[2]*layerUH_IP[1,1]))                          # Overland flow IP
        
  #Updating the saturation Layers. The "boxes are shifted on step ahead"
    for Lst in 1:Lty
      if(Lst==1)
         LayersP = LayerUpdate(ddist[Lst,1:NoL],outx[Lst], LayersP, layerUH_P, nodaysvector[Lst,1:NoL], NoL)
    end       
      if(Lst==2)
        LayersIP = LayerUpdate(ddist[Lst,1:NoL],outx[Lst], LayersIP, layerUH_IP, nodaysvector[Lst,1:NoL], NoL)       
      end
    end  
   
  #summing up groundwater states after this timesteps water has left. antHorlag is just the slowest nodaysvector
    for Lst in 1:Lty
      if(Lst == 1)
         lyrs[Lst] = sum(LayersP)
         subsurface[Lst] = sum(LayersP[2:NoL,1:antHorlag[1]]) # gives the sum of layers after todays runoff has taken place
      end
      if(Lst==2)
         lyrs[Lst] = sum(LayersIP)
         subsurface[Lst] = sum(LayersIP[2:NoL,1:antHorlag[2]])
      end
    end
   #println("AntHorlag 1 =", antHorlag[1] )         #14.06 OK
   #println("AntHorlag 2 =", antHorlag[2] )         #14.06 OK
   #println("nodaysvector 1 =", nodaysvector[1,2])  #14.06 OK
   #println("nodaysvector 2 =", nodaysvector[2,2])  #14.06 OK

   #println("UHriver =", UHriver )         #14.06 OK

  #runoff
  qsimutxInit .= (((GDTInit/1000)*(area[1]+area[2]))/Timeresinsec) .* UHriver     #[m3/s] runof solely due to initial moisture conditions 
  #println("qsimutxInit ",i," =", round(qsimutxInit[1] , digits= 4)," ",round(qsimutxInit[2], digits= 4))
  qsimutxP .= (((GDT_P/1000)*area[1])/Timeresinsec) .* UHriver     #[m3/s]
  #println("qsimutxP ",i," =", round(qsimutxP[1] , digits= 4)," ",round(qsimutxP[2] , digits= 4))
  qsimutxIP .= (((GDT_IP/1000)*area[2])/Timeresinsec) .* UHriver   #[m3/s]  
  qsimutxOF .= (((GDT_OF/1000)*totarea)/Timeresinsec) .* UHriver   #[m3/s]  # Overland flow. Not a contribution, just the fraction OF of total runoff
  qsimutxOFP .= (((GDT_OFP/1000)*area[1])/Timeresinsec) .* UHriver   #[m3/s]  # Overland flow P. Not a contribution, just the fraction OF of total runoff
  qsimutxOFIP .= (((GDT_OFIP/1000)*area[2])/Timeresinsec) .* UHriver   #[m3/s]  # Overland flow IP. Not a contribution, just the fraction OF of total runoff        
  qsimutxSSFP .= (((SSF_P/1000)*area[1])/Timeresinsec) .* UHriver   #[m3/s]  # Subsurface flow P. Not a contribution, just the fraction OF of total runoff
  qsimutxSSFIP .= (((SSF_IP/1000)*area[2])/Timeresinsec) .* UHriver   #[m3/s]  # Subsurface flow IP. Not a contribution, just the fraction OF of total runoff  
  qsimutxSS_contribP = (((SS_contribP/1000)*area[1])/Timeresinsec) .* UHriver   #[m3/s]# Not a contribution, just to assess the fraction of SS vs P
  qsimutxSS_contribIP = (((SS_contribIP/1000)*area[2])/Timeresinsec) .* UHriver   #[m3/s] Not a contribution, just to assess the fraction of SS vs P
  qsimutxP_contribP = (((P_contribP/1000)*area[1])/Timeresinsec) .* UHriver   #[m3/s] Not a contribution, just to assess the fraction of SS vs P
  qsimutxP_contribIP = (((P_contribIP/1000)*area[2])/Timeresinsec) .* UHriver   #[m3/s]  Not a contribution, just to assess the fraction of SS vs P    
  #println(round(sum(qsimutxOF),digits=6)," ",round(sum(qsimutxOFP)+sum(qsimutxOFIP),digits=6))
        
  #Total response
  qsimutx[1:noDT] .= 0.0                                   # this is written anew for each timestep
  qsimutx[1:noDT] .= qsimutx[1:noDT] .+ qsimutxP[1:noDT]   #adding contribution from permeable areas
  qsimutx[1:noDT] .= qsimutx[1:noDT] .+ qsimutxIP[1:noDT]  #adding contribution from impermeable areas

  #Combined (p and IP) contribution from SS and from P (P this timestep)
  qsimutxSS_contrib .= qsimutxSS_contribP .+ qsimutxSS_contribIP
  qsimutxP_contrib .= qsimutxP_contribP .+ qsimutxP_contribIP
        
  #println("noDT = ", noDT)              # Ok 14.06.2023
  #println("UHriver = ", length(UHriver)) # Ok 14.06.2023, same as the above     
  
  #Updating the routing in river network from iniyil moisture conditions
  QRivxInit, QRDInit = RiverUpdate(noDT,QRivxInit,qsimutxInit)

  #Updating the routing in river network total)
  QRivx, QRD = RiverUpdate(noDT,QRivx,qsimutx) # During update, QRivx from previous timestpæi sshifted and todays contrib is added
                                                                  # RN is also a reservoir [mm], water from todays event is stored for future runoff in the RN, relevant for WB.
  qmm_state = -9999
  #Updating the routing in river network P
  QRivxP,QRDP = RiverUpdate(noDT,QRivxP,qsimutxP)
    
  #Updating the routing in river network IP
  QRivxIP,QRDIP = RiverUpdate(noDT,QRivxIP,qsimutxIP)
  
  #Updating the routing in river network OF
  QRivxOF,QRDOF = RiverUpdate(noDT,QRivxOF,qsimutxOF)
  
  #Updating the routing in river network OFP
  QRivxOFP,QRDOFP = RiverUpdate(noDT,QRivxOFP,qsimutxOFP)
        
  #Updating the routing in river network OFIP
  QRivxOFIP,QRDOFIP = RiverUpdate(noDT,QRivxOFIP,qsimutxOFIP)
  
  #Updating the routing in river network SSF_P
  QRivxSSF_P,QRDSSF_P = RiverUpdate(noDT,QRivxSSF_P,qsimutxSSFP)
        
  #Updating the routing in river network SSF_IP
  QRivxSSF_IP,QRDSSF_IP = RiverUpdate(noDT,QRivxSSF_IP,qsimutxSSFIP)  
  
  #Updating the routing in river network for SS_contrib
  QRivxSS_contrib, QRDSS_contrib = RiverUpdate(noDT,QRivxSS_contrib,qsimutxSS_contrib)  

  #Updating the routing in river network for P_contrib
  QRivxPrecip_contrib, QRDP_contrib = RiverUpdate(noDT,QRivxPrecip_contrib,qsimutxP_contrib)  
        
  Qmm = (QRD*Timeresinsec*1000/totarea)  #QRD in mm/timestep
  
  #Waterbalance calculations InitLayers plus precip = Stored water plus Discharge                                          
  # Stored water 1)Initial Layers minus this timeteps discharge 2)plus precip to be distributed in future
                                            #timesteps 3) water stored in channel
   #WBP = (SPinn1 + SPinn2 + RinnP) - (lyrs[1] + sum(QRivxP)*(Timeresinsec*1000/area[1]))   #mm  last term inludes discharge and strage in rivernetwork                                         
   #println("WBP = ", WBP, " Total inn P= ",SPinn1+ SPinn2)                                             
   #WBIP = (SIPinn1 + SIPinn2+ RinnIP) - (lyrs[2] + sum(QRivxIP)*(Timeresinsec*1000/area[2]))   #mm  last term inludes discharge and strage in rivernetwork                                         
   #println("WBIP = ", WBIP," Total inn IP= ",SIPinn1+ SIPinn2)                                             
                                                  
  #Assigning outdata to vector, one vector for each timestep
   simresult[i, 1:5] = [i,round(ext_precip[1],digits=6),round(QRD,digits=6),
                round(QRDP,digits = 6),round(QRDIP,digits = 6)]
   simresult[i, 6:9] = [round(lyrs[1],digits = 6),round(lyrs[2],digits = 6),round(totdef[1],digits=2), 
                round(totdef[2],digits=2)]
   simresult[i, 10:11] = [round(Qmm,digits=6),round(qmm_state,digits=6)]
   simresult[i, 12:14] = [round(subsurface[1],digits=6), round(subsurface[2],digits=6),round(QRDOF,digits=6)]
   simresult[i, 15:16] = [round(QRDOFP,digits=6),round(QRDOFIP,digits=6)]
   simresult[i, 17:18] = [round(QRDSSF_P,digits=6),round(QRDSSF_IP,digits=6)]         
   simresult[i, 19:21] = [round(grvstate[1],digits=6),shp1, shp2]
                                  
  qberegn[i] = QRD      #Simulated runoff for current timestep. Used for skillscores etc.
  lyrberegn[i] = lyrs[1]
  SSberegn[i] = round(QRDSS_contrib,digits=6)
  Pberegn[i] = round(QRDP_contrib,digits=6)
  Initberegn[i] = round(QRDInit,digits=6)
end # for days2, i 
    
   maxQ = maximum(simresult[1:days2,3])
   imax = findall(simresult[1:days2,3] .== maxQ)
   maxP = maximum(ext_precipvec[1:days2])
   Pwhere = findall(ext_precipvec[1:days2] .== maxP)
    
    simresult2[((l-1)*NumSim)+j,1:21] = simresult[imax[1],1:21]
    simresult2[((l-1)*NumSim)+j,22:24] = [Pwhere[1],round(ext_precipvec[Pwhere[1]],digits=6), round(cel, digits=6)] 
    simresult3[((l-1)*NumSim)+j,1:days2] = ext_precipvec[1:days2]
    simresult3[((l-1)*NumSim)+j,(days2+1):(days2+days2)]= qberegn[1:days2]
    simresult3[((l-1)*NumSim)+j,(days2+days2+1):(days2+days2+days2)]= lyrberegn[1:days2]
    simresult3[((l-1)*NumSim)+j,(days2+days2+days2+1):(days2+days2+days2+days2)]= Initberegn[1:days2] # SSberegn[1:days2]            # Runoff contribution fra SS
    simresult3[((l-1)*NumSim)+j,(days2+days2+days2+days2+1):(days2+days2+days2+days2+days2)]= qberegn[1:days2].-Initberegn[1:days2] # Runoff contribution fra Precip
    #simresult3[j,(days2+days2+days2+days2+days2+1):(days2+days2+days2+days2+days2+days2)]= Initberegn[1:days2]
    #Need to start with a dry cathment again
       
        end # for j,stochastisc rounds and progressbar  
 
end # for k, fixed precip values simulated Numsim times

#Writing to files
#wwater = 1.0*persons*(140.0 + 42.0)/(86400.0*1000.0)# norsk vann equivalent use in L/day add 30% leakage from water net
toptitles = ["TimestpMaxFld","ExtPrecipDur","maxQ","maxQ_P","maxQ_IP","SST_P","SST_IP","Def_P","Def_IP","Qmm", "XX",
    "SS_P","SS_IP","OF", "OF_P","OF_IP","SSF_P","SSF_IP","InitSS", "Shp1", "Shp2", "TimesetpMaxPrec","PrecipMax", "SScel"]
  toptitles2 = ["Prec1","Prec2","Prec3","Prec4","Prec5","Prec6","Prec7","Prec8","Prec9","Prec10",
  	"Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","SS1","SS2","SS3","SS4","SS5","SS6","SS7","SS8","SS9","SS10",
    "Qinit1","Qinit2","Qinit3","Qinit4","Qinit5","Qinit6","Qinit7","Qinit8","Qinit9","Qinit10","Q-Qinit1","Q-Qinit2",
    "Q-Qinit3","Q-Qinit4","Q-Qinit5","Q-Qinit6","Q-Qinit7","Q-Qinit8","Q-Qinit9","Q-Qinit10"]

  CSV.write(utfile,DataFrame(simresult2,:auto), delim = ';', header=toptitles) 
  CSV.write(utfile2,DataFrame(simresult3,:auto), delim = ';',header=toptitles2) 


t2 = time_ns()

println("Number of timesteps in RN= ", noDT)
println("Location par GDP = ",GPloc)
println("Scale par GDP= ", GPsc)
println("Shape par GPD= ", GPsh)

return 

end  #end of func


