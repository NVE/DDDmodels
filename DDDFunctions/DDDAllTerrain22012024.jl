#DDD model Module based for hydrological predictions in urban and rural areas
#
#----------------------------------------------------------------------------------
#     Script name:  DDDAllTerrain22012024.R
#     Author:     Thomas Skaugen
#     Revised: 22.1.2024
#              
# Purpose   : DDD (Distance Distribution Dynamics) hydrological model simulates
#               1) saturated soil water flow, states of saturated (2-D)) and unsaturates zones
#               2) snow accumulation, distribution and melt
#               3) runoff from permable and impermeable surface for a given catchment.
#               4) saves the states for each day where precip > 0, i. e where SCHADEX might insert simulated extreme precip
#               5) Uses EB for snowmelt and evapotranspiration
#               6) estimates delay due to lakes
#               7) Estimates Mean daily discahrge (MAD) from regression equations if fraction of glaciers is large
#               8) Calculates Soilmoisture as a function of distace from river netwrok
#
# Model input: 1) catchment physical/geographical description, distance distributions for landscape features and river network
#              2) weather data, temperature,precipitation 
#              3) model states variables (optional)
#
# Model output: simulated 1)saturated and unsaturated soil water
#                         2)snow water equivalent and snow distribution
#                         3)river discharge
#                         4)state files for each rainy day. Remember to include temperature for 5 days t-1, t=0, t+1, t+2, t+3 
#
# Running DDD:    It is possible to save model state variables and run the model
#                 starting from saved state variables
#----------------------------------------------------------------------------------

function DDDAllTerrain(Gpar, startsim, tprm, prm ,ptqfile,utfile,r2fil, modstate,savestate, kal, spinuptime)


DDA = 6  # number of landscape types with distance distribution
# DDA=1 Permeable (P) areas
# DDA=2 ImPermeable (IP) areas
# DDA=3 Wetlands 
# DDA=4 Lakes (effective lake fraction)
# DDA=5 Glaciers  NB, Glaciated area is considered as P area for runoff dynamics
# DDA=6 Rivernetwork
#  
Lty = 2  # Areas with a subsurface and a snow distribution Lty[1] includes DDA[1,3,4,5] DDA[5] has no aeal extention 

#-------------------------------------preprocess start------------------------------------------

#####################VARIABLES FOR 10 ELEVATION ZONES and for Parts ###########################################
isoil = zeros(Lty,10)# precipitation and snowmelt from the elevation zones
misoil = zeros(Lty)# summed isoil for all elevation sones
gwgt = zeros(10)               # weights for glaciermelt pr each elevation zone
swgt = zeros(Lty,10)# weights for input to soils pr each elevation zone NB bogs are not part of this yet
spd = zeros(Lty,10)
swe_h = zeros(Lty,10)   #SWE pr altitude level
wcd = zeros(Lty,10)
sca = zeros(Lty,10)
nsno = zeros(Lty,10)
alfa = zeros(Lty,10)
ny = zeros(Lty,10)
snowfree = zeros(Lty,10)  #Indicator variable: 1 if sca  < gca and 0 if sca > gca
#tempvec = zeros(Lty,10)   # temperature vector  to save in statefile
PR = zeros(Lty,10)
PS = zeros(Lty,10)
MW = zeros(Lty,10)
SWrad = zeros(Lty,10)
LA = zeros(Lty,10)
LT = zeros(Lty,10)
Pa = zeros(10)
MWGLAC= zeros(Lty,10)
MWGLAC1 = zeros(Lty,10)
gisoil = zeros(Lty,10)   # glacialmelt from the elevation zones

#Addition for Enegy balance snowmelt
snro = zeros(Lty,10)            #Snow density
melt = zeros(Lty,10)           #instead of degreeday melt (MW)
snowdepth = zeros(Lty,10)        #hmm..
sndens = zeros(Lty,10)           #hm ...
#################################Parameters as vectors#############################################
a0 = zeros(Lty) # scale parameter unit snow
n0 = zeros(Lty) # shape parameter unit snow
d = zeros(Lty)

######################    Parameters for Distance distributions    #######################################
Ltyfrac = zeros(DDA) # areal frac.tion (AF) of permeable, impermeable, wetlands, glaciers and lakes(effective). AF of RN is negligible 
Ltymax = zeros(DDA)    # includes all landscapetypes Lakes are caluclate in program
Ltymid = zeros(DDA)    # includes all landscapetypes Lakes are caluclate in program
Ltystd = zeros(DDA)    # includes RN and Glaciers others are set to zero
Ltyz = zeros(DDA)      # includes P,IP and Wetlands, others are set to zero   

#################################States as vectors#################################################
M = zeros(Lty)
ICap = zeros(Lty)       # Infiltration capacity foor P[1] and IP[2] in [mm/s]
CritFlux = zeros(Lty)   # critical flux for overland flow
taux = zeros(Lty)       # albedo reduction in AlbedoUEB
lyrs = zeros(Lty)       # subsurface including overland flow
subsurface = zeros(Lty) # subsurface exluding overland flow
area = zeros(DDA)     # area estimates for all landscapetype
totdef = zeros(Lty)
Def = zeros(Lty)
PotEvap = zeros(Lty)
tosoil = zeros(Lty)
gmltosoil = zeros(Lty)
outx = zeros(Lty)
sm = zeros(Lty)
ea = zeros(Lty)
ea_sm = zeros(Lty)
ea_S = zeros(Lty)
OF = zeros(Lty)
Snowdepth = zeros(Lty)
Albedo = zeros(Lty)
snittT = zeros(Lty)
MLT = zeros(Lty)
MLA = zeros(Lty)
MSWrad = zeros(Lty)
snomag = zeros(Lty)
middelsca = zeros(Lty)
snofritt = zeros(Lty)
wcs = zeros(Lty)
GIsoil = zeros(Lty)  # mean glacial melt 
scaobx = zeros(10)
hprecip = zeros(10)
htemp = zeros(10)
hfelt = zeros(10)

########################Reading parameters#########################################################
#Hyposgraphic curve and locations
hfelt[1] = prm.val[5]+(prm.val[6]-prm.val[5])/2.0
hfelt[2] = prm.val[6]+(prm.val[7]-prm.val[6])/2.0
hfelt[3] = prm.val[7]+(prm.val[8]-prm.val[7])/2.0
hfelt[4] = prm.val[8]+(prm.val[9]-prm.val[8])/2.0
hfelt[5] = prm.val[9]+(prm.val[10]-prm.val[9])/2.0
hfelt[6] = prm.val[10]+(prm.val[11]-prm.val[10])/2.0
hfelt[7] = prm.val[11]+(prm.val[12]-prm.val[11])/2.0
hfelt[8] = prm.val[12]+(prm.val[13]-prm.val[12])/2.0
hfelt[9] = prm.val[13]+(prm.val[14]-prm.val[13])/2.0
hfelt[10] = prm.val[14]+(prm.val[15]-prm.val[14])/2.0

phi = prm.val[16]*pi/180  #converts Lat degrees to radians lat
thi = prm.val[17]*pi/180  #converts Lon degrees to radians Lon 

#Meteorological corrections
pkorr = tprm[4] #prm.val[18]         # Met corrections rain
skorr = tprm[5] #prm.val[19]                  # Met corrections snow
u = tprm[1] #prm.val[20]             # mean vind velocity [m/s]

#Snow parameters
pro = tprm[2] #prm.val[21]           # [fraction] liquid water content in snow
TX = tprm[3] #prm.val[22]            # Threshold temp for rain/snow
a0[1:Lty] .= prm.val[23:24]          # Alfa null in snofall unit distribution P, IP
d[1:Lty] .= prm.val[25:26]           # Decorrelation length (Measured in n) for precipitation P, IP

# Misc.
Timeresinsec = prm.val[27]           # temporal resolution of simulations, in seconds
MAD = prm.val[28]                    # mean annual discharge
totarea = prm.val[29]                # Total area of catchment in m2
NoL = Int(prm.val[30])                # Number of layers in subsurface including OF level
R = prm.val[31]                  # soil moisture content/100 for field capacity of soils 

#Celerities, subsurface and conduits
GshInt = prm.val[32]        # shapeparameter saturation Capital lambda
GscInt = tprm[6] #prm.val[33]        # scaleparameter saturation Capital lambda
Gshape, Gscale = Big2SmallLambda(GshInt, GscInt) # Coverting integrated celerity to layers 
OFVP  = tprm[7]  #prm.val[34]         # Overland flow velocity P
OFVIP = tprm[8]  #prm.val[35]                  # Overland flow velocity IP
Lv = tprm[9]     #prm.val[36]                     # lake celrity [m/s] Can it be fixed?, 0.01 is popular
rv = tprm[10]    #prm.val[37]                     # celerity of water in river/conduits network    

#Distance distributions Lty
Ltyfrac[1:5]  .= prm.val[38:42]  # areal fraction (AF) of permeable, impermeable, wetlands, lakes(effective) and glaciers. Not RN 
Ltymax[1:3] .= prm.val[43:45]    # Max. dist permeable, impermeable, wetlands, lakes(effective), glaciers and RN 
Ltymax[5:DDA] .= prm.val[46:47]    # Max. dist permeable, impermeable, wetlands, lakes(effective), glaciers and RN 
Ltymid[1:3] .= prm.val[48:50]    # Mean permeable, impermeable, wetlands, lakes(effective), glaciers and RN 
Ltymid[5:DDA] .= prm.val[51:52]    # Mean permeable, impermeable, wetlands, lakes(effective), glaciers and RN 
Ltystd[5:DDA] .= prm.val[53:54]    # Std. dev. permeable, impermeable, wetlands, lakes(effective), glaciers and RN 
Ltyz[1:3] .= prm.val[55:57]      # frac zero dist from RN, permeable, impermeable, wetlands, lakes(effective), glaciers and RN
# Glacierfractions of elevation zones
g1    = prm.val[58]       # areal fraction of glaciers in first elevation zone
g2    = prm.val[59]
g3    = prm.val[60]
g4    = prm.val[61]
g5    = prm.val[62]
g6    = prm.val[63]
g7    = prm.val[64]
g8    = prm.val[65]
g9    = prm.val[66]
g10   = prm.val[67]

meandailyP = prm.val[68]   # daily mean value precipitation
meandailyT = prm.val[69]   # daily mean value temperature
snfjell = 100.0*prm.val[70]# in parameterfile as fraction for MAD estimation, NOT for impermeble surfaces
persons = prm.val[71]
    
#Hardcoded infiltration capacity
#ICap[1] = (1000)*(GshInt*GscInt*Ltymid[1]/Timeresinsec) # Infiltration capacity for Permeable[1] [mm/s]. Equals mean Subsurface velocity
#ICap[2] = ICap[1]*0.001  # Infiltration capacity for ImPermeable[2] in

#provided by parameterfile
ICap[1] = prm.val[72]/3600.0#  Infiltration capacity for Permeable[1] Input from parameterfile is in mm/hour. we here transform to [mm/s]and later on to [mm/Timeresinsec] 
ICap[2] = prm.val[73]/3600.0# # Infiltration capacity for ImPermeable[2]  previously as ICap[1]*0.002  

CritFlux[1:Lty] .= prm.val[74:75] # Critical flux estimates. Sensitivity  unknown....

################################### End of reading parameters ############################################

#Merging assumed Glacier river network [5] with observed river network [6]. We assume independent normal distributions for 5 an 6
Ltymax[6] = Ltymax[6] + Ltymax[5]
Ltymid[6] = Ltymid[6] + Ltymid[5]  #mean of flowlength distribution
Ltystd[6] = Ltystd[6] + Ltystd[5]

#Infiltration capacity in mm/s. Must multiply with Timeresinsec
for Lst in 1:Lty 
    ICap[Lst] = ICap[Lst] .* Timeresinsec
end

#if no impermeable areas. Every item associated with Lty=2 (IP areas) should be zero
if(Ltyfrac[2] == 0.0)
  Lty = 1
end 

#reading input file (ptq)
ptqinn = CSV.read(ptqfile, DataFrame, header=["yr","mnt","day","hr","min","p01","p02","p03","p04","p05","p06","p07","p08","p09","p10",
        "t01","t02","t03","t04","t05","t06","t07","t08","t09","t10","q"])    # reading ptq data, elevation zones
                                                                                                                                
days = Int(length(ptqinn.yr))    # Length of timeseries

MAD2 = exp(-4.59 + 1.135*log(meandailyP)+0.97*log(totarea/1000000)-0.0603*meandailyT + 0.053*log(snfjell)-0.00014*hfelt[5])  # R2 = 0.99, Area in km2, Snfjell in fraction   
if (Ltyfrac[5] > 0.05 ) #Fraction of glaciers
     MAD = MAD2    # needs to adjust MAD if glaciers are present since MAD is used to estimate the groundwater reservoir
end 
println("Mad= ", MAD)

#Constants
hson = 10                               # Elevation zones
unitsnow = 0.1                          # mean of unit snow [mm]
n0 = unitsnow*a0                        # shape parameter of unit snow
gtcel = 0.99                            # threshold qunatile for groundwater capacity -> Overland flow
CFR = 2.5*(Timeresinsec/86400)*0.0833   # Fixed as 1/12 of estimate of CX= 2.5 for 24 hours 
len = Int(5*(86400/Timeresinsec))       # number of timestepes to spin up the model. Recommended to use timesteps equal to a minimum of 5 days to estimate the snowpack temperature.)

startsim = startsim + len               # taking into account estimating snowpack temperature
tempstart = zeros(len,10)               # matrix for storing temperatures when starting from states
tempstart = ptqinn[(startsim-len+1):startsim,16:25]  # assigning temperature values
STempvec = zeros(len)              # Temperature vector for estimating snowpack temperature

Pa[1:10] = 101.3*((293 .- 0.0065.*hfelt[1:10])/293.0).^5.26     # Air pressure as a function of height Zhu et al. 2012 and Stoy et al, 2018
MPa = mean(Pa)

#vectors and matrixes
ddist = zeros(Lty,NoL)          # distribution of water in layers
ddistx = zeros(Lty,NoL)         # water present in layers
ddistxtemp = zeros(Lty,NoL)     # temporay of the above
Magkap = zeros(Lty,NoL)
aktMag = zeros(Lty,NoL)         # current water in Layers
k = zeros(Lty,NoL)              # subsurface celerities (including that of overland flow) 
qberegn = zeros(days)            # for calculation of skill scores etc. 

if(kal == 0)
    simresult = zeros(days,39)        # matrix which into the results are written
end 

#grwpointresult <-matrix(0.0, days,29)     # matrix which into the groundwater point results are written

# Areas and weights
area = Ltyfrac.* totarea                    # areas of differnt Lty, a vector
                                            # NOTE that area[1] (permeable) includes area[5] (glaciers)
# area = (1-(bogfrac+glacfrac))*totarea      #area not covered by wetlands or glaciers. These three landscape types have different distancedistributions
# area2 = (1-(bogfrac))*totarea              #area in which we have hillslope process and rivernetwork processes

#elevarea = (totarea/hson)                  # area pr elevationzone
elevarea = ((area[1]+area[3])/hson)        # P 8including glaciers) and Bogs area pr elevationzone
gca = [g1,g2,g3,g4,g5,g6,g7,g8,g9,g10]     # fraction of glaciated area per elevation zone 
soilca = 1 .- gca                          # fraction of non-glaciated area per elevation zone

# gwgt is the fraction of glaciated area pr elevation zone in relation to total glaciated area 
# i.e. fraction of total glaciared area located in each elevation zone
if area[5] > 0                          # if glaciated area > 0.
    gwgt = gca .* elevarea/area[5]  # sums to 1, corrected for areafraction later 
else
    gwgt .= 0
end

# Finds the fraction of soils (and bogs) in each elevation zone in relation to total soil (and bog) area
# fraction of total soil- and bogarea located in each elvation zone  
# Glacer runoff is wil potentially come as an addition only to P area (including bogs) 
swgt[1,1:hson] .= soilca .* elevarea/(area[1]+ area[3]-area[5])   # sums to 1, corrected for areafraction later,rain, snowmelt and from glaciers 
#NOTE again, area[1] inludes area[5]
swgt[2,1:hson] .= 0.1   # 1 .* elevarea/(totarea)          # sums to 1, corrected for areafraction later,rain and snowmelt is only input to IP area
#Ser OK ut for case area[2]==0  31.1.2024

#Subsurface celerities 
for Lst in 1:Lty
    k[Lst,1:NoL] = CeleritySubSurface(NoL, Gshape, Gscale, Ltymid[Lst], Timeresinsec) # Celerity of subsurface (and overland) flow
   if(Lst==1)
    k[Lst,1] = k[Lst,1]*2   #03.07.2025 Double the original profile estimate due to equifinality, ECCO result #OFVP  # P m/s Holden et al. WRR,2008,  # Overland flow celerites P
   end
   if(Lst==2)
    k[Lst,1] = k[Lst,1]*2   #03.07.2025 Double the original profile estimate due to equifinality, ECCO result #OFVIP # IP, Sedyowati et al. 2017      # Overland flow celerites IP
   end 
end

#Overlandflow celerities
bogspeed = k[1,1] # same as for overland flow Permeable areas
                                    
#Unit hydrographs for Wetlands
antBogsteps =Int(trunc(Ltymax[3]/(bogspeed*Timeresinsec))+1) # Number of timsteps to drain wetlands
UHbog = zeros(antBogsteps)                  # Vector for Bog unit hydrograph
if(area[3] > 0)
    UHbog = SingleUH(bogspeed,Timeresinsec, Ltymid[3], Ltymax[3], Ltyz[3])
end
if(area[3] <= 0)
    UHbog = 1.0
end

#Unit hydrographs for P and IP based on celerities and distance distributions,exponentially distributed
antHorlag = zeros(Int64,Lty)
for Lst in 1:Lty
    antHorlag[Lst] = Int(trunc(Ltymax[Lst]/(k[Lst,NoL]*Timeresinsec))+1)       # +1 includes the final day.  NB only in soils part
end

nodaysvector = zeros(Int64,Lty, NoL)                        # Number of timesteps to drain each layer
for Lst in 1:Lty 
    nodaysvector[Lst,1:NoL] = Int.(trunc.(Ltymax[Lst]./k[Lst,1:NoL]./Timeresinsec).+1)   # integer time steps
end

 if(area[1] >0.0)
  layerUH_P = zeros(NoL,antHorlag[1])
 end
 if(area[2] >0.0)
  layerUH_IP = zeros(NoL,antHorlag[2])  
 end 
 
for i in 1: NoL 
  for Lst in 1:Lty 
   if(Lst==1) 
    layerUH_P[i,1:nodaysvector[Lst,i]] .= SingleUH(k[Lst,i], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
   end
   if(Lst==2) 
    layerUH_IP[i,1:nodaysvector[Lst,i]] .= SingleUH(k[Lst,i], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
   end
  end
end

#UH for RIVER/Conduits, normally distributed #6
#number of time units in river routing, rounded up. If set equal to one day, time is actually less
UHriver, noDT, nodaysRiv = SingleNormalUH(rv,Timeresinsec,Ltymid[6],Ltystd[6],Ltymax[6]) #This subroutine extracts a UH normally distributed. 

#UH for LAKES, normally distributed #4 
#Lakerouting parameters [m] are estimated according to the generic ellipse shape of a lake. A=π*a*b where a is the length [m] of the long 
# axis and b the short axis.  If we let b= a/4 we have that a=√(4A/π). The area A is calculated as A= Catchment_area*Fraction_effective_lake

Ltymax[4] = sqrt(4*area[4]/pi) 
Ltymid[4] = Ltymax[4]/2
Ltystd[4] = Ltymid[4]/2
UHLake, noDTLake, nodaysLake = SingleNormalUH(Lv,Timeresinsec,Ltymid[4],Ltystd[4],Ltymax[4]) #This subroutine extracts a UH normally distributed. 

qsimutx = zeros(noDT)      # common vector for runoff from all landscapetypes [m3/s]
qsimutxP = zeros(noDT)     # vector for runoff permeabel i m3/s
qsimutxIP = zeros(noDT)    # vector for runoff impermeable i m3/s
qsimutxBog = zeros(noDT)   # vector for runoff from wetlands bogss
qsimutxOF = zeros(noDT)    # vector for Overland Flow
QRivx = zeros(noDT)        # all qsimutvectors are stored UHriver timesteps back. 
QRivxP = zeros(noDT)
QRivxIP = zeros(noDT)
QRivxBog = zeros(noDT)
QRivxOF = zeros(noDT)

#Groundwater layers; 2dim levels, 1 fastest, NoL slowest#
if(area[1] >0.0)
  LayersP = zeros(NoL,antHorlag[1])
 end
if(area[2] > 0.0)
  LayersIP = zeros(NoL,antHorlag[2]) 
else
  LayersIP =zeros(NoL,1)
end 

BogLayers = zeros(antBogsteps) # no vertical dimension
LakeLayers = zeros(nodaysLake) # no vertical dimension 
P_LakeLayers = zeros(nodaysLake) # no vertical dimension 
IP_LakeLayers = zeros(nodaysLake) # no vertical dimension 
Bog_LakeLayers = zeros(nodaysLake) # no vertical dimension 
OF_LakeLayers = zeros(nodaysLake) # no vertical dimension 

for Lst in 1: Lty
  if(area[Lst]>0.0)
  Magkap[Lst,1:NoL], M[Lst] = LayerEstimation(GshInt,GscInt,Timeresinsec,Ltymax[Lst],Ltymid[Lst], 
        Ltyfrac[Lst]*MAD, Ltyfrac[Lst]*totarea, NoL, gtcel) 
  end   
end
#Note that MAD is for the entire catchment. It is fractioned in the above subroutine

#Extracting areas and heights in pyramide-plot for estimating SS states at distances from or in the RN
Parea = area[1] + area[3]+ area[4]+area[5]   # areas where the DD for Pareas apply: Pareas, 
AreasP, delta_dP = PyrAreas(NoL,Parea,Ltymax[1],nodaysvector[1,1:NoL], layerUH_P, antHorlag[1]) # A matrix of areas[in m2] for each time-step box, 
                                                                 # Important: same dimensions as Layers
                                                                 # The area drained pr time-step box for each layer
if(area[2] > 0.0)
  AreasIP, delta_dIP = PyrAreas(NoL,area[2],Ltymax[2],nodaysvector[2,1:NoL], layerUH_IP, antHorlag[2]) # A matrix of areas[in m2] for each time-step box, 
end

RNwdelta_d = Ltymax[6]/noDT

mindistant = 25 #number distance interval for each montoing of saturation if you change this you have to edit the push! statement at the end of the run
satintP = Int(trunc(Ltymax[1]/mindistant)) #Number of distances from RN for which we will investigate the saturation 
satintIP = Int(trunc(Ltymax[2]/mindistant)) #Number of distances from RN for which we will investigate the saturation 
rnint = Int(trunc(Ltymax[6]/mindistant)) #Number of distances from outlet for which we will investigate the water in the RN
mindist = zeros(mindistant)

#Distributed groundwater Pareas
grwpointresult = DataFrame(Datetime = DateTime[]) # specifyinbg a result dataframe into the groundwater point results are written
for i in 1:mindistant
 mindist[i] = i*satintP  # Distances at which we will monitor the simulated saturation
 insertcols!(grwpointresult,(i+1), string("gw",i)=> 0.0) 
end

rndist = zeros(mindistant) # distances in RN
#Distributed Water in river network
rnpointresult = DataFrame(Datetime = DateTime[]) # specifying a result dataframe wher the rn point results are written
for i in 1:mindistant
 rndist[i] = i*rnint  # Distances at which we will monitor the simulated water in the RN
 insertcols!(rnpointresult,(i+1), string("rn",i)=> 0.0) 
end

#States
 if(modstate == 1) 
    tilstfile = string("\\\\nve.no\\fil\\h\\HM\\Interne Prosjekter\\HydSimOveralt\\Totalmodell\\utdata\\",
    "statefile_",catchment,"_20200116_1220_.hdf5")
    #Run with states. The states from the previous times step is loaded
 @load tilstfile sca spd wcd nsno alfa ny sm smbog LayersP LayersIP BogLayers QRivx QRivxP QRivxIP QRivxBog QRivxOF LakeLayers P_LakeLayers IP_LakeLayers Bog_LakeLayers OF_LakeLayers tempstart taux #Run with states. The states from the previous times step is loaded 
 
end   

if (modstate == 0)                   # Run with initial values  
    smbog  = 0.0                          # Soilmoisture Bogs   
    QRivx[1:noDT] .= qsimutxP[1:noDT]             #discharge in m3/s, contribution from P
    QRivx[1:noDT] .= QRivx[1:noDT] .+ qsimutxIP[1:noDT] #adding contribution from IP
    QRivx[1:noDT] .= QRivx[1:noDT] .+ qsimutxBog[1:noDT]#adding contribution from bogs
    QRivxBog[1:noDT] .= qsimutxBog[1:noDT] 
    QRivxP[1:noDT] .= qsimutxP[1:noDT] 
    QRivxIP[1:noDT] .= qsimutxIP[1:noDT]
    totdef[1:Lty] .= M[1:Lty]
end

############################################################################
#                                                                          #
#simulation starts here.............. days is no. time steps               #
############################################################################
for i in startsim:days

  dato = Dates.DateTime(ptqinn.yr[i],ptqinn.mnt[i], ptqinn.day[i], ptqinn.hr[i],ptqinn.min[i]) #date
  #println(i," ",dato)      
         
  hr = ptqinn.hr[i] 
  DN = Dates.dayofyear(dato)           #daynumber        
  
  #Reads Precipitation and temperature for each elevation zone   
  htemp = ptqinn[i,16:25]
  hprecip = ptqinn[i,6:15] 

  meanprecip = mean(hprecip)
  meantemp =  mean(htemp)
  tempstart = TempstartUpdate(tempstart,htemp, len) #Updating the tempstart with this timesteps temperature 
 
  for Lst in 1:Lty     # landscape types, one snow regime for each landscape type P and IP. The other Lty have no snow

    for idim in 1:hson # elevation zones  

      #STempvec = TemperatureVector(ptqinn[(i-len+1):i,16:25],idim, len) # i is always greater than len
      STempvec = TemperatureVector(tempstart, idim, len)
      CGLAC = 0.0 # dummy, has no role in this version, energy balance is used
      CX = 0.0    # dummy, has no role in this version, energy balance is used
      TS = 0.0    # dummy, has no role in this version, energy balance is used

      #Estmates the mean snowpack temperature
     snittT[Lst] = SnowpackTemp(STempvec) 

      #Calculates Rain vs Snow and Snow/Glaciermelt
      WaterIn = NedbEB(CFR,DN,hprecip[idim],htemp[idim],spd[Lst,idim],Timeresinsec,TX,pkorr,skorr,
                         hr,u,Pa[idim],CGLAC,CX,TS,taux[Lst],snittT[Lst],phi,thi)

      #From WaterIn per elevation zone:         
      PS[Lst,idim] = WaterIn[1]        # Precipitation as snow
      PR[Lst,idim] = WaterIn[2]        # Precipitation as rain      
      MW[Lst,idim] = WaterIn[3]        # Snow melt potential
      MWGLAC[Lst,idim] = WaterIn[4]    # glaciers 
      SWrad[Lst,idim] = WaterIn[5]     # Short wave radiation (net)
      LA[Lst,idim] = WaterIn[6]        # Long wave radiation atomspheric
      LT[Lst,idim] = WaterIn[7]        # Long wave radiation terrestrial
      SH = WaterIn[8]                  # Turbulent heat
      LE = WaterIn[9]                  # Latent heat
      GH = WaterIn[10]                 # Ground heat
      PH = WaterIn[11]                 # Precipitation heat
      CC = WaterIn[12]                 # Cold content
      Albedo[Lst] = WaterIn[13]        # Albedo
      taux[Lst] = WaterIn[14]          # Aging factor of snow for albedo calulations
      snittT[Lst] = WaterIn[15]        # Snowpack temperature
      Tss = WaterIn[16]                # Snow surface temperature
      RH =  WaterIn[18]                # Relative humidity
      Cl = WaterIn[17]                 # CloudCover
      MWGLAC1[Lst,idim]=WaterIn[19]	   # Potential glacial melt from EB

      gisoil[Lst,idim] = MWGLAC1[1,idim]  # runoff due to glacial melt. is later corrected for snowfree areas
       
      # calculates Snow accumulation and Snowdistribution      
      FromSnow = SnowGamma(PR[Lst,idim],PS[Lst,idim],MW[Lst,idim],sca[Lst,idim],spd[Lst,idim],wcd[Lst,idim],
                        pro,nsno[Lst,idim],alfa[Lst,idim],ny[Lst,idim],a0[Lst],n0[Lst],a0[Lst],d[Lst])
    
      isoil[Lst,idim] = FromSnow[1]  # precipitation and snowmelt
      spd[Lst,idim] = FromSnow[2]    # snow water equivalent
      wcd[Lst,idim] = FromSnow[3]    # liquid water in snow
      sca[Lst,idim] = FromSnow[4] 
      nsno[Lst,idim] = FromSnow[5] 
      alfa[Lst,idim] = FromSnow[6] 
      ny[Lst,idim] = FromSnow[7]
     
      if (spd[idim] > 10000)     
	    if(kal == 1)
          spd[idim]  = 8000.0
	      skorr= 0.9*skorr
          println("Skorr is reduced, you build snowtowers = trend in SWE")
	    end 
      end
      
      if (sca[1,idim] < gca[idim])   # only for Lst = 1 permeable aras changed 15.02 2024
        snowfree[1, idim] = 1.0
      else
        snowfree[1,idim] = 0.0       # test for glaciermelt,  no melt if snowcovered, snowfree[idim] =0.0  
      end

     #Snowdepths and snowdenities
     snowdepthx = snowdepth[Lst,idim]
     snomag[Lst] =  sca[Lst, idim]*spd[Lst,idim] #snowreservoir this timestep

     #Calling snowdepth routine
     snowdepth[Lst,idim],sndens[Lst,idim] = NewSnowSDEB(PR[Lst,idim],PS[Lst,idim],htemp[idim],TX,
            MW[Lst,idim],snomag[Lst],snomag[Lst],snowdepthx)
     
    end #for elevation zones
 
    MSWrad[Lst] = mean(SWrad[Lst,1:hson])
    MLA[Lst] = mean(LA[Lst,1:hson])
    MLT[Lst] = mean(LT[Lst,1:hson])
    GIsoil[Lst] = mean(gisoil[Lst,1:hson])

    #Snowreservoir
    snomag[Lst] = mean(sca[Lst,1:hson] .* spd[Lst,1:hson]) #mean catchment SWE ,must multiply with SCA to have arealvalues
    swe_h[Lst,1:hson] = sca[Lst,1:hson] .* spd[Lst,1:hson] #SWE pr. elevation zone
    middelsca[Lst] = sum(sca[Lst,1:hson])/hson
    snofritt[Lst] = 1-middelsca[Lst]
    wcs[Lst] = mean(wcd[Lst,1:hson])

    #Estimate isoil from all elevation zones
    misoil[Lst] = sum(isoil[Lst,1:hson].*swgt[Lst,1:hson])                                      #snowmelt and rain on soils. 
    # isoil is weighted by the fraction of soils
    # in relation to soil an bog area pr elevation band
  end #For landscape types snow and rain
    
    # The following two statements soly for Lst ==1 (P areas)   
    m_r_onglac = sum(isoil[1,1:hson].*gwgt)            # snowmelt and rain from glaciated area. 
    # isoil is weighted by the fraction of glaciers 
    # pr elevation band in relation to total glaciated area  
    # glacier melt (gisoil) in mm this timestep. 
    # m_r_onglac + outglac is total output from glacier
    outglac = sum(gisoil[1,1:hson].*gwgt[1:hson].*snowfree[1,1:hson])  # gwgt because it is going to be scaled by glacfrac later on

    if(area[1] > 0.0) # area[1] includes area[5]
      ddistx[1,1:NoL] = LayerCapacityUpdate(LayersP, nodaysvector[1,1:NoL], Magkap[1,1:NoL], NoL)
     end
     if(area[2] > 0.0)
      ddistx[2,1:NoL] = LayerCapacityUpdate(LayersIP, nodaysvector[2,1:NoL], Magkap[2,1:NoL], NoL)
     end 
   
     gmltosoil[1] = misoil[1]*(1-(Ltyfrac[5]))+Ltyfrac[5]*(outglac+m_r_onglac) # input from rain, snow and glaciers always > 0 Needs some work!! [5] is glaciers
     gmltosoil[2] = sum(isoil[2,1:hson] .* swgt[2,1:hson]) # same as misoil[2]

    #Updating the deficit (for all sub surface layers, NOT overland flow layer) before Evapotranspiration
    for Lst in 1:Lty
      totdef[Lst] = sum(ddistx[Lst,2:NoL])            # This may be negative if input outstrips deficits
      
      PotEvap[Lst] = PotentialEvapPT(meantemp, MSWrad[Lst], MLA[Lst], MLT[Lst], MPa)
      if(PotEvap[Lst] < 0.0)
        PotEvap[Lst] = 0.0
      end
      Def[Lst] = max(totdef[Lst],0)
    end
	
	# Assigning overland flow if infiltration capapcity is exceeded and reduces tosoil in the process
	for Lst in 1:Lty
      OF[Lst],tosoil[Lst] = OFICap(gmltosoil[Lst],ICap[Lst])
	  #if (meanprecip > 0.0)
	   # println("precip=", meanprecip)
       # println("OF=",OF)
	   # println("tosoil=", tosoil)
	   # println(i," ",dato)
	  #end
    end 	  
   
    #calculating evapotranspiration from soilmoisture reservoir
    for Lst in 1:Lty
     sm[Lst], ea_sm[Lst], ea_S[Lst] = UnsaturatedEvapEB(tosoil[Lst], meantemp, sm[Lst], M[Lst], Def[Lst], PotEvap[Lst])
                  # updated sm
                  # evapotranspiration to be drawn from the soilmoisture
                  # evapotranspiration to be drawn from Layers  
    end
    
    #calculating additional (ea_S) evapotranspiration directly from Layers
    if(area[1] >0.0)
      LayersP, ea_S[1] = LayerEvap(LayersP, nodaysvector[1,1:NoL], ea_S[1],layerUH_P,NoL)
      ea[1] = ea_sm[1] + ea_S[1]
     end
     if(area[2] > 0.0)
      LayersIP, ea_S[2] = LayerEvap(LayersIP, nodaysvector[2,1:NoL], ea_S[2],layerUH_IP, NoL)
      ea[2] = ea_sm[2] + ea_S[2] # adding evapotranspiration from sm (first) and Layers (second)
     end 
       
    #Estimating current capacity in Layers after Evapotranspiration
    if(area[1] >0.0)
      ddistx[1,1:NoL] = LayerCapacityUpdate(LayersP, nodaysvector[1,1:NoL], Magkap[1,1:NoL], NoL)
     end
     if(area[2] > 0.0)
      ddistx[2,1:NoL] = LayerCapacityUpdate(LayersIP, nodaysvector[2,1:NoL], Magkap[2,1:NoL], NoL)
     end 

    #Updating the deficit (for all sub surface layers, NOT overland flow layer)
    for Lst in 1:Lty       
      totdef[Lst] = sum(ddistx[Lst,2:NoL])            # This may be negative if input outstrips deficits
      Def[Lst] = max(totdef[Lst],0)
    end 
     
  #Estimate isoil from all elevation zones and Lty
   for Lst in 1:Lty   
    # State of unsaturated zone, it can etiher be zero (complete saturation) or positive, the deficit. If negative we have overland flow 
    smlast = sm[Lst]
    # Call for the soilwater routine  calculates soil water (unsaturated) and evapotransipration 
    outx[Lst], sm[Lst] = UnsaturatedExEvap(tosoil[Lst], sm[Lst], R, Def[Lst])
    #if(round((tosoil[Lst]),digits=3)-round((outx[Lst]+(sm[Lst]-smlast)),digits = 3)!= 0.0)  # balance OK 20.06.2023     
    #  println(i,"hi, tosoil and sm not in balance", Lst)
    #  println(tosoil[Lst])
    # println(outx[Lst]+(sm[Lst]-smlast))
    #end   
   end
    
  #From Wetlands/Bogs
   smboglast = smbog
   Fromwetlands = WetlandsEB(tosoil[1], meantemp, middelsca[1], smbog, M[1], PotEvap[1])# takes on Permeable values, Lst =1
   outbog = Fromwetlands[1]
   smbog = Fromwetlands[2]
   eabog = Fromwetlands[3]

   #if(round((tosoil[1]),digits= 3)-round((outbog+(smbog+eabog-smboglast)),digits=3)!= 0.0)  # balance so far   
   #  println(i," hei after WetlandsEB")
   #  println(tosoil[1])
   #  println(outbog+(smbog+eabog-smboglast))
   #  println(sum(BogLayers))  #Boglayers is the routing of runoff generated from bogs, NOT the reservoir
   # end
  
  #Estimating input to Layers according to capacity
  for Lst in 1:Lty  
   ddist[Lst,1:NoL] = GrvInputDistributionICap(outx[Lst], NoL, ddistx[Lst,1:NoL], ddist[Lst,1:NoL], ICap[Lst])
  end
  
  #We need to include OF due to exceeeded infiltration capacity and UPDATE dist and outx. ddist will distribute outx taking into account both sat and inf excess
  for Lst in 1:Lty  
   if(OF[Lst]>0.0) # recall OF is only due to inf excess
        ddist[Lst,1] = ddist[Lst,1] + (OF[Lst]/(OF[Lst]+outx[Lst])) # will include water from exceedance of infiltration capacity
		ddist[Lst,2:NoL] = (1-ddist[Lst,1]).*ddist[Lst,2:NoL]  #updating the rest of ddist if OF
	    outx[Lst]= outx[Lst]+ OF[Lst]
   end
  end
  
  #Calculating UH for overland flow Layer #1
  OFDDD = 1
  if (OFDDD ==1)  
    for Lst in 1:Lty
     if(ddist[Lst,1]*outx[Lst] > 0)# overland flow for LST =1 (P) is confirmed for this event
      
      if(Lst==1)
        layerUH_P[1, 1:nodaysvector[1,1]] = OverlandFlowDynamicDD(k[Lst,1:NoL],ddist[Lst,1:NoL],outx[Lst], layerUH_P,
                    nodaysvector[Lst,1:NoL],NoL, Ltymid[Lst], CritFlux[Lst], Timeresinsec)
       end              
       if(Lst==2)
         if(area[2] > 0.0)
           layerUH_IP[1, 1:nodaysvector[2,1]] = OverlandFlowDynamicDD(k[Lst,],ddist[Lst,],outx[Lst], layerUH_IP,
                    nodaysvector[Lst,1:NoL],NoL, Ltymid[Lst], CritFlux[Lst], Timeresinsec)
         end
       end
     end
    end
  end

  ##Waterbalance calculations
  #RinnP =sum(QRivxP[2:noDT])*(Timeresinsec*1000/area[1]) # P water stored in RN from last TimeSep
  #RinnIP= sum(QRivxIP[2:noDT])*(Timeresinsec*1000/area[2]) # IP water stored in RN from last TimeStep
        
  #SPinn = sum(LayersP) + sum(ddist[1,1:NoL] .* outx[1] .* layerUH_P[1:NoL,1:antHorlag[1]])#Alle lag i P pluss dagens input før runoff
  ##SPinn1 = sum(LayersP[1,1:nodaysvector[1,1]]) + sum(ddist[1,1] .* outx[1] .* layerUH_P[1,1:nodaysvector[1,1]])#lag1
  ##SPinn2 = sum(LayersP[2,1:nodaysvector[1,2]]) + sum(ddist[1,2] .* outx[1] .* layerUH_P[2,1:nodaysvector[1,2]])#lag2 
  #SIPinn = sum(LayersIP) + sum(ddist[2,1:NoL] .* outx[2] .* layerUH_IP[1:NoL,1:antHorlag[2]])#Alle lag i IP pluss dagens input før runoff
  ##SIPinn1 = sum(LayersIP[1,1:nodaysvector[2,1]]) + sum(ddist[2,1] .* outx[2] .* layerUH_IP[1,1:nodaysvector[2,1]])
  ##SIPinn2  = sum(LayersIP[2,1:nodaysvector[2,2]]) + sum(ddist[2,2] .* outx[2] .* layerUH_IP[2,1:nodaysvector[2,2]])
   
  #GDT_P = sum(LayersP[1:NoL,1]) + sum(ddist[1,1:NoL] .* outx[1] .* layerUH_P[1:NoL,1])      #groundwater to be discharged into the rivernetwork + this timesteps contribution
  #GDT_IP = sum(LayersIP[1:NoL,1]) + sum(ddist[2,1:NoL] .* outx[2] .* layerUH_IP[1:NoL,1])   #groundwater to be discharged into the rivernetwork 
  if(area[1] >0.0)
    GDT_P = sum(LayersP[1:NoL,1]) + sum(ddist[1,1:NoL] .* outx[1] .* layerUH_P[1:NoL,1])      #groundwater to be discharged into the rivernetwork + this timesteps contribution
  end
  if(area[2] > 0.0)
    GDT_IP = sum(LayersIP[1:NoL,1]) + sum(ddist[2,1:NoL] .* outx[2] .* layerUH_IP[1:NoL,1])   #groundwater to be discharged into the rivernetwork 
  end 
  if(area[3] > 0.0)
   GDT_Bog =  BogLayers[1]+outbog*UHbog[1]         #bogwater to be discharged into the rivernetwork + this timesteps contribution
  end
  
  #Overland flow
  GDT_OF = Ltyfrac[1]*(LayersP[1,1]+(ddist[1,1]*outx[1]*layerUH_P[1,1]))
  if(area[2] > 0.0)
    GDT_OF  = GDT_OF +Ltyfrac[2]*(LayersIP[1,1]+(ddist[2,1]*outx[2]*layerUH_IP[1,1]))
  end  
  if(area[3] > 0.0)
    GDT_OF  = GDT_OF + Ltyfrac[3]*(BogLayers[1]+outbog*UHbog[1])# Overland flow part
  end         
            
  #Updating the saturation Layers
    for Lst in 1:Lty
      if(Lst==1)
         LayersP = LayerUpdate(ddist[Lst,1:NoL],outx[Lst], LayersP, layerUH_P, nodaysvector[Lst,1:NoL], NoL)
    end       
      if(Lst==2)
        LayersIP = LayerUpdate(ddist[Lst,1:NoL],outx[Lst], LayersIP, layerUH_IP, nodaysvector[Lst,1:NoL], NoL)       
      end
    end  

  BogLayers = BogLayerUpdate(outbog, BogLayers, UHbog, antBogsteps)#
   
  #summing up groundwater states
    for Lst in 1:Lty
      if(Lst == 1)
         lyrs[Lst] = sum(LayersP)                           # Todays precip is included, but ex todays runoff
         subsurface[Lst] = sum(LayersP[2:NoL,1:antHorlag[1]]) # gives the sum of layers after todays runoff has taken place
      end
      if(Lst==2)
         lyrs[Lst] = sum(LayersIP)
         subsurface[Lst] = sum(LayersIP[2:NoL,1:antHorlag[2]])
      end
    end
  
  #summing up bogstates
  boglyrs = sum(BogLayers)
        
  #reinstate original overland flow layer          
  for Lst in 1:Lty 
    if(Lst==1) 
     layerUH_P[1,1:nodaysvector[Lst,1]] .= SingleUH(k[Lst,1], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
    end
    if(Lst==2) 
     layerUH_IP[1,1:nodaysvector[Lst,1]] .= SingleUH(k[Lst,1], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
    end
  end

  GrWP = GrWPoint(nodaysvector[1,1:NoL], i, NoL, AreasP, LayersP, Parea, mindist, delta_dP) # returns a vector of mm groundwater at different points perpendicular to the RN
  
  #runoff this timestep
  if(area[1] > 0.0)
   qsimutxP .= (((GDT_P/1000)*area[1])/Timeresinsec) .* UHriver     #[m3/s]
  end
  if(area[2] > 0.0)
   qsimutxIP .= (((GDT_IP/1000)*area[2])/Timeresinsec) .* UHriver   #[m3/s]
  end
  if(area[3] > 0.0)
   qsimutxBog .= (((GDT_Bog/1000)*area[3])/Timeresinsec) .* UHriver #[m3/s]
  end
  qsimutxOF .= (((GDT_OF/1000)*totarea)/Timeresinsec) .* UHriver   #[m3/s]  # Overland flow. NOT a contribution, just extracted the fraction OF of total runoff
  
  #Total response
  qsimutx[1:noDT] .= 0.0  
  if(area[1] > 0.0)                                 # this is written anew for each timestep
   qsimutx[1:noDT] .= qsimutx[1:noDT] .+ qsimutxP[1:noDT]   #adding contribution from permeable areas
  end
  if(area[2] > 0.0)
   qsimutx[1:noDT] .= qsimutx[1:noDT] .+ qsimutxIP[1:noDT]  #adding contribution from impermeable areas
  end
  if(area[3] > 0.0)
   qsimutx[1:noDT] .= qsimutx[1:noDT] .+ qsimutxBog[1:noDT] #adding contribution from bogs
  end

  #Routing the contributions in the RN
  #Updating the routing in river network total)
  QRivx, QRD = RiverUpdate(noDT,QRivx,qsimutx)
  qmm_state = (sum(QRivx)-QRD)*(Timeresinsec*1000/totarea)       # This is also a reservoir [mm], water from todays event is stored for future runoff in the RN, relevant for WB.
                                                                  # Minus QRD since this will be handled in the WB equations
  #Updating the routing in river network Bogs
  QRivxBog,QRDBog = RiverUpdate(noDT,QRivxBog,qsimutxBog) 
    
  #Updating the routing in river network P
  QRivxP,QRDP = RiverUpdate(noDT,QRivxP,qsimutxP)  # OK OK 5.2.2024
    
  #Updating the routing in river network IP
  QRivxIP,QRDIP = RiverUpdate(noDT,QRivxIP,qsimutxIP) # OK 5.2.2024
  
  #Updating the routing in river network OF
  QRivxOF,QRDOF = RiverUpdate(noDT,QRivxOF,qsimutxOF)

  RnP = RiverPoint(noDT, QRivx, rndist, Ltymax[6]) # returns a vector of m3/s water in the RN at different distance from the outlet
  
  #Qmm = (QRD*Timeresinsec*1000/totarea)  #QRD in mm/day
  #Routing the contributions in the Lake
  Qm3s = LakeLayers[1] + QRD*UHLake[1]            #Lakewater to be discharged from outlet + this timesteps contribution, i.e. catchment response in m3/s
  LakeLayers = BogLayerUpdate(QRD, LakeLayers, UHLake, nodaysLake) # Same routine updating as for Bogs
  Qmm = (Qm3s*Timeresinsec*1000/totarea)  #GDT_Lake in mm/Timeresinsec

  P_Qm3s = P_LakeLayers[1] + QRDP*UHLake[1]            #Lakewater to be discharged from outlet + this timesteps contribution, i.e. catchemnt response 
  P_LakeLayers = BogLayerUpdate(QRDP, P_LakeLayers, UHLake, nodaysLake) # Same routine updating as for Bogs
  
  IP_Qm3s = IP_LakeLayers[1] + QRDIP*UHLake[1]            #Lakewater to be discharged from outlet + this timesteps contribution, i.e. catchemnt response 
  IP_LakeLayers = BogLayerUpdate(QRDIP, IP_LakeLayers, UHLake, nodaysLake) # Same routine updating as for Bogs
  
  Bog_Qm3s = Bog_LakeLayers[1] + QRDBog*UHLake[1]            #Lakewater to be discharged from outlet + this timesteps contribution, i.e. catchemnt response 
  Bog_LakeLayers = BogLayerUpdate(QRDBog, Bog_LakeLayers, UHLake, nodaysLake) # Same routine updating as for Bogs
  
  OF_Qm3s = OF_LakeLayers[1] + QRDOF*UHLake[1]            #not a contribution, just sttratifying the runoff 
  OF_LakeLayers = BogLayerUpdate(QRDOF, OF_LakeLayers, UHLake, nodaysLake) # Same routine updating as for Bogs
   
   #WBP = (SPinn + RinnP) - (lyrs[1] + sum(QRivxP)*(Timeresinsec*1000/area[1]))    #mm  last term inludes discharge and storage in rivernetwork                                         
   #println("WBP = ", WBP, " Total inn P= ",SPinn)                                             
   #WBIP = (SIPinn + RinnIP) - (lyrs[2] + sum(QRivxIP)*(Timeresinsec*1000/area[2]))   #mm  last term inludes discharge and storage in rivernetwork                                         
   #println("WBIP = ", WBIP," Total inn IP= ",SIPinn)  
##  Water balance before and after runoff calculations ok 20.06.2023. Also check that precip-(evap + sm + outx) == 0
   
  if (kal==0)
  #Assigning outdata to vector, one vector for each timestep
   simresult[i, 1:5] = [ptqinn.yr[i],ptqinn.mnt[i],ptqinn.day[i],ptqinn.hr[i],ptqinn.min[i]]
   simresult[i, 6:9] = [round(meanprecip,digits= 3),round(meantemp,digits = 3), ptqinn.q[i],round(Qm3s,digits=6)]
   simresult[i, 10:13] = [round(P_Qm3s,digits = 6),round(IP_Qm3s,digits = 6), round(Bog_Qm3s,digits = 6),round(middelsca[1],digits=3)]
   simresult[i, 14:16] = [round(snomag[1],digits = 6), round(lyrs[1],digits = 6),round(lyrs[2],digits = 6)]
   simresult[i, 17:20] = [round(totdef[1],digits=2), round(totdef[2],digits=2),round(sm[1],digits=6),round(sm[2],digits=6)]
   simresult[i, 21:24] = [round(ea[1],digits=6),round(ea[2],digits=6),round(Qmm,digits=6),round(smbog,digits=6)]
   simresult[i, 25:28] = [round(eabog,digits=6),round(qmm_state,digits=6), round(boglyrs,digits=6), round(middelsca[2],digits=3)]
   simresult[i, 29:32] = [round(snomag[2],digits= 6), round(wcs[1],digits=6),round(wcs[2],digits=6),round(subsurface[1],digits=6)]
   simresult[i, 33:34] = [round(subsurface[2],digits=6),round(OF_Qm3s,digits=6)]
   simresult[i, 35:36] = [round(outglac,digits=2), round(m_r_onglac,digits=2)]
   simresult[i, 37:39] = [round(GIsoil[1], digits=4),round(misoil[1], digits=4),round(snittT[1], digits=4)]
  end
  #grwpointresult[i, 1:29] <- c(ptqinn$yr[i],ptqinn$mnt[i],ptqinn$day[i],ptqinn$hr[i], round(GrWPoint[1:25],2))
                                
 
  push!(grwpointresult,(dato, GrWP[1],GrWP[2],GrWP[3],GrWP[4],GrWP[5],GrWP[6],GrWP[7],GrWP[8],GrWP[9],GrWP[10],GrWP[11],GrWP[12],GrWP[13],
  GrWP[14],GrWP[15],GrWP[16],GrWP[17],GrWP[18],GrWP[19],GrWP[20],GrWP[21],GrWP[22],GrWP[23],GrWP[24],GrWP[25]))  #writing groundwater soilmoisure to a DF [mm]
  
  push!(rnpointresult,(dato, RnP[1],RnP[2],RnP[3],RnP[4],RnP[5],RnP[6],RnP[7],RnP[8],RnP[9],RnP[10],RnP[11],RnP[12],RnP[13],
  RnP[14],RnP[15],RnP[16],RnP[17],RnP[18],RnP[19],RnP[20],RnP[21],RnP[22],RnP[23],RnP[24],RnP[25]))  #writing water in RN to a DF [m3/s]
  
  qberegn[i] = Qm3s      #Simulated runoff for current timestep. Used for skillscores etc.

  

#---------------------------------------------------------------------------------------------------

# Updating the  temperature matrix for estimating snowpack temperature
 tempstart = ptqinn[(i-len+1):i,16:25]  # assigning temperature values, i is always larger than len

# Saving states, SWE, sm, Layers, etc in a JLD2 file
if(i == (38165 + len) && savestate == 1) #  This number is ONE timestep less than startsim, i.e. the statefile is for the timestep before startsim 
  timestamp = Dates.format(dato, "yyyymmdd_HHMM")   
 
    statefile = string("\\\\nve.no\\fil\\h\\HM\\Interne Prosjekter\\HydSimOveralt\\Totalmodell\\utdata\\",
          "statefile_",catchment,"_",timestamp,"_.hdf5") 
    @save statefile sca spd wcd nsno alfa ny sm smbog LayersP LayersIP BogLayers QRivx QRivxP QRivxIP QRivxBog QRivxOF LakeLayers P_LakeLayers IP_LakeLayers Bog_LakeLayers OF_LakeLayers tempstart taux       
end
  #---------------------------------------------------------------------------------------------------

end # end of loop for number of timesteps in timeseries
 
println("nodaysLake=",nodaysLake)
println("nodaysRiver=",nodaysRiv)

 skillstart = Int(spinuptime*(86400/Timeresinsec))
 days2 = days
 
wwater = 1.0*persons*(140.0 + 42.0)/(86400.0*1000.0)# norsk vann equivalent use in L/day add 30% leakage from water net

 meansim = mean(wwater .+ qberegn[skillstart:days2])
 meanobs = mean(ptqinn.q[skillstart:days2])

# Computing skillscores NSE, kge, bias
 KGE, beta = (kge((wwater .+ qberegn[skillstart:days2]),ptqinn.q[skillstart:days2]))
 NSE = (nse((wwater .+ qberegn[skillstart:days2]),ptqinn.q[skillstart:days2]))
 bias = beta #(meansim/meanobs)
#bias = (meansim/meanobs) 

#Recall, tprm is the parameter vector
#                                       u,        pro,        TX,      PCritFlux,         OFP,        GshInt,   GscInt, persons
dfy = DataFrame(A=NSE,B=KGE,C=bias,D=tprm[1],E=tprm[2],F=tprm[3],G=tprm[4], H=tprm[5], II =tprm[6],IJ=tprm[7],
    JJ=tprm[8], KJ= tprm[9], LJ= tprm[10])
CSV.write(r2fil,DataFrame(dfy), delim = ';', header=false, append = true)


toptitles = ["Yr","Mnt","Day", "Hr","Min","Precip","Temp","Qobs","Qsim","Q_P","Q_IP","Q_Bogs","SCA_P","SWE_P","SS+_P", "SS+_IP","SSDef_P","SSDef_IP",
        "SM_P","SM_IP","Ea_P","Ea_IP","Qmm","SMBog","EaBog","qmm_state","Boglyrs","SCA_IP", "SWE_IP","WCS_P", "WCS_IP","SS_P", "SS_IP",
        "Q_OF", "outglac", "r_sm_outglac", "gisoil", "misoil","snittT[1]"]

 if (kal==0)
  CSV.write(utfile,DataFrame(simresult,:auto), delim = ';', header=toptitles)                
  println("M= ", round(M[1],digits=2)," ",round(M[2],digits =2))
  println("k_P= ", round(k[1,1],digits=6)," ",round(k[1,2],digits=6)," ",round(k[1,3],digits=6)," ",round(k[1,4],digits=6)," ",round(k[1,5],digits=6))
  if (area[2] >0.0)
   println("k_IP= ", round(k[2,1],digits=6)," ",round(k[2,2],digits=6)," ",round(k[2,3],digits=6)," ",round(k[2,4],digits=6)," ",round(k[2,5],digits=6))
  end
  println("Mean(Qsim)= ",meansim)
 end           

return ptqinn.q, qberegn, KGE,NSE,bias

end  #end of func
