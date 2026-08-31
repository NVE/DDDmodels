#Computes a single Normally distributed unit hydrograph (UH) based on the distance distributions and velocity
#Typically for Normal distance distributions, Rivers, Lakes

function SingleNormalUH(Lv,Timeresinsec,meanD,stdD,maxD) 

#using Distributions   

#meanD: scalar float = 129.84
#stdD: scalar float = 64.8
#maxD: scalar float = 538
#Timeresinsec: scalar Integer for example = 86400
#LV: scalar float, NB enters as an array float in other subroutines
    
nodaysLake = trunc(Int, maxD/Lv/Timeresinsec) + 1
if (nodaysLake > 1)     
  timeres = [0: 1: (nodaysLake-1);]                      # temporal resolution on riverrouting hydrogram, gives a sequance of number from 0 to 
  midLakescl = meanD/Lv/Timeresinsec #Actually mean response time (MRT) in days
  stdLakescl = stdD/Lv/Timeresinsec  #Actually standard deviation of response time (MRT) in days
  UHLake = zeros(nodaysLake)
  Lakenorm = Normal(midLakescl,stdLakescl) 
  UHLake = pdf.(Lakenorm,timeres)                      # makes  the pdf
  sum(UHLake)                                          # unclear purpose
  UHLake .= UHLake/sum(UHLake)                         # scale to give unity sum = 1
  noDT = Int(length(UHLake))
else
    UHLake = 1.0                                       # implies no Lake routing
    noDT = 1
end 

return UHLake, noDT, nodaysLake                                          #Unit hydrograph, i.e. weights distributing runoff in time
end

#nodaysRiv = Int(trunc(maxFl/rv/Timeresinsec)+1)
#if (nodaysRiv > 1)     
#  timeres = [0: 1: (nodaysRiv-1);]                      # temporal resolution on riverrouting hydrogram, gives a sequance of number from 0 to 
#  midFlscl = midFl/rv/Timeresinsec
#  stdFlscl = stdFl/rv/Timeresinsec
#  UHriver = zeros(nodaysRiv)
#  Rivnorm = Normal(midFlscl,stdFlscl) 
#  UHriver = pdf.(Rivnorm,timeres)                        # makes  the pdf
#  sum(UHriver)
#  UHriver .= UHriver/sum(UHriver)                         # scale to give unity sum = 1
#  noDT = Int(length(UHriver))
#else
#    UHriver = 1.0                      # implies no river routing
#    noDT = 1
#end
 
