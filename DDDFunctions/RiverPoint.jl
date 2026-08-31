#    Function RiverPoint
#
#------------------------------------------------------------------------
#     Description: Estimates water at a point in the river network
#
#     Author:  Thomas Skaugen
#     Revised: 05.08.2023
#--------------------------------------------------------------------------

function RiverPoint(vectorlengde,innQRivx, rnw_distvec, maxFl)

# Here we will adress the river UH at a certian distance from the outlet and extract the water amount 
# The RN UH is normally distributed implying that most water is at the center of the UH. This
# is at it should be since this "most water" indicates most number of river reaches and hence most erosion
# potential (most water exposed to river beds and -sides)
# delta_d is how wide are each box in the UH
# rnw_distvec the distance away from the outlet we want estimated

ant_rnw = length(rnw_distvec)       # number of points int he Rn we want estimated
rnw = zeros(ant_rnw) # vector for storing groundwater values for the current Layer 
rnwdelta_d = maxFl/vectorlengde

for i in 1:(ant_rnw-1) #to be certain we do not go beyond the matrix 
   j = trunc(Int, rnw_distvec[i]/rnwdelta_d) + 1
   rnw[i] = round(innQRivx[j], digits=3)
   #println("rnw: ",rnw[i])
end
             
return rnw

end
