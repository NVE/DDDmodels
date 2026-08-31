#    Function PyrAreas
#
#------------------------------------------------------------------------
#     Description: Estimates the areas asscociated with each layer at ecach time step due to celerity
#                  A postprocessing routine
#
#     Author:  Thomas Skaugen
#     Revised: 17.12.2019
# See Skaugen and Mengistu(2016)for the rationale behind this subroutine
# Need to calculate the ratios. 
# the parameter gamma of the exponential distance distribution is 1/midDL
# From Skaugen and Mengistu, 2016 we have that the ratio Kappa is related to 
# gamma= -log(kappa)/delta_d, where delta_d is the distance interval of choice
# The ratio kappa is then kappa= exp(-gamma*delta_d)=exp(-delta_d/midDL) 
#--------------------------------------------------------------------------


function PyrAreas(NoL,totarea,maxDl,nodaysvector, layerUH, antHorlag)

# antHorlag: 1 dim array Integer
# layerUH: 2 dim array float
# nodaysvector: 1 dim array integer
# NoL: scalar, integer    
# maxDl: scalar, integer
# totarea: scalar, float
   
  # total area distributed pr Layer
  Areas = zeros(antHorlag, NoL)                      # matrix(0.0, ncol=antHorlag,nrow=NoL)
  @inbounds for i in 1:NoL
    @simd for j in 1:nodaysvector[i]
       Areas[j,i] = totarea * layerUH[j,i]
    end
  end
  delta_d = maxDl ./ nodaysvector  # in meters (height as in the pyramid plots) pr. time-step box

#To be used in GRW_point subroutne: Layers__mm <- Layers*totarea/areas ######################
return Areas, delta_d
end

