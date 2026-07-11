#    Function tempstartUpdate
#
#------------------------------------------------------------------------
#     Description:  Updates the temperature matrix for snowpacktemperature estimation
#
#     Author: Thomas Skaugen
#     Revised: 19.02.2024
#--------------------------------------------------------------------------

function TempstartUpdate!(tempmatrix::Matrix{Float64}, temps::Vector{Float64}, len::Int)
# tempmatrix: timeseries of length 5 days 8length dep on temporal resolution) for 10 elevation zones
# temps: this timesteps temperatures
      @views tempmatrix[:,1:len-1] .= tempmatrix[:,2:len] # move the level of the matrix one timestep ahead
      tempmatrix[:,len] .= temps
end                
