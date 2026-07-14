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
    nrows = size(tempmatrix, 1)
    @inbounds for j in 1:len-1
        for i in 1:nrows
            tempmatrix[i,j] = tempmatrix[i,j+1] # move the level of the matrix one timestep ahead
        end
    end
    @inbounds for i in 1:nrows
        tempmatrix[i,len] = temps[i]
    end
end                
