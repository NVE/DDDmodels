"""
    nse(q_sim, q_obs)
Compute Nash-Sutcliffe efficiency
by Jan Magnusson. Revised by Anne Ellekjær Stavang and Thomas Skaugen
"""

function nse(q_sim, q_obs)

    ikeep = q_obs .> 0 # revidering
    q_obs = q_obs[ikeep]
    q_sim = q_sim[ikeep]

    if all(isnan, q_sim) || all(isnan, q_obs)

        nse = NaN
    
    else
        ikeep = .!isnan.(q_obs) .& .!isnan.(q_sim) 

        q_sim .= q_sim[ikeep]
        q_obs .= q_obs[ikeep]

        nse = 1 - sum((q_sim .- q_obs).^2) / sum((q_obs .- mean(q_obs)).^2)

    end

    return nse

end


