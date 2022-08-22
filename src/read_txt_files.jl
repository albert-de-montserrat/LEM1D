# START FUNCTION -> READ SEA LEVEL CURVE =======================================
function get_sea_level(sea_level_file::String)
    if sea_level_file == "Spratt2016-450"
        # LOAD ====================================================================
        curves0 = readdlm("Spratt2016-450ka.txt")
        # CORRECT-1 ===============================================================
        # get rid of NaN
        id = @. isnan(curves0[:, 2]) != 1
        sealevel = curves0[id, :]
        # DECLARE =================================================================
        nn = size(sealevel, 1)
        age = Array{Float64}(undef, nn, 1)
        sea_level = Array{Float64}(undef, nn, 1)
        correction = sealevel[1, 2]
        @simd for ii in 1:nn
            @inbounds age[ii] = sealevel[ii, 1] # already in [yrs]
            @inbounds sea_level[ii] = sealevel[ii, 2] - correction
        end
        # spratt450 columns: # Sea Level / meters above present day / climate reconstructions / Scaled first principal component of seven sea level reconstructions (0-430 ka),N

    elseif sea_level_file == "Spratt2016-800"
        # LOAD ====================================================================
        curves0 = readdlm("Spratt2016-800ka.txt")
        # CORRECT-1 ===============================================================
        # get rid of NaN
        id = @. isnan(curves0[:, 2]) != 1
        sealevel = curves0[id, :]
        # DECLARE =================================================================
        nn = size(sealevel, 1)
        age = Array{Float64}(undef, nn, 1)
        sea_level = Array{Float64}(undef, nn, 1)
        correction = sealevel[1, 2]
        @simd for ii in 1:nn
            @inbounds age[ii] = sealevel[ii, 1] # already in [yrs]
            @inbounds sea_level[ii] = sealevel[ii, 2] - correction
        end

    elseif sea_level_file == "Siddall2003-379"
        # LOAD ====================================================================
        curves0 = readdlm("Siddall2003-379ka.txt")
        # CORRECT-1 ===============================================================
        # get rid of NaN
        id = @. isnan(curves0[:, 2]) != 1
        sealevel = curves0[id, :]
        # DECLARE =================================================================
        nn = size(sealevel, 1)
        age = Array{Float64}(undef, nn, 1)
        sea_level = Array{Float64}(undef, nn, 1)
        correction = sealevel[1, 2]
        @simd for ii in 1:nn
            @inbounds age[ii] = sealevel[ii, 1] # already in [yrs]
            @inbounds sea_level[ii] = sealevel[ii, 2] - correction
        end

    elseif sea_level_file == "Bintanja3"
        # LOAD ====================================================================
        curves0 = readdlm("Bintanja-3Ma.txt")
        # CORRECT-1 ===============================================================
        # get rid of NaN
        id = @. isnan(curves0[:, 2]) != 1
        sealevel = curves0[id, :]
        # DECLARE =================================================================
        nn = size(sealevel, 1)
        age = Array{Float64}(undef, nn, 1)
        sea_level = Array{Float64}(undef, nn, 1)
        correction = sealevel[1, 2]
        @simd for ii in 1:nn
            @inbounds age[ii] = sealevel[ii, 1] # already in [yrs]
            @inbounds sea_level[ii] = sealevel[ii, 2] - correction
        end

    elseif sea_level_file == "Bintanja6"
        # LOAD ====================================================================
        curves0 = readdlm("Bintanja-6Ma.txt")
        # CORRECT-1 ===============================================================
        # get rid of NaN
        id = @. isnan(curves0[:, 2]) != 1
        sealevel = curves0[id, :]
        # DECLARE =================================================================
        nn = size(sealevel, 1)
        age = Array{Float64}(undef, nn, 1)
        sea_level = Array{Float64}(undef, nn, 1)
        correction = sealevel[1, 2]
        @simd for ii in 1:nn
            @inbounds age[ii] = sealevel[ii, 1] # already in [yrs]
            @inbounds sea_level[ii] = sealevel[ii, 2] - correction
        end

    elseif sea_level_file == "Bintanja3x2"
        # LOAD ====================================================================
        curves0 = readdlm("Bintanja-3x2.txt")
        # CORRECT-1 ===============================================================
        # get rid of NaN
        id = @. isnan(curves0[:, 2]) != 1
        sealevel = curves0[id, :]
        # DECLARE =================================================================
        nn = size(sealevel, 1)
        age = Array{Float64}(undef, nn, 1)
        sea_level = Array{Float64}(undef, nn, 1)
        correction = sealevel[1, 2]
        @simd for ii in 1:nn
            @inbounds age[ii] = sealevel[ii, 1] # already in [yrs]
            @inbounds sea_level[ii] = sealevel[ii, 2] - correction
        end
    end

    # OUTPUT
    return sea_level, age
end ## END READ SEA LEVEL FUNCTION

## START FUNCTION -> READ UPLIFT ===============================================
function read_uplift(fname)
    num = readdlm(fname)
    # -- correct time serie
    t = @. num[:, 1] - num[1, 1] # [year]
    # -- mm/year to m/year
    uplift = @. num[:, 2] / 1e3 # [mm/year] -> [m/year]

    return t, uplift
end ## END READ SEA LEVEL FUNCTION
