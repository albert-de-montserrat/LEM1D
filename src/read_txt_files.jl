
# START FUNCTION -> READ SEA LEVEL CURVE =======================================
function get_sea_level(sea_level_file::String)

    if sea_level_file == "Spratt2016-450kyrs"
        # LOAD ====================================================================
        path                    = "/home/albert/Dropbox/Riverini/DevelopNoTouch/sea-level/";
        sl                      = string(path,"Spratt2016-MLfile.txt");
        curves0                 = readdlm(sl);
        # CORRECT-1 ===============================================================
        # get rid of NaN
        id                      = @. isnan(curves0[:,2]) != 1
        sealevel                = curves0[id,:]
        # DECLARE =================================================================
        nn                      = size(sealevel,1)
        age                     = Array{Float64}(undef,nn,1)
        sea_level               = Array{Float64}(undef,nn,1)
        correction              = sealevel[1,2]
        @simd for ii = 1:nn
            @inbounds age[ii]             = sealevel[ii,1] * 1e3; # now in [yrs], originally in [kyrs]
            @inbounds sea_level[ii]       = sealevel[ii,2] - correction; # Sea Level,,,meters above present day,,climate reconstructions,,Scaled first principal component of seven sea level reconstructions (0-430 ka),N
        end

        # age                     = @. curves[:,1] * 1e3; # now in [yrs], originally in [kyrs]
        # sea_level               = @. curves[:,2] - curves[1,2]; # Sea Level,,,meters above present day,,climate reconstructions,,Scaled first principal component of seven sea level reconstructions (0-430 ka),N
        # CORRECT-2 ===============================================================
        # t=0 -> h = 0
        # sea_level               = @. sea_level - sea_level[1];

    elseif  sea_level_file == "Bintanja"
        age         = 1
        sea_level   = 1

    elseif  sea_level_file == "Spratt2016-800kyrs"
         # LOAD ====================================================================
         path                    = "/home/albert/Dropbox/Riverini/DevelopNoTouch/sea-level/";
         sl                      = string(path,"Spratt2016-800kyr.txt");
         curves0                 = readdlm(sl);
         # CORRECT-1 ===============================================================
         # get rid of NaN
         id                      = @. isnan(curves0[:,2]) != 1
         sealevel                = curves0[id,:]
         # DECLARE =================================================================
         nn                      = size(sealevel,1)
         age                     = Array{Float64}(undef,nn,1)
         sea_level               = Array{Float64}(undef,nn,1)
         correction              = sealevel[1,2]
         @simd for ii = 1:nn
             @inbounds age[ii]             = sealevel[ii,1] ; # now in [yrs], originally in [kyrs]
             @inbounds sea_level[ii]       = sealevel[ii,2] - correction; # Sea Level,,,meters above present day,,climate reconstructions,,Scaled first principal component of seven sea level reconstructions (0-430 ka),N
         end

    end

    # OUTPUT
    return sea_level, age

end ## END READ SEA LEVEL FUNCTION

## START FUNCTION -> READ UPLIFT ===============================================
function read_uplift()
    num     = readdlm("/Users/albert/Dropbox/Riverini/DevelopNoTouch/uplift_rates/Topo_Vy_Forearc.txt") ;
    # -- correct time serie
    t       = @. num[:,1] - num[1,1]; # [year]
    # -- mm/year to m/year
    uplift  = @. num[:,3] / 1e3; # [mm/year] -> [m.year]

    return t, uplift

end ## END READ SEA LEVEL FUNCTION
