abstract type SeaLevel end
abstract type Bintanja   end
abstract type Spratt450  end
abstract type Spratt800  end

struct SeaLevelOpts{T} 
    path::String
end

sea_level_curve(sealevel::SeaLevelOpts{Spratt450}) = string(sealevel.path,"Spratt2016-450kyr.txt");
sea_level_curve(sealevel::SeaLevelOpts{Spratt800}) = string(sealevel.path,"Spratt2016-800kyr.txt");
sea_level_curve(sealevel::SeaLevelOpts{Bintanja}) = string(sealevel.path,"Bintanja-3Ma.txt");

sea_level_starting_time(sea_level_file::SeaLevelOpts{Bintanja}) = 3e6
sea_level_starting_time(sea_level_file::SeaLevelOpts{Spratt800}) = 800e3
sea_level_starting_time(sea_level_file::SeaLevelOpts{Spratt450}) = 450e3

function get_sea_level(sealevelopts::SeaLevelOpts)

    sl = sea_level_curve(sealevelopts)
    curves0 = readdlm(sl)::Matrix{Float64}
    id = @. isnan(curves0[:,2]) != 1 # get rid of NaNs
    sealevel = curves0[id,:]
    nn = size(sealevel,1)
    age = Vector{Float64}(undef,nn)
    sea_level = Vector{Float64}(undef,nn)
    correction = sealevel[1,2]
    @simd for ii = 1:nn
        @inbounds sea_level[ii] = sealevel[ii,2] - correction; # Sea Level,,,meters above present day,,climate reconstructions,,Scaled first principal component of seven sea level reconstructions (0-430 ka),N
        @inbounds age[ii] = sealevel[ii,1]; # now in [yrs], originally in [kyrs]
    end
    
    return sea_level, age

end

function sea_level_corrections(sea_lvl_curve, sea_age, starting_time)
    sea_age = sea_age[end:-1:1] # reverse age array -> 0 = oldest age
    id_time = sea_age .>= maximum(sea_age) .- starting_time
    sea_age = sea_age[id_time]
    sea_lvl_curve = sea_lvl_curve[id_time]
    fsea = interpolate((reverse(sea_age),), reverse(sea_lvl_curve), Gridded(Linear()))
    return sea_lvl_curve, sea_age, fsea
end

function read_uplift()
    fname = "/home/albert/Dropbox/Riverini/DevelopNoTouch/uplift_rates/accre1_VertVel_870-871.txt"
    num = readdlm(fname)
    # correct time serie
    t = @. num[:,1] - num[1,1]; # [year]
    # mm/year to m/year
    uplift = @. num[:,2] / 1e3; # [mm/year] -> [m.year]

    return t, uplift

end