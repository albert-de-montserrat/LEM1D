###############################################################
#  FUNCTION -> find_shore_id                                    #
###############################################################
function find_shore_id(h_sea::Float64,sol::Array{Float64},n::Int64)
    for ii = 1:n-1
         # @inbounds if sol[ii] > h_sea && sol[ii+1] <= h_sea
        @inbounds if sol[ii] <= h_sea
            return ii
            break
        end
    end
end ## END find_shore_id FUNCTION ##############################

###############################################################
#  FUNCTION -> find_shore_id                                  #
###############################################################
# function find_shore_id(h_sea::Float64,sol::CuArray{Float32},n::Int64)
#     numblocks = ceil(Int, n/256)
#     CUDA.@sync begin
#         @cuda threads=256 blocks=numblocks idx=CUfind_shore_id(h_sea,sol,n)
#     end
#     return idx 
# end

function CUfind_shore_id(h_sea,sol,n)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for ii = index:stride:n
        # @inbounds if sol[ii] > h_sea && sol[ii+1] <= h_sea
        @inbounds if sol[ii] <= h_sea
            return ii
            break
        end
    end
end ## END find_shore_id FUNCTION ##############################

###############################################################
#  FUNCTION -> localmaxmin                                    #
###############################################################
function localmaxmin(y::Array{Float64,1},nn::Int64)
    dummy_min = Array{Int32}(undef,nn,1)
    dummy_max = Array{Int32}(undef,nn,1)
    dummy_ids = Array{Int32}(undef,nn,1)
    fill!(dummy_min,0)
    fill!(dummy_max,0)
    fill!(dummy_ids,0)

    # ---- Check 1st node
    cmin,cmax,cids = 1,1,1
    if y[2] > y[1] # local minima
        dummy_min[1] = 1
        cmin        += 1
        dummy_ids[1] = 1
        cids        += 1
    else # local maxima
        dummy_max[1] = 1
        cmax        += 1
        dummy_ids[1] = 1
        cids        += 1
    end

    # ---- Check 2:end-1 node
    @inbounds @simd for ii = 2:nn-2
        # global cmax,cmin,cids
        if y[ii+1] > y[ii] && y[ii-1] > y[ii] # local minima
            dummy_min[cmin] = ii
            cmin    += 1
            dummy_ids[cids] = ii
            cids    += 1
        elseif y[ii+1] < y[ii] && y[ii-1] < y[ii] # local maxima
            dummy_max[cmax] = ii
            cmax    += 1
            dummy_ids[cids] = ii
            cids    += 1
        end
    end

    # ---- Check last node
    if y[nn] < y[nn-1] # local minima
        dummy_min[cmin] = nn
        cmin        += 1
        dummy_ids[cids] = nn
        cids        += 1
    else # local maxima
        dummy_max[cmax] = nn
        cmax        += 1
        dummy_ids[cids] = nn
        cids        += 1
    end

    # ---- Remove zeros
    idx_min = Array{Int32}(undef, cmin-1,1)
    idx_max = Array{Int32}(undef, cmax-1,1)
    idx_ids = Array{Int32}(undef, cids-1,1)
    @inbounds for ii = 1:cmin-1
        idx_min[ii] = dummy_min[ii]
    end
    @inbounds for ii = 1:cmax-1
        idx_max[ii] = dummy_max[ii]
    end
    @inbounds for ii = 1:cids-1
        idx_ids[ii] = dummy_ids[ii]
    end

    return idx_min,idx_max,idx_ids
end ## END localmaxmin FUNCTION ##############################

###############################################################
#  FUNCTION -> segdist                                        #
###############################################################
function segdist(idx_ids::Array{Int32,2},x::Vector{Float64},z::Vector{Float64})
    nx        = length(x)
    nsegments = length(idx_ids) - 1
    Lsegments = Array{Float64}(undef,nsegments,1)
    dist2high = Array{Float64}(undef,nx,1)

    @inbounds for ii = 1:nsegments
        id1 = idx_ids[ii];
        id2 = idx_ids[ii+1];

        for jj = id1:id2

            if  z[id2] > z[id1]
                xhigh = x[id2]
            else
                xhigh = x[id1]
            end

            dist2high[jj] = abs(x[jj] - xhigh)

        end
    end

    return dist2high

end ## END segdist FUNCTION ##################################
