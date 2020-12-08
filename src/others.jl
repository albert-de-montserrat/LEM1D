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
end ## END localmaxmin FUNCTION

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

end ## END segdist FUNCTION

function firstderivative(x,z) # 1st derivative
    dz = (z[2] - z[1])*2
    ∇z = similar(x)

    @avx for i ∈ 2:length(x)-1 
        ∇z[i] = (x[i+1] - x[i-1]) / dz
    end

    ∇z[1]   = (x[2] - x[1]) / (dz)
    ∇z[end] = (x[end] - x[end-1]) / (dz)
    
    return ∇z
end

function secondderivative(x,z) # 2nd derivative
    dz = 2*(z[2] - z[1])^2
    Δz = similar(x)

    @avx for i ∈ 2:length(x)-1 
        Δz[i] = (x[i+1] - 2*x[i] + x[i-1]) / dz
    end

    Δz[1]   = Δz[2]
    Δz[end] = Δz[end-1]

    return Δz
end

function erodedvolume(z0,z,x)
    δz = z - z0
    return integrate(x,δz,TrapezoidalEvenFast())
end

function riverlength(x,z,nnod)        
    L = fill(0.0,nnod)        
    @inbounds @fastmath for i ∈ 2:nnod
        L[i] = L[i-1] + sqrt( (x[i]-x[i-1])^2  + (z[i]-z[i-1])^2) 
    end
    return L
end

# @inline function riverlength!(L,x,z,nnod)                   
#     @inbounds @fastmath for i ∈ 2:nnod
#         L[i] = L[i-1] + sqrt( (x[i]-x[i-1])^2  + (z[i]-z[i-1])^2) 
#     end
# end

# @inline function updater!(r, Kr, kappa_a, L, h, m, dx, x,z , nn)
#     riverlength!(L,x,z,nn)
#     r .= @. -Kr * (kappa_a^m) * (L^(h * m)) / dx 
# end

@inline function updater!(r, Kr, kappa_a, L, h, m, dx, x,z , nn)
    Threads.@threads for i ∈ 2:nn
        L[i] = L[i-1] + sqrt( dx  + (z[i]-z[i-1])*(z[i]-z[i-1])) 
        # r[i] = -Kr * (kappa_a^m) * (L[i]^(h * m)) / dx 
        r[i] = -Kr * pow32(kappa_a,m) * pow32(L[i],h * m) / dx 

    end
end

@inline function pow32(a::Float64,b::Number)
    a = convert(Float32,a)^convert(Float32,b)
    return convert(Float64,a)
end