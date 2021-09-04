±(a::Number, b::Number) = (a-b, a+b)

function find_shore_id(h_sea, River::Profile, Terrace::Profile, TerracePhysics, id::NamedTuple)
    # unpack 
    zr, zt, n = River.z, Terrace.z, Terrace.nnod
    # indices
    ir, it = 0, 0
    midpoint = n÷2 
    mr, mt = id.terrace.i1, id.river.i1
    # We start searching from the center of the profile and start moving with ±1 steps to both sides
    for j = 1:midpoint

        idx = ifelse(
            j == 1,
            mr ± 0, # ± 0 to make it type stable
            mr ± j
        )

        if ir == 0
            for i in idx
                # check river
                @inbounds ir = ifelse(
                    zr[i] <= h_sea,
                    i,
                    0
                )
            end
        end
       
        idx = ifelse(
            j == 1,
            mt ± 0, # ± 0 to make it type stable
            mt ± j
        )

        if it == 0
            for i in idx
                # check terrace
                @inbounds it = ifelse(
                    zt[i] <= h_sea,
                    i,
                    0
                )
            end
        end

        if (ir>0) && (it>0)
            infected_nodes_terrace = infected_indices(
                Terrace, TerracePhysics, it
            )
            infected_nodes_river = infected_indices(
                River, TerracePhysics, it
            )
            return (terrace = infected_nodes_terrace, river = infected_nodes_river)
        end

    end


end

function infected_indices(
    profile::Profile, TerracePhysics::TerraceParameters, id_shore::Int64
    )
    
    idx, i1 = id_shore, id_shore
    z = profile.z
    h_wb = 3/2*TerracePhysics.h_wb # influence of sea bed erosion (x1.5 as a safety buffer zone)
    zshore = z[idx] # height of sea level
    zᵢ = z[idx+1]

    while (zshore - zᵢ) < h_wb
        idx += 1
        @inbounds zᵢ = z[idx]
    end

    (i1 = i1, i2 = idx)
    
end

function find_shore_id(h_sea::Float64,sol::Array{Float64},n::Int64)
    for ii = 1:n-1
         # @inbounds if sol[ii] > h_sea && sol[ii+1] <= h_sea
        @inbounds if sol[ii] <= h_sea
            return ii
            break
        end
    end
end

function find_shore_id_terrace(h_sea::Float64,sol::Array{Float64},n::Int64)
    i1=1
    for ii = 1:n-1
         # @inbounds if sol[ii] > h_sea && sol[ii+1] <= h_sea
        @inbounds if sol[ii] <= h_sea
            i1 =  ii
        end
        @inbounds if h_sea-sol[ii] > 500
            return i1, ii
        end
    end
end