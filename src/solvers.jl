#-------------------------------------------------------------------------
function terracenewton!(TerracePhysics::TerraceParameters,dt::Float64,
    sol::Array{Float64},h_sea::Float64,i1::Int64,sol_n,sol_old)
    #=
        Add derivation of newton-raphson
    =#
    P, h_wb, βz = TerracePhysics.P0, TerracePhysics.h_wb, TerracePhysics.βz
    cte1 = βz * P * dt
    cte2 = 4.0/h_wb
    cte3 = cte1 * cte2
    # cte4 = Array{Float64}(undef, Threads.nthreads())
    copyto!(sol_old, sol)
    residual = 1.0
    while residual > 1e-4
        copyto!(sol_n,sol)
        # @turbo for ii = i1:length(sol)
        #     cte4    = exp(-cte2*(h_sea-sol[ii]))
        #     sol[ii]-= (sol[ii] + cte1 * cte4 - sol_old[ii]) * (1/ (1.0 + cte3 * cte4))            
        # end       
        @tturbo for ii = i1:length(sol)
            sol[ii]-= (sol[ii] + cte1 * exp(-cte2*(h_sea-sol[ii])) - sol_old[ii]) * (1/ (1.0 + cte3 * exp(-cte2*(h_sea-sol[ii]))))            
        end       
        residual = norm2(sol,sol_n)
    end
end

#-------------------------------------------------------------------------
function terracenewton!(βz::Float64,P::Float64,h_wb::Float64,dt::Float64,
    sol::Array{Float64},h_sea::Float64,i1::Int64,sol_n,sol_old)
    #=
        Add derivation of newton-raphson
    =#
    cte1        = βz * P * dt
    cte2        = 4.0/h_wb
    cte3        = cte1 * cte2
    # sol_n      .= copy(sol) .+1e3   
    copyto!(sol_old,sol)
    residual    = 1.0
    while residual > 1e-4
        copyto!(sol_n,sol)
        @turbo for ii = i1:length(sol)
            cte4    = exp(-cte2*(h_sea-sol[ii]))
            sol[ii]-= (sol[ii] + cte1 * cte4 - sol_old[ii]) * (1/ (1.0 + cte3 * cte4))
        end       
        residual = norm2(sol,sol_n)
    end
end

# CLIFF RETREAT ================================================================
function cliffretreat!(Terrace,h_sea,P0,h_wb,βx,Δt)
    shelf_depth         = 10; # max depth of continental shelf
    bool1               = @. (abs(h_sea-Terrace.z)<shelf_depth)
    bool2               = @. h_sea>=Terrace.z
    id_shelf            = [bool1[kk] == true && bool2[kk] == true ? true : false for kk = 1:length(bool1)]            
    z_sea_terrace       = Terrace.z[id_shelf];
    h2sea_terrace       = @. abs(h_sea - z_sea_terrace);
    
    # -- Cliff retreat
    a                   = Terrace.x[id_shelf]
    b                   = @. P0*exp(-4*h2sea_terrace/h_wb)
    integral            = integrate(a,b,TrapezoidalEvenFast())
    # integral            = integrate1D(a,b)
    dx_cliff            = (P0 - integral) * Δt * βx 
    # -- Modify geometry
    x_terrace               = deepcopy(Terrace.x)
    id_cliff                = findlast(x->x>h_sea,Terrace.z)  # find clif position
    x_terrace[1:id_cliff] .-= abs(dx_cliff)
    # iz_correct              = id_cliff + 1
    idx                     = id_cliff-25:id_cliff+25
    Terrace.z[idx]         .= interpolate((x_terrace[idx],), Terrace.z[idx], Gridded(Linear()))(Terrace.x[idx])
    Terrace.z               = interpolate((x_terrace,), Terrace.z, Gridded(Linear()))(Terrace.x)
    Terrace.z[id_cliff]     = Terrace.z[id_cliff+1]

end
# ==============================================================================

@inline integrate1D(x,f) = sum(f*(x[2]-x[1])) 

@inline function fixterracex!(x,id)
    x0 = x[id]
    for i ∈ 1:id-1
        @inbounds if x[i] >= x0
            x[i] -= x0 - 0.01*i
        end 
    end    
end