abstract type Model end
abstract type TerraceOnly <:Model end
abstract type TerraceKnickpoint <:Model end
abstract type TerraceDiffusion <:Model end
abstract type TerraceKnickpointDiffusion <:Model end

struct ModelOpts{T} end
    
function solve(solver::Type{ModelOpts{TerraceOnly}},
    TerracePhysics, RiverPhysics, RiverArrays, HillSlopePhysics, SubmarinePhysics,
    Δt, River, Terrace, h_sea, id_shore, buffer,
    ScratchTerrace_FEM, Scratch_FEM, to)

    @timeit to "newton terrace" terracenewton2!(
        TerracePhysics,
        (Δt * yr),
        Terrace.z,
        h_sea,
        id_shore.terrace)

    return Terrace, River,to
end

function solve(solver::Type{ModelOpts{TerraceKnickpoint}},
    TerracePhysics, RiverPhysics, RiverArrays, HillSlopePhysics, SubmarinePhysics,
    Δt, River, Terrace, h_sea, id_shore, buffer,
    ScratchTerrace_FEM, Scratch_FEM, to)

    # RIVER PROFILE
    @timeit to "newton river" terracenewton2!(
        TerracePhysics, 
        Δt * yr, 
        River.z, 
        h_sea, 
        id_shore.river, 
        )

    # TERRACE PROFILE            
    @timeit to "newton terrace" terracenewton2!(
        TerracePhysics,
        (Δt * yr),
        Terrace.z,
        h_sea,
        id_shore.terrace
        )     

    # KNICKPOINT MIGRATION ADVECTION -> FLUX LIMITER TVD 
    @timeit to "update r" updater!(RiverPhysics, RiverArrays, River)
    @timeit to "TVD" TVD!(River, buffer, RiverPhysics, RiverArrays, Δt, id_shore.river.i1)

    return Terrace, River,to
end

function solve(solver::Type{ModelOpts{TerraceKnickpointDiffusion}},
    TerracePhysics, RiverPhysics, RiverArrays, HillSlopePhysics, SubmarinePhysics,
    Δt, River, Terrace, h_sea, id_shore, buffer,
    ScratchTerrace_FEM, Scratch_FEM, to)

    # -- RIVER PROFILE
    @timeit to "newton river" terracenewton2!(
        TerracePhysics, 
        Δt * yr, 
        River.z, 
        h_sea, 
        id_shore.river, 
    )

    # -- TERRACE PROFILE            
    @timeit to "newton terrace" terracenewton2!(
        TerracePhysics,
        (Δt * yr),
        Terrace.z,
        h_sea,
        id_shore.terrace
    )   
                                               
    # KNICKPOINT MIGRATION ADVECTION -> FLUX LIMITER TVD
    @timeit to "update r" updater!(RiverPhysics, RiverArrays, River)
    @timeit to "TVD" TVD!(River, buffer, RiverPhysics, RiverArrays, Δt, id_shore.river.i1)

    # Diffusion over Δt (FEM)
    @timeit to "FEM River" River.z = femsolver(Δt, Δx, nn, River, h_sea,
                                               Scratch_FEM,        
                                               HillSlopePhysics,
                                               SubmarinePhysics)

    # TERRACE 
    # Diffusion over Δt
    @timeit to "FEM Terrace" Terrace.z = femsolver(
            Δt, 
            Terrace, 
            h_sea,
            ScratchTerrace_FEM,
            Scratch_FEM,
            HillSlopePhysics
    )

    return Terrace, River,to                                               
end

function solve(solver::Type{ModelOpts{TerraceDiffusion}},
    TerracePhysics, RiverPhysics, RiverArrays, HillSlopePhysics, SubmarinePhysics,
    Δt, River, Terrace, h_sea, id_shore, buffer, 
    ScratchTerrace_FEM, Scratch_FEM, to)

    # -- RIVER PROFILE
    @timeit to "newton river" terracenewton2!(
        TerracePhysics, 
        Δt * yr, 
        River.z, 
        h_sea, 
        id_shore.river, 
    )
    
    # TERRACE
    @timeit to "FEM Terrace" Terrace.z = femsolver(
            Δt, 
            Terrace, 
            h_sea,
            ScratchTerrace_FEM,
            Scratch_FEM,
            HillSlopePhysics
    )
            
    return Terrace, River,to                                               
end

function old_solver!(Terrace, TerracePhysics, h_sea, id_shore, dt)
    P, h_wb, βz = TerracePhysics.P0, TerracePhysics.h_wb, TerracePhysics.βz
    for i in id_shore.terrace.i1: id_shore.terrace.i2
        @inbounds Terrace.z[i] += - βz .* P .* exp(-4*(h_sea-Terrace.z[i])/h_wb) * dt
    end
end

#-------------------------------------------------------------------------
function terracenewton!(TerracePhysics::TerraceParameters,dt::Float64,
    sol::Array{Float64},h_sea::Float64,infected_nodes::NamedTuple,sol_n,sol_old)
    #=
        Add derivation of newton-raphson
    =#
    P, h_wb, βz = TerracePhysics.P0, TerracePhysics.h_wb, TerracePhysics.βz
    cte1 = βz * P * dt
    cte2 = 4.0/h_wb
    cte3 = cte1 * cte2
    copyto!(sol_old, sol)
    residual = 1.0
    i1, i2 = infected_nodes.i1, infected_nodes.i2
    while residual > 1e-8
        copyto!(sol_n,sol)    
        @tturbo for ii in i1:i2
            sol[ii] -= (sol[ii] + cte1 * exp(-cte2*(h_sea-sol[ii])) - sol_old[ii]) * (1/ (1 + cte3 * exp(-cte2*(h_sea-sol[ii]))))
        end       
        residual = norm2(sol, sol_n, infected_nodes)
    end
end

function _newton(sol,h_sea, cte1, cte2)
    h0, h, hprevious = sol, sol, sol
    residual = 1.0
    its = 0
    while its < 21 && residual > 1e-8
        cte3 = cte1 * exp(-cte2*(h_sea-h))
        up = h - hprevious + cte3
        down = 1/(1 - cte3)
        h -= up * down
        residual = abs(h-h0)/abs(h0)
        h0 = h
        its += 1
    end
    h
end

function terracenewton2!(TerracePhysics::TerraceParameters,dt::Float64,
    sol::Array{Float64},h_sea::Float64,infected_nodes::NamedTuple)
    #=
        Add derivation of newton-raphson
    =#
    P, h_wb, βz = TerracePhysics.P0, TerracePhysics.h_wb, TerracePhysics.βz
    cte1 = βz * P * dt
    cte2 = 4.0/h_wb
    i1, i2 = infected_nodes.i1, infected_nodes.i2
    for ii in i1:i2
        @inbounds sol[ii] = _newton(sol[ii] ,h_sea, cte1, cte2)
    end
    
end

#-------------------------------------------------------------------------
function terracenewton!(βz::Float64,P::Float64,h_wb::Float64,dt::Float64,
    sol::Array{Float64},h_sea::Float64,i1::Int64,sol_n,sol_old)
    #=
        Add derivation of newton-raphson
    =#
    cte1 = βz * P * dt
    cte2 = 4.0/h_wb
    cte3 = cte1 * cte2
    copyto!(sol_old, sol)
    residual    = 1.0
    while residual > 1e-3
        copyto!(sol_n,sol)
        @turbo for ii = i1:length(sol)
            cte4    = exp(-cte2*(h_sea-sol[ii]))
            sol[ii]-= (sol[ii] + cte1 * cte4 - sol_old[ii]) * (1/ (1.0 + cte3 * cte4))
        end       
        residual = norm2(sol, sol_n)/norm(sol_n)
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
