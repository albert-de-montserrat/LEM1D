abstract type Model end
abstract type TerraceOnly <:Model end
abstract type TerraceKnickpoint <:Model end
abstract type TerraceDiffusion <:Model end
abstract type TerraceKnickpointDiffusion <:Model end

struct ModelOpts{T} end
    
function solve(solver::Type{ModelOpts{TerraceOnly}},
    TerracePhysics, RiverPhysics, RiverArrays, HillSlopePhysics, SubmarinePhysics,
    Δt, River, Terrace, h_sea, id_shore, buffer1, buffer2,
    ScratchTerrace_FEM, Scratch_FEM, to)

    # -- TERRACE PROFILE
    zt = Terrace.z
    # @timeit to "newton terrace" terracenewton!(
    #     TerracePhysics,
    #     (Δt * yr),
    #     zt,
    #     h_sea,
    #     id_shore.terrace, 
    #     buffer1,
    #     buffer2)

    @timeit to "newton terrace" terracenewton2!(
        TerracePhysics,
        (Δt * yr),
        zt,
        h_sea,
        id_shore.terrace)

    return Terrace, River,to
end

function solve(solver::Type{ModelOpts{TerraceKnickpoint}},
    TerracePhysics, RiverPhysics, RiverArrays, HillSlopePhysics, SubmarinePhysics,
    Δt, River, Terrace, h_sea, id_shore, buffer1, buffer2,
    ScratchTerrace_FEM, Scratch_FEM, to)

    # -- RIVER PROFILE
    @timeit to "newton river" terracenewton!(
        TerracePhysics, 
        Δt * yr, 
        River.z, 
        h_sea, 
        id_shore.river + 1, 
        buffer1, 
        buffer2)

    # -- TERRACE PROFILE            
    @timeit to "newton terrace" terracenewton!(
        TerracePhysics,
        Δt * yr,
        Terrace.z,
        h_sea,
        id_shore.terrace, 
        buffer1,
        buffer2)            

    # KNICKPOINT MIGRATION ADVECTION -> FLUX LIMITER TVD 
    @timeit to "update r" updater!(RiverPhysics, RiverArrays, River)
    @timeit to "TVD" TVD!(River, buffer1, RiverPhysics, RiverArrays, Δt, id_shore.river)

    return Terrace, River,to
end

function solve(solver::Type{ModelOpts{TerraceKnickpointDiffusion}},
    TerracePhysics, RiverPhysics, RiverArrays, HillSlopePhysics, SubmarinePhysics,
    Δt, River, Terrace, h_sea, id_shore, buffer1, buffer2,
    ScratchTerrace_FEM, Scratch_FEM, to)

   # -- RIVER PROFILE
    @timeit to "newton river" terracenewton!(
        TerracePhysics, 
        Δt * yr, 
        River.z, 
        h_sea, 
        id_shore.river + 1, 
        buffer1, 
        buffer2)

    # -- TERRACE PROFILE            
    @timeit to "newton terrace" terracenewton!(
        TerracePhysics,
        Δt * yr,
        Terrace.z,
        h_sea,
        id_shore.terrace, 
        buffer1,
        buffer2)
                                               
    # KNICKPOINT MIGRATION ADVECTION -> FLUX LIMITER TVD
    @timeit to "update r" updater!(RiverPhysics, RiverArrays, River)
    @timeit to "TVD" TVD!(River, buffer1, RiverPhysics, RiverArrays, Δt, id_shore.river)

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
    Δt, River, Terrace, h_sea, id_shore, buffer1, buffer2,
    ScratchTerrace_FEM, Scratch_FEM, to)

    # -- TERRACE PROFILE            
    @timeit to "newton terrace" terracenewton!(
        TerracePhysics,
        (Δt * yr),
        Terrace.z,
        h_sea,
        (id_shore.terrace), 
        buffer1,
        buffer2)
    
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