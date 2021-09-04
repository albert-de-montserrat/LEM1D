import Pkg; Pkg.activate(".")
using LEM1D
using TimerOutputs, LoopVectorization

Δx = 30.0
Δt = 100
U0=0.0001
βz= 2e-6
sea_level_file = Bintanja
path_sea_level = "/home/albert/Dropbox/Riverini/Pluto/"
solver_opt = TerraceKnickpoint

function main(;
    Δx = 30.0, Δt = 100, U0=0.0001, βz= 2e-6, sea_level_file = Bintanja, 
    path_sea_level = "/home/albert/Dropbox/Riverini/Pluto/",
    solver_opt = TerraceOnly,
    slope = -1/50
    )

    uplift_type = "constant" # "constant" or "armel"
    solver = ModelOpts{solver_opt}
        # Options:
        #  - TerraceOnly
        #  - TerraceKnickpoint
        #  - TerraceDiffusion
        #  - TerraceKnickpointDiffusion
    sea_level_file  = sea_level_file 
        # no need to use  "" anymore. Options:
        #  - Spratt2016-450
        #  - Spratt2016-800
        #  - Bintanja
    sealevelopts = SeaLevelOpts{sea_level_file}(path_sea_level)
    # =========================================================================
    
    # INITIAL GEOMETRY ========================================================
    slope = slope 
    # slope = -tand(5)#1/50
    elementype = "linear"
    x, z0, _ ,nn, nel, e2n = mesher(Δx,slope,elementype)
    nnodel = elementype == "linear" ? 2 : 3
    z1 = copy(z0)
    River = Profile{RiverProfile}(x, z0, e2n, nn, Δx)
    Terrace = Profile{TerraceProfile}(x, z1, e2n, nn, Δx)
    # =========================================================================

    # FEM =====================================================================
    Scratch_FEM = sparsitymatrix(e2n, nel, nnodel)
    ScratchTerrace_FEM = ScratchTerraceFEM{Float64,Int64}(nn, Scratch_FEM.K)
    # =========================================================================

    # PHYSICAL PARAMETERS - Subaerial diffusion ===============================    
    # --- Rivers
    T = vanAlbada
    RiverPhysics = KnickPointsParameters{Float64, T}(
        Kr = 1e-6, # whipple and tucker 1999_ dimensional coeff of erosion
        h = 1.92, # whipple and tucker 1999
        kappa_a = 4.6071, # whipple and tucker 1999_ area-length coeff
        m = 0.5,
        n = 1.0,
        φ = FluxLimiter{T},
    )

    RiverArrays = KnickPointsArrays{typeof(x)}(River, nn, RiverPhysics)
    # --- Terrace
    D = 4.4e-4 # Fernandes and Dietrich, 1997 (Hillslope evolution...)
    # =========================================================================

    # PHYSICAL PARAMETERS - below sea level (from Malatesta for wave erosion) =
    TerracePhysics = TerraceParameters{Float64}(
        βz = βz, 
        h_wb = 15.0,
        P0 = 5e-5, # shallowest
        P_off = 5e-2, # offshore
    )
    # =========================================================================

    # PHYSICAL PARAMETERS - diffusion equation ================================
    HillSlopePhysics = HillSlopeParameters{Float64}(   
        κa = 5e-5, # check miguels paper
        diffusion = 0.1, # subaerial hill diffusion (K, m2/yr)
        alpha_dif = 1.0, # precipitation rate (m/yr)
    )
    # =========================================================================

    # PHYSICAL PARAMETERS - submarine diffusion (from Miguelito) ==============
    SubmarinePhysics = SubmarineParameters{Float64}(
        K_s = 1e2, # submarine diffusion coefficient
        λ = 5e-4, # submarine diffusion decay coefficient
    )
    # =========================================================================
    
    # SEA LEVEL ===============================================================
    starting_time  = sea_level_starting_time(sealevelopts)
    
    # SEA LEVEL ===============================================================
    sea_lvl_curve, sea_age  = get_sea_level(sealevelopts)
    sea_lvl_curve, sea_age, fsea = sea_level_corrections(sea_lvl_curve, sea_age, starting_time)
    # =========================================================================

    # # UPLIFT RATES FROM ARMEL =================================================
    # t_uplift = Vector{Float64}(undef,20)
    # uplift_background = Vector{Float64}(undef, 20)
    # if uplift_type == "armel"
    #     t_uplift, uplift_background = read_uplift()
    #     uplift_background = @. uplift_background * 2
        
    #     # -- define time interval
    #     t0          = 3e6        # starting time in [kyr]
    #     tf          = t0 + starting_time    # final time in [kyr]
    #     # -- take time interval from time and uplift rate series
    #     it0 = t_uplift .>= t0
    #     uplift_background = view(uplift_background, it0)
    #     t_uplift    = view(t_uplift, it0)
    #     t_uplift    = @. t_uplift - t_uplift[1]
    #     # =========================================================================
    #     fU = interpolate((t_uplift,), uplift_background, Gridded(Linear()))
    # end

    # SOLVER ==================================================================
    Δt = Δt
    break_time = maximum(sea_age.-Δt)
    t = 0.0    # initialise time

    # SOME PREALLOCATIONS    
    nit = Int64(floor(break_time / Δt))
    to = TimerOutput()
    t1 = time()
    buffer = similar(z0)
    id_shore = (
        terrace = (i1 = Terrace.nnod÷2, i2 = 1), 
        river = (i1 = River.nnod÷2, i2 = 1)
    )

    terrace_age = zeros(nn)

    @timeit to "main loop" for it = 1:nit

        # GET SEA LVL AND UPLIFT RATE =========================================
        @timeit to "interpolate sea level" h_sea = fsea(t)

        # -- find uplift rate at current time
        # @timeit to "interpolate uplift" U = (uplift_type == "constant" ? U0 : fU(t))::Float64
        @timeit to "interpolate uplift" U = U0

        # GET NODE BELOW SEA LEVEL===========================================
        @timeit to "id_shore" id_shore = find_shore_id(h_sea, River, Terrace, TerracePhysics, id_shore)
        
        # SOLVE EQUATIONS =====================================================
        @timeit to "solver" Terrace, River, to = solve(solver,
                TerracePhysics, RiverPhysics, RiverArrays, HillSlopePhysics, SubmarinePhysics,
                Δt, River, Terrace, h_sea, id_shore, buffer,
                ScratchTerrace_FEM, Scratch_FEM, to)

        # BACKGROUND UPLIFT ===================================================
        @timeit to "add uplift" uplift!(Terrace.z, River.z, U*Δt)
       
        # ERODED VOLUME =======================================================
        # push!(erosion, erodedvolume(sol0,sol,x))
        # erosionOut[it] = erodedvolume(sol0,River.z,x)

        t +=  Δt

        # TERRACE AGE =========================================================
        @turbo for i in id_shore.terrace.i1:length(x)
            terrace_age[i] = t
        end

        if t > break_time
            break
        end

    end

    t2 = time()
    ftime = t2-t1
    println(to)
    println("\n Finished after $ftime seconds ")

    return Terrace, River, terrace_age

end

function runall(δt, δx)
    # to  = TimerOutput()
    for dt in float(50:50:500), dx in δx
        
        Terrace,  = main(Δx = dx, Δt = dt)
        fname = string("output_dt_old/Potato_dt_",dt,"_dx_", dx,".h5")
        
        h5open(fname, "w") do file
                g           = create_group(file, "Terrace")       # create a group
                p           = create_group(file, "Parameters")  # create a group
                g["x"]      = Terrace.x                     # create a scalar dataset inside the group
                g["z"]      = Terrace.z                     # create a scalar dataset inside the group
                p["dt"]     = dt
                p["dx"]     = dx
        end
    end
    # to
end