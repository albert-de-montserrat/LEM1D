# include("C:/Users/crosetto/Documents/TerracesJulia/LEM1D-main/src/FEM.jl")
# using Pkg; Pkg.activate(".")
using LEM1D
using LinearAlgebra, Interpolations, LoopVectorization, Printf

function main(slope, βz, h_wb, recurrence_t, uplift)
    # =========================================================================
    uplift_type = "constant"   # "constant" or "armel"
    sea_level_file = "Spratt2016-450"
    # options:
    #  - Spratt2016-450
    #  - Spratt2016-800
    #  - Bintanja3
    #  - Bintanja6
    #  - Bintanja3x2
    U0 = 0.002
    # =========================================================================

    # INITIAL GEOMETRY ========================================================
    # profile = "slope";
    # options:
    #   "slope"         -> linear profile
    #   "equilibrium"   -> from excel: exponential decay + sea floor topo
    Δx = 5.0 # [m] distance between points
    elementype = "linear"
    L = 200e3
    Terrace = TerraceProfile(L, Δx, slope; intercept=2e3)
    x, z0, nn = deepcopy(Terrace.x), deepcopy(Terrace.z), Terrace.nnod
    z1 = copy(z0)
    # =========================================================================

    # UPLIFTS =================================================================
    eq_ctr = 0 # of eqs
    eq_time = 0 # time counter for megathrust
    # =========================================================================

    # PHYSICAL PARAMETERS - below sea level (from Malatesta for wave erosion) =
    P0 = 5e-5 # shallowest
    # =========================================================================

    # TIME CONSTANTS ==========================================================
    yr = 60.0 * 60.0 * 24.0 * 365.0
    # =========================================================================

    # SEA LEVEL ===============================================================
    # take the las X years of the sea lvl variations curve
    if sea_level_file === "Spratt2016-450"
        starting_time = 430e3
    elseif sea_level_file === "Spratt2016-800"
        starting_time = 798e3
    elseif sea_level_file === "Siddall2003-379"
        starting_time = 378.96e3
    elseif sea_level_file === "Bintanja3"
        starting_time = 3e6
    elseif sea_level_file === "Bintanja6" || "Bintanja3x2"
        starting_time = 6e6
    end
    sea_lvl_curve, sea_age = get_sea_level(sea_level_file)
    sea_age = sea_age[end:-1:1] # reverse age array -> 0 = oldest age
    id_time = sea_age .>= maximum(sea_age) .- starting_time
    sea_age = sea_age[id_time]
    sea_lvl_curve = sea_lvl_curve[id_time]
    fsea = interpolate((reverse(sea_age),), reverse(sea_lvl_curve), Gridded(Linear()))
    # =========================================================================

    # UPLIFT RATES FROM ARMEL =================================================
    if uplift_type === "armel"
        fname = fuplift
        t_uplift, uplift_background = read_uplift(fname)
        uplift_background = @. uplift_background    # * 2 when Armel's U is too low   
        # -- define time interval
        t0 = 0        # starting time in [kyr]
        tf = t0 + starting_time    # final time in [kyr] tsh = screenshot time
        # -- take time interval from time and uplift rate series
        it0 = t_uplift .>= t0
        uplift_background = view(uplift_background, it0)
        t_uplift = view(t_uplift, it0)
        t_uplift = @. t_uplift - t_uplift[1]
        # ======================================================================
        fU = interpolate((t_uplift,), uplift_background, Gridded(Linear()))
    end
    # =========================================================================

    # SOLVER ==================================================================
    Δt = 50.0 # time step (years)
    t_plot = 0.0
    t = sea_age[end]     # initialise time
    id_shore = 1

    # SOME PREALLOCATIONS    
    nit = Int64(floor(maximum(sea_age .- Δt) / Δt))
    tOut = Vector{Float64}(undef, nit)
    buffer1, buffer2 = deepcopy(Terrace.z), deepcopy(Terrace.z)
    t1 = time()
    terrace_age = zeros(nn)
    reoccupation_time = zeros(nn, 3)
    reoccupation_id = fill(:aerial, nn)
    reoccupation_id[Terrace.z .≤ fsea(t)] .= :submarine

    break_time = maximum(sea_age .- Δt)
    println("Starting solver... \n")
    U = U0
    for it in 1:nit

        # GET SEA LVL AND UPLIFT RATE =========================================
        h_sea = fsea(t)
        # -- find uplift rate at current time
        if uplift_type === "armel"
            # -- interpolate
            U = fU(t)
        end
        # =====================================================================

        # GET COORDS BELOW SEA LEVEL===========================================
        id_shore_terrace = find_shore_id(h_sea, Terrace.z, nn)
        # =====================================================================

        # TERRACES SOLVER =====================================================
        terracenewton!(
            βz,
            P0,
            h_wb,
            Δt * yr,
            Terrace.z,
            h_sea,
            id_shore_terrace + 1,
            nn,
            buffer1,
            buffer2,
        )
        # =====================================================================

        # MEGATHRUST ==========================================================
        eq_time += Δt # time counter for megathrust
        if eq_time ≥ recurrence_t
            Terrace.z .+= uplift
            eq_time = 0 # reset eq timer
            eq_ctr += 1 # eq counter
        end
        # =====================================================================
        
        # BACKGROUND UPLIFT ===================================================
        Terrace.z .+= U * Δt
        # =====================================================================

        # REOCUPATION==========================================================
        update_reoccupation!(reoccupation_time, reoccupation_id, t, h_sea, Terrace.z)
        # =====================================================================

        t += Δt
        tOut[it] = t

        # TERRACE AGE =========================================================
        for i in id_shore_terrace:length(x)
            terrace_age[i] = t
        end
        # =====================================================================

        t > break_time && break
    end

    @show eq_ctr
    t2 = time()
    ftime = t2 - t1
    println("\n Finished after $ftime seconds ")
    println("Final time = ", t)

    return Terrace, tOut, terrace_age, reoccupation_time
end

deg2slope(deg) = tand(deg)
degs = -10
slope = deg2slope(degs)
h_wb = 15.0
βz = 2.0e-6
uplift = 0
recurrence_t = 0

main(slope, βz, h_wb, recurrence_t, uplift)