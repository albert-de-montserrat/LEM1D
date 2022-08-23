using LEM1D
using LinearAlgebra, Interpolations, LoopVectorization, Printf

function main(slope, βz, h_wb, recurrence_t, uplift)
    # =========================================================================
    uplift_type = "constant"
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
    Δx = 5.0 # [m] resolution
    L = 200e3
    Terrace = TerraceProfile(L, Δx, slope; intercept=2e3)
    x, z0, nn = deepcopy(Terrace.x), deepcopy(Terrace.z), Terrace.nnod
    z1 = copy(z0)
    # =========================================================================

    # UPLIFTS =================================================================
    eq_ctr = 0 # of eqs
    eq_time = 0 # time counter for megathrust
    # =========================================================================

    # PHYSICAL PARAMETERS =====================================================
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

    # SOLVER ==================================================================
    Δt = 50.0 # time step (years)
    t = sea_age[end]     # initialise time

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
    U = U0
    println("Starting solver... \n")
    for it in 1:nit

        # GET SEA LEVEL =======================================================
        h_sea = fsea(t)
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