using DelimitedFiles: Threads
using NumericalIntegration: Threads
include("Mixer.jl")

Δx, U0, βz =  30.0, 0.0001, 2e-6

function main(Δx, U0, βz, sea_level_file)
    # =========================================================================
    ievol = 1 # 1 = plot evolution step-by-step, 0 = only result
    igif = 0 # 1 = make gif of evolution
    yrs_Ma = 1 # 1 = use Ma, 0 = yrs since start
    isave = 0 # 1 = save evolution in txt file; 0 = don"t save
    uplift_type = "constant" # "constant" or "armel"
    solver = ModelOpts{TerraceOnly}
        # Options:
        #  - TerraceOnly
        #  - TerraceKnickpoint
        #  - TerraceDiffusion
        #  - TerraceKnickpointDiffusion
    path_sea_level = "/home/albert/Dropbox/Riverini/Pluto/"
    sea_level_file  = Bintanja 
        # no need to use  "" anymore. Options:
        #  - Spratt2016-450
        #  - Spratt2016-800
        #  - Bintanja
    sealevelopts = SeaLevelOpts{sea_level_file}(path_sea_level)
       
    dt_evol         = 1000          # plot evolution every X [yrs]
    # U0              = 0.002
    # =========================================================================
    
    # INITIAL GEOMETRY ========================================================
    slope = -1/50
    elementype = "linear"
    x,z0,dx2,nn,nel,e2n = mesher(Δx,slope,elementype)
    nnodel = elementype == "linear" ? 2 : 3
    z1 = copy(z0)
    River = Profile{RiverProfile}(x, z0, e2n, nn, Δx)
    Terrace = Profile{TerraceProfile}(x, z1, e2n, nn, Δx)
    # =========================================================================

    # FEM =====================================================================
    Scratch_FEM = sparsitymatrix(e2n, nel, nnodel)
    ScratchTerrace_FEM = ScratchTerraceFEM{Float64,Int64}(nn, Scratch_FEM.K)
    # =========================================================================

    # UPLIFTS =================================================================
    km                       = 1e3;
    eq_ctr                   = 0;            # # of eqs
    eq_time                  = 0;            # time counter for megathrust
    recurrence_t             = 150;          # [yr] megathrust eq recurrence time
    uplift                   = 0.2;          # [m] uplift per megathrust cycle
    range_onshore_pulse      = [0 1e3]*km;   # range where onshore uplift is active: 1st) closest to shoreline
    range_offshore_pulse	 = [0 1e3]*km;   # range where onshore uplift is active: 1st) closest to shoreline
    # =========================================================================

    # PHYSICAL PARAMETERS - Subaerial diffusion ===============================    
    # --- Rivers
    # K = 3.14e-7  # whipple and tucker 1999_ dimensional coeff of erosion
    T = vanAlbada
    RiverPhysics = KnickPointsParameters{Float64, T}(
        Kr = 5e-7, # whipple and tucker 1999_ dimensional coeff of erosion
        h = 1.92, # whipple and tucker 1999
        kappa_a = 4.6071, # whipple and tucker 1999_ area-length coeff
        m = 0.5,
        n = 1.1,
        φ = FluxLimiter{T},
    )
    RiverArrays = KnickPointsArrays{typeof(x)}(River, nn, RiverPhysics)
    # --- Terrace
    D           = 4.4e-4 # Fernandes and Dietrich, 1997 (Hillslope evolution...)
    # =========================================================================
    

    # PHYSICAL PARAMETERS - below sea level (from Malatesta for wave erosion) =
    βx = 1e-8
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
    Δt          = 1e2
    break_time  = maximum(sea_age.-Δt)
    t_plot      = 0.0
    t           = 0.0    # initialise time

    # SOME PREALLOCATIONS    
    nit         = Int64(floor(break_time / Δt))
    to          = TimerOutput()
    t1          = time()
    buffer1, buffer2 = similar(z0),similar(z0)
    id_shore = (terrace = Terrace.nnod÷2, river = River.nnod÷2)

    @timeit to "main loop" for it = 1:nit

        # GET SEA LVL AND UPLIFT RATE =========================================
        @timeit to "interpolate sea level" h_sea = fsea(t)

        # -- find uplift rate at current time
        # @timeit to "interpolate uplift" U = (uplift_type == "constant" ? U0 : fU(t))::Float64
        @timeit to "interpolate uplift" U = U0

        # GET NODE BELOW SEA LEVEL===========================================
        @timeit to "id_shore" id_shore = find_shore_id(h_sea, River, Terrace, id_shore)
        
        # SOLVE EQUATIONS =====================================================
        @timeit to "solver" Terrace, River, to = solve(solver,
                TerracePhysics, RiverPhysics, RiverArrays, HillSlopePhysics, SubmarinePhysics,
                Δt, River, Terrace, h_sea, id_shore, buffer1, buffer2,
                ScratchTerrace_FEM, Scratch_FEM, to)
          
        # BACKGROUND UPLIFT ===================================================
        @timeit to "add uplift" uplift!(Terrace.z, River.z, U*Δt)

        # # TERRACES SOLVER := MALATESTA ========================================
        # if imalatesta == 1
        #     # -- RIVER PROFILE
        #     @timeit to "newton river" terracenewton!(TerracePhysics, Δt * yr,River.z, h_sea, id_shore + 1, nn,
        #                                              buffer1, buffer2)
            
        #     # # -- TERRACE PROFILE            
        #     @timeit to "newton terrace" terracenewton!(TerracePhysics,Δt * yr,Terrace.z,h_sea,id_shore_terrace+1, nn,
        #                                                buffer1,buffer2)
        #                                                to
        #     # println(t)
        #     # @timeit to "cliff retreat" cliffretreat!(Terrace,h_sea,P0,h_wb,βx,Δt*yr)
        # end
        # # =====================================================================

        # # # TERRACE ==========================================================
        # # # Diffusion over Δt
        # # @timeit to "FEM Terrace" Terrace.z = femsolver(Δt,Δx,e2n,nel,nn, Terrace, h_sea,
        # #             K,M,iMC,jMC,
        # #             κa, alpha_dif, diffusion,
        # #             K_s,λ) 
        # # # ==================================================================

        # if solver == "FW" 
        #     # EXPLICIT ADVECTION -> FLUX LIMITER TVD =============================
        #     @timeit to "update r" updater!(r, Kr, kappa_a, L, h, m, Δx^2, x, River.z, nn)
        #     @timeit to "TVD"  TVD!(River.z,buffer1,n,φ,r,rf,Δt, Δx,id_shore)

        # elseif solver == "diffusion+FW"

        #     # DIFFUSION + ADVECTION ==============================================
        #     # Advection over Δt/2
        #     updater!(r, Kr, kappa_a, L, h, m, Δx, x, River.z, nn)
        #     @timeit to "stream power law" TVD!(River.z,buffer1,n,φ,r,rf,Δt/2, Δx,id_shore)

        #     # Diffusion over Δt (FEM)
        #     @timeit to "FEM River" River.z = femsolver(Δt,Δx,e2n,nel,nn, River, h_sea,
        #                                                K,M,iMC,jMC,        
        #                                                κa, alpha_dif, diffusion,
        #                                                K_s,λ)
              
        #     # Advection over Δt/2
        #     updater!(r, Kr, kappa_a, L, h, m, Δx, x, River.z, nn)
        #     @timeit to "stream power law" TVD!(River.z,buffer1,n,φ,r,rf,Δt/2, Δx,id_shore)
        #     # ==================================================================

        #     # TERRACE ==========================================================
        #     # Diffusion over Δt
        #     @timeit to "FEM Terrace" Terrace.z = femsolver(Δt,Δx,e2n,nel,nn, Terrace, h_sea,
        #                                                    K,M,iMC,jMC,
        #                                                    κa, alpha_dif, diffusion,
        #                                                    K_s,λ) 
        #     # femsolver(Δt,Δx,e2n,nel,nn, Terrace, h_sea,
        #     #     K,M,iMC,jMC,
        #     #     D) 
        #     # ==================================================================
        # else
        #     nothing
        # end
        # # ===================================================================

        # # MEGATHRUST ==========================================================
        # eq_time += Δt;   # time counter for megathrust
        # if eq_time >= recurrence_t
        #     # Terrace
        #     x_intersect = x[find_shore_id(h_sea, Terrace.z, nn)]
        #     id_eq = (x_intersect - range_onshore_pulse[2]) .<= x.<=(x_intersect - range_onshore_pulse[1])
        #     Terrace.z[id_eq] .+= uplift
        #     # River
        #     x_intersect = x[find_shore_id(h_sea, River.z, nn)]
        #     id_eq = (x_intersect - range_onshore_pulse[2]) .<= x.<=(x_intersect - range_onshore_pulse[1])
        #     River.z[id_eq] .+= uplift            
            
        #     eq_time = 0 # reset eq timer
        #     eq_ctr += 1 # eq counter
        # end
        # # =====================================================================

        # # BACKGROUND UPLIFT ===================================================
        # @timeit to "add uplift" uplift!(Terrace.z, River.z, U*Δt)
        # # =====================================================================

        # ERODED VOLUME =======================================================
        # push!(erosion, erodedvolume(sol0,sol,x))
        # erosionOut[it] = erodedvolume(sol0,River.z,x)
        # =====================================================================

        t +=  Δt

        if t > break_time
            break
        end

    end

    t2 = time()
    ftime = t2-t1
    println(to)
    println("\n Finished after $ftime seconds ")

    return Terrace, River

end


function runall(U,βz)
    to  = TimerOutput()
    ifile = 1 
    for i=1:length(U), j=1:length(βz)
        
        Terrace, River = main(20.0,U[i],βz[j]);
        fname = string("output/Potato_",ifile,".h5")
        
        h5open(fname, "w") do file
                g           = g_create(file, "Terrace")       # create a group
                r           = g_create(file, "Rivers")       # create a group
                p           = g_create(file, "Parameters")  # create a group
                g["x"]      = Terrace.x                     # create a scalar dataset inside the group
                g["z"]      = Terrace.z                     # create a scalar dataset inside the group
                r["z"]      = River.z                       # create a scalar dataset inside the group
                p["U"]      = U[i]
                p["beta_z"] = βz[j]
        end
        ifile+=1
    end
    to
end

U   = reverse(range(0.0,0.0002,length = 2))
βz  = reverse(range(1.3e-5, 1e-6,length = 2))

U   = 0
βz  = 1.3e-6
@time runall(U,βz)

Δx, U0, βz =  5.0, 0.0001, 2e-6
Terrace, River = main(5.0, 0.0001, 2e-6, Bintanja)
Terrace, River = main(30.0, 0.0001, 2e-6, "Spratt2016-450kyrs")
 
function foo()
    for _ in 1:50
        main(30.0, 0.0001, 2e-6, "Spratt2016-450kyrs");
    end
end
# options:
        #  - Spratt2016-450kyrs
        #  - Spratt2016-800kyrs
        #  - Bintanja
