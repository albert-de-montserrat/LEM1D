using BenchmarkTools,LinearAlgebra, SparseArrays, Interpolations, DelimitedFiles, NumericalIntegration,LoopVectorization,StaticArrays
using TimerOutputs, HDF5, CUDA
using Plots
# using CSV,DataFrames
include("src/read_txt_files.jl")
include("src/river_functions.jl")
include("src/solvers.jl")
include("src/others.jl")
include("src/FEM.jl")
# using Plots,
# gr()

function main(Δx,U0,βz)
    # =========================================================================
    iparallel       = 0            # 1 = parallel solvers
    ievol           = 1            # 1 = plot evolution step-by-step, 0 = only result
    igif            = 0            # 1 = make gif of evolution
    yrs_Ma          = 1            # 1 = use Ma, 0 = yrs since start
    isave           = 0            # 1 = save evolution in txt file; 0 = don"t save
    imalatesta      = 1            # 1 = malatesta on
    uplift_type     = "constant"   # "constant" or "armel"
    solver          = "FW"         # solver scheme
    # solver          = "diffusion+FW"; # solver scheme
    output_path     = "./"
    output_name     = "test.txt"
    sea_level_file  = "Spratt2016-800kyrs"
        # options:
        #  - Spratt2016-450kyrs
        #  - Spratt2016-800kyrs
        #  - Bintanja
    dt_evol         = 1000          # plot evolution every X [yrs]
    # U0              = 0.002
    # =========================================================================
    
    # INITIAL GEOMETRY ========================================================
    # profile = "slope";
    # options:
    #   "slope"         -> linear profile
    #   "equilibrium"   -> from excel: exponential decay + sea floor topo
    # Δx      =  1 # [m] distance between points
    slope                   = -1.0 / 50.0
    elementype              = "linear"
    x,z0,dx2,nn,nel,e2n     = mesher(Δx,slope,elementype)
    elementype == "linear" ? nnodel=2 : nnodel=3
    z1                      = copy(z0)
    River                   = Profile{RiverProfile}(x,z0)
    Terrace                 = Profile{TerraceProfile}(x,z1)
    # River                   = CuProfile{RiverProfile}(CuArray(x),CuArray(z0))
    # Terrace                 = CuProfile{TerraceProfile}(CuArray(x),CuArray(z1))
    # =========================================================================

    # FEM =====================================================================
    # nel,e2n                 = connectivitymatrix(nn)
    K,M,iMC,jMC             = sparsitymatrix(e2n,nel,nnodel)
    # =========================================================================

    # UPLIFTS =================================================================
    km                      = 1e3;
    eq_ctr                  = 0;            # # of eqs
    eq_time                 = 0;            # time counter for megathrust
    recurrence_t            = 150;          # [yr] megathrust eq recurrence time
    uplift                  = 0.2;          # [m] uplift per megathrust cycle
    range_onshore_pulse     = [0 1e3]*km;   # range where onshore uplift is active: 1st) closest to shoreline
    range_offshore_pulse	= [0 1e3]*km;   # range where onshore uplift is active: 1st) closest to shoreline
    # =========================================================================

    # PHYSICAL PARAMETERS - Subaerial diffusion ===============================    
    # --- Rivers
    # K = 3.14e-7  # whipple and tucker 1999_ dimensional coeff of erosion
    Kr          = 5e-7    # whipple and tucker 1999_ dimensional coeff of erosion
    h           = 1.92    # whipple and tucker 1999
    kappa_a     = 4.6071  # whipple and tucker 1999_ area-length coeff
    m           = 0.5
    n           = 1.1
    L           = riverlength(x,z0,nn)
    r           = @. -Kr * (kappa_a^m) * (L^(h * m)) / Δx # EQUATION    
    # Buffers
    rf          = zero(z0)
    φ           = selectfluxlimiter("vanAlbada")
    # --- Terrace
    D           = 4.4e-4 # Fernandes and Dietrich, 1997 (Hillslope evolution...)
    # =========================================================================
    
    # PHYSICAL PARAMETERS - below sea level (from Malatesta for wave erosion) =
    # βz  = 2.5e-6 # 1.3e-5 or 7.5e-6
    βx          = 1e-8
    # βz          = 1e-6 # 1.3e-5 or 7.5e-6
    h_wb        = 15.0
    P0          = 5e-5 # shallowest
    P_off       = 5e-2 # offshore
    # =========================================================================

    # PHYSICAL PARAMETERS - diffusion equation ================================
    κa          = 5e-5        # check miguels paper
    diffusion   = 0.25        # subaerial hill diffusion (K, m2/yr)
    alpha_dif   = 1.0         # precipitation rate (m/yr)
    # =========================================================================

    # PHYSICAL PARAMETERS - submarine diffusion (from Miguelito) ==============
    K_s         = 1e2         # submarine diffusion coefficient
    λ           = 5e-4        # submarine diffusion decay coefficient
    # =========================================================================

    # TIME CONSTANTS ==========================================================
    yr          = 60.0 * 60.0 * 24.0 * 365.0
    day         = 60.0 * 60.0 * 24.0    
    # =========================================================================

    # SEA LEVEL ===============================================================
    starting_time           = 850e3                    # take the las X years of the sea lvl variations curve
    sea_lvl_curve, sea_age  = get_sea_level(sea_level_file)
    sea_age                 = sea_age[end:-1:1]        # reverse age array -> 0 = oldest age
    id_time                 = sea_age .> maximum(sea_age) .- starting_time
    sea_age                 = sea_age[id_time]
    sea_age                 = sea_age .- minimum(sea_age)
    sea_lvl_curve           = sea_lvl_curve[id_time]
    fsea                    = interpolate((reverse(sea_age),), sea_lvl_curve, Gridded(Linear()))
    # =========================================================================

    # UPLIFT RATES FROM ARMEL =================================================
    if uplift_type == "armel"
        t_uplift, uplift_background = read_uplift()
        uplift_background = @. uplift_background * 2
        # -- define time interval
        t0          = 38.9e6        # starting time in [kyr]
        tf          = t0 + 450e3    # final time in [kyr]
        # -- take time interval from time and uplift rate series
        uplift_background =
            uplift_background[@. (t_uplift <= tf) == (t_uplift >= t0)]
        t_uplift    = t_uplift[@. (t_uplift <= tf) == (t_uplift >= t0)]
        t_uplift    = @. t_uplift - t_uplift[1]
    end
    # =========================================================================

    # SOLVER ==================================================================
    Δt          = 1e2
    # sol         = copy(z0)   # initial solution = initial profile
    # terrace     = copy(sol)
    t_plot      = 0.0
    t           = 0.0    # initialise time
    id_shore    = 1
    # erosion     = Float64[]

    # SOME PREALLOCATIONS    
    # h2sea_terrace = zeros(nn)
    nit         = Int64(floor(maximum(sea_age .- Δt) / Δt))
    erosionOut  = Vector{Float64}(undef,nit)
    river_z     = Array{Float64}(undef,nn,nit)
    terrace_z   = similar(river_z)
    to          = TimerOutput()
    t1          = time()
    buffer1,buffer2 = similar(z0),similar(z0)

    # @sync begin
    #     Threads.@spawn River = river(River,βz,P0,h_wb,
    #                     n,φ,r,rf, Kr, kappa_a, L, h, m,
    #                     sea_age,U0,fsea,
    #                     K,M,iMC,jMC,D,
    #                     e2n,nel,
    #                     Δx,Δt, yr
    #                     )

    #     Terrace = terrace(Terrace,βz,P0,h_wb,
    #                     sea_age,U0,fsea, 
    #                     K,M,iMC,jMC,D,
    #                     e2n,nel,
    #                     Δx,Δt, yr
    #                     )

    # end


# while t < maximum(sea_age .- Δt)
    for it = 1:nit

        # GET SEA LVL AND UPLIFT RATE =========================================
        h_sea = fsea(t)

        # -- find uplift rate at current time
        # if uplift_type == "armel"
        #     # -- interpolate
        #     U = interpolate((t_uplift,), uplift_background, Gridded(Linear()))(t)
        # elseif uplift_type == "constant"
            U = U0
        # end
        # =====================================================================

        # GET COORDS BELOW SEA LEVEL===========================================
        id_shore            = find_shore_id(h_sea, River.z, nn)
        id_shore_terrace    = find_shore_id(h_sea, Terrace.z, nn)
        # =====================================================================

        # TERRACES SOLVER := MALATESTA ========================================
        if imalatesta == 1
            # -- RIVER PROFILE
            @timeit to "newton river" terracenewton!(βz,P0,h_wb,Δt * yr,River.z,h_sea,id_shore + 1,nn,
                                                     buffer1,buffer2)
            # # -- TERRACE PROFILE            
            @timeit to "newton terrace" terracenewton!( βz,P0,h_wb,Δt * yr,Terrace.z,h_sea,id_shore_terrace+1,nn,
                                                        buffer1,buffer2)
         
            # println(t)
            # @timeit to "cliff retreat" cliffretreat!(Terrace,h_sea,P0,h_wb,βx,Δt*yr)
        end
        # =====================================================================

        # TERRACE ==========================================================
        # Diffusion over Δt
        # @timeit to "FEM Terrace" Terrace.z .= femsolver(Δt,Δx,e2n,nel,nn, Terrace, h_sea,
        #     K,M,iMC,jMC,
        #     D) 
        # ==================================================================

        if solver == "FW" 
            # EXPLICIT ADVECTION -> FLUX LIMITER TVD =============================
            @timeit to "update r" updater!(r, Kr, kappa_a, L, h, m, Δx^2, x, River.z, nn)
            @timeit to "TVD"  TVD!(River.z,buffer1,n,φ,r,rf,Δt, Δx,id_shore)

        elseif solver == "diffusion+FW"

            # DIFFUSION + ADVECTION ==============================================
            # Advection over Δt/2
            updater!(r, Kr, kappa_a, L, h, m, Δx, x, River.z, nn)
            @timeit to "stream power law" TVD!(River.z,buffer1,n,φ,r,rf,Δt, Δx,id_shore)

            # Diffusion over Δt (FEM)
            @timeit to "FEM River" River.z = femsolver(Δt,Δx,e2n,nel,nn, River, h_sea,
                K,M,iMC,jMC,
                κa, alpha_dif, diffusion,
                K_s,λ) 
              
            # Advection over Δt/2
            updater!(r, Kr, kappa_a, L, h, m, Δx, x, River.z, nn)
            @timeit to "stream power law" TVD!(River.z,buffer1,n,φ,r,rf,Δt, Δx,id_shore)
            # ==================================================================

            # TERRACE ==========================================================
            # Diffusion over Δt
            @timeit to "FEM Terrace" Terrace.z = femsolver(Δt,Δx,e2n,nel,nn, Terrace, h_sea,
                K,M,iMC,jMC,
                D) 
            # ==================================================================
          
        end
        # ===================================================================

        # MEGATHRUST ==========================================================
        # eq_time = eq_time + Δt;   # time counter for megathrust
        # if eq_time >= recurrence_t
        #     x_intersect                     = interpolate((sol,), x, Gridded(Linear()))(h_sea);
        #     id_onshore                      = x<=(x_intersect - range_onshore_pulse[1] & ...
        #         x>=(x_intersect - range_onshore_pulse[2]);
        #     id_offshore                     = x>=(x_intersect + range_offshore_pulse[1]) & ...
        #         x<=(x_intersect + range_offshore_pulse[2]);
        #     eq_ctr                          = eq_ctr + 1;
        #     terrace[id_onshore|id_offshore]	= terrace[id_onshore | id_offshore] + uplift;
        #     sol[id_onshore | id_offshore]	= sol[id_onshore | id_offshore] + uplift;
        #     eq_time                         = 0;
        # end
        # =====================================================================

        # # BACKGROUND UPLIFT ===================================================
        # Terrace.z  .+= U * Δt
        # River.z    .+= U * Δt
        # # =====================================================================
        
        # ERODED VOLUME =======================================================
        # push!(erosion, erodedvolume(sol0,sol,x))
        # erosionOut[it] = erodedvolume(sol0,River.z,x)
        # =====================================================================

        # # FILL OUTPUT PROFILES ================================================
        # @views begin
        #     river_z[:,it]      = River.z
        #     terrace_z[:,it]    = Terrace.z
        # end
        # # =====================================================================
        
        t +=  Δt

        if t > maximum(sea_age.-Δt) #|| t > t_uplift(end)
            break
        end

    end

    t2 = time()
    ftime = t2-t1
    println(to)
    println("\n Finished after $ftime seconds ")

    # return x, river_z, terrace_z, erosionOut,to
    return Terrace, River

end


function runall(U,βz)
    to  = TimerOutput()
    ifile = 1 
    for i=1:length(U), j=1:length(βz)
        
        Terrace, River = main(20.0,U[i],βz[j]);
        fname = string("output/Potato_",ifile,".h5")
        
        h5open(fname, "w") do file
                g           = g_create(file, "River")       # create a group
                p           = g_create(file, "Parameters")  # create a group
                g["x"]      = Terrace.x                     # create a scalar dataset inside the group
                g["z"]      = Terrace.z                     # create a scalar dataset inside the group
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

Δx,U0,βz =  30.0,0.0,2e-5
main(30.0,0.0,2e-5);