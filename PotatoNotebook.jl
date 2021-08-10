### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ cc42df60-7ab2-11eb-16c2-4d21ab90084d
begin
	using Markdown
	using InteractiveUtils
	using BenchmarkTools,LinearAlgebra, SparseArrays, Interpolations, DelimitedFiles, NumericalIntegration,LoopVectorization,StaticArrays
	using TimerOutputs, HDF5, CUDA
	using Plots
end

# ╔═╡ 4a6e9172-7ab3-11eb-35e5-f52c7b6a0574
begin
	# using CSV,DataFrames
	include("/home/albert/Documents/LEM/src/read_txt_files.jl")
	include("/home/albert/Documents/LEM/src/river_functions.jl")
	include("/home/albert/Documents/LEM/src/solvers.jl")
	include("/home/albert/Documents/LEM/src/others.jl")
	include("/home/albert/Documents/LEM/src/FEM.jl")
end;

# ╔═╡ c095df4e-7ab4-11eb-0d31-c33fe91ff3b2
md"# Main function"

# ╔═╡ 4e811130-7ab4-11eb-1c78-57910763f899

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
    sea_level_file  = "Spratt2016-450kyrs"
    # sea_level_file  = "Spratt2016-800kyrs"
        # options:
        #  - Spratt2016-450kyrs
        #  - Spratt2016-800kyrs
        #  - Bintanja
    dt_evol         = 1000          # plot evolution every X [yrs]
    # =========================================================================
    
    # INITIAL GEOMETRY ========================================================
    # profile = "slope";
    # options:
    #   "slope"         -> linear profile
    #   "equilibrium"   -> from excel: exponential decay + sea floor topo
    # Δx      =  1 # [m] distance between points
    slope                   = -1/50 
    elementype              = "linear"
    x,z0,dx2,nn,nel,e2n     = mesher(Δx,slope,elementype)
    nnodel = elementype == "linear" ? 2 : 3
    z1                      = copy(z0)
    River                   = Profile{RiverProfile}(x,z0)
    Terrace                 = Profile{TerraceProfile}(x,z1)
    # =========================================================================

    # FEM =====================================================================
    K,M,iMC,jMC = sparsitymatrix(e2n,nel,nnodel)
    # =========================================================================

    # UPLIFTS =================================================================
    km = 1e3;
    eq_ctr = 0;            # # of eqs
    eq_time = 0;            # time counter for megathrust
    recurrence_t = 150; # [yr] megathrust eq recurrence time
    uplift = 0.2; # [m] uplift per megathrust cycle
    range_onshore_pulse = [0 1e3]*km;   # range where onshore uplift is active: 1st) closest to shoreline
    range_offshore_pulse = [0 1e3]*km;   # range where onshore uplift is active: 1st) closest to shoreline
    # =========================================================================

    # PHYSICAL PARAMETERS - Subaerial diffusion ===============================    
    # --- Rivers
    # K = 3.14e-7  # whipple and tucker 1999_ dimensional coeff of erosion
    Kr          = 5e-4    # whipple and tucker 1999_ dimensional coeff of erosion
    h           = 1.92    # whipple and tucker 1999
    kappa_a     = 4.6071  # whipple and tucker 1999_ area-length coeff
    m           = 0.2
    n           = 1.0
    L           = riverlength(x,z0,nn)
    r           = @. -Kr * (kappa_a^m) * (L^(h * m)) / Δx # EQUATION    
    # Buffers
    rf          = zero(z0)
    φ           = selectfluxlimiter("vanAlbada")
    # --- Terrace
    D           = 4.4e-4 # Fernandes and Dietrich, 1997 (Hillslope evolution...)
    # =========================================================================
    
    # PHYSICAL PARAMETERS - below sea level (from Malatesta for wave erosion) =
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
    starting_time = 3e6  # take the las X years of the sea lvl variations curve
    sea_lvl_curve, sea_age  = get_sea_level(sea_level_file)
    sea_age = sea_age[end:-1:1]        # reverse age array -> 0 = oldest age
    id_time = sea_age .> maximum(sea_age) .- starting_time
    sea_age = sea_age[id_time]
    sea_age = sea_age .- minimum(sea_age)
    sea_lvl_curve = sea_lvl_curve[id_time]
    fsea = interpolate((reverse(sea_age),), sea_lvl_curve, Gridded(Linear()))
    # =========================================================================

    # UPLIFT RATES FROM ARMEL =================================================
    if uplift_type == "armel"
        t_uplift, uplift_background = read_uplift()
        uplift_background = @. uplift_background * 2
        # -- define time interval
        t0          = 6e6        # starting time in [kyr]
        tf          = t0 + 800e3    # final time in [kyr]
        # -- take time interval from time and uplift rate series
        it0 = t_uplift .>= t0
        uplift_background = view(uplift_background, it0)
        t_uplift    = view(t_uplift, it0)
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
    nit         = Int64(floor(maximum(sea_age .- Δt) / Δt))
    # erosionOut  = Vector{Float64}(undef,nit)
    hseaOut     = Vector{Float64}(undef,nit)
    river_z     = Array{Float64}(undef,nn,nit)
    terrace_z   = similar(river_z)
    to          = TimerOutput()
    t1          = time()
    buffer1, buffer2 = similar(z0),similar(z0)

    for it = 1:nit

        # GET SEA LVL AND UPLIFT RATE =========================================
        h_sea = fsea(t)
	    hseaOut[it] = h_sea
        # -- find uplift rate at current time
        if uplift_type == "armel"
            # -- interpolate
            U = interpolate((t_uplift,), uplift_background, Gridded(Linear()))(t)
        elseif uplift_type == "constant"
            U = U0
        end
        # =====================================================================

        # GET COORDS BELOW SEA LEVEL===========================================
        id_shore            = find_shore_id(h_sea, River.z, nn)
        id_shore_terrace    = find_shore_id(h_sea, Terrace.z, nn)
        # id_shore_terrace,i2 = find_shore_id_terrace(h_sea, Terrace.z, nn)
        # =====================================================================

        # TERRACES SOLVER := MALATESTA ========================================
        if imalatesta == 1
            # -- RIVER PROFILE
            @timeit to "newton river" terracenewton!(βz,P0,h_wb,Δt * yr,River.z, h_sea, id_shore + 1, nn, buffer1,buffer2)
            # # -- TERRACE PROFILE            
            @timeit to "newton terrace" terracenewton!( βz,P0,h_wb,Δt * yr,Terrace.z, h_sea, id_shore_terrace+1, nn, buffer1, buffer2)
         
            # println(t)
            # @timeit to "cliff retreat" cliffretreat!(Terrace,h_sea,P0,h_wb,βx,Δt*yr)
        end
        # =====================================================================

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
              
            # # Advection over Δt/2
            # updater!(r, Kr, kappa_a, L, h, m, Δx, x, River.z, nn)
            # @timeit to "stream power law" TVD!(River.z,buffer1,n,φ,r,rf,Δt, Δx,id_shore)
            # ==================================================================

            # TERRACE ==========================================================
            # Diffusion over Δt
            @timeit to "FEM Terrace"  Terrace.z = femsolver(Δt,Δx,e2n,nel,nn, Terrace, 					   h_sea, K, M, iMC, jMC, κa, alpha_dif, diffusion, K_s, λ)
            # ==================================================================
        else
            nothing
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
        Terrace.z .+= U * Δt
        River.z .+=  U * Δt
        # # =====================================================================
        
        # ERODED VOLUME =======================================================
        # push!(erosion, erodedvolume(sol0,sol,x))
        # erosionOut[it] = erodedvolume(sol0,River.z,x)
        # =====================================================================

        # FILL OUTPUT PROFILES ================================================
        @views begin
            @inbounds river_z[:,it]      = River.z
            @inbounds terrace_z[:,it]    = Terrace.z
        end
        # =====================================================================
        
        t +=  Δt
		
		# if t > 100e3
		# 	break
		# end

        if t > maximum(sea_age.-Δt) #|| t > t_uplift(end)
            break
        end

    end

    # return x, river_z, terrace_z, erosionOut,to
    return Terrace, River, river_z, terrace_z, hseaOut, to

end

# ╔═╡ a2516ab2-7ab4-11eb-109c-81032adf42ed
md"# Input parameters"

# ╔═╡ 64ae3dfc-7ab4-11eb-0029-a3bb8e4aa170
begin
	Δx = 5.0 # lateral resolution of the terrace/river profile
	U0 = 0.0004 # passive uplift
	βz = 5e-6 # vertical incision (power? rate?)
end;

# ╔═╡ b4c49c6c-7ab4-11eb-28b6-b55de512a167
md"# Run model"

# ╔═╡ 65cefaa2-7ab4-11eb-2b2a-8f94f3bd4da2
Terrace, River, river_z, terrace_z, hseaOut, to = main(Δx, U0, βz);

# ╔═╡ 0c9107c2-7ab5-11eb-2201-af66cca99e08
md"# Plots"

# ╔═╡ 14386c90-7ab8-11eb-22d1-8ba1f2393e60
to

# ╔═╡ 1a8e60d4-7ab5-11eb-3fa0-bb41326eee86
begin
	hline([hseaOut[end]],
		color=:blue,
		label="sea")
	
	plot!(Terrace.x/1e3, River.z,
		color=:red,
		label="river")
	
	plot!(Terrace.x/1e3, Terrace.z,
		color=:green,
		label="terrace")

	ylims!(-150, 200)
	xlims!(95, 120)
end

# ╔═╡ 5a25a8fc-985d-42ed-ac1e-cad6ca49ca41
River.z

# ╔═╡ 838077a4-7abc-11eb-36e7-05060956e8e5
anim = @animate for it in 1:15:length(hseaOut)
	# t2 = time[it]/1e3
	
	hline([hseaOut[it]],
		color=:blue,
		label="sea level")
	
	plot!(Terrace.x/1e3, river_z[:,it],
		color=:red,
		label="river")
	
	plot!(Terrace.x/1e3, terrace_z[:,it],
		color=:green,
		label="terrace")

	ylims!(-150, 200)
	xlims!(95, 110)
end;

# ╔═╡ 03a4476e-7abd-11eb-0a9f-dbe1bdd74782
gif(anim, "anim_fps15.gif", fps = 1000)

# ╔═╡ Cell order:
# ╠═cc42df60-7ab2-11eb-16c2-4d21ab90084d
# ╟─4a6e9172-7ab3-11eb-35e5-f52c7b6a0574
# ╟─c095df4e-7ab4-11eb-0d31-c33fe91ff3b2
# ╠═4e811130-7ab4-11eb-1c78-57910763f899
# ╟─a2516ab2-7ab4-11eb-109c-81032adf42ed
# ╠═64ae3dfc-7ab4-11eb-0029-a3bb8e4aa170
# ╟─b4c49c6c-7ab4-11eb-28b6-b55de512a167
# ╠═65cefaa2-7ab4-11eb-2b2a-8f94f3bd4da2
# ╟─0c9107c2-7ab5-11eb-2201-af66cca99e08
# ╠═14386c90-7ab8-11eb-22d1-8ba1f2393e60
# ╠═1a8e60d4-7ab5-11eb-3fa0-bb41326eee86
# ╠═5a25a8fc-985d-42ed-ac1e-cad6ca49ca41
# ╠═838077a4-7abc-11eb-36e7-05060956e8e5
# ╠═03a4476e-7abd-11eb-0a9f-dbe1bdd74782
