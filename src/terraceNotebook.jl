### A Pluto.jl notebook ###
# v0.12.6

using Markdown
using InteractiveUtils

# ╔═╡ 75f53f00-0ae0-11eb-2523-e988a6f2e181
begin
	using BenchmarkTools,LinearAlgebra, Interpolations, DelimitedFiles, Plots, LoopVectorization, NumericalIntegration
	include("read_txt_files.jl")
	include("river_functions.jl")
	include("solvers.jl")
	include("others.jl")
end

# ╔═╡ 8a912dc0-0ae2-11eb-03f9-0f451ad8b10c
begin # some flags
	iparallel       = 0            # 1 = parallel solvers
	ievol           = 1            # 1 = plot evolution step-by-step, 0 = only 
	igif            = 0            # 1 = make gif of evolution
	yrs_Ma          = 1            # 1 = use Ma, 0 = yrs since start
    isave           = 0            # 1 = save evolution in txt file; 0 = don"t 
	imalatesta      = 1            # 1 = malatesta on
    uplift_type     = "constant"   # "constant" or "armel"
	solver          = "FW"         # solver scheme
	output_path     = "./"
    output_name     = "test.txt"
	sea_level_file  = "Spratt2016-800kyrs"
	        # options:
	        #  - Spratt2016-450kyrs
	        #  - Spratt2016-800kyrs
	        #  - Bintanja
	dt_evol         = 1000          # plot evolution every X [yrs]
    U0              = 0.001 		# passive uplift
	yr  			= 60.0 * 60.0 * 24.0 * 365.0
end

# ╔═╡ c55008b2-0ae2-11eb-05fe-35f20d380ea4
begin # INITIAL GEOMETRY ========================================================
	# profile = "slope";
	# options:
    #   "slope"         -> linear profile
	#   "equilibrium"   -> from excel: exponential decay + sea floor topo
    dx      =  20.0 # [m] distance between points
	slope   = -1.0 / 50.0
    x       = collect(0.0:dx:200e3)
    dx2     = dx * dx
	nn      = length(x)
	sol0    = @. slope * x + 2e3    # 2e3 is the intercept
	terrace0= copy(sol0)    # 2e3 is the intercept
end

# ╔═╡ 09195ea4-0ae3-11eb-08e4-c99929247d28
begin # UPLIFTS =================================================================
	km                      = 1e3;
    eq_ctr                  = 0;            # # of eqs
	eq_time                 = 0;            # time counter for megathrust
    recurrence_t            = 150;          # [yr] megathrust eq recurrence time
	uplift                  = 0.2;          # [m] uplift per megathrust cycle
    range_onshore_pulse     = [0 1e3]*km;   # range where onshore uplift is active: 1st) closest to shoreline
	range_offshore_pulse	= [0 1e3]*km;   # range where onshore uplift is active: 1st) closest to shoreline
end

# ╔═╡ 26c54580-0ae3-11eb-3263-f176fad834ee
begin # PHYSICAL PARAMETERS - stream power law ====================================
	L = abs(x[1] - x[end])
    # K = 3.14e-7  # whipple and tucker 1999_ dimensional coeff of erosion
    K = 5e-6  # whipple and tucker 1999_ dimensional coeff of erosion
    h = 1.92  # whipple and tucker 1999
    kappa_a = 4.6071  # whipple and tucker 1999_ area-length coeff
    m = 0.5
    n = 1
    r = @. -K * (kappa_a^m) * (x^(h * m)) / dx 
end

# ╔═╡ 457170e4-0ae3-11eb-149e-6d96ff1ea287
begin # PHYSICAL PARAMETERS - terrace formation (from Malatesta for wave erosion) =
	# beta_z  = 2.5e-6 # 1.3e-5 or 7.5e-6
    beta_z  = 1e-6 # 1.3e-5 or 7.5e-6
	h_wb    = 100.0
    P0      = 5e-5 # shallowest
	P_off   = 5e-2 # offshore
end

# ╔═╡ 5740103c-0ae3-11eb-31d5-d795e058982b
begin # PHYSICAL PARAMETERS - diffusion equation ================================
	kappa     = 5e-5        # check miguels paper
    diffusion = 1e-3        # subaerial hill diffusion (K, m2/yr)
	alpha_dif = 1.0         # precipitation rate (m/yr)
end

# ╔═╡ 97929f92-0ae3-11eb-119b-71390eb9e1d3
begin # SEA LEVEL ===============================================================
	starting_time           = 450e3                    # take the las X years of the sea lvl variations curve
    sea_lvl_curve, sea_age  = get_sea_level(sea_level_file)
	sea_age                 = sea_age[end:-1:1]        # reverse age array -> 0 = oldest age
    id_time                 = sea_age .> maximum(sea_age) .- starting_time	  
	sea_age                 = sea_age[id_time]
    sea_age                 = sea_age .- minimum(sea_age)
	sea_lvl_curve           = sea_lvl_curve[id_time]	
end

# ╔═╡ c602f700-0ae3-11eb-3005-9d5af1882906
begin # UPLIFT RATES FROM ARMEL =================================================
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
	else
		uplift_background,t_uplift=1,1
	end
end

# ╔═╡ e865764a-0ae7-11eb-200c-e7fc36fd90d3
begin
	# PRE-SOLVER ==================================================================
	dt          = 1e2 			# time step
	t_plot      = 0.0
	t           = 0.0    		# initialise time
	id_shore    = 1	
end

# ╔═╡ 35efe2a8-0ae4-11eb-393f-3be330f265da
function mainsolver(sea_age,sea_lvl_curve,t,dt,yr,
		x,sol,terrace,
		t_uplift,uplift_background,U0,uplift_type,
		beta_z,P0,h_wb,imalatesta,
		n,r,dx)

	nn 		 	= length(x)
	tOut 	 	= Float64[]
	hseaOut  	= Float64[]
	# erosionOut	= Float64[]
	riverOut  	= Array{Float64,1}[]
	terraceOut 	= Array{Float64,1}[]
	sol0 		= copy(sol)
	
	while t < maximum(sea_age .- dt)

        # GET SEA LVL AND UPLIFT RATE =========================================
        h_sea = interpolate((reverse(sea_age),), sea_lvl_curve, Gridded(Linear()))(t)

        # -- find uplift rate at current time
        if uplift_type == "armel"
            # -- interpolate
            U = interpolate((t_uplift,), uplift_background, Gridded(Linear()))(t)
        elseif uplift_type == "constant"
            U = U0
        end
        # =====================================================================

        # GET COORDS BELOW SEA LEVEL===========================================
        id_shore        = find_shore_id(h_sea, sol, nn)
        # =====================================================================

        # SOLVER BELOW SEA LEVEL = MALATESTA ==================================
        if imalatesta == 1
            # -- RIVER PROFILE
            # terracenewton!(beta_z,P0,h_wb,dt * yr,sol,h_sea,id_shore + 1,nn)
           
            # -- TERRACE PROFILE
            id_shore_terrace  = find_shore_id(h_sea, terrace, nn)
            terracenewton!(beta_z,P0,h_wb,dt * yr,terrace,h_sea,id_shore_terrace+1,nn)
		end
        # =====================================================================
        
        # SOLVER ABOVE SEA LEVEL ==============================================
        solve_advection!(sol, n, r, dx, dt, 1, id_shore)        
        # =====================================================================
		
		# # ERODED VOLUME =======================================================
		# if t > 0.0
		# 	erosion = erodedvolume(riverOut[end],sol,x,id_shore)
		# else
		# 	erosion = erodedvolume(sol0,sol,x,id_shore)
		# end
		# # =====================================================================

        # BACKGROUND UPLIFT ===================================================
        terrace  .+= U * dt
        sol      .+= U * dt
        # =====================================================================
        	 
		t += dt

		push!(tOut, t)
		push!(hseaOut, h_sea)
		# push!(erosionOut, erosion)
		push!(riverOut, copy(sol))
		push!(terraceOut, copy(terrace))

        # if t > 1e3 #|| t > t_uplift(end)
        #     break
		# end
		
	end

	return riverOut, terraceOut, tOut, hseaOut
end

# ╔═╡ 68d55cd8-0af4-11eb-0dc7-bbd95dd10d45
timefinder(time::Vector{Float64},t::Number) = argmin(abs.(time.-float(t)))

# ╔═╡ 23b62cde-0ae5-11eb-05b1-e10c2d478e23
river_z,terrace_z, time, h_sea = mainsolver(sea_age,sea_lvl_curve,t,dt,yr,
		x,sol0,terrace0,
		t_uplift,uplift_background,U0,uplift_type,
		beta_z,P0,h_wb,imalatesta,
		n,r,dx)

# ╔═╡ 16d9a0de-0ae9-11eb-1d1f-07ed9a98cc31
begin	
	t1 = time[end]/1e3
	
	hline([h_sea[end]],
		color=:blue)
	
	plot!(x/1e3, river_z[end],
		color=:red,
		label="river")
	
	plot!(x/1e3, terrace_z[end],
		color=:green,
		label="terrace")

	xlims!(75,150)
	ylims!(-200,400)
	title!("time $t1 kyrs")
end

# ╔═╡ 7aa3489e-0be4-11eb-1860-ab0c09bee274
begin
	function firstderivative(z,x) # 1st derivative
	    dz = (z[2] - z[1])*2
	    ∇z = similar(x)
	
	    @inbounds for i ∈ 2:length(x)-1 
	        ∇z[i] = (x[i+1] - x[i-1]) / dz
	    end
	
	    ∇z[1]   = (x[2] - x[1]) / (dz)
	    ∇z[end] = (x[end] - x[end-1]) / (dz)
	    
	    return ∇z
	end
	
	function secondderivative(z,x) # 2nd derivative
	    dz = 2*(z[2] - z[1])^2
	    Δz = similar(x)
	
	    @inbounds for i ∈ 2:length(x)-1 
	        Δz[i] = (x[i+1] - 2*x[i] + x[i-1]) / dz
	    end
	
	    Δz[1]   = Δz[2]
	    Δz[end] = Δz[end-1]
	
	    return Δz
	end
end

# ╔═╡ 396c043c-0be5-11eb-0b49-2be9189080da
begin
	dz=firstderivative(x,river_z[end])
	dz_terrace=firstderivative(x,terrace_z[end])
	dzz=secondderivative(x,river_z[end])
	dzz_terrace=secondderivative(x,terrace_z[end])
end

# ╔═╡ 85730772-0be5-11eb-2861-8bb8e90bc525
begin	
	plot(x/1e3, dz_terrace,
		color=:blue,
		label="terrace")

	plot!(x/1e3, dz,
		color=:red,
		label="river")

	xlims!(75,150)
	ylims!(-0.1,0)
	title!("time $t1 kyrs")
end

# ╔═╡ de4afc50-0be6-11eb-0f3d-c343541ee0f7
begin	
	plot(x/1e3, dzz_terrace./maximum(dzz_terrace),
		color=:blue,
		label="terrace")

	plot!(x/1e3, dzz./maximum(dzz),
		color=:red,
		label="river")

	xlims!(90,130)
	ylims!(-1,1)
	title!("time $t1 kyrs")
end

# ╔═╡ a5b0d85a-0af7-11eb-0d86-d372e97c147e
begin
	it = timefinder(time, 100e3)
	t2 = time[it]/1e3
	
	hline([h_sea[it]],
		color=:blue)
	
	plot!(x/1e3, river_z[it],
		color=:red,
		label="river")
	
	plot!(x/1e3, terrace_z[it],
		color=:green,
		label="terrace")

	xlims!(75,175)
	ylims!(-200,1e3)
	title!("time $t2 kyrs")
end

# ╔═╡ 5407d816-0b4f-11eb-0199-d5d5a2cd192a
# function instanterosion(river_z,x)
# 	volume = Vector{Float64}(undef,length(river_z))

# 	Threads.@threads for i ∈ 2:length(river_z)
# 		volume[i] = erodedvolume(river_z[i-1],river_z[i],x)
# 	end

# 	volume[1] = volume[2]
	
# 	return volume
	
# end

# ╔═╡ 666ec31c-0b4e-11eb-11dc-4518cca1359c
# begin
# 	plot(time/1e3, erosionOut,
# 		color=:red)
	
# 	title!("eroded volume")
# end

# ╔═╡ d6e99588-0afc-11eb-141c-cbd637033019
# make animation (it's quite slow, consider comenting it down)
anim = @animate for it ∈ 1:5:length(time)
	t2 = time[it]/1e3
	
	hline([h_sea[it]],
		color=:blue,
		label="sea_level")
	
	plot!(x/1e3, river_z[it],
		color=:red,
		label="river")
	
	plot!(x/1e3, terrace_z[it],
		color=:green,
		label="terrace")

	xlims!(75,175)
	ylims!(-200,1e3)
	title!("time $t2 kyrs")
end 

# ╔═╡ dfc32998-0afd-11eb-0d19-4d8729c6870f
gif(anim, "anim_fps15.gif", fps = 1000)

# ╔═╡ Cell order:
# ╠═75f53f00-0ae0-11eb-2523-e988a6f2e181
# ╠═8a912dc0-0ae2-11eb-03f9-0f451ad8b10c
# ╠═c55008b2-0ae2-11eb-05fe-35f20d380ea4
# ╟─09195ea4-0ae3-11eb-08e4-c99929247d28
# ╠═26c54580-0ae3-11eb-3263-f176fad834ee
# ╠═457170e4-0ae3-11eb-149e-6d96ff1ea287
# ╟─5740103c-0ae3-11eb-31d5-d795e058982b
# ╟─97929f92-0ae3-11eb-119b-71390eb9e1d3
# ╟─c602f700-0ae3-11eb-3005-9d5af1882906
# ╟─e865764a-0ae7-11eb-200c-e7fc36fd90d3
# ╟─35efe2a8-0ae4-11eb-393f-3be330f265da
# ╟─68d55cd8-0af4-11eb-0dc7-bbd95dd10d45
# ╟─23b62cde-0ae5-11eb-05b1-e10c2d478e23
# ╟─16d9a0de-0ae9-11eb-1d1f-07ed9a98cc31
# ╟─7aa3489e-0be4-11eb-1860-ab0c09bee274
# ╟─396c043c-0be5-11eb-0b49-2be9189080da
# ╟─85730772-0be5-11eb-2861-8bb8e90bc525
# ╟─de4afc50-0be6-11eb-0f3d-c343541ee0f7
# ╠═a5b0d85a-0af7-11eb-0d86-d372e97c147e
# ╟─5407d816-0b4f-11eb-0199-d5d5a2cd192a
# ╟─666ec31c-0b4e-11eb-11dc-4518cca1359c
# ╠═d6e99588-0afc-11eb-141c-cbd637033019
# ╠═dfc32998-0afd-11eb-0d19-4d8729c6870f
