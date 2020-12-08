function terrace(Terrace,βz,P0,h_wb,
    sea_age,U0,fsea,
    K,M,iMC,jMC,D,
    e2n,nel,
    Δx,Δt, yr)

    buffer1,buffer2 = similar(Terrace.z),similar(Terrace.z)
    nit             = Int64(floor(maximum(sea_age .- Δt) / Δt))
    nn              = length(buffer1)
    t               = 0.0
    for it = 1:nit

        # GET SEA LVL AND UPLIFT RATE =========================================
        h_sea   = fsea(t)
        U0
        # =====================================================================

        # GET COORDS BELOW SEA LEVEL===========================================
        id_shore_terrace    = find_shore_id(h_sea, Terrace.z, nn)
        # =====================================================================

        # TERRACES SOLVER := MALATESTA ========================================
        terracenewton!( βz,P0,h_wb,Δt * yr,Terrace.z,h_sea,id_shore_terrace+1,nn,
                        buffer1,buffer2)
        # =====================================================================

        # Diffusion =======================================================
        # Terrace.z = femsolver(Δt,Δx,e2n,nel,nn, Terrace, h_sea,
        #     K,M,iMC,jMC,
        #     D) 
        # ==================================================================

        t +=  Δt

        if t > maximum(sea_age.-Δt) #|| t > t_uplift(end)
            break
        end
        
    end
    return Terrace
end

#-------------------------------------------------------------------------
function river(River,βz,P0,h_wb,
    n,φ,r,rf, Kr, kappa_a, L, h, m,
    sea_age,U0,fsea,
    K,M,iMC,jMC,D,
    e2n,nel,
    Δx,Δt, yr)

    buffer1,buffer2 = similar(River.z),similar(River.z)
    nit             = Int64(floor(maximum(sea_age .- Δt) / Δt))
    nn              = length(buffer1)
    t               = 0
    for it = 1:nit

        # GET SEA LVL AND UPLIFT RATE =========================================
        h_sea   = fsea(t)
        U       = U0
        # =====================================================================

        # GET COORDS BELOW SEA LEVEL===========================================
        id_shore  = find_shore_id(h_sea, River.z, nn)
        # =====================================================================

        # TERRACES SOLVER := MALATESTA ========================================
        terracenewton!( βz,P0,h_wb,Δt * yr,River.z,h_sea,id_shore+1,nn,
                        buffer1,buffer2)
        # =====================================================================

        # EXPLICIT ADVECTION -> FLUX LIMITER TVD =============================
        updater!(r, Kr, kappa_a, L, h, m, Δx^2, River.x, River.z, nn)
        TVD!(River.z,buffer1,n,φ,r,rf,Δt, Δx,id_shore)
        # =====================================================================

        t +=  Δt

        if t > maximum(sea_age.-Δt) #|| t > t_uplift(end)
            break
        end
        
    end
    return River
end
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
function terracenewton!(βz::Float64,P::Float64,h_wb::Float64,dt::Float64,
    sol::Array{Float64},h_sea::Float64,i1::Int64,i2::Int64,sol_n,sol_old)
    """
    Add derivation of newton-raphson
    """

    # sol_old    .= copy(sol)
    cte1        = βz * P * dt
    cte2        = 4.0/h_wb
    cte3        = cte1 * cte2
    # sol_n      .= copy(sol) .+1e3   
    copyprofile!(sol_old,sol)
    residual    = 1.0
    while residual > 1e-4
        copyprofile!(sol_n,sol)
        @avx for ii = i1:i2
            cte4    = exp(-cte2*(h_sea-sol[ii]))
            sol[ii]-= (sol[ii] + cte1 * cte4 - sol_old[ii]) / (1.0 + cte3 * cte4)
        end       
        # sol_n      .= copy(sol)
        residual = norm2(sol,sol_n)
    end
end

function terracenewton!(βz::Float64,P::Float64,h_wb::Float64,dt::Float64,
    sol::CuArray{Float32,1},h_sea::Float64,i1::Int64,i2::Int64)

    sol_old     = copy(sol)
    cte1        = βz * P * dt
    cte2        = 4.0/h_wb
    cte3        = cte1 * cte2
    sol_n       = copy(sol) .+1e3   

    while norm(sol.-sol_n) < 1e-4       
        sol_n   = copy(sol)
        # cte4    = exp(-cte2*(h_sea-sol))
        # sol    -= (sol + cte1 * cte4 - sol_old) / (1.0 + cte3 * cte4)      
        cte4        = exp(-cte2*(h_sea-sol[i1:i2]))
        sol[i1:i2] .-= (sol[i1:i2] + cte1 * cte4 - sol_old[i1:i2]) / (1.0 + cte3 * cte4)
    end
end

# SOLVE WAVE EROSION ===========================================================
function solve_terrace!(beta_z::Float64,P::Float64,h_wb::Float64,dt::Float64,
    sol::Array{Float64},h_sea::Float64,i1::Int64,i2::Int64)
    @inbounds @fastmath @simd  for ii = i1:i2
        sol[ii] +=  - beta_z * P* exp(-4*(h_sea-sol[ii])/h_wb) * dt
    end

    return sol
end
# ==============================================================================

# PARALLEL - SOLVE WAVE EROSION ================================================
function parallel_solve_terrace!(beta_z::Float64,P::Float64,h_wb::Float64,dt::Float64,
    sol::Array{Float64},h_sea::Float64,i1::Int64,i2::Int64)
    i1 = id_shore[1]+1
    i2 = nn
    Threads.@threads for ii = i1:i2
        @inbounds sol[ii] -=  beta_z * P * exp(-4*(h_sea-sol[ii])/h_wb) * dt
    end

    # return sol
end
# ==============================================================================

# SOLVE ADVECTION ==============================================================
function solve_advection!(sol::Array{Float64,1},n::Int64,r::Array{Float64,1},
    dx::Float64,dt::Float64,i1::Int64,i2::Int64)

    @inbounds @fastmath @simd for ii ∈ i1:i2
        sol[ii] += dt * r[ii] * 1.0 * ((abs(sol[ii+1] - sol[ii])/dx))
    end

    return sol
end
# ==============================================================================

# SOLVE ADVECTION ==============================================================
function solve_advection_parallel!(sol::Array{Float64,1},n::Int64,r::Array{Float64,2},
    dx::Float64,dt::Float64,i1::Int64,i2::Int64,bc_id::Array{Int32,2})

    dummy   = Array{Int32}(undef,length(sol)-length(bc_id),1)
    c       = 1

    @inbounds for ii = 1:length(bc_id)-1
        for jj = bc_id[ii]+1 : bc_id[ii+1]-1
            dummy[c] = jj
            c += 1
        end
    end

    Threads.@threads  for ii in dummy
        @inbounds @fastmath sol[ii] += dt * r[ii] * 1 * ((abs(sol[ii+1] - sol[ii])/dx)^n)
    end

    return sol
end
# ==============================================================================

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
        if x[i] >= x0
            x[i] -= x0 - 0.01*i
        end 
    end    
end

#==============================================================================
FLUX LIMITER ADVECTION FUNCTIONS
==============================================================================#

# ====================================================================
@inline function rfactor!(rf,sol)
    """
    This is Wikipedias definition, paper has a slightly different one
    """
    for i ∈ 1:(length(sol)-2)
        # s1 =  sol[i+1] - sol[i]
        # s1 == 0.0 ? rf[i]= 1.0 : rf[i] = (sol[i] - sol[i-1]) / s1
        s1 =  sol[i+1] - sol[i]
        s1 == 0.0 ? rf[i]= 1.0 : rf[i] = (sol[i+2] - sol[i+1]) / s1

    end
    rf[end]     = 1.0
    rf[end-1]   = 1.0

end

@inline minmod(rf)      = max(0.,min(1.,rf))
@inline superbee(rf)    = max(0.,min(1.,2rf),min(rf,2.))
@inline vanAlbada(rf)   = (rf^2+rf)/(rf^2+1)
@inline vanAlbada2(rf)  = (2rf)/(rf^2+1)
@inline vanLeer(rf)     = (rf+abs(rf))/(abs(rf)+1)
@inline ospre(rf)       = 1.5(rf^2+rf)/(rf^2+rf+1)
@inline hquick(rf)      = 2.0(rf+abs(rf))/(rf+2)
@inline koren(rf)       = max(0,min(2rf,min((1+2rf)/3,2)))

function selectfluxlimiter(type)
    if type == "minmod"
        f = minmod
    elseif type == "superbee"
        f = superbee
    elseif type == "vanAlbada"
        f = vanAlbada
    elseif type == "vanAlbada2"
        f = vanAlbada2
    elseif type == "vanLeer"
        f = vanLeer
    elseif type == "ospre"
        f = ospre
    elseif type == "hquick"
        f = hquick
    elseif type == "koren"
        f = koren
    end    
end

@inline function TVD!(sol,sol0,n,φ,r,rf,dt,dx,ilast)
    rfactor!(rf,sol)
    copyprofile!(sol0,sol)
    @inbounds for i ∈ 2:ilast
        # allocations
        soli        = sol0[i]
        solplus     = sol0[i+1]
        solminus    = sol0[i-1]
        ri          = r[i]
        # solve
        α₀,α₁       = 0.5 * (1+sign(ri))*ri, 0.5 * (1-sign(ri))*ri
        Frlo        = α₀*soli + α₁*solplus    # f right low
        Fllo        = α₀*solminus + α₁*soli   # f left  low
        Frhi        = 0.5*ri*(soli + solplus) - (dt/2/dx)*(solplus-soli)*ri^2
        Flhi        = 0.5*ri*(soli + solminus) - (dt/2/dx)*(soli-solminus)*ri^2
        Fᵣ          = Frlo+φ(rf[i]) * (Frhi - Frlo )
        Fₗ          = Fllo+φ(rf[i-1]) * (Flhi - Fllo )
        # non linear term
        # multiplier  = ((sol0[i] - sol0[i+1])/dx)^(n-1)
        multiplier  = pow32((sol0[i] - sol0[i+1])/dx, n-1) # |-> speeds up the whole solver by 60%
        # update solution
        sol[i]     -= dt*(Fᵣ-Fₗ)/dx * multiplier

    end
end

@inline function copyprofile!(dest,profile)
    for i ∈ axes(profile,1)
        @inbounds dest[i] = profile[i]
    end
end

