fluxlimiter(::Type{FluxLimiter{minmod}}, rf)      = max(0.,min(1.,rf))
fluxlimiter(::Type{FluxLimiter{superbee}}, rf)    = max(0.,min(1.,2rf),min(rf,2.))
fluxlimiter(::Type{FluxLimiter{vanAlbada}}, rf)   = (rf^2+rf)/(rf^2+1)
fluxlimiter(::Type{FluxLimiter{vanAlbada2}}, rf)  = (2rf)/(rf^2+1)
fluxlimiter(::Type{FluxLimiter{vanLeer}}, rf)     = (rf+abs(rf))/(abs(rf)+1)
fluxlimiter(::Type{FluxLimiter{ospre}}, rf)       = 1.5(rf^2+rf)/(rf^2+rf+1)
fluxlimiter(::Type{FluxLimiter{hquick}}, rf)      = 2.0(rf+abs(rf))/(rf+2)
fluxlimiter(::Type{FluxLimiter{koren}}, rf)       = max(0,min(2rf,min((1+2rf)/3,2)))

@inline function rfactor!(rf, sol, ilast)
    """
    This is Wikipedias definition, paper has a slightly different one
    """
    @inbounds Threads.@threads for i in 1:ilast
       s1 =  sol[i+1] - sol[i]
       rf[i] = s1 == 0.0 ? 1.0 : (sol[i+2] - sol[i+1]) * (1/s1)
	end
    @inbounds rf[ilast] = 1.0
    @inbounds rf[ilast-1] = 1.0
end

@inline function updater!(RiverPhysics, RiverArrays, River)
    # unpack
    Kr, kappa_a, h, m = RiverPhysics.Kr, RiverPhysics.kappa_a,  RiverPhysics.h, RiverPhysics.m
    L, r = RiverArrays.L, RiverArrays.r
    dx, z, nn = River.Δx^2, River.z, River.nn
    # constant values
    cte1 = -Kr*(kappa_a^m)/dx
    cte2 = h*m
    # tight loop
    Threads.@threads for i in 2:nn
        @inbounds L[i] = L[i-1] + pow32(dx  + (z[i]-z[i-1])^2, 0.5) 
        @inbounds r[i] = cte1 * pow32(L[i], cte2)
    end
end

function TVD!(River, sol0, RiverPhysics::KnickPointsParameters, RiverArrays::KnickPointsArrays, dt, ilast)
    # Add reference to paper of Govers and Schwangart
    
    n, φ = RiverPhysics.n, RiverPhysics.φ
    r, rf = RiverArrays.r, RiverArrays.rf
    sol, dx = River.z, River.Δx

    rfactor!(rf, sol, ilast)
    copyto!(sol0, sol)
    # cache some constant values to avoid redundance
    invdx = 1/dx
    invdxdt = dt*invdx
    half_dt = dt*0.5

    for i in 2:ilast
        # cache
        @inbounds soli, solplus, solminus = sol0[i], sol0[i+1], sol0[i-1]
        @inbounds ri = r[i]
        ri_half = ri*0.5
        dtdx = half_dt*invdx*ri*ri
        # solve
        α₀, α₁ = ri_half * (1+sign(ri)), ri_half*(1-sign(ri))
        Frlo = α₀*soli + α₁*solplus    # f right low
        Fllo = α₀*solminus + α₁*soli   # f left  low
        Frhi = ri_half*(soli + solplus) - dtdx*(solplus-soli)
        Flhi = ri_half*(soli + solminus) - dtdx*(soli-solminus)
        @inbounds Fᵣ = Frlo+fluxlimiter(φ, rf[i]) * (Frhi - Frlo )
        @inbounds Fₗ = Fllo+fluxlimiter(φ, rf[i-1]) * (Flhi - Fllo )
        # non linear term
        @inbounds multiplier  = pow32((sol0[i] - sol0[i+1])*invdx, n-1) # |-> pow32 speeds up the whole solver by 60%
        # update solution
        @inbounds sol[i] -= invdxdt*(Fᵣ-Fₗ) * multiplier
    end
end