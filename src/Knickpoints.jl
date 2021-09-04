fluxlimiter(::Type{FluxLimiter{minmod}}, rf)      = max(0.,min(1.,rf))
fluxlimiter(::Type{FluxLimiter{superbee}}, rf)    = max(0.,min(1.,2rf),min(rf,2.))
fluxlimiter(::Type{FluxLimiter{vanAlbada}}, rf)   = (rf^2+rf)/(rf^2+1)
fluxlimiter(::Type{FluxLimiter{vanAlbada2}}, rf)  = (2rf)/(rf^2+1)
fluxlimiter(::Type{FluxLimiter{vanLeer}}, rf)     = (rf+abs(rf))/(abs(rf)+1)
fluxlimiter(::Type{FluxLimiter{ospre}}, rf)       = 1.5(rf^2+rf)/(rf^2+rf+1)
fluxlimiter(::Type{FluxLimiter{hquick}}, rf)      = 2.0(rf+abs(rf))/(rf+2)
fluxlimiter(::Type{FluxLimiter{koren}}, rf)       = max(0,min(2rf,min((1+2rf)/3,2)))


@inline function updater!(RiverPhysics, RiverArrays, River)
    # unpack
    Kr, kappa_a, h, m = RiverPhysics.Kr, RiverPhysics.kappa_a,  RiverPhysics.h, RiverPhysics.m
    L, r = RiverArrays.L, RiverArrays.r
    dx, z, x = River.Δx^2, River.z, River.x

    # constant values
    cte1 = -Kr*(kappa_a^m)/dx
    cte2 = h*m
    # update river length
    river_length!(L, River)
    @tturbo for i ∈ 2:length(x)
        r[i] = cte1 * L[i]^(cte2) 
    end
end

function river_length!(L, River)
    x, z, dx = River.x, River.z, River.Δx^2
    @tturbo for i in 2:length(x)
        L[i] = √(dx  + (z[i]-z[i-1])^2) 
    end
    @inbounds @fastmath for i in 2:length(x)
        L[i] += L[i-1]
    end
end

@inline function rfactor!(rf, sol, ilast)
    """
    Wikipedia definition, paper has a slightly different one
    """
    @inbounds for i in 1:ilast
        s1 =  sol[i+1] - sol[i]
        if s1 == 0.0
            rf[i] = 1.0
        elseif s1 <= 0.0
            rf[i] = 0.0
        else
            rf[i] = (sol[i+2] - sol[i+1]) / s1
        end
	end
    @inbounds rf[ilast] = 1.0
    @inbounds rf[ilast-1] = 1.0
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

    @fastmath for i in 2:ilast
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
        @inbounds multiplier  = ((sol0[i] - sol0[i+1])*invdx)^(n-1) # |-> pow32 speeds up the whole solver by 60%
        # update solution
        @inbounds sol[i] = max(sol[i] - invdxdt*(Fᵣ-Fₗ) * multiplier, sol[i+1])
        # @inbounds sol[i] -= invdxdt*(Fᵣ-Fₗ) * multiplier
    end
end