abstract type LimiterFunction end 
abstract type minmod <:LimiterFunction end 
abstract type superbee <:LimiterFunction end 
abstract type vanAlbada <:LimiterFunction end 
abstract type vanAlbada2 <:LimiterFunction end 
abstract type vanLeer <:LimiterFunction end 
abstract type ospre <:LimiterFunction end 
abstract type hquick <:LimiterFunction end 
abstract type koren <:LimiterFunction end 

struct FluxLimiter{T<:LimiterFunction} end

struct KnickPointsParameters{T, N}
    Kr::T
    h::T
    kappa_a::T
    m::T
    n::T
    φ::Type{FluxLimiter{N}}

    function KnickPointsParameters{T, N}(;
        Kr = 5e-7,    # whipple and tucker 1999_ dimensional coeff of erosion
        h = 1.92,    # whipple and tucker 1999
        kappa_a = 4.6071,  # whipple and tucker 1999_ area-length coeff
        m = 1.0,
        n = 1.1,
        φ = FluxLimiter{vanAlbada},
        ) where {T,N}

        new{T, N}(
            Kr,
            h,
            kappa_a,
            m,
            n,
            φ,
        )

    end

end

struct KnickPointsArrays{T}
    L::T
    r::T
    rf::T
    
    function KnickPointsArrays{T}(River, nn, RiverPhysics) where T
        x, z, Δx = River.x, River.x, River.Δx
        L = riverlength(x, z, nn)
        r = @. -RiverPhysics.Kr * (RiverPhysics.kappa_a^RiverPhysics.m) * (L^(RiverPhysics.h * RiverPhysics.m)) / Δx # EQUATION    
        rf = zero(z)
        new{T}(L,r,rf)
    end

end

struct TerraceParameters{T}
    βz::T
    h_wb::T
    P0::T
    P_off::T

    function TerraceParameters{T}(;
        βz = 1e-6,
        h_wb = 15.0,
        P0 = 5e-5,
        P_off = 5e-2,
    ) where T
        new{T}(
            βz,
            h_wb,
            P0,
            P_off,
        )
    end
end
   
struct HillSlopeParameters{T}
    κa::T
    diffusion::T
    alpha_dif::T

    function HillSlopeParameters{T}(;
        κa = 5e-5,
        diffusion = 0.1 ,
        alpha_dif = 1.0 ,
    ) where T
        new{T}(
            κa,
            diffusion,
            alpha_dif,
        )
    end

end

struct SubmarineParameters{T}
    K_s::T
    λ::T

    function SubmarineParameters{T}(;
        K_s = 1e2,
        λ = 5e-4,
    ) where T
        new{T}(
            K_s,
            λ,
        )
    end

end