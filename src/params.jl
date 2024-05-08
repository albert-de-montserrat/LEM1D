abstract type AbstractParams end

struct TerraceParams{T}
    h_wb::T
    βz::T
    βx::T
    Poff::T
    P0::T
    # buffer arrays
    sol_n::Vector{T}
    sol_old::Vector{T}

    function TerraceParams(;
        h_wb::T = 150e0,
        βz::T   = 1.0e-5,
        βx::T   = 5e-5,
        Poff::T = 5e-2,
        P0::T   = 5e-7,
        n       = 10
    ) where T
        
        new{T}(h_wb::T, βz, βx, Poff, P0, zeros(n), zeros(n))
    end
end