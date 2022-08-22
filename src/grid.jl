struct TerraceProfile{T,I,F}
    x::T
    z::T
    nnod::I
    Δx::F

    function TerraceProfile(L::T, Δx::T, slope; intercept=2e3) where T
        x       = collect(zero(T):Δx:L)
        z       = @. slope * x + intercept
        nn      = length(x)
        new{Vector{T},Int64,T}(x,z,nn,Δx)
    end

end