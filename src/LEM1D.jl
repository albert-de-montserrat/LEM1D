module LEM1D
    using LinearAlgebra,
        Interpolations, DelimitedFiles, LoopVectorization, Printf
    
    using ForwardDiff

    @inline function f_df(f::F, x::R) where {F,R<:Real}
        T = typeof(ForwardDiff.Tag(f, R))
        df = f(ForwardDiff.Dual{T}(x, one(x)))
        return df.value, first(df.partials)
    end    

    include("read_txt_files.jl")
    export get_sea_level, read_uplift

    include("grid.jl")
    export TerraceProfile
    
    include("terraces.jl")
    export terrace_vertical_erosion!, terrace_retreat!, find_shore_id, terrace_vertical_explicit!
    
    include("others.jl")
    export variable_recurrence_time_fixed, poisson_time_serie, update_reoccupation!
end