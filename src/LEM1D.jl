module LEM1D
    using LinearAlgebra,
        Interpolations, DelimitedFiles, LoopVectorization, Printf
    
    include("read_txt_files.jl")
    export get_sea_level, read_uplift

    include("grid.jl")
    export TerraceProfile
    
    include("terraces.jl")
    export terracenewton!, terrace_retreat!, find_shore_id
    
    include("others.jl")
    export variable_recurrence_time_fixed, poisson_time_serie, update_reoccupation!
end