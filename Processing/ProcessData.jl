using HDF5, Makie
import LinearAlgebra: norm

function getparameters(folder) 
    files   = filter(x->endswith(x, ".h5"), readdir(folder))   
    
    fid     = h5open(joinpath(folder,files[110]), "r")
    zref    = read(fid, "River/z")
    zref  .-= zref[1]

    U       = Array{Float64}(undef,length(files))
    βz      = similar(U)
    rms     = similar(U)
    @inbounds for i ∈ 1:length(files)
        fid     = h5open(joinpath(folder,files[i]), "r")
        U[i]    = read(fid, "Parameters/U")
        βz[i]   = read(fid, "Parameters/beta_z")

        z       = read(fid, "River/z")
        z     .-= z[1]
        rms[i]  = norm(zref.-z)
    end

    return U,βz,rms
end

U,βz,rms = getparameters(folder)

R = rms' .* ones(length(U))

scatter(U,rms)
scatter(βz,rms)
heatmap(U,βz,R)


fid     = h5open("output/Test_1.h5", "r")
z    = read(fid, "River/z")
x   = read(fid, "River/x")
