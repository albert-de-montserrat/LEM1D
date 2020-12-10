# module moduleFEM

# using StaticArrays,SparseArrays,LinearAlgebra

# export meshgen1d,matrixassembly,dirichlet

abstract type RiverProfile end
abstract type TerraceProfile end

mutable struct Profile{T}
    x::Vector{Float64}
    z::Vector{Float64}
end

# mutable struct CuProfile{T}
#     x::CuArray{Float32,1}
#     z::CuArray{Float32,1}
# end

function femsolver(Δt,Δx,e2n,nel,nnod, River::Profile{RiverProfile},hsea,
    K,M,iMC,jMC,
    κa, alpha_dif, diffusion,
    K_s,λ)

    x           = River.x
    z           = River.z

    # Diffusion equation ---------------------------
    isea        = z.>hsea
    imarine     = isea.==false
    L           = riverlength(x[isea],z[isea],sum(isea))
    κ_aerial    = @. diffusion + κa * alpha_dif * L
    # ----------------------------------------------

    hw          = @. abs(hsea-z)
    κ           = @. K_s*exp(-λ*hw)
    κ[isea]    .= κ_aerial
    α           = 1/2
        #  alpha = 1/2 --> Crank Nicholson (IMPLICIT)
        #  alpha = 1   --> backward difference (FULL IMPLICIT) 
        #  alpha = 0   --> forward difference (EXPLICIT) // not stable for now
    TOL         = 1e-3
    bc_rhs      = fill(0.0,nnod)
    T0          = copy(z)
    T           = copy(z)
    Tdummy      = copy(T)

    bcdof       = [1, nnod]
    bcval       = [z[1], z[nnod]]
    ifree       = 2:nnod-1

    for it = 1:40
        # ASSEMBLE MATRICES ---------------------------------------------------
        fill!(K,0.0)
        fill!(M,0.0)
        K,M         = matrixassembly(e2n,nnod, nel, κ, Δx, K, M, iMC,jMC)
        # K,M         = matrixassembly_threaded(e2n,nnod, nel, κ, Δx, K, M, iMC,jMC)
        KK          = @. Δt*α*K + M

        # DIRICHLET BOUNDARY CONDITIONS ----------------------------------------
        fill!(bc_rhs,0.0)
        KK,bc_rhs   = dirichlet(KK,bc_rhs,bcdof,bcval)
    
        # RIGHT HAND SIDE VECTOR -----------------------------------------------
        rhs         = (@. - Δt*(1-α)*K + M) * T0 - bc_rhs
        # rhs[bcdof] .= bcval    
        
        # SOLVE SYSTEM ---------------------------------------------------------        
        # T           = KK\rhs
        
        if @isdefined F
            factorize!(F,KK[ifree,ifree])
        else
            F = factorize(KK[ifree,ifree])
        end
        T[ifree] = F\rhs[ifree]
        
        # # UPDATE SUBMARINE HEIGHT ------------------------------------------------------
        @views begin
            hw          = hsea.-T[imarine]
            κ[imarine] .= @. K_s*exp(-λ*hw)
            κ[isea]    .= κ_aerial
        end
        if it > 1 && norm2(T,Tdummy) < TOL;break;end
        Tdummy      = copy(T)
    end

    return T

end

function femsolver(Δt,Δx,e2n,nel, nn,Terrace::Profile{TerraceProfile},hsea,
    K,M,iMC,jMC,
    D)

    isea        = 1:find_shore_id(hsea, Terrace.z, nn)
    x           = Terrace.x
    z           = Terrace.z
    nnod        = length(x)
    
    κ           = fill(0.0,nnod)
    κ[isea]    .= D
    α           = 1/2
        #  alpha = 1/2 --> Crank Nicholson (IMPLICIT)
        #  alpha = 1   --> backward difference (FULL IMPLICIT) 
        #  alpha = 0   --> forward difference (EXPLICIT) // not stable for now
    TOL         = 1e-3
    bc_rhs      = fill(0.0,nnod)
    bcdof       = [1, nnod]
    bcval       = [z[1], z[nnod]]
    ifree       = 2:nnod-1

    fill!(K,0.0)
    fill!(M,0.0)
    K,M         = matrixassembly(e2n,nnod, nel, κ, Δx, K, M, iMC,jMC)
    # K,M         = matrixassembly_threaded(e2n,nnod, nel, κ, Δx, K, M, iMC,jMC)
    KK          = @.(Δt*α*K + M)

    # DIRICHLET BOUNDARY CONDITIONS ----------------------------------------
    fill!(bc_rhs,0.0)
    KK,bc_rhs   = dirichlet(KK,bc_rhs,bcdof,bcval)

    # RIGHT HAND SIDE VECTOR -----------------------------------------------
    rhs         = (@. - Δt*(1-α)*K + M) * z - bc_rhs
    rhs[bcdof] .= bcval    
    
    # SOLVE SYSTEM ---------------------------------------------------------        
    z[isea]    .= KK[isea,isea]\view(rhs,isea)

    return z

end

 ## START meshgen1d FUNCTION ###################################################
function meshgen1d(x0::Float64,z0::Float64,xf::Float64,zf::Float64,
    nnod::Int64)

    x    = Array{Float64}(undef,nnod,1)
    z    = Array{Float64}(undef,nnod,1)
    dx   = (xf-x0) / (nnod-1)
    dz   = (zf-z0) / (nnod-1)
    x[1] = x0;
    z[1] = z0;

    # Create coordinate arrays
    @avx for ii = 2:nnod
        x[ii] = x[ii-1] + dx
        z[ii] = z[ii-1] + dz
    end

    # Create connectivity matrix
    nel = nnod - 1
    e2n = Array{Int64}(undef,2,nel)
    @simd for ii = 1:nel
        e2n[1,ii] = ii
        e2n[2,ii] = ii+1
    end

    # Output
    return x,z,nel,e2n

end #### END meshgen1d FUNCTION ################################################

## START meshgen1d FUNCTION ###################################################
function connectivitymatrix(nnod::Int64)

    # Create connectivity matrix
    nel = nnod - 1
    e2n = Array{Int64}(undef,2,nel)
    @simd for ii = 1:nel
        e2n[1,ii] = ii
        e2n[2,ii] = ii+1
    end

    # Output
    return nel,e2n

end #### END meshgen1d FUNCTION ################################################


## START shapefunctions FUNCTION ###############################################
function shapefunctions(element_order::Int64)

   if element_order == 1
       # local ip coordinates
       ipx = [sqrt(1/3) -sqrt(1/3)]
       w   = 1

    #    N1   = Vector{Float64}(undef,2)
    #    N2   = Vector{Float64}(undef,2)
    #    @inbounds for ii = 1 : 2
    #        N1[ii] = 0.5 * (1 - ipx[ii])
    #        N2[ii] = 0.5 * (1 + ipx[ii])
    #    end
    # ∇N   = [[-1 1]/2,[-1 1]/2]
       N1   = @SVector [0.5 * (1 - ipx[i]) for i ∈ 1:2]
       N2   = @SVector [0.5 * (1 + ipx[i]) for i ∈ 1:2]
       N    = [N1',N2']
       ∇N   = @SVector [-1/2, 1/2]
       ∇N   = [∇N',∇N']
   end

   # Output
   return w,N,∇N

end #### END shapefunctions FUNCTION ###########################################

#### START matrixassembly FUNCTION #############################################
function matrixassembly(e2n::Array{Int64,2},nnod::Int64,
    nel::Int64,dx::Float64)

    # Random variables
    # nnodel  = size(e2n,1)
    nnodel  = 2
    matel   = nnodel * nnodel
    nip     = nnodel    
    κ       = fill(1.0,nnod)
    detJ    = dx/2; # when using  N
    invJ    = 2/dx; # when using ∇N
    kappaEl = Array{Float64}(undef,nnodel,1)
    qEl     = Array{Float64}(undef,nnodel,1)
    q       = fill(0.0,nnod)
    q[1]    = 1.0

    # Calculate K and M bandwith
    bandwith = nnodel + 1
    nfill    = bandwith * nnod - 2 # number of nonzeros

    # Get shape functions
    weight,N,∇N  = shapefunctions(1)

    # Allocate FEM matrices
    # K      = Array{Float64,2}(undef,nnod,nnod) # conductivity matrix
    # M      = Array{Float64,2}(undef,nnod,nnod) # mass matrix
    # Q      = Array{Float64,1}(undef,nnod,1)    # flux vector (i.e. Neumann BC)

    # Allocate FEM element matrices
    Kel     = Array{Float64}(undef,nnodel,nnodel); fill!(Kel,0.0)
    Mel     = Array{Float64}(undef,nnodel,nnodel); fill!(Mel,0.0)
    Qel     = Array{Float64}(undef,nnodel,1); fill!(Qel,0.0)

    # Allocate FEM sparse matrices
    # K       = Array{Float64}(undef,nfill,1) # conductivity matrix
    # M       = Array{Float64}(undef,nfill,1) # mass matrix

    # Temporal allocations
    Kv      = Array{Float64}(undef,nel,nnodel*nnodel);  # storage for conductivity matrix coefficients
    Mv      = Array{Float64}(undef,nel,nnodel*nnodel);  # storage for mass matrix coefficients
    Q_all   = Array{Float64}(undef,nel,nnodel);         # storage for mass matrix coefficients
    iMC     = Array{Int64}(undef,nel,nnodel*nnodel);    # storage for the matrix i indices
    jMC     = Array{Int64}(undef,nel,nnodel*nnodel);    # storage for the matrix j indices

    # Other allocations
    I1      = SVector(1,1)'
    ij      = Array{Int64}(undef,nnodel,1)
    rows    = Array{Int64}(undef,nnodel,nnodel)
    cols    = Array{Int64}(undef,nnodel,nnodel)
    kappaEl = Array{Float64}(undef,2,1) # zeros(MVector{2})
    qEl     = Array{Float64}(undef,2,1) # zeros(MVector{2})

    # Assembly loop
    # Threads.@threads for iel = 1 : nel
    @inbounds for iel = 1 : nel
        # global Kel, Mel,Qel

        # Conductivity and Sediment flux at the nodes of i-th element
        @avx for ii = 1: nnodel
            ij[ii]      = e2n[ii,iel]
            kappaEl[ii] = κ[ij[ii]]
            qEl[ii]     = q[ij[ii]]
        end

        # Numerical integration loop
        for ip = 1 : nip
            # kappa at i-th integration point
            KappaEl  = kappaEl' * N[ip]';
            # kappa at i-th integration point
            QEl      = qEl' * N[ip]';
            # Conductivity elemental matrix
            Kel     += ∇N[ip]' * KappaEl .* ∇N[ip] .* weight .* invJ
            # Mass elemental matrix
            Mel     += N[ip]' *  N[ip] .* weight .* detJ
            # Sediment input force vector
            Qel     += N[ip]' .* qEl;
        end

        rows          .= ij * I1
        cols          .= rows'
        Kv[iel,:]      = Kel[:]
        Mv[iel,:]      = Mel[:]
        Q_all[iel,:]   = Qel[:]
        iMC[iel,:]     = rows[:]
        jMC[iel,:]     = cols[:]

        # c = 1
        # for i ∈ 1:nnodel,j ∈ 1:nnodel
        #     # rows[c]     = EL2NOD[j,iel]    # row location in global matrices
        #     # cols[c]     = EL2NOD[i,iel]    # row location in global matrices
        #     Kv[c,iel]   = Kel[j,i]
        #     Mv[c,iel]   = Mel[j,i]
        #     Q_all[c,iel]= Qel[j,i]
        #     iMC[c,iel]  = e2n[j,iel]
        #     jMC[c,iel]  = e2n[i,iel]
        #     c       +=1
        # end

        # Reset element matrices
        fill!(Kel , 0.0)
        fill!(Mel , 0.0)
        fill!(Qel , 0.0)

    end

    M   = sparse(iMC[:],jMC[:],Mv[:])
    K   = sparse(iMC[:],jMC[:],Kv[:])
    # Q   = accumarray(e2n[:], Q_all[:])
    Q   = accumQ(e2n[:], Q_all[:] , nnod)

    return K,M,Q

end #### END matrixassembly FUNCTION ###########################################

#### START matrixassembly FUNCTION #############################################
function matrixassembly(e2n::Array{Int64,2},nnod::Int64,
    nel::Int64,κ,dx::Float64,K,M,iMC,jMC)

    # Random variables
    # nnodel  = size(e2n,1)
    nnodel  = 2
    matel   = nnodel * nnodel
    nip     = nnodel    
    # κ       = fill(1.0,nnod)
    detJ    = dx/2; # when using  N
    invJ    = 2/dx; # when using ∇N
    kappaEl = Array{Float64}(undef,nnodel,1)
    # qEl     = Array{Float64}(undef,nnodel,1)
    # q       = fill(0.0,nnod)
    # q[1]    = 1.0

    # Get shape functions
    weight,N,∇N  = shapefunctions(1)

    # Allocate FEM element matrices
    Kel     = Array{Float64}(undef,nnodel,nnodel); fill!(Kel,0.0)
    # Mel     = Array{Float64}(undef,nnodel,nnodel); fill!(Mel,0.0)
    # Qel     = Array{Float64}(undef,nnodel,1); fill!(Qel,0.0)

    # Temporal allocations
    # Q_all   = Array{Float64}(undef,nel,nnodel);         # storage for mass matrix coefficients

    # Other allocations
    kappaEl = Array{Float64}(undef,2,1) # zeros(MVector{2})
    # qEl     = Array{Float64}(undef,2,1) # zeros(MVector{2})
    KappaEl = Array{Float64}(undef,1,1) # zeros(MVector{2})
    # QEl     = Array{Float64}(undef,1,1) # zeros(MVector{2})
    ij      = Array{Int64}(undef,nnodel,1)
    
    A       = [∇N[ip]' * ∇N[ip] * weight * invJ for ip = 1:2]
    B       = [N[ip]'  * N[ip]  * weight * detJ for ip = 1:2]
    Mel     = B[1]+B[2]

    # Assembly loop
    @inbounds for iel = 1 : nel

        # Conductivity and Sediment flux at the nodes of i-th element
        # @simd for ii = 1: nnodel
        #     ij[ii]      = e2n[ii,iel]
        #     kappaEl[ii] = κ[ij[ii]]
        #     # qEl[ii]     = q[ij[ii]]
        # end
        
        kappaEl = @SVector [κ[e2n[ii,iel]] for ii = 1 : 2] 

        # Numerical integration loop
        @views for ip = 1 : nip
            # kappa at i-th integration point
            # mul!(KappaEl, N[ip], kappaEl)
            # # kappa at i-th integration point
            # mul!(QEl, qEl', N[ip]')
            # Conductivity elemental matrix
            # @. Kel  += A[ip] * KappaEl[1]     
            @. Kel  += A[ip] * (N[ip] * kappaEl)
            # Mass elemental matrix
            # @. Mel  += B[ip]
            # # Sediment input force vector
            # @. Qel  += A[ip]' .* qEl
        end

        for ii ∈ 1:nnodel*nnodel
            i       = iMC[ii,iel]
            j       = jMC[ii,iel]
            K[i,j] += Kel[ii]
            M[i,j] += Mel[ii]
            # @show Mel[ii]
            # if ii < nnodel +1
            #     Q_all[iel,ii] = Qel[ii]
            # end
        end
       
        # Reset element matrices
        fill!(Kel , 0.0)
        # fill!(Mel , 0.0)
        # fill!(Qel , 0.0)

    end

    # Q   = accumQ(e2n[:], Q_all[:] , nnod)
    
    return K,M

end #### END matrixassembly FUNCTION ###########################################

#### START matrixassembly FUNCTION #############################################
function matrixassembly_threaded(e2n::Array{Int64,2},nnod::Int64,
    nel::Int64,κ,dx::Float64,K,M,iMC,jMC)

    # Random variables
    nnodel  = 2
    matel   = nnodel * nnodel
    nip     = nnodel    
    detJ    = dx/2; # when using  N
    invJ    = 2/dx; # when using dN
    kappaEl = [Array{Float64}(undef,nnodel,1) for i ∈ 1:Threads.nthreads()]
    # qEl     = Array{Float64}(undef,nnodel,1)
    # q       = fill(0.0,nnod)
    # q[1]    = 1.0

    # Get shape functions
    weight,N,dN  = shapefunctions(1)

    # Allocate FEM element matrices
    Kel     = [fill(0.0,nnodel,nnodel)  for i ∈ 1:Threads.nthreads()]
    # Qel     = Array{Float64}(undef,nnodel,1); fill!(Qel,0.0)

    # Temporal allocations
    # Q_all   = Array{Float64}(undef,nel,nnodel);         # storage for mass matrix coefficients

    # Other allocations
    kappaEl = [Vector{Float64}(undef,2)  for i ∈ 1:Threads.nthreads()]
    KappaEl = [Array{Float64}(undef,1)  for i ∈ 1:Threads.nthreads()]
    ij      = Array{Int64}(undef,nnodel,1)
    
    A       = [dN[ip]' * dN[ip] * weight * invJ for ip = 1:2]
    B       = [N[ip]'  * N[ip]  * weight * detJ for ip = 1:2]
    Mel     = B[1]+B[2]
    ndel    = nnodel*nnodel

    # Assembly loop -> add blocks 
    Threads.@threads for iel = 1 : nel

        kappaEl[Threads.threadid()] = @SVector [κ[e2n[ii,iel]] for ii = 1: 2]  # nnodel = 2

            # # Numerical integration loop
        @views for ip = 1 : nip
            # kappa at i-th integration point
            mul!(KappaEl[Threads.threadid()], N[ip],kappaEl[Threads.threadid()])
            # Conductivity elemental matrix
            @. (Kel[Threads.threadid()]  += A[ip] * KappaEl[Threads.threadid()])
        end
        
        @inbounds for ii ∈ 1:ndel
            i       = @views iMC[ii,iel]
            j       = @views jMC[ii,iel]
            K[i,j] += @views Kel[Threads.threadid()][ii]
            M[i,j] += Mel[ii]
        end
           
        # Reset element matrices
        fill!(Kel[Threads.threadid()], 0.0)

    end

    return K,M

end #### END matrixassembly FUNCTION ###########################################

#### START matrixassembly FUNCTION #############################################
function sparsitymatrix(e2n::Array{Int64,2},nel::Int64,nnodel::Int64)

    # nnodel  = 2

    # Temporal allocations
    iMC     = Array{Int64}(undef,nel,nnodel*nnodel);    # storage for the matrix i indices
    jMC     = Array{Int64}(undef,nel,nnodel*nnodel);    # storage for the matrix j indices

    # Other allocations
    ij      = @MVector ones(Int32,nnodel)
    I1      = (@SVector ones(Int32,nnodel))'
    # I1      = SVector(1,nnodel-1)'
    # ij      = MVector(1,nnodel-1)

    # Assembly loop
    @inbounds for iel = 1 : nel
        
        # Conductivity and Sediment flux at the nodes of i-th element
        @avx for ii = 1: nnodel
            ij[ii]      = e2n[ii,iel]            
        end
      
        rows          = ij * I1
        iMC[iel,:]    = vec(rows)
        jMC[iel,:]    = vec(rows')

    end

    M   = sparse(vec(iMC),vec(jMC),rand(length(jMC)))
    K   = copy(M)

    return K,M,iMC',jMC'

end #### END matrixassembly FUNCTION ###########################################

@inline function accumQ(idx, V, nnod)
    Q   = fill(0.0,nnod)
    @inbounds @simd for ii = 1 : nnod
        Q[idx[ii]] += V[ii]
    end
    return Q
end

#### START dirichlet FUNCTION ##################################################
function dirichlet(KK,bc_rhs,bcdof,bcval)

    @inbounds for i = 1:length(bcdof)
        bc_rhs                 += Array(KK[:,bcdof[i]] * bcval[i])
        KK[bcdof[i],:]         .= 0;
        KK[:,bcdof[i]]         .= 0;
        KK[bcdof[i],bcdof[i]]   = 1;
    end

    return KK,bc_rhs

end #### END dirichlet FUNCTION ################################################

## START accumarray FUNCTION ###################################################
@inline function accumarray(subs,M)
    C   = Vector{Float64}(undef,n)
    @avx for i ∈ 1:length(subs)
        idx     = subs[i]
        C[idx] += M[i]
    end

end #### END accumarray FUNCTION ###############################################

## START norm2 FUNCTION ########################################################
@inline function norm2(A,B)
    norm   = 0.0
    @inbounds for i ∈ 1:length(A)
         norm += (A[i] - B[i])^2
    end
    return sqrt(norm)

end #### END norm2 FUNCTION ####################################################

# end #### END module FUNCTION ###################################################

function mesher(Δx,slope,type::String)
    if type == "linear"
        x,z,dx2,nn,nel,e2n = linearmesh(Δx,slope)
    elseif type == "quadratic"
        x,z,dx2,nn,nel,e2n = quadraticmesh(Δx,slope)
    end
    
    return x,z,dx2,nn,nel,e2n
end

function linearmesh(Δx,slope)
    # -- Mesh
    x       = collect(0:Δx:200e3)
    z       = @. slope * x + 2e3    # 2e3 is the intercept
    dx2     = Δx * Δx
    nn      = length(x)

    # -- Connectivity matrix
    nel = nn - 1
    e2n = Array{Int64}(undef,2,nel)
    @simd for ii = 1:nel
        e2n[1,ii] = ii
        e2n[2,ii] = ii+1
    end
    
    # Output
    return x,z,dx2,nn,nel,e2n

end

function quadraticmesh(Δx,slope)
    # -- Mesh
    x       = collect(0:Δx/2:200e3)
    z       = @. slope * x + 2e3    # 2e3 is the intercept
    dx2     = Δx * Δx
    nn      = length(x)

    # -- Connectivity matrix
    nel     = Int(nn/2 - 1)
    e2n     = Array{Int64}(undef,3,nel)
    idx     = 1:2:nn
    @simd for ii = 1:nel        
        e2n[1,ii] = idx[ii]
        e2n[2,ii] = idx[ii]+1
        e2n[3,ii] = idx[ii]+2
    end
    
    # -- Output
    return x,z,dx2,nn,nel,e2n

end