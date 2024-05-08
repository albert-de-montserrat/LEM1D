function terrace_vertical_erosion!(
    terrace::TerraceProfile,
    params::TerraceParams,
    dt::Float64,
    h_sea,
    i1::Int64,
)

    (;h_wb, βz, P0, sol_n, sol_old) = params
    sol  = terrace.z
    i2   = terrace.nnod
    cte1 = βz * P0 * dt
    cte2 = 4.0 / h_wb
    copyprofile!(sol_old, sol)
    residual = Inf

    while residual > 1e-8
        copyprofile!(sol_n, sol)
        @tturbo for ii in i1:i2
            sol[ii] -= f(cte1, cte2, h_sea, sol[ii], sol_old[ii])
        end
        residual = norm2(sol, sol_n)
    end

end

function terrace_vertical_erosion!(
    βz,
    P,
    h_wb,
    dt,
    sol,
    h_sea,
    i1::Int64,
    i2::Int64,
    sol_n,
    sol_old,
)
    cte1 = βz * P * dt
    cte2 = 4.0 / h_wb
    copyprofile!(sol_old, sol)
    residual = Inf

    while residual > 1e-8
        copyprofile!(sol_n, sol)
        @tturbo for ii in i1:i2
            sol[ii] -= f(cte1, cte2, h_sea, sol[ii], sol_old[ii])
        end
        residual = norm2(sol, sol_n)
    end

end

@inline _f(cte1, cte2, h_sea, x, x_old) = muladd(cte1, exp(-cte2 * (h_sea - x)), x) - x_old

@inline function f(cte1, cte2, h_sea, x, x_old)
    x, dx = f_df(x->_f(cte1, cte2, h_sea, x, x_old), x)
    return x * inv(dx)
end

function terrace_vertical_explicit!(
    βz,
    P,
    h_wb,
    dt,
    sol,
    h_sea,
    i1::Int64,
    i2::Int64,
)
    cte1 = βz * P * dt
    cte2 = 4.0 / h_wb

    for ii in i1:i2
        sol[ii] -= cte1 * exp(-cte2 * (h_sea - sol[ii]))
    end
end

function terrace_vertical_explicit(
    βz,
    P,
    h_wb,
    dt,
    sol,
    h_sea,
)
    cte1 = βz * P * dt
    cte2 = 4.0 / h_wb

    return @. sol - cte1 * exp(-cte2 * (h_sea - sol))
end

function terrace_retreat!(
    profile,
    params::TerraceParams,
    dt::Float64,
    h_sea::Float64,
    i1::Int64,
)

    @inline sea_floor(h) = h_sea - h

    (;h_wb, βx, P0, Poff) = params
    
    dx = (profile.x[2] - profile.x[1])
    _h_wb = inv(h_wb)
    int = 0.0
    @inbounds @fastmath for i in i1:profile.nnod-1
        h = (h_sea - profile.z[i])
        h > 200.0 && break # lets define the shelf at depths of > 200m
        int += exp( -2.0 * (sea_floor(profile.z[i]) + sea_floor(profile.z[i+1])) * _h_wb) 
        # it should be multiplied by 4, but we use two becase we take the mean
        # of two consecutive sea floor elevation elements
    end
    Δx = -dt * βx * (Poff - P0 * int * dx) # amount of retreat
    
    # erode cliff
    cliff_x  = profile.x[i1] + Δx # new x-position
    cliff_z0 = profile.z[i1]
    @inbounds for i in i1:-1:1
        profile.x[i] ≤ cliff_x && break
        profile.z[i] = cliff_z0
    end

end

@inline function copyprofile!(dest, profile)
    @tturbo for i in eachindex(profile)
        dest[i] = profile[i]
    end
end

function find_shore_id(h_sea::T,sol::Vector{T},n::Int64) where T
    for ii = 1:n-1
        @inbounds if sol[ii] ≤ h_sea
            return ii
            break
        end
    end
end
