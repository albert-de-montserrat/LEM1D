# function terrace_vertical_erosion!(
#     βz::Float64,
#     P::Float64,
#     h_wb::Float64,
#     dt::Float64,
#     sol::Vector{Float64},
#     h_sea::Float64,
#     i1::Int64,
#     i2::Int64,
#     sol_n,
#     sol_old,
# )
#     cte1 = βz * P * dt
#     cte2 = 4.0 / h_wb
#     cte3 = cte1 * cte2
#     copyprofile!(sol_old, sol)
#     residual = Inf
#     while residual > 1e-8
#         copyprofile!(sol_n, sol)
#         @turbo for ii in i1:i2
#             cte4 = exp(-cte2 * (h_sea - sol[ii]))
#             sol[ii] -= (sol[ii] + cte1 * cte4 - sol_old[ii]) / (1 + cte3 * cte4)
#         end
#         residual = norm2(sol, sol_n)
#     end
# end

function terrace_vertical_erosion!(
    βz::Float64,
    P::Float64,
    h_wb::Float64,
    dt::Float64,
    sol::Vector{Float64},
    h_sea::Float64,
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
    βz::Float64,
    P::Float64,
    h_wb::Float64,
    dt::Float64,
    sol::Vector{Float64},
    h_sea::Float64,
    i1::Int64,
    i2::Int64,
)
    cte1 = βz * P * dt
    cte2 = 4.0 / h_wb

    @turbo for ii in i1:i2
        sol[ii] -= cte1 * exp(-cte2 * (h_sea - sol[ii]))
    end
end

function terrace_retreat!(
    profile,
    βx::Float64,
    Poff::Float64,
    P0::Float64,
    h_wb::Float64,
    dt::Float64,
    h_sea::Float64,
    i1::Int64,
)
    @inline sea_floor(h) = h_sea - h

    dx = (profile.x[2] - profile.x[1])
    _h_wb = inv(h_wb)
    int = 0.0
    @inbounds for i in i1:profile.nnod-1
        h = (h_sea - profile.z[i])
        h > 200.0 && break # lets define the shelf at depths of > 200m
        int += exp( -2.0 * (sea_floor(profile.z[i]) + sea_floor(profile.z[i+1])) * _h_wb) 
        # it should be multiplied by 4, but we use two becase we take the mean
        # of two consecutive sea floor elevation elements
    end
    Δx = -dt * βx * (Poff - P0 * int * dx) # amount of retreat
    
    # erode cliff
    # cliff_x0 = profile.x[i1] # original x-location
    cliff_x  = profile.x[i1] + Δx # new x-position
    cliff_z0 = profile.z[i1]
    @inbounds for i in i1:-1:1
        profile.x[i] ≤ cliff_x && break
        profile.z[i] = cliff_z0
    end

    # println("Retreat of $(abs(Δx))m. Original x = $(profile.x[i1]), new x = $(cliff_x) ")
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