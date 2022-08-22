function terracenewton!(
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
    cte3 = cte1 * cte2
    copyprofile!(sol_old, sol)
    residual = Inf
    while residual > 1e-6
        copyprofile!(sol_n, sol)
        @turbo for ii in i1:i2
            cte4 = exp(-cte2 * (h_sea - sol[ii]))
            sol[ii] -= (sol[ii] + cte1 * cte4 - sol_old[ii]) / (1.0 + cte3 * cte4)
        end
        residual = norm2(sol, sol_n)
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