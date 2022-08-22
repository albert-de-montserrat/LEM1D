@inline function norm2(A::Vector{T}, B::Vector{T}) where T
    norm = zero(T)
    @tturbo for i in eachindex(A)
         norm += (A[i] - B[i])^2
    end
    √(norm)
end 

@inline function norm2(A::Vector{T}, B::Vector{T}, indices::NamedTuple) where T
    norm = zero(T)
    i1, i2 = indices.i1, indices.i2
    @tturbo for i in i1:i2
         norm += (A[i] - B[i])^2
    end
    √(norm)
end 

function variable_recurrence_time_fixed(
    time_intervals, variable_rec_time, background_rec_time, t
)
    for i in eachindex(time_intervals)
        if time_intervals[i][1] ≤ t ≤ time_intervals[i][2]
            return variable_rec_time[i]
        end
    end
    return background_rec_time
end

function poisson_time_serie(ti, tf, recurrence_t)
    time_lentgh = abs(tf - ti)
    eq_times = [-1]
    eq_distribution = Poisson(1 / recurrence_t)
    eq_times = rand(eq_distribution, Int(time_lentgh))
    eq_times = findall(i -> i != 0, eq_times) .+ Int(ti)
    return eq_times
end

function update_reoccupation!(reoccupation_time, reoccupation_id, t, h_sea, z)
    # Threads.@threads 
    @inbounds for i in eachindex(z)
        if z[i] > h_sea
            reoccupation_id[i] = :aerial
        else
            if reoccupation_id[i] === :reoccupied
                reoccupation_time[i, 2] = t

            elseif reoccupation_id[i] === :aerial
                reoccupation_id[i] = :reoccupied
                @views reoccupation_time[i, 1:2] .= t
                reoccupation_time[i, 3] += 1
            end
        end
    end
end
