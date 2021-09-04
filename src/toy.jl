
function river_length(River)
    x, z, dx = River.x, River.z, River.Δx^2
    L = zero(x)
    @tturbo for i in 2:length(x)
        L[i] = √(dx  + (z[i]-z[i-1])^2) 
    end
    cumsum(L)
end

function lazycumsum(A::Vector{T}, idx) where T
    s = zero(T)
    @turbo for i in 1:idx
        s += A[i]
    end
    s
end

