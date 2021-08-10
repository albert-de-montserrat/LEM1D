function partialsort_cutoff!(x,cutoff)
    iswap = 1
    for i in eachindex(x)
      if x[i] <= cutoff
        if iswap != i
          x[iswap], x[i] = x[i], x[iswap]
        end
        iswap = iswap + 1
      end
    end
    return iswap - 1
end

x = rand(20)
y=deepcopy(x)
idx = partialsort_cutoff!(x,0.5)

@benchmark partialsort_cutoff!($y, $0.5) setup=(y=rand(20)) evals=1
