using BenchmarkTools, StructuredLight

suite = BenchmarkGroup(["hg"])

rs = LinRange(-5, 5, 256)
ψ = hg(rs, rs, m=2, n=5)

for m ∈ 1:1, n ∈ 1:1
    suite["hg"]["m=$m; n=$n"] = @benchmarkable hg($rs, $rs; m=$m, n=$n)
end

tune!(suite)
result = run(suite)

result["hg"]["m=1; n=1"]