using Domaca01: Tridiagonalna, ZgornjaTridiagonalna, SpodnjaTridiagonalna
using Domaca01: tridiag, lu, inv_lastni


A = [4.0 1 -2 2 ; 1 2 0 1 ; -2 0 3 -2 ; 2 1 -2 -1]
H, Q = tridiag(A)

# Get smallest eigenvalue
λ, v = inv_lastni(A, 0)

# Check
test = (A * v ./ v .- λ)
println(isapprox(test, [0,0,0,0], atol=1e-6))
test