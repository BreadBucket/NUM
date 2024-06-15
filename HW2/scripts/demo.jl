using Domaca02
using Printf


println("Normal distribution:")
for x in (-0.99:0.33:1)
	@printf("  Φ(%5.2f) = %.6f\n", x, normal(x))
end
println()

# Distance between two faces of A and B is 1
A = ((0,0,0), (1,1,1))
B = ((2,0,0), (3,1,1))
f = cubeForce(A, B, 10)

println("Force between cubes A and B is:\n  $(f).")
