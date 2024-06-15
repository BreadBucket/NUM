using Domaca02
using Test

@testset "Normal" begin
	@test isapprox(normal(-2), 0.02275, rtol=1e-3)
	@test isapprox(normal(-1), 0.15866, rtol=1e-3)
	@test normal(0) ≈ 0.5
	@test isapprox(normal(1),  0.84134, rtol=1e-3)
	@test isapprox(normal(2),  0.97725, rtol=1e-3)
end
