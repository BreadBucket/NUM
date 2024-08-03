using Domaca03
using Test


@testset "nihalo 0" begin
	t = 0
	l = 1.2
	pos_0 = pi/2
	
	Y = nihalo(l, t, pos_0, 0, 10);
	@test isapprox(Y, pos_0, rtol=1e-6)
end


@testset "nihalo & nihalo_path" begin
	t = 10
	l = 1.2
	pos_0 = pi/2
	
	x, y = nihalo_path(l, t, pos_0, 0, t*1000);
	Y = nihalo(l, t, pos_0, 0, t*1000);
	
	@test x[1] == 0
	@test y[1] == pos_0
	@test isapprox(x[end], t, rtol=1e-6)
	@test isapprox(y[end], Y, rtol=1e-6)
end
