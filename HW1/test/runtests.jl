using Domaca01: Tridiagonalna, ZgornjaTridiagonalna, SpodnjaTridiagonalna
using Domaca01: lu, inv_lastni
using LinearAlgebra
using Test


################################################################


@testset "Triagonalna matrika" begin
	T = Tridiagonalna([1,2], [3,4,5], [6,7])
	M = [3 6 0; 1 4 7 ; 0 2 5]
	
	@test T[1,1] == M[1,1]
	@test T[1,2] == M[1,2]
	@test T[1,3] == M[1,3]
	@test T[2,1] == M[2,1]
	@test T[2,2] == M[2,2]
	@test T[2,3] == M[2,3]
	@test T[3,1] == M[3,1]
	@test T[3,2] == M[3,2]
	@test T[3,3] == M[3,3]
	
	T[2,2] = M[2,2] = 10
	T[2,1] = M[2,1] = 12
	T[1,2] = M[1,2] = 99
	
	@test T[1,1] == M[1,1]
	@test T[1,2] == M[1,2]
	@test T[1,3] == M[1,3]
	@test T[2,1] == M[2,1]
	@test T[2,2] == M[2,2]
	@test T[2,3] == M[2,3]
	@test T[3,1] == M[3,1]
	@test T[3,2] == M[3,2]
	@test T[3,3] == M[3,3]
end


@testset "ZgornjaTridiagonalna matrika" begin
	T = ZgornjaTridiagonalna([3,4,5], [6,7])
	M = [3 6 0; 0 4 7 ; 0 0 5]
	
	@test T[1,1] == M[1,1]
	@test T[1,2] == M[1,2]
	@test T[1,3] == M[1,3]
	@test T[2,1] == M[2,1]
	@test T[2,2] == M[2,2]
	@test T[2,3] == M[2,3]
	@test T[3,1] == M[3,1]
	@test T[3,2] == M[3,2]
	@test T[3,3] == M[3,3]
	
	T[2,2] = M[2,2] = 10
	T[1,2] = M[1,2] = 99
	
	@test T[1,1] == M[1,1]
	@test T[1,2] == M[1,2]
	@test T[1,3] == M[1,3]
	@test T[2,1] == M[2,1]
	@test T[2,2] == M[2,2]
	@test T[2,3] == M[2,3]
	@test T[3,1] == M[3,1]
	@test T[3,2] == M[3,2]
	@test T[3,3] == M[3,3]
end


@testset "SpodnjaTridiagonalna matrika" begin
	T = SpodnjaTridiagonalna([1,2], [3,4,5])
	M = [3 0 0; 1 4 0 ; 0 2 5]
	
	@test T[1,1] == M[1,1]
	@test T[1,2] == M[1,2]
	@test T[1,3] == M[1,3]
	@test T[2,1] == M[2,1]
	@test T[2,2] == M[2,2]
	@test T[2,3] == M[2,3]
	@test T[3,1] == M[3,1]
	@test T[3,2] == M[3,2]
	@test T[3,3] == M[3,3]
	
	T[2,2] = M[2,2] = 10
	T[2,1] = M[2,1] = 12
	
	@test T[1,1] == M[1,1]
	@test T[1,2] == M[1,2]
	@test T[1,3] == M[1,3]
	@test T[2,1] == M[2,1]
	@test T[2,2] == M[2,2]
	@test T[2,3] == M[2,3]
	@test T[3,1] == M[3,1]
	@test T[3,2] == M[3,2]
	@test T[3,3] == M[3,3]
end


################################################################


function eq(A, B)
	if (size(A) != size(B))
		return false
	end
	
	for y=firstindex(A,1):lastindex(A,1)
	for x=firstindex(A,2):lastindex(A,2)
		if (A[y,x] != B[y,x])
			return false
		end
	end
	end
	
	return true
end


################################################################


@testset "L*U" begin
	ML = [1 0 0 0 ; 2 3 0 0 ; 0 4 5 0 ; 0 0 6 7]
	MU = [1 2 0 0 ; 0 3 4 0 ; 0 0 5 6 ; 0 0 0 7]
	M  = ML*MU
	
	L = SpodnjaTridiagonalna([2,4,6], [1,3,5,7])
	U = ZgornjaTridiagonalna([1,3,5,7], [2,4,6])
	A = L*U
	
	@test eq(M, A)
end


@testset "LU razcep" begin
	M = [2 3 0 ; -6 -11 13 ; 0 0 -8]
	A = Tridiagonalna([-6, 0], [2, -11, -8], [3, 13])
	L, U = lu(A)
	LU = L*U
	
	@test eq(M, A)
	@test eq(LU, A)
end


################################################################


@testset "inv_lastni" begin
	A = [4.0 1 -2 2 ; 1 2 0 1 ; -2 0 3 -2 ; 2 1 -2 -1]
	λ0,v0 = inv_lastni(A, -10)
	λ1,v1 = inv_lastni(A, 0)
	λ2,v2 = inv_lastni(A, 3)
	λ3,v3 = inv_lastni(A, 10)
	
	@test λ0 ≈ -2.197516977439423
	@test A*v0 ≈ λ0*v0
	@test λ1 ≈  1.084364463773226
	@test A*v1 ≈ λ1*v1
	@test λ2 ≈  2.2685314064312427
	@test A*v2 ≈ λ2*v2
	@test λ3 ≈  6.844621107234968
	@test A*v3 ≈ λ3*v3
end


################################################################

