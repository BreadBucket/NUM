module Domaca01
import Base: copy
import Base: getindex, setindex!, firstindex, lastindex, size
import Base: *, \, -
using LinearAlgebra


################################################################


struct Tridiagonalna
	l # elements on lower diagonal
	d # elements on main diagonal diagonal
	u # elements on upper diagonal
end

firstindex(T::Tridiagonalna) = 1
firstindex(T::Tridiagonalna, ax) = 1

lastindex(T::Tridiagonalna) = length(T.d)
lastindex(T::Tridiagonalna, ax) = (ax == 1 || ax == 2) ? length(T.d) : 0

size(T::Tridiagonalna) = (length(T.d), length(T.d))
size(T::Tridiagonalna, ax) = (ax == 1) ? length(T.d) : (ax == 2) ? length(T.d) : 0

"""
	e = T[i,j]
	Return element from i-th row and j-th column.
"""
function getindex(T::Tridiagonalna, i::Int, j::Int)
	if j == i
		return T.d[j]
	elseif j+1 == i
		return T.l[j]
	elseif j-1 == i
		return T.u[i]
	else
		return zero(T.d[1])
	end
end

"""
	T[i,j] = e
	Set element at i-th row and j-th column.
"""
function setindex!(T::Tridiagonalna, e, i::Int, j::Int)
	if j == i
		T.d[j] = e
	elseif j+1 == i
		T.l[j] = e
	elseif j-1 == i
		T.u[i] = e
	end
end

"""
	Multiply matrix T with vector x.
"""
function *(T::Tridiagonalna, x::Vector)
	y = zero(x)
	n = length(T.d)
	
	y[1] = T[1,1] * x[1] + T[1,2] * x[2]  # prva dva
	for i=2:n-1
		y[i] = T[i,i-1] * x[i-1] + T[i,i] * x[i] + T[i,i+1] * x[i+1]  # po tri elemente
	end
	y[n] = T[n, n-1] * x[n-1] + T[n,n] * x[n]  # zadnje dva
	
	return y
end

"""
	Creare copy of T
"""
function copy(T::Tridiagonalna)
	return Tridiagonalna(copy(T.l), copy(T.d), copy(T.u))
end


################################################################


struct ZgornjaTridiagonalna
	d	# diagonal
	u	# top diagonal
end

firstindex(T::ZgornjaTridiagonalna) = 1
firstindex(T::ZgornjaTridiagonalna, ax) = 1

lastindex(T::ZgornjaTridiagonalna) = length(T.d)
lastindex(T::ZgornjaTridiagonalna, ax) = (ax == 1 || ax == 2) ? length(T.d) : 0

size(T::ZgornjaTridiagonalna) = (length(T.d), length(T.d))
size(T::ZgornjaTridiagonalna, ax) = (ax == 1) ? length(T.d) : (ax == 2) ? length(T.d) : 0

"""
	e = T[i,j]
	Return element from i-th row and j-th column.
"""
function getindex(T::ZgornjaTridiagonalna, i::Int, j::Int)
	if j == i
		return T.d[i]
	elseif j-1 == i
		return T.u[i]
	else
		return zero(T.d[1])
	end
end

"""
	T[i,j] = e
	Set element at i-th row and j-th column.
"""
function setindex!(T::ZgornjaTridiagonalna, e, i::Int, j::Int)
	if j == i
		T.d[i] = e
	elseif j-1 == i
		T.u[i] = e
	end
end


"""
	Solve L*x = v
"""
function \(U::ZgornjaTridiagonalna, v)
	n = length(v)
	v = copy(v)
	
	if (U[n,n] != 0)
		v[n] /= U[n,n]
	else
		v[n] = 0
	end
	
	for i=n-1:-1:1
		v[i] -= U[i, i+1]*v[i+1]
		if (U[i,i] != 0)
			v[i] /= U[i,i]
		elseif (v[i] != 0)
			throw("Unsolveable.")
		end
	end
	
	return v
end


################################################################


struct SpodnjaTridiagonalna
	l	# lower diagonal
	d	# diagonal
end

firstindex(T::SpodnjaTridiagonalna) = 1
firstindex(T::SpodnjaTridiagonalna, ax) = 1

lastindex(T::SpodnjaTridiagonalna) = length(T.d)
lastindex(T::SpodnjaTridiagonalna, ax) = (ax == 1 || ax == 2) ? length(T.d) : 0

size(T::SpodnjaTridiagonalna) = (length(T.d), length(T.d))
size(T::SpodnjaTridiagonalna, ax) = (ax == 1) ? length(T.d) : (ax == 2) ? length(T.d) : 0

"""
	e = T[i,j]
	Return element from i-th row and j-th column.
"""
function getindex(T::SpodnjaTridiagonalna, i::Int, j::Int)
	if j == i
		return T.d[i]
	elseif j+1 == i
		return T.l[j]
	else
		return zero(T.d[1])
	end
end

"""
	T[i,j] = e
	Set element at i-th row and j-th column.
"""
function setindex!(T::SpodnjaTridiagonalna, e, i::Int, j::Int)
	if j == i
		T.d[i] = e
	elseif j+1 == i
		T.l[j] = e
	end
end


"""
	Solve L*x = v
"""
function \(L::SpodnjaTridiagonalna, v)
	n = length(v)
	v = copy(v)
	
	if (L[1,1] != 0)
		v[1] /= L[1,1]
	else
		v[1] = 0
	end
	
	for i=2:n
		v[i] -= L[i, i-1]*v[i-1]
		if (L[i,i] != 0)
			v[i] /= L[i,i]
		elseif (v[i] != 0)
			throw("Unsolveable.")
		end
	end
	
	return v
end


###############################################################


"""
	Multiply L*U
"""
function *(L::SpodnjaTridiagonalna, U::ZgornjaTridiagonalna)
	n = size(L)[1]
	A = Tridiagonalna(zeros(n-1), zeros(n), zeros(n-1))
	
	for y=1:n
		for x=max(1,y-1):min(n,y+1)
			for i=max(1,y-1):y
				A[y,x] += L[y,i] * U[i,x]
			end
		end
	end
	
	return A
end


###############################################################


"""
	Calculate matrices L and U such that: A = L*U
"""
function lu(A::Tridiagonalna)
	n = size(A)[1]
	L = SpodnjaTridiagonalna(zeros(n-1), ones(n))
	U = ZgornjaTridiagonalna(copy(A.d), copy(A.u))
	
	for i=2:n
		j = i-1
		L[i,j] = l = A[i,j] / U[i-1,j]
		U[i,j+1] -= U[i-1,j+1] * l
	end
	
	return L, U
end


"""
	Calculate Householder transformation.
"""
function tridiag(A)
	n = size(A, 1)
	if (n < 2) || (size(A, 2) < 2)
		throw("Velikos matrike mora biti vsaj 3x3.")
	end
	
	sign(n) = (n < 0) ? -1 : 1
	
	A = copy(A)
	tQ, Q = zeros(n,n), zeros(n,n) + I	# Temporary Q and final Q
	H = Tridiagonalna(zeros(n-1), zeros(n), zeros(n-1))
	
	for k=1:n-2
		# Calculate alpha
		a = 0
		for i=k+1:n
			a += A[i,k] * A[i,k]
		end
		a = -sign(A[k+1,k]) * sqrt(a)
		
		
		# Skip sqrt and simplify multiplication
		# r = √(a * (a - A[k+1,k]) / 2)
		rr = (a - A[k+1,k]) * a		# (r*2)^2 / 2
		
		# Define `v`; avoid instantiation
		# v[1...k]   = 0
		# v[k+1]     = (A[k+1,k] - a)/(r*2)
		# v[k+2...n] = A[k+1,k]/(r*2)
		vk1 = A[k+1,k] - a	# v[k+1] but simplified in the whole context
		
		
		# Avoid sqrt() and array `v` instantiation with manual multiplication.
		# Q = I - 2⋅v⋅v^T
		
		# Prepare identity part of Q by cleaning up previous Q
		# [ 1       ]
		# [   1 0 0 ]
		# [   0 ? ? ]
		# [   0 ? ? ]
		tQ[k,k] = 1
		for i=k+1:n
			tQ[k,i] = 0
			tQ[i,k] = 0
		end
		
		# Fill remaining part of Q with `I - 2⋅v⋅v^T`
		for y=k+1:n
		for	x=k+1:n
			yv = (y == k+1) ? vk1 : A[y,k]
			xv = (x == k+1) ? vk1 : A[x,k]
			vvT2_yx = yv*xv/rr	# (2⋅v⋅v^T)[y,x]; 2 is included in `rr`
			
			# I - (2⋅v⋅v^T)
			if (x == y)
				tQ[y,x] = 1 - vvT2_yx
			else
				tQ[y,x] = -vvT2_yx
			end
			
		end
		end
		
		A = tQ*A*tQ	# Q == Q^T
		Q[k:n, k+1:n] = Q[k:n, k+1:n] * tQ[k+1:n, k+1:n]
		
		# Copy to H
		H[k,k] = A[k,k]
		H[k+1,k] = A[k+1,k]
		H[k,k+1] = A[k,k+1]
	end
	
	# Copy to H (bottom right 2x2)
	H[n-1,n-1] = A[n-1,n-1]
	H[n-1,n] = A[n-1,n]
	H[n,n-1] = A[n,n-1]
	H[n,n] = A[n,n]
	
	return H, Q
end


###############################################################


"""
    Calculate eigen-value and eigen-vector
     of a matrix `A` from an initial approximation `l`.
"""
function inv_lastni(A, l, maxit::Int=100, tol=1e-10)
	H, Q = tridiag(A)
	
	# H - λI
	Hλ = copy(H)
	for i=1:size(Hλ, 1)
		Hλ[i,i] -= l
	end
	
	L, U = lu(Hλ)
	x = rand(size(A, 2))
	λ = l
	
	# Series of approximations
	for i=1:maxit
		y = U \ (L \ x)
		
		# Select abs-maximal λ
		max = argmax(abs.(y))
		λ = (x[max] / y[max]) + l
		x = y / y[max]
		
		if norm(H*x - λ*x, Inf) < tol
			return λ, Q*x
		end
		
	end
	
	throw("Did not converge.")
end


###############################################################


end # module Domaca01
