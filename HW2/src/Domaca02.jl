module Domaca02
export romberg, trapezoidSegmented, simpsonSegmented
export force, cubeForce, normal


################################################################


"""
	Calculate integral using trapezoid method.
		f: integrated function
		a: lower bound of integration
		b: upper bound of integration
"""
function trapezoid(f, a, b)
	return ((b-a)*(f(a)+f(b)))/2
end


"""
	Calculate integral using segmented trapezoid method.
		f: integrated function
		a: lower bound of integration
		b: upper bound of integration
		h: step size
"""
function trapezoidSegmented(f, a, b, h)
	sum = 0
	
	if b < a
		_a = a; a = b; b = _a
	end
	
	p = a
	while p < b
		p2 = p+h
		if p2 > b
			break
		end
		sum += f(p) + f(p2)
		p = p2
	end
	
	sum = sum*h/2
	
	# Add area of final strip (`(b-a)%h ≠ 0`)
	if p < b
		sum += trapezoid(f, p, b)
	end
	
	return sum
end


"""
	Calculate integral using rombergs method.
		f:     integrated function
		a:     lower bound of integration
		b:     upper bound of integration
		tol:   error tolerance
		maxit: maximum number of divisions
"""
function romberg(f, a, b; tol=1e-10, maxit::Int=50)
	h = (b-a)
	R = [ trapezoidSegmented(f, a, b, h) ]
	
	for k in 2:maxit
		h /= 2
		r = R[1]	# Virtual R[row-1,col-1]
		
		push!(R, R[end])	# Assumption: capacity not tight
		R[1] = trapezoidSegmented(f, a, b, h)
		
		# Calc row
		m4 = 1
		for i in 2:k
			m4 *= 4
			_r = R[i]
			R[i] = ((m4 * R[i-1]) - r) / (m4 - 1)
			r = _r
		end
		
		# Error tolerance reached
		if abs(R[end] - r) <= tol
			return R[end]
		end
		
	end
	
	throw("Maximum allowed iterations reached without convergence.")
end


"""
	Calculate normal distribution up to `x`: Φ(x) = P(X ≤ x)
	Integral is calculated using the adaptive simpsons method.
"""
function normal(x; tol=1e-10)
	if x == 0
		return 0.5
	end
		
	m = 1/sqrt(2*pi)
	f(t) = exp(-t*t/2)
	
	if x > 0
		return 0.5 + m*romberg(f, 0, x, tol=tol)
	else
		return 0.5 - m*romberg(f, x, 0, tol=tol)
	end
		
end


################################################################


"""
	Calculate integral using segmented simpsons method.
		f: integrated function
		a: lower bound of integration
		b: upper bound of integration
		n: number of segments
	Number of invocations of `f` is (2*n + 1).
"""
function simpsonSegmented(f, a, b, n::Int)
	l = b - a
	h = (l/n)/2.0
	hh = h+h
	
	p = a
	sum = f(p)
	
	for _ in 1:n-1
		sum = sum .+ (f(p+h).*4) .+ (f(p+hh).*2)
		p += hh
	end
	
	sum = sum .+ (f(p+h).*4) .+ f(p+hh)
	return sum .* (h/3)
end


################################################################


"""
	Force between two different points with constants set to 1.
"""
function force(x1, y1, z1, x2, y2, z2)
	x = x2 - x1
	y = y2 - y1
	z = z2 - z1
	l = x*x + y*y + z*z
	
	if l == 0
		throw("Null vector: ($x1,$y1,$z1) -- ($x2,$y2,$z2)")
	end
	
	return (x/l, z/l, y/l)
end

"""
	Force between two cubes with constants set to 1.
	The force is calculated using the segmented simpsons method.
		cubeA: tuple of lower and uper bounds of the first cube
		cubeB: tuple of lower and uper bounds of the second cube
		n:     number of divisions in the simpsons method
		
"""
function cubeForce(cubeA::Tuple{Tuple,Tuple}, cubeB::Tuple{Tuple,Tuple}, n::Int)
	Ax1, Ay1, Az1 = cubeA[1]
	Ax2, Ay2, Az2 = cubeA[2]
	Bx1, By1, Bz1 = cubeB[1]
	Bx2, By2, Bz2 = cubeB[2]
	
	return simpsonSegmented(
		x1 -> simpsonSegmented(
		y1 -> simpsonSegmented(
		z1 -> simpsonSegmented(
		x2 -> simpsonSegmented(
		y2 -> simpsonSegmented(
		z2 -> (
			force(x1, y1, z1, x2, y2, z2)
		), Bz1, Bz2, n
		), By1, By2, n
		), Bx1, Bx2, n
		), Az1, Az2, n
		), Ay1, Ay2, n
		), Ax1, Ax2, n
	)
end


################################################################
end # module Domaca02
