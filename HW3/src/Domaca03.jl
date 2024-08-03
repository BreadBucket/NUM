module Domaca03
export nihalo, nihalo_path


################################################################


"""
	Calculate pendulum position (angle offset) after some time.
		l:   pendulum length
		t:   time in seconds
		θ0:  initial offset
		Δθ0: initial angular velocity
		n:   simulation granularity (t/n)
"""
function nihalo(l, t, θ0, Δθ0, n)
	if (t < 0)
		throw("t must be non-negative.")
	elseif (l <= 0)
		throw("l must be positive.")
	elseif (n <= 0)
		throw("n must be positive.")
	end
	
	g = 9.80665 / l
	h = t/n
	
	ϕ = θ0
	v = Δθ0
	
	function f(pos, vel)
		return [vel, -g * sin(pos)]
	end
	
	for i in 1:n
		k1 = f(ϕ, v) * h
		k2 = f(ϕ + k1[1]/2, v + k1[2]/2) * h
		k3 = f(ϕ + k2[1]/2, v + k2[2]/2) * h
		k4 = f(ϕ + k3[1],   v + k3[2]) * h
		
		ϕ += (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6
		v += (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6
	end
	
	return ϕ
end


"""
	Calculate pendulum positions and their timestamps durring the simulation
		l:   pendulum length
		t:   time in seconds
		θ0:  initial offset
		Δθ0: initial angular velocity
		n:   simulation granularity (t/n)
"""
function nihalo_path(l, t, θ0, Δθ0, n)
	if (t < 0)
		throw("t must be positive.")
	elseif (l <= 0)
		throw("l must be positive.")
	elseif (n <= 0)
		throw("n must be positive.")
	end
	
	g = 9.80665 / l
	h = t/n
	
	ϕ = θ0
	v = Δθ0
	
	x = zeros(n+1)
	y = zeros(n+1)
	x[1] = 0
	y[1] = ϕ
	
	function f(pos, vel)
		return [vel, -g * sin(pos)]
	end
	
	for i in 1:n
		k1 = f(ϕ, v) * h
		k2 = f(ϕ + k1[1]/2, v + k1[2]/2) * h
		k3 = f(ϕ + k2[1]/2, v + k2[2]/2) * h
		k4 = f(ϕ + k3[1],   v + k3[2]) * h
		
		ϕ += (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6
		v += (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6
		
		x[i+1] = h*i
		y[i+1] = ϕ
	end
	
	return (x, y)
end


################################################################
end # module Domaca03
