using Domaca03
using Plots
using Printf


t = 10
l = 1.2
pos_0 = pi/2

print("Starting position: ")
print(pos_0 * 180/pi)
println("°")

x, y = nihalo_path(l, t, pos_0, 0, t*1000);
# Y = nihalo(l, t, pos_0, 0, t*1000);	# Y == y[end]

@printf "Final position after %g seconds: " t;
print(y[end] * 180/pi)
println("°")

plot(x, y * 180/pi, label="pendulum", legend=:topright)
plot!(x, sin.(x) * 90, label="sin(x)*90")
plot!(x, sin.(x * 2.5) * 90, label="sin(x*2.5)*90")
