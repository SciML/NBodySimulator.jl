@@ -0,0 +1,104 @@
using Plots
t = t1:τ:result.solution.t[end-1]
cc = get_position.(result, t)
pl = plot()
for i = 1:length(result.simulation.system.bodies)

    y = zeros(length(t))
    for j=1:length(t)
        positions = @view cc[j]
        if result.simulation.boundary_conditions isa PeriodicBoundaryConditions
            L = result.simulation.boundary_conditions[2]
            map!(x ->  x -= L * floor(x / L), positions, positions)
        end
        y[j]=cc[j][:,i][3]
    end
    plot!(pl, t, y, ylim = [0, L], xlim = [t1:t[end]])
end


#-------------Velocity

using Plots
t = t1:τ:result.solution.t[end-1]
cc = get_velocity.(result, t)
pl = plot()
for i = 1:length(result.simulation.system.bodies)

    y = zeros(length(t))
    for j=1:length(t)
        y[j]=cc[j][:,i][3]
    end
    plot!(pl, t, y, ylim = [0, L], xlim = [t1:t[end]])
end
display(pl)


#---------Temperature
using Plots
t = t1:τ:result.solution.t[end-1]
temper = temperature.(result, t)
plot(t, temper, xlabel="t, s", ylabel = "T, °K", label="Temperature, °K" )

#--------Energy
using Plots
t = t1:τ:result.solution.t[end-2]
e_kin = kinetic_energy.(result, t)/(kb*T)
e_pot = potential_energy.(result, t)/(kb*T)
#e_tot = total_energy.(result, t)/(kb*T)
e_tot = e_kin+e_pot
plot(t, [signif.(e_tot,4), signif.(e_kin,4), signif.(e_pot,4)], label=["Total energy","Kinetic energy","Potential energy"], legend = :right, ylabel="E/kT", xlabel="t, s", formatter=:scientific)
plot(d["t"], [signif.(d["e_tot"],4), signif.(d["e_kin"],4), signif.(d["e_pot"],4)], label=["Total energy","Kinetic energy","Potential energy"], legend = :right, ylabel="E/kT", xlabel="t, s", formatter=:scientific)


_e_tot=signif.(e_tot,4)
_e_pot=signif.(e_pot,4)
_e_kin=signif.(e_kin,4)
plot(t, [_e_tot, _e_kin, _e_pot], label=["Total energy","Kinetic energy","Potential energy"], legend = :right, ylabel="Energy, J", xlabel="t, s", formatter=:scientific)

#----------------------------------------
using DifferentialEquations
function f!(du,u,p,t)
    du[3] = -u[4]
    du[4] = -u[3]
    du[1] = dot([du[3]-du[1], du[4]-du[2]], [u[3]-u[1], u[4]-u[2]])
    du[2] = sqrt((u[3]-u[1])^2+(u[4]-u[2])^2)-1

end

# x1, y1, x2, y2
u0 = [0.0, 0.5, 0.0, -0.5]
M = [0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1]
tspan = (0.0, 10.0)

prob = ODEProblem(f!,u0,tspan, mass_matrix=M)
sol = solve(prob, Rosenbrock23(), dt = 0.01);

using Plots
u = sol(2)
scatter([u[1], u[3]], [u[2], u[4]], ylim = [0, 1.0], xlim = [-1.0, 0.0])

#-----------rdf--------------
using Plots
import GR
(rs, grf) = rdf(result)
plot(rs/σ, grf, xlim=[0, 0.4999L/σ], label=["Radial distribution function"],ylabel="g(r)", xlabel="r/σ")
plot(d["rs"]/σ, d["grf"], xlim=[0, 0.4999L/σ], label=["Radial distribution function"],ylabel="g(r)", xlabel="r/σ")

using Plots
import GR
(rs, grf) = rdf(result)
plot(rs/1e-10, grf, xlim=[0, 0.4999L/1e-10], label=["Radial distribution function"],ylabel="g(r)", xlabel="r, Å")
plot(d["rs"]/1e-10, d["grf"], xlim=[0, 0.4999L/1e-10], label=["Radial distribution function"],ylabel="g(r)", xlabel="r, Å")
plot(d["rs"], d["grf"], label=["Radial distribution function"],ylabel="g(r)", xlabel="r, Å")
plot(rs, grf, label=["Radial distribution function"],ylabel="g(r)", xlabel="r, Å")

#-------------msd-------------------
using Plots
import GR
plot(ts, dr2/1e-20, label=["Mean square distance"],ylabel="MSD, Å^2", xlabel="t, s")
plot(d["ts"], d["dr2"]/1e-20, label=["Mean square distance"],ylabel="MSD, Å^2", xlabel="t, s")
plot(ts, 100*dr2, label=["Mean square distance"],ylabel="MSD, Å^2", xlabel="t, s")
#----------aus file-------------
e_kin = d["e_kin"];
e_tot = d["e_tot"];
e_pot = d["e_pot"];



#---------------------------------

T0 = d["T0"]
temper = d["temper"]
thermostat = d["thermostat"]
t = d["t"]


#---------------------------------
Ok, I just updated AbstractPlotting! If you do `Pkg.add("Makie"); Pkg.checkout.(("Makie", "AbstractPlotting")))` the following should work
https://gist.github.com/SimonDanisch/8c59e79e6eaf8267796334bc2e6daa25#file-molecule_recipe-jl


function calculate!(c,a,b)
    @. c = a*b
end

function test()
    n = 100
    a = rand(n,n)
    b = rand(n,n)
    c = zero(a)
    m=70
    calculate!(c[:,1:m], a[:,1:m], b[:,1:m])
end