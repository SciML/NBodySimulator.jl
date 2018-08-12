abstract type Thermostat
end

struct NullThermostat <: Thermostat
end


struct AndersenThermostat{tType <: Real, νType <: Real} <: Thermostat
    T::tType
    ν::νType
end

struct BerendsenThermostat{τType <: Real,tType <: Real} <: Thermostat
    T::tType
    τ::τType
    γ::τType
end

function BerendsenThermostat(T::Real, τ::Real)
    BerendsenThermostat(T, τ, 0.5 / τ)
end

function berendsen_acceleration!(dv, v, ms, kb, N, Nc, p::BerendsenThermostat)
    T = md_temperature(v, ms, kb, N, Nc)
    if inv(T) == Inf
        @. dv += p.γ * v
    else
        @. dv += p.γ * (p.T / T - 1) * v
    end
end

# N - number of particles
# Nc - number of constraints
function md_temperature(vs, ms, kb, N, Nc)
    e_kin = sum(dot(ms, vec(sum(vs.^2, dims=1))))
    temperature = e_kin / (kb * (3 * N - Nc))
    return temperature
end

struct NoseHooverThermostat{tType <: Real,τType <: Real} <: Thermostat
    T::tType
    τ::τType
end

function nosehoover_acceleration!(dv, u, v, ms, kb, N, Nc, ζind, p::NoseHooverThermostat)
    @. dv -= u[ζind] * v
    @. dv[:,end] = 0
    T = md_temperature(v[:, 1:N], ms, kb, N, Nc)
    ndf = 3 * N - Nc 
    v[ζind] = inv(p.τ)^2 * (T/p.T - (ndf+1)/ndf)
    #v[ζind] = inv(p.Q) * ( sum(dot(ms, vec(sum(v[:,1:N].^2, 1)))) - (3 * N - Nc) * kb * p.T)
end

struct LangevinThermostat{tType <: Real, gType <: Real} <: Thermostat
    T::tType
    γ::gType
end
