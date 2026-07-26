"""
Abstract supertype for thermostat configurations.

Thermostats control a simulation's temperature. Use one of the concrete exported
thermostat types when constructing an [`NBodySimulation`](@ref).
"""
abstract type Thermostat end
"""
No thermostat.
"""
struct NullThermostat <: Thermostat
end

"""
    AndersenThermostat(T, ν)

Andersen thermostat targeting temperature `T` with collision frequency `ν`.

# Arguments

- `T`: Target temperature.
- `ν`: Collision frequency.

# Fields

- `T`: Target temperature.
- `ν`: Collision frequency.

# Examples

```julia
using NBodySimulator

thermostat = AndersenThermostat(300.0, 1.0)
```
"""
struct AndersenThermostat{tType <: Real, νType <: Real} <: Thermostat
    T::tType
    ν::νType
end

"""
    BerendsenThermostat(T, τ)

Berendsen thermostat targeting temperature `T` with coupling time `τ`.

# Arguments

- `T`: Target temperature.
- `τ`: Temperature-coupling time.

# Fields

- `T`: Target temperature.
- `τ`: Temperature-coupling time.
- `γ`: Rescaling coefficient derived from `τ`.

# Examples

```julia
using NBodySimulator

thermostat = BerendsenThermostat(300.0, 0.1)
```
"""
struct BerendsenThermostat{τType <: Real, tType <: Real} <: Thermostat
    T::tType
    τ::τType
    γ::τType
end

function BerendsenThermostat(T::Real, τ::Real)
    return BerendsenThermostat(T, τ, 0.5 / τ)
end

function berendsen_acceleration!(dv, v, ms, kb, N, Nc, p::BerendsenThermostat)
    T = md_temperature(v, ms, kb, N, Nc)
    return if inv(T) == Inf
        @. dv += p.γ * v
    else
        @. dv += p.γ * (p.T / T - 1) * v
    end
end

# N - number of particles
# Nc - number of constraints
function md_temperature(vs, ms, kb, N, Nc)
    e_kin = sum(dot(ms, vec(sum(vs .^ 2, dims = 1))))
    temperature = e_kin / (kb * (3 * N - Nc))
    return temperature
end

"""
    NoseHooverThermostat(T, τ)

Nose-Hoover thermostat targeting temperature `T` with relaxation time `τ`.

# Arguments

- `T`: Target temperature.
- `τ`: Relaxation time.

# Fields

- `T`: Target temperature.
- `τ`: Relaxation time.

# Examples

```julia
using NBodySimulator

thermostat = NoseHooverThermostat(300.0, 0.1)
```
"""
struct NoseHooverThermostat{tType <: Real, τType <: Real} <: Thermostat
    T::tType
    τ::τType
end

function nosehoover_acceleration!(dv, u, v, ms, kb, N, Nc, ζind, p::NoseHooverThermostat)
    @. dv -= u[ζind] * v
    @. dv[:, end] = 0
    T = md_temperature(v[:, 1:N], ms, kb, N, Nc)
    ndf = 3 * N - Nc
    return v[ζind] = inv(p.τ)^2 * (T / p.T - (ndf + 1) / ndf)
    #v[ζind] = inv(p.Q) * ( sum(dot(ms, vec(sum(v[:,1:N].^2, 1)))) - (3 * N - Nc) * kb * p.T)
end

"""
    LangevinThermostat(T, γ)

Langevin thermostat targeting temperature `T` with friction coefficient `γ`.

# Arguments

- `T`: Target temperature.
- `γ`: Friction coefficient.

# Fields

- `T`: Target temperature.
- `γ`: Friction coefficient.

# Examples

```julia
using NBodySimulator

thermostat = LangevinThermostat(300.0, 1.0)
```
"""
struct LangevinThermostat{tType <: Real, gType <: Real} <: Thermostat
    T::tType
    γ::gType
end
