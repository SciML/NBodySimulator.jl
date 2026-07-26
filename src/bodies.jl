"""
    Body

Abstract supertype for particle definitions used by N-body systems.

`Body` implementations provide an initial position, velocity, and mass. The built-in
[`MassBody`](@ref), [`ChargedParticle`](@ref), and [`MagneticParticle`](@ref) types
cover the standard interaction models.
"""
abstract type Body end

"""
    MassBody(r, v, m)

Particle with position `r`, velocity `v`, and mass `m`.

# Arguments

- `r`: Three-dimensional initial position.
- `v`: Three-dimensional initial velocity.
- `m`: Particle mass.

# Fields

- `r`: Initial position.
- `v`: Initial velocity.
- `m`: Mass.

# Examples

```julia
using NBodySimulator, StaticArrays

body = MassBody(SVector(0.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), 1.0)
```
"""
struct MassBody{cType <: Real, mType <: Real} <: Body
    r::SVector{3, cType}
    v::SVector{3, cType}
    m::mType
end

"""
    ChargedParticle(r, v, m, q)

Particle with position `r`, velocity `v`, mass `m`, and charge `q`.

# Arguments

- `r`: Three-dimensional initial position.
- `v`: Three-dimensional initial velocity.
- `m`: Particle mass.
- `q`: Electric charge.

# Fields

- `r`, `v`, `m`: Position, velocity, and mass.
- `q`: Electric charge.

# Examples

```julia
using NBodySimulator, StaticArrays

particle = ChargedParticle(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1.0, 1.0)
```
"""
struct ChargedParticle{cType <: Real, mType <: Real, qType <: Real} <: Body
    r::SVector{3, cType}
    v::SVector{3, cType}
    m::mType
    q::qType
end

"""
    MagneticParticle(r, v, m, mm)

Particle with position `r`, velocity `v`, mass `m`, and magnetic moment `mm`.

# Arguments

- `r`: Three-dimensional initial position.
- `v`: Three-dimensional initial velocity.
- `m`: Particle mass.
- `mm`: Three-dimensional magnetic moment.

# Fields

- `r`, `v`, `m`: Position, velocity, and mass.
- `mm`: Magnetic moment.

# Examples

```julia
using NBodySimulator, StaticArrays

particle = MagneticParticle(
    SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1.0, SVector(0.0, 0.0, 1.0)
)
```
"""
struct MagneticParticle{cType <: Real, mType <: Real, mmType <: Real} <: Body
    r::SVector{3, cType}
    v::SVector{3, cType}
    m::mType
    mm::SVector{3, mmType}
end

struct WaterMolecule <: Body
    O::MassBody
    H1::MassBody
    H2::MassBody
end
"""
    generate_bodies_in_cell_nodes(n, m, v_dev, L; rng = MersenneTwister(n))

Generate `n` equal-mass bodies on a cubic-cell grid with normally distributed
velocities.

# Arguments

- `n`: Number of bodies to generate.
- `m`: Mass assigned to every body.
- `v_dev`: Standard deviation of each velocity component.
- `L`: Side length of the cubic cell.

# Keyword Arguments

- `rng`: Random-number generator used for velocities.

# Examples

```julia
using NBodySimulator

bodies = generate_bodies_in_cell_nodes(8, 1.0, 0.1, 4.0)
```
"""
function generate_bodies_in_cell_nodes(
        n::Int, m::Real, v_dev::Real, L::Real;
        rng = MersenneTwister(n)
    )
    velocities = v_dev * randn(rng, Float64, (3, n))
    bodies = MassBody[]

    count = 1
    dL = L / (ceil(n^(1 / 3)))
    for x in (dL / 2):dL:L, y in (dL / 2):dL:L, z in (dL / 2):dL:L
        if count > n
            break
        end
        r = SVector(x, y, z)
        v = SVector{3}(velocities[:, count])
        body = MassBody(r, v, m)
        push!(bodies, body)
        count += 1
    end
    return bodies
end
