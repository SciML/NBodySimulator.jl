abstract type BoundaryConditions end

"""
    PeriodicBoundaryConditions(boundary)

Periodic boundary conditions over explicit lower and upper bounds for each coordinate.

# Arguments

- `boundary`: Six-element vector `(xmin, xmax, ymin, ymax, zmin, zmax)`.
- `L`: A scalar side length, expanded to the interval `[0, L]` on each axis.

# Fields

- `boundary`: Lower and upper bounds for the three Cartesian coordinates.

# Examples

```julia
using NBodySimulator

boundary = PeriodicBoundaryConditions(10.0)
```
"""
struct PeriodicBoundaryConditions{cType <: Real} <: BoundaryConditions
    boundary::SVector{6, cType}
end

PeriodicBoundaryConditions(L::Real) = PeriodicBoundaryConditions(SVector(0, L, 0, L, 0, L))

# Explicit `show` so the printed representation is independent of the module the
# IO context resolves to (the default struct printer qualifies the type name with
# `NBodySimulator.` when the binding is not visible from the printing module, e.g.
# inside an isolated `@safetestset` module).
function Base.show(io::IO, bc::PeriodicBoundaryConditions{cType}) where {cType}
    print(io, "PeriodicBoundaryConditions{", cType, "}(")
    show(io, bc.boundary)
    return print(io, ")")
end

function Base.iterate(pbc::PeriodicBoundaryConditions, state = 1)
    state > length(pbc.boundary) && return nothing
    return pbc.boundary[state], state + 1
end

function Base.getindex(pbc::PeriodicBoundaryConditions, i::Integer)
    1 <= i <= length(pbc.boundary) || throw(BoundsError(pbc, i))
    return pbc.boundary[i]
end

"""
    InfiniteBox()

Boundary condition representing an unbounded three-dimensional domain.

# Fields

- `boundary`: Six infinite lower and upper bounds.

# Examples

```julia
using NBodySimulator

boundary = InfiniteBox()
```
"""
struct InfiniteBox{cType <: Real} <: BoundaryConditions
    boundary::SVector{6, <:cType}
end

InfiniteBox() = InfiniteBox(SVector(-Inf, Inf, -Inf, Inf, -Inf, Inf))

function Base.show(io::IO, bc::InfiniteBox{cType}) where {cType}
    print(io, "InfiniteBox{", cType, "}(")
    show(io, bc.boundary)
    return print(io, ")")
end

"""
    CubicPeriodicBoundaryConditions(L)

Periodic cubic boundary condition with side length `L`.

# Arguments

- `L`: Cubic-cell side length.

# Fields

- `L`: Cubic-cell side length.

# Examples

```julia
using NBodySimulator

boundary = CubicPeriodicBoundaryConditions(10.0)
```
"""
struct CubicPeriodicBoundaryConditions{cType <: Real} <: BoundaryConditions
    L::cType
end

function Base.show(io::IO, bc::CubicPeriodicBoundaryConditions{cType}) where {cType}
    print(io, "CubicPeriodicBoundaryConditions{", cType, "}(")
    show(io, bc.L)
    return print(io, ")")
end

function get_interparticle_distance(ri, rj, pbc::PeriodicBoundaryConditions)
    rij = ri - rj
    x, y, z = rij
    while x < pbc[1]
        x += pbc[2] - pbc[1]
    end
    while x >= pbc[2]
        x -= pbc[2] - pbc[1]
    end
    while y < pbc[3]
        y += pbc[4] - pbc[3]
    end
    while y >= pbc[4]
        y -= pbc[4] - pbc[3]
    end
    while z < pbc[5]
        z += pbc[6] - pbc[5]
    end
    while z >= pbc[6]
        z -= pbc[6] - pbc[5]
    end
    rij = @SVector [x, y, z]
    r2 = rij[1]^2 + rij[2]^2 + rij[3]^2
    r = sqrt(r2)
    return (rij, r, r2)
end

function get_interparticle_distance(ri, rj, bc::CubicPeriodicBoundaryConditions)
    rij = ri - rj
    x, y, z = rij
    size = bc.L
    radius = 0.5 * size
    while x >= radius
        x -= size
    end
    while x < -radius
        x += size
    end
    while y >= radius
        y -= size
    end
    while y < -radius
        y += size
    end
    while z >= radius
        z -= size
    end
    while z < -radius
        z += size
    end
    rij = @SVector [x, y, z]
    r2 = rij[1]^2 + rij[2]^2 + rij[3]^2
    r = sqrt(r2)
    return (rij, r, r2)
end

function get_interparticle_distance(ri, rj, ::BoundaryConditions)
    rij = ri - rj
    r2 = rij[1]^2 + rij[2]^2 + rij[3]^2
    r = sqrt(r2)
    return (rij, r, r2)
end
