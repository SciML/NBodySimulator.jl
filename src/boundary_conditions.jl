abstract type BoundaryConditions
end

struct PeriodicBoundaryConditions{cType <: Real} <: BoundaryConditions
    boundary::SVector{6,cType}
end

PeriodicBoundaryConditions(L::Real) = PeriodicBoundaryConditions(SVector(0, L, 0, L, 0, L))

function Base.iterate(pbc::PeriodicBoundaryConditions,state=1)
  state > length(pbc.boundary) && return nothing
  pbc.boundary[state], state + 1
end

function Base.getindex(pbc::PeriodicBoundaryConditions, i::Integer)
    1 <= i <= length(pbc.boundary) || throw(BoundsError(pbc, i))
    pbc.boundary[i]
end

struct InfiniteBox{cType <: Real} <: BoundaryConditions
    boundary::SVector{6,<:cType}
end

InfiniteBox() = InfiniteBox(SVector(-Inf, Inf, -Inf, Inf, -Inf, Inf))

struct CubicPeriodicBoundaryConditions{cType <: Real} <: BoundaryConditions
    L::cType
end

function get_interparticle_distance(ri, rj, pbc::PeriodicBoundaryConditions)
    rij = ri - rj
    x, y, z = rij
    while x  < pbc[1]   x += pbc[2]-pbc[1] end
    while x >= pbc[2]   x -= pbc[2]-pbc[1] end
    while y  < pbc[3]   y += pbc[4]-pbc[3] end
    while y >= pbc[4]   y -= pbc[4]-pbc[3] end
    while z  < pbc[5]   z += pbc[6]-pbc[5] end
    while z >= pbc[6]   z -= pbc[6]-pbc[5] end
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
    while x >= radius    x -= size end
    while x < -radius    x += size end
    while y >= radius    y -= size end
    while y < -radius    y += size end
    while z >= radius    z -= size end
    while z < -radius    z += size end
    rij = @SVector [x, y, z]
    r2 = rij[1]^2 + rij[2]^2 + rij[3]^2
    r = sqrt(r2)
    return (rij, r, r2)
end

function get_interparticle_distance(ri, rj, ::BoundaryConditions)
    rij = ri - rj
    r2 = rij[1]^2 + rij[2]^2 + rij[3]^2
    r = sqrt(r2)
    (rij, r, r2)
end
