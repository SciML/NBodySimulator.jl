abstract type BoundaryConditions
end

struct PeriodicBoundaryConditions{cType <: Real} <: BoundaryConditions
    boundary::SVector{6,cType}
end

PeriodicBoundaryConditions(L::Real) = PeriodicBoundaryConditions(SVector(0, L, 0, L, 0, L))

Base.start(::PeriodicBoundaryConditions) = 1

Base.done(pbc::PeriodicBoundaryConditions, state) = state > length(pbc.boundary)

function Base.next(pbc::PeriodicBoundaryConditions, state) 
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

function apply_boundary_conditions!(ri, rj, pbc::PeriodicBoundaryConditions, R2)
    rij = @MVector [Inf, Inf, Inf]
    success = false
    rij2 = 0

    for dx in [0,-pbc[2],pbc[2]], dy in [0,-pbc[4],pbc[4]], dz in [0,-pbc[6],pbc[6]]
        rij = @MVector [ri[1] - rj[1] + dx, ri[2] - rj[2] + dy, ri[3] - rj[3] + dz]
        for x in (1, 2, 3) 
            rij[x] -= (pbc[2x] - pbc[2x - 1]) * div(rij[x], (pbc[2x] - pbc[2x - 1]))
        end
        rij2 = dot(rij, rij)
        if  rij2 < R2
            success = true
            break
        end
    end
    return (rij, rij2, success)
end

function apply_boundary_conditions!(ri, rj, pbc::CubicPeriodicBoundaryConditions, R2)
    rij = @MVector [Inf, Inf, Inf]
    success = false
    rij2 = 0

    shifts = @SVector [0,-pbc.L,pbc.L]
    for dx in shifts, dy in shifts, dz in shifts
        rij = @MVector [ri[1] - rj[1] + dx, ri[2] - rj[2] + dy, ri[3] - rj[3] + dz]
        for x in (1, 2, 3) 
            rij[x] -= pbc.L * div(rij[x], pbc.L)
        end
        rij2 = dot(rij, rij)
        if  rij2 < R2
            success = true
            break
        end
    end
    return (rij, rij2, success)
end

function apply_boundary_conditions!(ri, rj, pbc::BoundaryConditions, R2)
    rij = ri - rj
    (rij, dot(rij,rij),true)
end