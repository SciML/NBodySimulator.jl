"""
The potentials or force field determines the interaction of particles and, therefore, their acceleration.
"""
abstract type PotentialParameters end

struct LennardJonesParameters{pType <: Real} <: PotentialParameters
    ϵ::pType
    σ::pType
    R::pType
    σ2::pType
    R2::pType
end

function LennardJonesParameters(ϵ::Real, σ::Real, R::Real)
    LennardJonesParameters(ϵ, σ, R, σ^2, R^2)
end

function LennardJonesParameters()
    LennardJonesParameters(1.0, 1.0, 2.5)
end

function Base.show(stream::IO, pp::LennardJonesParameters)
    println(stream, "Lennard-Jones:")
    print(stream, "\tϵ:")
    show(stream, pp.ϵ)
    println(stream)
    print(stream, "\tσ:")
    show(stream, pp.σ)
    println(stream)
    print(stream, "\tR:")
    show(stream, pp.R)
    println(stream)
end

struct GravitationalParameters{gType <: Real} <: PotentialParameters
    G::gType
end

function GravitationalParameters()
    GravitationalParameters(6.67408e-11)
end

function Base.show(stream::IO, pp::GravitationalParameters)
    println(stream, "Gravitational:")
    print(stream, "\tG:")
    show(stream, pp.G)
    println(stream)
end

struct ElectrostaticParameters{pType <: Real} <: PotentialParameters
    k::pType
    R::pType
    R2::pType
end

function ElectrostaticParameters()
    ElectrostaticParameters(9e9, Inf, Inf)
end

function ElectrostaticParameters(k::Real)
    ElectrostaticParameters(k, Inf, Inf)
end

function ElectrostaticParameters(k::Real, R::Real)
    ElectrostaticParameters(k, R, R^2)
end

function Base.show(stream::IO, pp::ElectrostaticParameters)
    println(stream, "Electrostatic:")
    print(stream, "\tk:")
    show(stream, pp.k)
    println(stream)
end

struct MagnetostaticParameters{mType <: Real} <: PotentialParameters
    μ_4π::mType
end

function MagnetostaticParameters()
    MagnetostaticParameters(1e-7)
end

function Base.show(stream::IO, pp::MagnetostaticParameters)
    println(stream, "Magnetostatic:")
    print(stream, "\tμ/4π:")
    show(stream, pp.μ_4π)
    println(stream)
end

struct SPCFwParameters{pType <: Real} <: PotentialParameters
    rOH::pType
    aHOH::pType
    kb::pType
    ka::pType
end

function pairwise_lennard_jones_acceleration!(dv,
        rs,
        i::Integer,
        indxs::Vector{<:Integer},
        ms::Vector{<:Real},
        p::LennardJonesParameters,
        pbc::BoundaryConditions)
    T = eltype(rs)
    force1, force2, force3 = zero(T), zero(T), zero(T)
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]

    for j in indxs
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            (rij, r, rij_2) = get_interparticle_distance(ri, rj, pbc)

            if rij_2 < p.R2
                σ_rij_6 = (p.σ2 / rij_2)^3
                σ_rij_12 = σ_rij_6^2
                factor = (2 * σ_rij_12 - σ_rij_6) / rij_2
                force1 += factor * rij[1]
                force2 += factor * rij[2]
                force3 += factor * rij[3]
            end
        end
    end
    coeff = 24 * p.ϵ / ms[i]
    dv[1] += coeff * force1
    dv[2] += coeff * force2
    dv[3] += coeff * force3
end

function pairwise_electrostatic_acceleration!(dv,
        rs,
        i::Integer,
        n::Integer,
        qs::Vector{<:Real},
        ms::Vector{<:Real},
        exclude::Dict{Int, Vector{Int}},
        p::ElectrostaticParameters,
        pbc::BoundaryConditions)
    T = eltype(rs)
    force1, force2, force3 = zero(T), zero(T), zero(T)
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    @inbounds for j in 1:n
        if !in(j, exclude[i])
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij, r, r2 = get_interparticle_distance(ri, rj, pbc)
            if r2 < p.R2
                factor = qs[j] / (r * r2)
                force1 += factor * rij[1]
                force2 += factor * rij[2]
                force3 += factor * rij[3]
            end
        end
    end
    coeff = p.k * qs[i] / ms[i]
    dv[1] += coeff * force1
    dv[2] += coeff * force2
    dv[3] += coeff * force3
end

function gravitational_acceleration!(dv,
        rs,
        i::Integer,
        n::Integer,
        bodies::Vector{<:MassBody},
        p::GravitationalParameters)
    T = eltype(rs)
    accel1, accel2, accel3 = zero(T), zero(T), zero(T)
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    @inbounds for j in 1:n
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            factor = -p.G * bodies[j].m / norm(rij)^3
            accel1 += factor * rij[1]
            accel2 += factor * rij[2]
            accel3 += factor * rij[3]
        end
    end

    dv[1] += accel1
    dv[2] += accel2
    dv[3] += accel3
end

function magnetostatic_dipdip_acceleration!(dv,
        rs,
        i::Integer,
        n::Integer,
        bodies::Vector{<:MagneticParticle},
        p::MagnetostaticParameters)
    T = eltype(rs)
    force1, force2, force3 = zero(T), zero(T), zero(T)
    mi = bodies[i].mm
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    @inbounds for j in 1:n
        if j != i
            mj = bodies[j].mm
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            rij4 = dot(rij, rij)^2
            r = rij / norm(rij)
            mir = dot(mi, r)
            mij = dot(mj, r)
            contrib = (mi * mij + mj * mir + r * dot(mi, mj) - 5 * r * mir * mij) / rij4
            force1 += contrib[1]
            force2 += contrib[2]
            force3 += contrib[3]
        end
    end

    coeff = 3 * p.μ_4π / bodies[i].m
    dv[1] += coeff * force1
    dv[2] += coeff * force2
    dv[3] += coeff * force3
end

function harmonic_bond_potential_acceleration!(dv,
        rs,
        i::Integer,
        ms::Vector{<:Real},
        neighbouhoods,
        p::SPCFwParameters)
    T = eltype(rs)
    force1, force2, force3 = zero(T), zero(T), zero(T)
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    @inbounds for (j, k) in neighbouhoods[i]
        rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
        rij = ri - rj
        r = norm(rij)
        d = r - p.rOH
        factor = -d * k / r
        force1 += factor * rij[1]
        force2 += factor * rij[2]
        force3 += factor * rij[3]
    end

    coeff = one(T) / ms[i]
    dv[1] += coeff * force1
    dv[2] += coeff * force2
    dv[3] += coeff * force3
end

function valence_angle_potential_acceleration!(dv,
        rs,
        a::Integer,
        b::Integer,
        c::Integer,
        ms::Vector{<:Real},
        p::SPCFwParameters)
    ra = @SVector [rs[1, a], rs[2, a], rs[3, a]]
    rb = @SVector [rs[1, b], rs[2, b], rs[3, b]]
    rc = @SVector [rs[1, c], rs[2, c], rs[3, c]]

    rba = ra - rb
    rbc = rc - rb
    rcb = rb - rc

    rbaXbc = cross(rba, rbc)
    pa = normalize(cross(rba, rbaXbc))
    pc = normalize(cross(rcb, rbaXbc))

    cosine = dot(rba, rbc) / (norm(rba) * norm(rbc))

    if cosine > 1
        cosine = 1
    elseif cosine < -1
        cosine = -1
    end
    aHOH = acos(cosine)

    force = -p.ka * (aHOH - p.aHOH)
    force_a = pa * force / norm(rba)
    force_c = pc * force / norm(rbc)
    force_b = -(force_a + force_c)

    @. dv[:, a] += force_a / ms[a]
    @. dv[:, b] += force_b / ms[b]
    @. dv[:, c] += force_c / ms[c]
end
