abstract type PotentialParameters
end

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
    print(stream, "\tϵ:"); show(stream, pp.ϵ); println(stream)
    print(stream, "\tσ:"); show(stream, pp.σ); println(stream)
    print(stream, "\tR:"); show(stream, pp.R); println(stream)
end

struct GravitationalParameters{gType <: Real} <: PotentialParameters
    G::gType
end

function GravitationalParameters()
    GravitationalParameters(6.67408e-11)
end

function Base.show(stream::IO, pp::GravitationalParameters)
    println(stream, "Gravitational:")
    print(stream, "\tG:"); show(stream, pp.G);
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
    print(stream, "\tk:"); show(stream, pp.k);
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
    print(stream, "\tμ/4π:"); show(stream, pp.μ_4π);
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

    force = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]

    for j ∈ indxs
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            (rij, r, rij_2) = get_interparticle_distance(ri, rj, pbc)

            if rij_2 < p.R2
                σ_rij_6 = (p.σ2 / rij_2)^3
                σ_rij_12 = σ_rij_6^2
                force += (2 * σ_rij_12 - σ_rij_6 ) * rij / rij_2
            end
        end
    end
    @. dv +=  24 * p.ϵ * force / ms[i]
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

    force = @SVector [0.0, 0.0, 0.0]
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    @inbounds for j = 1:n
        if !in(j, exclude[i])
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij, r, r2 = get_interparticle_distance(ri, rj, pbc)
            if r2 < p.R2
                force += qs[j] * rij / (r*r2)
            end
        end
    end
    @. dv += p.k * qs[i] * force / ms[i]
end


function gravitational_acceleration!(dv,
    rs,
    i::Integer,
    n::Integer,
    bodies::Vector{<:MassBody},
    p::GravitationalParameters)

    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    @inbounds for j = 1:n
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            accel -= p.G * bodies[j].m * rij / norm(rij)^3
        end
    end

    @. dv += accel
end

function magnetostatic_dipdip_acceleration!(dv,
    rs,
    i::Integer,
    n::Integer,
    bodies::Vector{<:MagneticParticle},
    p::MagnetostaticParameters)

    force = @SVector [0.0, 0.0, 0.0];
    mi = bodies[i].mm
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    @inbounds for j = 1:n
        if j != i
            mj = bodies[j].mm
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            rij4 = dot(rij, rij)^2
            r =  rij / norm(rij)
            mir = dot(mi, r)
            mij = dot(mj, r)
            force += (mi * mij + mj * mir + r * dot(mi, mj) - 5 * r * mir * mij) / rij4
        end
    end

    @. dv += 3 * p.μ_4π * force / bodies[i].m
end

function harmonic_bond_potential_acceleration!(dv,
    rs,
    i::Integer,
    ms::Vector{<:Real},
    neighbouhoods::Dict{Int,Vector{Tuple{Int,Float64}}},
    p::SPCFwParameters)

    force = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    @inbounds for (j, k) in neighbouhoods[i]
        rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
        rij = ri - rj
        r = norm(rij)
        d = r - p.rOH
        force -= d * k * rij / r
    end

    @. dv += force / ms[i]
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

    if cosine>1
        cosine = 1
    elseif cosine<-1
        cosine =-1
    end
    aHOH = acos(cosine)

    force = - p.ka * (aHOH - p.aHOH)
    force_a = pa * force / norm(rba)
    force_c = pc * force / norm(rbc)
    force_b = -(force_a + force_c)

    @. dv[:,a] += force_a / ms[a]
    @. dv[:,b] += force_b / ms[b]
    @. dv[:,c] += force_c / ms[c]
end
