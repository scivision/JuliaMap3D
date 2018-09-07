#!/usr/bin/env julia
using Test

struct Ellipsoid
    a::Float64  # semi-major axis [m]
    f::Float64  # flattening
    b::Float64  # semi-minor axis [m]
end

function ecef2geodetic(x, y, z, ell::Ellipsoid)
    #=
    convert ECEF (meters) to geodetic coordinates

    input
    -----
    x,y,z  [meters] target ECEF location                             [0,Infinity)
    ell    reference ellipsoid

    output
    ------
    lat,lon   (degrees/radians)
    alt  (meters)

    Algorithm is based on
    http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
    This algorithm provides a converging solution to the latitude equation
    in terms of the parametric or reduced latitude form (v)
    This algorithm provides a uniform solution over all latitudes as it does
    not involve division by cos(phi) or sin(phi)
    =#

    ea = ell.a
    eb = ell.b
    rad = hypot(x, y)
# Constant required for Latitude equation
    rho = atan(eb * z, ea * rad)
# Constant required for latitude equation
    c = (ea^2 - eb^2) / hypot(ea * rad, eb * z)
# Starter for the Newtons Iteration Method
    vnew = atan(ea * z, eb * rad)
# Initializing the parametric latitude
    v = 0
    for i = 1:5
        v = deepcopy(vnew)
# %% Newtons Method for computing iterations
        vnew = v - ((2 * sin(v - rho) - c * sin(2 * v)) /
                    (2 * (cos(v - rho) - c * cos(2 * v))))

        if all(v ≈ vnew)
            break
        end
    end
# %% Computing latitude from the root of the latitude equation
    lat = atan(ea * tan(vnew), eb)
    # by inspection
    lon = atan(y, x)

    alt = (((rad - ea * cos(vnew)) * cos(lat)) +
           ((z - eb * sin(vnew)) * sin(lat)))

    lat > -pi / 2 && lat < pi / 2 || throw(DomainError("-90 <= lat <= 90"))
    lon > -pi && lon < 2 * pi || throw(DomainError("-180 <= lat <= 360"))

    return [rad2deg(lat), rad2deg(lon), alt]

end

function geodetic2ecef(lat, lon, alt, ell::Ellipsoid)
    #=
    Point
    input:
    -----
    lat, lon (degrees)
    alt (altitude, meters)    [0, Infinity)
    ell    reference ellipsoid
    
    output: ECEF x,y,z (meters)
    =#

    lat = deg2rad(lat)
    lon = deg2rad(lon)
    
    lat > -pi / 2 && lat < pi / 2 || throw(DomainError("-90 <= lat <= 90"))

    lon > -pi && lon < 2 * pi || throw(DomainError("-180 <= lat <= 360"))

    alt > 0 || throw(DomainError("altitude in  [0, Infinity)"))
    # radius of curvature of the prime vertical section
    N = get_radius_normal(lat, ell)
    # Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic
    # coordinates.
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (ell.b / ell.a)^2 + alt) * sin(lat)

    return [x, y, z]
end

function get_radius_normal(lat_radians, ell)
    # Compute normal radius of planetary body

    a = ell.a
    b = ell.b

    return a^2 / sqrt(a^2 * cos(lat_radians)^2 + b^2 * sin(lat_radians)^2)
end


if PROGRAM_FILE == splitdir(@__FILE__)[end]
  lla0 = [42., -82., 200.]
  ell = Ellipsoid(6378137., 1 / 298.2572235630, 6378137. * (1 - 1 / 298.2572235630))
  xyz0 = geodetic2ecef(lla0..., ell)
  
  @test isapprox(xyz0,[660.6753e3,-4700.9487e3,4245.738e3], rtol=1e-4)
  
  lla1 = ecef2geodetic(xyz0..., ell)
  @test isapprox(lla1, lla0)
  
end
