#!/usr/bin/env julia
using Test

struct Ellipsoid
    a::Float64  # semi-major axis [m]
    f::Float64  # flattening
    b::Float64  # semi-minor axis [m]
end

function ecef2geodetic(x::Real, y::Real, z::Real, ell::Ellipsoid)
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

    based on:
    You, Rey-Jer. (2000). Transformation of Cartesian to Geodetic Coordinates without Iterations.
    Journal of Surveying Engineering. doi: 10.1061/(ASCE)0733-9453
    =#

    r = sqrt(x^2 + y^2 + z^2)

    E = sqrt(ell.a^2 - ell.b^2)

    # eqn. 4a
    u = sqrt(0.5 * (r^2 - E^2) + 0.5 * sqrt((r^2 - E^2)^2 + 4 * E^2 * z^2))

    Q = hypot(x, y)

    huE = hypot(u, E)

    # eqn. 4b
    Beta = atan(huE / u * z / hypot(x, y))

    # eqn. 13
    eps = ((ell.b * u - ell.a * huE + E^2) * sin(Beta)) / (ell.a * huE * 1 / cos(Beta) - E^2 * cos(Beta))

    Beta += eps
# %% final output
    lat = atan(ell.a / ell.b * tan(Beta))

    lon = atan(y, x)

    # eqn. 7
    alt = sqrt((z - ell.b * sin(Beta))^2 + (Q - ell.a * cos(Beta))^2)


    lat > -pi / 2 && lat < pi / 2 || throw(DomainError("-90 <= lat <= 90"))
    lon > -pi && lon < 2 * pi || throw(DomainError("-180 <= lat <= 360"))

    return (rad2deg(lat), rad2deg(lon), alt)

end

function geodetic2ecef(lat::Real, lon::Real, alt::Real, ell::Ellipsoid)
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

    return (x, y, z)
end


function enu2uvw(east::Real, north::Real, up::Real,
                 lat0::Real, lon0::Real)

    lat0 = deg2rad(lat0)
    lon0 = deg2rad(lon0)

    t = cos(lat0) * up - sin(lat0) * north
    w = sin(lat0) * up + cos(lat0) * north

    u = cos(lon0) * t - sin(lon0) * east
    v = sin(lon0) * t + cos(lon0) * east

    return (u, v, w)
end


function get_radius_normal(lat_radians, ell)
    # Compute normal radius of planetary body

    a = ell.a
    b = ell.b

    return a^2 / sqrt(a^2 * cos(lat_radians)^2 + b^2 * sin(lat_radians)^2)
end


if basename(PROGRAM_FILE) == basename(@__FILE__)

  if length(ARGS) == 3
  
  else
    lla0 = (42., -82., 200.)
    ell = Ellipsoid(6378137., 1 / 298.2572235630, 6378137. * (1 - 1 / 298.2572235630))
    
    x, y, z = geodetic2ecef(lla0..., ell)
    
    @test isapprox(x, 660675.2518247)
    @test isapprox(y, -4700948.68316)
    @test isapprox(z, 4245737.66222)
    
    lat1, lon1, alt1 = ecef2geodetic(x,y,z, ell)
    
    @test isapprox(lat1, lla0[1])
    @test isapprox(lon1, lla0[2])
    @test isapprox(alt1, lla0[3])
    
  end
end
