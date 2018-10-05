#!/usr/bin/env julia

include("ecef.jl")

function ecef2enu(x::Real, y::Real, z::Real, lat0::Real, lon0::Real, h0::Real, ell::Ellipsoid)
    #=
    input
    -----
    x,y,z  [meters] target ECEF location                             [0,Infinity)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid

    output:
    -------
    e,n,u   East, North, Up [m]
    =#
    
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell)

    return uvw2enu(x - x0, y - y0, z - z0, lat0, lon0)
end


function uvw2enu(u::Real, v::Real, w::Real, lat0::Real, lon0::Real)
    
    lat0 = deg2rad(lat0)
    lon0 = deg2rad(lon0)

    t = cos(lon0) * u + sin(lon0) * v
    East = -sin(lon0) * u + cos(lon0) * v
    Up = cos(lat0) * t + sin(lat0) * w
    North = -sin(lat0) * t + cos(lat0) * w

    return (East, North, Up)
end


function aer2enu(az::Real, el::Real, srange::Real)
    #=
    input:
    ------
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    deg    degrees input/output  (False: radians in/out)

    output:
    -------
    e,n,u   East, North, Up [m]

    =#
    
    el = deg2rad(el)
    az = deg2rad(az)

    r = srange * cos(el)

    return (r * sin(az), r * cos(az), srange * sin(el))
end


if basename(PROGRAM_FILE) == basename(@__FILE__)

  if length(ARGS) == 6
  
  else
    ell = Ellipsoid(6378137., 1 / 298.2572235630, 6378137. * (1 - 1 / 298.2572235630))
    
    lla0 = (42., -82., 200.)
    aer0 = (33., 70., 1000.)
    xyz = (660930.1927610816, -4.701424222957011e6, 4.246579604632881e6)
    east0, north0, up0 = (186.277521, 286.842228, 939.692621)
    
    east, north, up = ecef2enu(xyz..., lla0..., ell)
    
    @test isapprox(east, east0)
    @test isapprox(north, north0)
    @test isapprox(up, up0)
    
  end
end
