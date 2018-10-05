#!/usr/bin/env julia

include("ecef.jl")
include("enu.jl")

function aer2ecef(az::Real, el::Real, srange::Real,
                  lat0::Real, lon0::Real, alt0::Real,
                  ell::Ellipsoid)
    #=
    convert target azimuth, elevation, range (meters) from observer at lat0,lon0,alt0 to ECEF coordinates.

     Input:
     -----
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

     output: ECEF x,y,z  [meters]

    if you specify NaN for srange, return value z will be NaN
    =#
    
    # Origin of the local system in geocentric coordinates.
    x0, y0, z0 = geodetic2ecef(lat0, lon0, alt0, ell)
    # Convert Local Spherical AER to ENU
    e1, n1, u1 = aer2enu(az, el, srange)
    # Rotating ENU to ECEF
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0)
    # Origin + offset from origin equals position in ECEF
    return (x0 + dx, y0 + dy, z0 + dz)
end


if basename(PROGRAM_FILE) == basename(@__FILE__)
    lla0 = (42., -82., 200.)
    aer0 = (33., 70., 1000.)    
    ell = Ellipsoid(6378137., 1 / 298.2572235630, 6378137. * (1 - 1 / 298.2572235630))
    
    println(aer2ecef(aer0..., lla0..., ell))
end
