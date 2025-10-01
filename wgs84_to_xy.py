#!/usr/bin/env python3
# wgs84_to_xy.py
# Convert WGS-84 lat/lon (deg) to local EN (x,y in meters) on stdout.
# Reads from stdin: one "lat lon" per line. Lines starting with '#' and
# blank lines are passed through unchanged.
#
# Origin (deg): +44.789915759, -95.942430803

import sys
import math


ORIGIN_LAT_DEG =  44.789915759
ORIGIN_LON_DEG = -95.942430803


def _geodetic_to_ecef(lat_rad, lon_rad, h_m=0.0):
    """Convert WGS-84 geodetic to ECEF (meters)."""
    a = 6378137.0                # semi-major axis (m)
    f = 1.0 / 298.257223563      # flattening
    e2 = f * (2.0 - f)           # first eccentricity squared

    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)

    N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)

    x = (N + h_m) * cos_lat * cos_lon
    y = (N + h_m) * cos_lat * sin_lon
    z = (N * (1.0 - e2) + h_m) * sin_lat
    return x, y, z


def _ecef_to_enu(dx, dy, dz, lat0_rad, lon0_rad):
    """Rotate ECEF delta (m) into local ENU at origin (radians)."""
    sin_lat0 = math.sin(lat0_rad)
    cos_lat0 = math.cos(lat0_rad)
    sin_lon0 = math.sin(lon0_rad)
    cos_lon0 = math.cos(lon0_rad)

    e = -sin_lon0 * dx + cos_lon0 * dy
    n = (-sin_lat0 * cos_lon0 * dx
         -sin_lat0 * sin_lon0 * dy
         +cos_lat0 * dz)
    return e, n  # up (u) omitted


def main():
    lat0_rad = math.radians(ORIGIN_LAT_DEG)
    lon0_rad = math.radians(ORIGIN_LON_DEG)

    x0, y0, z0 = _geodetic_to_ecef(lat0_rad, lon0_rad, 0.0)

    for raw in sys.stdin:
        line = raw.rstrip("\n")

        if not line:
            print()
            continue

        if line.lstrip().startswith('#'):
            print(line)
            continue

        parts = line.split()
        if len(parts) < 2:
            print(f"# skipped: {line}")
            continue

        try:
            lat_deg = float(parts[0])
            lon_deg = float(parts[1])
        except ValueError:
            print(f"# skipped: {line}")
            continue

        lat_rad = math.radians(lat_deg)
        lon_rad = math.radians(lon_deg)

        x, y, z = _geodetic_to_ecef(lat_rad, lon_rad, 0.0)
        dx, dy, dz = (x - x0), (y - y0), (z - z0)

        e, n = _ecef_to_enu(dx, dy, dz, lat0_rad, lon0_rad)
        print(f"{e:.6f} {n:.6f}")


if __name__ == "__main__":
    main()
