#!/usr/bin/env python3
# Convert WGS-84 lat/lon[/h] to local ENU (East, North) XY coordinates in meters.
# Origin is either the first coordinate in the file or explicitly provided.

import argparse
import math
import sys
from typing import List, Tuple, Optional

# WGS-84 ellipsoid constants
WGS84_A = 6378137.0            # semi-major axis (m)
WGS84_F = 1.0 / 298.257223563  # flattening
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)  # eccentricity squared

def deg2rad(d: float) -> float:
    return d * math.pi / 180.0

def llh_to_ecef(lat_deg: float, lon_deg: float, h_m: float = 0.0) -> Tuple[float, float, float]:
    """Geodetic (lat, lon, h) -> ECEF (X, Y, Z)."""
    lat = deg2rad(lat_deg)
    lon = deg2rad(lon_deg)
    s = math.sin(lat)
    c = math.cos(lat)
    chi = math.sqrt(1.0 - WGS84_E2 * s * s)
    N = WGS84_A / chi
    X = (N + h_m) * c * math.cos(lon)
    Y = (N + h_m) * c * math.sin(lon)
    Z = (N * (1.0 - WGS84_E2) + h_m) * s
    return X, Y, Z

def enu_rotation(lat0_deg: float, lon0_deg: float) -> List[List[float]]:
    """Rotation matrix R so [e n u]^T = R * ([X Y Z]^T - r0)."""
    lat0 = deg2rad(lat0_deg)
    lon0 = deg2rad(lon0_deg)
    sL = math.sin(lat0); cL = math.cos(lat0)
    sO = math.sin(lon0); cO = math.cos(lon0)
    # East, North, Up unit vectors in ECEF
    rE = (-sO,         cO,       0.0)
    rN = (-sL*cO,     -sL*sO,    cL)
    rU = ( cL*cO,      cL*sO,    sL)
    return [list(rE), list(rN), list(rU)]

def mat_vec_mul3(R: List[List[float]], v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2],
        R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2],
        R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2],
    )

def parse_points(path: str) -> List[Tuple[float, float, float]]:
    pts: List[Tuple[float, float, float]] = []
    with open(path, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                raise ValueError(f"{path}:{ln}: expected 'lat lon [h]'")
            try:
                lat = float(parts[0])
                lon = float(parts[1])
                h = float(parts[2]) if len(parts) >= 3 else 0.0
            except ValueError:
                raise ValueError(f"{path}:{ln}: invalid number")
            if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
                raise ValueError(f"{path}:{ln}: lat/lon out of range")
            pts.append((lat, lon, h))
    if not pts:
        raise ValueError(f"{path}: no valid coordinates found")
    return pts

def write_xy(path: str, xy: List[Tuple[float, float]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for x, y in xy:
            f.write(f"{x:.6f} {y:.6f}\n")

def convert_ll_to_local_xy(
    pts_llh: List[Tuple[float, float, float]],
    origin: Optional[Tuple[float, float, float]] = None
) -> Tuple[List[Tuple[float, float]], Tuple[float, float, float]]:
    """Convert list of (lat,lon,h) to ENU XY (east,north)."""
    if origin is None:
        lat0, lon0, h0 = pts_llh[0]
    else:
        lat0, lon0, *rest = origin
        h0 = rest[0] if rest else 0.0

    X0, Y0, Z0 = llh_to_ecef(lat0, lon0, h0)
    R = enu_rotation(lat0, lon0)

    xy: List[Tuple[float, float]] = []
    for (lat, lon, h) in pts_llh:
        X, Y, Z = llh_to_ecef(lat, lon, h)
        dX, dY, dZ = X - X0, Y - Y0, Z - Z0
        e, n, _u = mat_vec_mul3(R, (dX, dY, dZ))
        xy.append((e, n))  # east, north
    return xy, (lat0, lon0, h0)

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert WGS-84 lat lon [h] to local ENU XY (meters)."
    )
    ap.add_argument("input", help="input file: 'lat lon' or 'lat lon h' per line")
    ap.add_argument("output", help="output file: 'x y' per line (meters)")
    ap.add_argument("--origin", nargs=2, type=float, metavar=("LAT0", "LON0"),
                    help="origin (lat lon) in degrees; default: first point")
    ap.add_argument("--origin-h", type=float, default=0.0,
                    help="origin height (m), only used with --origin")
    args = ap.parse_args()

    try:
        pts = parse_points(args.input)
        origin = (args.origin[0], args.origin[1], args.origin_h) if args.origin else None
        xy, (lat0, lon0, h0) = convert_ll_to_local_xy(pts, origin)
        write_xy(args.output, xy)
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        return 1

    sys.stdout.write(
        f"Converted {len(xy)} points with origin "
        f"lat0={lat0:.8f}, lon0={lon0:.8f}, h0={h0:.3f} m\n"
        f"Wrote: {args.output}\n"
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())