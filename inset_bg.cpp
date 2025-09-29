// inset_bg.cpp
//
// Inset a polygon inward using Boost.Geometry buffer().
// Input:  one point per line, "x y" (whitespace separated), not necessarily closed.
// Output: outer ring of the inset polygon, one "x y" per line (not closed).
//
// Usage:
//   inset_bg <input.txt> <offset> <output.txt> [--join miter|round] [--miter 4.0] [--circle 16]
//
// Notes:
// - <offset> is a positive inward distance (same units as input). We apply -offset in buffer().
// - If the inset collapses (offset too large), program throws.
// - We pick the largest resulting polygon if buffer produces multiple parts.

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/simplify.hpp>
#include <boost/geometry/algorithms/unique.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace bg = boost::geometry;

using Pt       = bg::model::d2::point_xy<double>;
using Polygon  = bg::model::polygon<Pt, /*Clockwise=*/true, /*Closed=*/true>;
using MP       = bg::model::multi_polygon<Polygon>;

namespace { // TU-local helpers

// ---- Hard-coded simplify settings ----
constexpr bool  kSimplifyInput  = true;
constexpr double kSimplifyInputTol  = 0.01; // meters
constexpr bool  kSimplifyOutput = true;
constexpr double kSimplifyOutputTol = 0.01; // meters
// -------------------------------------

struct Args {
  std::string in_path;
  std::string out_path;
  double offset{};             // positive inward distance
  std::string join{"miter"};   // "miter" | "round"
  double miter_limit{4.0};
  int circle_pts{16};
};

double ParseDouble(const std::string& s) {
  std::size_t pos{};
  double v = std::stod(s, &pos);
  if (pos != s.size()) throw std::runtime_error("invalid floating-point: \"" + s + "\"");
  return v;
}
int ParseInt(const std::string& s) {
  std::size_t pos{};
  int v = std::stoi(s, &pos);
  if (pos != s.size()) throw std::runtime_error("invalid integer: \"" + s + "\"");
  return v;
}

Args ParseArgs(int argc, char** argv) {
  if (argc < 4) {
    throw std::runtime_error(
      "usage:\n"
      "  inset_bg <input.txt> <offset> <output.txt>\n"
      "           [--join miter|round] [--miter 4.0] [--circle 16]\n"
      "notes:\n"
      "  <offset> is a positive inward distance; program applies -offset for buffer()\n"
    );
  }
  Args a;
  a.in_path  = std::string(argv[1]);
  a.offset   = ParseDouble(argv[2]); // positive inward
  a.out_path = std::string(argv[3]);

  for (int i = 4; i < argc; ++i) {
    std::string key = argv[i];
    if (key == "--join" && i + 1 < argc) {
      a.join = std::string(argv[++i]);
      if (a.join != "miter" && a.join != "round") {
        throw std::runtime_error("invalid --join (must be miter|round)");
      }
    } else if (key == "--miter" && i + 1 < argc) {
      a.miter_limit = ParseDouble(std::string(argv[++i]));
      if (a.miter_limit <= 0.0) throw std::runtime_error("--miter must be > 0");
    } else if (key == "--circle" && i + 1 < argc) {
      a.circle_pts = ParseInt(std::string(argv[++i]));
      if (a.circle_pts < 4) throw std::runtime_error("--circle must be >= 4");
    } else {
      throw std::runtime_error("unknown or incomplete arg: " + key);
    }
  }

  if (a.offset <= 0.0) throw std::runtime_error("<offset> must be > 0 (inward distance)");
  return a;
}

std::vector<Pt> ReadPoints(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open input: " + path);

  std::vector<Pt> pts;
  pts.reserve(1024);

  double x{}, y{};
  while (in >> x >> y) {
    pts.emplace_back(x, y);
  }
  if (pts.size() < 3) throw std::runtime_error("not enough points in: " + path);
  return pts;
}

Polygon MakePolygonFromPoints(const std::vector<Pt>& pts) {
  Polygon poly;
  auto& ring = poly.outer();

  ring.clear();
  ring.reserve(pts.size() + 1);
  for (const auto& p : pts) ring.push_back(p);

  if (ring.size() >= 2 && !bg::equals(ring.front(), ring.back())) {
    ring.push_back(ring.front()); // close
  }

  bg::correct(poly); // enforce CW + closed for our Polygon typedef
  return poly;
}

const Polygon& BiggestByArea(const MP& mp) {
  if (mp.empty()) throw std::runtime_error("buffer produced no polygons (collapsed?)");
  const Polygon* best = &mp.front();
  auto best_area = std::abs(bg::area(*best));
  for (const auto& poly : mp) {
    auto area = std::abs(bg::area(poly));
    if (area > best_area) {
      best = &poly;
      best_area = area;
    }
  }
  return *best;
}

void WriteOuterClosed(const Polygon& poly, const std::string& path) {
  const auto& ring = poly.outer();
  if (ring.size() < 4) {
    throw std::runtime_error("inset ring degenerate (fewer than 3 unique points)");
  }

  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open output: " + path);

  out.setf(std::ios::fixed);
  out << std::setprecision(6);

  // Write all ring points
  for (std::size_t i = 0; i < ring.size(); ++i) {
    out << ring[i].x() << ' ' << ring[i].y() << '\n';
  }

  // If Boost somehow gave us an open ring, explicitly close it in the file.
  if (ring.front().x() != ring.back().x() || ring.front().y() != ring.back().y()) {
    out << ring.front().x() << ' ' << ring.front().y() << '\n';
  }
}

} // anonymous

int main(int argc, char** argv) {
  try {
    const auto args    = ParseArgs(argc, argv);
    const auto raw_pts = ReadPoints(args.in_path);
    const auto poly_in = MakePolygonFromPoints(raw_pts);

    // Optional input simplify (Douglas–Peucker)
    Polygon work = poly_in;
    if (kSimplifyInput && kSimplifyInputTol > 0.0) {
      Polygon simp;
      bg::simplify(work, simp, kSimplifyInputTol);
      bg::correct(simp);
      work = std::move(simp);
    }

    // Buffer strategies (inset)
    bg::strategy::buffer::distance_symmetric<double> distance(-args.offset); // inward
    bg::strategy::buffer::side_straight              side;
    bg::strategy::buffer::end_flat                   end;
    bg::strategy::buffer::point_circle               circle(args.circle_pts);

    MP result;
    if (args.join == "miter") {
      bg::strategy::buffer::join_miter join(args.miter_limit);
      bg::buffer(work, result, distance, side, join, end, circle);
    } else { // "round"
      bg::strategy::buffer::join_round join;
      bg::buffer(work, result, distance, side, join, end, circle);
    }

    // Pick the biggest output polygon
    const auto& best0 = BiggestByArea(result);
    Polygon best = best0;

    // Optional output simplify
    if (kSimplifyOutput && kSimplifyOutputTol > 0.0) {
      Polygon simp_out;
      bg::simplify(best, simp_out, kSimplifyOutputTol);
      bg::correct(simp_out);
      bg::unique(simp_out);  // drop exact duplicates, if any
      best = std::move(simp_out);
    }

    WriteOuterClosed(best, args.out_path);

    std::cout << "Inset wrote " << (best.outer().size() - 1)
              << " points to: " << args.out_path << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}
