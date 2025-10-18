#include "BoundarySwaths.hpp"

#include <mp-units/systems/si/unit_symbols.h>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/algorithms/correct.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <gsl-lite/gsl-lite.hpp>
namespace gsl = gsl_lite;

#include <algorithm>
#include <string>
#include <vector>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <exception>
#include <type_traits>
#include <cstdlib>

template<typename T>
void DbgOut(const char* str, const T& x)
  { std::cerr << str << '=' << x << std::endl; }

#define DBG(x) DbgOut(BOOST_STRINGIZE(x), (x))

namespace fs  = std::filesystem;
namespace ggl = boost::geometry;

namespace tjg {

using Ring = ggl::model::ring<Pt>;

} // tjg

namespace {

std::ostream& operator<<(std::ostream& os, const tjg::Pt& p)
  { return os << p.x << ' ' << p.y; }

std::istream& operator>>(std::istream& os, tjg::Pt& p)
  { return os >> p.x >> p.y; }

#ifdef TJG_NOT_USED
std::ostream& operator<<(std::ostream& os, const tjg::Ring& r) {
  for (const auto& p: r)
    os << p << '\n';
  return os;
}
#endif

template<class Geo>
requires (!std::is_same_v<Geo, tjg::Ring>)
void EnsureValid(const Geo& geo) {
  auto failure = ggl::validity_failure_type{};
  if (ggl::is_valid(geo, failure)) [[likely]]
    return;
  auto msg = std::string{"Invalid geometry: "};
  msg += ggl::validity_failure_type_message(failure);
  throw std::runtime_error{msg};
} // EnsureValid

tjg::Polygon ReadPolygon(const fs::path& path) {
  auto in = std::ifstream{path, std::ios::binary};
  if (!in) throw std::runtime_error{"cannot open input: " + path.string()};
  auto poly = tjg::Polygon{};
  auto* ring = &poly.outer();
  auto pts = std::vector<tjg::Pt>{};
  auto line = std::string{};
  auto pt = tjg::Pt{};
  while (getline(in, line)) {
    if (line.empty())
      continue;
    if (line != "# inner") {
      auto iss = std::istringstream{line};
      iss >> pt;
      ring->push_back(pt);
      continue;
    }
    // Start a new inner ring.
    poly.inners().emplace_back();
    ring = &poly.inners().back();
  }
  in.close();
  ggl::correct(poly);
  EnsureValid(poly);
  return poly;
} // ReadPolygon

// Writes all the swaths associated with a single polygon.
void WriteSwath(const fs::path& path, const tjg::MultiPath& swath) {
  if (swath.empty())
    return;
  const auto ext = path.extension();
  int num = 0;
  for (const auto& pth: swath) {
    auto path2 = path;
    auto ext2  = fs::path{std::to_string(++num)};
    ext2 += ext;
    path2.replace_extension(ext2);
    auto out = std::ofstream{path2, std::ios::binary};
    if (!out) throw std::runtime_error{"cannot open output: " + path2.string()};
    out << std::fixed << std::setprecision(2);
    for (const auto& pt: pth) out << pt << '\n';
    out.close();
  }
} // WriteSwath

// Writes all swaths associated with a field.
void WriteSwaths(const fs::path& path, const std::vector<tjg::MultiPath>& swaths)
{
  if (swaths.empty())
    return;
  const auto ext = path.extension();
  auto num = 0;
  for (const auto& swath: swaths) {
    auto path2 = path;
    auto ext2  = fs::path{std::to_string(++num)};
    ext2 += ext;
    path2.replace_extension(ext2);
    WriteSwath(path2, swath);
  }
} // WriteSwaths

} // anonymous

int main(int argc, const char* argv[]) {
  using namespace mp_units::si::unit_symbols;
  try {
    auto arg0 = fs::path{argv[0]};

    if (argc != 4) {
      std::cerr << "usage:\n" << arg0.string()
        << "  <input.xy> <offset_m> <output.xy>\n";
      return EXIT_FAILURE;
    }

    const auto in_path  = fs::path{argv[1]};
    const auto offset   = std::stod(argv[2]) * m; // inward distance
    const auto out_path = fs::path{argv[3]};

    if (offset <= 0.0 * m) throw std::runtime_error{"<offset_m> must be > 0"};

    // ---- Read ORIGINAL perimeter and build polygon
    auto poly_in = ReadPolygon(in_path);
    auto swaths = tjg::BoundarySwaths(poly_in, offset);
    WriteSwaths(out_path, swaths);
  }
  catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} // main
