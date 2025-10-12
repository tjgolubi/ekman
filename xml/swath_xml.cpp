// swatch_xml.cpp
// Parse ISOXML (ISO 11783-10) Guidance swaths: GGP/GPN/LSG/PNT
// C++23, pugixml, exceptions, filesystem paths, fstreams, "almost always auto"

#include <pugixml.hpp>

#include <gsl-lite/gsl-lite.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <optional>
#include <stdexcept>
#include <charconv>
#include <limits>
#include <cmath>
#include <cctype>

namespace fs = std::filesystem;
namespace gsl = gsl_lite;

namespace {

using std::string;
using std::string_view;

using Metres  = double;
using Degrees = double;

//--------------- Small helpers (duplicated for now) ------------------

[[noreturn]] void Fail(string_view msg) {
  throw std::runtime_error{std::string{msg}};
}

string Sanitize(string_view s) {
  auto out = string{};
  out.reserve(s.size());
  for (auto ch : s) {
    auto c = static_cast<unsigned char>(ch);
    auto ok = std::isalnum(c) || ch == '_' || ch == '-' || ch == ' ';
    out.push_back(ok ? ch : '_');
  }
  return out;
} // Sanitize

string RequireAttr(const pugi::xml_node& node, gsl::czstring attrName) {
  if (auto a = node.attribute(attrName); a && a.value()[0] != '\0')
    return string{a.value()};
  Fail(std::string{"missing required attribute '"} + attrName + "' on <"
       + node.name() + ">");
} // RequireAttr

std::optional<string>
OptAttr(const pugi::xml_node& node, gsl::czstring attrName) {
  if (auto a = node.attribute(attrName); a && a.value()[0] != '\0')
    return string{a.value()};
  return std::nullopt;
}

double RequireDouble(const pugi::xml_node& node, gsl::czstring attrName) {
  constexpr auto nan = std::numeric_limits<double>::quiet_NaN();
  auto v = node.attribute(attrName).as_double(nan);
  if (std::isfinite(v))
    return v;
  Fail(std::string{"attribute '"} + attrName + "' on <" + node.name()
       + "> is not a finite number");
} // RequireDouble

std::optional<double>
OptionalDouble(const pugi::xml_node& node, gsl::czstring attrName) {
  constexpr auto nan = std::numeric_limits<double>::quiet_NaN();
  auto a = node.attribute(attrName);
  if (a) {
    auto v = a.as_double(nan);
    if (std::isfinite(v))
      return v;
  }
  return std::nullopt;
} // OptionalDouble

std::vector<char> ReadAll(const fs::path& file) {
  auto ifs = std::ifstream{file, std::ios::binary};
  if (!ifs) Fail("cannot open file: " + file.string());
  ifs.seekg(0, std::ios::end);
  auto size = ifs.tellg();
  if (size < 0) Fail("failed to stat file: " + file.string());
  ifs.seekg(0, std::ios::beg);
  auto buf = std::vector<char>(static_cast<std::size_t>(size));
  if (!ifs.read(buf.data(), static_cast<std::streamsize>(buf.size())))
    Fail("failed to read file: " + file.string());
  return buf;
} // ReadAll

//----------------------- Domain models (lightweight) -----------------------

struct Customer {
  string id;
  string name;
};

struct Farm {
  string id;
  string name;
  std::optional<string> customer_ref;
};

struct PointLLA {
  Degrees lat_deg=0;
  Degrees lon_deg=0;
  std::optional<Metres> alt_m;
};

enum class LineStringType : int {
  PolygonExterior = 1,
  PolygonInterior = 2,
  TramLine        = 3,
  SamplingRoute   = 4,
  GuidancePattern = 5,
  Drainage        = 6,
  Fence           = 7,
  Flag            = 8,
  Obstacle        = 9
};

struct LineString {
  LineStringType type;
  std::optional<string> designator;
  std::optional<string> id; // LSG@F
  std::vector<PointLLA> points;
};

enum class GuidanceType : int { AB=1, APlus=2, Curve=3, Pivot=4, Spiral=5 };

struct GuidancePattern {
  string id;                               // GPN@A
  std::optional<string> designator;        // GPN@B
  std::optional<GuidanceType> kind;        // GPN@C
  std::optional<Degrees> heading_deg;      // GPN@G
  std::optional<LineString> geometry;      // LSG with A=5
};

struct GuidanceGroup {
  string id;                         // GGP@A
  std::optional<string> name;        // GGP@B
  std::vector<GuidancePattern> patterns;
};

struct Partfield {
  string id;                         // PFD@A
  std::optional<string> code;        // PFD@B
  string designator;                 // PFD@C
  std::optional<string> customer_ref;// PFD@E -> CTR@A
  std::optional<string> farm_ref;    // PFD@F -> FRM@A
  std::vector<GuidanceGroup> groups; // GGP children
};

//-------------------------- Parsing --------------------------

std::unordered_map<string, Customer> ParseCustomers(const pugi::xml_node& root) {
  auto map = std::unordered_map<string, Customer>{};
  for (auto ctr : root.children("CTR")) {
    auto id = RequireAttr(ctr, "A");
    auto name = Sanitize(RequireAttr(ctr, "B"));
    map.emplace(id, Customer{ id, name });
  }
  return map;
} // ParseCustomers

std::unordered_map<string, Farm> ParseFarms(const pugi::xml_node& root) {
  auto map = std::unordered_map<string, Farm>{};
  for (auto frm : root.children("FRM")) {
    auto id = RequireAttr(frm, "A");
    auto name = Sanitize(RequireAttr(frm, "B"));
    auto cid = OptAttr(frm, "I"); // CustomerIdRef
    map.emplace(id, Farm{ id, name, cid });
  }
  return map;
} // ParseFarms

LineStringType ToLineStringType(int v) {
  switch (v) {
    case 1: return LineStringType::PolygonExterior;
    case 2: return LineStringType::PolygonInterior;
    case 3: return LineStringType::TramLine;
    case 4: return LineStringType::SamplingRoute;
    case 5: return LineStringType::GuidancePattern;
    case 6: return LineStringType::Drainage;
    case 7: return LineStringType::Fence;
    case 8: return LineStringType::Flag;
    case 9: return LineStringType::Obstacle;
    default: Fail("unknown LSG@A (LineStringType) value: " + std::to_string(v));
  }
} // ToLineStringType

std::optional<GuidanceType> ToGuidanceTypeOpt(const pugi::xml_node& gpn) {
  if (auto v = OptAttr(gpn, "C")) {
    auto iv = 0;
    auto res = std::from_chars(v->data(), v->data() + v->size(), iv);
    if (res.ec == std::errc{}) {
      switch (iv) {
        case 1: return GuidanceType::AB;
        case 2: return GuidanceType::APlus;
        case 3: return GuidanceType::Curve;
        case 4: return GuidanceType::Pivot;
        case 5: return GuidanceType::Spiral;
        default: return std::nullopt;
      }
    }
  }
  return std::nullopt;
} // ToGuidanceTypeOpt

LineString ParseLSG(const pugi::xml_node& lsg) {
  auto type_i = static_cast<int>(RequireDouble(lsg, "A")); // enum stored as NMTOKEN digits
  auto type = ToLineStringType(type_i);
  auto designator = OptAttr(lsg, "B");
  auto id = OptAttr(lsg, "F");
  auto pts = std::vector<PointLLA>{};
  pts.reserve(64);
  for (auto pnt : lsg.children("PNT")) {
    // Per schema, C=lat, D=lon; E=optional altitude
    auto lat = RequireDouble(pnt, "C");
    auto lon = RequireDouble(pnt, "D");
    auto alt = OptionalDouble(pnt, "E");
    pts.push_back(PointLLA{lat, lon, alt});
  }

#if 0
  if (pts.size() < 2)
    Fail("LSG must contain at least two PNT points to form a linestring");
#endif
  return LineString{type, designator, id, std::move(pts)};
} // ParseLSG

GuidancePattern ParseGPN(const pugi::xml_node& gpn) {
  auto id   = RequireAttr(gpn, "A");
  auto name = OptAttr(gpn, "B");
  auto kind = ToGuidanceTypeOpt(gpn);
  auto hdg  = OptionalDouble(gpn, "G");

  // Expect exactly one LSG (schema allows one); we’ll take the first if multiple exist
  std::optional<LineString> geom;
  for (auto lsg : gpn.children("LSG")) {
    auto ls = ParseLSG(lsg);
    geom = std::move(ls);
    break;
  }

  return GuidancePattern{id, name, kind, hdg, std::move(geom)};
} // ParseGPN

GuidanceGroup ParseGGP(const pugi::xml_node& ggp) {
  auto id = RequireAttr(ggp, "A");
  auto name = OptAttr(ggp, "B");
  auto patterns = std::vector<GuidancePattern>{};
  patterns.reserve(8);
  for (auto gpn : ggp.children("GPN"))
    patterns.push_back(ParseGPN(gpn));
  if (patterns.empty()) {
    // Schema permits zero? (We’ll tolerate but warn via exception?)
  }
  return GuidanceGroup{ id, name, std::move(patterns) };
} // ParseGGP

Partfield ParsePFD(const pugi::xml_node& pfd) {
  auto id   = RequireAttr(pfd, "A");
  auto code = OptAttr(pfd, "B");
  auto des  = RequireAttr(pfd, "C");
  auto ctr  = OptAttr(pfd, "E"); // CustomerIdRef
  auto frm  = OptAttr(pfd, "F"); // FarmIdRef

  auto groups = std::vector<GuidanceGroup>{};
  for (auto ggp : pfd.children("GGP"))
    groups.push_back(ParseGGP(ggp));

  return Partfield{ id, code, Sanitize(des), ctr, frm, std::move(groups) };
} // ParsePFD

struct Index {
  std::unordered_map<string, Customer> customers;
  std::unordered_map<string, Farm>     farms;
  std::unordered_map<string, Partfield> partfields;
};

Index BuildIndex(const pugi::xml_node& root) {
  auto out = Index{};
  out.customers  = ParseCustomers(root);
  out.farms      = ParseFarms(root);
  for (auto pfd : root.children("PFD")) {
    auto pf = ParsePFD(pfd);
    auto id = pf.id;
    out.partfields.emplace(id, std::move(pf));
  }
  return out;
} // BuildIndex

//-------------------------- Reporting --------------------------

void PrintSummary(const Index& idx) {
  auto totalGGP = std::size_t{0};
  auto totalGPN = std::size_t{0};
  auto totalLSG = std::size_t{0};
  auto totalPNT = std::size_t{0};

  // Reference checks (customer/farm)
  for (const auto& [id, pf] : idx.partfields) {
    if (pf.customer_ref && !idx.customers.contains(*pf.customer_ref)) {
      std::cerr << "warning: PFD " << id << " references missing CTR "
                << *pf.customer_ref << '\n';
    }
    if (pf.farm_ref && !idx.farms.contains(*pf.farm_ref)) {
      std::cerr << "warning: PFD " << id << " references missing FRM "
                << *pf.farm_ref << '\n';
    }
  }

  std::cout << "Customers: " << idx.customers.size() << '\n';
  for (const auto& [id, c] : idx.customers)
    std::cout << "  CTR  : " << c.name << " (id=" << id << ")\n";
  std::cout << "Farms: " << idx.farms.size() << '\n';
  for (const auto& [id, f] : idx.farms) {
    std::cout << "  FRM  : " << f.name << " (id=" << id << ", CTR="
              << (f.customer_ref ? *f.customer_ref : "?") << ")\n";
  }

  std::cout << "Fields: " << idx.partfields.size() << '\n';
  for (const auto& [pfdId, pf] : idx.partfields) {
    std::cout << "  PFD " << pfdId << " : \"" << pf.designator << "\""
              << " (CTR=" << (pf.customer_ref ? *pf.customer_ref : "?")
              << ", FRM=" << (pf.farm_ref ? *pf.farm_ref : "?") << ")\n";

    for (const auto& ggp : pf.groups) {
      ++totalGGP;
      std::cout << "    GGP " << ggp.id
                << (ggp.name ? " \"" + *ggp.name + "\"" : "") << '\n';
      for (const auto& gpn : ggp.patterns) {
        ++totalGPN;
        std::cout << "      GPN " << gpn.id;
        if (gpn.designator) std::cout << " \"" << *gpn.designator << "\"";
        if (gpn.kind) std::cout << " kind=" << static_cast<int>(*gpn.kind);
        if (gpn.heading_deg) std::cout << " heading=" << *gpn.heading_deg << " deg";
        if (gpn.geometry) {
          ++totalLSG;
          auto& ls = *gpn.geometry;
          std::cout << " LSG(type=" << static_cast<int>(ls.type)
                    << ", pts=" << ls.points.size() << ")";
          totalPNT += ls.points.size();
        }
        std::cout << '\n';
      }
    }
  }

  std::cout << "Totals:\n";
  std::cout << "  GuidanceGroups:  " << totalGGP << '\n';
  std::cout << "  Patterns:        " << totalGPN << '\n';
  std::cout << "  LineStrings:     " << totalLSG << '\n';
  std::cout << "  Points:          " << totalPNT << '\n';
} // PrintSummary

//-------------------------- Main --------------------------

int run(int argc, gsl::czstring argv[]) {
  auto path = (argc >= 2) ? fs::path{argv[1]} : fs::path{"swaths.xml"};
  if (!fs::exists(path))
    Fail("file not found: " + path.string());

  auto buffer = ReadAll(path);

  auto doc = pugi::xml_document{};
  auto parseResult =
                doc.load_buffer(buffer.data(), buffer.size(), pugi::parse_full);
  if (!parseResult)
    Fail(std::string{"XML parse error: "} + parseResult.description());

  auto root = doc.child("ISO11783_TaskData");
  if (!root) {
    // Some exports use ISO11783_TaskFile; tolerate if root not found
    root = doc.first_child();
    if (!root) Fail("no root element found");
  }

  auto idx = BuildIndex(root);
  PrintSummary(idx);
  return EXIT_SUCCESS;
} // run

} // local

int main(int argc, gsl::czstring argv[]) {
  try {
    return run(argc, argv);
  }
  catch (const std::exception& e) {
    std::cerr << "std::exception: " << e.what() << '\n';
  }
  catch (...) {
    std::cerr << "unknown exception\n";
  }
  return EXIT_FAILURE;
} // main
