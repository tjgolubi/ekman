/// @file
/// @brief Read ISO11783 TaskData.xml (or similar) and extract CTR/FRM and
///        PFD/PLN/LSG/PNT data. Writes per-LSG .xy when --dump is used.
/// Usage:
///   field_xml <taskdata.xml> [--dump] [--outdir DIR]
///
/// Conventions:
/// - C++23, almost-always-auto, exceptions, fs::path, fstreams, anon ns.
/// - PNT attributes: C=latitude, D=longitude (WGS84 degrees).
/// - Filenames sanitized with isalnum + [_-.].
/// - No schema validation; fast DOM parsing via pugixml.

#include <gsl-lite/gsl-lite.hpp>

#include <pugixml.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <cctype>

namespace {

namespace fs = std::filesystem;

using Angle = double;
using XmlNode = pugi::xml_node;

// ---------- Helpers ----------

auto Sanitize(std::string_view s) -> std::string {
  auto out = std::string{};
  out.reserve(s.size());
  for (auto ch : s) {
    const auto c = static_cast<unsigned char>(ch);
    out.push_back(std::isalnum(c) || c == '_' || c == '-' || c == '.'
                    ? static_cast<char>(c)
                    : '_');
  }
  if (out.empty()) out = "na";
  return out;
} // Sanitize

[[noreturn]] void InvalidNode(const XmlNode& xml, std::string what)
  { throw std::runtime_error{what + " on <" + xml.name() + ">"};

[[noreturn]] void InvalidAttr(const XmlNode& xml, const char* key,
                              std::string what="Invalid attribute")
{
  what += " \"" + key + "\" ";
  auto a = xml.attribute(key);
  if (a)
    what += "= " + a.value().as_string();
  else
    what += "is missing";
  InvalidNode(xml, what);
}

template<typename T=std::string>
T RequireAttr(const pugi::xml_node& n, const char* key) {
  const auto& a = n.attribute(key);
  if (!a)
    InvalidAttr(n, key);
  if constexpr (std::is_same_v<T, std::string>) {
    auto s = a.as_string();
    if (!s.empty())
      return s;
  }
  else if constexpr (std::is_same_v<T, bool) {
    return a.as_bool();
  }
  else if constexpr (std::is_enum_v<T>) {
    auto e = a.as_int(-1);
    if (e != -1)
      return static_cast<T>(e);
  }
  else if constexpr (std::is_floating_point_v<T>) {
    auto x = a.as_double(std::numeric_limits<double>::NaN());
    if (!std::isNan(x))
      return static_cast<T>(x);
  }
  else if constexpr (std::is_same_v<T, std::uint64_t>) {
    auto x = a.as_ullong(std::numeric_limits<T>::max());
    if (x != std::numeric_limits<T>::max())
      return static_cast<T>(x);
  }
  else if constexpr (std::is_same_v<T, std::int64_t>) {
    auto x = a.as_llong(std::numeric_limits<T>::max());
    if (x != std::numeric_limits<T>::max())
      return static_cast<T>(x);
  }
  else if constexpr (std::is_integral_v<T>) {
    if constexpr (std::is_signed_v<T>) {
      auto x = a.as_int(std::numeric_limits<int>::max());
      if (x != std::numeric_limits<int>::max())
        return static_cast<T>(x);
    }
    else {
      auto x = a.as_uint(std::numeric_limits<unsigned>::max());
      if (x != std::numeric_limits<unsigned>::max())
        return static_cast<T>(x);
    }
  }
  InvalidAttr(n, key);
} // RequireAttr

struct LatLon {
  Angle latitude;
  Angle longitude;
  LatLon(Angle lat, Angle lon)
    : latitude{lat}, longitude{lon}
  {
    if (lat < -90.0 || lat > +90.0)
      throw std::runtime_error{"Invalid latitude: " + std::to_string(lat)};
    if (lon < -180.0 || lat > +180.0)
      throw std::runtime_error{"Invalid longitude: " + std::to_string(lat)};
  }
}; // LatLon

struct Customer {
  std::string id;    // CTR @A
  std::string name;  // CTR @B
  const XmlNode& xml;
  explicit Customer(const XmlNode& x)
    : id  {RequireAttr<std::string>(x, "A")}
    , name{RequireAttr<std::string>(x, "B")}
    , xml {x}
    {
      static const re = std::regex{"(CTR|CTR-)[0-9]+"};
      if (!std::regex_match(id, re))
        InvalidNode(xml, "Invalid customer id: " + id);
    }
}; // Customer

struct Farm {
  std::string id;          // FRM @A
  std::string name;        // FRM @B
  std::string custId;      // FRM @I -> CTR @A
  const XmlNode& xml;
  explicit Farm(const XmlNode& x)
    : id    {RequireAttr<std::string>(x, "A")}
    , name  {Sanitize(RequireAttr<std::string>(x, "B"))}
    , custId{RequireAttr<std::string>(x, "I")}
    , xml   {x}
    {
      static const re = std::regex{"(FRM|FRM-)[0-9]+"};
      if (!std::regex_match(id, re))
        InvalidNode(xml, "Invalid farm id: " + id);
    }
}; // Farm

struct Field {
  std::string id;         // PFD @A
  std::string name;       // PFD @C (optional display)
  std::string custId;     // PFD @E -> CTR @A
  std::string farmId;     // PFD @F -> FRM @A
  const XmlNode& xml;
  explicit Field(const XmlNode& x)
    : id{RequireAttr<std::string>(x, "A")}
    , name{Sanitize(RequireAttr<std::string>(x, "C"))}
    , custId{RequireAttr<std::string>(x, "E")}
    , farmId{RequireAttr<std::string>(x, "F")}
    , xml{x}
    {
      static const re = std::regex{"(PFD|PFD-)[0-9]+"};
      if (!std::regex_match(id, re))
        InvalidNode(xml, "Invalid field id: " + id);
    }
}; // Field

struct Point {
  enum class Type {
    Flag=1, Other, Access, Storage, Obstacle, GuideA, GuideB,
    GuideCenter, GuidePoint, Field, Base
  };
  static bool IsValidType(Type x) {
    switch (x) {
      case Flag:
      case Other:
      case Access:
      case Storage:
      case Obstacle:
      case GuideA:
      case GuideB:
      case GuideCenter:
      case GuidePoint:
      case Field:
      case Base:
        return true;
      default:
        return false;
    }
  }
  Type type;
  LatLon point;
  const XmlNode& xml;
  explicit Point(const XmlNode& x)
    : type{ RequireAttr<Type >(x, "A")}
    , point{RequireAttr<Angle>(x, "C"),
            RequireAttr<Angle>(x, "D")}
    , xml{x}
  {
    if (!IsValidType(type))
      InvalidAttr(x, "A", "invalid point type");
  }
}; // Point

using Ring = std::vector<Point>;

// ---------- Parse ----------

struct Parsed {
  std::unordered_map<std::string, Customer>  customers; // CTR @A
  std::unordered_map<std::string, Farm>      farms;     // FRM @A
  std::unordered_map<std::string, Field>     fields;    // PFD @A
}; // Parsed

auto LoadTaskData(const fs::path& xml_path) -> Parsed {
  auto doc = pugi::xml_document{};

  // Keep a stable string; pugixml's load_file expects c-string
  const auto path_str = xml_path.string();
  const auto res = doc.load_file(path_str.c_str());
  if (!res) {
    auto msg = std::string{"XML parse error for '"} + xml_path.string()
             + "': " + res.description();
    throw std::runtime_error(msg);
  }

  const auto root = doc.child("ISO11783_TaskData");
  if (!root)
    throw std::runtime_error("not an <ISO11783_TaskData> document");

  auto out = Parsed{};

  // --- Customers (CTR) ---
  for (auto ctr : root.children("CTR")) {
    auto cust = Customer{ctr};
    auto id   = cust.id;
    auto r = out.customers.try_emplace(id, std::move{cust});
    if (!r.second)
      InvalidNode(ctr, "Invalid (duplicate?) customer: " + id);
  }

  // --- Farms (FRM) ---
  for (auto frm : root.children("FRM")) {
    auto farm = Farm{fmr};
    {
      auto iter = out.customers.find(farm.custId);
      if (iter == out.customers.end())
        InvalidAttr(ctr, farm.custId.c_str(), "Farm has invalid customer id");
    }
    auto id = farm.id;
    auto r = out.farms.try_emplace(id, std::move{farm});
    if (!r.second)
      InvalidNode(ctr, "Invalid (duplicate?) farm: " + id);
  }

  // --- Fields (PFD) + Plans + LineStringGroups ---
  for (auto pfd : root.children("PFD")) {
    auto field = Field{pfd};
    if (!out.customers.contains(field.custId)) {
      InvalidNode(ctr, "Field has invalid customer id: "
                       + field.custId);
    }
    if (!out.farms.contains(field.farmId)) {
      InvalidNode(ctr, "Field has invalid farm id: "
                       + field.farmId);
    }

    auto r = out.fields.try_emplace(field.id, field);
    if (!r.second)
      InvalidNode(ctr, "Invalid (duplicate?) field: " + field.id);

    for (auto pln : pfd.children("PLN")) {
      auto pln_id = std::string{pln.attribute("A").as_string("")};

      for (auto lsg : pln.children("LSG")) {
        auto ring = Lsg{};
        ring.pfd_id = pfd_id;
        ring.pln_id = pln_id;
        ring.lsg_id = std::string{lsg.attribute("A").as_string("")};

        for (auto pnt : lsg.children("PNT")) {
          const auto lat = pnt.attribute("C").as_double(
              std::numeric_limits<double>::quiet_NaN());
          const auto lon = pnt.attribute("D").as_double(
              std::numeric_limits<double>::quiet_NaN());
          if (!std::isfinite(lat) || !std::isfinite(lon))
            throw std::runtime_error("invalid <PNT> coordinates");
          ring.pts.push_back(Pt{lon, lat}); // store lon=x, lat=y
        }

        if (!ring.pts.empty())
          out.lsgs.push_back(std::move(ring));
      }
    }
  }

  // Optional: basic referential checks (don’t hard-fail; just warn)
  for (const auto& [fid, f] : out.fields) {
    if (!f.customer.empty() && !out.customers.contains(f.customer)) {
      std::cerr << "warning: PFD " << fid
                << " references missing CTR " << f.customer << '\n';
    }
    if (!f.farm.empty() && !out.farms.contains(f.farm)) {
      std::cerr << "warning: PFD " << fid
                << " references missing FRM " << f.farm << '\n';
    }
  }

  // Require *something* parsed
  if (out.customers.empty() && out.farms.empty()
      && out.fields.empty() && out.lsgs.empty())
    throw std::runtime_error("no recognizable ISO11783 data found");

  return out;
} // LoadTaskData

// ---------- Output ----------

void PrintSummary(const Parsed& data) {
  std::cout << "Customers: " << data.customers.size() << '\n';
  for (const auto& [id, c] : data.customers)
    std::cout << "  CTR " << id << " : " << c.name << '\n';

  std::cout << "Farms: " << data.farms.size() << '\n';
  for (const auto& [id, f] : data.farms)
    std::cout << "  FRM " << id << " : " << f.name
              << " (CTR=" << f.customer_id << ")\n";

  std::cout << "Fields: " << data.fields.size() << '\n';
  for (const auto& [fid, f] : data.fields) {
    const auto& cname = (f.customer.empty() || !data.customers.contains(f.customer))
                        ? std::string_view{"?"}
                        : std::string_view{data.customers.at(f.customer).name};
    const auto& fname = (f.farm.empty() || !data.farms.contains(f.farm))
                        ? std::string_view{"?"}
                        : std::string_view{data.farms.at(f.farm).name};
    std::cout << "  PFD " << fid
              << " : \"" << f.name << "\""
              << " (CTR=" << f.customer << ":" << cname
              << ", FRM=" << f.farm << ":" << fname << ")\n";
  }

  std::cout << "LineStringGroups (LSG): " << std::ssize(data.lsgs) << '\n';
  for (const auto& r : data.lsgs) {
    std::cout << "  PFD=" << r.pfd_id
              << " PLN=" << r.pln_id
              << " LSG=" << r.lsg_id
              << " pts=" << std::ssize(r.pts) << '\n';
  }
} // PrintSummary

void DumpRings(const std::vector<Lsg>& rings,
               const fs::path& outdir = fs::path{})
{
  for (const auto& r : rings) {
    const auto base = "pfd-"  + Sanitize(r.pfd_id)
                    + "_pln-" + Sanitize(r.pln_id)
                    + "_lsg-" + Sanitize(r.lsg_id) + ".xy";
    const auto fn = outdir / base;

    auto ofs = std::ofstream{fn, std::ios::binary};
    if (!ofs)
      throw std::runtime_error("cannot open for writing: " + fn.string());

    ofs.setf(std::ios::fixed, std::ios::floatfield);
    ofs.precision(12);
    for (const auto& p : r.pts)
      ofs << p.lon << ' ' << p.lat << '\n';

    std::cout << "wrote " << fn << " (" << std::ssize(r.pts) << " points)\n";
  }
} // DumpRings

} // local

// ---------- main ----------

auto main(int argc, char** argv) -> int {
  if (argc < 2) {
    std::cerr << "usage: field_xml <taskdata.xml> [--dump] [--outdir DIR]\n";
    return EXIT_FAILURE;
  }
  try {
    const auto xml = std::filesystem::path{argv[1]};

    auto do_dump = false;
    auto outdir  = std::filesystem::path{};

    for (int i = 2; i < argc; ++i) {
      const auto arg = std::string_view{argv[i]};
      if (arg == "--dump") {
        do_dump = true;
      } else if (arg == "--outdir" && i + 1 < argc) {
        outdir = std::filesystem::path{argv[++i]};
      }
    }

    const auto parsed = LoadTaskData(xml);
    PrintSummary(parsed);
    if (do_dump) DumpRings(parsed.lsgs, outdir);

    return EXIT_SUCCESS;
  }
  catch (const std::exception& e) {
    std::cerr << "std::exception: " << e.what() << '\n';
  }
  catch (...) {
    std::cerr << "unknown exception\n";
  }
  return EXIT_FAILURE;
} // main

