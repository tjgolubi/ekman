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
#include <unordered_map>
#include <regex>
#include <string>
#include <string_view>
#include <vector>
#include <optional>
#include <iterator>
#include <functional>
#include <concepts>
#include <utility>
#include <stdexcept>
#include <limits>
#include <type_traits>
#include <cstring>
#include <cmath>
#include <cctype>

namespace {

namespace fs = std::filesystem;
namespace gsl = gsl_lite;

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

template<typename T>
concept Enum = std::is_enum_v<T>;

template<Enum E> constexpr E enum_cast(int value);

[[noreturn]] void InvalidNode(const XmlNode& xml, std::string what) {
  what.reserve(what.size() + 8 + std::strlen(xml.name()));
  what += " on <";
  what += xml.name();
  what += ">";
  throw std::runtime_error{what};
} // InvalidNode

[[noreturn]] void InvalidAttr(const XmlNode& xml, const char* key,
                              std::string what="Invalid attribute")
{
  what.reserve(what.size() + std::strlen(key) + 40);
  what += " \"";
  what += key;
  what += "\" ";
  auto a = xml.attribute(key);
  if (a) {
    what += "= ";
    what += a.as_string();
  }
  else {
    what += "is missing";
  }
  InvalidNode(xml, what);
} // InvalidAttr

template<typename T>
std::optional<T> GetAttr(const XmlNode& x, const char* key) {
  const auto& a = x.attribute(key);
  if (!a)
    return std::nullopt;
  if constexpr (std::is_same_v<T, std::string>) {
    return a.as_string();
  }
  else if constexpr (std::is_same_v<T, bool>) {
    return a.as_bool();
  }
  else if constexpr (std::is_enum_v<T>) {
    auto e = a.as_int(-1);
    if (e != -1)
      return enum_cast<T>(e);
  }
  else if constexpr (std::is_floating_point_v<T>) {
    auto x = a.as_double(std::numeric_limits<double>::infinity());
    if (std::isfinite(x))
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
  else {
    static_assert(false, "GetAttr: invalid type");
  }
  InvalidAttr(x, key);
} // GetAttr

template<typename T=std::string>
T RequireAttr(const pugi::xml_node& x, const char* key) {
  const auto& a = x.attribute(key);
  if (!a)
    InvalidAttr(x, key);
  auto v = GetAttr<T>(x, key);
  if (v)
    return *v;
  if constexpr (std::is_same_v<T, std::string>) {
    auto s = a.as_string();
    if (s[0] != '\0')
      return s;
  }
  InvalidAttr(x, key);
} // RequireAttr

struct LatLon {
  Angle latitude;
  Angle longitude;
  LatLon(Angle lat, Angle lon)
    : latitude{lat}, longitude{lon}
  {
    if (lat < -90.0 || lat > +90.0)
      throw std::runtime_error{"Invalid latitude: " + std::to_string(lat)};
    if (lon < -180.0 || lon > +180.0)
      throw std::runtime_error{"Invalid longitude: " + std::to_string(lon)};
  }
}; // LatLon

struct Customer {
  std::string id;    // CTR @A
  std::string name;  // CTR @B
  const XmlNode* xml = nullptr;
  Customer() = default;
  explicit Customer(const XmlNode& x)
    : id  {RequireAttr<std::string>(x, "A")}
    , name{RequireAttr<std::string>(x, "B")}
    , xml {&x}
    {
      static const auto re = std::regex{"(CTR|CTR-)[0-9]+"};
      if (!std::regex_match(id, re))
        InvalidNode(x, "Invalid customer id: " + id);
    }
}; // Customer

struct Farm {
  std::string id;          // FRM @A
  std::string name;        // FRM @B
  std::string custId;      // FRM @I -> CTR @A
  const XmlNode* xml = nullptr;
  Farm() = default;
  explicit Farm(const XmlNode& x)
    : id    {RequireAttr<std::string>(x, "A")}
    , name  {Sanitize(RequireAttr<std::string>(x, "B"))}
    , custId{RequireAttr<std::string>(x, "I")}
    , xml   {&x}
    {
      static const auto re = std::regex{"(FRM|FRM-)[0-9]+"};
      if (!std::regex_match(id, re))
        InvalidNode(x, "Invalid farm id: " + id);
    }
}; // Farm

struct Point {
  enum class Type {
    Flag=1, Other, Access, Storage, Obstacle, GuideA, GuideB,
    GuideCenter, GuidePoint, Field, Base
  };
  static constexpr const char* Name(Type x) {
    switch (x) {
      case Type::Flag:        return "Type::Flag";
      case Type::Other:       return "Type::Other";
      case Type::Access:      return "Type::Access";
      case Type::Storage:     return "Type::Storage";
      case Type::Obstacle:    return "Type::Obstacle";
      case Type::GuideA:      return "Type::GuideA";
      case Type::GuideB:      return "Type::GuideB";
      case Type::GuideCenter: return "Type::GuideCenter";
      case Type::GuidePoint:  return "Type::GuidePoint";
      case Type::Field:       return "Type::Field";
      case Type::Base:        return "Type::Base";
      default: return nullptr;
    }
  } // Name(Type)
  Type type;
  LatLon point;
  const XmlNode* xml;
  using TypeValidator = std::function<bool(Type)>;
  static constexpr bool AnyType(Type) { return true; }
  Point() = default;
  explicit Point(const XmlNode& x, TypeValidator validator=AnyType)
    : type{ RequireAttr<Type >(x, "A")}
    , point{RequireAttr<Angle>(x, "C"),
            RequireAttr<Angle>(x, "D")}
    , xml{&x}
    {
      if (!validator(type)) {
        auto msg = std::string{"Point: invalid type in this context: "}
                 + Name(type);
        throw std::range_error{msg};
      }
    }
}; // Point

template<> constexpr Point::Type enum_cast<Point::Type>(int x) {
  auto e = static_cast<Point::Type>(x);
  if (Point::Name(e)) [[likely]]
    return e;
  throw std::range_error{"Invalid Point::Type: " + std::to_string(x)};
} // enum_cast

struct LineString { // LSG
  enum class Type {
    Exterior=1,
    Interior,
    TramLine,
    Sampling,
    Guidance,
    Drainage,
    Fence,
    Flag,
    Obstacle
  };
  static constexpr const char* Name(Type x) {
    switch (x) {
      case Type::Exterior:  return "Type::Exterior";
      case Type::Interior:  return "Type::Interior";
      case Type::TramLine:  return "Type::TramLine";
      case Type::Sampling:  return "Type::Sampling";
      case Type::Guidance:  return "Type::Guidance";
      case Type::Drainage:  return "Type::Drainage";
      case Type::Fence:     return "Type::Fence";
      case Type::Flag:      return "Type::Flag";
      case Type::Obstacle:  return "Type::Obstacle";
      default: return nullptr;
    }
  } // Name(Type)
  Type type;
  std::vector<Point> points;
  const XmlNode* xml = nullptr;

  bool empty() const { return points.empty(); }
  std::size_t size() const { return points.size(); }

  LineString() = default;
  explicit LineString(const XmlNode& x,
                      Point::TypeValidator point_validator=Point::AnyType)
    : type{RequireAttr<Type>(x, "A")}
    , xml{&x}
  {
    for (const auto& pnt: x.children("PNT"))
      points.emplace_back(Point{pnt, point_validator});
  }
}; // LineString

template<> constexpr LineString::Type enum_cast<LineString::Type>(int x) {
  auto e = static_cast<LineString::Type>(x);
  if (LineString::Name(e)) [[likely]]
    return e;
  throw std::range_error{"Invalid LineString type: " + std::to_string(x)};
}; // enum_cast

struct Polygon { // PLN
  enum class Type {
    Boundary=1,
    Treatment,
    Water,
    Building,
    Road,
    Obstacle,
    Flag,
    Other,
    Field,
    Headland,
    Buffer,
    Windbreak
  }; // Type

  static constexpr const char* Name(Type x) {
    switch (x) {
      case Type::Boundary:  return "Type::Boundary";
      case Type::Treatment: return "Type::Treatment";
      case Type::Water:     return "Type::Water";
      case Type::Building:  return "Type::Building";
      case Type::Road:      return "Type::Road";
      case Type::Obstacle:  return "Type::Obstacle";
      case Type::Flag:      return "Type::Flag";
      case Type::Other:     return "Type::Other";
      case Type::Field:     return "Type::Field";
      case Type::Headland:  return "Type::Headland";
      case Type::Buffer:    return "Type::Buffer";
      case Type::Windbreak: return "Type::Windbreak";
      default: return nullptr;
    }
  } // Name(Type)

  Type type;
  LineString outer;
  std::vector<LineString> inners;
  const XmlNode* xml = nullptr;

  Polygon() = default;
  explicit Polygon(const XmlNode& x)
    : type{RequireAttr<Type>(x, "A")}
    , xml{&x}
  {
    auto point_validator = [](Point::Type t) -> bool {
      return (t == Point::Type::Field);
    };
    for (const auto& lsg: x.children("LSG")) {
      auto ring = LineString{lsg, point_validator};
      switch (ring.type) {
        case LineString::Type::Exterior:
          if (!outer.empty())
            throw std::runtime_error{"Polygon: multiple exterior rings"};
          outer = std::move(ring);
          break;
        case LineString::Type::Interior:
          inners.emplace_back(std::move(ring));
          break;
        default: {
          auto msg = std::string{"Polygon: unexpected LineString type: "}
                   + LineString::Name(ring.type);
          throw std::runtime_error{msg};
        }
      }
    }
    if (outer.empty())
      throw std::runtime_error{"Polygon: missing exterior ring"};
    if (std::ssize(outer) < 4)
      throw std::runtime_error{"Polygon: exterior ring too small"};
    for (const auto& r: inners) {
      if (std::ssize(r) < 4)
        throw std::runtime_error{"Polygon: inter ring too small"};
    }
  }
}; // Polygon

template<> Polygon::Type enum_cast<Polygon::Type>(int x) {
  auto e = static_cast<Polygon::Type>(x);
  if (Polygon::Name(e)) [[likely]]
    return e;
  throw std::range_error{"Invalid Polygon type: " + std::to_string(x)};
} // enum_cast

struct Field {
  std::string id;         // PFD @A
  std::string name;       // PFD @C (optional display)
  std::string custId;     // PFD @E -> CTR @A
  std::string farmId;     // PFD @F -> FRM @A
  std::vector<Polygon> parts;
  const XmlNode* xml = nullptr;
  Field() = default;
  explicit Field(const XmlNode& x)
    : id{RequireAttr<std::string>(x, "A")}
    , name{Sanitize(RequireAttr<std::string>(x, "C"))}
    , custId{RequireAttr<std::string>(x, "E")}
    , farmId{RequireAttr<std::string>(x, "F")}
    , xml{&x}
  {
    static const auto re = std::regex{"(PFD|PFD-)[0-9]+"};
    if (!std::regex_match(id, re))
      InvalidNode(x, "Invalid field id: " + id);
    for (const auto& pln: x.children("PLN"))
      parts.emplace_back(Polygon{pln});
  }
}; // Field

// ---------- Parse ----------

struct TaskData {
  std::unordered_map<std::string, Customer>  customers; // CTR @A
  std::unordered_map<std::string, Farm>      farms;     // FRM @A
  std::unordered_map<std::string, Field>     fields;    // PFD @A
}; // TaskData

auto LoadTaskData(const fs::path& xml_path) -> TaskData {
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

  auto out = TaskData{};

  // --- Customers (CTR) ---
  for (auto ctr : root.children("CTR")) {
    auto cust = Customer{ctr};
    auto id   = cust.id;
    auto r = out.customers.try_emplace(id, std::move(cust));
    if (!r.second)
      InvalidNode(ctr, "Invalid (duplicate?) customer: " + id);
  }

  // --- Farms (FRM) ---
  for (auto frm : root.children("FRM")) {
    auto farm = Farm{frm};
    {
      if (!out.customers.contains(farm.custId))
        InvalidNode(frm, "Farm has invalid customer id: " + farm.custId);
    }
    auto id = farm.id;
    auto r = out.farms.try_emplace(id, std::move(farm));
    if (!r.second)
      InvalidNode(frm, "Invalid (duplicate?) farm: " + id);
  }

  // --- Fields (PFD) + Plans + LineStringGroups ---
  for (auto pfd : root.children("PFD")) {
    auto field = Field{pfd};
    if (!out.customers.contains(field.custId))
      InvalidNode(pfd, "Field has invalid customer id: " + field.custId);
    if (!out.farms.contains(field.farmId))
      InvalidNode(pfd, "Field has invalid farm id: " + field.farmId);

    auto id = field.id;
    auto r = out.fields.try_emplace(id, std::move(field));
    if (!r.second)
      InvalidNode(pfd, "Invalid (duplicate?) field: " + id);
  }

  // Optional: basic referential checks (don’t hard-fail; just warn)
  for (const auto& [fid, f] : out.fields) {
    if (!f.custId.empty() && !out.customers.contains(f.custId)) {
      std::cerr << "warning: PFD " << fid
                << " references missing CTR " << f.custId << '\n';
    }
    if (!f.farmId.empty() && !out.farms.contains(f.farmId)) {
      std::cerr << "warning: PFD " << fid
                << " references missing FRM " << f.farmId << '\n';
    }
  }

  // Require *something* parsed
  if (out.customers.empty() && out.farms.empty() && out.fields.empty())
    throw std::runtime_error("no recognizable ISO11783 data found");

  return out;
} // LoadTaskData

// ---------- Output ----------

void PrintSummary(const TaskData& data) {
  std::cout << "Customers: " << std::ssize(data.customers) << '\n';
  for (const auto& [id, c] : data.customers)
    std::cout << "  CTR " << id << " : " << c.name << '\n';

  std::cout << "Farms: " << std::ssize(data.farms) << '\n';
  for (const auto& [id, f] : data.farms) {
    std::cout << "  FRM " << id << " : " << f.name
              << " (CTR=" << f.custId << ")\n";
  }

  std::cout << "Fields: " << std::ssize(data.fields) << '\n';
  for (const auto& [fid, f] : data.fields) {
    auto custName = std::string_view{"?"};
    auto custIter = data.customers.find(f.custId);
    if (custIter != data.customers.end())
      custName = custIter->second.name;
    auto farmName = std::string_view{"?"};
    auto farmIter = data.farms.find(f.farmId);
    if (farmIter != data.farms.end())
      farmName = farmIter->second.name;

    std::cout << "  PFD " << fid
              << " : \"" << f.name << "\""
              << " (CTR=" << f.custId << ":" << custName
              << ", FRM=" << f.farmId << ":" << farmName << ")\n";
    if (f.parts.empty())
      std::cout << "    No polygons";
    for (const auto& p: f.parts) {
      std::cout << "    Polygon: " << Polygon::Name(p.type)
                << " outer=" << std::ssize(p.outer) << " pts";
      auto innerCount = std::ssize(p.inners);
      if (innerCount > 0) {
        std::cout << ", " << innerCount << " inner rings";
        auto sum = std::size_t{0};
        for (const auto& r: p.inners)
          sum += std::ssize(r);
        std::cout << ", total=" << sum << " points";
      }
      std::cout << std::endl;

    }
  }

} // PrintSummary

} // local

// ---------- main ----------

auto main(int argc, char** argv) -> int {
  if (argc < 2) {
    std::cerr << "usage: field_xml <taskdata.xml>\n";
    return EXIT_FAILURE;
  }
  try {
    const auto xml = fs::path{argv[1]};

    const auto parsed = LoadTaskData(xml);
    PrintSummary(parsed);

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

